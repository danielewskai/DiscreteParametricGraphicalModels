# Packages ----------------------------------------------------
library(spatstat.data)
library(spatstat.geom)
library(dplyr)
library(igraph)


# Data --------------------------------------------------------
data(amacrine)
dataset_name <- "amacrine_on_C4_2x2"

pp <- amacrine[marks(amacrine) == "on"]
pp <- unmark(pp)

plot(pp, main = "Amacrine ON cells")


# Normalize the observation window to [0, 1]^2 ----------------
W <- Window(pp)

pp_unit <- ppp(
  x = (pp$x - W$xrange[1]) / diff(W$xrange),
  y = (pp$y - W$yrange[1]) / diff(W$yrange),
  window = owin(c(0, 1), c(0, 1))
)

plot(pp_unit, main = "Amacrine ON cells, normalized")


# Build binary observations ----------------------------------
#
# We divide the full window into 12 x 12 blocks.
# Then each block is divided locally into 2 x 2 cells:
#
#   v3  v4
#   v1  v2
#
# Each block gives one binary vector X_i in {0,1}^4.
n_blocks_x <- 12
n_blocks_y <- 12

make_2x2_observations <- function(pp_unit, n_blocks_x, n_blocks_y) {
  block_width <- 1 / n_blocks_x
  block_height <- 1 / n_blocks_y
  
  X <- matrix(
    0,
    nrow = n_blocks_x * n_blocks_y,
    ncol = 4
  )
  
  block_info <- data.frame(
    block_id = seq_len(n_blocks_x * n_blocks_y),
    bx = NA_integer_,
    by = NA_integer_,
    n_points = NA_integer_,
    n_occupied_regions = NA_integer_
  )
  
  k <- 1
  
  for (bx in seq_len(n_blocks_x)) {
    for (by in seq_len(n_blocks_y)) {
      
      x0 <- (bx - 1) * block_width
      x1 <- bx * block_width
      y0 <- (by - 1) * block_height
      y1 <- by * block_height
      
      inside_x <- if (bx < n_blocks_x) {
        pp_unit$x >= x0 & pp_unit$x < x1
      } else {
        pp_unit$x >= x0 & pp_unit$x <= x1
      }
      
      inside_y <- if (by < n_blocks_y) {
        pp_unit$y >= y0 & pp_unit$y < y1
      } else {
        pp_unit$y >= y0 & pp_unit$y <= y1
      }
      
      inside <- inside_x & inside_y
      
      x_local <- (pp_unit$x[inside] - x0) / block_width
      y_local <- (pp_unit$y[inside] - y0) / block_height
      
      if (length(x_local) > 0) {
        ix <- pmin(2, floor(2 * x_local) + 1)
        iy <- pmin(2, floor(2 * y_local) + 1)
        
        # Numbering:
        #
        #   v3  v4
        #   v1  v2
        #
        region <- ifelse(
          ix == 1 & iy == 1, 1,
          ifelse(
            ix == 2 & iy == 1, 2,
            ifelse(ix == 1 & iy == 2, 3, 4)
          )
        )
        
        X[k, unique(region)] <- 1L
      }
      
      block_info$bx[k] <- bx
      block_info$by[k] <- by
      block_info$n_points[k] <- length(x_local)
      block_info$n_occupied_regions[k] <- sum(X[k, ])
      
      k <- k + 1
    }
  }
  
  colnames(X) <- paste0("v", 1:4)
  
  list(
    X = X,
    block_info = block_info
  )
}

obs <- make_2x2_observations(
  pp_unit = pp_unit,
  n_blocks_x = n_blocks_x,
  n_blocks_y = n_blocks_y
)

X <- obs$X
block_info <- obs$block_info

dim(X)
head(X)
table(rowSums(X))


# Graph C4 ----------------------------------------------------
#
# Local layout:
#
#   v3  v4
#   v1  v2
#
# Edges are side-neighbor relations:
# v1-v2, v2-v4, v4-v3, v3-v1.
edges <- rbind(
  c("v1", "v2"),
  c("v2", "v4"),
  c("v4", "v3"),
  c("v3", "v1")
)

colnames(edges) <- c("a", "b")

coords <- data.frame(
  vertex = paste0("v", 1:4),
  x = c(0.25, 0.75, 0.75, 0.25),
  y = c(0.25, 0.25, 0.75, 0.75)
)

g <- graph_from_edgelist(edges, directed = FALSE)

is_chordal(g)$chordal
is_tree(g)

plot(
  g,
  layout = as.matrix(coords[, c("x", "y")]),
  vertex.label = coords$vertex,
  main = "C4 graph"
)


# Support of the hard-core model ------------------------------
is_admissible <- function(x, edges) {
  all(apply(edges, 1, function(e) {
    !(x[e[1]] == 1 && x[e[2]] == 1)
  }))
}

admissible <- apply(X, 1, is_admissible, edges = edges)

table(admissible)
mean(admissible)

X_adm <- X[admissible, , drop = FALSE]

all_states <- as.matrix(expand.grid(rep(list(0:1), 4)))
colnames(all_states) <- paste0("v", 1:4)

support <- all_states[
  apply(all_states, 1, is_admissible, edges = edges),
  ,
  drop = FALSE
]

support
apply(support, 1, paste0, collapse = "")


# Observed counts on the support ------------------------------
state_id <- function(mat) {
  apply(mat, 1, paste0, collapse = "")
}

support_states <- state_id(support)
observed_states <- state_id(X_adm)

observed_counts <- setNames(integer(length(support_states)), support_states)
tab <- table(observed_states)

observed_counts[names(tab)] <- as.integer(tab)

observed_counts
sum(observed_counts)


# Fit mult_G(1, y) --------------------------------------------
#
# We use theta = log(y).
# For x in the support:
#
#   P_theta(X = x) = exp(theta^T x) / delta(theta)

loglik <- function(theta, support, counts) {
  eta <- as.vector(support %*% theta)
  
  # log-sum-exp for numerical stability
  m <- max(eta)
  log_delta <- m + log(sum(exp(eta - m)))
  
  sum(counts * eta) - sum(counts) * log_delta
}

fit_model <- function(counts, support) {
  fit <- optim(
    par = rep(0, ncol(support)),
    fn = function(theta) -loglik(theta, support, counts),
    method = "BFGS",
    control = list(maxit = 10000)
  )
  
  theta_hat <- fit$par
  names(theta_hat) <- colnames(support)
  
  eta <- as.vector(support %*% theta_hat)
  m <- max(eta)
  
  prob_hat <- exp(eta - m)
  prob_hat <- prob_hat / sum(prob_hat)
  
  expected <- sum(counts) * prob_hat
  
  T_LR <- 2 * sum(
    ifelse(counts > 0, counts * log(counts / expected), 0)
  )
  
  list(
    theta_hat = theta_hat,
    y_hat = exp(theta_hat),
    prob_hat = prob_hat,
    expected = expected,
    T_LR = T_LR,
    convergence = fit$convergence
  )
}

fit_obs <- fit_model(
  counts = as.numeric(observed_counts),
  support = support
)

fit_obs$convergence
fit_obs$T_LR
fit_obs$y_hat


# GOF table ---------------------------------------------------
decode_state <- function(s) {
  bits <- as.integer(strsplit(s, "")[[1]])
  active <- paste0("v", which(bits == 1))
  
  if (length(active) == 0) {
    return("none")
  }
  
  paste(active, collapse = " + ")
}

gof_table <- data.frame(
  state = names(observed_counts),
  observed = as.numeric(observed_counts),
  expected = fit_obs$expected,
  fitted_probability = fit_obs$prob_hat,
  profile = vapply(names(observed_counts), decode_state, character(1)),
  row.names = NULL
) |>
  mutate(
    signed_pearson_resid = (observed - expected) / sqrt(expected),
    abs_pearson_resid = abs(signed_pearson_resid),
    direction = case_when(
      observed > expected ~ "overrepresented",
      observed < expected ~ "underrepresented",
      TRUE ~ "matched"
    )
  )

gof_table |>
  arrange(desc(abs_pearson_resid))


# Parametric bootstrap ----------------------------------------
# We refit the model in every bootstrap sample.

B <- 500
set.seed(123)

T_boot <- numeric(B)

for (b in seq_len(B)) {
  boot_counts <- as.vector(
    rmultinom(
      n = 1,
      size = sum(observed_counts),
      prob = fit_obs$prob_hat
    )
  )
  
  fit_boot <- fit_model(
    counts = boot_counts,
    support = support
  )
  
  T_boot[b] <- fit_boot$T_LR
  
  if (b %% 100 == 0) {
    message("bootstrap ", b, " / ", B)
  }
}

p_value_bootstrap <- (1 + sum(T_boot >= fit_obs$T_LR)) / (B + 1)

fit_obs$T_LR
p_value_bootstrap

hist(
  T_boot,
  breaks = 40,
  main = "Bootstrap distribution of LR statistic",
  xlab = expression(T[LR])
)

abline(v = fit_obs$T_LR, lty = 2, lwd = 2)


# Final summary -----------------------------------------------

summary_result <- list(
  dataset = dataset_name,
  n_points = pp$n,
  n_blocks_x = n_blocks_x,
  n_blocks_y = n_blocks_y,
  n_observations = nrow(X),
  local_discretization = "2x2 local grid",
  graph = "C4 cycle on 2x2 side-neighbor grid",
  graph_is_decomposable = is_chordal(g)$chordal,
  p_vertices = ncol(X),
  n_edges = nrow(edges),
  n_total_patterns = nrow(X),
  n_admissible_patterns = nrow(X_adm),
  admissible_fraction = mean(admissible),
  n_support = nrow(support),
  n_observed_states = sum(observed_counts > 0),
  observed_counts = observed_counts,
  T_LR = fit_obs$T_LR,
  p_value_bootstrap = p_value_bootstrap,
  bootstrap_B = B,
  y_hat = fit_obs$y_hat,
  row_sum_table_all = table(rowSums(X)),
  row_sum_table_admissible = table(rowSums(X_adm))
)

summary_result

gof_table |>
  arrange(desc(abs_pearson_resid))