# Graphical Count Models: Real-Data Experiments

This repository contains code for real-data experiments illustrating goodness-of-fit procedures for graphical count models, with a focus on the case $r=1$ of the graphical multinomial model $\mathrm{mult}_G(1,\underline{y})$.

The repository contains two real-data examples:

1. a Rydberg-atom experiment, where each observation is a binary measurement shot;
2. a spatial point-pattern experiment based on retinal amacrine cells, where binary observations are obtained by discretizing local spatial blocks.

Both examples use the same support-restricted goodness-of-fit methodology.

---

## Data used in this repository

### 1. Rydberg-atom data

The first experiment uses binary measurement data from a Rydberg-atom experiment.

In this setting, each graph vertex represents an atom. Each measurement shot records whether each atom is observed in the Rydberg-excited state or not. Therefore, every raw observation is already a binary vector

```math
\underline{X}^{(i)} \in \{0,1\}^{V}.
```

The graph $G=(V,E)$ represents pairwise incompatibility constraints between atoms. If two vertices are connected by an edge, then the two corresponding atoms should not be simultaneously excited in an ideal graph-constrained measurement.

The code analyzes selected four-vertex graph configurations:

| Experiment | Graph | Description |
|---|---:|---|
| `fig2d` | $P_4$ | path on four vertices |
| `fig2j` | $S_4$ | star on four vertices |
| `fig2e` | $C_4$ | cycle on four vertices |
| `fig2k` | PAN4 | paw graph on four vertices |

For each graph, the analysis first retains only observations that satisfy the graph constraints and then fits $\mathrm{mult}_G(1,\underline{y})$ on this admissible part of the sample.

The corresponding code is contained in the Jupyter notebook:

    Rydberg_code.ipynb

or in the corresponding notebook file in this repository.

---

### 2. Amacrine-cell point-pattern data

The second experiment uses the `amacrine` dataset from the R package `spatstat.data`.

This dataset contains spatial locations of displaced amacrine cells in the retina of a rabbit. The point pattern is marked by cell type, with two types usually denoted as `"on"` and `"off"`. In the main analysis, only the `"on"` cells are used.

Unlike the Rydberg data, the raw amacrine observations are not binary vectors. They are spatial point locations. Therefore, the code first transforms the point pattern into binary local occupancy vectors.

The observation window is normalized to $[0,1]^2$ and divided into $12 \times 12$ local blocks. Each local block is then divided into a $2 \times 2$ grid.

Each block gives one binary vector

$$
\underline{X}^{(i)}=(X^{(i)}_1,X^{(i)}_2,X^{(i)}_3,X^{(i)}_4)\in\{0,1\}^4,
$$

where $X^{(i)}_j=1$ if at least one `"on"` amacrine cell falls in subregion $v_j$ of block $i$, and $X^{(i)}_j=0$ otherwise.

The graph used in the amacrine experiment is the cycle $C_4$ on the local $2\times2$ grid:

$$
v_1-v_2,\qquad
v_2-v_4,\qquad
v_4-v_3,\qquad
v_3-v_1.
$$

Thus, admissible configurations are those in which no two side-adjacent subregions are simultaneously occupied.

The corresponding code is contained in the R script:

    amacrine_on_C4.R

or in the corresponding R file in this repository.

---

## Statistical pipeline

Both experiments use the same statistical pipeline:

1. construct binary observations $\underline{X}^{(1)},\ldots,\underline{X}^{(n)}\in\{0,1\}^V$;
2. fix a graph $G=(V,E)$;
3. retain only observations belonging to the graph-constrained support;
4. fit $\mathrm{mult}_G(1,\underline{y})$ by maximum likelihood;
5. compute a likelihood-ratio goodness-of-fit statistic;
6. obtain a p-value by parametric bootstrap, refitting the model in each bootstrap sample.

The activity parameters $\underline{y}$ are estimated from the data. Therefore, the bootstrap refits the model in every bootstrap sample.

---

## Model

Let $G=(V,E)$ be a graph. In the case $r=1$, the support of the model consists of binary vectors satisfying the graph constraints:

```math
\mathcal{N}_{G,1}
=
\left\{
\underline{x}\in\{0,1\}^{V}
:
x_u x_v=0
\ \forall\, uv\in E
\right\}.
```

Equivalently, no two adjacent vertices of $G$ can be simultaneously active.

The fitted model is

$$
\mathbb{P}_{\underline{y}}(\underline{X}=\underline{x})=\frac{\prod_{v\in V} y_v^{x_v}}{\delta_G(\underline{y})},\qquad\underline{x}\in\mathcal{N}_{G,1},
$$

where

$$
\delta_G(\underline{y})=\sum_{\underline{z}\in\mathcal{N}_{G,1}}\prod_{v\in V} y_v^{z_v}.
$$

Equivalently, with $\theta_v=\log y_v$,

$$
\mathbb{P}_{\underline{\theta}}(\underline{X}=\underline{x})=\frac{\exp(\underline{\theta}^{\top}\underline{x})}{\sum_{\underline{z}\in\mathcal{N}_{G,1}}\exp(\underline{\theta}^{\top}\underline{z})}.
$$

In the experiments in this repository, the support $\mathcal{N}_{G,1}$ is small enough to be enumerated exactly. Therefore, the likelihood and the normalizing constant are computed directly by summing over all admissible configurations.

---

## Goodness-of-fit procedure

For each experiment, binary observations are first filtered to the graph-constrained support:

$$
\underline{X}^{(i)}\in\mathcal{N}_{G,1}.
$$

Let $O_{\underline{x}}$ be the observed count of configuration $\underline{x}$ among the admissible observations. After fitting the model, the expected count is

$$
E_{\underline{x}}=n_{\mathrm{adm}}\widehat{\mathbb{P}}(\underline{X}=\underline{x}),
$$

where $n_{\mathrm{adm}}$ is the number of admissible observations.

The likelihood-ratio statistic is

$$
T_{\mathrm{LR}}=2\sum_{\underline{x}\in\mathcal{N}_{G,1}}O_{\underline{x}}\log\frac{O_{\underline{x}}}{E_{\underline{x}}},
$$

with the convention that terms with $O_{\underline{x}}=0$ are equal to zero.

The p-value is computed by parametric bootstrap:

1. simulate bootstrap counts from the fitted model;
2. refit the model in the bootstrap sample;
3. recompute $T_{\mathrm{LR}}$;
4. compare the observed statistic with the bootstrap distribution.

The bootstrap p-value is computed as

```math
p_{\mathrm{val}}
=
\frac{
1+
\left|
\left\{
b : T_b \ge T_{\mathrm{obs}}
\right\}
\right|
}{
B+1
}.
```

where $B$ is the number of bootstrap replications.

---

## How to run the experiments

### Rydberg experiment

Open the notebook

    Rydberg_code.ipynb

and run all cells.

The notebook loads the data, defines the graph, filters observations to the graph-constrained support, fits $\mathrm{mult}_G(1,\underline{y})$, computes $T_{\mathrm{LR}}$, and evaluates the p-value by parametric bootstrap with refitting.

### Amacrine experiment

Open R or RStudio and run:

    source("amacrine_on_C4.R")

The script returns the main object

    summary_result

and the goodness-of-fit table

    gof_table

The output includes:

- the number of total observations;
- the number and fraction of admissible observations;
- the graph-constrained support;
- observed counts on the support;
- fitted activity parameters $\widehat{\underline{y}}$;
- the likelihood-ratio statistic $T_{\mathrm{LR}}$;
- the bootstrap p-value;
- expected counts and Pearson residuals for admissible configurations.
