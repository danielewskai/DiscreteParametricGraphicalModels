# Graphical Count Models: Real-Data Experiments

This repository contains code for real-data experiments illustrating goodness-of-fit procedures for graphical count models, with a focus on the case \(r=1\) of the graphical multinomial model \(\mathrm{mult}_G(1,y)\).

The repository contains two types of experiments:

1. a Rydberg-atom experiment, where each observation is a binary measurement shot;
2. a spatial point-pattern experiment based on retinal amacrine cells, where binary observations are obtained by discretizing local spatial blocks.

Both experiments use the same statistical pipeline:

1. construct binary observations \(X^{(1)},\ldots,X^{(n)}\in\{0,1\}^V\);
2. fix a graph \(G=(V,E)\);
3. retain only observations belonging to the graph-constrained support;
4. fit \(\mathrm{mult}_G(1,y)\) by maximum likelihood;
5. compute a likelihood-ratio goodness-of-fit statistic;
6. obtain a p-value by parametric bootstrap, refitting the model in each bootstrap sample.

The activity parameters \(y\) are estimated from the data. Therefore, the bootstrap refits the model in every bootstrap sample.

---

## Model

Let \(G=(V,E)\) be a graph. In the case \(r=1\), the support of the model consists of binary vectors satisfying the graph constraints:

\[
\mathcal N_{G,1}=\{x\in\{0,1\}^{V}:x_u x_v=0\text{ for every }uv\in E\}.
\]

Equivalently, no two adjacent vertices of \(G\) can be simultaneously active.

The fitted model is

\[
\mathbb P_\underline{y}(\underline{X}=\underline{x})=\frac{\prod_{v\in V} y_v^{x_v}}{\delta_G(\underline{y})},\qquad x\in\mathcal N_{G,1},
\]

where

\[
\delta_G(\underline{y})=\sum_{z\in\mathcal N_{G,1}}\prod_{v\in V} y_v^{z_v}.
\]

Equivalently, with \(\theta_v=\log y_v\),

\[
\mathbb P_\underline{\theta}(\underline{X}=\underline{x})=\frac{\exp(\underline{\theta}^\top \underline{x})}{\sum_{\underline{z}\in\mathcal N_{G,1}}\exp(\underline{\theta}^\top \underline{z})}.
\]

In the experiments in this repository, the support \(\mathcal N_{G,1}\) is small enough to be enumerated exactly. Therefore, the likelihood and the normalizing constant are computed directly by summing over all admissible configurations.

---

## Goodness-of-fit procedure

For each experiment, binary observations are first filtered to the graph-constrained support:

\[
X^{(i)}\in\mathcal N_{G,1}.
\]

Let \(O_underline{x}\) be the observed count of configuration \(x\) among the admissible observations. After fitting the model, the expected count is

\[
E_\underline{x}=n_{\mathrm{adm}}\widehat{\mathbb P}(\underline{X}=\underline{x}),
\]

where \(n_{\mathrm{adm}}\) is the number of admissible observations.

The likelihood-ratio statistic is

\[
T_{\mathrm{LR}}=2\sum_{\underline{x}\in\mathcal N_{G,1}}O_\underline{x}\log\frac{O_\underline{x}}{E_\underline{x}},
\]

with the convention that terms with \(O_\underline{x}=0\) are equal to zero.

The p-value is computed by parametric bootstrap:

1. simulate bootstrap counts from the fitted model;
2. refit the model in the bootstrap sample;
3. recompute \(T_{\mathrm{LR}}\);
4. compare the observed statistic with the bootstrap distribution.

The bootstrap p-value is computed as

\[
p_val=\frac{1+\#\{T_b\ge T_{\mathrm{obs}}\}}{B+1},
\]

where \(B\) is the number of bootstrap replications.
