# WinSampleSize 🏆

Sample size calculation and power evaluation for hierarchical composite endpoints using a series of derived win statistics while accounting for the correlation structure among component endpoints.

## Overview 🔍

`WinSampleSize` is an R research codebase for two-arm clinical trial design with hierarchical composite endpoints (HCEs). The workflow combines:

- copula-based generation of correlated component outcomes
- hierarchical pairwise win/loss comparison across endpoint levels
- estimation of win and loss probabilities under `H0` and `HA`
- estimation of `Xi` variance components for win-statistic inference
- theoretical sample size and power calculation for multiple win-based estimands
- empirical validation through Monte Carlo simulation

The current implementation supports the following derived win statistics:

- Net Benefit (`NB`)
- Win Ratio (`WR`)
- Win Odds (`WO`)
- DOOR-style probability scale (`DOOR`)

## What The Code Does 🧠

For a prespecified endpoint hierarchy, the code compares each treatment subject against each control subject. A pair is first evaluated on the highest-priority endpoint. If that comparison is unresolved (a tie), the pair moves to the next endpoint, and so on. This produces:

- overall win probability: `tau_w = P(Treatment wins)`
- overall loss probability: `tau_l = P(Treatment loses)`
- tie probability: `1 - tau_w - tau_l`
- endpoint-specific win/loss contributions across the hierarchy
- `Xi` variance components used in asymptotic variance formulas

The main design idea is:

1. Specify endpoint distributions under `H0` and `HA`.
2. Specify the dependence structure among endpoints.
3. Estimate `tau` and `Xi` under `H0` and `HA` from large simulated super-populations.
4. Compute the required sample size for a target power.
5. Optionally validate theoretical power with empirical power simulations.

## Supported Endpoint Types 🧩

The code currently supports five endpoint types:

- `survival`
- `ordinal`
- `binary`
- `continuous`
- `count`

Direction of benefit in the current implementation:

- `survival`: larger observed event-free time is better
- `ordinal`: larger category is better
- `binary`: larger value is better
- `continuous`: larger value is better, subject to a clinically meaningful `threshold`
- `count`: smaller value is better

Important: the order of endpoints in the list defines the clinical hierarchy.

## Repository Structure 📁

- `DynSampleGener.R`: data generation from marginal endpoint models and a copula dependence structure
- `DynWinVarEstFUNC.R`: win/loss kernel construction, `Xi` estimation, sample size, and power functions
- `Simulation Examples/SimulationEngine.R`: adaptive Monte Carlo engine used by scenario scripts
- `Simulation Examples/Scenario1 (Cont+Cont).R`: continuous + continuous example
- `Simulation Examples/Scenario2 (Cont+Bin).R`: continuous + binary example
- `Simulation Examples/Scenario3 (Surv+Cont).R`: survival + continuous example
- `Simulation Examples/Scenario4 (Bin+Cont).R`: binary + continuous example
- `Heart-FID Case1.R`: case study motivated by HEART-FID style settings
- `ShinyApp Examples/`: prototype scripts for interactive examples

Some scripts in the repository are exploratory or legacy. The main maintained workflow is centered on `DynSampleGener.R`, `DynWinVarEstFUNC.R`, and `Simulation Examples/SimulationEngine.R`.

## Requirements ⚙️

This repository is organized as sourceable R scripts rather than an installable R package.

Core packages used across the project:

```r
install.packages(c("copula", "survival", "Matrix", "pbapply", "dplyr", "mvtnorm"))
```

`parallel` is part of base R and is used for multi-core simulation.

## Endpoint Specification 🧪

Each endpoint is defined as a list. Examples:

```r
list(type = "survival", dist = "Exponential", params = list(lambda = 0.10))
list(type = "ordinal", prob = c(0.2, 0.3, 0.5))
list(type = "binary", prob = 0.35)
list(type = "continuous", params = list(mu = 30, sigma = 15), threshold = 5)
list(type = "count", params = list(lambda = 0.60))
```

Notes:

- For `Gaussian` copulas, `copula_param` must be a positive-definite correlation matrix with dimension equal to the number of endpoints.
- For `Clayton`, `Frank`, and `Gumbel`, `copula_param` is interpreted as Kendall's tau.
- `Follow_up.Time` affects survival endpoints through administrative censoring.

## Core Functions 🛠️

### Data generation

- `Generating_Sample(endpoints, copula_type, copula_param, Follow_up.Time, N.Super)`
- `CALC.Observed.Corr(data, endpoints)`
- `CALC.Observed.CorrMul(data, endpoints)`

### Win-statistic building blocks

- `Calc.Kernal.Matrix(Group.Treat, Group.Control, endpoints)`
- `Calc.Xi(Win_Kernal, Loss_Kernal)`

### Variance formulas

- `Var_NB(m, n, Xi)`
- `Var_logWR(m, n, Xi, tau_w, tau_l)`
- `Var_logWO(m, n, Xi, tau_w, tau_l)`
- `Var_DOOR(m, n, Xi)`

### Design and performance

- `Calc.SampleSize(tau_w.HA, tau_l.HA, Xi.HA, Xi.H0, tau_w.H0, tau_l.H0, alpha, beta, Sample.rho, Metric)`
- `Calc.TheoPower(tau_w.HA, tau_l.HA, Xi.HA, Xi.H0, tau_w.H0, tau_l.H0, alpha, m, Sample.rho, Metric)`
- `Calc.AttPower(RUNNING, alpha, m, n, copula_type, copula_param, endpoints.Ctrl, endpoints.Trt, Follow_up.Time, useParallel, numCores)`

`Sample.rho` is the allocation ratio `n / m`, where `m` is the treatment-arm sample size and `n` is the control-arm sample size.

## Minimal Example 💡

```r
source("DynSampleGener.R")
source("DynWinVarEstFUNC.R")

endpoints.HA <- list(
  list(type = "continuous", params = list(mu = 4, sigma = 10), threshold = 8),
  list(type = "continuous", params = list(mu = 36, sigma = 15), threshold = 6)
)

endpoints.H0 <- list(
  list(type = "continuous", params = list(mu = 3, sigma = 12), threshold = 8),
  list(type = "continuous", params = list(mu = 30, sigma = 15), threshold = 6)
)

CORR <- matrix(c(1, 0.4,
                 0.4, 1), nrow = 2)

M <- 8000
N <- 8000
Follow_up.Time <- 200

PopData_H0 <- Generating_Sample(
  endpoints = endpoints.H0,
  copula_type = "Gaussian",
  copula_param = CORR,
  N.Super = M + N,
  Follow_up.Time = Follow_up.Time
)

Pop.Treat.HA <- Generating_Sample(
  endpoints = endpoints.HA,
  copula_type = "Gaussian",
  copula_param = CORR,
  N.Super = M,
  Follow_up.Time = Follow_up.Time
)

Pop.Treat.H0 <- PopData_H0[1:M, ]
Pop.Control.H0 <- PopData_H0[(M + 1):(M + N), ]

Kernel.H0 <- Calc.Kernal.Matrix(
  Group.Treat = Pop.Treat.H0,
  Group.Control = Pop.Control.H0,
  endpoints = endpoints.H0
)

Kernel.HA <- Calc.Kernal.Matrix(
  Group.Treat = Pop.Treat.HA,
  Group.Control = Pop.Control.H0,
  endpoints = endpoints.HA
)

Xi.H0 <- Calc.Xi(Kernel.H0$Win_Kernal, Kernel.H0$Loss_Kernal)
Xi.HA <- Calc.Xi(Kernel.HA$Win_Kernal, Kernel.HA$Loss_Kernal)

ss_wr <- Calc.SampleSize(
  tau_w.HA = Kernel.HA$tau_w,
  tau_l.HA = Kernel.HA$tau_l,
  Xi.HA = Xi.HA,
  Xi.H0 = Xi.H0,
  tau_w.H0 = Kernel.H0$tau_w,
  tau_l.H0 = Kernel.H0$tau_l,
  alpha = 0.05,
  beta = 0.15,
  Sample.rho = 1,
  Metric = "WR"
)

power_wr <- Calc.TheoPower(
  tau_w.HA = Kernel.HA$tau_w,
  tau_l.HA = Kernel.HA$tau_l,
  Xi.HA = Xi.HA,
  Xi.H0 = Xi.H0,
  tau_w.H0 = Kernel.H0$tau_w,
  tau_l.H0 = Kernel.H0$tau_l,
  alpha = 0.05,
  m = ss_wr$m.sample,
  Sample.rho = 1,
  Metric = "WR"
)

ss_wr
power_wr
```

In practice, design quantities should usually be averaged over many Monte Carlo replicates instead of relying on a single large draw. The scenario scripts in `Simulation Examples/` automate that process.

## Running The Included Simulations ▶️

Example scripts can be run directly with `Rscript`:

```bash
Rscript "Simulation Examples/Scenario1 (Cont+Cont).R"
Rscript "Simulation Examples/Scenario2 (Cont+Bin).R"
Rscript "Simulation Examples/Scenario3 (Surv+Cont).R"
Rscript "Simulation Examples/Scenario4 (Bin+Cont).R"
Rscript "Heart-FID Case1.R"
```

The scenario scripts use the adaptive engine in `Simulation Examples/SimulationEngine.R` to:

- estimate `tau` and `Xi` until Monte Carlo precision targets are met
- determine a fixed baseline sample size, typically from `WR`
- evaluate theoretical power, empirical power, and type I error
- export summary tables and convergence plots

## Derived Measures 📐

The code works with the following quantities:

- `NB = tau_w - tau_l`
- `WR = tau_w / tau_l`
- `WO = (1 + NB) / (1 - NB)`
- `DOOR = 0.5 * (1 + NB)`

These quantities are all evaluated through asymptotic variance formulas built from the estimated `Xi` components.

## Practical Notes 📝

- For multi-endpoint HCEs, the dependence structure can materially affect both the win probabilities and the required sample size.
- Continuous endpoints use threshold-based comparisons, so the design is sensitive to the clinically meaningful margin you specify.
- The codebase is research-oriented and script-driven; it is best used as a reproducible framework rather than as a polished R package.
- Prototype files in `ShinyApp Examples/` may lag behind the current core scripts.

## Citation And Context 📚

This repository was built to support methodological work on sample size determination for hierarchical composite endpoints using correlated win statistics, including case-study style analyses such as the HEART-FID setting.
