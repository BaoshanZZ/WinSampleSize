# HEART-FID Example Organization

This folder is the maintained HEART-FID workflow. The example is split into two branches:

## FORSS branch

Final scripts:

- `Calibrate_HeartFID_Copula.R`: calibrates the latent Gaussian copula so simulated observed associations match the HEART-FID Harrell-C/Kendall targets.
- `Run_HeartFID_Independent_vs_Calibrated.R`: runs the maintained FORSS analysis for the independent, direct-input, and calibrated correlation scenarios.
- `Build_HeartFID_Independent_vs_Calibrated_Tables.R`: builds the final manuscript tables and result text from the maintained summary CSV.

Final outputs:

- `HeartFID_output/Summary_HEARTFID_Independent_vs_Calibrated.csv`
- `HeartFID_output/Table.HeartFID.Power.Main.tex`
- `HeartFID_output/Table.HeartFID.Power.ThreeScenarios.tex`
- `HeartFID_output/Table.HeartFID.Decomposition.ThreeScenarios.tex`
- `HeartFID_output/Table.HeartFID.Computation.tex`
- `HeartFID_output/Results.HeartFID.Main.tex`

Support/cache outputs:

- `HeartFID_output/calibration/`
- `HeartFID_output/intermediate/`
- `HeartFID_output/convergence_plots/`

Current observed association targets:

- `D1--D2`: `2 * Harrell C(D1, HF hospitalization count) - 1 = -0.2211`
- `D1--D3`: `2 * Harrell C(D1, 6MWD change) - 1 = 0.5173`
- `D2--D3`: `Kendall tau-b(HF hospitalization count, 6MWD change) = -0.1032`

The direct-input scenario uses these same raw endpoint-scale targets as latent
Gaussian correlations.

Current adaptive Monte Carlo stopping settings in the maintained run:

- `EPSILON_tau = 1e-3`
- `EPSILON_xi = 1e-4`

## Barnhart branch

Reference scripts copied from the Barnhart implementation are kept in:

- `Barnhart_reference/`

These are reference materials for comparison, not the maintained FORSS run scripts.

## Deleted as redundant/test files

The old adaptive/test HEART-FID files were removed because they were superseded by the maintained workflow above:

- root-level `Heart-FID Case1.R`
- root-level `Indep_T0.R`
- `Sim2.R`
- output-directory copies `Sim1.R`, `Sim_3.R`, and `Sim4.R`
- old adaptive output files: `Summary_HEARTFID_Adaptive.csv`, `precalc_baseline.rds`, `result_*.rds`, and `legacy/`
- duplicate table names: `Table.HeartFID.*.TwoScenarios.tex`

Run order:

```bash
Rscript RealStudy_HeartFID/Run_HeartFID_Independent_vs_Calibrated.R
Rscript RealStudy_HeartFID/Build_HeartFID_Independent_vs_Calibrated_Tables.R
```
