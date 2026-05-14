# Simulation Examples

This folder mixes driver scripts, one-off comparison scripts, and generated outputs.
The current standardized scenario runs use `N_sp = 2000` and `B_MAX = 3000`.

Primary driver scripts:
- `Scenario1 (Cont+Cont).R`
- `Scenario2 (Cont+Bin).R`
- `Scenario3 (Surv+Cont).R`
- `Scenario4 (Bin+Cont).R`

Current standardized summary outputs:
- `Summary.Scenario1.Nsp2000.csv`
- `Summary.Scenario2.Nsp2000.csv`
- `Summary.Scenario3.Nsp2000.csv`
- `Summary.Scenario4.Nsp2000.csv`

Supporting analysis scripts:
- `Power_Compare_Plot.R`
- `Build_Computation_Table.R`
- `Build_Decomposition_Table.R`
- `Build_Operating_Characteristic_Tables.R`
- `Build_EqualVariance_Sensitivity_Appendix.R`
- `Build_Nsp_Representative_Comparison.R`
- `Barnhart_Examples.R`

Generated artifacts:
- `convergence_plots/` contains adaptive Monte Carlo convergence PDFs
- `Comparison.Nsp8000_vs_2000.rho0.csv` contains the representative Nsp comparison summary
- `Summary.Power.Compare.Nsp2000.csv` stores theory-versus-empirical power gaps
- `Table.*.tex` stores table exports based on the standardized Nsp2000 runs

Legacy files:
- Older summary files without the `.Nsp2000` suffix are retained for archival or direct comparison with earlier runs.
