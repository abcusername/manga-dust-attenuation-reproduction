# Stage 1 Reproduction Notes

## Completed Work

- inspected the paper structure and identified the radial-scaling part as the first target for reproduction
- used a MaNGA MAPS file for one-galaxy testing
- extracted R/Re, line fluxes, and sigma_Ha from the FITS data
- constructed E(B-V)_gas from Halpha/Hbeta
- applied simplified quality-control cuts
- computed radial-bin Spearman correlations
- exported comparison plots for N2S2, O3O2, and sigma_Ha

## Current Limitation

This is still an initial reproduction:
- the analysis currently uses a simplified BPT-cut scheme
- only a subset of the paper’s variables and figures has been reproduced
- numerical summary tables still need to be exported

## Next Step

- export numerical rho tables
- reproduce more of the paper’s figure logic
- improve consistency with the original selection criteria
- compare trends more explicitly with the published paper
