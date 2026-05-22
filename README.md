# SINDyMJO-RMM-ENSO

This repository contains the R scripts used in an undergraduate thesis on data-driven modelling of the Madden-Julian Oscillation (MJO) under different ENSO backgrounds. The workflow applies Sparse Identification of Nonlinear Dynamics (SINDy) to the Real-time Multivariate MJO (RMM) index and compares the inferred low-dimensional dynamical structures between boreal winter and summer.

## Repository Contents

- `SindyMJO_RMM_ONI.R`: main analysis workflow. It reads RMM and ONI data, constructs winter and summer samples, applies smoothing and velocity estimation, filters active MJO days, runs SINDy, selects models using corrected AIC, classifies model structures, and extracts periods from linear models.
- `MJO_polar.R`: helper function that converts `RMM1` and `RMM2` into MJO amplitude and phase.

Raw data files and generated outputs are not included in this repository. They should be downloaded or regenerated locally.

## Data Sources

The workflow uses public climate index datasets:

- RMM index: Australian Bureau of Meteorology, Real-time Multivariate MJO Index  
  https://www.bom.gov.au/climate/mjo/graphics/rmm.74toRealtime.txt
- ONI index: NOAA Climate Prediction Center, Oceanic Niño Index  
  https://cpc.ncep.noaa.gov/data/indices/oni.ascii.txt

## Method Summary

The analysis treats `RMM1` and `RMM2` as a two-dimensional phase-space representation of the MJO. The main steps are:

1. Read RMM and ONI data.
2. Match daily RMM records with ENSO background states.
3. Construct boreal winter and summer samples.
4. Apply a 9-point centered moving average to `RMM1` and `RMM2`.
5. Estimate velocities using finite-difference schemes.
6. Select active MJO days using an amplitude threshold.
7. Perform repeated random sampling and SINDy sparse regression.
8. Select candidate models using corrected AIC.
9. Classify model structures according to nonzero coefficient masks.
10. Estimate oscillation periods from purely linear models.

## R Environment

Install missing packages before running the workflow. The scripts were developed in R and use the following packages:

- `ggplot2`
- `reshape2`
- `ggpubr`
- `patchwork`
- `R1magic`
- `latex2exp`
- `viridis`

## Usage

From the repository directory:

```bash
Rscript SindyMJO_RMM_ONI.R [rmm.txt] [oni.txt] [out_dir] [y0] [y1] [Nreal] [Nsamp] [amp] [dpi] [years_mode] [oni_thr] [seed]
```

Main arguments:

- `rmm.txt`: path to the RMM input file.
- `oni.txt`: path to the ONI input file.
- `out_dir`: output directory.
- `y0`, `y1`: analysis years.
- `Nreal`: number of random realizations.
- `Nsamp`: sample size for each realization.
- `amp`: MJO amplitude threshold.
- `years_mode`: `all`, `neutros`, `ninos`, or `ninas`.
- `oni_thr`: ONI threshold for ENSO classification.
- `seed`: random seed.

## Outputs

The script writes figures and CSV tables to the specified output directory. Typical outputs include coefficient statistics, model-structure masks, model-class frequencies, representative model coefficients, and linear-model period distributions.


## Notes

The workflow is adapted from the SINDy-MJO code of Diaz, Barreiro, and Rubido and modified for RMM input and ONI-based ENSO classification. Users should cite the original data providers and relevant methodological references when reusing the workflow.
