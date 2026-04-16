# Project Summary

## Project Title

Initial reproduction of radial dust-attenuation scaling relations in MaNGA galaxies

## Motivation

I prepared this project for discussion with Professor Cheng Li during the Tsinghua astronomy forum. The work is based on one of his published papers on radial dust attenuation scaling relations in MaNGA galaxies.

## Reference Paper

Li, N. & Li, C. (2024)  
**Estimating Dust Attenuation From Galactic Spectra. III. Radial variations of dust attenuation scaling relations in MaNGA galaxies**

## Current Focus

The current stage focuses on reproducing part of the radial behavior of **E(B-V)_gas**, especially its correlation with:

- log([NII]/[SII])
- log([OIII]/[OII])
- log(sigma_Ha)

across different radial bins.

## What I Have Done

### 1. Read MaNGA MAPS data
I used a MaNGA MAPS FITS file and extracted the quantities needed for emission-line analysis, including:
- [OII]
- Hβ
- [OIII]5008
- Hα
- [NII]6585
- [SII]6718+6732
- sigma_Ha
- elliptical radius (R/Re)

### 2. Apply quality cuts
I implemented:
- inverse-variance based error estimation
- S/N cuts for major lines
- simplified mask filtering
- a simplified BPT-based exclusion scheme

### 3. Compute E(B-V)_gas
Using the Balmer decrement Hα/Hβ, I estimated E(B-V)_gas for valid spaxels.

### 4. Analyze radial correlations
I divided the data into radial bins and computed the Spearman correlation coefficient between E(B-V)_gas and:
- log([OIII]/[OII])
- log([NII]/[SII])
- log(sigma_Ha)

### 5. Export comparison plots
I generated comparison plots for the radial dependence of these correlations under different selection strategies, including:
- raw selection
- BPT-cut selection
- S/N-cut selection

## Current Status

This repository represents an initial reproduction rather than a full replication of the paper.

The current emphasis is on:
- understanding the MaNGA MAPS data structure
- reproducing radial-bin correlation behavior
- testing the impact of selection cuts
- comparing my trends with the published results

## Next Steps

Planned next steps include:
- exporting numerical summary tables of rho in each radial bin
- reproducing more figures from the paper
- refining the BPT selection to better match the paper methodology
- extending the analysis to E(B-V)_star and additional properties
