# MaNGA Dust Attenuation Reproduction

This repository contains an initial reproduction of part of the analysis in the paper:

**Estimating Dust Attenuation From Galactic Spectra. III. Radial variations of dust attenuation scaling relations in MaNGA galaxies**

The current stage focuses on the radial behavior of **E(B-V)_gas** and its correlation with selected emission-line and kinematic properties in MaNGA spaxels.

## Project Goal

The goal of this project is to reproduce, at an initial stage, part of the radial dust-attenuation analysis presented in the paper, especially the relation between **E(B-V)_gas** and:

- log([NII]/[SII])
- log([OIII]/[OII])
- log(sigma_Ha)

across different radial bins in units of **R/Re**.

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
I generated figure panels for the radial dependence of these correlations under different selection strategies, including:
- raw selection
- BPT-cut selection
- S/N-cut selection

## Repository Structure

```text
manga-dust-attenuation-reproduction/
├── docs/      # project notes and stage summaries
├── src/       # analysis scripts
├── figures/   # exported radial-correlation plots
├── results/   # summary outputs
