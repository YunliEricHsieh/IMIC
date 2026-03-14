# IMIC: Integration of Metatranscriptomes Into Community GEMs

[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fismejo%2Fwraf109-blue.svg)](https://doi.org/10.1093/ismejo/wraf109)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b-orange)](https://www.mathworks.com/products/matlab.html)
[![R](https://img.shields.io/badge/R-4.4+-blue)](https://www.r-project.org/)

**IMIC** (Integration of Metatranscriptomes Into Community GEMs) is a method that integrates **metatranscriptomic data** into microbial community genome-scale metabolic models (GEMs). It dramatically improves predictions of individual species growth rates and metabolic interactions compared to abundance-only approaches (MICOM).

This repository contains **all code, models, data, and scripts** needed to reproduce the results from the paper:

> **Integration of metatranscriptomics data improves the predictive capacity of microbial community metabolic models**  
> Yunli Eric Hsieh, Kshitij Tandon, Heroen Verbruggen, Zoran Nikoloski  
> *The ISME Journal*, Volume 19, Issue 1, wraf109 (2025)  
> [DOI: 10.1093/ismejo/wraf109](https://doi.org/10.1093/ismejo/wraf109)

## What is IMIC?

Traditional community GEMs rely only on species abundance. **IMIC** uses metatranscriptomic profiles (TPM values) to impose **condition-specific flux constraints** via GPR rules with a relaxation parameter β and a single balancing factor λ (automatically determined via sensitivity analysis).  
The result: growth rates that correlate strongly with both relative and absolute abundances, plus superior prediction of metabolite concentration changes and inter-species interactions.

## Key Features

- Automated construction of compartmentalized community GEMs (supports CarveMe, gapseq, KBase, etc.)
- IMIC flux constraint integration (one tunable parameter λ)
- Direct comparison with CoCo-GEM and MICOM
- Sensitivity analysis, flux-sum analysis, and flux variability analysis
- Full reproduction of the paper's ganjang fermentation community and synthetic consortia case studies
- Ready-to-run R scripts that generate all paper figures

## Dependencies

- **MATLAB** (tested on R2023b)
- **COBRA Toolbox** v3.4
- **COMMIT** v1.1.2
- **Gurobi Optimizer** v11.0 (recommended; academic license free)
- **R** (tested on 4.4.3) – for figure generation

### Installation (5 minutes)

1. Clone the repository:
   ```bash
   git clone https://github.com/YunliEricHsieh/IMIC.git
   cd IMIC

2. Install COBRA Toolbox and COMMIT (follow official guides).
3. Install and link Gurobi to MATLAB.
4. Add the following folders to your MATLAB path:
   ```matlab
   addpath(genpath('Program/matlab'));
   ```

## Reproducing the Paper Results
All scripts are in Program/matlab/. Run them in MATLAB:
```matlab
% 1. Build the community model
build_community_model

% 2. Predict community growth rates (IMIC vs CoCo-GEM)
calculate_community_growth

% 3. Run sensitivity analysis for the parameter λ
sensitivity_analysis
```
- study_case.m – Main case-study analysis
- R scripts in Program/R/ (Fig2.R, Fig3.R, …, FigS*.R) generate all paper figures once the MATLAB tables are produced.

## Repository Structure
```
IMIC/
├── Program/
│   ├── matlab/               # Core code (model building, IMIC implementation, comparisons)
│   │   ├── model_generation/
│   │   ├── IMIC/             # Core IMIC functions
│   │   ├── CoCo-GEM/, MICOM/, COMMIT/
│   │   └── *.m               # Main entry-point scripts
│   └── R/                    # Figure generation (ggplot2-based)
├── models/                   # Draft and consensus GEMs (.mat)
├── study_case/               # Input data for case studies
├── table/                    # Generated result tables
├── media/                    # Culture media definitions (.mat)
├── LICENSE
└── README.md
```

## Citation
Please cite the paper if you use IMIC:
```
@article{hsieh2025imic,
  title   = {Integration of metatranscriptomics data improves the predictive capacity of microbial community metabolic models},
  author  = {Hsieh, Yunli Eric and Tandon, Kshitij and Verbruggen, Heroen and Nikoloski, Zoran},
  journal = {The ISME Journal},
  volume  = {19},
  number  = {1},
  pages   = {wraf109},
  year    = {2025},
  doi     = {10.1093/ismejo/wraf109}
}
```
