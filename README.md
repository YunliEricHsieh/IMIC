# **Integration of metatranscriptomics data improves the predictive capacity of microbial community metabolic models**

## Description

This repository contains code to produce the results presented in the manuscript: [Integration of metatranscriptomics data improves the predictive capacity of microbial community metabolic models](https://doi.org/10.1093/ismejo/wraf109).


## Dependencies
- [MATLAB](https://www.mathworks.com/products/matlab.html) (tested on v2023b)
- [COBRA](https://github.com/opencobra/cobratoolbox/tree/master) (tested onv3.4)
- [COMMIT](https://github.com/pwendering/COMMIT) (tested on v1.1.2)
- [R](https://www.r-project.org/) (tested on 4.4.3)
- [Gurobi solver](https://support.gurobi.com/hc/en-us/articles/4534161999889-How-do-I-install-Gurobi-Optimizer) (tested on version 11.0.0)

## Workflow for the publication results
- The script `IMIC/Program/matlab/model_generation/build_community_model.m` contains the workflow for building the model
- The script `Program/matlab/calculate_community_growth.m` performs the prediction of community growth rates for the IMIC and CoCo-GEM approaches
