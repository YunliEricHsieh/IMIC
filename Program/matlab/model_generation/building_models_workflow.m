% load options
options

% KBase, CarveMe, and gapseq draft models were created outside of
% MATLAB and are available upon request

%% Convert draft models to a common format that can be used for evaluation
% % Convert and translate the draft models to MNXref namespace
CarveMe_model_conversion

gapseq_model_conversion

KBase_model_conversion

%% Remove all exchange and biomass reactions from the models
remove_biomass_and_exchange_rxns

%% Add metabolite formulae
add_formulae_to_models

%% Merge models from different approaches
merge_metabolic_models

%% add a universal prokaryotic biomass to all models
add_universal_biomass_for_each_tools

% add a universal prokaryotic biomass to consensus model
add_universal_biomass

%% Iterative gap filling (COMMIT)
run_iterative_gap_filling

%% Building community models in one big model
build_community_model
