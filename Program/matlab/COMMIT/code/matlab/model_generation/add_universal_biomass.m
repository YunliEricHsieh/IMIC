% Add the universal prokaryotic biomass reaction from CarveMe models to all
% models, dependent on Gram stain information (Machado et al., 2018, Nucleic Acid Res.;
% Xavier et al., 2017, Metab. Eng.)
opts = delimitedTextImportOptions("NumVariables", 1);
opts.DataLines = [2, Inf];
biomass_rxn = readtable("~/IMIC/Program/matlab/COMMIT/data/gap-filling/Biomass_reaction.csv", opts);
biomass_rxn = table2cell(biomass_rxn);

biomass_name = 'universal_Biomass_reaction';

methods = {'consensus'};
modelTopDir =  '~/IMIC/models';
%modelTopDir =  '~/IMIC/study_case/models';

%% Merged models
disp('Consensus')
for i=1:numel(methods)
    workspace = fullfile(modelTopDir, strcat(methods{i}, '_draft_models.mat'));
    load(workspace)
    for j=1:numel(merged_models)
        model = merged_models{j};
        
        model = addReaction(model, 'BIOMASS_Reaction',...
            'reactionName', biomass_name',...
            'reactionFormula',  biomass_rxn{:},...
            'reversible', false,...
            'lowerBound', 0,...
            'upperBound', 1000, ...
            'objectiveCoef', 1,...
            'printLevel', 0);

        % add BIOMASS as a product to the biomass reaction
        model = changeRxnMets(model, 'BIOMASS[c]', 'BIOMASS[c]',...
            'BIOMASS_Reaction', 1);
        model = removeRxns(model, 'EX_BIOMASS_c');
        % add an exchange reaction for BIOMASS
        [~,model]=evalc('addSinkReactions(model, ''BIOMASS[c]'', 0, 1000);');
        
        merged_models{j} = model;
    end
    workspace = fullfile(modelTopDir, strcat(methods{i}, '_draft_models_biomass.mat'));
    save(workspace, 'merged_models')
    clear merged_models
end
