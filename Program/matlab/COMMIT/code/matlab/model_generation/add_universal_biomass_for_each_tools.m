% Add the universal prokaryotic biomass reaction from CarveMe models to all
% models, dependent on Gram stain information (Machado et al., 2018, Nucleic Acid Res.;
% Xavier et al., 2017, Metab. Eng.)
opts = delimitedTextImportOptions("NumVariables", 1);
opts.DataLines = [2, Inf];
biomass_rxn = readtable("~/IMIC/Program/matlab/COMMIT/data/gap-filling/Biomass_reaction.csv", opts);
biomass_rxn = table2cell(biomass_rxn);

biomass_name = 'universal_Biomass_reaction';

methods = {'carveme','gapseq','kbase'};
modelTopDir =  '~/IMIC/models';

for k=1:numel(methods)
    disp(methods{k})
    workspace = fullfile(modelTopDir, methods{k}, strcat(methods{k}, '_models_metFormulas.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};

        model = addReaction(model, 'BIOMASS_Reaction', ...
            'reactionName', biomass_name',...
            'reactionFormula',  biomass_rxn{:},...
            'reversible', false,...
            'lowerBound', 0,...
            'upperBound', 1000, ...
            'objectiveCoef', 1,...
            'printLevel', 0);

        if strcmp(methods{k},'carveme')
            model = addMetabolite(model, 'BIOMASS[c]');
        end

        % add BIOMASS as a product to the biomass reaction
        model = changeRxnMets(model, 'BIOMASS[c]', 'BIOMASS[c]',...
            'BIOMASS_Reaction', 1);
        model = removeRxns(model, 'EX_BIOMASS_c');
        % add an exchange reaction for BIOMASS
        [~,model]=evalc('addSinkReactions(model, ''BIOMASS[c]'', 0, 1000);');
        
        models{j} = model;
    end
    
    workspace = fullfile(modelTopDir, strcat(methods{k}, '_draft_models_biomass.mat'));
    save(workspace, 'models')
    clear models
end
