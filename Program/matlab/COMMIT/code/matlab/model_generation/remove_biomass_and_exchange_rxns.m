% remove exchange reactions and biomass reactions
options; clear

modelTopDir = '~/IMIC/models';
%modelTopDir = '~/IMIC/study_case/models';
methods = {'carveme','gapseq','kbase'};

for i = 1:numel(methods)
    disp(methods{i})
    workspace = fullfile(modelTopDir, methods{i}, strcat(methods{i},'_models.mat'));
    load(workspace)

    for j=1:numel(models)
        model = models{j};
        ec = model.EC;
        idx_bio = strcmp(model.rxns, 'BIOMASS_Reaction');
        
        % find exchange reactions
        exchange_rxns = logical(findExchangeReactions(model))';
        
        to_remove = model.rxns(logical(exchange_rxns+idx_bio));
        model = removeRxns(model, to_remove, 'metFlag', false);
        
        % deal with EC numbers separately since there can be errors if the
        % number of metabolites and reactions is equal
        model.EC = ec(~logical(exchange_rxns + idx_bio));
        
        models{j} = model;
    end
     disp('done')

     save_workspace = fullfile(modelTopDir, methods{i}, strcat(methods{i},'_models_no_medium_no_biomass.mat'));
     save(save_workspace,'models')
     clear models
end
