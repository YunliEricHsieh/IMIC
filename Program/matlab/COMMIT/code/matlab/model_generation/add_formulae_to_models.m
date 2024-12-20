% add formulae to models

formulaTab = readtable('~/IMIC/Program/matalb/COMMIT/data/tables/MNXref/MNXref_MET_FORMULAE.csv',...
    'ReadVariableNames', true);
methods = {'carveme','gapseq','kbase'};
modelTopDir = '~/IMIC/study_case/models';

for j = 1:numel(methods)
    disp(methods{j})
    workspace = fullfile(modelTopDir, methods{j}, strcat(methods{j}, '_models_no_medium_no_biomass.mat'));
    load(workspace)
        
    for k=1:numel(models)
        models{k} = addMetFormulae(models{k}, formulaTab);
    end

    save_workspace = fullfile(modelTopDir, methods{j},  strcat(methods{j},'_models_metFormulas.mat'));
    save(save_workspace, 'models')
    clear models
end
