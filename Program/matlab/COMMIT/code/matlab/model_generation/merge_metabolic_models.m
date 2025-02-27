% merge draft metabolic models from different approaches
options
clearvars -except dbFile ncpu

% set up parallel pool
c = parcluster;
c.NumWorkers = ncpu;
delete(gcp('nocreate'))
P = parpool(c);

methods = {'carveme', 'gapseq', 'kbase'};
modelTopDir = '~/IMIC/models';
%modelTopDir = '~/IMIC/study_case/models';

disp('-------------------------------------------------------------------')
disp('START')
disp('-------------------------------------------------------------------')

disp('loading the universal database')
load(dbFile)

disp('-------------------------------------------------------------------')
disp('Loading and collecting models from the different approaches...')
for j=1:numel(methods)
    workspace = fullfile(modelTopDir, methods{j},...
    strcat(methods{j}, '_models_metFormulas.mat'))

    load(workspace);
    eval([methods{j}, '= models;']);
end
disp('------------------------------')
merged_models = {};

for j=1:numel(models)
    id = models{j}.id;
    fprintf('Model #%d (%s)\n', j, id)
    fprintf('Collecting %d models...\n', numel(methods))
        
    % create models variable for every OTU in each habitat
    models_to_merge = {};
    for k=1:numel(methods)
        eval(['models_to_merge = vertcat(models_to_merge,',...
            methods{k}, '{j});'])
    end
    disp('------------------------------')
        
    % run merging function
    merged_models{j} = mergeModels(models_to_merge, dbModel_MNXref_balanced);
    wo_del = sum(cellfun(@(x)numel(x.rxns), models_to_merge));
    w_del = numel(merged_models{j}.rxns);
    fprintf('Number of deleted (merged) reactions:\t%d\n', wo_del-w_del)
    disp('------------------------------')
end
    
workspace = fullfile(modelTopDir,'consensus_draft_models');
disp('saving workspace')
save(workspace, 'merged_models')
clear models
disp('-------------------------------------------------------------------')
