modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';
methods = {'consensus','carveme','gapseq','kbase'};
timepoint = {'20d', '40d', '60d', '90d', '180d'};

ncpu = 14;
delete(gcp('nocreate'));
parpool(ncpu);

%% testing with the opt lambda from each methods
lambda = [12, 20, 15, 25];

for i = 1:numel(methods)
    model_workspace = fullfile(modelDir, [methods{i},'_com.mat']);
    load(model_workspace)
    
    results = {};
    
    for j = 1:numel(timepoint)
        disp('----------------------------------------------------------------------')
        fprintf('\n################# %s\n\n', timepoint{j})
        transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{j}, '.csv']);
        transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);

        disp('Calculate Community Growth Rate')
        com_solution = IMIC(com_model, transcriptTable, lambda(i));

        maximum_value = -com_solution.objval;

        results_table = var_analysis(com_model, transcriptTable,lambda(i), maximum_value);

        results = [results; results_table];
    end

    results = cell2table(vertcat(results{:}));
    results.Properties.VariableNames = {'ModelID','Max_growth_rate','Min_growth_rate'};
    writetable(results, fullfile(tablesDir, 'variability_analysis',[methods{i},'_variability_analysis_with_opt_lambda.csv']));
end

%% testing with the same lambda
lambda = 12;

for i = 1:numel(methods)
    model_workspace = fullfile(modelDir, [methods{i},'_com.mat']);
    load(model_workspace)
    
    results = {};
    
    for j = 1:numel(timepoint)
        disp('----------------------------------------------------------------------')
        fprintf('\n################# %s\n\n', timepoint{j})
        transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{j}, '.csv']);
        transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);

        disp('Calculate Community Growth Rate')
        com_solution = IMIC(com_model, transcriptTable, lambda);

        maximum_value = -com_solution.objval;

        results_table = var_analysis(com_model, transcriptTable,lambda, maximum_value);

        results = [results; results_table];
    end

    results = cell2table(vertcat(results{:}));
    results.Properties.VariableNames = {'ModelID','Max_growth_rate','Min_growth_rate'};
    writetable(results, fullfile(tablesDir, 'variability_analysis',[methods{i},'_variability_analysis.csv']));
end

