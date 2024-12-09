modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';
methods = {'consensus','carveme','gapseq','kbase'};
timepoint = {'20d', '40d', '60d', '90d', '180d'};

<<<<<<< HEAD
lambda = 12;

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);
=======
lambda = [12, 14, 15, 50];
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7

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
<<<<<<< HEAD
        com_solution = IMIC(com_model, transcriptTable, lambda);

        maximum_value = -com_solution.objval;

        results_table = var_analysis(com_model, transcriptTable,lambda, maximum_value);
=======
        com_solution = IMIC(com_model, transcriptTable, lambda(i));

        maximum_value = -com_solution.objval;

        results_table = var_analysis(com_model, transcriptTable,lambda(i), maximum_value);
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7

        results = [results; results_table];
    end

    results = cell2table(results);
    results.Properties.VariableNames = {'ModelID','Max_growth_rate','Min_growth_rate'};
<<<<<<< HEAD
    writetable(results, fullfile(tablesDir, 'variability_analysis',[methods{i},'_variability_analysis.csv']));
=======
    writetable(results, fullfile(tablesDir, 'variability_analysis',[methods{i},'_variability_analysis_with_opt_lambda.csv']));
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7
end




