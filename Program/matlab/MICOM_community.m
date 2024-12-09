modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';

methods = {'consensus','carveme','gapseq','kbase'};
timepoint = {'20d', '40d', '60d', '90d', '180d'};

% set up the trade-off parameter
alpha = 0.5;

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);
environment = getEnvironment();

for i = 1:numel(methods)
    % load models
    model_workspace = fullfile(modelDir, [methods{i},'_com.mat']);
    load(model_workspace)

    micom_cell = cell(numel(timepoint),1);
    micom_ab1_cell = cell(numel(timepoint),1);

    parfor j = 1:numel(timepoint)
        restoreEnvironment(environment);
        changeCobraSolver('gurobi','all');

        disp('----------------------------------------------------------------------')
        fprintf('\n################# %s\n\n', timepoint{j})
        
        abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_',timepoint{j}, '.csv']); 
        abTable = readtable(abFile, 'ReadVariableNames', true);

        disp('Calculate Max Growth Rate with coco')
        micom_solution = MICOM_max_growth(com_model, abTable);
        micom_max_growth = -micom_solution.objval;

        disp('Calculate Community Growth Rate with coco')
        com_solution = MICOM(com_model, micom_max_growth, abTable, alpha);

        % MICOM without abundance info
        com_solution1 = MICOM_ab1(com_model, 1);

        micom_cell{j} = com_solution.x(contains(com_model.rxns,'BIOMASS_R'));
        micom_ab1_cell{j} = com_solution1.x(contains(com_model.rxns,'BIOMASS_R'));

    end

    micom_results_cells = vertcat(micom_cell{:});
    micom_ab1_results_cells = vertcat(micom_ab1_cell{:});

    results = table(micom_results_cells, micom_ab1_results_cells);
    results.Properties.VariableNames = {'MICOM','MICOM_ab1'};
    writetable(results, fullfile(tablesDir,'predicted_growth',[methods{i},'_MICOM_results_table.csv']));
end