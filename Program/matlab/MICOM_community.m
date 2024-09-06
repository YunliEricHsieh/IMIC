modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';

timepoint = {'20d', '40d', '60d', '90d', '180d'};

% set up the trade-off parameter
alpha = 0.5;

% load models
model_workspace = fullfile(modelDir, 'consensus_com.mat');
load(model_workspace)

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);

micom_cell = cell(numel(timepoint),1);

parfor j = 1:numel(timepoint)
    disp('----------------------------------------------------------------------')
    fprintf('\n################# %s\n\n', timepoint{j})
        
    abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_',timepoint{j}, '.csv']); 
    abTable = readtable(abFile, 'ReadVariableNames', true);

    disp('Calculate Max Growth Rate with coco')
    micom_solution = MICOM_max_growth(com_model, abTable);
    micom_max_growth = -micom_solution.objval;

    disp('Calculate Community Growth Rate with coco')
    com_solution = MICOM(com_model, coco_max_growth, abTable, alpha);

    micom_cell{j} = com_solution.x(contains(com_model.rxns,'BIOMASS_R'));

end

micom_results_cells = vertcat(micom_cell{:});

results = table(micom_results_cells);
results.Properties.VariableNames = {'MICOM'};
writetable(results, fullfile(tablesDir, 'predicted_growth','MICOM_results_table.csv'));