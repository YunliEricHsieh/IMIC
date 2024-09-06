modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';
methods = {'consensus','carveme','gapseq','kbase'};
timepoint = {'20d', '40d', '60d', '90d', '180d'};

% set up the parameters
delta = 10; % CoCo
gamma = 30; % CoCo
alpha = 0.5; % CoCo
lambda = 12; % IMIC

null_v = zeros(14,1);
null_v(:,1) = NaN;

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);

for i = 1:numel(methods)
    model_workspace = fullfile(modelDir, [methods{i},'_com.mat']);
    load(model_workspace)
    
    coco_cell = cell(numel(timepoint), 1);
    imic_cell = cell(numel(timepoint), 1);
    
    parfor j = 1:numel(timepoint)
        disp('----------------------------------------------------------------------')
        fprintf('\n################# %s\n\n', timepoint{j})
        
        abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_',timepoint{j}, '.csv']);
        transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{j}, '.csv']);
    
        abTable = readtable(abFile, 'ReadVariableNames', true);
        transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);
        
        disp('Calculate Max Growth Rate with coco')
        micom_solution = MICOM_max_growth(com_model, abTable);
        coco_max_growth = -micom_solution.objval;


        disp('Calculate Community Growth Rate with coco')
        coco_com_solution = coco(com_model, coco_max_growth, ...
            transcriptTable, abTable, alpha, delta, gamma);

        if isfield(coco_com_solution, 'x') && ~contains(coco_com_solution.status, 'NUMERIC')
            coco_cell{j} = coco_com_solution.x(contains(com_model.rxns,'BIOMASS_R'));
        else
            coco_cell{j} = null_v;
        end
        
        disp('Calculate Community Growth Rate')
        com_solution = IMIC(com_model, transcriptTable, lambda);
                
        if isfield(com_solution, 'x') && ~contains(com_solution.status, 'NUMERIC')
            imic_cell{j} = com_solution.x(contains(com_model.rxns,'BIOMASS_R'));
        else
            imic_cell{j} = null_v;
        end
                
    end
    
    coco_results_cells = vertcat(coco_cell{:});
    imic_results_cells = vertcat(imic_cell{:});

    results = table(coco_results_cells,imic_results_cells);
    results.Properties.VariableNames = {'CoCo','IMIC'};
    writetable(results, fullfile(tablesDir, 'predicted_growth', [methods{i},'_results_table.csv']));

    clear results
end

