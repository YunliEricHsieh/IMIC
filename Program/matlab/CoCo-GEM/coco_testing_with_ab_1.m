modelDir =  '~/IMIC/models';
methods = {'consensus','carveme','gapseq','kbase'};

tablesDir = '~/IMIC/table';
timepoint = {'20d', '40d', '60d', '90d', '180d'};

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);

alpha = 1;
delta = [10,20,30,40,50,60,70,80,90,100];
gamma = [10,20,30,40,50,60,70,80,90,100];

null_v = NaN(14,1);

for l = 1:numel(methods)
    model_workspace = fullfile(modelDir, [methods{l},'_com.mat']);
    load(model_workspace)

    results = [];

    for i = 1:numel(timepoint)
        disp('----------------------------------------------------------------------')
        fprintf('\n################# %s\n\n', timepoint{i});
    
        abFile = fullfile(tablesDir,'abundance_table', ['relative_ab_', timepoint{i}, '.csv']);
        transcriptFile = fullfile(tablesDir,'abundance_table', ['rna_ab_and_geneID_', timepoint{i}, '.csv']);
    
        abTable = readtable(abFile, 'ReadVariableNames', true);
        transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);
    
        disp('Calculate Max Growth Rate with coco')
        coco_solution = coco_max_growth_ab1(com_model, abTable);
        coco_max_growth = -coco_solution.objval;

        % Calculate the total combinations of delta and gamma
        totalCombinations = numel(delta) * numel(gamma);
        tmp_results = zeros(14, totalCombinations); % Adjusted to preallocate the array

        parfor idx = 1:totalCombinations
            [j, k] = ind2sub([numel(delta), numel(gamma)], idx); % Convert linear index to subscript
        
            fprintf('Calculating for delta = %d and gamma = %d...\n', delta(k), gamma(j));
            coco_com_solution = coco_ab1(com_model, coco_max_growth, ...
                transcriptTable, abTable, alpha, delta(k), gamma(j));

            if isfield(coco_com_solution, 'x') && ~contains(coco_com_solution.status, 'NUMERIC')
                tmp_results(:, idx) = coco_com_solution.x(contains(com_model.rxns,'BIOMASS_R')); 
            else
                tmp_results(:, idx) = null_v;
            end
        end

        % Flatten and concatenate the tmp_results to the main results
        results = [results; tmp_results]; 
    end

    writematrix(results,  fullfile(tablesDir, 'parameter_test', [methods{l},'_coco_test_with_abundance_1.csv']));

end
