modelDir =  '~/IMIC/models';
<<<<<<< HEAD
model_workspace = fullfile(modelDir, 'consensus_com.mat');
load(model_workspace)
=======
methods = {'consensus','carveme','gapseq','kbase'};
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7

tablesDir = '~/IMIC/table';
timepoint = {'20d', '40d', '60d', '90d', '180d'};

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);

% testing different parameter value
<<<<<<< HEAD
alpha = [0.3,0.5,0.7,0.9];
=======
%alpha = [0.3,0.5,0.7,0.9];
alpha = [0.5];
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7
delta = [10,20,30,40,50,60,70,80,90,100];
gamma = [10,20,30,40,50,60,70,80,90,100];

null_v = NaN(14,1);

<<<<<<< HEAD
for l = 1:numel(alpha)

    results = [];

    for i = 1:numel(timepoint)
        disp('----------------------------------------------------------------------')
        fprintf('\n################# %s\n\n', timepoint{i});
    
        abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_', timepoint{i}, '.csv']);
        transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{i}, '.csv']);
    
        abTable = readtable(abFile, 'ReadVariableNames', true);
        transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);
    
        disp('Calculate Max Growth Rate with coco')
        coc_solution = MICOM_max_growth(com_model, abTable);
        coco_max_growth = -coc_solution.objval;
    
        % Calculate the total combinations of delta and gamma
        totalCombinations = numel(delta) * numel(gamma);
        tmp_results = zeros(14, totalCombinations); % Adjusted to preallocate the array

        parfor idx = 1:totalCombinations
            [j, k] = ind2sub([numel(delta), numel(gamma)], idx); % Convert linear index to subscript
        
            fprintf('Calculating for delta = %d and gamma = %d...\n', delta(k), gamma(j));
            coco_com_solution = coco(com_model, coco_max_growth, ...
                transcriptTable, abTable, alpha(l), delta(k), gamma(j));

            if isfield(coco_com_solution, 'x') && ~contains(coco_com_solution.status, 'NUMERIC')
                tmp_results(:, idx) = coco_com_solution.x(contains(com_model.rxns,'BIOMASS_R'));
            else
                tmp_results(:, idx) = null_v;
            end
        end

        % Flatten and concatenate the tmp_results to the main results
        results = [results; tmp_results]; 
    end

    alphas = {'0.3','0.5','0.7','0.9'};
    % Save the final results
    writematrix(results, fullfile(tablesDir, 'parameter_test', ['coco_test_alpha_', alphas{l},'.csv']));

=======
for m = 1:numel(methods)
    model_workspace = fullfile(modelDir, [methods{m},'_com.mat']);
    load(model_workspace)

    results = [];

    for l = 1:numel(alpha)

        for i = 1:numel(timepoint)
            disp('----------------------------------------------------------------------')
            fprintf('\n################# %s\n\n', timepoint{i});
    
            abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_', timepoint{i}, '.csv']);
            transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{i}, '.csv']);
    
            abTable = readtable(abFile, 'ReadVariableNames', true);
            transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);
    
            disp('Calculate Max Growth Rate with coco')
            coc_solution = MICOM_max_growth(com_model, abTable);
            coco_max_growth = -coc_solution.objval;
    
            % Calculate the total combinations of delta and gamma
            totalCombinations = numel(delta) * numel(gamma);
            tmp_results = zeros(14, totalCombinations); % Adjusted to preallocate the array

            parfor idx = 1:totalCombinations
                [j, k] = ind2sub([numel(delta), numel(gamma)], idx); % Convert linear index to subscript
        
                fprintf('Calculating for delta = %d and gamma = %d...\n', delta(k), gamma(j));
                coco_com_solution = coco(com_model, coco_max_growth, ...
                    transcriptTable, abTable, alpha(l), delta(k), gamma(j));

                if isfield(coco_com_solution, 'x') && ~contains(coco_com_solution.status, 'NUMERIC')
                    tmp_results(:, idx) = coco_com_solution.x(contains(com_model.rxns,'BIOMASS_R'));
                else
                    tmp_results(:, idx) = null_v;
                end
            end

            % Flatten and concatenate the tmp_results to the main results
            results = [results; tmp_results]; 
        end

        %alphas = {'0.3','0.5','0.7','0.9'};
        alphas = {'0.5'};
        % Save the final results
        writematrix(results, fullfile(tablesDir, 'parameter_test', [methods{m},'_coco_test_alpha_', alphas{l},'.csv']));
    end
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7
end
