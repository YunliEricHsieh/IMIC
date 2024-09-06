%%% Sensitivity analysis %%%
modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';
timepoint = {'20d', '40d', '60d', '90d', '180d'};
methods = {'consensus','carveme','gapseq','kbase'};

% lambda for testing
lambda = [0.1,0.5,1,2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,22,24,25,50,75,100];

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);

for k = 1:numel(methods)
    % load models
    model_workspace = fullfile(modelDir, [methods{k},'_com.mat']);
    load(model_workspace)

    results_cells = cell(numel(lambda), 1);

    parfor i = 1:numel(lambda)
        
            solution = zeros(1,numel(timepoint));

            for j = 1:numel(timepoint)
                disp('----------------------------------------------------------------------')
                fprintf('\n################# %s\n\n', timepoint{j})
    
                transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{j}, '.csv']);
                transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);

                com_solution = IMIC(com_model, transcriptTable, lambda(i));

                solution(j) = -com_solution.objval;
            end

            results_cells{i} = solution';
    end

    results_table = vertcat(results_cells{:});
    Lambda = [repmat({'L = 0.1'},5,1);repmat({'L = 0.5'},5,1);repmat({'L = 1'},5,1);repmat({'L = 2'},5,1);
        repmat({'L = 3'},5,1);repmat({'L = 4'},5,1);repmat({'L = 5'},5,1);repmat({'L = 6'},5,1);
        repmat({'L = 7'},5,1);repmat({'L = 8'},5,1);repmat({'L = 9'},5,1);repmat({'L = 10'},5,1);
        repmat({'L = 12'},5,1);repmat({'L = 14'},5,1);repmat({'L = 15'},5,1);repmat({'L = 16'},5,1);
        repmat({'L = 18'},5,1);repmat({'L = 20'},5,1);repmat({'L = 22'},5,1);repmat({'L = 24'},5,1);
        repmat({'L = 25'},5,1);repmat({'L = 50'},5,1);repmat({'L = 75'},5,1);repmat({'L = 100'},5,1)];

    results = table(Lambda,results_table);
    results.Properties.VariableNames = {'Lambda','Solutions'};
    writetable(results, fullfile(tablesDir, 'sensitivity_result', [methods{k},'_sensitivity_analysis.csv']));
end