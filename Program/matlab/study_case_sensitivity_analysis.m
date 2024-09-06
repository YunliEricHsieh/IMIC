modelDir =  '~/IMIC/study_case/models';
tablesDir = '~/IMIC/study_case/table';
methods = {'consensus'};
timepoint = {'0h_1', '0h_2', '0h_3', '4h_1', '4h_2','4h_3','8h_1', '8h_2','8h_3','24h_1', '24h_2','24h_3'};
ratio = {'1-1_','1000-1_','1-1000_'};
output = {'equal','more','less'};

% lambda for testing
lambda = [0.1,0.5,1,2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,22,24,25,50,75,100];

% load models
model_workspace = fullfile(modelDir, [methods{1},'_com.mat']);
load(model_workspace)

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);

for k = 1:numel(ratio)

    results_cells = cell(numel(lambda), 1);

    parfor i = 1:numel(lambda)
        
            solution = zeros(1,numel(timepoint));

            for j = 1:numel(timepoint)
                disp('----------------------------------------------------------------------')
                fprintf('\n################# %s\n\n', timepoint{j})
    
                transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_',ratio{k}, timepoint{j}, '.csv']);
                transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);

                com_solution = IMIC(com_model, transcriptTable, lambda(i));

                solution(j) = -com_solution.objval;
            end

            results_cells{i} = solution';
    end

    results_table = vertcat(results_cells{:});
    Lambda = [repmat({'L = 0.1'},12,1);repmat({'L = 0.5'},12,1);repmat({'L = 1'},12,1);repmat({'L = 2'},12,1);
        repmat({'L = 3'},12,1);repmat({'L = 4'},12,1);repmat({'L = 5'},12,1);repmat({'L = 6'},12,1);
        repmat({'L = 7'},12,1);repmat({'L = 8'},12,1);repmat({'L = 9'},12,1);repmat({'L = 10'},12,1);
        repmat({'L = 12'},12,1);repmat({'L = 14'},12,1);repmat({'L = 15'},12,1);repmat({'L = 16'},12,1);
        repmat({'L = 18'},12,1);repmat({'L = 20'},12,1);repmat({'L = 22'},12,1);repmat({'L = 24'},12,1);
        repmat({'L = 25'},12,1);repmat({'L = 50'},12,1);repmat({'L = 75'},12,1);repmat({'L = 100'},12,1)];
    Time = repmat([repmat({'0h'},3,1);repmat({'4h'},3,1);repmat({'8h'},3,1);repmat({'24h'},3,1)],24,1);

    results = table(Lambda,results_table, Time);
    results.Properties.VariableNames = {'Lambda','Solutions','Time'};
    writetable(results, fullfile(tablesDir, 'sensitivity_result', [output{k},'_sensitivity_analysis.csv']));
end



