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
        
            solution_cell = cell(numel(timepoint),1);

            for j = 1:numel(timepoint)
                disp('----------------------------------------------------------------------')
                fprintf('\n################# %s\n\n', timepoint{j})
    
                transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_',ratio{k}, timepoint{j}, '.csv']);
                transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);

                com_solution = IMIC(com_model, transcriptTable, lambda(i));

                solution_cell{j} = com_solution.x(contains(com_model.rxns,'BIOMASS_R'));
            end

            results_cells{i} = vertcat(solution_cell{:});
    end
    
    results_table = vertcat(results_cells{:});
    Species = repmat({'EC';'PP'},288,1);
    Time = repmat([repmat({'0h'},6,1);repmat({'4h'},6,1);repmat({'8h'},6,1);repmat({'24h'},6,1)],24,1);
    Lambda = [repmat({'L = 0.1'},24,1);repmat({'L = 0.5'},24,1);repmat({'L = 1'},24,1);repmat({'L = 2'},24,1);
        repmat({'L = 3'},24,1);repmat({'L = 4'},24,1);repmat({'L = 5'},24,1);repmat({'L = 6'},24,1);
        repmat({'L = 7'},24,1);repmat({'L = 8'},24,1);repmat({'L = 9'},24,1);repmat({'L = 10'},24,1);
        repmat({'L = 12'},24,1);repmat({'L = 14'},24,1);repmat({'L = 15'},24,1);repmat({'L = 16'},24,1);
        repmat({'L = 18'},24,1);repmat({'L = 20'},24,1);repmat({'L = 22'},24,1);repmat({'L = 24'},24,1);
        repmat({'L = 25'},24,1);repmat({'L = 50'},24,1);repmat({'L = 75'},24,1);repmat({'L = 100'},24,1)];

    results = table(Species,results_table,Time,Lambda);
    results.Properties.VariableNames = {'Species','Growth_rate','Time','Lambda'};
    writetable(results, fullfile(tablesDir, 'parameter_test', [output{k},'_lambda_testing.csv']));

end



