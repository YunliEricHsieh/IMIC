modelDir =  '~/IMIC/study_case/models';
tablesDir = '~/IMIC/study_case/table';
methods = {'consensus'};
timepoint = {'0h_1', '0h_2', '0h_3', '4h_1', '4h_2','4h_3','8h_1', '8h_2','8h_3','24h_1', '24h_2','24h_3'};
ratio = {'1-1_', '1000-1_', '1-1000_'};
output = {'equal', 'more', 'less'};

lambda = 9;

null_v = zeros(2,1);
null_v(:,1) = NaN;

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);

model_workspace = fullfile(modelDir, [methods{1},'_com.mat']);
load(model_workspace)

for i = 1:numel(ratio)

    results = [];
    
    for j = 1:numel(timepoint)
        disp('----------------------------------------------------------------------')
        fprintf('\n################# %s\n\n', timepoint{j})
        
        tmp_results = [];
        transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', ratio{i}, timepoint{j}, '.csv']);
        transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);
            
        disp('Calculate Community Growth Rate')
        com_solution = IMIC(com_model, transcriptTable, lambda);

        if contains(com_solution.status, 'NUMERIC') 
            tmp_results = [tmp_results, null_v];
        else
            tmp_results = [tmp_results, com_solution.x(contains(com_model.rxns,'BIOMASS_R'))];
        end
        
        results = [results; tmp_results];

        clear tmp_results
    end
    
    Species = repmat({'PP';'EC'},12,1);
    Time = [repmat({'0h'},6,1);repmat({'4h'},6,1);repmat({'8h'},6,1);repmat({'24h'},6,1)]; 

    results_table = table(Species,results,Time);
    results_table.Properties.VariableNames = {'Species','Growth_rate','Timepoint'};
    writetable(results_table, fullfile(tablesDir, 'predicted_growth', [output{i},'_results_table.csv']));

end

