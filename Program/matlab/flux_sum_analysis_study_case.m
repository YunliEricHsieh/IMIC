%% Flux sum analysis for each metabolite in [e]
modelDir =  '~/IMIC/study_case/models';
tablesDir = '~/IMIC/study_case/table';

timepoint = {'0h_1', '0h_2', '0h_3', '4h_1', '4h_2','4h_3','8h_1', '8h_2','8h_3','24h_1', '24h_2','24h_3'};
ratio = {'1-1_','1000-1_','1-1000_'};
output = {'equal','more','less'};

% load models
model_workspace = fullfile(modelDir, 'consensus_com.mat');
load(model_workspace)

lambda = 6;

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);

tar_met_IDs = [];

for k = 1:numel(ratio)

    flux_sum_min = [];
    metabolite_ID = {};
  
    parfor i = 1:numel(timepoint)

        disp('----------------------------------------------------------------------')
        fprintf('\n################# %s\n\n', timepoint{i})
    
        transcriptFile = fullfile(tablesDir, 'abundance_table',['rna_ab_and_geneID_',ratio{k}, timepoint{i}, '.csv']);
        transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);

        disp('Calculate Community Growth Rate')
        com_solution = IMIC(com_model, transcriptTable, lambda);

        maximum_value = -com_solution.objval;

        % flux sum analysis
        [min_flux_sum, mets_ID] = flux_sum(com_model, transcriptTable, lambda, maximum_value, tar_met_IDs);

        flux_sum_min(:, i) = min_flux_sum;
        metabolite_ID{i} = mets_ID;
    end
    
    out_transport = vertcat(metabolite_ID{1});
    results = table(out_transport,flux_sum_min);
    writetable(results, fullfile(tablesDir, 'flux_sum_analysis', [output{k},'_flux_sum.csv']));

end

