modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';

timepoint = {'20d', '40d', '60d', '90d', '180d'};

% read the community model
model_workspace = fullfile(modelDir, 'consensus_com.mat');
load(model_workspace);

% set up the parameters
delta = 30; % CoCo
gamma = 100; % CoCo
alpha = 0.3; % CoCo
lambda = 12; % IMIC

ncpu = 5;
delete(gcp('nocreate'));
parpool(ncpu);
environment = getEnvironment();

% MICOM
flux_sum_micom_ab1 = [];
metabolite_ID_micom_ab1 = {};
flux_sum_micom = [];
metabolite_ID_micom = {};

% CoCo-GEM
flux_sum_coco_ab1 = [];
metabolite_ID_coco_ab1 = {};
flux_sum_coco = [];
metabolite_ID_coco = {};

% IMIC
flux_sum = [];
metabolite_ID = {};

tar_met_IDs = [];
  
parfor i = 1:numel(timepoint)

    restoreEnvironment(environment);
    changeCobraSolver('gurobi','all');

    disp('----------------------------------------------------------------------')
    fprintf('\n################# %s\n\n', timepoint{i})

    abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_',timepoint{i}, '.csv']);
    transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{i}, '.csv']);
    
    abTable = readtable(abFile, 'ReadVariableNames', true);
    transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);
    
    % coco ab1 flux sum
    coco_ab1_solution = coco_max_growth_ab1(com_model, abTable);
    coco_ab1_max_growth = -coco_ab1_solution.objval;

    [coco_ab1_flux_sum_value, coco_ab1_mets_ID] = coco_ab1_flux_sum(com_model, coco_ab1_max_growth, ...
        transcriptTable, abTable, alpha, delta, gamma, tar_met_IDs);
    
    flux_sum_coco_ab1(:, i) = coco_ab1_flux_sum_value;
    metabolite_ID_coco_ab1{i} = coco_ab1_mets_ID;
    
    % coco flux sum
    coc_solution = MICOM_max_growth(com_model, abTable);
    coco_max_growth = -coc_solution.objval;
    
    [coco_flux_sum_value, coco_mets_ID] = coco_flux_sum(com_model, coco_max_growth, ...
        transcriptTable, abTable, alpha, delta, gamma, tar_met_IDs);

    flux_sum_coco(:, i) = coco_flux_sum_value;
    metabolite_ID_coco{i} = coco_mets_ID;

    % MICOM ab1 flux sum
    [micom_ab1_flux_sum_value, micom_ab1_mets_ID] = MICOM_ab1_flux_sum(com_model, alpha, tar_met_IDs);

    flux_sum_micom_ab1(:, i) = micom_ab1_flux_sum_value;
    metabolite_ID_micom_ab1{i} = micom_ab1_mets_ID;

    % MICOM flux sum
    micom_solution = MICOM_max_growth(com_model, abTable);
    micom_max_growth = -micom_solution.objval;

    [micom_flux_sum_value, micom_mets_ID]  = MICOM_flux_sum(model, micom_max_growth, abTable, alpha, tar_met_IDs);

    flux_sum_micom(:, i) = micom_flux_sum_value;
    metabolite_ID_micom{i} = micom_mets_ID;

    % IMIC
    [flux_sum_value, mets_ID] = IMIC_flux_sum(com_model, transcriptTable, lambda, tar_met_IDs);

    flux_sum(:, i) = flux_sum_value;
    metabolite_ID{i} = mets_ID;

end

% write the table after calculating for one time point
out_transport = vertcat(metabolite_ID{1});
results = table(out_transport,flux_sum);
writetable(results, fullfile(tablesDir,'flux_sum_analysis','14_MAG_media_mets_flux.csv'));
    
out_transport = vertcat(metabolite_ID_coco{1});
coco_results = table(out_transport,flux_sum_coco);
writetable(coco_results, fullfile(tablesDir,'flux_sum_analysis','coco_14_MAG_media_mets_flux.csv'));

out_transport = vertcat(metabolite_ID_coco_ab1{1});
coco_ab1_results = table(out_transport,flux_sum_coco_ab1);
writetable(coco_ab1_results, fullfile(tablesDir,'flux_sum_analysis','coco_ab1_14_MAG_media_mets_flux.csv'));    

out_transport = vertcat(metabolite_ID_micom{1});
micom_results = table(out_transport,flux_sum_micom);
writetable(micom_results, fullfile(tablesDir,'flux_sum_analysis','MICOM_14_MAG_media_mets_flux.csv'));
    
out_transport = vertcat(metabolite_ID_micom_ab1{1});
micom_ab1_results = table(out_transport,flux_sum_micom_ab1);
writetable(micom_ab1_results, fullfile(tablesDir,'flux_sum_analysis','MICOM_ab1_14_MAG_media_mets_flux.csv'));
    


