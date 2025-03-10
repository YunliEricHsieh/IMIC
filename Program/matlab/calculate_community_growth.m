modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';
methods = {'consensus','carveme','gapseq','kbase'};
timepoint = {'20d', '40d', '60d', '90d', '180d'};

% set up the parameters
delta = [30, 100, 50, 10]; % CoCo for each methods
gamma = [100, 10, 10, 10]; % CoCo for each methods
alpha = [0.3, 0.5, 0.5, 0.5]; % CoCo for each methods
lambda = 12; % IMIC

null_v = zeros(14,1);
null_v(:,1) = NaN;

ncpu = 5;
delete(gcp('nocreate'));
parpool(ncpu);

for i = 1:numel(methods)
    model_workspace = fullfile(modelDir, [methods{i},'_com.mat']);
    load(model_workspace)
    
    coco_cell = cell(numel(timepoint), 1);
    imic_cell = cell(numel(timepoint), 1);

    % create the abundance cell and genome ID cell
    ab_cell = cell(numel(timepoint), 1);
    rep_cell = cell(numel(timepoint), 1);
    genome_ID = cell(numel(timepoint), 1);
    
    parfor j = 1:numel(timepoint)
        disp('----------------------------------------------------------------------')
        fprintf('\n################# %s\n\n', timepoint{j})
        
        abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_',timepoint{j}, '.csv']);
        transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{j}, '.csv']);
    
        abTable = readtable(abFile, 'ReadVariableNames', true);
        transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);

        % check if 'replication_rate' is double or not
        if iscell(abTable.replication_rate)
            abTable.replication_rate = cellfun(@str2double, abTable.replication_rate);
        end
        
        disp('Calculate Max Growth Rate with coco')
        micom_solution = MICOM_max_growth(com_model, abTable);
        coco_max_growth = -micom_solution.objval;


        disp('Calculate Community Growth Rate with coco')
        coco_com_solution = coco(com_model, coco_max_growth, ...
            transcriptTable, abTable, alpha(i), delta(i), gamma(i));

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

        % find the biomass reactions
        bio_rxn = com_model.rxns(contains(com_model.rxns, 'BIOMASS_R'));
        ab_table = [];
        rep_table = [];
        genome_table = [];

        for k = 1:numel(bio_rxn)
            num = extractAfter(bio_rxn{k}, 'BIOMASS_Reaction_');
            ab = abTable.relative_ab(find(contains(abTable.Genome,['KG',num,'_genomic'])));
            rep = abTable.replication_rate(find(contains(abTable.Genome,['KG',num,'_genomic'])));
            ID = abTable.Genome(find(contains(abTable.Genome,['KG',num,'_genomic'])));

            ab_table = [ab_table; ab];
            rep_table = [rep_table; rep];
            genome_table = [genome_table; ID];
        end
        
        ab_cell{j,1} = ab_table;
        rep_cell{j,1} = rep_table;
        genome_ID{j,1} = genome_table;
                
    end
    
    coco_results_cells = vertcat(coco_cell{:});
    imic_results_cells = vertcat(imic_cell{:});

    ab_info = vertcat(ab_cell{:});
    rep_info = vertcat(rep_cell{:});
    ID_info = vertcat(genome_ID{:});

    results = table(ID_info, ab_info, rep_info, coco_results_cells, imic_results_cells);
    results.Properties.VariableNames = {'ID', 'MAG_ab', 'Replication_rate', 'CoCo', 'IMIC'};
    writetable(results, fullfile(tablesDir, 'predicted_growth', [methods{i},'_results_table.csv']));

    clear results
end

