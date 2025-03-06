modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';

methods = {'consensus','carveme','gapseq','kbase'};
timepoint = {'20d', '40d', '60d', '90d', '180d'};

% set up the trade-off parameter
alpha = [0.3,0.5,0.7,0.9, 1];

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);
environment = getEnvironment();

null_v = NaN(14,1);

for i = 1:numel(methods)
    % load models
    model_workspace = fullfile(modelDir, [methods{i},'_com.mat']);
    load(model_workspace)

    for l = 1:numel(alpha)
        micom_cell = cell(numel(timepoint),1);
        micom_ab1_cell = cell(numel(timepoint),1);

        % create the abundance cell and genome ID cell
        ab_cell = cell(numel(timepoint), 1);
        rep_cell = cell(numel(timepoint), 1);
        genome_ID = cell(numel(timepoint), 1);

        parfor j = 1:numel(timepoint)
            restoreEnvironment(environment);
            changeCobraSolver('gurobi','all');

            disp('----------------------------------------------------------------------')
            fprintf('\n################# %s\n\n', timepoint{j})
        
            abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_',timepoint{j}, '.csv']); 
            abTable = readtable(abFile, 'ReadVariableNames', true);

            % check if 'replication_rate' is double or not
            if iscell(abTable.replication_rate)
                abTable.replication_rate = cellfun(@str2double, abTable.replication_rate);
            end

            disp('Calculate Max Growth Rate with coco')
            micom_solution = MICOM_max_growth(com_model, abTable);
            micom_max_growth = -micom_solution.objval;

            disp('Calculate Community Growth Rate with coco')
            com_solution = MICOM(com_model, micom_max_growth, abTable, alpha(l));

            % MICOM without abundance info
            com_solution1 = MICOM_ab1(com_model, 1);

            if isfield(com_solution, 'x') && ~contains(com_solution.status, 'NUMERIC')
                micom_cell{j} = com_solution.x(contains(com_model.rxns,'BIOMASS_R'));
            else
                micom_cell{j} = null_v;
            end

            if isfield(com_solution1, 'x') && ~contains(com_solution1.status, 'NUMERIC')
                micom_ab1_cell{j} = com_solution1.x(contains(com_model.rxns,'BIOMASS_R'));
            else
                micom_ab1_cell{j} = null_v;
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

        micom_results_cells = vertcat(micom_cell{:});
        micom_ab1_results_cells = vertcat(micom_ab1_cell{:});

        ab_info = vertcat(ab_cell{:});
        rep_info = vertcat(rep_cell{:});
        ID_info = vertcat(genome_ID{:});
    
        alphas = {'0.3','0.5','0.7','0.9', '1'};

        results = table(ID_info, ab_info, rep_info, micom_results_cells, micom_ab1_results_cells);
        results.Properties.VariableNames = {'ID', 'MAG_ab', 'Replication_rate', 'MICOM', 'MICOM_ab1'};
        writetable(results, fullfile(tablesDir,'predicted_growth',[methods{i},'_MICOM_results_table_', alphas{l},'.csv']));
    end
end