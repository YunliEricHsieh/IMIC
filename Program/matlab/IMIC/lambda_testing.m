modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';
timepoint = {'20d', '40d', '60d', '90d', '180d'};
methods = {'consensus','carveme','gapseq','kbase'};

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);
parfevalOnAll(@maxNumCompThreads, 0, 2); 

% testing the different value of lambda
lambda = [0.1,0.5,1,2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,22,24,25,50,75,100];

null_v = zeros(14,1);
null_v(:,1) = NaN;

for i = 1:numel(methods)
    model_workspace = fullfile(modelDir, [methods{i},'_com.mat']);
    load(model_workspace)
    
    results = [];
    
    for j = 1:numel(timepoint)
        disp('----------------------------------------------------------------------')
        fprintf('\n################# %s\n\n', timepoint{j})
    
        tmp_results = [];
        transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{j}, '.csv']);
        transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);

        % integrate the abundance info into the final result tables
        abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_',timepoint{j}, '.csv']); 
        abTable = readtable(abFile, 'ReadVariableNames', true);

        % check if 'replication_rate' is double or not
        if iscell(abTable.replication_rate)
            abTable.replication_rate = cellfun(@str2double, abTable.replication_rate);
        end
    
        parfor k = 1:numel(lambda)
            disp('Calculate Community Growth Rate')
            com_solution = IMIC(com_model, transcriptTable, lambda(k));

            if contains(com_solution.status, 'NUMERIC') 
                tmp_results = [tmp_results, null_v];
            else
                tmp_results = [tmp_results, com_solution.x(contains(com_model.rxns,'BIOMASS_R'))];
            end
        end
    
        results = [results; tmp_results];

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

    ab_info = vertcat(ab_cell{:});
    rep_info = vertcat(rep_cell{:});
    ID_info = vertcat(genome_ID{:});

    table1 = table(ID_info, ab_info, rep_info);
    table1.Properties.VariableNames = {'ID', 'MAG_ab', 'Replication_rate'};

    table2 = array2table(results);
    table2.Properties.VariableNames = {'L = 0,1','L = 0.5','L = 1','L = 2','L = 3','L = 4','L = 5','L = 6', ...
        'L = 7','L = 8','L = 9','L = 10', 'L = 12', ...
        'L = 14', 'L = 15', 'L = 16', 'L = 18', 'L = 20', 'L = 22', ...
        'L = 24', 'L = 25', 'L = 50', 'L = 75', 'L = 100'};

    results_table = [table1, table2];
    writetable(results_table, fullfile(tablesDir, 'parameter_test', [methods{i},'_lambda_test_resutls.csv']));
    clear results_table results tmp_results
end

