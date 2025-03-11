modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';

timepoint = {'20d', '40d', '60d', '90d', '180d'};

model_workspace = fullfile(modelDir, 'iterative/consensus.mat');
load(model_workspace, 'GF');

% find all the common reactions across all MAGs
comm_rxns = intersect(GF{1}.rxns,GF{2}.rxns);

for i = 3:14
    comm_rxns = intersect(comm_rxns,GF{i}.rxns);
end

% read the community model
model_workspace = fullfile(modelDir, 'consensus_com.mat');
load(model_workspace);

% sellect the reaction with GPR rules
rxn_list = [];

for i = 1:numel(comm_rxns)
   search_list = com_model.grRules(contains(com_model.rxns,comm_rxns{i}));

   if sum(~cellfun(@(x) ischar(x) && isempty(x), search_list)) > 0
       rxn_list = [rxn_list; comm_rxns(i)];
   end
end

% calculate gene expression abundance of reaction by using GPR rules
flux_cor = zeros(numel(rxn_list),5);
flux_cor_p = zeros(numel(rxn_list),5);
rxns = com_model.rxns;
rxns(find(contains(rxns, '_r'))) = {'RE'};
rxn_flux_cell = cell(numel(rxn_list), 5);
ab_cell = cell(numel(rxn_list), 5);
expression_value_cell = cell(numel(rxn_list),5);

for i = 1:numel(timepoint)
    abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_',timepoint{i}, '.csv']);
    transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{i}, '.csv']);
    
    abTable = readtable(abFile, 'ReadVariableNames', true);
    transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);

    com_solution = IMIC(com_model, transcriptTable, 12);

    % calculate f(g)
    express_value = zeros(numel(com_model.rxns),1);
    
    for j = 1:numel(com_model.rxns)
        if ~isempty(com_model.rules{j}) && ~contains(com_model.grRules{j},'spontaneous')
            % split the genes by 'or'
            gpr_split = strsplit(com_model.rules{j},'|');
            gpr_split_TPM = [];

            for k = numel(gpr_split)
                list = extractBetween(gpr_split{k},'x(',')');
                list_TPM = [];

                for l = 1:numel(list)
                    gene_id = com_model.genes(str2num(list{l}));
                    list_TPM(l) = transcriptTable.TPM(strcmp(transcriptTable.Geneid,gene_id));
                end
            
                % find min value
                gpr_split_TPM(k) = min(list_TPM);
            end

            % sum the value
            express_value(j) = sum(gpr_split_TPM);
        end
    end

    % find the flux value and f(g) for common reaction, which having GPR
    % rule
    for j = 1:numel(rxn_list)
        expression_value_cell{j,i} = express_value(contains(rxns, rxn_list{j}));

        rxn_flux_cell{j,i} = com_solution.x(contains(rxns, rxn_list{j}));

        tem_rxn = rxns(contains(rxns, rxn_list{j}));

        ab_table = [];

        for k = 1:numel(tem_rxn)
            num = extractAfter(tem_rxn{k}, rxn_list{j});
            ab = abTable.relative_ab(find(contains(abTable.Genome,['KG',num,'_genomic'])));
            ab_table = [ab_table; ab];
        end
        
        ab_cell{j,i} = ab_table;
    end
    
    % do the Spearman correlation analysis
    for j = 1:numel(rxn_list)
        [flux_cor(j,i), flux_cor_p(j,i)] = corr(rxn_flux_cell{j,i}, ab_cell{j,i},'type','Spearman');
    end
end

% find the reaction which having high correlation between flux and relative
% abundance
tmp_rxn = cell(5,1);

for i = 1:numel(timepoint)

    % correct p-value for different time points
    [c_pvalues, c_alpha, h] = fdr_BH(flux_cor_p(:,i), 0.05);

    % find h = 1 (H0 is rejected) and correlation coefficient > 0.7 
    tmp_rxn{i} = rxn_list(find(h == 1 & flux_cor(:,i) > 0.7));
end

% find the intersection over all timepoint (except for timepoint 4)
cor_reaction = intersect(tmp_rxn{1}, tmp_rxn{2});
cor_reaction = intersect(cor_reaction, tmp_rxn{3});
cor_reaction = intersect(cor_reaction, tmp_rxn{5});

% save the list of core reactions
rxn_table = table(cor_reaction);
rxn_table.Properties.VariableNames = {'Core_reactions'};
writetable(rxn_table, fullfile(tablesDir, 'Core_reactions_list.csv'));

% find the fluxes of core reactions in each timepoint
ids = find(contains(rxn_list, cor_reaction));

rxn_flux_1 = rxn_flux_cell(ids,1);
rxn_flux_2 = rxn_flux_cell(ids,2);
rxn_flux_3 = rxn_flux_cell(ids,3);
rxn_flux_4 = rxn_flux_cell(ids,4);
rxn_flux_5 = rxn_flux_cell(ids,5);

% find the abundance of MAG in each timepoint
ab_info_1 = ab_cell(ids,1);
ab_info_2 = ab_cell(ids,2);
ab_info_3 = ab_cell(ids,3);
ab_info_4 = ab_cell(ids,4);
ab_info_5 = ab_cell(ids,5);

% find the correaction factor, f(g), of gene expression for the core reactions
expression_1 = expression_value_cell(ids,1);
expression_2 = expression_value_cell(ids,2);
expression_3 = expression_value_cell(ids,3);
expression_4 = expression_value_cell(ids,4);
expression_5 = expression_value_cell(ids,5);

RxnID = repelem(cor_reaction, 14);

% save the table
result_table = table(vertcat(RxnID{:}),vertcat(ab_info_1{:}),vertcat(ab_info_2{:}),vertcat(ab_info_3{:}),vertcat(ab_info_4{:}),vertcat(ab_info_5{:}), ...
    vertcat(rxn_flux_1{:}),vertcat(rxn_flux_2{:}),vertcat(rxn_flux_3{:}),vertcat(rxn_flux_4{:}),vertcat(rxn_flux_5{:}), ...
    vertcat(expression_1{:}),vertcat(expression_2{:}),vertcat(expression_3{:}),vertcat(expression_4{:}),vertcat(expression_5{:}));
result_table.Properties.VariableNames = {'RxnID','Ab_info_1','Ab_info_2','Ab_info_3','Ab_info_4','Ab_info_5', ...
    'rxn_flux_1','rxn_flux_2','rxn_flux_3','rxn_flux_4','rxn_flux_5',...
    'expression_1','expression_2','expression_3','expression_4','expression_5'};

writetable(result_table, fullfile(tablesDir,'key_reaction','flux_value_of_key_rxns.csv'));

%% testing the impact of these reaction on community growth
ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);

vari_sum = zeros(5,numel(cor_reaction));

for i = 1:numel(timepoint)
    
    % load metatranscriptomic data
    transcriptFile = fullfile(tablesDir, 'abundance_table', ['rna_ab_and_geneID_', timepoint{i}, '.csv']);
    transcriptTable = readtable(transcriptFile, 'ReadVariableNames', true);

    % IMIC
    ori_solution = IMIC(com_model, transcriptTable, 12);
    sum_value = sum(ori_solution.x(contains(com_model.rxns, 'BIOMASS_R')));

    parfor j = 1:numel(cor_reaction)

        % find the IDs of knocking out reactions
        knock_out = com_model.rxns(contains(com_model.rxns, cor_reaction{j}));
        k_model = removeRxns(com_model, knock_out);

        % IMIC
        k_solution = IMIC(k_model, transcriptTable, 12);

        % calculate the variability of sum of individual growth rate
        tmp_sum = sum(k_solution.x(contains(k_model.rxns, 'BIOMASS_R')));
        vari_sum(i,j) = (sum_value - tmp_sum)/sum_value;
    end

end

writematrix(vari_sum, fullfile(tablesDir, 'key_reaction','the impact of core reaction on community growth rate.csv'));