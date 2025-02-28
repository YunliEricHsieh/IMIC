function com_solution = coco_ab1(model, max_growth, transcriptTable, abTable, alpha, delta, gamma)
% Function for doing CoCo-GEM wiht abundance equal 1.
% Input:
%       model:                  community metabolic models as struct object
%       max_growth:             the maxmium community growth rate calculate
%                               from MICOM_max_growth function
%       transcriptTable:        table containing the abundance of
%                               transcriptomic data from each organism
%       abTable:                table containing the relative abundance of
%                               each organism
%       alpha:                  user-specified trade-off
%       delta:                  user-specified parameter
%       gamma:                  user-specified parameter
% Output:
%       com_solution:           Solution, returned as a real vector or real 
%                               array.


% find model ID number
num = extractBetween(abTable.Genome,'KG','_');

% calculate expression fold change
MIN_EXP_SCALE = 10;

% correct null expression values using 1/MIN_EXP_SCALE of minimum expression value for each gene
log_express = transcriptTable.log_normalized;
null_express = log_express == 0;
log_express(null_express) = max(log_express);
log_express(null_express) = min(log_express)/MIN_EXP_SCALE;

fold_change = log_express/mean(log_express);

% transcript ratio
gene_list = unique(extractBefore(transcriptTable.Geneid,'_'));

for i = 1:numel(gene_list)
    count_sum(i) = sum(transcriptTable.all_99_bam(contains(transcriptTable.Geneid, gene_list{i})));
end

max_count_sum = max(count_sum);

alpha_matrix = count_sum/max_count_sum;

count_log_matrix = log(count_sum + 1);

% define lower and upper bounds:
ub = model.ub;
lb = model.lb;

% calculate gene set expression by fold change data
for i = 1:numel(model.rxns)
    % reactions have GPR
    if ~isempty(model.rules{i}) && ~contains(model.grRules{i},'spontaneous')

        % split the genes by 'or'
        gpr_split = strsplit(model.rules{i},'|');
        gpr_split_fold_change = [];

        for j = numel(gpr_split)
            list = extractBetween(gpr_split{j},'x(',')');
            list_fold_change = [];

            for k = 1:numel(list)
                gene_id = model.genes(str2num(list{k}));
                list_fold_change(k) = fold_change(strcmp(transcriptTable.Geneid,gene_id));
            end
            
            % find min for 'and'
            gpr_split_fold_change(j) = min(list_fold_change);
        end

        % find max for 'or'
        rxn_expr = max(gpr_split_fold_change);

        % new ub for reactions
        deltas_up = delta * count_log_matrix(contains(gene_list, extractBefore(gene_id,'_')));
        alphas = alpha_matrix(contains(gene_list, extractBefore(gene_id,'_')));

        if rxn_expr >= 1
            factors = alphas * ...
                (1 + gamma * alphas * log(rxn_expr)) + (1 - alphas);
        else
            factors = alphas / ...
                (1 + gamma * alphas * abs(log(rxn_expr))) + (1 - alphas);
        end
           
        ub(i) = deltas_up * factors;
    end
end

% export reactions
export = model.rxns(contains(model.rxns, 'export_'));

for i = 1:numel(export)
    
    t = char(export{i});
    inx = strfind(t,'_');
    n = t(inx(2)+1:end);

    ex = find(strcmp(model.rxns,export{i}));
    ub(ex) = ub(ex) * 1; % abundance = 1
    
    clear ex t inx n 
end

% define the objective:
f = zeros(size(model.S,2),1);

for i = 1:numel(num)
    bio = find(ismember(model.rxns,['BIOMASS_Reaction_',num{i}]));
    % set up minimum growth rate
    lb(bio) = 0.001;
    clear bio
end

% equality constraints:
% LHS
Aeq = model.S;
% RHS
beq = zeros(size(model.S,1),1);

% constraint of max growth rate
ab = zeros(size(Aeq,2),1)';
b = zeros(1,1);

b = -alpha * max_growth;

for i = 1:numel(num)
    bio = find(strcmp(model.rxns,['BIOMASS_Reaction_',num{i}]));
    ab(bio) = -1; % abundance = 1
    clear bio
end

A_ineq = ab;
b_ineq = b;

% positive semidefinite matrix
H = zeros(size(Aeq,2));
biomass_rxns = find(contains(model.rxns,'BIOMASS_Reaction'));

for i = 1:numel(biomass_rxns)
    H(biomass_rxns(i),biomass_rxns(i)) = 1;
end

params.Threads = 1;

% using Gurobi solver
problem.Q = sparse(H);
problem.A = [Aeq; A_ineq];
problem.rhs = [beq; b_ineq];
problem.obj = f;
problem.lb = lb;
problem.ub = ub;
problem.sense = [repelem('=',size(beq,1),1); repelem('<',size(b_ineq,1),1)];
problem.vtype = repelem('C',size(Aeq,2),1);
problem.modelsense = 'min';

com_solution = gurobi(problem, params);

%com_solution = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub)
end