 function results_table = var_analysis(model, transcriptTable, lambda, maximum_value)
% Function for doing the variability analysis.
% Input:
%       model:                  community metabolic models as struct object
%       transcriptTable:        table containing the abundance of
%                               transcriptomic data from each organism
%       lambda:                 balancing factor
%       maximum_value:          maximum community growth rate calculated
%                               from IMIC
% Output:
%       results_tabel:          max/min growth rate value for each MAGs

% find model ID number
MAG_ID = model.rxns(contains(model.rxns,'BIOMASS_R'));

% find M : the max value of f_i(g)
express_value = [];
count = 0;

for i = 1:numel(model.rxns)
    if ~isempty(model.rules{i}) && ~contains(model.grRules{i},'spontaneous')
        % split the genes by 'or'
        gpr_split = strsplit(model.rules{i},'|');
        gpr_split_TPM = [];

        for j = numel(gpr_split)
            list = extractBetween(gpr_split{j},'x(',')');
            list_TPM = [];

            for k = 1:numel(list)
                gene_id = model.genes(str2num(list{k}));
                list_TPM(k) = transcriptTable.TPM(strcmp(transcriptTable.Geneid,gene_id));
            end
            
            % find min value
            gpr_split_TPM(j) = min(list_TPM);
        end

        % sum the value
        rxn_transcript_ab = sum(gpr_split_TPM);

        count = count + 1;
        express_value(count) =  rxn_transcript_ab;
    end
end

% define lower and upper bounds for v, B+ and B-:
ub = [model.ub; ones(size(express_value,2)*2,1)*1000];
lb = [model.lb; zeros(size(express_value,2)*2,1)];

% equality constraints:
% LHS
a = [zeros(size(model.S,2),1); -ones(size(express_value,2)*2,1)*lambda]';

bio = find(contains(model.rxns,'BIOMASS_Reaction_'));
a(bio) = 1;

Aeq = [model.S, zeros(size(model.S,1),size(express_value,2)*2); a];

% RHS
b = maximum_value;
beq = [zeros(size(model.S,1),1);b];

% inequality constraints:
n_col = size(Aeq,2);
n_row = size(express_value,2)*3;

A_ineq = zeros(n_row, n_col);
b_ineq = zeros(size(express_value,2)*3,1);

row_count = 0;

for i = 1:numel(model.rxns)
    
    % reactions have GPR
    if ~isempty(model.rules{i}) && ~contains(model.grRules{i},'spontaneous')

        % inequality constraints
        row_count = row_count + 1;

         % V - Vm*B+ + Vm*B- < Vm*f/M
        A_ineq(row_count, [i, size(model.S,2)+row_count, size(model.S,2)+size(express_value,2)+row_count]) = [1, -1000, 1000];
        b_ineq(row_count) = ub(i) * express_value(row_count) / max(express_value);

        % B+ - B- < 1-f/M
        A_ineq(size(express_value,2)+row_count, [size(model.S,2)+row_count, size(model.S,2)+size(express_value,2)+row_count]) = [1, -1];
        b_ineq(size(express_value,2)+row_count) = 1 - (express_value(row_count) / max(express_value));

        % -B+ + B- < f/M
        A_ineq(size(express_value,2)*2+row_count, [size(model.S,2)+row_count, size(model.S,2)+size(express_value,2)+row_count]) = [-1, 1];
        b_ineq(size(express_value,2)*2+row_count) = (express_value(row_count) / max(express_value));
    end
end

problem.A = [Aeq; A_ineq];
problem.rhs = [beq; b_ineq];
problem.lb = lb;
problem.ub = ub;
problem.sense = [repelem('=',size(beq,1),1); repelem('<',size(b_ineq,1),1)];
problem.vtype = repelem('C',size(Aeq,2),1);

results_table = {};
for i = 1:numel(MAG_ID)
    num = extractAfter(MAG_ID{i},'BIOMASS_Reaction_');
    dis = ['Calculating the maximum growth rate for model : KG_', num];
    disp(dis)

    % define the objective:
    f = [zeros(size(model.S,2),1); zeros(size(express_value,2)*2,1)];
    bio = find(ismember(model.rxns,['BIOMASS_Reaction_',num]));

    % using Gurobi solver
    % maximum growth rate
    f(bio) = -1;
    problem.obj = f;
    solution_max = gurobi(problem);

    if contains(solution_max.status, 'INFEASIBLE')
        solution_max.objval = 'NA';
    else
        solution_max.objval = -solution_max.objval;
    end

    % minimum growth rate
    f(bio) = 1;
    problem.obj = f;
    solution_min = gurobi(problem);

    if contains(solution_min.status, 'INFEASIBLE')
        solution_min.objval = 'NA';
    end

    results_t = {['KG_',num], solution_max.objval, solution_min.objval};
    results_table = [results_table; results_t];

    clear f bio num solution_max solution_min results_t
end
end