function [flux_sum_value, mets_ID] = IMIC_flux_sum(model, transcriptTable, lambda, tar_met_IDs)
% Function for doing IMIC and flux sum analysis.
% Input:
%       model:                  community metabolic models as struct object
%       transcriptTable:        table containing the abundance of
%                               transcriptomic data from each organism
%       lambda:                 balancing factor
%       tar_met_IDs:            the list of metabolites
% Output:
%       flux_sum_value:         flux sum value of each targeted metabolite 
%       mets_ID:                the list of tested metabolite

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

if ~isempty(tar_met_IDs)
    % use the provided metabolite IDs
    tar_mets = tar_met_IDs;
else
    % find the metabolites in '[e]'
    tar_mets = model.mets(contains(model.mets, '[e]'));
end

% define lower and upper bounds for v, B+, B-, and flux sum variable:
ub = [model.ub; ones(size(express_value,2)*2,1)*1000; ones(numel(tar_mets),1)*10^9];
lb = [model.lb; zeros(size(express_value,2)*2,1); zeros(numel(tar_mets),1)];

% equality constraints:
% LHS
% flux sum constraint:
tar_mets_num = find(ismember(model.mets, tar_mets));
fs = [-0.5*abs(model.S(tar_mets_num, :)), zeros(numel(tar_mets), size(express_value,2)*2), eye(numel(tar_mets))];

Aeq = [model.S, zeros(size(model.S,1),size(express_value,2)*2+numel(tar_mets)); fs];

% RHS
beq = [zeros(size(model.S,1),1); zeros(numel(tar_mets),1)];

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

% define the objective:
f = [zeros(size(model.S,2),1); ones(size(express_value,2)*2,1)*lambda; zeros(numel(tar_mets),1)];

bio = find(contains(model.rxns,'BIOMASS_Reaction'));
f(bio) = -1;

% using Gurobi solver
problem.A = [Aeq; A_ineq];
problem.rhs = [beq; b_ineq];
problem.obj = f;
problem.lb = lb;
problem.ub = ub;
problem.sense = [repelem('=',size(beq,1),1); repelem('<',size(b_ineq,1),1)];
problem.vtype = repelem('C',size(Aeq,2),1);

solution = gurobi(problem);

mets_ID = tar_mets;
flux_sum_value = solution.x(size(Aeq,2)-numel(tar_mets)+1:end);

end
