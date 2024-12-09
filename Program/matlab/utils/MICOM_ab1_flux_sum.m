function [flux_sum_value, mets_ID] = MICOM_ab1_flux_sum(model, alpha, tar_met_IDs)
% Function for the second setp of MICOM.
% Input:
%       model:                  community metabolic models as struct object
%       abTable:                table containing the relative abundance of
%                               each organism
%       alpha:                  user-specified trade-off
%       tar_met_IDs:            the list of metabolites
% Output:
%       flux_sum_value:         flux sum value of each targeted metabolite 
%       mets_ID:                the list of tested metabolite

if ~isempty(tar_met_IDs)
    % use the provided metabolite IDs
    tar_mets = tar_met_IDs;
else
    % find the metabolites in '[e]'
    tar_mets = model.mets(contains(model.mets, '[e]'));
end

% calculate maximum community growth rate
bio = find(contains(model.rxns, 'BIOMASS_R'));
model.c(bio) = 1;
solution = optimizeCbModel(model);
max_growth = solution.f;

% define lower and upper bounds:
ub = [model.ub; ones(numel(tar_mets),1)*10^9];
lb = [model.lb; zeros(numel(tar_mets),1)];

% define the objective:
f = [zeros(size(model.S,2),1); zeros(numel(tar_mets),1)];

% set up minimum growth rate
lb(bio) = 0.001;

% equality constraints:
% LHS
% flux sum constraint:
tar_mets_num = find(ismember(model.mets, tar_mets));
fs = [-0.5*abs(model.S(tar_mets_num, :)), eye(numel(tar_mets))];

Aeq = [model.S, zeros(size(model.S,1),numel(tar_mets)); fs];

% RHS
beq = [zeros(size(model.S,1),1); zeros(numel(tar_mets),1)];

% constraint of max growth rate
ab = zeros(size(Aeq,2),1)';
b = zeros(1,1);

% -u_c < -alpha * u_max
ab(bio) = -1; % set all the abundance equal 1
b = -alpha * max_growth;

A_ineq = ab;
b_ineq = b;

% positive semidefinite matrix
H = zeros(size(Aeq,2));
biomass_rxns = find(contains(model.rxns,'BIOMASS_Reaction'));

for i = 1:numel(biomass_rxns)
    H(biomass_rxns(i),biomass_rxns(i)) = 1;
end

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

com_solution = gurobi(problem);

mets_ID = tar_mets;
flux_sum_value = com_solution.x(size(Aeq,2)-numel(tar_mets)+1:end);

end