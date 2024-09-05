modelDir =  '~/IMIC/models';

% load models
model_workspace = fullfile(modelDir, 'consensus_com.mat');
load(model_workspace)

% calculate maximum community growth rate
bio = find(contains(com_model.rxns, 'BIOMASS_R'));
com_model.c(bio) = 1;
solution = optimizeCbModel(com_model);
max_com_growth = solution.f;

% define lower and upper bounds:
ub = com_model.ub;
lb = com_model.lb;

% define the objective:
f = zeros(size(com_model.S,2),1);

% equality constraints:
% LHS
Aeq = com_model.S;
% RHS
beq = zeros(size(com_model.S,1),1);

% constraint of max growth rate
ab = zeros(size(Aeq,2),1)';
b = zeros(1,1);

% -u_c < -alpha * u_max
ab(bio) = -1; % set all the abundance equal 1
b = -1 * max_com_growth;

A_ineq = ab;
b_ineq = b;

% positive semidefinite matrix
H = zeros(size(Aeq,2));
biomass_rxns = find(contains(com_model.rxns,'BIOMASS_Reaction'));

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
