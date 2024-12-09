function com_solution = MICOM_ab1(model, alpha)
% Function for the second setp of MICOM.
% Input:
%       model:                  community metabolic models as struct object
%       abTable:                table containing the relative abundance of
%                               each organism
%       alpha:                  user-specified trade-off
% Output:
%       solution:               Solution, returned as a real vector or real 
%                               array.

% calculate maximum community growth rate
bio = find(contains(model.rxns, 'BIOMASS_R'));
model.c(bio) = 1;
solution = optimizeCbModel(model);
max_growth = solution.f;

% define lower and upper bounds:
ub = model.ub;
lb = model.lb;

% define the objective:
f = zeros(size(model.S,2),1);

% set up minimum growth rate
lb(bio) = 0.001;

% equality constraints:
% LHS
Aeq = model.S;
% RHS
beq = zeros(size(model.S,1),1);

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

%com_solution = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub)
end