function com_solution = MICOM(model, max_growth, abTable, alpha)
% Function for the second setp of MICOM.
% Input:
%       model:                  community metabolic models as struct object
%       max_growth:             the maximum connumity growth rate
%                               calculated from the first step of MICOM
%       abTable:                table containing the relative abundance of
%                               each organism
%       alpha:                  user-specified trade-off
% Output:
%       solution:               Solution, returned as a real vector or real 
%                               array.


% find model ID number
num = extractBetween(abTable.Genome,'KG','_');

% define lower and upper bounds:
ub = model.ub;
lb = model.lb;

% export reactions
export = model.rxns(contains(model.rxns, 'export_'));

for i = 1:numel(export)
    
    t = char(export{i});
    inx = strfind(t,'_');
    n = t(inx(2)+1:end);

    ex = find(strcmp(model.rxns,export{i}));
    ub(ex) = ub(ex) * abTable.relative_ab(find(contains(abTable.Genome,['KG',n,'_genomic'])));
    
    clear ex t inx n 
end

% define the objective:
f = zeros(size(model.S,2),1);

bio = find(ismember(model.rxns,'BIOMASS_Reaction_'));

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

b = -alpha * max_growth;

for i = 1:numel(num)
    bio = find(strcmp(model.rxns,['BIOMASS_Reaction_',num{i}]));
    ab(bio) = -abTable.relative_ab(find(contains(abTable.Genome,['KG',num{i},'_genomic'])));
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