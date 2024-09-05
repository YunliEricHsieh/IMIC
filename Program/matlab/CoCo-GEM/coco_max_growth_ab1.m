function solution = coco_max_growth_ab1(model, abTable)
% Function for finding the maximum community growth rate with abundance equal 1.
% Input:
%       model:                  community metabolic models as struct object
%       abTable:                table containing the relative abundance of
%                               each organism
% Output:
%       solution:               Solution, returned as a real vector or real 
%                               array.

% find model number
num = extractBetween(abTable.Genome,'KG','_');

% define lower and upper bounds:
ub = model.ub;
lb = model.lb;

% add constraint for export reactions
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
    f(bio) = -1 * 1; % abundance = 1
    % set up minimum growth rate
    lb(bio) = 0.001;
    clear bio
end

% equality constraints:
% LHS
Aeq = model.S;
% RHS
beq = zeros(size(model.S,1),1);

% using Gurobi solver
problem.A = Aeq;
problem.rhs = beq;
problem.obj = f;
problem.lb = lb;
problem.ub = ub;
problem.sense = repelem('=',size(beq,1),1);
problem.vtype = repelem('C',size(Aeq,2),1);

solution = gurobi(problem);

% using linprog
%[solution, max_growth, xflag] = linprog(f, [], [], Aeq, beq, lb, ub);

end