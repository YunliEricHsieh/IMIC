<<<<<<< HEAD
function [min_flux_sum, mets_ID] = MICOM_flux_sum(model, max_growth, micom_obj_value, abTable, alpha, tar_met_IDs)
% Function for the flux sum analysis.
=======
function [flux_sum_value, mets_ID]  = MICOM_flux_sum(model, max_growth, abTable, alpha, tar_met_IDs)
% Function for the second setp of MICOM.
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7
% Input:
%       model:                  community metabolic models as struct object
%       max_growth:             the maximum connumity growth rate
%                               calculated from the first step of MICOM
<<<<<<< HEAD
%       micom_obj_value:        calculated from the second step of MICOM
%                               for quadratic constraint
%       abTable:                table containing the relative abundance of
%                               each organism
%       alpha:                  user-specified trade-off
%       tar_met_IDs:            the list of imported metabolites
% Output:
%       min_flux_sum:           minimum flux sum value of each imported metabolite 
%       mets_ID:                the list of tested metabolite


=======
%       abTable:                table containing the relative abundance of
%                               each organism
%       alpha:                  user-specified trade-off
%       tar_met_IDs:            the list of metabolites
% Output:
%       flux_sum_value:         flux sum value of each targeted metabolite 
%       mets_ID:                the list of tested metabolite

>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7
if ~isempty(tar_met_IDs)
    % use the provided metabolite IDs
    tar_mets = tar_met_IDs;
else
    % find the metabolites in '[e]'
    tar_mets = model.mets(contains(model.mets, '[e]'));
end

% find model ID number
num = extractBetween(abTable.Genome,'KG','_');

<<<<<<< HEAD
% define lower and upper bounds for v and flux sum variable:
=======
% define lower and upper bounds:
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7
ub = [model.ub; ones(numel(tar_mets),1)*10^9];
lb = [model.lb; zeros(numel(tar_mets),1)];

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

<<<<<<< HEAD
% set up minimum growth rate
bio = find(ismember(model.rxns,'BIOMASS_Reaction_'));
=======
% define the objective:
f = [zeros(size(model.S,2),1); zeros(numel(tar_mets),1)];

bio = find(ismember(model.rxns,'BIOMASS_Reaction_'));

% set up minimum growth rate
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7
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
<<<<<<< HEAD
problem.A = [Aeq; A_ineq];
problem.rhs = [beq; b_ineq];
problem.quadcon.Qc = sparse(H);
problem.quadcon.rhs = [micom_obj_value]; % quadratic constraint
problem.quadcon.q = zeros(1, size(Aeq,2));
problem.quadcon.sense = '=';

=======
problem.Q = sparse(H);
problem.A = [Aeq; A_ineq];
problem.rhs = [beq; b_ineq];
problem.obj = f;
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7
problem.lb = lb;
problem.ub = ub;
problem.sense = [repelem('=',size(beq,1),1); repelem('<',size(b_ineq,1),1)];
problem.vtype = repelem('C',size(Aeq,2),1);
<<<<<<< HEAD

mets_ID = tar_mets;
min_flux_sum = [];

% minimum flux sum analysis for each metabolite in [e]
for i = 1:numel(tar_mets)

    % define the objective:
    f = zeros(size(Aeq,2),1);
    num = size(model.S,2) + i;

    % using Gurobi solver
    % minimum flux
    f(num) = 1;
    problem.obj = f;
    solution_min = gurobi(problem);

    if isfield(solution_min, 'x') && ~contains(solution_min.status, 'NUMERIC')
        min_flux_sum = [min_flux_sum; solution_min.objval];
    else 
        % the model couldn't find feasible solution 
        min_flux_sum = [min_flux_sum; nan];
    end  
 
    clear f num solution_min
end
=======
problem.modelsense = 'min';

com_solution = gurobi(problem);

mets_ID = tar_mets;
flux_sum_value = com_solution.x(size(Aeq,2)-numel(tar_mets)+1:end);
>>>>>>> 03da4f62cc99fdbd6ffdf79926f3f674a8b2cfd7

end