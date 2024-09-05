function changed_model = convertgapseqModel(model, translate, translationDB)
% Convert a metabolic model reconstructed using the KBase web service to a form
% that can be used with the FastGapFilling function. For this, field names
% will be changed, added or deleted.
% Input:
%           struct model:               metabolic model that has been
%                                       reconstructed with KBase
%           logical translate:          if true, model ids will be translated
%                                       to MNXref ids
% Outputd:
%           struct changed_model:       updated KBase model

% Specify fields to delete
delete = {'proteins', 'osenseStr', 'geneNames', 'metFormulas', 'modelVersion',...
    'modelNotes'};

% Specify field names that have to be changed
original_names = {'modelID','rxnECNumbers'};
new_names = {'id','EC'};

% Check if all fields that should be changed or deleted are present
check_fields = all(cellfun(@(x)isfield(model, x), [delete, original_names]));
if ~check_fields
    error('The input model does not have all fields that are usually contained in a gapseq reconstruction')
end

% Change model fields
for i=1:numel(original_names)
   tmp = model.(original_names{i});
   model.(new_names{i}) = tmp;
   model = rmfield(model, original_names{i});
end
clear tmp

% Remove '0' from the model compartments
model.comps = regexprep(model.comps, '0', '');
model.compNames = regexprep(model.compNames, '0', '');

% Remove '0' from the reaction ids, reaction names and metabolite ids
model.rxns = regexprep(model.rxns, '0$', '');

% change the reaction ids of export reactions by translating the exported
% metabolite ids to the MNXref namespace
ex_rxns = cellfun(@(x)contains(x, 'EX_'), model.rxns);
tmp_num = find(cellfun(@(x)contains(x, 'EX_'), model.rxns));
tmp_suffix = cellfun(@(x)regexp(x, '_.$', 'match'), model.rxns(ex_rxns));
tmp_mets = regexp(model.rxns(ex_rxns), 'cpd[0-9]*', 'match');
tmp_mets = translateIDs(tmp_mets, 'met', translationDB.metTab, 'ModelSEED', 'MNXref');

for i = 1:numel(tmp_mets)
    if ~isempty(tmp_mets{i})
        model.rxns(tmp_num(i)) = cellstr(strcat('EX_', tmp_mets{i}, tmp_suffix{i}));
    end
end
clear ex_rxns tmp_suffix tmp_mets tmp_num

% change the reaction ids of sink reactions by translating the exported
% metabolite ids to the MNXref namespace
sink_rxns = cellfun(@(x)contains(x, 'DM_'), model.rxns);
tmp_num = find(cellfun(@(x)contains(x, 'DM_'), model.rxns));
tmp_suffix = cellfun(@(x)regexp(x, '_.$', 'match'), model.rxns(sink_rxns));
tmp_mets = regexp(model.rxns(sink_rxns), 'cpd[0-9]*', 'match');
tmp_mets = translateIDs(tmp_mets, 'met', translationDB.metTab, 'ModelSEED', 'MNXref');

for i = 1:numel(tmp_mets)
    if ~isempty(tmp_mets{i})
        model.rxns(tmp_num(i)) = cellstr(strcat('EX_', tmp_mets{i}, tmp_suffix{i}));
    end
end
clear sink_rxns tmp_suffix tmp_mets

model.mets = regexprep(model.mets, '0\]', '\]');

% Remove the periplasmal compartment
model.comps = {'c'; 'e'};
model.compNames = {'cytosol'; 'extracellular space'};

% Find if there are reactions belonging to periplasmal compartment or not
if find(contains(model.rxns, '_p'))
    error('The model contains the reactions belonging to periplasmal compartment.')
end

% move all metabolites remaining in the periplasm to the extra-cellular
% compartment (change in every reaction and delete the old ids)
old_met_idx = find(contains(model.mets, '[p]'));
old_mets = model.mets(contains(model.mets, '[p]'));
new_mets = regexprep(old_mets, '\[p\]$', '\[e\]');
for i=1:numel(old_mets)
    rxns_met = findRxnsFromMets(model, old_mets{i});
    for j=1:numel(rxns_met)
        if ~ismember(new_mets{i}, model.mets)
            name = model.metNames{old_met_idx(i)};
            formula = model.metFormulas(old_met_idx(i));
            %note = model.metNotes(old_met_idx(i));
            model = addMetabolite(model, new_mets{i}, 'metName', name);
            model.metFormulas(end) = formula;
            %model.metNotes(end) = note;
        end
        model = changeRxnMets(model, old_mets{i}, new_mets{i}, rxns_met{j});
    end
end
model = removeMetabolites(model, old_mets);

% add 'csense' field to add the type of constraint for each metabolite
model.csense = repmat('E', size(model.S, 1), 1);

% Change the id of the biomass reaction
model.c(find(model.c == 1)) = 0;
model.c(find(contains(model.rxns,'bio'))) = 1;
model.rxns{logical(model.c)} = 'BIOMASS_Reaction';

% change '_' to '.' in the gene names
model.genes = cellfun(@(x)regexprep(x, '_', '\.'), model.genes,...
    'UniformOutput', false);
model.genes = cellfun(@(x)regexprep(x, '\.', '\|', 'once'), model.genes,...
            'UniformOutput', false);

% Remove fields from the gapseq model
model = rmfield(model, delete);

% re-order the field names by size:
field_size = cellfun(@(x)numel(model.(x)), fieldnames(model));
[~,field_order]= sort(field_size);
idx_first = find(ismember(fieldnames(model), {'id', 'description'}));
field_order(ismember(field_order, idx_first)) = [];
field_order = vertcat(idx_first, field_order);
model = orderfields(model, field_order);

if translate
    model = translateModel(model, 'ModelSEED', 'MNXref', translationDB);
end

% unify duplicate metabolite IDs
model = removeDuplicateMets(model);

% unify duplicate reaction IDs
tmp = load(fullfile('data', 'gap-filling', 'database',...
    'Universal-model-MNXref-balanced.mat'));
dbModel = tmp.dbModel_MNXref_balanced;
model = removeDuplicateRxns(model, dbModel);

changed_model = model;

end
