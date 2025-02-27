% Convert all models to a uniform standard
%% CarveMe
% Tables
tablesDir = '~/IMIC/Program/matlab/COMMIT/data/tables';
topDir = '~/IMIC';

ecTransFile = fullfile(tablesDir, 'corrected-EC-numbers.csv');
coFactorFile = fullfile(tablesDir, 'cofactors_from_KEGG.csv');
formulaeFile = fullfile(tablesDir, 'MNXref', 'MNXref_MET_FORMULAE.csv');
geneIDFile = fullfile(topDir, 'table/abundance_table','rna_ab_and_geneID_20d.csv')
%geneIDFile = fullfile(topDir, 'table/abundance_table','rna_ab_and_geneID_1-1_0h_1.csv');

ID_list = readtable(geneIDFile, 'ReadVariableNames',true);
formulaTab = readtable(formulaeFile, 'ReadVariableNames', true, 'Delimiter', '\t');
ecTranslationTable = readtable(ecTransFile, 'ReadVariableNames', false);
coFactorsTab = readtable(coFactorFile, 'Delimiter', '\t', 'ReadVariableNames', false);
translationDB = loadTranslationDB;

fprintf('################# Converting CarveMe draft models #################\n')
  
CarveMeModelDir = fullfile(topDir, 'models/carveme/draft');
DraftWorkspace = fullfile(topDir, 'models/carveme/draft/carveme_models_draft.mat');
ModelWorkspace = fullfile(topDir, 'models/carveme/carveme_models.mat');
CarveMeModels = dir(fullfile(CarveMeModelDir, '*.xml'));
CarveMeModels = {CarveMeModels.name};
%names = cellfun(@(x)regexp(x, strcat('KG*[0-9]*[^\.]*'), 'match'),CarveMeModels);
names = cellfun(@(x)regexp(x, '\w*.xml', 'match'),CarveMeModels);

n = numel(CarveMeModels);
    
for j=1:n
    fprintf('Reading file #%d/%d: %s\n', j, n, CarveMeModels{j})
    eval(strcat(names{j},'=','readCbModel(fullfile(CarveMeModelDir,CarveMeModels{j}));'))
end
    
% Collect the single models into a cell array 'models'
eval(strcat('models = {', strjoin(names, ';'),'};'))

% make the stoichiometric matridces sparse to be able to save the models
% properly
for j=1:n
    models{j}.S = sparse(models{j}.S);
end
    
save(DraftWorkspace, 'models')

% convert the model to the the standard format
fprintf('\nConverting to standard format and translating to MNXref namespace\n')
for j=1:numel(models)
    fprintf('\nModel %d/%d:\n', num2str(j))
    models{j} = convertCarveMeModel(models{j}, true, translationDB);
    models{j}.EC = correctEC(models{j}.EC, ecTranslationTable);
    % add a formula
    models{j} = addMetFormulae(models{j}, formulaTab);
    % change gene ID
    models{j} = carveme_gene_IDs(models{j},ID_list);
end

save(ModelWorkspace, 'models')
