% Convert all models to a uniform standard
%% gapseq
% Tables
tablesDir = '~/IMIC/Program/matlab/COMMIT/data/tables';
topDir = '~/IMIC';

formulaeFile = fullfile(tablesDir, 'MNXref', 'MNXref_MET_FORMULAE.csv');
geneIDFile = fullfile(topDir, 'table/abundance_table','rna_ab_and_geneID_20d.csv');
%geneIDFile = fullfile(topDir, 'table','rna_ab_and_geneID_1-1_0h_1.csv');

ID_list = readtable(geneIDFile, 'ReadVariableNames',true);
formulaTab = readtable(formulaeFile, 'ReadVariableNames', true, 'Delimiter', '\t');
translationDB = loadTranslationDB;

fprintf('################# Converting gapseq draft models #################\n')
  
gapseqModelDir = fullfile(topDir, 'models/gapseq/draft');
DraftWorkspace = fullfile(topDir, 'models/gapseq/draft/gapseq_models_draft.mat');
ModelWorkspace = fullfile(topDir, 'models/gapseq/gapseq_models.mat');
gapseqModels = dir(fullfile(gapseqModelDir, '*.xml'));
gapseqModels = {gapseqModels.name};
%names = cellfun(@(x)regexp(x, strcat('cds_KG*[0-9]*[^\-]*'), 'match'),gapseqModels);
names = cellfun(@(x)regexp(x, strcat('cds_*[\w]*'), 'match'),gapseqModels);

n = numel(gapseqModels);
    
for j=1:n
    fprintf('Reading file #%d/%d: %s\n', j, n, gapseqModels{j})
    eval(strcat(names{j},'=','readCbModel(fullfile(gapseqModelDir,gapseqModels{j}));'))
end

% Collect the single models into a cell array 'models'
eval(strcat('models = {', strjoin(names, ';'),'};'))

% make the stoichiometric matridces sparse to be able to save the models
% properly
for j=1:n
    models{j}.S = sparse(models{j}.S);
    % Change model ID
    %models{j}.modelID =  strcat('KG_',num2str(sscanf(extractAfter(models{j}.modelID,'_'),'KG%d')));
    models{j}.modelID = extractAfter(models{j}.modelID, '_');
end
    
save(DraftWorkspace, 'models')

% convert the model to the the standard format
fprintf('\nConverting to standard format and translating to MNXref namespace\n')
for j=1:numel(models)
    fprintf('\nModel %d/%d:\n', num2str(j))
    models{j} = convertgapseqModel(models{j}, true, translationDB);
    % add a formula
    models{j} = addMetFormulae(models{j}, formulaTab);
    % change gene ID
    models{j} = gapseq_gene_IDs(models{j},ID_list);
end
    
save(ModelWorkspace, 'models');

