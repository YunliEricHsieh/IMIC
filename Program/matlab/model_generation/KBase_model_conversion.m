% Convert all models to a uniform standard
%% KBase
% Tables
tablesDir = '~/IMIC/Program/matlab/COMMIT/data/tables';
topDir = '~/IMIC';

ModelSEEDRxnEC = fullfile(tablesDir, 'ModelSEED_RXN_EC.csv');
formulaeFile = fullfile(tablesDir, 'MNXref', 'MNXref_MET_FORMULAE.csv');
ecTransFile = fullfile(tablesDir, 'corrected-EC-numbers.csv');

formulaTab = readtable(formulaeFile, 'ReadVariableNames', true, 'Delimiter', '\t');
ecRxnTable = readtable(ModelSEEDRxnEC, 'ReadVariableNames', false);
ecTranslationTable = readtable(ecTransFile, 'ReadVariableNames', false);
translationDB = loadTranslationDB;

fprintf('################# Converting KBase draft models #################\n')
  
KBaseModelDir = fullfile(topDir, 'models/kbase/draft');
DraftWorkspace = fullfile(topDir, 'models/kbase/draft/kbase_models_draft.mat');
ModelWorkspace = fullfile(topDir, 'models/kbase/kbase_models.mat');
KBaseModels = dir(fullfile(KBaseModelDir, '*.xml'));
KBaseModels = {KBaseModels.name};
%names = cellfun(@(x)regexp(x, strcat('KG*[0-9]*[^\.]*'), 'match'),KBaseModels);
names = cellfun(@(x)regexp(x, '\w*.xml', 'match'),KBaseModels);

n = numel(KBaseModels);
    
for j=1:n
    fprintf('Reading file #%d/%d: %s\n', j, n, KBaseModels{j})
    eval(strcat(names{j},'=','readCbModel(fullfile(KBaseModelDir,KBaseModels{j}));'))
end
    
% Collect the single models into a cell array 'models'
eval(strcat('models = {', strjoin(names, ';'),'};'))

% make the stoichiometric matridces sparse to be able to save the models
% properly
for j=1:n
    models{j}.S = sparse(models{j}.S);
    % Add E.C. numbers
    % Retrieve ec numbers from ModelSEED reaction database
    ec = ecNumbersFromModelSEED(models{j}.rxns, ecRxnTable);
    % Correct the EC numbers as they could be outdated
    models{j}.EC =  correctEC(ec, ecTranslationTable);
    % Change gene names
    models{j}.genes =  erase(models{j}.genes, 'gene-');
    models{j}.genes =  erase(models{j}.genes, '.CDS');
    % Change model ID
    %models{j}.modelID =  strcat('KG_',num2str(sscanf(extractBefore(models{j}.description,'_'),'KG%d')));
    models{j}.modelID = extractBefore(models{j}.description,'_');
end
    
save(DraftWorkspace, 'models')
 
% convert the model to the the standard format
fprintf('\nConverting to standard format and translating to MNXref namespace\n')
for j=1:numel(models)
    fprintf('\nModel %d/%d:\n', num2str(j))
    models{j} = convertKBaseModel(models{j}, true, translationDB);
    % add a formula
    models{j} = addMetFormulae(models{j}, formulaTab);
end
    
save(ModelWorkspace, 'models');

