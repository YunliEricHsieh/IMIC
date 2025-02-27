disp('----------------------------------------------------------------------')
disp('                       Iterative gap filling                          ')
disp('----------------------------------------------------------------------')

% load options script 
disp('Loading required data...')

options

% medium
load(mediumFile, 'medium')
%load(LB_mediumFile, 'LB_medium')

% Universal database
load(dbFile)

methods = {'carveme','gapseq','kbase','consensus'};
%methods = {'consensus'};
modelDir =  '~/IMIC/models';
%modelDir =  '~/IMIC/study_case/models';

% start parallel pool
delete(gcp('nocreate'))
parpool(ncpu);

for j = 1:numel(methods)

    disp('----------------------------------------------------------------------')
    fprintf('%s models for OTU composition in %s dataset\n', methods{j})
    
    % load model cell array
    model_workspace = fullfile(modelDir, [methods{j}, '_draft_models_biomass.mat']);
    
    load(model_workspace)
    
    if ~exist('merged_models', 'var')
        merged_models = models; clear models
    end
     
    % load OTU abundances to take the subset found in the current study
    if exist('otuDir', 'var') && ~isempty(otuDir)
        otu_file = fullfile(otuDir, habitat{j}, 'otutab.txt');
        tab_merged = readAbundancesFromFile(otu_file);
        order = tab_merged;
        
        model_ids = cellfun(@(x)strtok(x.id), merged_models,...
            'UniformOutput', false);
        merged_models = merged_models(ismember(model_ids, tab_merged.Properties.RowNames));
        clear tab_merged
    elseif exist('otuDir', 'var') && ~isempty(otuDir) && ~exist('otuDir', 'dir')
        error('Given path to OTU abundance directory does not exist: %s', otuDir)
    end
    
    % create cell arrays for individual auxotrophic media
    % get the list of complementary media
    if exist('mediaDir', 'var') && ~isempty(mediaDir) && exist('mediaDir', 'dir')
        auxo_media = dir(fullfile(mediaDir, '*.tsv'));
        auxo_media = fullfile({auxo_media.folder}, {auxo_media.name})';
        auxo_media = auxo_media(contains(auxo_media, strcat(model_ids, '.tsv')));
        
        for k=1:numel(auxo_media)
            tmp_tab = importdata(auxo_media{k});
            tmp_tab = tmp_tab.textdata(2:end,1);
            auxo_media{k} = strcat(translateIDs(tmp_tab, 'met', [], 'ModelSEED',...
                'MNXref', false), '[e]');
        end
    elseif exist('mediaDir', 'var') && ~isempty(mediaDir) && ~exist('mediaDir', 'dir')
        error('Given path to auxotrophic media directory does not exist: %s', mediaDir)
    else
        auxo_media = repelem({medium}, numel(merged_models), 1);
        %auxo_media = repelem({LB_medium}, numel(merged_models), 1);
    end
    
    % write single .mat files for every model in a subdirectory
    tmpModelDir = fullfile(modelDir, 'tmp_models', [methods{j}, tmp_spec]);
    [s, fileList] = createModelDir(merged_models, tmpModelDir);
    
    clear merged_models
    
    if s
        
        fileList = fullfile(tmpModelDir, fileList);
        
        % run iterative gap filling
        [GF, EX, gf_order, solutions, exc, gf, bio] = ...
            iterativeGapFilling(fileList, medium, auxo_media,...
            dbModel_MNXref_balanced, weights, epsilon, include_sink,...
            order, iterations, seq_sim_workspace);
        % [GF, EX, gf_order, solutions, exc, gf, bio] = ...
        %     iterativeGapFilling(fileList, LB_medium, auxo_media,...
        %     dbModel_MNXref_balanced, weights, epsilon, include_sink,...
        %     order, iterations, seq_sim_workspace);

        outDir = fullfile('~/IMIC/models/iterative');
        % outDir = fullfile('~/IMIC/study_case/models/iterative');

        save(fullfile(outDir, methods{j}),...
            'GF', 'EX', 'gf_order', 'solutions', 'exc', 'gf', 'bio')
        
        % write results to file
        opt = find(all(solutions==gf_order, 2));
        V = {exc, gf, bio, opt(1)};
        vars = {'exc', 'gf', 'bio', 'opt'};
        
        for l=1:numel(vars)
            writetable(array2table(V{l}, 'VariableNames', sprintfc('M_%d', 1:size(V{l}, 2))),...
                fullfile(outDir, [methods{j}, '-', vars{l}, '.txt']),...
                'WriteVariableNames', true, 'WriteRowNames', false,...
                'Delimiter', '\t')
        end
        
        EX = translateIDs(strtok(EX, '[')', 'met', [], 'MNXref', 'NAMES');
        writetable(cell2table(EX, 'VariableNames', {'exported'}),...
           [outDir, '/Exchanged_mets_', methods{j},'.txt'])
        
        clear GF EX gf_order solutions exc gf bio
        
        rmdir(tmpModelDir, 's')
        
    end
end
