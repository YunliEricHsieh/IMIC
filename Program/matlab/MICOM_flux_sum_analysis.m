%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["EMPTY", "EmptyEquation"];
opts.VariableTypes = ["string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["EMPTY", "EmptyEquation"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["EMPTY", "EmptyEquation"], "EmptyFieldRule", "auto");

% Import the data
rxn_description = readtable('~/IMIC/Program/matlab/COMMIT/data/tables/rxn_desciptions.csv', opts);

clear opts

%% FVA analysis for each transport reaction
modelDir =  '~/IMIC/models';
tablesDir = '~/IMIC/table';

timepoint = {'20d', '40d', '60d', '90d', '180d'};

% read the community model
model_workspace = fullfile(modelDir, 'consensus_com.mat');
load(model_workspace);

alpha = 0.3;

% find transport reactions according to the rxn description
keywords = {'exchange', 'transport', 'import', 'export'};

rxn_ex = contains(rxn_description.EmptyEquation,keywords{1});
rxn_tr = contains(rxn_description.EmptyEquation,keywords{2});
rxn_im = contains(rxn_description.EmptyEquation,keywords{3});
rxn_exp = contains(rxn_description.EmptyEquation,keywords{4});

any = rxn_ex | rxn_tr | rxn_im | rxn_exp;
all_transport_rxn = rxn_description(any,:);

tmp_transport = {};
import_mets = {};

% find the import reactions existing in the model
for j = 1:numel(com_model.rxns)

    if sum(ismember(all_transport_rxn.EMPTY, extractBefore(com_model.rxns{j},'_')))>0

        % find the common transport reactions over all MAGs
        ID = all_transport_rxn.EMPTY(ismember(all_transport_rxn.EMPTY, extractBefore(com_model.rxns{j},'_')));
        rxn_ID = com_model.rxns(contains(extractBefore(com_model.rxns,'_'), ID));

        % remove the backward reactions
        rxn_ID = rxn_ID(~contains(rxn_ID, '_r'));

        if numel(rxn_ID) == 14
            
            % choose the reactions which transport the metabolite from [e]
            % to [c]
            reactionFormula = printRxnFormula(com_model, com_model.rxns{j}, false); 
            parts = strsplit(reactionFormula{1}, '->');
            lhs = strtrim(parts{1});
    
            % Check for '[e]' in left hand side
            if contains(lhs, '[e]') 
                tmp_transport = [tmp_transport; com_model.rxns{j}];
                tmp_mets = regexp(lhs, '\w+\[e]', 'match');
                import_mets = [import_mets, tmp_mets];
            end
        end
    end
end
    
% find the unique metabolites IDs
import_mets = unique(import_mets)';

% remove the proton from the list
import_mets = import_mets(~cellfun(@(x) isequal(x, 'MNXM1[e]'), import_mets));

flux_sum_min = [];
metabolite_ID = {};

ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);

for i = 1:numel(timepoint)

    disp('----------------------------------------------------------------------')
    fprintf('\n################# %s\n\n', timepoint{i})
    
    abFile = fullfile(tablesDir, 'abundance_table', ['relative_ab_',timepoint{i}, '.csv']); 
    abTable = readtable(abFile, 'ReadVariableNames', true);

    disp('Calculate Max Growth Rate')
    micom_solution = MICOM_max_growth(com_model, abTable);
    max_growth = -micom_solution.objval;

    disp('Calculate Community Growth Rate with MICOM')
    com_solution = MICOM(com_model, max_growth, abTable, alpha);

    micom_obj_value = com_solution.objval;

    % FVA for transport reactions
    [min_flux_sum, mets_ID]  = MICOM_min_flux_sum(com_model, max_growth, micom_obj_value, abTable, alpha, import_mets);

    flux_sum_min(:, i) = min_flux_sum;
    metabolite_ID{i} = mets_ID;

end

imp_metabolite = vertcat(metabolite_ID{1});
results = table(imp_metabolite,flux_sum_min);
writetable(results, fullfile(tablesDir,'flux_sum_analysis','MICOM_flux_sum.csv'));

