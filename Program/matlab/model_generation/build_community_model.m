modelDir = '~/IMIC/models/iterative';
modelOut = '~/IMIC/models';
methods = {'consensus'};

% medium that has been used for gap filling
% case 1
%mediumFile = 'media/M9-medium-anaerobic.mat';
% case 2
LB_mediumFile = '~/IMIC/media/LB-medium.mat';

for k = 1:numel(methods)
    methods{k}
    model_workspace = fullfile(modelDir, [methods{k}, '.mat']);
    load(model_workspace)

%    load(mediumFile, 'medium')
    load(LB_mediumFile, 'LB_medium')

    all_medium = LB_medium;
%    all_medium = medium;
    all_medium = extractBefore(all_medium,'[e]');

    % set up an initial model. We choose the first gap-filled model
    com_model = GF{gf_order(1)};

    % find model ID number
    % case 1
%    num = extractAfter(com_model.id,'KG_');
    % case 2
    if contains(com_model.id, 'Ecoli')
        num = '1';
    else
        num = '2';
    end

    % assign reaction IDs with compartment number
    com_model.rxns = replace(com_model.rxns,'_c',['_c',num]);
    com_model.rxns = replace(com_model.rxns,'_e',['_e',num]);
    com_model.rxns = replace(com_model.rxns,'BIOMASS_Reaction',['BIOMASS_Reaction_',num]);
    com_model.rxns = replace(com_model.rxns,'sink_BIOMASS[c]',['sink_BIOMASS[c',num,']']);
    com_model.rxns(contains(com_model.rxns, 'export_')) = ...
        append(com_model.rxns(contains(com_model.rxns, 'export_')),'_', num);

    % assign metabolites IDs with compartment number except for the metabolites
    % in extracellular space
    com_model.mets = replace(com_model.mets,'[c',['[c_',num]);

    com_model.description = 'community_model';
    com_model.id = 'community_model';

    clear num
    
    % remove the initial model from gf_order list
    gf_order(1) = [];

    % merge models into a big model (distinguish different organism with different compartment)
    for i = gf_order
    
        model = GF{i};
        disp(['Model id: ', model.id])

        % find model ID number
        % case 1
%        num = extractAfter(model.id,'KG_');
        % case 2
        if contains(model.id, 'Ecoli')
            num = '1';
        else
            num = '2';
        end

        % assign reaction IDs with compartment number
        model.rxns = replace(model.rxns,'_c',['_c',num]);
        model.rxns = replace(model.rxns,'_e',['_e',num]);
        model.rxns = replace(model.rxns,'BIOMASS_Reaction',['BIOMASS_Reaction_',num]);
        model.rxns = replace(model.rxns,'sink_BIOMASS[c]',['sink_BIOMASS[c',num,']']);
        model.rxns(contains(model.rxns, 'export_')) = ...
            append(model.rxns(contains(model.rxns, 'export_')),'_', num);

        % assign metabolites IDs with compartment number
        model.mets = replace(model.mets,'[c',['[c_',num]);

        % add the reactions into community model
        for j = 1:numel(model.rxns)

            % reaction formula
            rxn_formula = char(printRxnFormula(model,model.rxns(j)));
            
            % add reaction into community model
            com_model = addReaction(com_model,model.rxns{j},'reactionFormula',rxn_formula,... 
                                'geneRule', model.grRules{j});

            % add EC number into the community model
            if ~ismissing(model.EC(j))

                % find the reaction ID
                rxnID = model.rxns{j};

                % find the index of reaction in the community model
                index = find(ismember(com_model.rxns, rxnID));

                % add EC number into the community model
                com_model.EC(index) = model.EC(j);
            end

        end
    
        com_model.comps(i) = {['c_',num]};
        com_model.compNames(i) = {['cytosol_',num]};

        clear num model
    end

    % add 'extracellular space' into compNames
    com_model.comps(numel(GF)+1) = {'e'};
    com_model.compNames(numel(GF)+1) = {'extracellular space'};

    % Now, all the models in the community have same extracellular space.
    % Therefore, we can remove EX and sink reactions, which were used to
    % connect exchanged metabolites between models in the community in the COMMIT.
    exported = com_model.rxns(contains(com_model.rxns, 'export_'));
    exported = extractBetween(exported, '_','_');

    % keep the metabolites provided from medium
    exported_met = setdiff(exported, all_medium);

    EX = com_model.rxns(find(contains(com_model.rxns,'EX_')));
    sink = com_model.rxns(find(contains(com_model.rxns,'sink_')));

    for i = 1:numel(exported_met)
        if sum(contains(sink,[exported_met{i},'[e]'])) > 0
            sink_ID = sink(contains(sink,[exported_met{i},'[e]']));
            com_model = removeRxns(com_model,sink_ID);
        elseif sum(contains(EX,[exported_met{i},'[e]'])) > 0
            EX_ID = sink(contains(EX,[exported_met{i},'[e]']));
            com_model = removeRxns(com_model,EX_ID);
        end
    end

    save(fullfile(modelOut, [methods{k},'_com.mat']), 'com_model')
end