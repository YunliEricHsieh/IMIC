modelDir = '~/IMIC/study_case/models/iterative';
modelOut = '~/IMIC/study_case/models';
% methods = {'carveme','gapseq','kbase','consensus'}; % case 1
methods = {'consensus'}; % case 2

% medium that has been used for gap filling
% case 1
% mediumFile = 'media/M9-medium-anaerobic.mat';
% case 2
LB_mediumFile = '~/IMIC/media/LB-medium.mat';

for k = 1:numel(methods)
    methods{k}
    model_workspace = fullfile(modelDir, [methods{k}, '.mat']);
    load(model_workspace)

    % load(mediumFile, 'medium') % case 1
    load(LB_mediumFile, 'LB_medium') % case 2

    all_medium = LB_medium; % case 2
    % all_medium = medium; % case 1
    all_medium = extractBefore(all_medium,'[e]');

    % set up an initial model.
    % com_model = GF{gf_order(1)}; % in case 1, we choose the first gap-filled model

    % in case 2, we choose the first P. putida model
    if ismember(GF{1}.id, 'Pputida')
        com_model = GF{1};
        gf_order(1) = [];
    else
        com_model = GF{2};
        gf_order(2) = [];
    end

    % find model ID number
    % case 1
    % num = extractAfter(com_model.id,'KG_');
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
    % gf_order(1) = []; % case 1

    % merge models into a big model (distinguish different organism with different compartment)
    for i = gf_order
    
        model = GF{i};
        disp(['Model id: ', model.id])

        % find model ID number
        % case 1
        % num = extractAfter(model.id,'KG_');
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
            rxn_formula = char(printRxnFormula(model,model.rxns(j), false));


            if sum(ismember(com_model.rxns, model.rxns{j})) == 0
                % add reaction into community model
                disp(['Add reaction: ', model.rxns{j}])

                com_model = addReaction(com_model,model.rxns{j},'reactionFormula',rxn_formula);

                % add EC number into the community model
                % find the reaction ID
                rxnID = model.rxns{j};

                % find the index of reaction in the community model
                index = find(ismember(com_model.rxns, rxnID));

                % add EC number into the community model
                com_model.EC(index) = model.EC(j);

                % add GPR rules into the community model
                com_model.grRules(index) = model.grRules(j);

            end
        end

        % add genes into the community model
        com_model.genes = unique([com_model.genes; model.genes]);
       
        com_model.comps(i) = {['c_',num]};
        com_model.compNames(i) = {['cytosol_',num]};

        clear num model
    end

    % remove rules field
    com_model = rmfield(com_model, 'rules');

    % create rules field
    com_model = generateRules(com_model);

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
