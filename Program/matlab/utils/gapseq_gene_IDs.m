function model = gapseq_gene_IDs(model,ID_list)
% Convert gene IDs in gapseq model to a form that 
% is used in transcriptome data.
% Input:
%           struct model:               metabolic model that has been
%                                       reconstructed with gapseq
% Output:
%           struct changed_model:       updated gapseq model

gene = model.genes(contains(model.genes, 'cds'));
gene = erase(gene, 'gp|lcl.');

% separate contig ID and gene/protein IDs
gene = split(gene, '.cds.');

% extract gene/protein IDs
for i = 1:length(gene)
    t = char(gene{i,2});
    inx = strfind(t,'.');

    %gene{i,2} = t(1:inx(2)-1);
    gene{i,2} = t(1:inx(1)-1);
end

Geneid = replace(ID_list.Geneid,'_','.');
Chr = replace(ID_list.Chr,'_','.');
proteinID = replace(ID_list.proteinID,'_','.');

IDs = {};

for i = 1:length(gene)
    % check duplicated protein ID
    % have protein ID and the duplication does not exit
    if sum(contains(proteinID,gene{i,2})) == 1
        IDs{i} = ID_list.Geneid(contains(proteinID,gene{i,2}));

    % does not have protein ID, using gene ID instead
    elseif sum(contains(proteinID,gene{i,2})) == 0
        IDs{i} =  ID_list.Geneid(contains(Geneid,gene{i,2}));

    % the duplication does exit, check with contig ID
    % using cross comparison to find the exact gene name
    else
        warning('There are duplicated protein IDs.')
        IDs{i} = ID_list.Geneid(intersect(find((contains(proteinID,gene{i,2}))), ...
            find(contains(Chr, gene{i,1}))));

    end

    IDs = [IDs{:}]';
end

for i = 1:length(gene)
    model.genes(intersect(find(contains(model.genes,gene{i,2})), ...
        find(contains(model.genes,gene{i,1})))) = IDs(i);
end
