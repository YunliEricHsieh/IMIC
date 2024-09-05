% Options for setting

%~~~~~~~~~~~~~ General settings ~~~~~~~~~~~~~%

% number of parallel workers
ncpu = 20;

% experiments to be analyzed
methods = {'carveme','kbase','gapseq','consensus'};

% path to top directory
topDir = '~/IMIC';

% medium file
mediumFile = 'media/M9-medium-anaerobic.mat';
LB_mediumFile = 'media/LB-medium-anaerobic.mat';

% directory for auxotrophy media (files must be of form 'model id'.tsv)
% optional: Remove the following line or set mediaDir='' if you do not want
% to use auxotrophic media to fill gaps in the first model of each
% iteration. The minimal medium given above will be used instead.
% directory for auxotrophy media (files must be of form 'model id'.tsv)
mediaDir = '';

% top directory for OTU abundance (subdirectories must be habitat >
% experiment > otutab.txt
% optional: remove the following line or set otuDir='' if you do not want
% to use OTU abundances to select at subset of models found the experiment
% and just use all models
otuDir = '';

% threshold for the biomass reaction
epsilon = 1E-3;

% number of iterations for the iterative approach
iterations = 100;

% print level
verbose = true;

% provide an order for iterative gap filling (optional)
order = [];

% whether or not sink reactions should be included in the gap-filling
% objective
include_sink = false;

% change to working directory
cd(topDir)

% add all COMMIT matlab scripts to path
addpath(genpath('~/IMIC/Program/matlab/COMMIT/code/matlab'))

%~~~~~~~~~~~~~ Model workspace and output directory ~~~~~~~~~~~~~%

tmp_spec = '_all';

%~~~~~~~~~~~~~ Gap-filling resources ~~~~~~~~~~~~~%

% file for the database
dbFile = '~/IMIC/Program/matlab/COMMIT/data/gap-filling/database/Universal-model-MNXref-balanced.mat';

% path to sequence similarity workspace
seq_sim_workspace = '';


%~~~~~~~~~~~~~ Reaction-type-specific weights ~~~~~~~~~~~~~%

% reactions with sequence evidence (E <= 1E-6) 
weights.sequence = 25;
% transport reactions
weights.transport = 100;
% transport reactions operating on highly-permeable metabolites
weights.permeable = 50;
% general weight for metabolic reactions
weights.metabolic = 50;
% weight for making a reaction reversible
weights.reverse = 25;
% weight for exchange reactions from the database
weights.exchange = 100000;
% weight for allowed exchange / uptake reactions (from previous solutions)
weights.uptake = 1;
% weight for sink reactions if they should be determined in the gap-filling program
weights.sink = 0;
% weights for reactions already contained in the model
weights.model = 0;

%~~~~~~~~~~~~~ metabolite classification resources ~~~~~~~~~~~~~%

% ChEBI ontology DAG workspace
chebiOntologyWS = '~/IMIC/Program/matlab/COMMIT/data/metabolite-classification/ontologyGraph.mat';

% tab-separated ontology file
ontologyFile = '~/IMIC/Program/matlab/COMMIT/data/metabolite-classification/chebi_ontology.csv';

% brite Hierarchy file
briteFile = '~/IMIC/Program/matlab/COMMIT/data/metabolite-classification/briteHierarchy_ext.csv';
