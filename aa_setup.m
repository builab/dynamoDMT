%%%%%%%%%%%%%%%%%%%%
% Set up the project
% dynamoDMT v0.1
%%%%%%%%%%%%%%%%%%%%
% Import tomogram & create catalog using the GUI
% For asymetric wedge, you might want to check it using this guide
% https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Indicating_the_missing_wedge
% 

%%%%%% Instructions %%%%%%
% Step 1: Change the following paths to the right directories.
%%%% modDir contains the files with the tomograms and path references the root directory for analysis
% Step 2: Indicate the delimiter that references all mod files of interest.
%%%% CU428*/tip_CP.mod
% Step 3: Indicate the suffix that needs to be removed in the various .mod files.
%%%% In the same way, indicate the suffix that needs to be added to reference the mrc files.
% Step 4: Change the pixel size
%%%% 

%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

%%%%%%% Variables subject to change %%%%%%%%

modDir = '/london/data0/20220404_TetraCU428_Tip_TS/ts/';
prjPath = '/london/data0/20220404_TetraCU428_Tip_TS/ts/tip_CP_dPhi/';
modFileDelimiter = 'CU428*/tip_CP.mod';
stringToBeRemoved = '/tip_CP.mod';
recSuffix = '_rec';
apixel = '8.48';

%%%%%%%%% Do not change anything under here %%%%%%%%%%

pathToModelScript = fullfile(sprintf('%sdynamoDMT', prjPath), 'createModTxt.sh');
cmdStr = [pathToModelScript ' ' modDir ' ' modFileDelimiter ' ' prjPath];
system(cmdStr);
modelfile = sprintf('%smodfiles.txt', prjPath);
catalogs = sprintf('%scatalogs', prjPath);
listOfTomograms = sprintf('%scatalogs/listOfTomograms', prjPath);
c001Dir = sprintf('%scatalogs/c001', prjPath);
modelDir = sprintf('%smodels/', prjPath);
mkdir(catalogs);
mkdir(listOfTomograms);
mkdir(c001Dir);
mkdir(modelDir);

catPath = sprintf('%scatalogs/c001', prjPath);
docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath); % This has to be consistent with the vll after catalog creation
vllFilePath = sprintf('%scatalogs/tomograms.vll', prjPath);

%% Populate the necessary files into a single location. Use this directory
% to create .vll and .doc files

% $./vllAndDocScript.sh modDir modelDestination modelfile stringTobeRemoved docFilePath vllFilePath 

pathToModelScript = fullfile(sprintf('%sdynamoDMT', prjPath), 'vllAndDocScript.sh');

% v0.2b change
cmdStr = [pathToModelScript ' ' modDir ' ' listOfTomograms ' ' modelfile ' ' stringToBeRemoved ' ' docFilePath ' ' vllFilePath ' ' recSuffix '.mrc ' apixel];
system(cmdStr);

% Create new catalogue from vll file, delete old one
dcm('c', catPath, 'fromvll', vllFilePath, 'delete_old', 1)

% Create text files using model2point
% Imod coordinate should be in text file, clicking along the filament (no direction needed)
% model2point -Contour imodModel.mod imodModel.txt

pathToModelScript = fullfile(sprintf('%sdynamoDMT', prjPath), 'model2pointscript.sh');

cmdStr = [pathToModelScript ' ' modDir ' ' modelDir ' ' modelfile ' ' stringToBeRemoved];
system(cmdStr);
