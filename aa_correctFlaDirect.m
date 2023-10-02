%%%%%%%%%%%%%%%%%%%%
% Correct the flagella direction to be compatible with the old AA & avoid coneflip later
% Can be incorporated into aa_setup if needed. It overwrites the model file from aa_setup
% dynamoDMT v0.2b
%%%%%%%%%%%%%%%%%%%%

%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

%%%%%%% Variables subject to change %%%%%%%%

modDir = '/london/data0/20220404_TetraCU428_Tip_TS/ts/'; % Tomo folder
prjPath = '/london/data0/20220404_TetraCU428_Tip_TS/ts/doublet/';
modelDir = sprintf('%smodels/', prjPath);
modFileName = 'doublet.mod';
modFixFileName = 'doublet_fix.mod';
flaDirectFile = sprintf('%sflaDirect.txt', prjPath); % Use it when above option is 1
docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath); % This has to be consistent with the vll after catalog creation
recSuffix = '_rec'; % Without the .mrc
pathToModelScript = fullfile(sprintf('%sdynamoDMT', prjPath), 'correctModDirect.sh');


%%%%%%%%% Do not change anything under here %%%%%%%%%%
modelfile = sprintf('%smodfiles.txt', prjPath);
catalogs = sprintf('%scatalogs', prjPath);
listOfTomograms = sprintf('%scatalogs/listOfTomograms', prjPath);
c001Dir = sprintf('%scatalogs/c001', prjPath);
modelDir = sprintf('%smodels/', prjPath);

% Read the flaDirect
fileFD = fopen(flaDirectFile); FD = textscan(fileFD,'%s %d'); fclose(fileFD);
tomoList = FD{1,1}'; % get tomogram List
tomoDirect = FD{1,2}'; % get total number of tomograms
flaDirect = containers.Map(tomoList, tomoDirect);

% Read tomo list
fileID = fopen(docFilePath); D = textscan(fileID,'%d %s'); fclose(fileID);
tomoID = D{1,1}'; % get tomogram ID
nTomo = length(D{1,2}); % get total number of tomograms

%% Loop through all tomograms
for idx = 1:nTomo
    tomo = D{1,2}{idx,1};
    [tomoPath,tomoName,ext] = fileparts(tomo);
    % Modify specific to name
    tomoName = strrep(tomoName, recSuffix, ''); % Remove the rec part of the name
    if isKey(flaDirect, tomoName) > 0
    	disp(flaDirect(tomoName));
		cmdStr = [pathToModelScript ' ' modDir tomoName '/' modFileName ' ' num2str(flaDirect(tomoName)) ' ' modDir tomoName '/' modFixFileName];
		disp(cmdStr);
		system(cmdStr);
		% Model2Point
		cmdStr2 = ['model2point -Contour ' modDir tomoName '/' modFixFileName ' ' modelDir tomoName '.txt'];
		disp(cmdStr2);
		system(cmdStr2);
    else
    	disp([tomoName ' does not have a flaDirect value']);
    end
end
