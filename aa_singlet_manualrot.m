%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to manually assign rot angle for singlets
% dynamoDMT v0.2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/london/data0/20221128_TetraCU428Membrane_26k_TS/singlet/';

%%%%%%%%

%% Input
docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath);
filamentListFile = sprintf('%sfilamentList.csv', prjPath);
modelDir = sprintf('%smodels', prjPath);

narot = [174.73 -163 -126 -65 -18 -30.6 42 97 120];

% Read the list of filament to work with
filamentList = readcell(filamentListFile, 'Delimiter', ',');

%% Crop & generate initial average
for idx = 1:length(filamentList)
  tableName = [modelDir '/' filamentList{idx} '.tbl'];
  % Back up
  tableBackup = [modelDir '/' filamentList{idx} '_orig.tbl'];
  copyfile(tableName, tableBackup);
  disp(['Reading ' filamentList{idx}]);
  tImport = dread(tableName);
  tImport(:, 9) = narot(idx);
  dwrite(tImport, tableName); 
end
