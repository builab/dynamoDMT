%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to generate initial average
% dynamoDMT v0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/london/data0/20220404_TetraCU428_Tip_TS/ts/tip_complex/';

%%%%%%%%

% Input
docFilePath = sprintf('%scatalogs/tomograms.doc', prjDir);
filamentListFile = sprintf('%sfilamentList.csv', prjDir);
modelDir = sprintf('%smodels', prjDir);
particleDir = sprintf('%sparticles', prjDir);
boxSize = 168; % Extracted subvolume size
mw = 12; % Number of parallel workers to run
lowpass = 24; % Filter the initial average to 60 (50?) Angstrom

% Change by Huy
filamentList = readcell(filamentListFile, 'Delimiter', ',');

% Updated filamentList with tip complex
filamentListOneFile = 'filamentListOne.csv';
filamentListOne = {};

% Crop & generate initial average
for idx = 1:length(filamentList)
  tableName = [modelDir '/' filamentList{idx} '.tbl'];
  targetFolder = [particleDir '/' filamentList{idx}];
  disp(['Reading ' filamentList{idx}]);
  tImport = dread(tableName);
  %tImport = tImport(1, :);		

  
  % Cropping subtomogram out
  dtcrop(docFilePath, tImport, targetFolder, boxSize, 'mw', mw); % mw = number of workers to run
  
  % Plotting (might not be optimum since plotting everything here)
  dtplot(tImport, 'pf', 'oriented_positions');
  
end

% Merge only the list with particles_00001
count = 0
tCombine = [];
mkdir([particleDir '/one']);
for idx = 1:length(filamentList)
  tableName = [particleDir '/' filamentList{idx} '/crop.tbl'];
  tOri = dread(tableName);
  if isfile([particleDir '/' filamentList{idx} '/particle_000001.em'])
	disp([filamentList{idx} ' is valid']);
	copyfile([particleDir '/' filamentList{idx} '/particle_000001.em'], [particleDir '/one/particle_' sprintf('%06d.em', count + 1)]);
  	filamentListOne{count + 1, 1} = filamentList{idx};
  	tCombine = [tCombine; tOri(1, :)];
  	tCombine (:, 1) = 1:size(tCombine, 1)';
  	dwrite(tOri(1,:), [particleDir '/' filamentList{idx} '/crop_one.tbl']);
  	count = count + 1;
  end
end
writecell(filamentListOne, filamentListOneFile);
dwrite(tCombine, [particleDir '/one/crop.tbl'])

oa = daverage([particleDir '/one'], 't', tCombine, 'fc', 1, 'mw', mw);
dwrite(dynamo_bandpass(oa.average, [1 lowpass]), 'init_template.em');

