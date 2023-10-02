%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to align repick average of each doublet with a reference
% and transform all the alignment to an updated table.
% dynamoDMT v0.2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This might be changed to accommodate the intraAlnRepick
%
%%%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/london/data0/20220404_TetraCU428_Tip_TS/ts/tip_CP_dPhi/';

%% Input
pixelSize = 8.48;
boxSize = 96;
filamentRepickListFile = 'filamentRepickList.csv';
particleDir = sprintf('%sparticles_repick', prjPath);
alnDir = sprintf('%sintraAln_repick', prjPath);
previewDir =[particleDir '/preview']; % created from previously
mw = 12; % Number of parallel workers to run
gpu = [0:5]; % Alignment using gpu
initRefFile = 'reference_all.em'; % Use the best reference that you have. If not, the updated ref from previous step
coneFlip = 0; % Search for polarity. 1 is yes. Recommended to pick with polarity and set to 0
alnLowpass = 40; % Angstrom
avgLowpass = 30; % Angstrom
zshift_limit = 6; % ~4nm shift limit in pixel for 8 nm repeat, 8nm shift for 16-nm repeat
newRefFile = 'reference_repick.em';


%%
filamentList = readcell(filamentRepickListFile, 'Delimiter', ',');
noFilament = length(filamentList);
template = dread(initRefFile);
newTemplate = zeros(boxSize, boxSize, boxSize);


alnLowpassPix = round(pixelSize/alnLowpass*boxSize);
mkdir(previewDir)
%v0.2b
% Need to go into alnDir to read the intraAln project
cd(alnDir)

%% Calculate the alignment of the filamentAverage to the initial reference
% transform the corresponding table for all particles
for idx = 1:noFilament
	% v0.2b
	%aPath = ([particleDir '/' filamentList{idx} '/template.em']); % Read the path of the alignment project average
	%tPath = ([particleDir '/' filamentList{idx} '/crop.tbl']); 
	aPath = ddb([filamentList{idx} ':a']); % Read the path of the alignment project average
	tPath = ddb([filamentList{idx} ':rt']);
	filamentAvg = dread(aPath);
    	disp(filamentList{idx})
	if coneFlip > 0
  		sal = dalign(dynamo_bandpass(filamentAvg,[1 alnLowpassPix]), dynamo_bandpass(template,[1 alnLowpassPix]),'cr',10,'cs',5,'ir',360,'is',10,'dim', boxSize, 'limm',1,'lim',[5,5,zshift_limit],'rf',5,'rff',2, 'cone_flip', 1); % cone_flip
	else
		% For repick with initial angle, restrict ir search
		%sal = dalign(dynamo_bandpass(filamentAvg,[1 alnLowpassPix]), dynamo_bandpass(template,[1 alnLowpassPix]),'cr',10,'cs',5,'ir',360,'is',10,'dim', boxSize, 'limm',1,'lim',[5,5,zshift_limit],'rf',5,'rff',2); % cone_flip
  		% For normal repick without set initial angle
		sal = dalign(dynamo_bandpass(filamentAvg,[1 alnLowpassPix]), dynamo_bandpass(template,[1 alnLowpassPix]),'cr',10,'cs',5,'ir',360,'is',5,'dim',boxSize, 'limm',1,'lim',[5, 5, zshift_limit],'rf',5,'rff',2); % no cone_flip
	end
	
	%dview(sal.aligned_particle);
	% Write out preview
	newTemplate = newTemplate + sal.aligned_particle;
	filt_aligned_particle = dynamo_bandpass(sal.aligned_particle, [1 round(pixelSize/avgLowpass*boxSize)]);
	img = sum(filt_aligned_particle(:,:,floor(boxSize/2) - 10: floor(boxSize/2) + 10), 3);
	imwrite(mat2gray(img), [previewDir '/' filamentList{idx} '_aln.png'])
	
	% Read last table from alignment
	tFilament = dread(tPath);
	% Read last transformation & applied to table
	tFilament_ali = dynamo_table_rigid(tFilament, sal.Tp);
	% Write table
	dwrite(tFilament_ali, [particleDir '/' filamentList{idx} '/aligned.tbl']);
	% Write aligned intraAvg
	dwrite(sal.aligned_particle, [alnDir '/avg/' filamentList{idx} '_aln.em']);

end

cd ..

%% Calculate average
newTemplate = newTemplate/noFilament;
dwrite(newTemplate, newRefFile);
