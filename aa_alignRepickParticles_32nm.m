%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to align repick particles
% and transform all the 32-nm alignment to an updated table.
% dynamoDMT v0.11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/london/data0/20220404_TetraCU428_Tip_TS/ts/base_CP/';

% Input
filamentListFile = 'filamentList.csv';
particleDir = sprintf('%sparticles_repick_32nm', prjPath);
mw = 12; % Number of parallel workers to run
gpu = [0:5]; % Alignment using gpu
template_name = 'ref_base_32nm.em';
tableFileName = 'merged_particles_repick_32nm.tbl'; % merged particles table all
starFileName = 'merged_particles_repick_32nm.star'; % star file name for merged particles
tableOutFileName = 'merged_particles_repick_32nm_align.tbl'; % merged particles table all
pAlnAll = 'pAlnRepickParticles32nm';
refMask = 'masks/mask_cp_base_bin4.em';
lowpass = 40; % 30Angstrom filter in fourier pixel
zshift_limit = 10; % ~8 nm


filamentList = readcell(filamentListFile, 'Delimiter', ',');
noFilament = length(filamentList);


template = dread(template_name);
 
% Combine all the particles into one table
% create table array
targetFolder = {};
tableName ={};

for idx = 1:noFilament
	targetFolder{idx} = [particleDir '/' filamentList{idx}];
	tableName{idx} = [particleDir '/' filamentList{idx} '/aligned.tbl'];
end


% create ParticleListFile object (this object only exists temporarily in matlab)
plfClean = dpkdata.containers.ParticleListFile.mergeDataFolders(targetFolder,'tables',tableName);

% create and write the .star file
plfClean.writeFile(starFileName)

% create merged table
tMergedClean = plfClean.metadata.table.getClassicalTable();
dwrite(tMergedClean,tableFileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform alignment of all particles with the ref
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dcp.new(pAlnAll,'t', tableFileName, 'd', targetFolder{1}, 'template', template_name, 'masks','default','show',0, 'forceOverwrite',1);
dvput(pAlnAll,'data',starFileName)
dvput(pAlnAll,'file_mask',refMask)

% set alignment parameters
dvput(pAlnAll,'ite', [2 2]);
dvput(pAlnAll,'dim', [144 144]);
dvput(pAlnAll,'low', [40 40]);
dvput(pAlnAll,'cr', [15 6]);
dvput(pAlnAll,'cs', [5 2]);
dvput(pAlnAll,'ir', [15 6]);
dvput(pAlnAll,'is', [5 2]);
dvput(pAlnAll,'rf', [5 5]);
dvput(pAlnAll,'rff', [2 2]);
dvput(pAlnAll,'lim', [zshift_limit zshift_limit]);
dvput(pAlnAll,'limm',[1 2]);
dvput(pAlnAll,'sym', 'c1'); % 
    
% set computational parameters
dvput(pAlnAll,'dst','matlab_gpu','cores',1,'mwa',mw);
dvput(pAlnAll,'gpus',gpu);

% check/unfold/run
dvrun(pAlnAll,'check',true,'unfold',true);

aPath = ddb([pAlnAll ':a']);
a = dread(aPath);
tPath = ddb([pAlnAll ':t:ite=last']); % This makes convertion to Relion better
dwrite(dread(tPath), tableOutFileName);
dwrite(dynamo_bandpass(a,[1 lowpass])*(-1),['result_alnRepickParticles_32nm_INVERTED_all.em']);
