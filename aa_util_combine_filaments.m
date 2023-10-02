%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utility script to combine many filament tables into one
% dynamoDMT v0.2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Before Running Script %%%%%%%%%%
%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/london/data0/20221128_TetraCU428Membrane_26k_TS/doublet_16nm/';

%% Input
filamentListFile = 'filamentList.csv';
particleDir = sprintf('%sparticles_segment', prjPath);
tableFileName = 'merged_particles_segment.tbl'; % merged particles table all

%%
filamentList = readcell(filamentListFile, 'Delimiter', ',');
noFilament = length(filamentList);
 
% Combine all the particles into one table
% create table array
targetFolder = {};
tableName ={};

for idx = 1:noFilament
	targetFolder{idx} = [particleDir '/' filamentList{idx}];
	tableName{idx} = [particleDir '/' filamentList{idx} '/crop.tbl'];
end


% create ParticleListFile object (this object only exists temporarily in matlab)
plfClean = dpkdata.containers.ParticleListFile.mergeDataFolders(targetFolder,'tables',tableName);

% create and write the .star file
%plfClean.writeFile(starFileName)

% create merged table
tMergedClean = plfClean.metadata.table.getClassicalTable();
dwrite(tMergedClean,tableFileName)