%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to go through the filament List and exclude any filament with less than X particles
% by making a new filament list and back up the old one
% This avoids later script crashed or waste of processing time
% dynamoDMT v0.2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/london/data0/20220404_TetraCU428_Tip_TS/ts/tip_CP_dPhi/';

%%%%%%%%

%% Input
minPartNo = 10; % Minimum number of particles, 4 is ok, better 10 to get better alignment
useRepick = 1; % Working with repick particles instead of initial particles
boxSize = 96; % Original extracted subvolume size

%%
if useRepick > 0
	particleDir = sprintf('%sparticles_repick', prjPath);
	filamentListFile = sprintf('%sfilamentRepickList.csv', prjPath);
else
	particleDir = sprintf('%sparticles', prjPath);
	filamentListFile = sprintf('%sfilamentList.csv', prjPath);
end


%% Reading the list file
filamentList = readcell(filamentListFile, 'Delimiter', ',');

filamentListNew = {};

%% Looping through every filament
for idx = 1:length(filamentList)
    tableName = [particleDir '/' filamentList{idx} '/crop.tbl'];
    try
    	tOri = dread(tableName);
    	if size(tOri, 1) < minPartNo
            disp(['Skip ' filamentList{idx} ': less than' num2str(minPartNo) ' particles'])
        	continue
    	end
    catch
    	disp(['Skip ' filamentList{idx} ': due to error'])
    end
 	% If it is cropping out
  	filamentListNew{end + 1, 1} = filamentList{idx};
end

%% 0.2b Writing new list
if size(filamentListNew, 1) < size(filamentList, 1)
	% Backup old filamentList & write new one
	copyfile(filamentListFile, [filamentListFile '.bak']);
	writecell(filamentListNew, filamentListFile);
end