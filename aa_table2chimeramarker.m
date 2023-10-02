%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to move the origin of the subtomo for new table, cmm file & imod file
% dynamoDMT v0.2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Box center = BoxSize/2 + 1
% Shift vector = NewOrigin - Box center

%%%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/london/data0/20221128_TetraCU428Membrane_26k_TS/singlet/';

%%%%%%%%
%% Input
docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath);
modelDir = sprintf('%smodels', prjPath);
recSuffix = '_rec'; % Without the .mrc
tableAlnFileName = 'merged_particles_align.tbl'; % merge particles after particle alignment for robust 
shiftVector = [0 0 0]; % Shift vector in pixel, measure from the map just like relion (-6 for Atubule, 10 for B-tubule)

outName = 'doublet_center'; % Output will be name outName + tomoName + file extension
radius = 6;


% Read the list of filament to work with
tOri = dread(tableAlnFileName);

% Transform table
T = dynamo_rigid('shifts', -shiftVector); % create transformation
tOri_shift = dynamo_table_rigid(tOri, T); % transform table

% Reset the origin
tOri_adjusted = tOri_shift;

tOri_adjusted(:, 24:26) = tOri_shift(:, 24:26) + round(tOri_shift(:, 4:6));
tOri_adjusted(:, 4:6) = tOri_shift(:, 4:6) - round(tOri_shift(:, 4:6));

fileID = fopen(docFilePath); D = textscan(fileID,'%d %s'); fclose(fileID);
tomoID = D{1,1}'; % get tomogram ID
nTomo = length(D{1,2}); % get total number of tomograms

for idx = 1:nTomo
    tomo = D{1,2}{idx,1};
    tomono = D{1,1}(idx);
    [tomoPath,tomoName,ext] = fileparts(tomo);
    tomoName = strrep(tomoName, recSuffix, ''); % Remove the rec part of the name
    tTomo = tOri_adjusted(tOri_adjusted(:,20) == tomono, :);
    if isempty(tTomo) == 1
        continue;
    end
	dwrite(tTomo, [modelDir '/' outName '_' tomoName '.tbl']);
	dynamo_table2chimeramarker([modelDir '/' outName '_' tomoName '.cmm'], tTomo, radius);
	dlmwrite([modelDir '/' outName '_' tomoName '.txt'], tTomo(:, 23:26), 'delimiter', ' ');
end

% You can do 'point2model singlet.txt singlet.mod'