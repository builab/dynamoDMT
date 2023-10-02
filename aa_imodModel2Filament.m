%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to convert IMOD model to filament torsion model
% dynamoDMT v0.2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using GUI https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Filament_model
% Imod coordinate should be in text file, clicking along the filament (no direction needed)
% model2point -Contour imodModel.mod imodModel.txt
% NOTE: If the filament twist (microtubule/CP), we need to define subunits_dphi to describe the torsion.
% however, it might be related to the polarity of the filament (- or + sign).
% NOTE: filament number to Column 23


%%%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/london/data0/20220404_TetraCU428_Tip_TS/ts/tip_CP_dPhi/';

%%%%%%%%

%% Input
docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath);
modelDir = sprintf('%smodels', prjPath);
c001Dir = sprintf('%scatalogs/c001', prjPath);
recSuffix = '_rec'; % The suffix path without .mrc
pixelSize = 8.48; % Angstrom per pixel
periodicity = 82.8; % Using 84.5 of doublet, 82.8 for CP tip, 86 for CP base
subunits_dphi = 0.72;  % For the tip CP 0.72, base CP 0.5, doublet 0
subunits_dz = periodicity/pixelSize; % in pixel repeating unit dz = 8.4 nm = 168 Angstrom/pixelSize
filamentListFile = sprintf('%sfilamentList.csv', prjPath);
minPartNo = 4; % Minimum particles number per Filament

% loop through all tomograms
fileID = fopen(docFilePath); D = textscan(fileID,'%d %s'); fclose(fileID);
tomoID = D{1,1}'; % get tomogram ID
nTomo = length(D{1,2}); % get total number of tomograms

filamentList = {};

%% Loop through tomograms
for idx = 1:nTomo
    tomo = D{1,2}{idx,1};
    [tomoPath,tomoName,ext] = fileparts(tomo);
    % Modify specific to name
    tomoName = strrep(tomoName, recSuffix, ''); % Remove the rec part of the name from IMOD
    imodModel = [modelDir '/' tomoName '.txt'];
    modelout = strrep(imodModel, '.txt', '.omd');
    
    allpoints = load(imodModel);
    
    m = {}; % Cell array contains all filament
    contour = unique(allpoints(:, 1));
    % Loop through filaments
    for i = 1:length(contour)
        filamentid = contour(i);
        points = allpoints(allpoints(:, 1) == filamentid, 2:4);
        m{i} = dmodels.filamentWithTorsion();
        m{i}.subunits_dphi = subunits_dphi;
        m{i}.subunits_dz = subunits_dz;
        m{i}.name = [tomoName '_' num2str(contour(i))];
        % Import coordinate
        m{i}.points = points;
        % Create backbone
        m{i}.backboneUpdate();
        % Update crop point (can change dz)
        m{i}.updateCrop();
        % Link to catalog
        m{i}.linkCatalogue(c001Dir, 'i', idx);
        m{i}.saveInCatalogue();
        
        % Testing this block
        t = m{i}.grepTable();
        
        % 0.2b addition
        t(:,23) = contour(i);
        if (size(t, 1) < minPartNo)
        	disp(['Skip ' tomoName ' Contour ' num2str(contour(i)) ' with less than ' num2str(minPartNo) ' particles'])
        	continue
        end
        % Add the good to the list
        filamentList{end + 1, 1} = [tomoName '_' num2str(contour(i))];
		dwrite(t, [modelDir '/' tomoName '_' num2str(contour(i)) '.tbl']);
        
        % Optional for visualization of table
        dtplot(t, 'pf', 'oriented_positions');
        view(-230,30);axis equal;
        hold on;
    end
    print([modelDir '/' tomoName] , '-dpng');
    close all;
    
    % Write the DynamoModel
    dwrite(m, modelout)

end

%% Write out list file
writecell(filamentList, filamentListFile);
