%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to convert IMOD model to filament torsion model
% dynamoDMT v0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using GUI https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Filament_model
% Imod coordinate should be in text file, clicking along the filament (no direction needed)
% model2point -Contour imodModel.mod imodModel.txt
% Write out the filament list/folder for further processing as well
% NOTE: Important to have tomogram number 
% NOTE: If the filament twist (microtubule/CP), we need to define subunits_dphi to describe the torsion.
% however, it might be related to the polarity of the filament (- or + sign).

%%%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/london/data0/20220404_TetraCU428_Tip_TS/ts/tip_complex/';

%%%%%%%%

% Input
docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath);
modelDir = sprintf('%smodels', prjPath);
c001Dir = sprintf('%scatalogs/c001', prjPath);
pixelsize = 8.48; % Angstrom per pixel
periodicity = 168; % Using 16-nm of doublet for DMT 
subunits_dphi = 0;  % For the tip CP
subunits_dz = periodicity/pixelsize; % in pixel repeating unit dz = 8.4 nm = 168 Angstrom/pixelSize
filamentListFile = sprintf('%sfilamentList.csv', prjPath);

% loop through all tomograms
fileID = fopen(docFilePath); D = textscan(fileID,'%d %s'); fclose(fileID);
tomoID = D{1,1}'; % get tomogram ID
nTomo = length(D{1,2}); % get total number of tomograms

filamentList = {};

% Loop through tomograms
for idx = 1:nTomo
    tomo = D{1,2}{idx,1};
    [tomoPath,tomoName,ext] = fileparts(tomo);
    % Modify specific to name
    tomoName = strrep(tomoName, '_rec', ''); % Remove the rec part of the name
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
        
        % Add to the list
        filamentList{end + 1, 1} = [tomoName '_' num2str(contour(i))];

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
        % Optional for visualization of table
        %dtplot(['particles/' tomoName '_' num2str(contour(i)) '/crop.tbl'], 'pf', 'oriented_positions');
    end
    
    % Write the DynamoModel
    dwrite(m, modelout)

end

% Write out list file
writecell(filamentList, filamentListFile);
