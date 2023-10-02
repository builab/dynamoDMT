%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to apply alignment parameters to repick filament with torsion model
% Should have same parameters as imodModel2Filament
% dynamoDMT v0.2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using GUI https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Filament_model
% Now also printing output
% NOTE: v0.2b add CC filter using median(cc) - 3*mad(cc) 
% NOTE: v0.2b add different contour (Done)
% NOTE: v0.2b eliminate duplicate less than dTh (Done)
% NOTE: v0.2b find nearest neighbour for angle
% 0.2b Now use col 23 as filament number
% merged_particles.tbl intentionally switch the polarity of particles
% Seems to work better running on merged_particles.tbl instead of merged_particles_aln.tbl
% Implement a nearest neighbour for angle assignnment, still doesn't work yet due to Dynamo careless angle interpolation.

%%%%%%%% Before Running Script %%%%%%%%%%
%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/london/data0/20220404_TetraCU428_Tip_TS/ts/tip_CP_dPhi/';


%% Input
docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath);
modelDir = sprintf('%smodels_repick', prjPath);
origParticleDir = sprintf('%sparticles', prjPath);
particleDir = sprintf('%sparticles_repick', prjPath);
c001Dir = sprintf('%scatalogs/c001', prjPath);
pixelSize = 8.48; % Angstrom per pixel
periodicity = 82.8; % 82.8 for tipCP, xx for baseCP, 169 doublet
boxSize = 96;
mw = 12;
subunits_dphi = 0.72;  % For the tip CP 0.72, baseCP 0.5, doublet 0
subunits_dz = periodicity/pixelSize; % in pixel repeating unit dz = 8.4 nm = 168 Angstrom/pixelSize
filamentRepickListFile = sprintf('%sfilamentRepickList.csv', prjPath);
tableAlnFileName = 'merged_particles_align.tbl'; % merge particles before particle alignment for robust but must be merged_particles_align to use doInitialAngle
avgLowpass = 40; % Angstrom
dTh = 30; % Distance Threshold in Angstrom
doExclude = 1; % Exclude particles too close
doOutlier = 1; % Exclude outlier using CC using MAD
doInitialAngle = 0; % Only turn on for axoneme case now, absolutely not for microtubule

%% loop through all tomograms
fileID = fopen(docFilePath); D = textscan(fileID,'%d %s'); fclose(fileID);
tomoID = D{1,1}'; % get tomogram ID
nTomo = length(D{1,2}); % get total number of tomograms

tAll = dread(tableAlnFileName);

filamentRepickList = {};

%% Loop through tomograms
for idx = 1:nTomo
    tomo = D{1,2}{idx,1};
    [tomoPath,tomoName,ext] = fileparts(tomo);
    tomono = D{1,1}(idx);
    % Modify specific to name
    tomoName = strrep(tomoName, '_rec', ''); % Remove the rec part of the name
    tTomo = tAll(tAll(:,20) == tomono, :);
    if isempty(tTomo) == 1
        continue;
    end
      
    modelout =   [modelDir '/' tomoName '.omd'];
    contour = unique(tTomo(:, 23));
    
    m = {}; % Cell array contains all filament
    
    for i = 1:length(contour)        
        tContour = tTomo(tTomo(:, 23) == contour(i), :);
        phi = median(tContour(:, 9)); % Same as AA  
    
        % v0.2b Important: this step invert the Y axis, doing for each contour might help to check for polarity
        if doExclude > 0
            tContourEx = dpktbl.exclusionPerVolume(tContour, dTh/pixelSize);
            % Make sure to sort as before by particles number for not inverting angle
            tContour = sortrows(tContourEx, 1);
            display(['Exclude ' num2str(length(tContour) - length(tContourEx)) ' particles due to proximity']);
        end
       
        if doOutlier > 0
            cc = tContour(:, 10);
            x = median(cc);
            y = mad(cc);
            tContour = tContour(cc > x - 3*y, :);
            display(['Contour ' num2str(contour(i)) ': Exclude ' num2str(sum(cc <= x - 3*y)) ' particles']);
        end
        
        if isempty(tContour) == 1
            continue;
        end
        points = tContour(:, 24:26) + tContour(:, 4:6);
        
        % v0.2b Check for polarity in the original crop file
       
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
        t = m{i}.grepTable();

        %v0.2b addition
        if isempty(t) == 1
          	warning(['Skip: ' tomoName  '_' num2str(contour(i)) 'does not have any particles!']);
        	continue;
        end

        t(:,23) = contour(i); % Additing contour number (filament)
        
        if doInitialAngle > 0
            if abs(subunits_dphi) > 0     % Twist
            	midIndex = floor(size(t, 1)/2);
            	t(:, 9) = t(:, 9) - t(midIndex, 9) + phi; 
            else
           		t(:, 9) = phi; % This works will in case of doublet, in case of tip/base cp, make the middle value to this and then same shift
           	end
        end
        
        % Check point
        dwrite(t, [modelDir '/' tomoName '_' num2str(contour(i)) '.tbl']);
        targetFolder = [particleDir '/'  tomoName '_' num2str(contour(i))];
        
        % Cropping subtomogram out
        % 0.2b
        try
       		dtcrop(docFilePath, t, targetFolder, boxSize, 'mw', mw);
        	tCrop = dread([targetFolder '/crop.tbl']);
        	oa_all = daverage(targetFolder, 't', tCrop, 'fc', 1, 'mw', mw);
        	dwrite(dynamo_bandpass(oa_all.average, [1 round(pixelSize/avgLowpass*boxSize)]), [targetFolder '/average.em']);
        	% Average the middle region again
        	if size(tCrop, 1) > 15
            	midIndex = floor(size(tCrop, 1)/2);
            	tCrop = tCrop(midIndex - 3: midIndex + 4, :);
        	end
        	oa = daverage(targetFolder, 't', tCrop, 'fc', 1, 'mw', mw);
        	dwrite(dynamo_bandpass(oa.average, [1 round(pixelSize/avgLowpass*boxSize)]), [targetFolder '/template.em']);
        
        	% Plotting save & close. dtplot seems to error if only 1 particles
        	if size(tCrop, 1) > 1
            	dtplot(tCrop, 'pf', 'oriented_positions');
            	view(-230, 30); axis equal;
            	print([targetFolder '/repick_' tomoName '_' num2str(contour(i))] , '-dpng');
            	close all
        	end
       	catch
  			warning(['Skip: ' tomoName  '_' num2str(contour(i)) 'does not have enough particles!'])
  			continue;
  		end	
  		% If it is cropping out
  		filamentRepickList{end + 1, 1} = [tomoName  '_' num2str(contour(i))];   
    end
    % Write the DynamoModel
    dwrite(m, modelout);
end

%% Write filament list out
writecell(filamentRepickList, filamentRepickListFile);


