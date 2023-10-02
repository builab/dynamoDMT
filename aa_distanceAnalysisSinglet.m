% Central Pair Distance Script 
% Draws 3D structures of the central pairs 
% This script uses distance2curve.m, shadedErrorBar.m
% Plots original distance graphs, each line representing an individual central pair
% The best input is the repick table file due to the regularity of the points
% Aggregates data into one mean curve
% Pixel resolution was converted with ratio 1 pixel to 1.011 nanometers 
% .eps graph parameters (ex. font, fontsize, etc.) can be found between lines 217-229 and are subject to change
% Exports the final overall mean distance of all chosen central pairs as an eps file
% Exports csv file with table format [x, y, y+std, y-std, std]
% Modify the code from Max Tong to make it easier
% Output is found in the output directory

% % % % % % % % % INSTRUCTIONS % % % % % % % % % 
% Ensure the input path is correct
% Input contains particles and output is project path (current directory)
% Tomo should contain list of all central pairs of interest, located in the same directory
% Table files in this case have to be in /london/data0/20220404_TetraCU428_Tip_TS/ts/cp_transition_analysis/particles/CU...../
% % % % % % % % % % % % % % % % % % % % % % % % %

%%%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

%%% Input
prjPath = '/london/data0/20221128_TetraCU428Membrane_26k_TS/singlet/';
% Adjust the input to particles or particles_repick (ideally particles_repick)
inputPath = sprintf('%sparticles/', prjPath);
outputPrefix = 'singlet'; % CP or singlet. Do it strictly

% Selected tomograms for plotting
%tomograms = ["CU428lowmag_11", "CU428lowmag_14", "CU428lowmag_22", "CU428lowmag_29"];
tomograms = ["CU428lowmag_07", "CU428lowmag_11", "CU428lowmag_14", "CU428lowmag_22", "CU428lowmag_29"];

% hard code here
if strcmp(outputPrefix, 'CP') == 1
	microtubuleList = [1]; % For CP
else 
	microtubuleList = [1:9]; %for A-tubule
end
pixelSize = 10.11; % Angstrom
outlierDist = 1000; % Angstrom for CP 500 is safe
fitGroup = 5;
plotTomo = 0; % Turn to 0 if no need

% Limit to which is a parallel to eliminate noise
distThresAngst = 320; % In Angstrom

addpath(genpath(prjPath));

%%%%%%%%% Do not change anything under here %%%%%%%%%%
numberOfTomo = length(tomograms);

midpointI = [];
mIndexes = [];
distanceM = {};
distanceIndex = 1;
    
for i = 1:numberOfTomo 
    %Parse through data and plot original cp in xyz space  
    disp(tomograms(i))   
    % For A-tubule, need also 9 to 1
    for microtubuleId = microtubuleList
    	disp(['Microtubule ' num2str(microtubuleId)])
    	% Important: project the short microtubule to long microtubule to avoid problem
    	% Perhaps, the better one is compare Y of the tip point
        dVectors = [];    
        tbl1 = dread([inputPath sprintf('%s', tomograms(i)) '_' num2str(microtubuleId) '/crop.tbl']);
        m1 = tbl1(:,4:6) + tbl1(:,24:26);
        
        % For singlet, also 9 to 1
        if microtubuleId == 9
        	tbl2 = dread([inputPath sprintf('%s', tomograms(i)) '_' num2str(1) '/crop.tbl']);
        else
        	tbl2 = dread([inputPath sprintf('%s', tomograms(i)) '_' num2str(microtubuleId+1) '/crop.tbl']);
        end

        m2 = tbl2(:,4:6) + tbl2(:,24:26);
 
 		% Check the orientation (base is up or down)
 		isBaseUp = 0;
 		if m1(1, 2) > m1(end, 2)
 			isBaseUp = 1;
 		end
 		if isBaseUp > 0
 			% Microtubule 2 is shorter
 			if m2(end, 2) > m1(end, 2)
 				mtemp = m1;
 				m1 = m2;
 				m2 = mtemp;
 			end
 		else
 			% Microtubule 2 is shorter
 			if m2(end, 2) < m1(end, 2)
 				mtemp = m1;
 				m1 = m2;
 				m2 = mtemp;
 			end
 		end	
 				   
        CS = cat(1,0,cumsum(sqrt(sum(diff(m1,[],1).^2,2))));
        dd = interp1(CS, m1, unique([CS(:)' linspace(0,CS(end),100)]),'pchip');

        %Plot c1, c2, and distances between particles


        CS2 = cat(1,0,cumsum(sqrt(sum(diff(m2,[],1).^2,2))));
        dd2 = interp1(CS2, m2, unique([CS2(:)' linspace(0,CS2(end),100)]),'pchip');
        
        if plotTomo > 0
        	figure('Name', tomograms(i)), hold on
        	plot3(m1(:,1),m1(:,2),m1(:,3),'.b-');
        	plot3(dd(:,1),dd(:,2),dd(:,3),'.r-');
        	axis image, view(3), legend({'Original','Interp. Spline'});
        	hold on
        	plot3(m2(:,1),m2(:,2),m2(:,3),'.b-');
        	plot3(dd2(:,1),dd2(:,2),dd2(:,3),'.r-');
        	axis image, view(3), legend({'Original','Interp. Spline'});
        end
        
        %Find the distance between each point
        curvexy = m2;
        mapxy = m1;
        [xy,distance,t] = distance2curve(curvexy,mapxy,'linear');
        
        for idx = 1:length(xy)
            pt1 = m1(idx,:);
            pt2 = xy(idx,:);
            if plotTomo > 0
            	plot3([pt1(1) pt2(1)],[pt1(2) pt2(2)],[pt1(3) pt2(3)]);
            end
            dVectors = [dVectors; pt2 - pt1];
        end
        
        %   Central pair rotational twist measurements
        mRot = [];
        for k = 1:length(dVectors)-1
            mRot = [mRot; vrrotvec(dVectors(k+1,:),dVectors(k,:))];
        end
        
        hold off
		% Cleaning
        distance(distance >= outlierDist/pixelSize) = [];
        x = (1:1:length(distance))';
        mean_arr = movmean(distance,3);
        median_arr = medfilt1(distance);
		% plot([1:1:length(distance)], mean_arr);
        distanceM{distanceIndex,1} = x;
        distanceM{distanceIndex,2} = distance;
        distanceM{distanceIndex,3} = mean_arr;
        distanceM{distanceIndex,5} = median_arr;
        
        % TODO This part of the code might be make a lot simpler with a distance threshold
		% Find the linear slope of 5 points (fitGroup) and find the slope with the most
		% negative position
		% For A-tubule, 5 points is better
        T = table(x,distance);
        [p,~,mu] = polyfit(T.x, median_arr, 5);
        distanceM{distanceIndex, 4} = polyval(p,x,[],mu);
        m=Inf;
        ss=[];
        index=1;
        % HUY: The end is very prone to error so, should exclude the end here
        for j=fitGroup:1:length(distance)-fitGroup
           p=polyfit(T.x(j-fitGroup+1:j,:), T.distance(j-fitGroup+1:j,:), 1);
           ss=[ss;p];
           if (m>p(1))
               index=idivide((2*j-fitGroup+1),int16(2))+1;
               m=p(1);
           end
        end
        if distanceM{distanceIndex, 2}(1) < distThresAngst/pixelSize
        	index = 1;
        end
        mIndexes = [mIndexes;index];
        
        %For simplicity midpoint = mIndexes
        midpointI = [midpointI; index];
        distanceIndex = distanceIndex + 1;

    end
end

% Manual part to filter data
%ind = find(midpointI > 13)
        
%   Central pair distance measurements
for index = 1:length(midpointI)
    distanceM{index,1} = (distanceM{index,1} - double(midpointI(index)))*pixelSize/10;
end
    
distanceM2 = {};
for i = 1:1:size(distanceM,1)
    distanceM2{i,1} = distanceM{i,2};
end
    
%Add zeros at the beginning of the shorter arrays
for i = 1:1:size(distanceM2,1)
   for j = 1:1:size(distanceM2,1)
        beg = numel(distanceM2{j,1}(1:midpointI(j)))-numel(distanceM2{i,1}(1:midpointI(i)));
        edd = numel(distanceM2{j,1}(midpointI(j):end))-numel(distanceM2{i,1}(midpointI(i):end));
        hed = max(0,beg);
        ed = max(0,edd);
        distanceM2{i,1}=[zeros(1, hed)'; distanceM2{i,1}; zeros(1, ed)'];
        midpointI(i) = midpointI(i) + hed;
    end
end
    
meanM = [];
standardD = [];
%Find mean but exclude values with 0s
for i = 1:1:length(distanceM2{1,1})
   cur = [];
   summation = 0;
   test = [];
   n = 0;
   for j = 1:1:size(distanceM2,1)
       if (distanceM2{j,1}(i) ~= 0)
           cur = [cur; distanceM2{j,1}(i)];
           summation = summation + distanceM2{j,1}(i);
           n = n + 1;
       end
   end
   meanM = [meanM; mean(cur)];
   standardD = [standardD; std(cur, "omitnan")*pixelSize/10];
end
    
startValue = 0;
endValue = 0;
for i = 1:size(distanceM,1)
    if (startValue > distanceM{i,1}(1))
   		startValue = distanceM{i,1}(1);
    end
end
    
for i = 1:size(distanceM,1)
    if (endValue < distanceM{i,1}(length(distanceM{i,1})))
        endValue = distanceM{i,1}(length(distanceM{i,1}));
    end
end
    
%   Enumerate starting from first point till last point and change the
%   index based on 1 pixel : 1.011nm
meanM = [(startValue:pixelSize/10:endValue)' meanM*pixelSize/10];
        
%Plot original
figure('Name', 'Original');
set(gca, 'DefaultLineLineWidth', 2)
set(gca, 'LineWidth', 2)

hold on;
for i = 1:size(distanceM,1)
    plot(distanceM{i,1}, distanceM{i,2});
end
hold off;
ylim([10 80])
xlim([0 150])
% Hard code
set(gcf, 'PaperPositionMode', 'auto');
hold off;
print(sprintf('%s/%sOriginalDistanceGraph', prjPath, outputPrefix),'-depsc2');

    
%Plot overall mean with standard d and generate eps file
gcf = figure('Name', 'Mean curve with std D');
hold on;
set(gca, 'DefaultLineLineWidth', 2)

shadedErrorBar(meanM(:,1), meanM(:,2), standardD(:,1), 'lineProps', '-b', 'transparent', 1)
set(gcf,'units','pixels','Position', [100 100 700 375]);
plotMidPoint = plot(meanM(midpointI,1), meanM(midpointI, 2));
%set(plotMidPoint, 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
%set(gca,'FontSize',13,'FontName','Arial');
%set(gca,'XTick', -30:5:40);
%set(gca,'XTickLabel',-30:5:40);
set(gca, 'LineWidth', 2)
ylim([10 80])
xlim([0 150])

% Hard code
set(gcf, 'PaperPositionMode', 'auto');
hold off;
print(sprintf('%s/%sDistanceGraph', prjPath, outputPrefix),'-depsc2');
    
%     %Plot overall moving mean with standard d
%     figure('Name', 'Moving Mean curve with std D');
%     hold on;
%     shadedErrorBar(meanM(:,1), movmean(meanM(:,2),3), standardD(:,1), 'lineProps', '-b', 'transparent', 1)
%     hold off;
    
%Plot overall median filter with standard d
%figure('Name', 'Median filter curve with std D');
%hold on;
%shadedErrorBar(meanM(:,1), medfilt1(meanM(:,2)), standardD(:,1), 'lineProps', '-b', 'transparent', 1)
%ylimit([15 50])
%xlim([-25 150])
%hold off;
    
%     %Plot overall mean with standard error
%     figure('Name', 'Mean curve with std Err');
%     hold on;
%     errorbar(meanM(:,1), meanM(:,2),standardD(:,1)/sqrt(size(meanM(:,2),2)));
%     plot(meanM(:,1), meanM(:,2));
%     hold off;
    
final = [meanM(:,1) meanM(:,2) meanM(:,2)+standardD(:,1) meanM(:,2)-standardD(:,1) standardD(:,1)];
csvPath = sprintf('%s/MeanSingletDistance.csv', prjPath);
csvwrite(csvPath, final);
    
fprintf('Complete! Data was saved as [x,y, y+std, y-std, std] to %s\n', csvPath);
fprintf('Mean figure was exported to %s/CPDistanceGraph.eps\n', prjPath);
fprintf('Mean body central pair distance: %d nanometers\n', mean(final(1:25,2)));
fprintf('Mean tip central pair distance: %d nanometers\n', mean(final(35:end,2)));
