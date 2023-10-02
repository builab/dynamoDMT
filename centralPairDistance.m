% Central Pair Distance Script 
% Draws 3D structures of the central pairs 
% This script uses distance2curve.m, shadedErrorBar.m, and extractTableData.sh 
% NOTE: extractTableData.sh needs read and write permissions $chmod 777 *
% Plots original distance graphs, each line representing an individual central pair
% Aggregates data into one mean curve
% Exports the final overall mean distance of all chosen central pairs as an eps file
% Exports csv file with table format [x, y, y+std, y-std, std]

% Output is found in the output directory

% % % % % % % % % INSTRUCTIONS % % % % % % % % % 
% Ensure the input path is correct
% Input contains particles and output is project path (current directory)
% CentralPair should contain list of all central pairs of interest, located in the same directory
% Table files in this case have to be in /london/data0/20220404_TetraCU428_Tip_TS/ts/cp_transition_analysis/particles/CU...../
% % % % % % % % % % % % % % % % % % % % % % % % %



inputPath = '/london/data0/20220404_TetraCU428_Tip_TS/ts/cp_transition_analysis/particles/';
CentralPair = ["CU428m_TS022" "CU428_TS013" "CU428_TS019" "CU428_TS218" "CU428_TS216" "CU428_TS229" ];
pixelSize = 0.848; % nm

%  "CU428_TS216" "CU428_TS229"



%%%%%%%%% Do not change anything under here %%%%%%%%%%
numberOfCentralPairs = length(CentralPair);
drawGraphs(numberOfCentralPairs, inputPath, CentralPair, pixelSize);

function drawGraphs(numberOfCentralPairs, inputPath, CentralPair, pixelSize)
    clc;
    outputPath = pwd;
    mkdir output;
    modDir1 = sprintf('%s/output/log1.txt', outputPath);
    modDir2 = sprintf('%s/output/log2.txt', outputPath);

    midpointI = [];
    mIndexes = [];
    distanceM = {};
    distanceIndex = 1;
    
    for i = 1:numberOfCentralPairs

        pathToModelScript = fullfile(sprintf('%s', outputPath), 'extractTableData.sh');
        cmdStr = [pathToModelScript ' ' modDir1 ' ' modDir2 ' ' inputPath ' ' sprintf('%s', CentralPair(i))];
        system(cmdStr);

        %Parse through data and plot original cp in xyz space
        m1 = readtable(modDir1);
        m1 = m1{:,:};
        CS = cat(1,0,cumsum(sqrt(sum(diff(m1,[],1).^2,2))));
        dd = interp1(CS, m1, unique([CS(:)' linspace(0,CS(end),100)]),'pchip');

        %Plot c1, c2, and distances between particles
        figure('Name', CentralPair(i)), hold on
        plot3(m1(:,1),m1(:,2),m1(:,3),'.b-');
        plot3(dd(:,1),dd(:,2),dd(:,3),'.r-');
        axis image, view(3), legend({'Original','Interp. Spline'});

        m2 = readtable(modDir2);
        m2 = m2{:,:};
        CS2 = cat(1,0,cumsum(sqrt(sum(diff(m2,[],1).^2,2))));
        dd = interp1(CS2, m2, unique([CS2(:)' linspace(0,CS2(end),100)]),'pchip');

        hold on
        plot3(m2(:,1),m2(:,2),m2(:,3),'.b-');
        plot3(dd(:,1),dd(:,2),dd(:,3),'.r-');
        axis image, view(3), legend({'Original','Interp. Spline'});
        
        %Find the distance between each point
        curvexy = m2;
        mapxy = m1;
        [xy,distance,t] = distance2curve(curvexy,mapxy,'linear');
        
        for idx = 1:length(xy)
            pt1 = m1(idx,:);
            pt2 = xy(idx,:);
            plot3([pt1(1) pt2(1)],[pt1(2) pt2(2)],[pt1(3) pt2(3)]);
        end
        
        hold off
%         h = figure();
        distance(distance >= 50) =[];
        x = (1:1:length(distance))';
        mean_arr = movmean(distance,3);
        median_arr = medfilt1(distance);
%         plot([1:1:length(distance)], mean_arr);
        distanceM{distanceIndex,1} = x;
        distanceM{distanceIndex,2} = distance;
        distanceM{distanceIndex,3} = mean_arr;
        distanceM{distanceIndex,5} = median_arr;
        
%       Find the linear slope of 10 points and find the slope with the most
%       negative position
        T = table(x,distance);
        [p,~,mu] = polyfit(T.x, median_arr, 5);
        distanceM{distanceIndex, 4} = polyval(p,x,[],mu);
        distanceIndex = distanceIndex + 1;
%         hold on
%         plot(x,f);
%         hold off
        m=Inf;
        ss=[];
        index=1;
        for j=10:1:length(distance)
           p=polyfit(T.x(j-9:j,:), T.distance(j-9:j,:), 1);
           ss=[ss;p];
           if (m>p(1))
               index=idivide((2*j-9),int16(2))+1;
               m=p(1);
           end
        end
        mIndexes = [mIndexes;index];  
        
        %Find the first 5 and last 5 points; average; then find midpoint
         s=mean(distance(1:10,:));
         temp=tail(T,10);
         E=mean(temp.distance);
         fMean= (s + E)/2;
         midpoint = 0;
         for k=1:1:length(distance)
             if (fMean<distance(k))
                 midpoint=k;
             end
         end
         midpointI= [midpointI;midpoint];
    end

    for index = 1:length(midpointI)
        distanceM{index,1} = (distanceM{index,1} - double(midpointI(index)))*.pixelSize;
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
        standardD = [standardD; std(cur, "omitnan")*pixelSize];
%         test = [test; summation/n];
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
%   index based on 1 pixel : 0.848nm
    meanM = [(startValue:pixelSize:endValue)' meanM*pixelSize];
    
    
    %Plot original
    figure('Name', 'Original');
    hold on;
    for i = 1:size(distanceM,1)
        plot(distanceM{i,1}, distanceM{i,2});
    end
    hold off;
    
%     %Plot mean
%     figure('Name', '3rd order moving mean');
%     hold on;
%     for i = 1:size(distanceM,1)
%         plot(distanceM{i,1}, distanceM{i,3});
%     end
%     hold off;
%     
%     %Plot median filter
%     figure('Name', '3rd order median filter');
%     hold on;
%     for i = 1:size(distanceM,1)
%         plot(distanceM{i,1}, distanceM{i,5});
%     end
%     hold off;
    
    %Plot overall mean with standard d and generate eps file
    gcf = figure('Name', 'Mean curve with std D');
    hold on;
    shadedErrorBar(meanM(:,1), meanM(:,2), standardD(:,1), 'lineProps', '-b', 'transparent', 1)
    set(gcf,'units','pixels','Position', [100 100 700 375]);
    plotMidPoint = plot(meanM(midpointI,1), meanM(midpointI, 2));
    set(plotMidPoint, 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
    set(gca,'FontSize',13,'FontName','Arial');
    set(gca,'XTick', -30:5:40);
    set(gca,'XTickLabel',-30:5:40);
    set(gcf, 'PaperPositionMode', 'auto');
    hold off;
    print(sprintf('%s/output/CPDistanceGraph', outputPath),'-depsc2');
    
%     %Plot overall moving mean with standard d
%     figure('Name', 'Moving Mean curve with std D');
%     hold on;
%     shadedErrorBar(meanM(:,1), movmean(meanM(:,2),3), standardD(:,1), 'lineProps', '-b', 'transparent', 1)
%     hold off;
    
    %Plot overall median filter with standard d
    figure('Name', 'Median filter curve with std D');
    hold on;
    shadedErrorBar(meanM(:,1), medfilt1(meanM(:,2)), standardD(:,1), 'lineProps', '-b', 'transparent', 1)
    hold off;
    
%     %Plot overall mean with standard error
%     figure('Name', 'Mean curve with std Err');
%     hold on;
%     errorbar(meanM(:,1), meanM(:,2),standardD(:,1)/sqrt(size(meanM(:,2),2)));
% %     plot(meanM(:,1), meanM(:,2));
%     hold off;
    
    final = [meanM(:,1) meanM(:,2) meanM(:,2)+standardD(:,1) meanM(:,2)-standardD(:,1) standardD(:,1)];
    csvPath = sprintf('%s/output/MeanCPDistance.csv', outputPath);
    csvwrite(csvPath, final);
    
    fprintf('Complete! Data was saved as [x,y, y+std, y-std, std] to %s\n', csvPath);
    fprintf('Mean figure was exported to %s/output/CPDistanceGraph.eps\n', outputPath);
    fprintf('Mean body central pair distance: %d nanometers\n', mean(final(1:25,2)));
    fprintf('Mean tip central pair distance: %d nanometers\n', mean(final(35:end,2)));
    
    delete(modDir1);
    delete(modDir2);
    
end
