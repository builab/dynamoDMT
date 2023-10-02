% % % % % % % % % % % % % % % % % % % DOCUMENTATION % % % % % % % % % % % % % % % % % % %
% Creates .txt files using $model2point -Contour .mod .txt
% Draws A tubules beginning from B tubule termination
% Using the CP as the normal vector, finds a plane to the base of the A tubule and calculates the distance between the base of the CP to the plane along the normal vector
% Reports the difference between largest and shortest A tubule distance
% Output is saved in a txt file
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % INSTRUCTIONS % % % % % % % % % % % % % % % % % % %
% Ensure the input path is correct, leads to directory with the .mod files
% tomograms should contain list of all tomograms of interest with A tubules picked with CP, located in the same directory
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
inputPath = '/storage2/Thibault/20220512_TetraCU428Membrane/isonet/corrected_tomos/ATubulesWithCP';
tomograms = ["CU428_TS013", "CU428_TS277", "CU428_TS196", "CU428_TS259", "CU428_TS268"];
% , "CU428_TS277", "CU428_TS196", "CU428_TS259", "CU428_TS268"
%%%%%%%%% Do not change anything under here %%%%%%%%%%
m = drawGraphs(inputPath, tomograms);
disp(m);
csvPath = [inputPath '/distanceOutput.csv'];
csvwrite(csvPath,m);
% z = [489.7963
%   145.5213
%   144.8277
%   798.4588
%   369.9988]
% mean(z)
function longestD = drawGraphs(inputPath, tomograms)
    clc;
    lBD = [];
    lTD = [];
    mBT = [];
    mTT = [];
    for j = 1:length(tomograms)
        cmdStr= ['model2point -Contour ' inputPath '/' sprintf('%s', tomograms(j)) '.mod ' inputPath '/' sprintf('%s', tomograms(j)) '.txt'];
        system(cmdStr);
        curPath = [inputPath '/' sprintf('%s', tomograms(j)) '.txt'];
        m1 = readtable(curPath);
        m1 = m1{:,:};
        CP = [m1(length(m1)-1,2) m1(length(m1)-1,3) m1(length(m1)-1,4)];
        CP = [CP; [m1(length(m1),2) m1(length(m1),3) m1(length(m1),4)]];
        figure('Name', sprintf('%s: A Tubules and CP', tomograms(j)))
        plot3([CP(1,1) CP(2,1)],[CP(1,2) CP(2,2)],[CP(1,3) CP(2,3)], '.b-');
        m1(end,:) = [];
        m1(end,:) = [];
        
        mDT = [];
        mBCP = [];
        mTCP = [];
        for i = 1:2:length(m1)
            ATubuleA = [m1(i,2) m1(i,3) m1(i,4)];
            ATubuleB = [m1(i+1,2) m1(i+1,3) m1(i+1,4)];
            ATubuleBase = [];
            ATubuleTip = [];
            dA = norm(CP(1,:) - ATubuleA);
            dB = norm(CP(1,:) - ATubuleB);
            if dA < dB
                ATubuleBase = ATubuleA;
                ATubuleTip = ATubuleB;
            else
                ATubuleBase = ATubuleB;
                ATubuleTip = ATubuleA;
            end
%           Find which one is actually the base of the A tubule
            hold on
            plot3([m1(i,2), m1(i+1,2)],[m1(i,3), m1(i+1,3)],[m1(i,4), m1(i+1,4)],'.b-');
            
%           Base
            [xy,distance,t] = distance2curve(CP,ATubuleBase,'linear');
            for idx = 1:size(xy)
                pt1 = ATubuleBase(idx,:);
                pt2 = xy(idx,:);
            end
            plot3([pt1(1) pt2(1)],[pt1(2) pt2(2)],[pt1(3) pt2(3)]);
            plot3([CP(1,1) pt2(1)],[CP(1,2) pt2(2)],[CP(1,3) pt2(3)], '.b-');
            
%           Tip
            [xyz,distance,t] = distance2curve(CP,ATubuleTip,'linear');
            for idx = 1:size(xyz)
                pt1 = ATubuleTip(idx,:);
                pt2 = xyz(idx,:);
            end
            plot3([pt1(1) pt2(1)],[pt1(2) pt2(2)],[pt1(3) pt2(3)]);
            plot3([CP(2,1) pt2(1)],[CP(2,2) pt2(2)],[CP(2,3) pt2(3)], '.b-');
            mBCP = [mBCP; xy(idx,:)];
            mTCP = [mTCP; xyz(idx,:)];
        end
        hold off;
        
%       Base Projection distances
        figure('Name', sprintf('%s: Distance between each base/tip projection', tomograms(j)));
        hold on;
        [mSBP,mD] = knnsearch(mBCP, CP(1,:),'K',length(mBCP),'Distance','euclidean');
        
        mBP = [];
        for i = 1:(length(mD)-1)
            cur = mD(i+1) - mD(i);
            mBP = [mBP; cur];
        end
        
        if length(mBP) <= 8
            mBP = [mBP; zeros(8-length(mBP))]
        end
        mBT = [mBT, mBP]
        plot(mBP)
        
%       Tip Projection distances
        [mSTP,mD] = knnsearch(mTCP, CP(2,:),'K',length(mTCP),'Distance','euclidean');
        mTP = [];
        for i = 1:(length(mD)-1)
            cur = mD(i+1) - mD(i);
            mTP = [mTP; cur];
        end
        if length(mTP) <= 8
            mTP = [mTP; zeros(8-length(mTP))]
        end
        mTT = [mTT, mTP];
        plot(mTP)
        hold off;
        
%       Range of Base
        longestB = 0;
        shortestB = intmax;
        for i = 1:size(mBCP)
            cur = pdist([CP(1,:);mBCP(i,:)], 'euclidean');
            if (cur > longestB)
                longestB = cur;
            end
            if (cur < shortestB)
                shortestB = cur;
            end
            mBP = [mBP; cur];
        end
        
%       Range of Tip
        longestT = 0;
        shortestT = intmax;
        for i = 1:size(mTCP)
            cur = pdist([CP(2,:);mTCP(i,:)], 'euclidean');
            if (cur > longestT)
                longestT = cur;
            end
            if (cur < shortestT)
                shortestT = cur;
            end
            mDT = [mDT; cur];
        end
        

%         figure('Name', sprintf('%s: Distance Graph', tomograms(j)));
%         hold on;
%         plot(mBP);
%         plot(mDT);
%         hold off;
        lBD = [lBD; (longestB - shortestB)];
        lTD = [lTD; (longestT - shortestT)];
        delete(curPath);
    end
    longestD = [lBD, lTD];
end

% % % % % % % % % % % % % % % % % % % % DOCUMENTATION % % % % % % % % % % % % % % % % % % % 
% % Creates .txt files using $model2point -Contour .mod .txt
% % Draws A tubules beginning from B tubule termination
% % Using the CP as the normal vector, finds a plane to the base of the A tubule and calculates the distance between the base of the CP to the plane along the normal vector
% % Reports the difference between largest and shortest A tubule distance
% % Output is saved in a txt file
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% % % % % % % % % % % % % % % % % % % % INSTRUCTIONS % % % % % % % % % % % % % % % % % % % 
% % Ensure the input path is correct, leads to directory with the .mod files
% % tomograms should contain list of all tomograms of interest with A tubules picked with CP, located in the same directory
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% 
% inputPath = '/storage2/Thibault/20220512_TetraCU428Membrane/isonet/corrected_tomos/ATubulesWithCP';
% tomograms = ["CU428_TS013"];
% 
% % , "CU428_TS277", "CU428_TS196", "CU428_TS259", "CU428_TS268"
% 
% %%%%%%%%% Do not change anything under here %%%%%%%%%%
% 
% m = drawGraphs(inputPath, tomograms);
% disp(m);
% csvPath = [inputPath '/distanceOutput.csv'];
% csvwrite(csvPath,m);
% 
% function longestD = drawGraphs(inputPath, tomograms)
%     clc;
%     lD = [];
%     for j = 1:length(tomograms)
%         cmdStr= ['model2point -Contour ' inputPath '/' sprintf('%s', tomograms(j)) '.mod ' inputPath '/' sprintf('%s', tomograms(j)) '.txt'];
%         system(cmdStr);
%         
%         curPath = [inputPath '/' sprintf('%s', tomograms(j)) '.txt'];
%         m1 = readtable(curPath);
%         m1 = m1{:,:};
%         CP = [m1(length(m1)-1,2) m1(length(m1)-1,3) m1(length(m1)-1,4)];
%         CP = [CP; [m1(length(m1),2) m1(length(m1),3) m1(length(m1),4)]];
%         figure('Name', sprintf('%s: A Tubules and CP', tomograms(j)))
%         plot3([CP(1,1) CP(2,1)],[CP(1,2) CP(2,2)],[CP(1,3) CP(2,3)], '.b-');
%         
%         m1(end,:) = [];
%         m1(end,:) = [];
%         
%         mD = [];
%         mPCP = [];
%         
%         for i = 1:2:length(m1)
%             ATubule = [m1(i,2) m1(i,3) m1(i,4)];
%     %         ATubule = [ATubule; [m1(i+1,2) m1(i+1,3) m1(i+1,4)]];
%             hold on
%             plot3([m1(i,2), m1(i+1,2)],[m1(i,3), m1(i+1,3)],[m1(i,4), m1(i+1,4)],'.b-');
%             [xy,distance,t] = distance2curve(CP,ATubule,'linear');
%             for idx = 1:size(xy)
%                 pt1 = ATubule(idx,:);
%                 pt2 = xy(idx,:);
%             plot3([pt1(1) pt2(1)],[pt1(2) pt2(2)],[pt1(3) pt2(3)]);
%             plot3([CP(1,1) pt2(1)],[CP(1,2) pt2(2)],[CP(1,3) pt2(3)], '.b-');
%             end
%             mPCP = [mPCP; xy(idx,:)];
%     %         mD = [mD; distance];
%         end
%         longest = 0;
%         shortest = intmax;
%         for i = 1:size(mPCP)
%             cur = pdist([CP(1,:);mPCP(i,:)], 'euclidean');
%             if (cur > longest)
%                 longest = cur;
%             end
%             if (cur < shortest)
%                 shortest = cur;
%             end
%             mD = [mD; cur];
%         end
%         hold off;
%         figure('Name', 'Distance Graph');
%         hold on;
%         plot(mD);
%         hold off;
%         lD = [lD; (longest - shortest)];
%         delete(curPath);
%     end
%     longestD = lD;
% end