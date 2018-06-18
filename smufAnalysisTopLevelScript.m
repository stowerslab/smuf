%% top-level script to process video data for scent-mark-to-female-urine (SMUF) behavior
% NOTE: requires Statistics and Image Processing toolboxes
% 
% This is the top level script to analyze urine marking. It calls "smufMakeEum" and "smufVideoToPixels"
% functions, groups data,, calculates statistics, and makes various plots. Steps are run independently
% by setting PROGRAM CONTROL options below, as explained in the documentation and comments.
% 
% Author: Jason Keller
% Date: June 2018
% 
% please cite: Keller, Stowers et al, Nature Neuroscience, 2018
%
%   INPUT: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [filepath]
%   [mouseNum]
%   [eum] - using 'smufMakeEum.m' - uses adaptive threshold to make map of urine pixels to search here
%   [saveBwMovie]
%

clear all; close all;
saveBwMovie = true;
 filepathRoot = 'C:\data\';

load 'C:\data\cageMask.mat' %created using makeMask.m
fontSz = 20;

%% PROGRAM CONTROL:
findROIs = 1; % first: manual step; do this only
findUrine = 0; % second: automated step; run and then QC output videos
combineMice = 0; % third, make groups to compare as necessary
readEthovision = 0; % optionally read input from Ethovision
plotRaster = 0; % optionally plot raster of urine spots detected
plotPixelsPerSpot = 0; % optionally plot pixels per spot
computeStats = 0; % fourth, compute stats for group
plotData = 0; % fifth, plot data - different plots are manually controlled by uncommenting them

%%% INPUT filenames for individual mice/videos (to access later, ex. 2nd entry is mouseNums{1,2}):
mouseNums = {'mouse1' 'mouse2'}; %all videos for TEST DAY 1
totalMice = size(mouseNums, 2);

% path to save summary data (for combining data across mice):
matNameSaveSummary = [filepathRoot, 'allSmufData_testDay1.mat'];

% for computing stats and plotting
groupsToCompare = {'testDay1', 'testDay2'};  % groups to compare (assumes "allSmufData_" is prefix)
    
stat1 = 'pixelsUrinatedPerFrame';
stat2 = 'pixelsUrinatedPerFrameNorm';  %normalized to maxium to facilitate comparing different mice/days
stat3 = 'spotsUrinatedPerFrame'; %tfor last frame this is totalSpots
stat4 = 'spotsUrinatedPerFrameNorm';
stat5 = 'spotsUrinatedPerFrame1stHalf';
stat6 = 'spotsUrinatedPerFrame2ndHalf';
stat7 = 'spotsUrinatedPerFrameNormMultDays'; %normalized for each mouse for all days comparing
stat8 = 'pixelsPerSpotPerFrame'; %may be useful looking at background vs. SMUF
stat9 = 'nosepokesStim'; %binary array over frames from Ethovision - 1 if nose in stimulus zone
stat10 = 'nosepokesTonic'; %binary array over frames from Ethovision - 1 if nose in tonic zone
stat11 = 'velocity'; % in cm/s, velocity of animal center point per frame, from Ethovision
stat12 = 'frameTimes';
stat13 = 'smufIndex'; % not used yet

smufStats = struct(stat1,[], stat2,[], stat3,[], stat4,[], stat5,[], stat6,[], stat7,[], stat8,[], stat9,[], stat10,[], stat11,[], stat12,[], stat13,[]);
numFramesAnalysis = 3601; %number of frames to analyze in videos
halfFramesAnalysis = (numFramesAnalysis-1)/2; %half number of frames to analyze in videos (should be integer) - corresponds to when stimulus is applied
frameRate = 1/15; %1/fps

%% do all manual steps for each mouse/video first (i.e. ROIs and thresholds if necessary):
if findROIs
    for k = 1:totalMice
        mouseNum = mouseNums{1,k};
        eum = smufMakeEum(filepathRoot, mouseNum, cageMask); %EUM = end urine map
        eumName = [filepathRoot, mouseNum, '_eum.mat'];
        save(eumName, 'eum'); %save data for each individual mouse
        clear eum;
    end
end

if findUrine
    % now do heavy lifting without supervision by calling "smufVideoToPixels":
    for k = 1:totalMice %#ok<*UNRCH>
        mouseNum = mouseNums{1,k};
        eumName = [filepathRoot, mouseNum, '_eum.mat'];
        load(eumName); 
        [pixelsUrinatedPerFrame, pixelsUrinatedPerFrameNorm, spotsUrinatedPerFrame, spotsUrinatedPerFrameNorm, pixelsPerSpotPerFrame, frameTimes, smufIndex] = smufVideoToPixels(filepathRoot, mouseNum, eum, saveBwMovie);
        matName = [filepathRoot, mouseNum, '_smufData.mat']; %save data for each individual mouse
        save(matName, 'pixelsUrinatedPerFrame', 'pixelsUrinatedPerFrameNorm', 'spotsUrinatedPerFrame', 'spotsUrinatedPerFrameNorm', 'pixelsPerSpotPerFrame', 'frameTimes', 'smufIndex');
    end
end

%% combine data across mice & compute pixelDiff / spotDiff raster:
if combineMice
    pixelsDiff = false(totalMice,numFramesAnalysis-1); %array of 0s and 1s for thresholded "urination events"
    spotsDiff = false(totalMice,numFramesAnalysis-1); %array of 0s and 1s for spots detected
    nosepokesStimAll = false(totalMice,numFramesAnalysis-1); %array of 0s and 1s for spots detected

    for k = 1:totalMice
        mouseNum = mouseNums{1,k};
        matName = [filepathRoot, mouseNum, '_smufData.mat'];
        load(matName);

        smufStats.pixelsUrinatedPerFrame(end+1,:) = pixelsUrinatedPerFrame; 
        smufStats.pixelsUrinatedPerFrameNorm(end+1,:) = pixelsUrinatedPerFrameNorm;
        smufStats.spotsUrinatedPerFrame(end+1,:) = spotsUrinatedPerFrame;
        smufStats.spotsUrinatedPerFrameNorm(end+1,:) = spotsUrinatedPerFrameNorm;
        smufStats.pixelsPerSpotPerFrame(end+1,:) = pixelsPerSpotPerFrame;
        smufStats.frameTimes(end+1,:) = frameTimes;
        smufStats.smufIndex(end+1,:) = 0;

        tempPixDiff = diff(pixelsUrinatedPerFrame);
        urineSpotThresh = 10; %set all "urination events" above thresh to '1' to make a raster (alternatively could use log scale):
        tempPixDiff(find(tempPixDiff<urineSpotThresh)) = false; %#ok<FNDSB>
        tempPixDiff(find(tempPixDiff>=urineSpotThresh)) = true; %#ok<FNDSB>
        pixelsDiff(k,:) = tempPixDiff;

        tempSpotDiff = diff(spotsUrinatedPerFrame);
        spotThresh = 1; %set all "spot events" above thresh to '1' to make a raster (alternatively could use log scale):
        tempSpotDiff(find(tempSpotDiff<spotThresh)) = false; %#ok<FNDSB>
        tempSpotDiff(find(tempSpotDiff>=spotThresh)) = true; %#ok<FNDSB>
        spotsDiff(k,:) = tempSpotDiff;   
        
        %read nosepoke data from Ethovision if necessary:
        if readEthovision
            [smufStats.nosepokesStim(end+1,:), smufStats.nosepokesTonic(end+1,:), smufStats.velocity(end+1,:)] = readNosepokeVelocity(filepathRoot, mouseNum, numFramesAnalysis);
            nosepokesStimAll(k,:) = smufStats.nosepokesStim(k,1:end-1);
        end
    end

    if plotRaster
%         hPixelRaster = figure; % plot raster of "urination events"
%         LineFormatVert.LineWidth = 3;
%         plotSpikeRaster(pixelsDiff, frameRate/60, 'PlotType','vertline',...
%                                    'VertSpikePosition',0,...
%                                    'VertSpikeHeight',0.8);
%         xlabel('time (min)', 'FontSize', fontSz);
%         ylabel('mice', 'FontSize', fontSz);
%         set(gca,'YTick',[]);
%         yLims = get(gca, 'YLim');
%         x = [2, 2];
%         y = [yLims(1) yLims(2)];
%         hold on;
%         plot(x, y, 'Color', [0 0.75 0.5], 'LineStyle', ':', 'LineWidth', 2); %dotted line for stimulus
%         hold off;
%         % plot2svg([filepathRoot, 'allSmufPixelsRaster.svg'], hPixelRaster);
        
        hSpotRaster = figure; % plot raster of "spot events"
        LineFormatSpots.LineWidth = 1;
        LineFormatSpots.Color = [0 0 0]; %mark urine spots with vertical black lines
        plotSpikeRaster(spotsDiff, frameRate/60, 'PlotType','vertline',...
                                   'LineFormat',LineFormatSpots,...
                                   'VertSpikePosition',0,...
                                   'VertSpikeHeight',0.8);
        hold on;
        LineFormatNose.Color = [1 0 0]; %mark nosepokes with horizontal red lines
        LineFormatNose.LineWidth = 2;
        plotSpikeRaster(nosepokesStimAll, frameRate/60, 'PlotType','vertline',...
                                   'LineFormat',LineFormatNose,...
                                   'VertSpikePosition',0,...
                                   'VertSpikeHeight',0.2);
        xlabel('time (min)', 'FontSize', fontSz);
        ylabel('mice', 'FontSize', fontSz);
        set(gca,'YTick',[]);
        yLims = get(gca, 'YLim');
        x = [2, 2];
        y = [yLims(1) yLims(2)];
        hold on;
        plot(x, y, 'Color', [0 0.75 0.5], 'LineStyle', ':', 'LineWidth', 2); %dotted line for stimulus
        hold off;
%         plot2svg([filepathRoot, 'spotRaster_', matNameSaveSummary, '.svg'], hSpotRaster);
    end
    % 
    if plotPixelsPerSpot %plot pixels per spot if desired
        hPixelsPerSpot = figure;
        hold on;
        for k = 1:totalMice
           plot(smufStats.frameTimes(k,:), smufStats.pixelsPerSpotPerFrame(k,:), 'o-');
        end
        xlabel('time (min)', 'FontSize', fontSz);
        ylabel('pixels per spot', 'FontSize', fontSz);
        hold off;
    end
    
    save(matNameSaveSummary, 'smufStats', 'pixelsDiff', 'spotsDiff');
end

%% Compute stats
if computeStats
    numGroupsToCompare = size(groupsToCompare, 2);
    smufStatsVec = repmat(smufStats, numGroupsToCompare, 1);

    meanPix = zeros(numFramesAnalysis, numGroupsToCompare);
    stdPix = zeros(numFramesAnalysis, numGroupsToCompare);
    semPix = zeros(numFramesAnalysis, numGroupsToCompare);
    meanciPix = zeros(numFramesAnalysis, 2, numGroupsToCompare);

    meanPixNorm = zeros(numFramesAnalysis, numGroupsToCompare);
    stdPixNorm = zeros(numFramesAnalysis, numGroupsToCompare);
    semPixNorm = zeros(numFramesAnalysis, numGroupsToCompare);
    meanciPixNorm = zeros(numFramesAnalysis, 2, numGroupsToCompare);

    meanSpots = zeros(numFramesAnalysis, numGroupsToCompare);
    stdSpots = zeros(numFramesAnalysis, numGroupsToCompare);
    semSpots = zeros(numFramesAnalysis, numGroupsToCompare);
    maxSpots = zeros(totalMice, 1);

    meanSpotsNorm = zeros(numFramesAnalysis, numGroupsToCompare);
    stdSpotsNorm = zeros(numFramesAnalysis, numGroupsToCompare);
    semSpotsNorm = zeros(numFramesAnalysis, numGroupsToCompare);
    
    meanSpotsNormMult = zeros(numFramesAnalysis, numGroupsToCompare);
    stdSpotsNormMult = zeros(numFramesAnalysis, numGroupsToCompare);
    semSpotsNormMult = zeros(numFramesAnalysis, numGroupsToCompare);
    
    meanSpots1stHalf = zeros(totalMice, 1);
    stdSpots1stHalf = zeros(totalMice, 1);
    semSpots1stHalf = zeros(totalMice, 1);
    
    meanSpots2ndHalf = zeros(totalMice, 1);
    stdSpots2ndHalf = zeros(totalMice, 1);
    semSpots2ndHalf = zeros(totalMice, 1);

    for i = 1:numGroupsToCompare
        load([filepathRoot, 'allSmufData_', groupsToCompare{1,i}, '.mat'])
        smufStatsVec(i) = smufStats; %create vector of structs for comparisons
        tempMaxSpots(:,1) = max(smufStats.spotsUrinatedPerFrame,[],2);  %#ok<SAGROW>
        maxSpots(:,1) = max(tempMaxSpots, maxSpots); % calculate to normalize by entire multi-day max # spots
      
        % calculate statistics:
        [meanPix(:,i), stdPix(:,i), semPix(:,i), meanciPix(:,:,i)] = grpstats(smufStats.pixelsUrinatedPerFrame,[],{'mean','std','sem','meanci'});
        [meanPixNorm(:,i), stdPixNorm(:,i), semPixNorm(:,i), meanciPixNorm(:,:,i)] = grpstats(smufStats.pixelsUrinatedPerFrameNorm,[],{'mean','std','sem','meanci'});
        [meanSpots(:,i),stdSpots(:,i),semSpots(:,i)] = grpstats(smufStats.spotsUrinatedPerFrame,[],{'mean','std','sem'});
        [meanSpotsNorm(:,i),stdSpotsNorm(:,i),semSpotsNorm(:,i)] = grpstats(smufStats.spotsUrinatedPerFrameNorm,[],{'mean','std','sem'}); 
        [meanSpots1stHalf(:,i),stdSpots1stHalf(:,i),semSpots1stHalf(:,i)] = grpstats(smufStatsVec(i).spotsUrinatedPerFrame1stHalf,[],{'mean','std','sem'});
        [meanSpots2ndHalf(:,i),stdSpots2ndHalf(:,i),semSpots2ndHalf(:,i)] = grpstats(smufStatsVec(i).spotsUrinatedPerFrame2ndHalf,[],{'mean','std','sem'});
    end
    
    for i = 1:numGroupsToCompare  %go back through one more time if normalizing by grand total of spots over multiple days
        for k = 1:totalMice
            smufStatsVec(i).spotsUrinatedPerFrameNormMultDays(k,:) = smufStatsVec(i).spotsUrinatedPerFrame(k,:)./maxSpots(k);
        end
        FriedData(:,i) = smufStatsVec(i).spotsUrinatedPerFrameNormMultDays(:,end); %save matrix for Friedman test later
        [meanSpotsNormMult(:,i),stdSpotsNormMult(:,i),semSpotsNormMult(:,i)] = grpstats(smufStatsVec(i).spotsUrinatedPerFrameNormMultDays,[],{'mean','std','sem'});
    end
    
end

%% plot data:
if plotData
    cmap = colormap('lines');
    close(1); % for some reason colormap() creates a new figure if one not already open, and this can't be suppressed unless axis object already exists

    % %plot urine pixels over time, NOT normalized %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     hPix = figure;
%     hold on;
%     title('SMUF behavior')
%     xlabel('time (min)', 'FontSize', 16);
%     ylabel('pixels with urine', 'FontSize', 16);
%     for i = 1:numGroupsToCompare
%         [hlA,hpA] = boundedline(smufStatsVec(i).frameTimes(1,:), meanPix(:,i), semPix(:,i), 'cmap', [cmap(i,:)],'alpha','y');
%     end
%     % meanSlope = meanPix(halfFramesAnalysis,groupNum)/halfFramesAnalysis;
%     % meanSlopeTrace = [1:1:numFramesAnalysis]*meanSlope; %plot continuation of mean slope from first 2 min
%     % plot(smufStatsVec(groupNum).frameTimes(1,:), meanSlopeTrace, 'k--')
%     axis tight;
%     set(gca, 'FontSize', 16);
%     x = [smufStatsVec(1).frameTimes(1,halfFramesAnalysis), smufStatsVec(1).frameTimes(1,halfFramesAnalysis)];
%     yLims = get(gca, 'YLim');
%     y = [yLims(1) yLims(2)];
%     plot(x, y, 'Color', [0 0.75 0.5], 'LineStyle', ':', 'LineWidth', 2); %dotted line for stimulus
%     hold off;
    % colormap(cmap(1:numGroupsToCompare,:))
    % cbar_axes = colorbar('location','EastOutside','YTick',[0.1 0.9],'YTickLabel',{[' Day 1'];[' Day 4']});
    
    % %plot urine pixels over time, NORMALIZED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hPixNorm = figure;
    % hold on;
    % title('SMUF behavior, normalized')
    % xlabel('time (min)', 'FontSize', 16);
    % ylabel('pixels with urine, normalized', 'FontSize', 16);
    % for i = 1:numGroupsToCompare
    %     [hlA,hpA] = boundedline(smufStatsVec(i).frameTimes(1,:), meanPixNorm(:,i), semPixNorm(:,i), 'cmap', [cmap(i,:)],'alpha','y');
    % end
    % axis tight;
    % set(gca, 'FontSize', 16);
    % x = [smufStatsVec(1).frameTimes(1,halfFramesAnalysis), smufStatsVec(1).frameTimes(1,halfFramesAnalysis)];
    % yLims = get(gca, 'YLim');
    % y = [yLims(1) yLims(2)];
    % plot(x, y, 'Color', [0 0.75 0.5], 'LineStyle', ':', 'LineWidth', 2); %dotted line for stimulus
    % hold off;

% %     %plot urine SPOTS over time, NOT normalized %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     hSpots = figure;
%     hold on;
%     xlabel('time (min)', 'FontSize', fontSz);
%     ylabel('# urine spots', 'FontSize', fontSz);
%     for i = 1:numGroupsToCompare
%         [hlA,hpA] = boundedline(smufStatsVec(i).frameTimes(1,:), meanSpots(:,i), semSpots(:,i), 'cmap', [cmap(i,:)],'alpha','y');
%     end
%     axis tight;
% %     axis([0 4.1 0 1.1]);
%     set(gca, 'FontSize', fontSz);
%     x = [smufStatsVec(1).frameTimes(1,halfFramesAnalysis), smufStatsVec(1).frameTimes(1,halfFramesAnalysis)];
%     yLims = get(gca, 'YLim');
%     y = [yLims(1) yLims(2)];
%     plot(x, y, 'Color', [0 0.75 0.5], 'LineStyle', ':', 'LineWidth', 2); %dotted line for stimulus
%     hold off;

    %plot urine SPOTS over time, NORMALIZED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hSpotsNorm = figure;
    hold on;
    xlabel('time (min)', 'FontSize', fontSz);
    ylabel('% max. urine spots', 'FontSize', fontSz);
    for i = 1:numGroupsToCompare
        [hlA,hpA] = boundedline(smufStatsVec(i).frameTimes(1,:), meanSpotsNormMult(:,i), semSpotsNormMult(:,i), 'cmap', [cmap(i,:)],'alpha','y');
    end
    axis tight;
    axis([0 4.1 0 1.1]);
    set(gca, 'FontSize', fontSz);
    x = [smufStatsVec(1).frameTimes(1,halfFramesAnalysis), smufStatsVec(1).frameTimes(1,halfFramesAnalysis)];
    yLims = get(gca, 'YLim');
    y = [yLims(1) yLims(2)];
    plot(x, y, 'Color', [0 0.75 0.5], 'LineStyle', ':', 'LineWidth', 2); %dotted line for stimulus
    hold off;

    % %plot individual mouse traces over time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % hFigRawUrineTraces = figure;
    % hold on;
    % cmap = colormap('lines');
    % lineStyles = {'-','-','-.','-.'};
    % legendLabels = {'Balb Male','B6 Male','Female','TrpC2 KO'};
    % plotHandles = zeros(numGroupsToCompare,1);
    % for i = 1:numGroupsToCompare    
    %     for j = 1:size(smufStatsVec(i).pixelsUrinatedPerFrame,1)
    %         plotHandles(i) = plot(smufStatsVec(i).frameTimes(1,:), smufStatsVec(i).pixelsUrinatedPerFrame(j,:), 'Color', cmap(i,:), 'LineWidth', 2, 'LineStyle', lineStyles{i});
    %     end
    % end
    % axis tight; 
    % legend(plotHandles, legendLabels)
    % set(gca, 'FontSize', 16);
    % x = [smufStatsVec(1).frameTimes(1,halfFramesAnalysis), smufStatsVec(1).frameTimes(1,halfFramesAnalysis)];
    % yLims = get(gca, 'YLim');
    % y = [yLims(1) yLims(2)];
    % plot(x, y, 'Color', [0 0.75 0.5], 'LineStyle', ':', 'LineWidth', 2); %dotted line for stimulus
    % xlabel('time (min)', 'FontSize', 16);
    % ylabel('pixels with urine', 'FontSize', 16);

%     % plot total # spots 0-2 min & 2-4 min (use 1st 'groupToCompare') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     hSpotsBarChart = figure;
%     hold on;
%     errorbar([1 2], [meanSpots1stHalf(1,1) meanSpots2ndHalf(1,1)], [semSpots1stHalf(1,1) semSpots2ndHalf(1,1)], 'rx', 'LineWidth', 2);
%     plot([1 2], [meanSpots1stHalf(1,1) meanSpots2ndHalf(1,1)], '*-', 'Color', 'k', 'MarkerSize', 5, 'LineWidth', 2);
%     hold off;
%     axis([0.5 2.5 0 meanSpots2ndHalf(1,1)+stdSpots2ndHalf(1,1)+0.2]);
%     set(gca, 'FontSize', fontSz);
%     set(gca,'XTick', [1 2], 'XTickLabel',{'before stim.';'after stim.'})
%     ylabel('# urine spots', 'FontSize', fontSz)
% %     ksResult1 = kstest(smufStatsVec(1).spotsUrinatedPerFrame1stHalf);
%     p2Min4MinWilcox = signrank(smufStatsVec(1).spotsUrinatedPerFrame1stHalf,smufStatsVec(1).spotsUrinatedPerFrame2ndHalf)


%     plot mean # spots (normalized, w/ SEM) over 4-day DREADD treatment, with overlaid individual traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     hDreadd = figure;
%     hold on;
%     for ip = 1:size(smufStatsVec(1).spotsUrinatedPerFrameNormMultDays,1)  %first plot individual mouse traces
%         indDreaddTrace = [smufStatsVec(1).spotsUrinatedPerFrameNormMultDays(ip,end) smufStatsVec(2).spotsUrinatedPerFrameNormMultDays(ip,end) smufStatsVec(3).spotsUrinatedPerFrameNormMultDays(ip,end) smufStatsVec(4).spotsUrinatedPerFrameNormMultDays(ip,end)];
%         plot(indDreaddTrace, 'Color', cmap(ip,:), 'LineWidth', 2.5, 'LineStyle', '--')
%     end
%     errorbar([1 2 3 4], [meanSpotsNormMult(end,1) meanSpotsNormMult(end,2) meanSpotsNormMult(end,3) meanSpotsNormMult(end,4)], [semSpotsNormMult(end,1) semSpotsNormMult(end,2) semSpotsNormMult(end,3) semSpotsNormMult(end,4)], 'Color', 'k', 'LineWidth', 4);
%     plot([1 2 3 4], [meanSpotsNormMult(end,1) meanSpotsNormMult(end,2) meanSpotsNormMult(end,3) meanSpotsNormMult(end,4)], 'o-', 'Color', 'k', 'MarkerSize', 6, 'LineWidth', 8);
%     hold off;
%     axis([0.5 4.5 0 max(max(meanSpotsNormMult))+max(max(semSpotsNormMult))+0.2]);
%     set(gca, 'FontSize', fontSz);
%     set(gca,'XTick', [1 2 3 4], 'XTickLabel',{'CNO';'saline';'CNO';'saline'})
%     ylabel('% max. urine spots', 'FontSize', fontSz)
%     [pFriedman,tblFriedman,statsFriedman] = friedman(FriedData,1);
%     cFriedman = multcompare(statsFriedman, 'estimate', 'friedman', 'ctype', 'dunn-sidak')

    %save figure to raster & vector format (NOTE plot2svg is from MATLAB central):
    % plot2svg([filepathRoot, 'allSmufPixelsOverTimeNorm.svg'], hSpotsNorm)
    % plot2svg([filepathRoot, '\plots\allSmufSpotsSem_CnoSaline.svg'], hDreadd);
end

