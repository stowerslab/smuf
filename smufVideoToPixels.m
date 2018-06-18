function [ pixelsUrinatedPerFrame, pixelsUrinatedPerFrameNorm, spotsUrinatedPerFrame, spotsUrinatedPerFrameNorm, pixelsPerSpotPerFrame, frameTimes, smufIndex] = smufVideoToPixels( filepath, mouseNum, eum, saveBwMovie )
%% smufVideoToPixels.m
% 
% This is the main function that takes a video and EUM and searches for new urine spots in each frame,
% using simple thresholding (after sharpen filter preprocessing) and connected components.
% 
% Author: Jason Keller
% Date: June 2018
% 
% please cite: Keller, Stowers et al, Nature Neuroscience, 2018
%
%   INPUT: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [filepath] - for contructing input video filename
%   [mouseNum] - for contructing input video filename
%   [eum] - using 'smufMakeEum.m' - uses adaptive threshold to make map of urine pixels to search here
%   [saveBwMovie] - whether or not to save threhsolded overlay movie for quality control
%
%   OUTPUT: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [pixelsUrinatedPerFrame] - 1D array (#frames) with urinated pixels detected in that frame
% [pixelsUrinatedPerFrameNorm] - 1D array (#frames) with normalized urinated pixels detected in that frame
% [spotsUrinatedPerFrame] - 1D array (#frames) with urinated spots detected in that frame
% [spotsUrinatedPerFrameNorm] - 1D array (#frames) with normalized urinated spots detected in that frame
% [pixelsPerSpotPerFrame] - 1D array (#frames) with urinated spots detected in that frame
% [frameTimes] - 1D array (#frames) with times based on frame rate
% [smufIndex] - not used for now
%    
% KNOWN ISSUES: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   (1) can never separate urinated pixels that overlap stimulus or each other with current method

%% read video & track mouse
filename = [filepath, mouseNum, '.mp4'];
display(['mouse ' mouseNum]);
smufObj = VideoReader(filename);

nFrames = smufObj.NumberOfFrames;
% NOTE: parameters below must be updated in both "smufMakeEum" and "smufVideoToPixels"
topCut = 5; %hard code to cut off black border for cropped video output (use makeMask.jpg to estimate)
bottomCut = 350;
leftCut = 24;
rightCut = 590;

vidHeight = bottomCut-topCut; %smufObj.Height;  
vidWidth = rightCut-leftCut; %smufObj.Width; 
% nTotalPixels = vidHeight*vidWidth;
frameRate = smufObj.FrameRate;  %in fps

%Preallocate the movie matrix
fullMov = zeros(vidHeight,vidWidth,3,nFrames,'uint8'); %need full RGB saved to output results on top
grMov = zeros(vidHeight,vidWidth,nFrames,'uint8'); %use green + red channels since they have most information
gPercent = 0.5; %red seems to counterbalance mouse shadow (light in green channel, dark in red channel)
rPercent = 1 - gPercent;
BW = false(vidHeight,vidWidth,nFrames); %this is a logical array for a binary image used to overlay urine detection results

urineDot = uint8([255 255 0]);  % for tracking urinated pixels later; use yellow for pee
%thresh parameters:
% NOTE: this should be updated to the same value both here and in smufMakeEum.m when changing conditions
thresh = 175;  %found using threshToolMod (see smufMakeEum function)

hWait = waitbar(1/nFrames,['Reading frames & thresholding ', mouseNum]); %set up progress bar
for k = 1:nFrames %Read one frame at a time:
    waitbar(k/nFrames, hWait); %update progress bar
    temp = read(smufObj,k); 
    fullMov(:,:,:,k) = temp(topCut:bottomCut-1, leftCut:rightCut-1,:); %read method returns a H-by-W-by-B-by-F matrix, video,where H is the image frame height, W isthe image frame width, B is the number of bandsin the image (for example, 3 for RGB), and F isthe number of frames read
    grMov(:,:,k) = squeeze(fullMov(:,:,1,k)*rPercent + fullMov(:,:,2,k)*gPercent);
    grMov(:,:,k) = imsharpen(grMov(:,:,k),'Radius',20,'Amount',3, 'Threshold', 0); %sharpen to accentuate urine edges; parameters found manually
    
    %THRESHOLD FOR URINE
    % adaptive thresholding not needed if sharpening image, but may be useful for uneven lighting conditions
    BW(:,:,k) = im2bw(grMov(:,:,k), thresh/255);
end

close(hWait);
clear temp grMov; %free up some memory

if ~saveBwMovie
    clear fullMov;  %free up some memory
% else %make thresholded movies as a sanity check for occlusions, noise, etc.
%     gBW = uint8(zeros(size(BW,1),size(BW,2),1,size(BW,3)));  %need to add empty dimension for videoWrite function
%     gBW(:,:,1,:) = uint8(BW);
%     a = find(BW);
%     gBW(a) = 2^8-1; %#ok<FNDSB>  %scale to max for writing to video
end

%% find urinated pixels & spots
% check for overlap of thresholded spots by spotOverlap%
% end product is to make arrays of [# pixels urinated, frames] and [# spots, frames] to show "urination events"
hWait = waitbar(1/nFrames,'Finding urinated pixels...'); %set up progress bar
pixelsUrinatedPerFrame = zeros(nFrames,1);
spotsUrinatedPerFrame = zeros(nFrames,1);
tempEum = eum; % for running spot count
searchPixels = eum; % for running pixel count
searchSpots = bwlabel(eum);  %potential urine spots are connected compenents of EUM; could use 'bwconncomp' and 'regionprops' here too
% figure; imshow(label2rgb(searchSpots));
totalSpots = max(max(searchSpots));
currentSpots = 0; %assume zero spots to start with
pixelsPerSpotPerFrame = zeros(nFrames,1);
totalUrinePixels = length(find(eum));
allSpotsUrinated = false(size(eum));  %keep running count for movie output (these should grow as searchSpots declines)
newPixelsUrinated = false(size(eum));
allPixelsUrinated = false(size(eum));  %keep running count for movie output (these should grow as searchPixels declines)

spotOverlap = 0.5; %percentage overlap needed before marking a spot as filled; rather arbitrary tradeoff between speed of detection and robustness

for k = 1 : nFrames
    currentFrame = BW(:,:,k);
    currentLabelFrame = bwlabel(BW(:,:,k)); %connected components
    %%% PIXELS:
    newPixelsUrinated(searchPixels) = currentFrame(searchPixels); %these are indeces within searchPixels that are above threshold
    searchPixels = searchPixels & ~newPixelsUrinated;  %remove urine pixels found from the search
    allPixelsUrinated = allPixelsUrinated | newPixelsUrinated;
    %%% SPOTS:
    pixelsPerSpotPerFrame(k) = 0;
    numSpotsRemaining = max(max(searchSpots));
    for spotK = 1:numSpotsRemaining %loop over all unfound spots
        spotIdxList = find(searchSpots==spotK); %total indeces for current spot being searched
        potentialSpotIdx = intersect(find(newPixelsUrinated), spotIdxList);  %overlap within current spot being searched
        labelsInSpot = unique(nonzeros(currentLabelFrame(potentialSpotIdx))); % find BWgl connected components inside current spot (s/b only 1 normally)
        labelK = 1;
        while labelK <= length(labelsInSpot) % now loop over those labels
            labelTotalIdxList = find(currentLabelFrame==labelsInSpot(labelK));
            labelOverlap = intersect(labelTotalIdxList, spotIdxList);
            if (length(labelOverlap) > spotOverlap*length(spotIdxList)) %if overlap is greater than spotOverlap%, consider the spot marked
                tempEum(spotIdxList) = 0; %delete the spot from the searchable area
                allSpotsUrinated(spotIdxList) = 1;
                searchSpots = bwlabel(tempEum); %redo connected components once a deletion happens
                currentSpots = totalSpots - max(max(searchSpots));
                pixelsPerSpotPerFrame(k) = length(spotIdxList);  %record total number of pixels in the spot found for this frame
                labelK = length(labelsInSpot) + 1; %exit while loop
            end
            labelK = labelK + 1;
        end
    end
    
    if saveBwMovie % turn marked pixels/spots yellow in current frame for quality control
        % if k > 100 % uncomment here if you want to turn on yellow pixels only after a delay
        tempR = fullMov(:,:,1,k); %pull out current frame to allow linear indexing
        tempG = fullMov(:,:,2,k);
        tempB = fullMov(:,:,3,k);
        tempR(allPixelsUrinated) = repmat(urineDot(1), length(find(allPixelsUrinated)), 1); %set movie pixels for QC (PIXELS)
        tempG(allPixelsUrinated) = repmat(urineDot(2), length(find(allPixelsUrinated)), 1);
        tempB(allPixelsUrinated) = repmat(urineDot(3), length(find(allPixelsUrinated)), 1);
%         tempR(allSpotsUrinated) = repmat(urineDot(1), length(find(allSpotsUrinated)), 1); %set movie pixels for QC (SPOTS)
%         tempG(allSpotsUrinated) = repmat(urineDot(2), length(find(allSpotsUrinated)), 1);
%         tempB(allSpotsUrinated) = repmat(urineDot(3), length(find(allSpotsUrinated)), 1);
        fullMov(:,:,1,k) = tempR;
        fullMov(:,:,2,k) = tempG;
        fullMov(:,:,3,k) = tempB;
        % end
    end
    
    pixelsUrinatedPerFrame(k) = totalUrinePixels - length(find(searchPixels));  %update running counts
    spotsUrinatedPerFrame(k) = currentSpots;
    waitbar(k/nFrames, hWait); %update progress bar
end
close(hWait);
%normalize pixel counts to max # pixels for each mouse, thus isolating timescale of response to compare across mice, whereas total # spots and total # pixels can give variability between mice
pixelsUrinatedPerFrameNorm = pixelsUrinatedPerFrame ./ max(pixelsUrinatedPerFrame);
spotsUrinatedPerFrameNorm = spotsUrinatedPerFrame ./ max(spotsUrinatedPerFrame);
smufIndex = 0; %not used for now, so just set to zero

frameTimes = 0:1/frameRate:1/frameRate*(nFrames-1); %in seconds
frameTimes = frameTimes/60; %in min

%% save movie for quality control
if saveBwMovie
%     fullMov = fullMov(:,:,:,800:7200);  %cut a certain segment by frame #'s if necessary
%     saveMpg4(filepath, mouseNum, gBW, '_ThreshOut')  %for saving thresholding
    saveMpg4(filepath, mouseNum, fullMov, '_YellowOut')  %for saving final urination overlay movie
end

%% other functions
function saveMpg4(filepath, mouseNum, movMatrix, addString)
    display(['Saving ' addString ' movie for ', mouseNum, '...'])
    bwSmufObj = VideoWriter([filepath, mouseNum, addString, '.mp4'], 'MPEG-4');
    bwSmufObj.FrameRate = 90;  %speed up for easier viewing
    open(bwSmufObj)
    writeVideo(bwSmufObj,movMatrix)
    close(bwSmufObj)
end

% keyboard % uncomment to examine variables before returning from function - type 'return' to exit 

end  %end function smufVideoToPixels



