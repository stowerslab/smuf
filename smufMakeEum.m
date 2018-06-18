function [eum] = smufMakeEum( filepath, mouseNum, cageMask )
%% smufMakeEum: function to read in SMUF video frame and create End Urine Map (EUM)
% 
% This function allows selection of ROIs to exclude from analysis (i.e. stimulus marks), since we do not spectrally or
% otherwise separate tonic water or semale stimuli. It is possible to build a setup to always make these in the same place, and
% thus automate this step, but it is also useful to have flexibility in where to pipette stimuli when an animal is moving
% around. Running only this step on all videos first allows the more time-consuming urine detection across
% all frames to be done unsupervised in batches.
% 
% Author: Jason Keller
% Date: June 2018
% 
% please cite: Keller, Stowers et al, Nature Neuroscience, 2018
%
%  INPUT:
%  [filepath]
%  [mouseNum]
%
%  OUTPUT:
%  [eum] - mask with final urine spots to be counted, cropped; an overlay image is saved for simple quality control check
% 

%% read video % select negative threshold to track mouse
filename = [filepath, mouseNum, '.mp4']; %I used MP4 but any file format that VideoReader can handle is acceptable
display(['making EUM, mouse ' mouseNum]);
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
% frameRate = smufObj.FrameRate;  %in fps;
numFramesAvg = 3;  %number of frames to avg; more is better for noise, less is better to capture last bit of urine

%thresh parameters, depends on lighting and camera conditions:
% NOTE: this should be updated to the same value both here and in smufVideoToPixels.m when changing conditions
thresh = 175;  %found using threshToolMod (needs uint8 input); 

%% Preallocate the movie matrix & read frames
fullMov = zeros(vidHeight,vidWidth,3,numFramesAvg,'uint8'); %need full RGB saved to output results on top
grMov = zeros(vidHeight,vidWidth,numFramesAvg,'uint8'); %use green + red channels since they have most information
gPercent = 0.5; %red seems to counterbalance mouse shadow (light in green channel, dark in red channel)
rPercent = 1 - gPercent;

for k = 1:numFramesAvg %use last 3 frames for EUM; read one frame at a time:
    frameNum = nFrames-numFramesAvg+k;
    tempS = read(smufObj,frameNum); 
    fullMov(:,:,:,k) = tempS(topCut:bottomCut-1, leftCut:rightCut-1,:); %read method returns a H-by-W-by-B-by-F matrix, video,where H is the image frame height, W isthe image frame width, B is the number of bandsin the image (for example, 3 for RGB), and F isthe number of frames read
    grMov(:,:,k) = squeeze(fullMov(:,:,1,k)*rPercent + fullMov(:,:,2,k)*gPercent);
    grMov(:,:,k) = imsharpen(grMov(:,:,k),'Radius',20,'Amount',3, 'Threshold', 0); %sharpen to accentuate urine edges; parameters found manually
end

% manually select exclusion pixels (e.g. stimulus, mouse, reflections)
lastFrames = uint8(mean(grMov,numFramesAvg));  %use mean over time instead of max projection

%% select exclusion spots 
hFig = figure;
hIm = imagesc(lastFrames);
axis off;
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure

hDoneControl = uicontrol(hFig, 'Style', 'pushbutton', 'String', 'DONE',...  %make button to press when done
                          'Callback', @doneCallback, 'Position', [50 50 100 20]);   %#ok<NASGU>
                      
% User selects stimulus ROIs with mouse (or any other exclusion pixels):
set(hFig, 'Name','Use mouse to circle urine exclusion areas, press DONE when finished');
exitStim = false;
exclude = logical(zeros(vidHeight,vidWidth)); %#ok<LOGL>
while(~exitStim)
    hStimRoi = imfreehand(gca);
    exclude = exclude | createMask(hStimRoi,hIm);  %OR the masks to combine
end
close(hFig); 
pause on; pause(0.1); pause off; %allow figure to close

%% thresholding:
% thresh = threshToolMod(lastFrames) %can uncomment here when first setting up
BW = im2bw(lastFrames, thresh/255);
eum = BW & cageMask(topCut:bottomCut-1, leftCut:rightCut-1) &~exclude;

%% make overlay image for QC:
urineDot = uint8([255 255 0]);
overlay = zeros(vidHeight,vidWidth,'uint8');
tempR = fullMov(:,:,1,end); %pull out current frame to allow linear indexing
tempG = fullMov(:,:,2,end);
tempB = fullMov(:,:,3,end);
tempR(eum) = repmat(urineDot(1), length(find(eum)), 1); %set movie pixels for QC
tempG(eum) = repmat(urineDot(2), length(find(eum)), 1);
tempB(eum) = repmat(urineDot(3), length(find(eum)), 1);
overlay(:,:,1) = tempR;
overlay(:,:,2) = tempG;
overlay(:,:,3) = tempB;
figure; imagesc(overlay)
filenameEum = [filepath, mouseNum, '_eumOverlay.jpg'];
imwrite(overlay, filenameEum);

%% callbacks
function doneCallback(hObject, ~)
    delete(hObject);
    exitStim = true;
end

end  %end function smufMakeEum



