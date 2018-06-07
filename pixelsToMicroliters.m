function [microliterArray] = pixelsToMicroliters(pixelArray)
%% pixelsToMicroliters.m
% 
% This function converts urine pixels detected on chromotography paper to microliter volume, based on 2nd-order polynomial fit
%
% Author: Jason Keller
% Date: June 2018
% 
% please cite: Keller, Stowers et al, Nature Neuroscience, 2018

% these are the a/b/c coefficients from the "urineCalibration" script (please replace for a new paper/light/camera setup)
% p = [0.0000089602 0.0254 0]; %for THICK paper from 'urineCalibration' script, but with y-intecept forced from -0.4450 to 0, since you can't have negative pixels
p = [0.0000015792 0.0179 0]; %for THIN paper

microliterArray = polyval(p, pixelArray);

end %end function