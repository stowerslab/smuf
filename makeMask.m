%% makeMask.m
% This script creates a mask from a single frame of a SMUF behavioral video,
% such that urinated pixels will only be found within the mask. It is helpful to
% exclude arena borders where reflections may cause noise.
% 
% Author: Jason Keller
% Date: June 2018
% 
% please cite: Keller, Stowers et al, Nature Neuroscience, 2018

cd 'C:\data\'
clear all; close all;
I = imread('makeMask.jpg'); %use JPEG created from single frame of behavior video
imshow(I)

hP = impoly();  %create mask manually by clicking mouse into a polygon, excluding cage walls / water spout, etc.
cageMask = hP.createMask(); 
save('cageMask.mat', 'cageMask') %save to MAT file for later use