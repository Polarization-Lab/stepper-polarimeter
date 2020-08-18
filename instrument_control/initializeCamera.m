function [vid] = initializeCamera(mode,ROIPosition,FramesPerTrigger)
%INITIALIZECAMERA Summary of this function goes here
%Kira Hart khart@optics.arizona.edu
%July 28 2020
%This script is designed to initialize the camera
%note for R2020 and beyond an additional add on is neededed
%see manual
%Bin size can be 1,2 or 4
% ROI POsition is [xstart,ystart,w,h]
%FramerPerTrigger is an interger

disp('Initializing Camera')
if mode == 1
    vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_FastMode'); %vid source for camera, 1x1 BIN
elseif mode == 2
    vid = videoinput('hamamatsu', 1, 'MONO16_BIN2x2_1024x1024_FastMode'); %2x2bin
elseif mode == 4
    vid = videoinput('hamamatsu', 1, 'MONO16_BIN4x4_512x512_FastMode'); %4x4bin
else
    disp('Invalid bin size, enter 1,2, or 4')
end

vid.FramesPerTrigger = FramesPerTrigger; 
vid.ROIPosition = ROIPosition;

end

