%Kira Hart khart@optics.arizona.edu
%July 28 2020
%This script is a test of stepper control functions

fclose('all');

clear all
close all
clc
instrreset %clear and reset any existing port communications

%%

%%Initialization
ROI =  [274 274 1500 1500];
mode = 1 ; % 1x1 binning
framespertrigger = 3;

COMmono = initializeMono();

xps = initializeMotor();

vid = initializeCamera(mode,ROI,framespertrigger);

src = getselectedsource(vid);
src.ExposureTimeControl = 'normal';
src.ExposureTime = 1; %Exposure time of Camera



