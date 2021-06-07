% Stepper Code
% Authors: James Heath
% Last Updated: 2021/01/01

% Description:
% This code performs calibration and produces sample results from the
% calculated calibration values. It can perform both Spectral and/or
% Amplitude calibrations before running the analysis

% Inputs: None
% Outputs: Saved files for calibration values, Mueller matrices, calculated
% diattenuation, linear diattenuation, polarizance, and depolarization index 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear variables
clear all; close all;clc;

%%
nSteps = 64; %StepperCal = 64, RGB950 = 40
LambdaList = linspace(450,750,31);
nLambda = length(LambdaList);

%image size. Needs to be more robust.
xBegROI = 746;
xEndROI = 755;
yBegROI = 746;
yEndROI = 755;

%Set Plots to 1 for plots to pop up and 0 for no plots
Plots = 1;

%Calibration file path
calibrationFilePath = 'D:\Measurements\Air_Calibrations\New_Diffuser\AirCalibration-06-Dec-2020.h5';
sampleFilePath = 'D:\Measurements\Air_Calibrations\New_Diffuser\Dichroic45-07-Dec-2020.h5';

%% Load air measurements

airMeasurement = LoadAirMeasurements(calibrationFilePath,LambdaList,nSteps);

%% Spectral Calibration

%Grab calibration data 
%Change '~' to 'W' to obtain W-matrix for spectral analysis
[caldata,~] = MM_SpectralCalibration(airMeasurement,xBegROI,xEndROI,yBegROI,yEndROI,nLambda,LambdaList,nSteps,Plots);

%% Amplitude Calibration
% Note: This code can take a significant amount of time to run
% Default image size is 1500x1500 <- do we want to change this?

W = MM_AmplitudeCalibration(airMeasurement,nLambda,LambdaList,nSteps,caldata,Plots);

%% Analysis
MM_Analysis(sampleFilePath,xBegROI,xEndROI,yBegROI,yEndROI,nLambda,LambdaList,nSteps,Amp,W,Plots)
