function [caldata,W] = MM_Calibration(xBegROI,xEndROI,yBegROI,yEndROI,nLambda,LambdaList,nSteps)
%% Load Calibration data
tic
for ii = 1:length(LambdaList) %for loading all lambda(s)
%Get h5 reference data
    Lambda = LambdaList(ii);
    [~,~,refVecs(ii,:)] = load_refdata(calibrationFilePath,Lambda,nSteps);

%Grab measured data from h5 file
    airMeasurement(ii,:,:,:) = ReadCalImages(calibrationFilePath,Lambda,refVecs(ii,:),nSteps);
end
toc
%% Set size of meandata (Maximum size is 10 x 10 due to memory)

airMeasurement = airMeasurement(1,:,yBegROI:yEndROI,xBegROI:xEndROI); %Grab middle 10 x 10 section


%% Fit Parameters

clc;

%Find the maximum value of air measurements
mx = max(airMeasurement,[],'all');
%Get amount of pixels
PixelCount = length(airMeasurement(1,1,:));

%Setup function for fitting (see below for input values):
%AirCalRoutine_MKeditLSQ(Amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,nSteps)
func = @(fitvals, inpvals) AirCalRoutine_MKeditLSQ(fitvals(1:PixelCount*nLambda),fitvals(PixelCount*nLambda+1:PixelCount*nLambda+nLambda),fitvals(PixelCount*nLambda+1+nLambda),fitvals(PixelCount*nLambda+2+nLambda:PixelCount*nLambda+1+2*nLambda),fitvals(PixelCount*nLambda+2+2*nLambda),fitvals(PixelCount*nLambda+3+2*nLambda),inpvals);

%Set bounds
lb = [mx*zeros(1,PixelCount*nLambda) -pi*ones(1,nLambda) -pi/2 -pi*ones(1,nLambda) -pi/2 -pi/2]; %Lower bound
ub = [4*mx*ones(1,PixelCount*nLambda) pi*ones(1,nLambda) pi/2 pi*ones(1,nLambda) pi/2 pi/2];     %Upper bound

%Starting fit guess (Amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP)
fitguess = [2*mx*ones([1,PixelCount*nLambda]), 2*pi/3*ones(1,nLambda), 0, 2*pi/3*ones(1,nLambda), 0, 0];
%% Curve fitting
%Least squares curve fitting and storing into caldata
tic
[fits,resnorm,res] = lsqcurvefit(func, fitguess, nSteps, airMeasurement,lb,ub);
caldata(:)=fits(:); %Grab fitted data
toc

%Indicators for code finishing
f = msgbox('Operation Completed');
load handel
sound(y,Fs)

%%
%RMSE of data
for ii = 1:nLambda
RMSE_var(ii) = std(res(ii,:)); %will need to fix size (1 x nLambda)
end
%Grab caldata for tables, plotting, and psgMM/psaMM
%Create table in radians
%% Separate fitted data into each desired fit value (angles in radians)
Amp = caldata(1:PixelCount*nLambda); %Amplitude
dg_rad = caldata(PixelCount*nLambda+1:PixelCount*nLambda+nLambda);             %PSG retardance
DelRad_g = caldata(PixelCount*nLambda+1+nLambda);                              %PSG theta
da_rad = caldata(PixelCount*nLambda+2+nLambda:PixelCount*nLambda+1+2*nLambda); %PSA retardance
DelRad_a = caldata(PixelCount*nLambda+2+2*nLambda);                            %PSA theta
LP_Rad = caldata(PixelCount*nLambda+3+2*nLambda);                              %Linear Polarizer theta


%% Plotting retardance fit data
%Steps for plotting
ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 

%Plot retardance curves (PSG and PSA delta)
figure(1)
plot(LambdaList,dg_rad,'*b')
xlabel('\lambda (nm)')
title('\delta_{g} radians');
figure(2)
plot(LambdaList,da_rad,'*b')
title('\delta_{a} radians');
xlabel('\lambda (nm)')
%% Plot RMSE
figure(3)
errorbar(LambdaList,RMSE_var);
title('RMSE vs Lambda');
xlabel('\lambda (nm)')
%% Plot amplitude Fits
close all;
%Plot Data
hold on;
plot(ThetaMotorGen,mean(squeeze(airMeasurement(1,:,:,:)),[2,3]),'*r'); %Measured
Irrad = AirCalRoutine_MKeditLSQ(Amp,dg_rad,DelRad_g,da_rad,DelRad_a,LP_Rad, 629); %Fits

%Plot fit
plot(0:0.01:2*pi,mean(fliplr(squeeze(Irrad(1,:,:,:))),[2,3]));%,'color',cRGB);
title(['Lambda =  ' num2str(LambdaList(1))]);
xlabel('PSG Rotation (Radians)')
ylabel('Camera Counts')

%Create table in degrees
dg_deg = dg_rad.*180/pi;
DelDeg_g = DelRad_g.*180/pi;
da_deg = da_rad.*180/pi;
DelDeg_a = DelRad_a.*180/pi;
LP_Deg = LP_Rad.*180/pi;

toc
%%
%Creating W Matrix
tic
[~,W] = AirCalRoutine_MKedit2(Amp,dg_rad,DelRad_g,da_rad,DelRad_a,LP_Rad,nSteps);
toc
end