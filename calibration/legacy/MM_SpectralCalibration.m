function [caldata,W] = MM_SpectralCalibration(airMeasurement, xBegROI,xEndROI,yBegROI,yEndROI,nLambda,LambdaList,nSteps,Plots)
%% Set size of meandata (Maximum size is 10 x 10 due to memory)

airMeasurement = airMeasurement(1,:,yBegROI:yEndROI,xBegROI:xEndROI); %Grab desired ROI

%% Fit Parameters 

disp('Setting up fit parameters')

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

disp('Fitting curves')

tic
[fits,resnorm,res] = lsqcurvefit(func, fitguess, nSteps, airMeasurement,lb,ub);
caldata(:)=fits(:); %Grab fitted data
toc

%Indicators for code finishing
disp('Curve fitting finished')
load handel
sound(y,Fs)

%%
%RMSE of data
for ii = 1:nLambda
RMSE_var(ii) = std(res(ii,:)); %will need to fix size (1 x nLambda)
end

%% Separate PSG and PSA retardance for plotting
dg_rad = caldata(PixelCount*nLambda+1:PixelCount*nLambda+nLambda);             %PSG retardance
da_rad = caldata(PixelCount*nLambda+2+nLambda:PixelCount*nLambda+1+2*nLambda); %PSA retardance

%% Plotting retardance fit data
%Steps for plotting
ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 
if (Plots == 1)
    
    disp('Plotting retardance curves and RMSE')
    
    %Plot retardance curves (PSG and PSA delta)
    figure(1)
    plot(LambdaList,dg_rad,'*b')
    xlabel('\lambda (nm)')
    title('\delta_{g} radians');
    figure(2)
    plot(LambdaList,da_rad,'*b')
    title('\delta_{a} radians');
    xlabel('\lambda (nm)')
    
    %Plot RMSE
    figure(3)
    errorbar(LambdaList,RMSE_var);
    title('RMSE vs Lambda');
    xlabel('\lambda (nm)')
end
%%
%Creating W Matrix

disp('Creating W Matrix')
tic
[~,W] = AirCalRoutine_MKedit2(Amp,dg_rad,DelRad_g,da_rad,DelRad_a,LP_Rad,nSteps);
toc

disp('Spectral Calibration finished')
end
