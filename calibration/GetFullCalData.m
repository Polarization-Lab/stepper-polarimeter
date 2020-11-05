
nSteps = 64; %StepperCal = 64, RGB950 = 40
LambdaList = linspace(450,750,4);
nLambda = length(LambdaList);

%image size. Needs to be more robust.
testimsz_1 = 1500;
testimsz_2 = 1500;

%%
%Calibration file path
calibrationFilePath = 'D:\Measurements\VVR_Analysis\AirCalibrations-15-Sep-2020.h5';
testFilePath = 'C:\Stepper Polarimeter\VVR_Analysis_Temp\VVR-16-Sep-2020.h5';

for ii = 1:length(LambdaList) %for loading all lambda(s)
%Get h5 reference data
    Lambda = LambdaList(ii);
%     refVecs(ii,:) = NormalizedReferenceData(calibrationFilePath,Lambda,nSteps);

%Grab measured data from h5 file
    meanfulldata(ii,:,:,:) = ReadDataImagesNoRef(calibrationFilePath,Lambda,nSteps);
%     meanfulldata(ii,:,:,:) = ReadFullCalImages(calibrationFilePath,Lambda,refVecs,nSteps); %needs to use dark field correction code at some point

end

%%
clc;
mx = max(meanfulldata,[],'all');
%Get amount of pixels
PixelCount = length(meanfulldata(1,1,:));

%meandataVecs = reshape(meandata,1,64,10,1

%Get a function for fitting the values
%Get spectrum data for "lambda colored" plots
% sRGB = spectrumRGB(Lambda, '1964_FULL');
% cRGB(:) = [sRGB(:,1) sRGB(:,2) sRGB(:,3)]; %Colors for graphing

%AirCalRoutine_MKedit2(Amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,nSteps)
func = @(fitvals, inpvals) AirCalRoutine_MKeditAmpCal(fitvals,inpvals(1),inpvals(2),inpvals(3),inpvals(4),inpvals(5),inpvals(6));
%Set bounds and known values, likely not optimal ask Meredith
% lb = [Amp ];
% ub = [3*mx*ones(1,PixelCount) pi*ones(1,PixelCount) pi*ones(1,PixelCount)];
% lb = mx*ones(1,PixelCount*nLambda);
% ub = 3*mx*ones(1,PixelCount*nLambda);
lb = 0;
ub = 3*mx;

%Input vector
% inpvec = [PixelCount, nSteps, ThetaMotorGen, ThetaMotorAna];

%Starting fit guess
% fitguess = [2*mx*ones([1,PixelCount]), 2*pi/3*ones([1,PixelCount]), 2*pi/3*ones([1,PixelCount])];, 2*pi/3, 0, 2*pi/3, 0, 0
fitguess = 2*mx;
%%

%constant across Lambda
% DelRad_g = 0.093174;
% DelRad_a = -0.041596;
% LP_Rad = 0.16073;
% 
% 
% dg_rad = [2.3471 2.2977 2.2434 2.1702 2.1280 2.0916 2.0450 1.9948 1.8772 1.8676 1.8914 1.8506 1.8142 1.8068 1.6952 1.6851 1.6444 1.5705 1.5951 1.6081 1.5555 1.5155 1.4630 1.5619 1.4011 1.4545 1.3613 1.4504 1.4149 1.3450];
% da_rad = [2.1675 2.1252 2.0909 2.0463 2.0107 1.9807 1.9509 1.9184 1.8163 1.8165 1.8448 1.8112 1.7840 1.7801 1.6734 1.6650 1.6264 1.5537 1.5787 1.5961 1.5414 1.5022 1.4497 1.5532 1.3878 1.4387 1.3482 1.4381 1.4008 1.3320];
% 
dg_rad_full = [dg_rad(1),dg_rad(11),dg_rad(21),dg_rad(31)];
da_rad_full = [da_rad(1),da_rad(11),da_rad(21),da_rad(31)];
%%
%Least squares curve fitting and storing into caldata
%Will only work for square images atm
clear ampcaldata
ampcaldata = zeros(nLambda,sqrt(PixelCount),sqrt(PixelCount));
for n = 1:nLambda
    inputvec = [dg_rad_full, DelRad_g, da_rad_full, DelRad_a, LP_Rad, nSteps];
    [fits,resnorm,res] = lsqcurvefit(func, fitguess, inputvec, meanfulldata,lb,ub);
    ampcaldata(:) = fits(:);
end
f = msgbox('Operation Completed');
% load handel
% sound(y,Fs)
%%
%Least squares curve fitting and storing into caldata
%Will only work for square images atm
clear ampcaldata
ampcaldata = zeros(nLambda,sqrt(PixelCount),sqrt(PixelCount));
for n = 1:nLambda
    inputvec = [dg_rad(n), DelRad_g, da_rad(n), DelRad_a, LP_Rad, nSteps];
    for ii = 1:length(1:sqrt(PixelCount))
        for jj = 1:length(1:sqrt(PixelCount))
            [fits,resnorm,res] = lsqcurvefit(func, fitguess, inputvec, meanfulldata(n,:,ii,jj),lb,ub);
            ampcaldata(n,ii,jj) = fits;
        end
    end
end


% load handel
% sound(y,Fs)
% caldata(:)=fits(:);
%%

imagesc(squeeze(ampcaldata(5,:,:)));colorbar;
title('Amplitude fit Lambda = 490 nm')

%%
%RMSE of data
for ii = 1:nLambda
RMSE_var(ii) = std(res(ii,:)); %will need to fix size (1 x nLambda)
end
%Grab caldata for tables, plotting, and psgMM/psaMM
%Create table in radians
% Amp = caldata(1:PixelCount*nLambda);
% dg_rad = caldata(PixelCount*nLambda+1:PixelCount*nLambda+nLambda); %change range
% DelRad_g = caldata(PixelCount*nLambda+1+nLambda);
% da_rad = caldata(PixelCount*nLambda+2+nLambda:PixelCount*nLambda+1+2*nLambda); %change range
% DelRad_a = caldata(PixelCount*nLambda+2+2*nLambda);
% LP_Rad = caldata(PixelCount*nLambda+3+2*nLambda);

%plot lambda vs del_a and del_g
%plot Amplitude vs meas# for nLambda
%print Theta_LP,PSA,PSG


ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 
%%
figure(1)
plot(LambdaList,dg_rad)
title('\delta_{g} radians');
figure(2)
plot(LambdaList,da_rad)
title('\delta_{a} radians');
figure(3)
errorbar(LambdaList,RMSE_var);
title('RMSE vs Lambda');
%%
%Plot Data
hold on;
plot(ThetaMotorGen,mean(squeeze(meandata(5,:,:,:)),[2,3]),'*r'); 
Irrad = AirCalRoutine_MKeditLSQ(Amp,dg_rad,DelRad_g,da_rad,DelRad_a,LP_Rad, 629);

%Plot fit
plot(0:0.01:2*pi,mean(squeeze(Irrad(5,:,:,:)),[2,3]));%,'color',cRGB);
title(['Lambda =  ' num2str(LambdaList(5))]);
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
Avg_Amplitude = mean(Amp,'all');
RMSE = RMSE_var;


T = table(LambdaList,Avg_Amplitude,dg_rad, DelRad_g, da_rad, DelRad_a, LP_Rad, RMSE)

% Avg_delta_g = mean(dg_deg,'all');
% Avg_Theta_g = mean(DelDeg_g,'all');
% Avg_delta_a = mean(da_deg,'all');
% Avg_Theta_a = mean(DelDeg_a,'all');
% Avg_LP_Theta = mean(LP_Deg,'all');

T = table(LambdaList,Avg_Amplitude,dg_deg, DelDeg_g, da_deg, DelDeg_a, LP_Deg, RMSE)
% 
% figure(2)
% imagesc(reshape(Amp,5,5));colormap(gca,'parula');hcb1=colorbar;
% sgtitle(['Amplitude ' num2str(Lambda) 'nm']