% Last Updated: 2020/12/14

%Set initial parameters

%Measurement steps
nSteps = 64; %StepperCal = 64, RGB950 = 40

%Wavlengths of interest
LambdaList = 630;%, 630, 650, 670];

%Store number of wavelengths
nLambda = length(LambdaList);

%image size for dynamic programming (not done yet)
ROI = 1500;

%%
%Calibration file path
calibrationFilePath = 'D:\Measurements\Air_Calibrations\New_Diffuser\AirCalibration-06-Dec-2020.h5';
testFilePath = 'D:\Measurements\Air_Calibrations\New_Diffuser\Dichroic45-07-Dec-2020.h5';

%% Load Images
clear meanfulldata
for ii = 1:nLambda %for loading all lambda(s)
%Get h5 reference data
    Lambda = LambdaList(ii);
    refVecs(ii,:) = NormalizedReferenceData(calibrationFilePath,Lambda,nSteps);
%Grab measured data from h5 file
    meanfulldata(ii,:,:,:) = ReadCalImages(calibrationFilePath,Lambda,refVecs(ii,:),nSteps);
end

%% Set Fit Parameters

%Find the maximum value of air measurements
mx = max(meanfulldata,[],'all');

%Get amount of pixels
PixelCount = length(meanfulldata(1,1,:));

%AirCalRoutine_MKedit2(Amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,nSteps)
func = @(fitvals, inpvals) AirCalRoutine_MKeditAmpCal(fitvals,inpvals(1),inpvals(2),inpvals(3),inpvals(4),inpvals(5),inpvals(6));

%Set bounds
lb = 0;             %Lower bound;
ub = 3*mx;        %Upper bound;

%Starting fit guess
fitguess =  mx;      %maximum value. Try using found amplitudes from spectral data: 1.062119705563137e+05
%% Set fit input values from Spectral Calibration

%values constant across Lambda (values in radians)
DelRad_g = -0.1825;     %PSG theta
DelRad_a = -0.0203;     %PSA theta
LP_Rad = -0.0610;       %Linear Polarizer Theta

%Full set of PSG and PSA retardances
dg_rad = [2.32821491720789,2.27923088133178,2.22232304280817,2.17410616636198,2.13309067300460,2.08198339906172,2.03090094594950,1.98507651881251,1.94578326429503,1.90544610463194,1.87293511987322,1.83200282613381,1.79245046412239,1.76172910227583,1.73727986513025,1.70695826521702,1.68115704522444,1.65592990218593,1.63405542419882,1.59480380716638,1.59024913139111,1.55487026365334,1.53312571366529,1.52411297116052,1.49549773451957,1.48652107775223,1.45656538537525,1.44059443821981,1.41146637680512,1.39566454487169,1.36910002336742];
da_rad = [2.21508374156057,2.16956964041893,2.12900440198022,2.09058916286353,2.04957236778996,2.00464758675349,1.96036373992076,1.92239496140048,1.88570322785251,1.85007659839389,1.81971725700290,1.78187920325098,1.74589814435504,1.71631647898904,1.69237436610493,1.66448524608816,1.63847463328053,1.61426210116849,1.59324181399593,1.55519471336503,1.54936618329363,1.51539823118199,1.49384390495641,1.48502994775830,1.45652115908295,1.44694983689127,1.41930418343572,1.40380613905006,1.37684807652435,1.36337808695621,1.33839227147499];

%Extract specific retardances for fitting amplitude
dg_rad_full = [dg_rad(1),dg_rad(10),dg_rad(13),dg_rad(19),dg_rad(21), dg_rad(23)];
da_rad_full = [da_rad(1),da_rad(10),da_rad(13),da_rad(19),da_rad(21), da_rad(23)];

%% Pixel by Pixel Fit
%Least squares curve fitting and storing into caldata
%Will only work for square images atm
tic
clear ampcaldata
% ampcaldata = zeros(nLambda,sqrt(PixelCount),sqrt(PixelCount)); %Reserve memory
for n = 1:nLambda %across all wavelengths
    inputvec = [dg_rad_full(6), DelRad_g, da_rad_full(6), DelRad_a, LP_Rad, nSteps]; %Input values vector 
    for ii = 1:length(1:sqrt(PixelCount))                                  %across all pixel rows
        for jj = 1:length(1:sqrt(PixelCount))                              %across all pixel columns
            [fits,resnorm,res] = lsqcurvefit(func, fitguess, inputvec, squeeze(meanfulldata(n,:,ii,jj)),lb,ub); %Least-squares curve fit routine
            ampcaldata(n,ii,jj) = fits;                                    %Grab fitted amplitude terms
            residualError(n,ii,jj) = resnorm;
            RMSE_var(n,ii,jj,:) = std(res(:));                             
        end
    end
end
toc

%Play Hallelujah when finished
load handel
sound(y,Fs)

%% Plotting fits

%Steps for plotting
ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 

%Plot Data
figure(4)
hold on;
plot(ThetaMotorGen,mean(squeeze(meanfulldata(4,:,:,:)),[2,3]),'*r'); %Plot measured data points
Irrad = AirCalRoutine_MKeditLSQ(ampcaldata(4,:,:),dg_rad_full,DelRad_g,da_rad_full,DelRad_a,LP_Rad, 256); %Plot fit

%Plot mean fit
range = linspace(0,2*pi,256);
plot(range,mean(squeeze(Irrad(1,:,:,:)),[2,3]));
title(['Lambda =  ' num2str(LambdaList(4))]);
xlabel('PSG Rotation (Radians)')
ylabel('Camera Counts')

%% Create W-matrix

[~,W_amp] = AirCalRoutine_MKedit2(ampcaldata,dg_rad(19),DelRad_g,da_rad(19),DelRad_a,LP_Rad,nSteps);
 
%%
for n = 1:nLambda 
    tic
        Wtest(n,:,:) = AirCalRoutine_MKeditAmpCalW(dg_rad_full(n),DelRad_g,da_rad_full(n),DelRad_a,LP_Rad,nSteps);
    toc
end
W = squeeze(Wtest);
%% Load Image data
clear imgdata;
for ii = 1:length(LambdaList) %for loading all lambda(s)
%Get h5 reference data
    Lambda = LambdaList(ii);
    [~,~,refVecsImgs(ii,:)] = load_refdata(testFilePath,Lambda,nSteps); %grab reference data
    imgdata(ii,:,:,:) = ReadDataImages(testFilePath, Lambda,refVecsImgs(ii,:),nSteps); %image data with reference
end


%% Create W^-1 and Mueller terms

tic
clear M
for n = 1:nLambda
    Winv = pinv(squeeze(W_amp(1,:,:))); %Create W^-1
    M(n,:,:) = Winv*(reshape(squeeze(imgdata(n,:,:,:)),nSteps,ROI*ROI)./repmat(ampcaldata(n,:),[nSteps 1])); %Create Mueller terms
end
toc

%% Rearrange Mueller terms into Mueller Matrix
clear mmVecs
for n = 1:nLambda
    mmVecs(n,:,:,:) = reshape(squeeze(M(n,:,:)),16,ROI,ROI); 
end

%% Save Mueller Matrix Images
for ii = 1:nLambda
    filename = ["670Amp45" num2str(LambdaList(nLambda))];
    DisplayMM(mmVecs,nLambda,LambdaList,filename)
end