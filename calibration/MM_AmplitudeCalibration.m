function W = AmplitudeCalibration(airMeasurement,nLambda,LambdaList,nSteps,caldata,Plots)
%% Set Fit Parameters

disp('Setting up fit parameters')

%Find the maximum value of air measurements
mx = max(airMeasurement,[],'all');

%Get amount of pixels
PixelCount = length(airMeasurement(1,1,:));

%AirCalRoutine_MKedit2(Amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,nSteps)
func = @(fitvals, inpvals) AirCalRoutine_MKeditAmpCal(fitvals,inpvals(1),inpvals(2),inpvals(3),inpvals(4),inpvals(5),inpvals(6));

%Set bounds
lb = 0;             %Lower bound;
ub = 3*mx;        %Upper bound;

%Starting fit guess
fitguess =  mx;      %maximum value. Try using found amplitudes from spectral data: 1.062119705563137e+05

%% Grab caldata values
% Separate fitted data into each desired fit value (angles in radians)
Amp = caldata(1:PixelCount*nLambda); %Amplitude
dg_rad = caldata(PixelCount*nLambda+1:PixelCount*nLambda+nLambda);             %PSG retardance
DelRad_g = caldata(PixelCount*nLambda+1+nLambda);                              %PSG theta
da_rad = caldata(PixelCount*nLambda+2+nLambda:PixelCount*nLambda+1+2*nLambda); %PSA retardance
DelRad_a = caldata(PixelCount*nLambda+2+2*nLambda);                            %PSA theta
LP_Rad = caldata(PixelCount*nLambda+3+2*nLambda);                              %Linear Polarizer theta
%% Set fit input values from Spectral Calibration

%Extract specific retardances for fitting amplitude
dg_rad_full = dg_rad(nLambda);
da_rad_full = da_rad(nLambda);

%% Pixel by Pixel Fit
%Least squares curve fitting and storing into caldata
%Will only work for square images atm

disp('Fitting amplitude curves')

tic
clear ampcaldata
% ampcaldata = zeros(nLambda,sqrt(PixelCount),sqrt(PixelCount)); %Reserve memory
 inputvec = [dg_rad_full, DelRad_g, da_rad_full, DelRad_a, LP_Rad, nSteps]; %Input values vector 
for ii = 1:length(1:sqrt(PixelCount))                                  %across all pixel rows
    for jj = 1:length(1:sqrt(PixelCount))                              %across all pixel columns
        [fits,resnorm,res] = lsqcurvefit(func, fitguess, inputvec, squeeze(airMeasurement(n,:,ii,jj)),lb,ub); %Least-squares curve fit routine
        ampcaldata(n,ii,jj) = fits;                                    %Grab fitted amplitude terms
        residualError(n,ii,jj) = resnorm;
        RMSE_var(n,ii,jj,:) = std(res(:));                             
    end
end

toc

%Play Hallelujah when finished
disp('Curve fitting finished')
load handel
sound(y,Fs)

%% Plotting fits

%Steps for plotting
ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 

%Plot Data
figure(4)
hold on;
plot(ThetaMotorGen,mean(squeeze(airMeasurement(nLambda,:,:,:)),[2,3]),'*r'); %Plot measured data points
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