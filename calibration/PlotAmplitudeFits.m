function PlotAmplitudeFits(ThetaMotorGen,airMeasurement,caldata,Lambda,nLambda,index)
%Inputs:
% ThetaMotorGen - data points for each step made by stepper
% airMeasurement - Air measurements loaded
% caldata - Calibration data
% Lambda - desired wavelength for plotting
% nLambda - number of lambda in calibration (used for indexing caldata)
% index - variable location of desired wavelength from LambdaList



% Separate fitted data into each desired fit value (angles in radians)
Amp = caldata(1:PixelCount*nLambda); %Amplitude
dg_rad = caldata(PixelCount*nLambda+1:PixelCount*nLambda+nLambda);             %PSG retardance
DelRad_g = caldata(PixelCount*nLambda+1+nLambda);                              %PSG theta
da_rad = caldata(PixelCount*nLambda+2+nLambda:PixelCount*nLambda+1+2*nLambda); %PSA retardance
DelRad_a = caldata(PixelCount*nLambda+2+2*nLambda);                            %PSA theta
LP_Rad = caldata(PixelCount*nLambda+3+2*nLambda);                              %Linear Polarizer theta

%Plot Data
hold on;
plot(ThetaMotorGen,mean(squeeze(airMeasurement(index,:,:,:)),[2,3]),'*r'); %Measured
Irrad = AirCalRoutine_MKeditLSQ(Amp,dg_rad,DelRad_g,da_rad,DelRad_a,LP_Rad, 629); %Fits, 629 points used for smooth line

%Plot fit
plot(0:0.01:2*pi,mean(fliplr(squeeze(Irrad(index,:,:,:))),[2,3]));%,'color',cRGB);
title(['Lambda =  ' num2str(Lambda)]);
xlabel('PSG Rotation (Radians)')
ylabel('Camera Counts')

%Create table in degrees
dg_deg = dg_rad.*180/pi;
DelDeg_g = DelRad_g.*180/pi;
da_deg = da_rad.*180/pi;
DelDeg_a = DelRad_a.*180/pi;
LP_Deg = LP_Rad.*180/pi;
