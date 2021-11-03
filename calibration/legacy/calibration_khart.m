%This is my attempt at a general calibration code using existing functions 
% I've written for the stepper. It is as general as possible 

%% Set Up
%file path, assume folder is in path 
fn = 'air-18-Aug-2020.h5';

%measurement info
wavelength = 550;
num_meas = 64;

%choose the ROI to calibrate
xstart = 800;
xcount = 10;
ystart = 800;
ycount = 10;
PixelCount = xcount * ycount;

%% Import Data
%import image data
images= ims;
%images = load_imagedata(fn, wavelength ,num_meas, xstart,xcount,ystart,ycount);
%ims =images;

mean_ims = mean(reshape(ims,[100,64]));

%import reference data
%[ref_mean, ref_std] = load_refdata(fn,wavelength,num_meas);
 

%% Reference Correction
corr = ref_mean/max(ref_mean);
for i = 1:num_meas
    images(i) = images(i) ./ corr(i);
end


%% Perform Fit
func = @(fitvals, inpvals) AirCalFit(fitvals(1:PixelCount),fitvals(PixelCount+1),fitvals(PixelCount+2),fitvals(PixelCount+3),fitvals(PixelCount+4),fitvals(PixelCount+5),inpvals(1));

%set limits
mx = max(images,[],'all') * ones(1,PixelCount);

%re
images = reshape(images,[64,10,10]);

%Set upper (ub) and lower (lb) boundaries for fitting values
lb = [mx  -pi -pi/2 -pi -pi/2 -pi/2];
ub = [2*mx pi  pi/2  pi  pi/2  pi/2];

options = optimset('MaxFunEvals',1000000); %Indicates how many iterations
opts = optimoptions('lsqcurvefit', 'Diagnostics','on', 'Display','iter-detailed');

jj = 1; %measurement number
%Setup Generator angles
ThetaMotorGen = (0:num_meas-1)*2*pi/num_meas; 
%Setting up analyzer and image size
ThetaMotorAna = 4.9*ThetaMotorGen;

%Least squares curve fitting 
init = [3*mx, 2*pi/3, 0, 2*pi/3, 0, 0]; %starting position
[fits,resnorm(jj),res] = lsqcurvefit(func,init, num_meas, images,lb,ub,opts);
caldata(jj,:)=fits(:);

[Irrad,~] = AirCalFit(caldata(1),caldata(2),caldata(3),caldata(4),caldata(5),caldata(6), 629);

%Plot each wavelength's fit
plot(ThetaMotorGen,mean_ims,'*r'); 
hold on;
plot(0:0.01:2*pi,Irrad);
xlabel('PSG Rotation (Radians)')
ylabel('Camera Counts')
