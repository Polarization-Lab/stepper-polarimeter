function []=GenMuellerVar2(fn_air,fn_sample)

%Function to read an h5 air measurement Stepper system parameters to create a
%W matrix. Then reconstruct the MM of a sample measurement from that W
%matrix. MM reconstruction is saved to same file as the sample measurement

close all; tic
addpath  /Users/meredith/Documents/MATLAB/stepper-polarimeter-master/calibration
addpath('/Users/meredith/Library/Mobile Documents/com~apple~CloudDocs/Documents/MATLAB/stepper-polarimeter-master/calibration')

%fn_air='/Volumes/StepperData/Measurements/Air_Calibrations/air_20210710_400-800-50nm/airMeas-400_800-50nm_TeleLens-10-Jul-2021.h5';
%fn_sample='/Volumes/StepperData/Measurements/Waveguides/Option2/Angle Measurements/Run6_07122021/sampleMeas-460_530_625-12-Jul-2021.h5';
PSG=h5readatt(fn_sample,'/images','PSG_positions');
nSteps = length(PSG);

waves=h5info(fn_sample,'/images/');
waves={waves.Groups.Name};
LambdaList=zeros(1,length(waves));
nLambda=length(LambdaList);
for l=1:nLambda
    extractAfter(waves(l),'wave');
    LambdaList(l)=str2double(char(extractAfter(waves(l),'wave')));
end

rect = GetROI(fn_sample,LambdaList(1)); 
rect(4) = rect(3); %make ROI square
nPixels=(rect(4)+1)^2; 

%Write MM to file by creating a new index
info = h5info(fn_sample);
index=length({info.Groups.Name})-1;
%h5create(fn_sample,strcat('/Mueller',num2str(index)),[nLambda,16,sqrt(nPixels),sqrt(nPixels)]);

%Movie Output
MRGB=VideoWriter([extractBefore(fn_sample,'.h5'),sprintf('_MM%d.mp4',index)],'MPEG-4');
MRGB.FrameRate = 1;
open(MRGB);set(0,'DefaultFigureVisible', 'off');

for n=1:nLambda
    W = squeeze(GenWMatrixFromStepperPars(fn_air,LambdaList(n))); 
    temp = reshape(ReadCalImages(fn_sample,LambdaList(n),ones(1,nSteps),nSteps,rect),nSteps,nPixels);
    MM=pinv(W)*temp;
    ViewMM(reshape(MM,[4 4 sqrt(nPixels) sqrt(nPixels)]));
    title(sprintf('$\\lambda$=%d[nm]',LambdaList(n)),'FontSize',30);
    f=getframe(gcf);
    writeVideo(MRGB,f);close all;
    %Append MM to h5 file
    h5create(fn_sample,strcat('/Mueller',num2str(index),'/wave',num2str(LambdaList(n))),[16,nPixels]);
    h5write(fn_sample,strcat('/Mueller',num2str(index),'/wave',num2str(LambdaList(n))),MM);

end
close(MRGB);close all;set(0,'DefaultFigureVisible', 'on');
h5writeatt(fn_sample,strcat('/Mueller',num2str(index),'/'),'Air_Measurement_file',fn_air);
h5writeatt(fn_sample,strcat('/Mueller',num2str(index),'/'),'ROI',rect);

