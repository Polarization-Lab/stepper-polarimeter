function []=GenMuellerVar2(fn_air,fn_sample)

%fn_air: h5 file air measurement from Stepper 
%fn_sample: h5 file sample measurement from Stepper
%fn_air in input to function GenWMatrixFromStepperPars that creates a W matrix (see EQ 7.3 In Polarized Light and Optical Systems) 
%W matrix reconstructs the Mueller matrix of fn_sample measurement
%Mueller matrix and singular values and vectors of coherency matrix (see EQ 14 https://doi.org/10.1364/OE.425295) are appended to fn_sample h5 file
%A movie of Mueller matrix over wavelength is saved in same folder as fn_sample

addpath  /Users/meredith/Documents/MATLAB/stepper-polarimeter-master/calibration
addpath('/Users/meredith/Library/Mobile Documents/com~apple~CloudDocs/Documents/MATLAB/stepper-polarimeter-master/calibration')
addpath('/Users/meredith/Documents/MATLAB/stepper-polarimeter-master')
%fn_air='/Volumes/StepperData/Measurements/Air_Calibrations/air_20210908_400-800-50/airmeas400-800-50-08-September-2021-09-Sep-2021.h5'
PSG=h5readatt(fn_sample,'/images','PSG_positions');
nSteps = length(PSG);

waves=h5info(fn_sample,'/images/');
waves={waves.Groups.Name};
LambdaList=zeros(1,length(waves));
nLambda=length(LambdaList);
for l=1:nLambda
    w=char(waves(l));
    LambdaList(l)=str2double(w((end-2):end));
end

rect = GetROI(fn_sample,LambdaList(1)); 
rect(4) = rect(3); %make ROI square
nPixels=(rect(4)+1)^2; 

%Write MM to file by creating a new index
info = h5info(fn_sample);
index=length({info.Groups.Name})-1;

%Movie Output
w=char(fn_sample);%clip out ".h5" part of filename
MRGB=VideoWriter([w(1:end-3),sprintf('_MM%d.mp4',index)],'MPEG-4');
MRGB.FrameRate = 1;
open(MRGB);set(0,'DefaultFigureVisible', 'off');

%Analysis for cohereny matrix
U = [1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 sqrt(-1) -sqrt(-1) 0];
PI = zeros(4,4,4,4);
for n=1:4
	for m=1:4
PI(n,m,:,:) = 0.5*U*(kron(reshape(U(n,:),2,2),conj(reshape(U(m,:),2,2))))*U';
	end
end
PI=reshape(PI,16,16);
C=zeros(nLambda,nPixels,16);

for n=1:nLambda
    W = squeeze(GenWMatrixFromStepperPars(fn_air,LambdaList(n))); 
    temp = reshape(ReadCalImages(fn_sample,LambdaList(n),ones(1,nSteps),nSteps,rect),nSteps,nPixels);
    MM=pinv(W)*temp;%16xnPixels matrix
    C(n,:,:)=0.25*MM'*PI;%nPixelsx16 matrix
    ViewMM(reshape(MM,[4 4 sqrt(nPixels) sqrt(nPixels)]));
    title(sprintf('$\\lambda$=%d[nm]',LambdaList(n)),'FontSize',30,'interpreter','latex');
    f=getframe(gcf);
    writeVideo(MRGB,f);close all;
    %Append MM to h5 file
    h5create(fn_sample,strcat('/Mueller',num2str(index),'/wave',num2str(LambdaList(n))),[16,nPixels]);
    h5write(fn_sample,strcat('/Mueller',num2str(index),'/wave',num2str(LambdaList(n))),MM);

end
close(MRGB);close all;set(0,'DefaultFigureVisible', 'on');
h5writeatt(fn_sample,strcat('/Mueller',num2str(index),'/'),'Air_Measurement_file',fn_air);
h5writeatt(fn_sample,strcat('/Mueller',num2str(index),'/'),'ROI',rect);

%Compute and Save Coherency Matrix SVD
UC=zeros(nPixels,4,4);
S=zeros(nPixels,4);

for n=1:nLambda
    for nn=1:nPixels
        [UC(nn,:,:),D] = svd(reshape(C(n,nn,:,:),4,4));%C is Hermitian
        S(nn,:)=diag(D);
    end
    h5create(fn_sample,strcat('/Real_CoherencySingularVectors',num2str(index),'/wave',num2str(LambdaList(n))),[nPixels,4,4]);
    h5write(fn_sample,strcat('/Real_CoherencySingularVectors',num2str(index),'/wave',num2str(LambdaList(n))),real(UC));
    h5create(fn_sample,strcat('/Imag_CoherencySingularVectors',num2str(index),'/wave',num2str(LambdaList(n))),[nPixels,4,4]);
    h5write(fn_sample,strcat('/Imag_CoherencySingularVectors',num2str(index),'/wave',num2str(LambdaList(n))),imag(UC));
    h5create(fn_sample,strcat('/CoherencySingularValues',num2str(index),'/wave',num2str(LambdaList(n))),[nPixels,4]);
    h5write(fn_sample,strcat('/CoherencySingularValues',num2str(index),'/wave',num2str(LambdaList(n))),S);
    clear UC;
    clear S;
end

