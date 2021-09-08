function [ROIdata,angle_hat]=GenFit2AirMeasVar2(fn,rep)
%For the input path to air measurement this function produces parameteric
%fit for rotation retarder polarimeter

%addpath  /Users/meredith/Documents/MATLAB/stepper-polarimeter-master/calibration
addpath('/Users/meredith/Library/Mobile Documents/com~apple~CloudDocs/Documents/MATLAB/stepper-polarimeter-master/calibration')

%fn='/Volumes/StepperData/Measurements/Air_Calibrations/air_20210710_400-800-50nm/airMeas-400_800-50nm_TeleLens-10-Jul-2021.h5';

PSG=h5readatt(fn,'/images','PSG_positions');
nSteps = length(PSG);%Number retarder rotations

waves=h5info(fn,'/images/');
waves={waves.Groups.Name};
LambdaList=zeros(1,length(waves));
nLambda=length(LambdaList);%Number of wavelengths measured
for l=1:nLambda
    extractAfter(waves(l),'wave');
    LambdaList(l)=str2double(char(extractAfter(waves(l),'wave')));
end

rect = GetROI(fn,LambdaList(1)); 
rect(3) = 3; %make 3 pixel wide ROI
rect(4) = rect(3); %make ROI square
nPixels=(rect(4)+1)^2; 
ROIdata=zeros(nLambda,nSteps,sqrt(nPixels),sqrt(nPixels));

for n=1:nLambda
	
    refVecs = NormalizedReferenceData(fn,LambdaList(n),nSteps);
    ROIdata(n,:,:,:) = ReadCalImages(fn,LambdaList(n),refVecs,nSteps,rect);

end


%Function handles
func_angles = @(fitvals, inpvals) GenAirMeasFromStepperParsV2(inpvals,fitvals,nSteps,LambdaList);
func_amp = @(fitvals, inpvals) GenAirMeasFromStepperParsV2(fitvals,inpvals,nSteps,LambdaList);
func_all = @(fitvals, inpvals) GenAirMeasFromStepperParsV2(fitvals(1:(nLambda*nPixels)),fitvals((nLambda*nPixels)+1:end),inpvals(1),inpvals(2:end));

%Angle initial values, upper and lower bounds
lb_angles = [-2*pi/3*(1/500)*2, 3*pi/3, -2*pi/3*(1/500)*2, 3*pi/3, -pi/2, -pi/2, -pi/2]; 
angle_hat = [-2*pi/3*(1/500), 4*pi/3, -2*pi/3*(1/500), 4*pi/3, 0, 0, 0];
ub_angles = [-2*pi/3*(1/500)*0.5, 10*pi/3, -2*pi/3*(1/500)*0.5, 10*pi/3, pi/2,  pi/2, pi/2];     

%Amplitude initial values, upper and lower bounds
mx = repmat(squeeze(max(max(max(ROIdata,[],4),[],3),[],2)),[1 nPixels]);
lb_amp = 0*mx(:)';
amp_hat = 2*mx(:)';
ub_amp = 4*mx(:)';

[all_hat] = lsqcurvefit(func_all, [amp_hat, angle_hat],[nSteps,LambdaList], ROIdata,[lb_amp,lb_angles],[ub_amp,ub_angles]);%fit all parameters
amp_hat=all_hat(1:(nLambda*nPixels));
angle_hat=all_hat((nLambda*nPixels)+1:end);

resnorm=zeros(rep,1);
res=zeros(rep,nLambda,nSteps*nPixels);
all_hat=zeros(rep,length(angle_hat)+length(amp_hat));

for x=1:rep
    
    [angle_hat] = lsqcurvefit(func_angles, angle_hat, amp_hat, ROIdata,lb_angles,ub_angles);                                %fit only angle parameters
    [amp_hat] = lsqcurvefit(func_amp, amp_hat, angle_hat, ROIdata,lb_amp,ub_amp);                               %fit only amplitude parameters
    [all_hat(x,:),resnorm(x),res(x,:,:)] = lsqcurvefit(func_all, [amp_hat, angle_hat],[nSteps,LambdaList], ROIdata,[lb_amp,lb_angles],[ub_amp,ub_angles]);%fit all parameters

    %Random guess
     angle_hat=lb_angles + (ub_angles-lb_angles).*rand(size(lb_angles));
     amp_hat=lb_amp + (ub_amp-lb_amp).*rand(size(lb_amp));

end
%Pick best solution of all reps 
[~,x]=min(resnorm);
all_hat=all_hat(x,:);
amp_hat=all_hat(1:(nLambda*nPixels));
angle_hat=all_hat((nLambda*nPixels)+1:end);
res_lambda=sqrt(mean(res(x,:,:).^2,2));

%Organize fit parameter output
start=1;
a_PSG = angle_hat(start); %PSG retardance magnitude slope and intercept
start=start+1;
b_PSG = angle_hat(start);               
start=start+1;
a_PSA = angle_hat(start);%PSA retardance magnitude slope and intercept
start=start+1;
b_PSA = angle_hat(start);               
start=start+1;
PSG_theta = angle_hat(start);  %PSG retarder orientaiton 
start=start+1;
PSA_theta = angle_hat(start);  %PSA retarder orientaiton 
start=start+1;
PSA_LP = angle_hat(start);     %PSA linear Polarizer theta 

%Write estimated system parameters to h5 file
%h5create(fn,'/cal',size(all_hat));
h5writeatt(fn,'/cal/','PSG_retardance_slope', a_PSG);
h5writeatt(fn,'/cal/','PSG_retardance_intercept', b_PSG);
h5writeatt(fn,'/cal/','PSG_theta', PSG_theta);
h5writeatt(fn,'/cal/','PSA_retardance_slope', a_PSA);
h5writeatt(fn,'/cal/','PSA_retardance_intercept', b_PSA);
h5writeatt(fn,'/cal/','PSA_theta', PSA_theta);
h5writeatt(fn,'/cal/','PSA_LP', PSA_LP);
h5writeatt(fn,'/cal/','ROI', rect);
h5writeatt(fn,'/cal/','Amp', amp_hat);
    
%Movie Output 
n_pixels=1:sqrt(nPixels);
m_pixels=1:sqrt(nPixels);

MRGB=VideoWriter([extractBefore(fn,'.h5'),'_RMSE_fit.mp4'],'MPEG-4');
MRGB.FrameRate = 1;
open(MRGB);

set(0,'DefaultFigureVisible', 'off');
for l=1:nLambda
    figure;plot((0:nSteps-1)*2*pi/nSteps,mean((ROIdata(l,:,n_pixels,m_pixels)),[3,4]),'*k'); hold on;%mean in ROI
    Irrad = GenAirMeasFromStepperParsV2(amp_hat,angle_hat,nSteps,LambdaList); 
    plot((0:nSteps-1)*2*pi/nSteps,mean((Irrad(l,:,n_pixels,m_pixels)),[3,4]));%mean in ROI
    title(['\lambda =  ' num2str(LambdaList(l)),' RMSE=',num2str(res_lambda(l),'%.2e')],'FontSize',30);
    xlabel('PSG Rotation Step [Radians]','FontSize',30)
    ylabel(sprintf('Avg Counts in %dx%d ROI',length(n_pixels),length(m_pixels)),'FontSize',30);
    xlim([0 (nSteps-1)*2*pi/nSteps]);xticks([0 2*pi*(nSteps-1)/nSteps]);xticklabels({'0',sprintf('2\\pi (%d-1)/%d',nSteps,nSteps)});
    a=ylim;yticks([a(1) a(2)]);set(gca,'FontSize',20);
    f=getframe(gcf);
    writeVideo(MRGB,f);close all;
end
%Retardance Magnitude Plot
figure;plot(LambdaList,a_PSG*LambdaList+b_PSG,'*-b'); hold on;
plot(LambdaList,a_PSA*LambdaList+b_PSA,'*-r');
xlabel('\lambda [nm]','FontSize',30);
ylabel('Retardance [rad]','FontSize',30);
h=legend('PSG','PSA');set(h,'FontSize',30);
f=getframe(gcf);
writeVideo(MRGB,f);close all;
close(MRGB);close all;set(0,'DefaultFigureVisible', 'on');
return
