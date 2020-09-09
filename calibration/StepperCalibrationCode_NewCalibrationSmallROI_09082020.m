% Stepper Calibration Code
% Authors: Lisa Li, James Heath, Kira Hart
% Last Updated: 2020/08/25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear variables
tic
clear all; close all;clc;
nSteps = 64; %StepperCal = 64, RGB950 = 40
LambdaList = linspace(450,740,30);
nLambda = length(LambdaList);

%image size. Needs to be more robust.
testimsz_1 = 1500;
testimsz_2 = 1500;
ROI = 5;

%Calibration file path
calibrationFilePath = 'Y:\Measurements\Dichroic_Analysis\air-01-Sep-2020.h5';
testFilePath = '\\765-polarizer.optics.arizona.edu\stepper\Stepper Polarimeter\Stepper Software Rebuild\TestData\Dichroic_Analysis\dichroic-45-07-Sep-2020.h5';

for ii = 1:length(LambdaList) %for loading all lambda(s)
%Get h5 reference data
    Lambda = LambdaList(ii);
    refVecs(ii,:) = NormalizedReferenceData(calibrationFilePath,Lambda,nSteps);

%Grab measured data from h5 file
    meandata(ii,:,:,:) = ReadCalImages(calibrationFilePath,Lambda,refVecs,nSteps); %needs to use dark field correction code at some point

end

%Find the max
mx = max(meandata,[],'all');

%%
clc;
%Get amount of pixels
PixelCount = length(meandata(1,1,:));

%meandataVecs = reshape(meandata,1,64,10,10);


%Get a function for fitting the values
%Get spectrum data for "lambda colored" plots
% sRGB = spectrumRGB(Lambda, '1964_FULL');
% cRGB(:) = [sRGB(:,1) sRGB(:,2) sRGB(:,3)]; %Colors for graphing

%AirCalRoutine_MKedit2(Amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,nSteps)
func = @(fitvals, inpvals) AirCalRoutine_MKeditLSQ(fitvals(1:PixelCount*nLambda),fitvals(PixelCount*nLambda+1:PixelCount*nLambda+nLambda),fitvals(PixelCount*nLambda+1+nLambda),fitvals(PixelCount*nLambda+2+nLambda:PixelCount*nLambda+1+2*nLambda),fitvals(PixelCount*nLambda+2+2*nLambda),fitvals(PixelCount*nLambda+3+2*nLambda),inpvals);

%Set bounds and known values, likely not optimal ask Meredith
% lb = [Amp ];
% ub = [3*mx*ones(1,PixelCount) pi*ones(1,PixelCount) pi*ones(1,PixelCount)];
lb = [mx*ones(1,PixelCount*nLambda) -pi*ones(1,nLambda) -pi/2 -pi*ones(1,nLambda) -pi/2 -pi/2];
ub = [3*mx*ones(1,PixelCount*nLambda) pi*ones(1,nLambda) pi/2 pi*ones(1,nLambda) pi/2 pi/2];


%Input vector
% inpvec = [PixelCount, nSteps, ThetaMotorGen, ThetaMotorAna];

%Starting fit guess
% fitguess = [2*mx*ones([1,PixelCount]), 2*pi/3*ones([1,PixelCount]), 2*pi/3*ones([1,PixelCount])];, 2*pi/3, 0, 2*pi/3, 0, 0
fitguess = [2*mx*ones([1,PixelCount*nLambda]), 2*pi/3*ones(1,nLambda), 0, 2*pi/3*ones(1,nLambda), 0, 0];
%%
%Least squares curve fitting and storing into caldata
[fits,resnorm,res] = lsqcurvefit(func, fitguess, nSteps, meandata,lb,ub);
caldata(:)=fits(:);
%%
%RMSE of data
for ii = 1:nLambda
RMSE_var(ii) = std(res(ii,:)); %will need to fix size (1 x nLambda)
end
%Grab caldata for tables, plotting, and psgMM/psaMM
%Create table in radians
Amp = caldata(1:PixelCount*nLambda);
dg_rad = caldata(PixelCount*nLambda+1:PixelCount*nLambda+nLambda); %change range
DelRad_g = caldata(PixelCount*nLambda+1+nLambda);
da_rad = caldata(PixelCount*nLambda+2+nLambda:PixelCount*nLambda+1+2*nLambda); %change range
DelRad_a = caldata(PixelCount*nLambda+2+2*nLambda);
LP_Rad = caldata(PixelCount*nLambda+3+2*nLambda);

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
plot(ThetaMotorGen,mean(squeeze(meandata(13,:,:,:)),[2,3]),'*r'); 
Irrad = AirCalRoutine_MKeditLSQ(Amp,dg_rad,DelRad_g,da_rad,DelRad_a,LP_Rad, 629);

%Plot fit
plot(0:0.01:2*pi,mean(squeeze(Irrad(13,:,:,:)),[2,3]));%,'color',cRGB);
title(['Lambda =  ' num2str(LambdaList(13))]);
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
% sgtitle(['Amplitude ' num2str(Lambda) 'nm'])
%%
%Creating W Matrix
tic
[~,W] = AirCalRoutine_MKedit2(Amp,dg_rad,DelRad_g,da_rad,DelRad_a,LP_Rad,nSteps);
toc
%% Generate Mueller Matrices
Lambda = 450;


refVecsImgs = NormalizedReferenceData(testFilePath, Lambda, nSteps);
% [imgdata, ROIimg] = ReadDataImages(testFilePath, Lambda,refVecsImgs,nSteps);
[imgdata, ROIimg] = ReadDataImagesNoRef(testFilePath, Lambda,nSteps);
% save('dichroic_575_45','imgdata','ROIimg');
f = msgbox('Operation Completed');
%%
load handel
sound(y,Fs)
%%
temp = reshape(Amp(:),30,5,5);

%%

tic
% nLambda=length(Lambda);
M = zeros(nLambda,16,1500*1500);
for n = 1:nLambda
    Winv = pinv(squeeze(W(n,:,:)));
    M(n,:,:) = Winv*(reshape(imgdata(:,:,:),nSteps,1500*1500)./repmat(Amp(:)',[nSteps 1]));
end

toc
%% .h5 ROI image1
tic
nLambda=length(LambdaList);
M = zeros(nLambda,16,ROI*ROI);
for n = 1:nLambda
    Winv = pinv(squeeze(W(n,:,:)));
    M(n,:,:) = Winv*(reshape(squeeze(ROIimg(n,:,:,:)),nSteps,ROI*ROI)./repmat(AmpArray(n,:),[nSteps 1]));
end
toc
%%
for n = 1:nLambda
    mmVecs(n,:,:,:) = reshape(squeeze(M(n,:,:)),16,ROI,ROI);
end
%% Normalize MM
for ii = 1:16
    NormMMVecs(ii,:,:) = squeeze(mmVecs(ii,:,:))./squeeze(mmVecs(1,:,:));
end

%% Create averaged scalar MM (Table form)
tic
for p=1:16
    avgMM(p) = mean(squeeze(mmVecs(p,:,:))/mean(squeeze(mmVecs(1,:,:))),'all');
end
NormAvgMM = reshape(avgMM,[4 4])'

for p=1:16
    TotalMM(p) = mean(squeeze(mmVecs(p,:,:)),'all');
end
TotalMM2 = reshape(TotalMM,[4 4])'
toc

%% Display image MMs
tic
close all;
mkdir(['Y:\Measurements\Dichroic_Analysis\MM_img\dichroic-45-small-ROI-' date]);
%mx=max(mmVecs(:));
mmNum = [ "00" "01" "02" "03" '10' '11' '12' '13' '20' '21' '22' '23' '30' '31' '32' '33'];

for nLambda = 1:30
    lims = [-abs(max(squeeze(mmVecs(nLambda,:)),[],'all')) abs(max(squeeze(mmVecs(nLambda,:)),[],'all'))];


    for p=1:16
        subplot(4,4,p)
        imshow(squeeze(mmVecs(nLambda,p,:,:)),lims,'Colormap',GWP);
    %     imshow(TotalMM(p),lims,'Colormap',GWP);
        title(strcat('M',mmNum(p)))
    end
    figure(1)
    t = (subplot(4,4,16).Position);
    colorbar('position', [t(1)+t(3) t(2) t(3)/3 t(4)*4.7] );
    sgtitle(['Dichroic 45' char(176) ' ' num2str(LambdaList(nLambda)) 'nm' ])
    F(nLambda) = getframe(gcf);
    saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\MM_img\dichroic-45-small-ROI-' date],['Lambda ' num2str(LambdaList(nLambda)) '_' date '.png']));
    close all
end
toc

%%
writerObj = VideoWriter('dichroicSmallROI.avi');
writerObj.FrameRate = 1;

open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj,frame);
end
close(writerObj);



%% Display Normalized MM
tic
close all;
%mx=max(mmVecs(:));

lims = [-abs(max(NormMMVecs,[],'all')) abs(max(NormMMVecs,[],'all'))];

mmNum = [ "00" "01" "02" "03" '10' '11' '12' '13' '20' '21' '22' '23' '30' '31' '32' '33'];
for p=1:16
    subplot(4,4,p)
    imshow(squeeze(NormMMVecs(p,:,:)),lims,'Colormap',GWP);
%     imshow(TotalMM(p),lims,'Colormap',GWP);
    title(strcat('M',mmNum(p)))
end
t = (subplot(4,4,16).Position);
colorbar('position', [t(1)+t(3) t(2) t(3)/3 t(4)*4.7] );
sgtitle(['Dichroic 45' char(176) ' ' num2str(Lambda) 'nm Norm' ])
toc

%% M_xy setup for easier calculations
for n = 1:nLambda
    Label = ['Dichroic 45' char(176)];
    WL = num2str(Lambda);

    M_00(n,:,:)= squeeze(mmVecs(n,1,:,:));
    M_01(n,:,:)= squeeze(mmVecs(n,2,:,:));
    M_02(n,:,:)= squeeze(mmVecs(n,3,:,:));
    M_03(n,:,:)= squeeze(mmVecs(n,4,:,:));
    M_10(n,:,:)= squeeze(mmVecs(n,5,:,:));
    M_11(n,:,:)= squeeze(mmVecs(n,6,:,:));
    M_12(n,:,:)= squeeze(mmVecs(n,7,:,:));
    M_13(n,:,:)= squeeze(mmVecs(n,8,:,:));
    M_20(n,:,:)= squeeze(mmVecs(n,9,:,:));
    M_21(n,:,:)= squeeze(mmVecs(n,10,:,:));
    M_22(n,:,:)= squeeze(mmVecs(n,11,:,:));
    M_23(n,:,:)= squeeze(mmVecs(n,12,:,:));
    M_30(n,:,:)= squeeze(mmVecs(n,13,:,:)); 
    M_31(n,:,:)= squeeze(mmVecs(n,14,:,:));
    M_32(n,:,:)= squeeze(mmVecs(n,15,:,:));
    M_33(n,:,:)= squeeze(mmVecs(n,16,:,:));
end

M_00(M_00 == 0) = NaN;
%% Average Transmission
mkdir(['Y:\Measurements\Dichroic_Analysis\MM_img\dichroic-45-small-ROI-' date]);
meanTransmission = mean(M_00,[2 3]);
plot(LambdaList,meanTransmission,'k-*')
title('Avg Transmission over \lambda')
xlabel('\lambda (nm)')
ylabel('Transmission')
saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\MM_img\dichroic-45-small-ROI-' date],'Avg_Transmission_Lambda.png'));
%% Diattenuation
close all;
mkdir(['Y:\Measurements\Dichroic_Analysis\Diattenuation\dichroic-45-small-ROI-' date]);

for n = 1:nLambda
%Equation 6.41 page 177 from Russell's book
% Dia = sqrt(mmVecs(2,:,:).^2 + mmVecs(3,:,:).^2 + mmVecs(4,:,:).^2)./mmVecs(1,:,:);
Dia(n,:,:) = sqrt(M_01(n,:,:).^2 + M_02(n,:,:).^2 + M_03(n,:,:).^2)./M_00(n,:,:);
%Plot Diattenuation map and histogram
figure(1)
subplot(1,2,1)
imshow(squeeze(Dia(n,:,:)),'colormap',parula);colorbar;
title(['2D Diattenuation per pixel ' Label num2str(LambdaList(n))]);
subplot(1,2,2)
histogram(squeeze(Dia(n,:,:)))
title(['Diattenuation Histogram (' num2str(LambdaList(n)) ')']);
saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\Diattenuation\dichroic-45-small-ROI-' date],['Lambda ' num2str(LambdaList(n)) '_' date '.png']));
close all;
end
%%
avgDia = mean(Dia,[2 3]);
plot(LambdaList,avgDia,'b*-');
title('Avg Diattenuation over \lambda')
xlabel('\lambda (nm)')
ylabel('Diattenuation')
saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\Diattenuation\dichroic-45-small-ROI-' date],['Avg_Dia_Lambda.png']));
%% Linear Diattenuation
mkdir(['Y:\Measurements\Dichroic_Analysis\Linear_Diattenuation\dichroic-45-small-ROI-' date]);

edges = linspace(-1.5,1.5,50);

for n = 1:nLambda
    %Equation 6.49 page 179 from Russell's book
    LinDia(n,:,:) = sqrt(M_01(n,:,:).^2 + M_02(n,:,:).^2)./M_00(n,:,:);
    %Plot LD map and histogram
    figure(2)
    subplot(1,2,1)
    imshow(squeeze(LinDia(n,:,:)),[0 1],'colormap',parula);colorbar;
    title(['2D Linear Diatten. per pixel ' Label num2str(LambdaList(n))]);
    subplot(1,2,2)
    histogram(squeeze(LinDia(n,:,:)))
    title(['Linear Diattenuation Histogram (' num2str(LambdaList(n)) ')']);
    saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\Linear_Diattenuation\dichroic-45-small-ROI-' date],['LD_Lambda_' num2str(LambdaList(n)) '_' date '.png']));
    close all;
end
%%
avgLD = mean(LinDia,[2 3]);
plot(LambdaList,avgLD,'r*-');
title('Avg Linear Diattenuation over \lambda')
xlabel('\lambda (nm)')
ylabel('Linear Diattenuation')
saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\Linear_Diattenuation\dichroic-45-small-ROI-' date],['Avg_LinDia_Lambda.png']));
%% Polarizance
mkdir(['Y:\Measurements\Dichroic_Analysis\Polarizance\dichroic-45-small-ROI-' date]);
for n = 1:nLambda
    %Equation 6.51 page 180 from Russell's book
    % Pol = sqrt(mmVecs(5,:,:).^2 + mmVecs(9,:,:).^2 + mmVecs(13,:,:).^2)./mmVecs(1,:,:);
    Pol(n,:,:) = sqrt(M_10(n,:,:).^2 + M_20(n,:,:).^2 + M_30(n,:,:).^2)./M_00(n,:,:);
    %Plot Polarizance map and histogram
    subplot(1,2,1)
    imshow(squeeze(Pol(n,:,:)),[0 1],'colormap',parula);colorbar;
    title(['2D Polarizance per pixel ' Label num2str(LambdaList(n))]);

    subplot(1,2,2)
    histogram(squeeze(Pol(n,:,:)))
%     imhist(squeeze(Pol))
    title(['Polarizance Histogram (' num2str(LambdaList(n)) ')']);
    saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\Polarizance\dichroic-45-small-ROI-' date],['Pol_Lambda_' num2str(LambdaList(n)) '_' date '.png']));
    close all;
end
%%
avgPol = mean(Pol,[2 3]);
plot(LambdaList,avgPol,'g*-');
title('Avg Polarizance over \lambda')
xlabel('\lambda (nm)')
ylabel('Polarizance')
saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\Polarizance\dichroic-45-small-ROI-' date],['Avg_Pol_Lambda.png']));
%% Retardation Map

tic
close all;

delta = acos(((M_00 + M_11 + M_22 + M_33)/2)-1);
del_H = delta/(2*sin(delta))*(M_23 - M_32);
del_45 = delta/(2*sin(delta))*(M_31 - M_13);
del_R = delta/(2*sin(delta))*(M_12 - M_21);

del_H = 0.5*(M_23 - M_32);
del_45 = 0.5*(M_31 - M_13);
del_R = 0.5*(M_12 - M_21);


%Equation(s) 6.32 from page 173 of Russell's book. Only works for halfwave
%retarder!
% del_H = pi*(sqrt((squeeze(mmVecs(5,:,:))+1)./2));
% del_H = pi*(sqrt((M_10 + 1)./2));
% % del_45 = pi*(sign(squeeze(mmVecs(6,:,:)).*sqrt((squeeze(mmVecs(10,:,:))+1)/2)));
% del_45 = pi*(sign(M_11.*sqrt((M_21 + 1)/2)));
% % del_R = pi*(sign(squeeze(mmVecs(7,:,:)).*sqrt((squeeze(mmVecs(16,:,:))+1)/2)));
% del_R = pi*(sign(M_12.*sqrt((M_33+1)/2)));
% 
% %Equation 6.27 from page 171 of Russell's book
delta_check = sqrt(del_H.^2 + del_45.^2 + del_R.^2);
edges = linspace(-1,1,200);
%Plot retardance
a=subplot(2,2,1)
histogram(delta_check,'BinEdges',edges);
title('\delta')
b=subplot(2,2,2)
histogram(del_H,'BinEdges',edges);
title('\delta_H')
c=subplot(2,2,3)
histogram(del_45,'BinEdges',edges);
title('\delta_{45}')
loc = subplot(2,2,4)
d = get(loc, 'Position');
histogram(del_R,'BinEdges',edges);
title('\delta_R')
% sgtitle([Label WL ' Retardance Maps'])
toc

%% Depolarization Index
close all;
mkdir(['Y:\Measurements\Dichroic_Analysis\DepolIndex\dichroic-45-small-ROI-' date]);
%Sum of squares variable
M_ij = zeros(nLambda,5,5);

for n = 1:nLambda
%create sum of squares and calculate DI. Equation 6.79 page 195 from Russell's book
    for ii = 1:16
        M_ij(n,:,:) = squeeze(mmVecs(n,ii,:,:)).^2 + squeeze(M_ij(n,:,:));
    end
end

for n= 1:nLambda
    % DepolIndex = sqrt(M_ij-squeeze(mmVecs(1,:,:)).^2)./(sqrt(3)*squeeze(mmVecs(1,:,:)));
    for ii = 1:16
        DepolIndex(n,:,:) = sqrt(squeeze(M_ij(n,:,:))-squeeze(M_00(n,:,:)).^2) ./ (sqrt(3)*squeeze(M_00(n,:,:)));
    end
    
    figure(1)
    %Plot DI map
    subplot(1,2,1)
    imshow(squeeze(DepolIndex(n,:,:)),[0 1],'colormap',parula);colorbar;
    title(['2D Depolarization Index per pixel ' Label num2str(LambdaList(n))]);

    subplot(1,2,2)
    histogram(squeeze(DepolIndex(n,:,:)))
    title(['Depolarization Index Histogram (' num2str(LambdaList(n)) ')']);
    saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\DepolIndex\dichroic-45-small-ROI-' date],['DI_Lambda_' num2str(LambdaList(n)) '_' date '.png']));
    close all;
end

%%
avgDI = mean(DepolIndex,[2 3]);
plot(LambdaList,avgDI,'m*-');
title('Depolarization Index over \lambda')
xlabel('\lambda (nm)')
ylabel('Depol. Index')
saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\DepolIndex\dichroic-45-small-ROI-' date],'Avg_DI_Lambda.png'));

