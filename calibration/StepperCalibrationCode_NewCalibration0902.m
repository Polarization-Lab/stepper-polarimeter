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


%Calibration file path
calibrationFilePath = 'D:\Measurements\Dichroic_Analysis\air-01-Sep-2020.h5';
testFilePath = 'D:\Measurements\Dichroic_Analysis\dichroic-45-07-Sep-2020.h5';

for ii = 1:length(LambdaList) %for loading all lambda(s)
%Get h5 reference data
    Lambda = LambdaList(ii);

    refVecs(ii,:) = NormalizedReferenceData(calibrationFilePath,Lambda,nSteps);

%Grab measured data from h5 file
    meandata(ii,:,:,:) = ReadCalImages(calibrationFilePath,Lambda,refVecs,nSteps); %needs to use dark field correction code at some point

end

%Find the max
mx = max(meandata,[],'all');

f = msgbox('Operation Completed');
%%
load handel
sound(y,Fs)

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
lb = [mx*zeros(1,PixelCount*nLambda) -pi*ones(1,nLambda) -pi/2 -pi*ones(1,nLambda) -pi/2 -pi/2];
ub = [3*mx*ones(1,PixelCount*nLambda) pi*ones(1,nLambda) pi/2 pi*ones(1,nLambda) pi/2 pi/2];


%Input vector
% inpvec = [PixelCount, nSteps, ThetaMotorGen, ThetaMotorAna];

%Starting fit guess
% fitguess = [2*mx*ones([1,PixelCount]), 2*pi/3*ones([1,PixelCount]), 2*pi/3*ones([1,PixelCount])];, 2*pi/3, 0, 2*pi/3, 0, 0
fitguess = [mx*ones([1,PixelCount*nLambda]), 2*pi/3*ones(1,nLambda), 0, 2*pi/3*ones(1,nLambda), 0, 0];
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
% sgtitle(['Amplitude ' num2str(Lambda) 'nm'])
%%
%Creating W Matrix
tic
[~,W] = AirCalRoutine_MKedit2(Amp,dg_rad,DelRad_g,da_rad,DelRad_a,LP_Rad,nSteps);
toc
%% Generate Mueller Matrices
% refVecsImgs = NormalizedReferenceData(testFilePath, Lambda, nSteps);
% [imgdata ROIimg] = ReadDataImages(testFilePath, Lambda,refVecsImgs,nSteps);

for ii = 1:length(LambdaList) %for loading all lambda(s)
%Get h5 reference data
    Lambda = LambdaList(ii);
%     refVecs(ii,:) = NormalizedReferenceData(calibrationFilePath,Lambda,nSteps);
    if Lambda == 470
        Lambda = 480;
        ii = ii + 1;
    end
%Grab measured data from h5 file
%     meandata(ii,:,:,:) = ReadCalImages(calibrationFilePath,Lambda,refVecs,nSteps); %needs to use dark field correction code at some point
    imgdata(ii,:,:,:)= ReadDataImagesNoRef(testFilePath,Lambda,nSteps);
end
% save('dichroic_575_45','imgdata','ROIimg');
f = msgbox('Operation Completed');
%%

tic
nLambda=length(LambdaList);
M = zeros(nLambda,16,1500*1500);
for n = 1:nLambda
    Winv = pinv(squeeze(W(n,:,:)));
    M(n,:,:) = Winv*(reshape(squeeze(imgdata(n,:,:,:)),nSteps,1500*1500)./repmat(ampcaldata(n,:),[nSteps 1]));
end
for n = 1:nLambda
mmVecs(n,:,:,:) = reshape(squeeze(M(n,:,:)),16,1500,1500);
end
toc
%% .h5 ROI image
tic
nLambda=length(Lambda);
M = zeros(nLambda,16,10*10);
for n = 1:nLambda
    Winv = pinv(squeeze(W(n,:,:)));
    M(n,:,:) = Winv*(reshape(ROIimg(:,:,:),nSteps,10*10)./repmat(Amp(:)',[nSteps 1]));
end
mmVecsROI = reshape(squeeze(M(1,:,:)),16,10,10);
mmVecs = mmVecsROI;
toc
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
mmNum = [ "00" "01" "02" "03" '10' '11' '12' '13' '20' '21' '22' '23' '30' '31' '32' '33'];
%mx=max(mmVecs(:));
for n = 1:nLambda
    lims = [-abs(max(squeeze(mmVecs(n,:)),[],'all')) abs(max(squeeze(mmVecs(n,:)),[],'all'))];


    for p=1:16
        subplot(4,4,p)
        imshow(squeeze(mmVecs(n,p,:,:)),lims,'Colormap',GWP);
    %     imshow(TotalMM(p),lims,'Colormap',GWP);
        title(strcat('M',mmNum(p)))
    end
    figure(1)
    t = (subplot(4,4,16).Position);
    colorbar('position', [t(1)+t(3) t(2) t(3)/3 t(4)*4.7] );
    sgtitle(['Dichroic 45' char(176) ' ' num2str(LambdaList(n)) 'nm' ])
    saveas(gcf,['FullROILambda ' num2str(LambdaList(n)) '.png']);
    close all;
end
toc
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
figure(n)
t = (subplot(4,4,16).Position);
colorbar('position', [t(1)+t(3) t(2) t(3)/3 t(4)*4.7] );
sgtitle(['Dichroic 45' char(176) ' ' num2str(Lambda) 'nm Norm' ])
toc

%% M_xy setup for easier calculations
Label = ['Dichroic 45' char(176)];
WL = num2str(Lambda);

M_00= squeeze(mmVecs(1,:,:));
M_01= squeeze(mmVecs(2,:,:));
M_02= squeeze(mmVecs(3,:,:));
M_03= squeeze(mmVecs(4,:,:));
M_10= squeeze(mmVecs(5,:,:));
M_11= squeeze(mmVecs(6,:,:));
M_12= squeeze(mmVecs(7,:,:));
M_13= squeeze(mmVecs(8,:,:));
M_20= squeeze(mmVecs(9,:,:));
M_21= squeeze(mmVecs(10,:,:));
M_22= squeeze(mmVecs(11,:,:));
M_23= squeeze(mmVecs(12,:,:));
M_30= squeeze(mmVecs(13,:,:)); 
M_31= squeeze(mmVecs(14,:,:));
M_32= squeeze(mmVecs(15,:,:));
M_33= squeeze(mmVecs(16,:,:));

%% Diattenuation
close all;

%Equation 6.41 page 177 from Russell's book
Dia = sqrt(mmVecs(2,:,:).^2 + mmVecs(3,:,:).^2 + mmVecs(4,:,:).^2)./mmVecs(1,:,:);

%Plot Diattenuation map and histogram
figure(1)
subplot(1,2,1)
imshow(squeeze(Dia(1,:,:)),[0 1],'colormap',parula);colorbar;
title(['2D Diattenuation per pixel ' Label WL]);
subplot(1,2,2)
imhist(squeeze(Dia))
title(['Diattenuation Histogram (' WL ')']);

%% Linear Diattenuation

%Equation 6.49 page 179 from Russell's book
LinDia = sqrt(mmVecs(2,:,:).^2 + mmVecs(3,:,:).^2)./mmVecs(1,:,:);

%Plot LD map and histogram
figure(2)
subplot(1,2,1)
imshow(squeeze(LinDia(1,:,:)),[0 1],'colormap',parula);colorbar;
title(['2D Linear Diatten. per pixel ' Label WL]);
subplot(1,2,2)
imhist(squeeze(LinDia))
title(['Linear Diattenuation Histogram (' WL ')']);
%% Polarizance
figure(3)

%Equation 6.51 page 180 from Russell's book
Pol = sqrt(mmVecs(5,:,:).^2 + mmVecs(9,:,:).^2 + mmVecs(13,:,:).^2)./mmVecs(1,:,:);

%Plot Polarizance map and histogram
subplot(1,2,1)
imshow(squeeze(Pol(1,:,:)),[0 1],'colormap',parula);colorbar;
title(['2D Polarizance per pixel ' Label WL]);

subplot(1,2,2)
imhist(squeeze(Pol))
title(['Polarizance Histogram (' WL ')']);
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
% delta_r = sqrt(del_Q.^2 + del_U.^2 + del_V.^2);
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

%Sum of squares variable
M_ij = zeros(10,10);

%create sum of squares and calculate DI. Equation 6.79 page 195 from Russell's book
for ii = 1:16
    M_ij = squeeze(mmVecs(ii,:,:)).^2 + M_ij;
end
DepolIndex = sqrt(M_ij-squeeze(mmVecs(1,:,:)).^2)./(sqrt(3)*squeeze(mmVecs(1,:,:)));

Mean_DI = mean(DepolIndex,'all')
%Plot DI map
subplot(1,2,1)
imshow(squeeze(DepolIndex(:,:)),[0 1],'colormap',parula);colorbar;
title(['2D Depolarization Index per pixel ' Label WL]);

subplot(1,2,2)
imhist(squeeze(DepolIndex))
title(['Depolarization Index Histogram (' WL ')']);
%% Jones Matrix Translation
%M = U.kronecker(J*,J).U^-1

