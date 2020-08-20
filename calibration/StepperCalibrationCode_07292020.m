% Stepper Calibration Code
% Version 1.0
% Original Authors: Lisa Li, James Heath
% Last Updated: 2020/07/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear variables
clear all; close all;clc;
i = sqrt(-1);
tic
test = 0;

%number of calibration trials. If there is only one, set = 1, etc
numCal = 1;


%kk cycles over the # trials within in the 'Lambda'nm_'date'_Exp-'second's
%folder. 
for kk = 1:numCal
    %setup calibration and test folders
    %For multiple trials: ['filepath\' num2str(kk)];
    %For single trial: 'filepath';
%     calibrationFolder = ['C:\Users\Jake.000\School\UASAL\Calibration\550nm\550nm_20200614_Exp-0.09s\' num2str(kk)]; %Multiple trials
    calibrationFolder = 'C:\Users\Jake.000\School\UASAL\Calibration\550nm\AirCal_20200715_550_Exp-0.085s'; %Single trial
    
    %Steps for each retarder measurement sequence
    nSteps = 64; %StepperCal = 64, RGB950 = 40

    %Open info file and display
    fid = fopen(strcat(calibrationFolder,'\info.txt'));
    info = textscan(fid, '%s', 'delimiter', '\n\r', 'whitespace', '');
    fclose(fid);
    info{:};
    temp = cat(1,info{:});
    AllInfo = strjoin(temp);
    AllInfo = convertCharsToStrings(AllInfo);

    %Grab starting, ending and lambda step from info.txt
    LambdaStart = str2num(string(extractBetween(AllInfo,"scan:","nm")));
    LambdaEnd = str2num(string(extractBetween(AllInfo,"nm --","nm")));
    delLambda = str2num(string(extractBetween(AllInfo,"size","nm")));

    %setup steps for for loop
    LambdaList = LambdaStart:delLambda:LambdaEnd;

    %create colors for graphing gorgeously in CIE 1964 format
    sRGB = spectrumRGB(LambdaList, '1964_FULL');

    %Setup Generator angles
    ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 
    %Setting up analyzer and image size
    ThetaMotorAna = 4.9*ThetaMotorGen;

    %import ref_readings dat file
    ref = importdata(strcat(calibrationFolder,'\ref_readings.dat'));

    hold off; 
    mx = 0;
    for jj = 1:length(LambdaList)

        %reads in binary file with: 'b' for big endian (default is little),
        %'float32' (default double)
%         inp =  fread(fopen(strcat(calibrationFolder,'\imageData\images',int2str(jj-1),'.bin'),'r'),'float32','b');
        inp = fopen(strcat(calibrationFolder,'\imageData\images',int2str(jj-1),'.bin'),'r');
        calimsz_1=fread(inp,1,'int32','b');
        calimsz_2=fread(inp,1,'int32','b');	
        inp = fread(fopen(strcat(calibrationFolder,'\imageData\images',int2str(jj-1),'.bin'),'r'),'single','b');

        for ii = 1:nSteps       
            %Grab binary image measurement, regularize it, and store the wanted
            %ROI into measdata variable
            measimg(ii,:,:) = reshape(inp(2*ii+(ii-1)*(calimsz_1*calimsz_2):2*ii+ii*(calimsz_1*calimsz_2)-1),[calimsz_1,calimsz_2])/ref(ii,jj);
            temp = measimg(ii,1295:1305,1015:1025); %measimg(ii,:,:);
            measdata(jj,ii) = mean(temp(:));
            new_mx = max(measimg,[],'all');
            if mx <= new_mx
                mx = new_mx;
            end
            test = test +1;
            
        end
        
        %Get a function for fitting the values
        %AirCalRoutine(amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,PSGangles,PSAangles)
        func = @(fitvals, inpvals) AirCalRoutine_MKedit(fitvals(1),fitvals(2),fitvals(3),fitvals(4),fitvals(5),fitvals(6),inpvals(1:nSteps),inpvals(65:128));
        cRGB(:) = [sRGB(:,1) sRGB(:,2) sRGB(:,3)]; %Colors for graphing

        %Set upper (ub) and lower (lb) boundaries for fitting values
        lb = [mx -pi -pi/2 -pi -pi/2 -pi/2];
        ub = [3*mx pi pi/2 pi pi/2 pi/2];

        %options = optimset('MaxFunEvals',1000000); %Indicates how many iterations
        opts = optimoptions('lsqcurvefit', 'Diagnostics','on', 'Display','iter-detailed');

        %Least squares curve fitting and storing into caldata
        [fits,resnorm(jj),res] = lsqcurvefit(func,[3*mx, 2*pi/3, 0, 2*pi/3, 0, 0], [ThetaMotorGen ThetaMotorAna], squeeze(measdata(jj,:)),lb,ub,opts);
        caldata(jj,:)=fits(:);

        [Irrad,~] = AirCalRoutine_MKedit(caldata(jj,1),caldata(jj,2),caldata(jj,3),caldata(jj,4),caldata(jj,5),caldata(jj,6), 0:0.01:2*pi, 4.9*(0:0.01:2*pi));
        %Plot each wavelength's fit
        plot(ThetaMotorGen,measdata(jj,:),'*r'); 
        hold on;
        plot(0:0.01:2*pi,Irrad,'color',cRGB);
        title(['Lambda =  ' num2str(LambdaList(jj))]);
        xlabel('PSG Rotation (Radians)')
        ylabel('Camera Counts')

    end
end
toc
%%

%RMSE of data
RMSE_var = std(res(:));

%Create table in radians
Lambda = LambdaList';
Amp = caldata(:,1);
dg_rad = caldata(:,2);
DelRad_g = caldata(:,3);
da_rad = caldata(:,4);
DelRad_a = caldata(:,5);
LP_Rad = caldata(:,6);
RMSE = RMSE_var;


table(Lambda,Amp, dg_rad,DelRad_g,da_rad,DelRad_a,LP_Rad, RMSE)

%Create table in degrees
dg_deg = dg_rad.*180/pi;
DelDeg_g = DelRad_g.*180/pi;
da_deg = da_rad.*180/pi;
DelDeg_a = DelRad_a.*180/pi;
LP_Deg = LP_Rad.*180/pi;
RMSE = RMSE_var;
PSG_Theta = DelDeg_g;
PSA_Theta = DelDeg_a;
table(Lambda,Amp, dg_deg,PSG_Theta,da_deg,PSA_Theta,LP_Deg, RMSE)

%%
%Creating W Matrix

tic
%Setup
temp = [0:360/nSteps:360].*pi/180;
psg_Rad = temp(1:nSteps);
psa_Rad = 4.9.*psg_Rad;

[~,W] = AirCalRoutine_MKedit(Amp,dg_rad,0,da_rad,0,LP_Rad,psg_Rad,psa_Rad);

W_Invs = pinv(W);
toc
%% Generate Mueller Matrices
clear testimg;
%Open test measurements
testFolder = 'C:\Users\Jake.000\School\UASAL\Calibration\550nm\VVR_20200624_550_Exp-0.09s';%'C:\Users\Jake.000\School\UASAL\Calibration\LP_Orientation-A\LP-A_20200622_550nm_Exp-0.09s';%'C:\Users\Jake.000\School\UASAL\Calibration\LCP_Orientation-A\LCP-A_20200623_450nm_Exp-0.28s';%'C:\Users\Jake.000\School\UASAL\Calibration\LP-Horiz_20200526_450nm-750nm-100nm'; %"\\765-Polarizer.optics.arizona.edu\StepperData\Measurements\Summer-2020_Calibration\air_20200522_450nm-750nm-100nm_test";

fid = fopen(strcat(testFolder,'\info.txt'));
rawinfo = textscan(fid, '%s', 'delimiter', '\n\r', 'whitespace', '');
fclose(fid);
rawinfo{:}

temp = cat(1,rawinfo{:});
AllInfo = strjoin(temp);
AllInfo = convertCharsToStrings(AllInfo);

%Grab starting, ending and lambda step from info.txt
LambdaStart = str2num(string(extractBetween(AllInfo,"scan:","nm")));
LambdaEnd = str2num(string(extractBetween(AllInfo,"nm --","nm")));
delLambda = str2num(string(extractBetween(AllInfo,"size","nm")));

ref = importdata(strcat(testFolder,'\ref_readings.dat'));
LambdaList = LambdaStart:delLambda:LambdaEnd;
nWaves = length(LambdaList);
for jj = 1:nWaves
    %reads in binary file with: 'b' for big endian (default is little),
    %'float32' (default double)
     inp = fopen(strcat(testFolder,'\imageData\images',int2str(jj-1),'.bin'),'r');
    testimsz_1=fread(inp,1,'int32','b');
	testimsz_2=fread(inp,1,'int32','b');	
    inp = fread(fopen(strcat(testFolder,'\imageData\images',int2str(jj-1),'.bin'),'r'),'single','b');
    
    for ii = 1:nSteps       
        %Grab binary image measurement, regularize it, and store the wanted
        %ROI into measdata variable
        testimg(ii, :, :) = reshape(inp(2*ii+(ii-1)*(testimsz_1*testimsz_2):2*ii+ii*(testimsz_1*testimsz_2)-1),[testimsz_1,testimsz_2])/ref(ii,jj);
    end
end

%% Create Mueller Matrix

mmVecs = W_Invs*reshape(testimg,64,testimsz_1*testimsz_2);
mmVecs = reshape(mmVecs,16,testimsz_1,testimsz_2);

%% Scalar Average Mueller Matrix

for p=1:16
    avgMM(p) = mean(squeeze(mmVecs(p,:,:))/mean(squeeze(mmVecs(1,:,:))),'all');
end
avgMM = reshape(avgMM,[4 4])'

for p=1:16
    TotalMM(p) = mean(squeeze(mmVecs(p,:,:)),'all');
end
TotalMM = reshape(TotalMM,[4 4])'


%% Display MMs
close all;

lims = [-abs(max(mmVecs,[],'all')) abs(max(mmVecs,[],'all'))];

mmNum = [ "00" "01" "02" "03" '10' '11' '12' '13' '20' '21' '22' '23' '30' '31' '32' '33'];
for p=1:16
    subplot(4,4,p)
    imshow(squeeze(mmVecs(p,:,:)),lims,'Colormap',GWP);
    title(strcat('M',mmNum(p)))
end
t = (subplot(4,4,16).Position);
colorbar('position', [t(1)+t(3) t(2) t(3)/3 t(4)*4.7] );
sgtitle('VVR2 550nm')

%% Diattenuation
close all;

%Equation 6.41 page 177 from Russell's book
Dia = sqrt(mmVecs(2,:,:).^2 + mmVecs(3,:,:).^2 + mmVecs(4,:,:).^2)./mmVecs(1,:,:);

%Plot Diattenuation map and histogram
figure(1)
subplot(1,2,1)
imshow(squeeze(Dia(1,:,:)),[0 1],'colormap',parula);colorbar;
title('2D Diattenuation per pixel Dichroic(0,0) 550nm');
subplot(1,2,2)
imhist(squeeze(Dia))
title('Diattenuation Histogram (750nm)');

%% Linear Diattenuation

%Equation 6.49 page 179 from Russell's book
LinDia = sqrt(mmVecs(2,:,:).^2 + mmVecs(3,:,:).^2)./mmVecs(1,:,:);

%Plot LD map and histogram
figure(2)
subplot(1,2,1)
imshow(squeeze(LinDia(1,:,:)),[0 1],'colormap',parula);colorbar;
title('2D Linear Diatten. per pixel Dichroic(45,0) 750nm');
subplot(1,2,2)
imhist(squeeze(LinDia))
title('Linear Diattenuation Histogram (750nm)');
%% Polarizance
figure(3)

%Equation 6.51 page 180 from Russell's book
Pol = sqrt(mmVecs(5,:,:).^2 + mmVecs(9,:,:).^2 + mmVecs(13,:,:).^2)./mmVecs(1,:,:);

%Plot Polarizance map and histogram
subplot(1,2,1)
imshow(squeeze(Pol(1,:,:)),[0 1],'colormap',parula);colorbar;
title('2D Polarizance per pixel Dichroic(45,0) 750nm');

subplot(1,2,2)
imhist(squeeze(Pol))
title('Polarizance Histogram (750nm)');
%% Retardation Map

tic
close all;

%Equation(s) 6.32 from page 173 of Russell's book. Only works for halfwave
%retarder!
del_Q = pi*(sqrt((squeeze(mmVecs(5,:,:))+1)./2));
del_U = pi*(sign(squeeze(mmVecs(6,:,:)).*sqrt((squeeze(mmVecs(10,:,:))+1)/2)));
del_V = pi*(sign(squeeze(mmVecs(7,:,:)).*sqrt((squeeze(mmVecs(16,:,:))+1)/2)));

%Equation 6.27 from page 171 of Russell's book
delta_r = sqrt(del_Q.^2 + del_U.^2 + del_V.^2);

%Plot retardance
a=subplot(2,2,1)
imshow(delta_r,[-pi pi],'colormap',parula); colorbar;
title('\delta')
b=subplot(2,2,2)
imshow(del_Q,[-pi pi],'colormap',parula); colorbar;
title('\delta_H')
c=subplot(2,2,3)
imshow(del_U,[-pi pi],'colormap',parula); colorbar;
title('\delta_{45}')
loc = subplot(2,2,4)
d = get(loc, 'Position');
imshow(del_V,[-pi pi],'colormap',parula); colorbar;
title('\delta_R')
sgtitle('Dichroic(45,0) 750nm Retardance Maps')
toc
%% Depolarization Index
close all;

%Sum of squares variable
M_ij = zeros(testimsz_1,testimsz_2);

%create sum of squares and calculate DI. Equation 6.79 page 195 from Russell's book
for ii = 1:16
    M_ij = squeeze(mmVecs(ii,:,:)).^2 + M_ij;
end
DepolIndex = sqrt(M_ij-squeeze(mmVecs(1,:,:)).^2)./(sqrt(3)*squeeze(mmVecs(1,:,:)));

%Plot DI map
imshow(squeeze(DepolIndex(:,:)),[0 1],'colormap',parula);colorbar;
title('2D Depolarization Index per pixel Dichroic(45,0) 750nm');

%% Mueller to Jones from Russell's book (6.12.3)

%Setting up easier variables for Mueller to Jones calculations
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

p = zeros(4,testimsz_1,testimsz_2); %Magnitude variable
phi = zeros(4,testimsz_1,testimsz_2); %Phase variable
JM = zeros(4,testimsz_1,testimsz_2); %Jones xx,xy,yx,yy

%Find amplitudes using Equation(s) 6.111 page 208 from Russell's book
% p(1,:,:) = sqrt(0.5*(squeeze(mmVecs(1,:,:))+squeeze(mmVecs(2,:,:))+squeeze(mmVecs(5,:,:))+squeeze(mmVecs(6,:,:))));
% p(2,:,:) = sqrt(0.5*(squeeze(mmVecs(1,:,:))-squeeze(mmVecs(2,:,:))+squeeze(mmVecs(5,:,:))-squeeze(mmVecs(6,:,:))));
% p(3,:,:) = sqrt(0.5*(squeeze(mmVecs(1,:,:))+squeeze(mmVecs(2,:,:))-squeeze(mmVecs(5,:,:))-squeeze(mmVecs(6,:,:))));
% p(4,:,:) = sqrt(0.5*(squeeze(mmVecs(1,:,:))-squeeze(mmVecs(2,:,:))-squeeze(mmVecs(5,:,:))+squeeze(mmVecs(6,:,:))));
p(1,:,:) = sqrt(0.5*(M_00 + M_01 + M_10 + M_11));
p(2,:,:) = sqrt(0.5*(M_00 - M_01 + M_10 - M_11));
p(3,:,:) = sqrt(0.5*(M_00 + M_01 - M_10 - M_11));
p(4,:,:) = sqrt(0.5*(M_00 - M_01 - M_10 + M_11));

%Relative phase Equation(s) 6.112 page 208 from Russell's book

% phi(2,:,:) = atan2((-squeeze(mmVecs(4,:,:))-squeeze(mmVecs(8,:,:))),(squeeze(mmVecs(3,:,:)) + squeeze(mmVecs(7,:,:))));
% phi(3,:,:) = atan2((squeeze(mmVecs(13,:,:))+squeeze(mmVecs(14,:,:))),(squeeze(mmVecs(9,:,:)) + squeeze(mmVecs(10,:,:))));
% phi(4,:,:) = atan2((squeeze(mmVecs(15,:,:))-squeeze(mmVecs(12,:,:))),(squeeze(mmVecs(11,:,:)) + squeeze(mmVecs(16,:,:))));
phi(2,:,:) = atan2((M_02 + M_12),(-M_03 - M_13));
phi(3,:,:) = atan2((M_20 + M_21),(M_30 + M_31));
phi(4,:,:) = atan2((M_22 + M_33),(M_32 - M_23));

%Jones Matrix subtracting phi_xx from equation 6.110 page 207
JM(1,:,:) = p(1,:,:).*exp(-i*-phi(4,:,:));
for ii = 2:4
    JM(ii,:,:) = p(ii,:,:).*exp(-i*phi(ii,:,:));
end


%%
close all;
figure(1)
subplot(2,4,1)
imshow(squeeze(real(JM(1,:,:))),[],'colormap',GWP);colorbar;
title('J_{xx}')
subplot(2,4,2)
imshow(squeeze(real(JM(2,:,:))),[],'colormap',GWP);colorbar;
title('J_{xy}')
subplot(2,4,5)
imshow(squeeze(real(JM(3,:,:))),[],'colormap',GWP);colorbar;
title('J_{yx}')
subplot(2,4,6)
imshow(squeeze(real(JM(4,:,:))),[],'colormap',GWP);colorbar;
title('J_{yy}')
subplot(2,4,3)
imshow(squeeze(angle(JM(1,:,:))),[],'colormap',GWP);colorbar;
title('\phi_{xx}')
subplot(2,4,4)
imshow(squeeze(angle(JM(2,:,:))),[],'colormap',GWP);colorbar;
title('\phi_{xy}')
subplot(2,4,7)
imshow(squeeze(angle(JM(3,:,:))),[],'colormap',GWP);colorbar;
title('\phi_{yx}')
subplot(2,4,8)
imshow(squeeze(angle(JM(4,:,:))),[],'colormap',GWP);colorbar;
title('\phi_{yy}')

