function [I,C,S]=ViewMMh5( fn_sample )

%This function accepts the file path to an h5 file produces by UA stepper
%polarimeter: https://github.com/Polarization-Lab/stepper-polarimeter/wiki
%Produces two PDF files of graphical output.

%OUTPUT VARIABLES 
%I: nLambdax16xROI_widthxROI_width the un-normalized Mueller matrix

%Example paths
% fn_sample1='/Volumes/StepperData/Measurements/Waveguides/Option2/SpaceMeasurements/Run3_09102021/OnAxis_Option2_530nm-10-Sep-2021.h5';
% fn_sample2='/Volumes/StepperData/Measurements/Waveguides/Option2/SpaceMeasurements/Off_Axis/-2deg/neg2deg_OffAxis_530nm-10-Sep-2021.h5';
% fn_sample3='/Volumes/StepperData/Measurements/Waveguides/Option2/SpaceMeasurements/Off_Axis/+2deg/pos2deg_OffAxis_530nm-10-Sep-2021.h5'

waves=h5info(fn_sample,'/images/');
waves={waves.Groups.Name};
LambdaList=zeros(1,length(waves));
nLambda=length(LambdaList);

for l=1:nLambda
    w=char(waves(l));
    LambdaList(l)=str2double(w((end-2):end));%grab last 3 characters for wavelength
end

%Select from index of Mueller matrices - minus 5 pickest latest version
info = h5info(fn_sample);
index=length({info.Groups.Name})-5;%Index keeps track on number of times intermediate images are processed into Mueller reconstruction. Indexing allows different calibration files or ROIs to be used on intermediate images.
index=4; %indexing for pos2deg_OffAxis_530nm-10-Sep-2021.h5 & neg2deg_OffAxis_530nm-10-Sep-2021.h
index=0; %indexing for OnAxis_Option2_530nm-10-Sep-2021.h5


rect=h5readatt(fn_sample,strcat('/Mueller',num2str(index),'/'),'ROI');
nPixels=(rect(4)+1)^2; 

I=zeros(nLambda,16,nPixels);
C=zeros(nLambda,nPixels,4,4);
S=zeros(nLambda,nPixels,4);
ID=zeros(4,4);
ID(1)=1;
Depol=zeros(nLambda,nPixels);
Pol_mag=zeros(nLambda,nPixels);
Di_mag=zeros(nLambda,nPixels);



for n=1:nLambda
     I(n,:,:)=h5read(fn_sample,strcat('/Mueller',num2str(index),'/wave',num2str(LambdaList(n))));
     %I = reshape(I,[4 4 sqrt(nPixels) sqrt(nPixels)]);
     %inverse fesnel mirror MM, uncoated aluminum at 45 degrees, 532 nm 
     %F=[1.089729264540746 0.03105722428292043 0 0;0.03105722428292043 1.0897292645407457 0 0;0 0 -1.064747642529994 -0.22990819155300934; 0 0 0.22990819155300934 -1.064747642529994];
     %I=pagemtimes(I,F);
     %I = reshape(I,[1 16 nPixels]);
     C(n,:,:,:)=h5read(fn_sample,strcat('/Real_CoherencySingularVectors',num2str(index),'/wave',num2str(LambdaList(n))));
     C(n,:,:,:)=squeeze(C(n,:,:,:))+1i*h5read(fn_sample,strcat('/Imag_CoherencySingularVectors',num2str(index),'/wave',num2str(LambdaList(n))));
     S(n,:,:)= h5read(fn_sample,strcat('/CoherencySingularValues',num2str(index),'/wave',num2str(LambdaList(n))));
     Inorm=squeeze(I(n,:,:)./repmat(I(n,1,:),[1 16]));%normalized MM
	Depol(n,:)=sqrt(sum((Inorm-repmat(ID(:),[1 nPixels])).^2))/sqrt(3);%L2 norm between ideal depolarizer and normalized MM
     Pol_mag(n,:)=sqrt(Inorm(2,:).^2+Inorm(3,:).^2+Inorm(4,:).^2);
     Di_mag(n,:)=sqrt(Inorm(5,:).^2+Inorm(9,:).^2+Inorm(13,:).^2);
end

%Plot truncates range to [0,1] for Polarizance Magnitude, Diattenuation Magnitude, and Depolarization Index
%Plotting example for n=nLambda

figure;
subplot(1,3,1);
imagesc(reshape(Pol_mag,sqrt(nPixels),sqrt(nPixels)),[0 1]);axis square; axis off;colormap jet;colorbar;
title('|Polarizance|','FontSize',15);
subplot(1,3,2);
imagesc(reshape(Di_mag,sqrt(nPixels),sqrt(nPixels)),[0 1]);axis square; axis off;colormap jet;colorbar;
title('|Diattenuation|','FontSize',15);
subplot(1,3,3);
imagesc(reshape(Depol,sqrt(nPixels),sqrt(nPixels)),[0 1]);axis square; axis off;colormap jet;colorbar;
title('Depolarization Index','FontSize',12);
%subplot(2,2,4);
%imagesc(reshape(I(1,1,:),sqrt(nPixels),sqrt(nPixels)));axis square; axis off;colormap jet;colorbar;
%title('M_{00}','FontSize',25);
w=char(fn_sample);
cmd=sprintf('print -dpdf %s',['''',w(1:end-3),sprintf('_MMpars%d_Lambda%d.pdf',index,LambdaList(n)),'''']);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gcf,'renderer','painters');drawnow;
eval(cmd); close all;

%Plot of Mueller Matrix

cmap = zeros(256,3);
% Red
cmap(128:184,1)= linspace(0,1,57);
cmap(184:256,1)= 1;
% Green
cmap(1:96,2) = linspace(1,0,96);
cmap(160:256,2)= linspace(0,1,97);
% Blue
cmap(1:72,3) = 1;
cmap(72:128,3)= linspace(1,0,57);

I=reshape(squeeze(I(n,:,:)),16,sqrt(nPixels),sqrt(nPixels));
figure;

    tempImg=[];
    for nn=1:4
        tempRow=[];
        for zz=1:4
            tempSubImg=squeeze(I(zz+(nn-1)*4,:,:));
            tempSubImg=[min(I(:))*zeros(sqrt(nPixels),1)';tempSubImg;min(I(:))*zeros(sqrt(nPixels),1)'];%Put borders on image
            tempSubImg=[min(I(:))*zeros(sqrt(nPixels)+2,1),tempSubImg,min(I(:))*zeros(sqrt(nPixels)+2,1)];
            tempRow=[tempRow;tempSubImg];
        end
        tempImg=[tempImg,tempRow];
	end
    t=prctile((abs(I(:))),90); %this adjusts the colorbar so highest value is in the 90th percentile. this helps with looking at more detail when max value doesn't appear in image that much 
    imagesc(tempImg,[-t,t]);axis square; axis off; colormap(cmap);colorbar;
    
cmd=sprintf('print -dpdf %s',['''',w(1:end-3),sprintf('_MM%d_Lambda%d.pdf',index,LambdaList(n)),'''']);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gcf,'renderer','painters');drawnow;
eval(cmd); close all;
    
%% 

figure;
S=reshape(S,sqrt(nPixels),sqrt(nPixels),4);
for nn=1:4
    subplot(2,2,nn)
    imagesc(S(:,:,nn));c=colorbar;axis off;axis equal;set(c,'FontSize',15);title(sprintf('\\lambda_{%d}',nn),'FontSize',25);
end
cmd=sprintf('print -dpdf %s',['''',w(1:end-3),sprintf('_EigenSpectrum%d_Lambda%d.pdf',index,LambdaList(n)),'''']);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gcf,'renderer','painters');drawnow;
eval(cmd); close all;

%Each Mueller-Jones matrix
%Analysis for cohereny matrix
U = [1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 sqrt(-1) -sqrt(-1) 0];
PI = zeros(4,4,4,4);
for nn=1:4
	for m=1:4
PI(nn,m,:,:) = 0.5*U*(kron(reshape(U(nn,:),2,2),conj(reshape(U(m,:),2,2))))*U';
	end
end
PI=reshape(PI,16,16);

%% scaling MM with eigenvalue weight -jaclyn

SumM = zeros(((sqrt(nPixels)+2)*4),((sqrt(nPixels)+2)*4)); %create an empty array to add the eigenmatrices to
%size is the size of image create in the loop below (with borders)

for z=1:4
    M=zeros(16,nPixels);
    for nn=1:nPixels
        temp=squeeze(C(n,nn,:,z));
        temp=temp*temp';
        M(:,nn)=reshape(temp,1,16)*PI';
    end
    M=reshape(M,16,sqrt(nPixels),sqrt(nPixels)).*permute((repmat(squeeze(S(:,:,z)),[1 1 16])),[3 1 2]);
    %M = reshape(M,[4 4 sqrt(nPixels) sqrt(nPixels)]);
    %inverse fesnel mirror MM, uncoated aluminum at 45 degrees, 532 nm 
    %F=[1.089729264540746 0.03105722428292043 0 0;0.03105722428292043 1.0897292645407457 0 0;0 0 -1.064747642529994 -0.22990819155300934; 0 0 0.22990819155300934 -1.064747642529994];
    %M = pagemtimes(M,F);
    %M = reshape(M,[16 sqrt(nPixels) sqrt(nPixels)]);
    figure;
    tempImg=[];
    for nn=1:4
        tempRow=[];
        for zz=1:4
            tempSubImg=squeeze(M(zz+(nn-1)*4,:,:));
            tempSubImg=[min(M(:))*zeros(sqrt(nPixels),1)';tempSubImg;min(M(:))*zeros(sqrt(nPixels),1)'];%Put borders on image
            tempSubImg=[min(M(:))*zeros(sqrt(nPixels)+2,1),tempSubImg,min(M(:))*zeros(sqrt(nPixels)+2,1)];
            tempRow=[tempRow;tempSubImg];
        end
        tempImg=[tempImg,tempRow];
	end
    t=prctile((abs(M(:))),90);
    imagesc(tempImg,[-t,t]);axis square; axis off; colormap(cmap);colorbar;
    SumM = SumM + tempImg; %adding eigen matrix to SumM each loop iteration 
   
cmd=sprintf('print -dpdf %s',['''',w(1:end-3),sprintf('_%dMM%d_EigenSCALEDLambda%d.pdf',z,index,LambdaList(n)),'''']);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gcf,'renderer','painters');drawnow;
eval(cmd); close all;

end

%plotting the sum of the eigenmatrices to compare to original MM 

t=prctile((abs(I(:))),90);
imagesc(SumM,[-t,t]);axis square; axis off; colormap(cmap);colorbar;
cmd=sprintf('print -dpdf %s',['''',w(1:end-3),sprintf('_%dMM%d_EigenSUMLambda%d.pdf',z,index,LambdaList(n)),'''']);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(gcf,'renderer','painters');drawnow;
eval(cmd); close all;


return
