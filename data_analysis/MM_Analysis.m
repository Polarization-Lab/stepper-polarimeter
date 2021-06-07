function MM_Analysis(samplefilepathname,xBegROI,xEndROI,yBegROI,yEndROI,nLambda,LambdaList,nSteps,Amp,W,SpecOrAmp)
%% Generate Mueller Matrices This is all analysis, different file!
for ii = 1:length(LambdaList) %for loading all lambda(s)
%Get h5 reference data
    Lambda = LambdaList(ii);
    [~,~,refVecsImgs(ii,:)] = load_refdata(samplefilepathname,Lambda,nSteps);
    [~, ROIimg(ii,:,:,:)] = ReadDataImages(samplefilepathname, Lambda,refVecsImgs(ii,:),nSteps);
%     imgdatanoref(ii,:,:,:) = ReadDataImagesNoRef(testFilePath, Lambda,nSteps);
% save('dichroic_575_45','imgdata','ROIimg');imgdata(ii,:,:,:)
end

%Message for finished operation and play Hallelujah
f = msgbox('Operation Completed');
load handel
sound(y,Fs)

%% Change ROI
%Grab ROI corresponding to Spectral or Amplitude analysis
%if SpecORrmp = 1, Spectral
%if SpecOrAmp = 0, Amplitude
if (SpecOrAmp == 1)
    
elseif (SpecOrAmp = 0)
end
%%
ampcaldata = reshape(Amp(:),31,10,10);

%%
imgdata = ROIimg;
tic
% nLambda=length(Lambda);
M = zeros(nLambda,16,10*10);
for n = 1:nLambda
    Winv = pinv(squeeze(W(n,:,:)));
    M(n,:,:) = Winv*(reshape(squeeze(ROIimg(n,:,:,:)),nSteps,ROIx*ROIy)./repmat(ampcaldata(n,:),[nSteps 1]));
end

toc
%% .h5 ROI image1
mmVecs = CreateMuellerMatrix(Amp,W,nLambda,10,ROIimg,nSteps);
%%
fn = strcat("Dichroic0_mmVecs_", date);
save(fn,"mmVecs");
%% Normalize MM
for ii = 1:16
    NormMMVecs(ii,:,:) = squeeze(mmVecs(ii,:,:))./squeeze(mmVecs(1,:,:));
end

%% Create averaged scalar MM (Table form)
AvgMMTable(mmVecs,670,LambdaList)

%% Display image MMs

filename = 'Dichroic45_12152020';
DisplayMM(mmVecs,nLambda,LambdaList,filename)


%% Display Normalized MM
DisplayNormMM(mmVecs,LambdaList,nLambda)

%% MM transmission curves 

TransmissionMMPlot(LambdaList,mmVecs)

%% M_xy setup for easier calculations
Label = ['Dichroic 0' char(176)];

for n = 1:nLambda
    Label = ['Dichroic 0' char(176)];

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

%% Diattenuation
DiattenuationMM(mmVecs,LambdaList,nLambda,Label)
%% Linear Diattenuation
LinDiattenuationMM(mmVecs,LambdaList,nLambda,Label)

%% Polarizance
Polarizance(mmVecs, LambdaList, nLambda, Label)
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
DepolIndex(mmVecs, LambdaList, nLambda, Label)
end