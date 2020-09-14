function DisplayNormMM(mmVecs,LambdaList,nLambda)
%One wavelength only. Needs work
for ii = 1:16
    NormMMVecs(ii,:,:) = squeeze(mmVecs(ii,:,:))./squeeze(mmVecs(1,:,:));
end
% Display Normalized MM
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
