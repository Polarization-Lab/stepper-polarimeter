function mmNorm = DisplayNormMM(mmVecs,LambdaList,filename)
%One wavelength only. Needs work
fn = filename;
for n = 1:length(LambdaList)
    for ii = 1:16
        NormMMVecs(ii,:,:) = squeeze(mmVecs(n,ii,:,:))./squeeze(mmVecs(n,1,:,:));
    end

    mmNorm(n,:,:,:) = NormMMVecs;

    mkdir(['D:\Measurements\Dichroic_Analysis\Norm_MM_img\Zero\' date]);
    % Display Normalized MM
    tic
    
    close all;
    %mx=max(mmVecs(:));
    
    figure('units','normalized','outerposition',[0 0 1 1])
    lims = [-abs(max(NormMMVecs,[],'all')) abs(max(NormMMVecs,[],'all'))];

    mmNum = [ "00" "01" "02" "03" '10' '11' '12' '13' '20' '21' '22' '23' '30' '31' '32' '33'];
    for p=1:16
        subplot(4,4,p)
        imshow(squeeze(NormMMVecs(p,:,:)),lims,'Colormap',GWP);
        text(5,5, [num2str(round(avgMMcomp(n,p)*100,2,'significant')) '%']);
    %     imshow(TotalMM(p),lims,'Colormap',GWP);
        title(strcat('M',mmNum(p)))
    end

    t = (subplot(4,4,16).Position);
    colorbar('position', [t(1)+t(3) t(2) t(3)/3 t(4)*4.7] );
    sgtitle(['Dichroic 45' char(176) ' ' num2str(LambdaList(n)) 'nm Norm' ])
    F(n) = getframe(gcf);
    saveas(gcf,fullfile(['D:\Measurements\Dichroic_Analysis\Norm_MM_img\Zero\' date],['Lambda ' num2str(LambdaList(n)) '_' date '.png']));
end

writerObj = VideoWriter(fn,'MPEG-4');
writerObj.FrameRate = 1;

open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj,frame);
end
close(writerObj)
toc
