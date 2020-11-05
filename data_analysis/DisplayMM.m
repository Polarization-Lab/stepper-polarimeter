function DisplayMM(mmVecs,nLambda,LambdaList,filename)

mkdir(['D:\Measurements\Dichroic_Analysis\MM_img\Zero\' date]);
fn = filename;
% Display image MMs
tic
close all;

%mx=max(mmVecs(:));
mmNum = [ "00" "01" "02" "03" '10' '11' '12' '13' '20' '21' '22' '23' '30' '31' '32' '33'];

for n = 1:nLambda
    lims = [-abs(max(squeeze(mmVecs(n,:)),[],'all')) abs(max(squeeze(mmVecs(n,:)),[],'all'))];
    lims = [-1 1];

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
    F(n) = getframe(gcf);
    saveas(gcf,fullfile(['D:\Measurements\Dichroic_Analysis\MM_img\Zero\' date],['Lambda ' num2str(LambdaList(n)) '_' date '.png']));
    close all
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
