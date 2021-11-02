for ii = 1:64
    figure(1);
    lims = [0 61000];
    imagesc(squeeze(imgdata(10,ii,:,:)),lims); colorbar;
    title(['Raw Sample 540 nm, Measurement = ' num2str(ii)])
    
    F(ii) = getframe(gcf);
    close all;
end

writerObj = VideoWriter('SampleImg','MPEG-4');
writerObj.FrameRate = 1;

open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj,frame);
end
close(writerObj)
