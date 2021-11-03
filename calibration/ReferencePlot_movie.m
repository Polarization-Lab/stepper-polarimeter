function intenvideo = ReferencePlot_movie(Directory,file,lambdas,num_meas)

intenvideo = VideoWriter(strcat(Directory,file,'-Referenceplotmovie.mp4'),'MPEG-4');
intenvideo.FrameRate = 32/15;
open(intenvideo);
set(0,'DefaultFigureVisible', 'off');
file = strcat(file,'.h5');


for jj = 1:length(lambdas)
for ii = 1:num_meas 
    group_name   = strcat('/images/wave',num2str(lambdas(jj)),'/meas',num2str(ii));
    refdetdata = h5read(strcat(Directory,file),strcat(group_name,'/refdata'));
    
    if refdetdata ==10.904
        plot(ii,refdetdata,'.r','MarkerSize',12)
    else
        plot(ii,refdetdata,'.k','MarkerSize',12)
    end
    title(strcat('\lambda =',num2str(lambdas(jj)),'[nm]'))
    hold on
    grid on
    xlabel('Measurement Number')
    ylabel('Ref Detector Voltage')
    axis([0 64 -11 11])
    
    frame = getframe(gcf);
    writeVideo(intenvideo,frame);
end 
hold off
end
set(0,'DefaultFigureVisible', 'on');


close(intenvideo);
