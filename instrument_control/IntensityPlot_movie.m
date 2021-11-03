function intenvideo = IntensityPlot_movie(Directory,file,lambdas,num_meas,Exps)
%Directoy - to h5 file where measurement data is held
%file  - name of .h5 file to be read, input without extension
%lambdas - single value or vector of wavelengths to show images for  

%Directory =  'D:\Measurements\Waveguides\option1\'; %change 
%file = 'SampleMeas-460_530_625-05-Jul-2021.h5'
%lambdas = 530;


intenvideo = VideoWriter(strcat(Directory,file,'-Intensityplotmovie.mp4'),'MPEG-4');
intenvideo.FrameRate = 32/15;
open(intenvideo);
set(0,'DefaultFigureVisible', 'off');
file = strcat(file,'.h5');

for jj = 1:length(lambdas)
for ii = 1:num_meas 
    group_name   = strcat('/images/wave',num2str(lambdas(jj)),'/meas',num2str(ii));
    images = h5read(strcat(Directory,file),strcat(group_name,'/imagedata'));
   
    PSG_steps = (ii-1)*2*pi/num_meas;
   
    if prctile(images(:),99)==2^16-1
        plot(PSG_steps,prctile(images(:),99),'.r','MarkerSize',12)
    else
        plot(PSG_steps,prctile(images(:),99),'.k','MarkerSize',12)
    end
    title(strcat('\lambda =',num2str(lambdas(jj)),'[nm]',' Exposure=',num2str(Exps(jj)),'[s]'))
    hold on
    grid on
    xlabel('PSG Rotation Step [Radians]')
    ylabel('max camera counts')
    axis([0 (num_meas-1)*2*pi/num_meas 0 70000])

    
    frame = getframe(gcf);
    writeVideo(intenvideo,frame);
end
hold off
end

set(0,'DefaultFigureVisible', 'on');


close(intenvideo);
