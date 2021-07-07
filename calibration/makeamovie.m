function video = makeamovie(Directory,file,lambdas)
%Directoy - to h5 file where measurement data is held
%file  - name of .h5 file to be read, input without extension
%lambdas - single value or vector of wavelengths to show images for  

%Directory =  'D:\Measurements\Waveguides\option1\'; %change 
%file = 'SampleMeas-460_530_625-05-Jul-2021.h5'
%lambdas = 530;


video = VideoWriter(strcat(Directory,file,'-Rawimagemoive'))
video.FrameRate = 32/15
open(video)
file = strcat(file,'.h5')

for jj = 1:length(lambdas)
for ii = 1:64 
    group_name   = strcat('/images/wave',num2str(lambdas(jj)),'/meas',num2str(ii));
    images = h5read(strcat(Directory,file),strcat(group_name,'/imagedata'));
    subplot(1,2,1)
    imagesc(images)
    colormap(gray)
    colorbar
    subplot(1,2,2)
    plot(ii,max(smooth(images)),'.k','MarkerSize',12)
    title(strcat('wavelength =',num2str(lambdas(jj)),'[nm]'))
    hold on
    grid on
    xlabel('measurment number')
    ylabel('max camera counts')
    axis([0 64 0 inf])
    frame = getframe(gcf)
    writeVideo(video,frame)
end
hold off
end


close(video)



