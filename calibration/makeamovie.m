function video = makeamovie(filepath,wavelength)
%filpath - to h5 file where measurement data is held 
%wavelength 
%nSteps for stepper should be 64 

filepath = '/Users/jaclynjohn/Documents/MATLAB/aircal.h5'; %change 
wavelength = 460;

video = VideoWriter('D:\Measurements\Waveguides\option1\530') %path where avi file will be saved
video.FrameRate = .1 
open(video)

for ii = 1:64
    group_name   = strcat('/images/wave',num2str(wavelength),'/meas',num2str(ii));
    images = h5read(filepath,strcat(group_name,'/imagedata'));
    %pause(0.5); %slow down frames
    imagesc(images)
    frame = getframe(gcf)
    writeVideo(video,frame)
end


close(video)



