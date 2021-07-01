function video = makeamovie(filepath,wavelength,nSteps)
%filpath - to h5 file where measurement data is held 
%wavelength 
%nSteps for stepper should be 64 

filepath = '/Users/jaclynjohn/Documents/MATLAB/aircal.h5';
wavelength = 460;
nSteps = 64;

video = VideoWriter('aircal')
open(video)

for ii = 1:nSteps 
    group_name   = strcat('/images/wave',num2str(wavelength),'/meas',num2str(ii));
    images = h5read(filepath,strcat(group_name,'/imagedata'));
    pause(0.5);
    imagesc(images)
    frame = getframe(gcf)
    writeVideo(video,frame)
end


close(video)


%to do: slow down frame rate

