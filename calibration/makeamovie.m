function video = makeamovie(filepath,wavelength,nSteps)
%filpath - to h5 file where measurement data is held 
%wavelength 
%nSteps for stepper should be 64 

for ii = 1:nSteps 
    group_name   = strcat('/images/wave',num2str(wavelength),'/meas',num2str(ii));
    images = h5read(filepath,strcat(group_name,'/imagedata'));
end


video = VideoWriter('aircal','Motion JPEG AVI')

open(video)


for k = 1:nSteps  
   imagesc(images)
   frame = getframe(gcf);
   writeVideo(video,frame);
end

%not done working on this
%will be used to look at all images from a measurement in a movie format
