function meandata = ReadCalImages(filepath,wavelength,ref,nSteps)

%Wavelength = scalar wavelength value
%filepath = string to .h5 file
%ref = 1x64 of mean reference fluxes normalized to first reference flux
%nSteps = number of measurements taken (64 usually)

fn = extractAfter(filepath,"Dichroic_Analysis\");
wv = wavelength;

for ii = 1:nSteps
    group_name   = strcat('/images/wave',num2str(wv),'/meas',num2str(ii));
    images = h5read(fn,strcat(group_name,'/imagedata'))/ref(ii);
%     mean_img = (images(1:100,1401:1500,1,1) + images(1:100,1401:1500,1,2) + images(1:100,1401:1500,1,3))/3;
    mean_img = (images(701:705,701:705,1,1) + images(701:705,701:705,1,2) + images(701:705,701:705,1,3))/3;
    meandata(ii,:,:) = mean_img;
end

