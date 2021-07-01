function meandata = ReadCalImages(fn,wavelength,ref,nSteps)

%Wavelength = scalar wavelength value
%filepath = string to .h5 file
%ref = 1x64 of mean reference fluxes normalized to first reference flux
%nSteps = number of measurements taken (64 usually)

fn = ['D:/Measurements/Air_Calibrations/BVO-UV_PSA-LP/BVO-UV_PSA-LPjune302021-30-Jun-2021.h5'];

%fn = extractAfter(filepath,"");

%wavelength = 460;

%nSteps = 64

for ii = 1:nSteps
    group_name   = strcat('/images/wave',num2str(wavelength),'/meas',num2str(ii));
    images = h5read(fn,strcat(group_name,'/imagedata'))/mean(ref);
    meandata(ii,:,:) = images;
end

