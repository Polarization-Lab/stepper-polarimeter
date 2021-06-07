function meandata = ReadCalImages(filepath,wavelength,ref,nSteps)

%Wavelength = scalar wavelength value
%filepath = string to .h5 file
%ref = 1x64 of mean reference fluxes normalized to first reference flux
%nSteps = number of measurements taken (64 usually)

fn = extractAfter(filepath,"Air_Calibrations\New_Diffuser\");
wv = wavelength;

for ii = 1:nSteps
    group_name   = strcat('/images/wave',num2str(wv),'/meas',num2str(ii));
    images = h5read(fn,strcat(group_name,'/imagedata'))/mean(ref);
    meandata(ii,:,:) = images;
end

