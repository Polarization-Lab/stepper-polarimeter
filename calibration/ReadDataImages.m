function [imgdata ROI_imgdata] = ReadDataImages(filepath, wavelength,ref,nSteps)

%fn = extractAfter(filepath,"\New_Diffuser\");

wv = wavelength;

for ii = 1:nSteps
    group_name   = strcat('/images/wave',num2str(wv),'/meas',num2str(ii));
    images       = h5read(fn,strcat(group_name,'/imagedata'))/mean(ref);
    ROI_imgdata(ii,:,:) = images(746:755,746:755); %ROI size image. Need to make this dynamic
    imgdata(ii,:,:) = images;                      %Full images
end

