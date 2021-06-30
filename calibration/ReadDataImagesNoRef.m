function imgdata = ReadDataImagesNoRef(filepath, wavelength, nSteps)

%fn = extractAfter(filepath,"Air_Calibrations\New_Diffuser\");

wv = wavelength;

for ii = 1:nSteps
    group_name   = strcat('/images/wave',num2str(wv),'/meas',num2str(ii));
    images       = h5read(filepath,strcat(group_name,'/imagedata'));%/mean(ref);
%     test_img     = (images(701:705,701:705,1,1) + images(701:705,701:705,1,2) + images(701:705,701:705,1,3))/3;
%     mean_img     = (images(:,:,1,1) + images(:,:,1,2) + images(:,:,1,3))/3;
   
%     ROI_imgdata(ii,:,:) = test_img;
    imgdata(ii,:,:) = images;
end

