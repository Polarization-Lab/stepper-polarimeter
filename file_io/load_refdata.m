function [ref_mean,ref_std,refVecs] = load_refdata(fn,wavelength,meas_num)
%LOAD_IMAGEDATA this function loads image data from the hdf5 file for a
%single wavelength. Will average over apporipriate number of frames and
%display associated metadata
%   fn is the filename
%   wavelength is the wavelength in nm
%   meas_num is the number of measurements in the file

%allocate space 
ref_mean = zeros(meas_num,1);
ref_std = zeros(meas_num,1);

for i = 1:meas_num
    group_name   = strcat('/images/wave',num2str(wavelength),'/meas',num2str(i));
    ref          = h5read(fn,strcat(group_name,'/refdata'));
    
    ref_mean(i)         = mean(ref);
    ref_std(i)          = std(ref);
end
refVecs = ref_mean/ref_mean(1);

end

