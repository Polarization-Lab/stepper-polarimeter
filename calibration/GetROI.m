function [rect] = GetROI(fn,wavelength)
    group_name   = strcat('/images/wave',num2str(wavelength),'/meas',num2str(1));
	norm = h5read(fn,strcat(group_name,'/imagedata'));
	[~,rect] = imcrop(norm/max(norm(:)));close all;
	rect=round(rect);
end
