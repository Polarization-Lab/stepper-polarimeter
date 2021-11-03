function [] = measurement(fp,fn,usernotes,num_meas,ROI,framesPerTrigger,wavelengths,exposures)
%MEASUREMENT this function execute a full measurement sequence
%for an input list of wavelengths and exposures

fclose('all');
close all
clc
instrreset

initialize
vid.ROIPosition = ROI;

%create HDF5 dataset
dat = date();
starttime = datestr(now);

name = strcat(fp,fn,'-',dat,'.h5');

if size(wavelengths) ~= size(exposures)
    disp('ERROR; exposure and wavelength vectors must be equal')
end

[~,len] = size(wavelengths);

%step through wavelengths
for i = 1:len

    wavelength=wavelengths(i);
    exposure = exposures(i);
    homeMotor(xps)
    wavelengthSweep(name,wavelength,exposure ,vid,num_meas, COMmono,COMdmm, xps, framesPerTrigger)
end

endtime = datestr(now);

[PSA, PSG] = generate_PSAG_angles(num_meas);

%write attibutes to directory
 h5writeatt(name,'/images/','start_time', starttime);
 h5writeatt(name,'/images/','end_time', endtime);
 h5writeatt(name,'/images/','user_notes', usernotes);
 h5writeatt(name,'/images/','PSG_positions', PSG); 
 h5writeatt(name,'/images/','PSA_positions', PSA); 

% close ports 
fclose('all');
close all
clc
instrreset

end

