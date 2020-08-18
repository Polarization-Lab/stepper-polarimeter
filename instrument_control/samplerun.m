%create HDF5 dataset

%writefilepath
fp = 'c:/Stepper Polarimeter/Stepper Software Rebuild/TestData/Dichroic/';
date = date();
starttime = datestr(now);

usernotes = 'taken by khart ; dichroic at 0 degrees, new exposures';

fn = 'dichroic0';
name = strcat(fp,fn,'-',date,'.h5');

num_meas = 64;

wavelength=450;
exposure = 0.35;
framesPerTrigger = 3;
homeMotor(xps)
wavelengthSweep(name,wavelength,exposure ,vid,num_meas, COMmono,COMdmm, xps, framesPerTrigger)

wavelength=550;
exposure = 0.15;
framesPerTrigger = 3;
homeMotor(xps)
wavelengthSweep(name,wavelength,exposure ,vid,num_meas, COMmono,COMdmm, xps, framesPerTrigger)

wavelength=650;
exposure = 0.17;
framesPerTrigger = 3;
homeMotor(xps)
wavelengthSweep(name,wavelength,exposure ,vid,num_meas, COMmono,COMdmm, xps, framesPerTrigger)

wavelength=750;
exposure = 0.25;
framesPerTrigger = 3;
homeMotor(xps)
wavelengthSweep(name,wavelength,exposure ,vid,num_meas, COMmono,COMdmm, xps, framesPerTrigger)

endtime = datestr(now);

[PSA, PSG] = generate_PSAG_angles(num_meas);

%write attibutes to directory
 h5writeatt(name,'/images/','start_time', starttime);
 h5writeatt(name,'/images/','end_time', endtime);
 h5writeatt(name,'/images/','user_notes', usernotes);
 h5writeatt(name,'/images/','PSG_positions', PSG); 
 h5writeatt(name,'/images/','PSA_positions', PSA); 

 % close ports 
fclose('all')
clear all
close all
clc
instrreset
