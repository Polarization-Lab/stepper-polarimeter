function [] = wavelengthSweep(fn,wavelength, exposure , vid, num_meas, COMmono,COMdmm, xps, framesPerTrigger)
%WAVELENGTHSWEEP Summary of this function goes here
%   fn  filename string path to .h5 fil.e
%   WAVELENGTH double wavelength of exposure [nm]
%   EXPOSURE double exposure time for this wavelength [sec]
%   PSG array with PSG positions
%   PSA array with PSA positions
%   VID video input object 
%   COMmono serial variable to monochromator 
%   COMdmm serial variable to dmm
%   xps XPS to motor 
%   framesPerTrigger number of frames per trigger

%load PSA/PSG Angles
[PSA, PSG] = generate_PSAG_angles(num_meas);
 
%set monochromator wavelength
changeWavelength(COMmono,wavelength)

%configure camera settings
src = getselectedsource(vid);
triggerconfig(vid, 'immediate');
vid.FramesPerTrigger = framesPerTrigger;
src.ExposureTimeControl = 'normal';
src.ExposureTime = exposure;

%given settings, find correct image dimension for saving
start(vid)
im = getdata(vid);
imdim = size(im);


%create wavelength fp
wavename = strcat('/images/wave',num2str(wavelength),'/'); %name for this wavelength group

%take darkfield image
NI_shutter(0)%close shutter
pause(1)

start(vid)
prepare_DMM(COMdmm,exposure)
pause(.2)
dark_im = getdata(vid);
trigger_DMM(COMdmm)
stop(vid)
NI_shutter(1) %open_shutter

pause(exposure)
dmmdata = read_DMM(COMdmm,exposure); %read DMM data from buffer


%write dark image to h5
imname = strcat('/images/wave',num2str(wavelength),'/darkdata');
dmmname= strcat('/images/wave',num2str(wavelength),'/darkref');

h5create(fn,imname,size(dark_im))%create image dataset
h5create(fn,dmmname,size(dmmdata))%create image dataset

h5write(fn,imname,dark_im);
h5write(fn,dmmname, dmmdata);



%set up measurement for number of PSG/PSA measurements
i=1;

%loop through positions
while i < num_meas +1
    g = PSG(i);
    a = PSA(i);
    
    movePSG(xps,g) %MOVE psg
    movePSA(xps,a) %MOVE PSA
    
    %define image group name
    meas_name = strcat(wavename,'meas',num2str(i),'/'); 
        
    start(vid)
    prepare_DMM(COMdmm,exposure);
    
    %get image data
    imdata = getdata(vid);
    trigger_DMM(COMdmm);
        
    %take images
    pause(exposure)
    dmmdata =read_DMM(COMdmm,exposure);
        
        
       
        
    %write to h5
    imname = strcat(meas_name,'imagedata');
    dmmname= strcat(meas_name,'refdata');

    h5create(fn,imname,size(imdata))%create image dataset
    h5create(fn,dmmname,size(dmmdata))%create image dataset

    h5write(fn,imname,imdata);
    h5write(fn,dmmname, dmmdata);

    clear imdata
    clear dmmdata

    stop(vid)
    i=i+1;
end

%write attibutes to directory
 h5writeatt(fn,wavename,'exposure_time', exposure);
 h5writeatt(fn,wavename,'frames_per_trigger', framesPerTrigger);
 h5writeatt(fn,wavename,'image_dimension', imdim);
end

