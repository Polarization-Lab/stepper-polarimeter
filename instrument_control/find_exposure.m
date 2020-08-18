function [PSAmax,PSGmax] = find_exposure(wavelength, nominalexposure,meas_num,COMmono,vid,xps)
%find_exposure this function helps find the PSA/PSG orientation with the
%maximum exposure for a given wavelength and then ups the exposure to fill
%the dynamic range of the camera

%set monochromator wavelength

changeWavelength(COMmono,wavelength)

%configure camera settings
src = getselectedsource(vid);
triggerconfig(vid, 'immediate');
vid.FramesPerTrigger = 1;
src.ExposureTimeControl = 'normal';
src.ExposureTime = nominalexposure;

preview(vid)

samp =20;

%load PSA/PSG Angles
[PSA, PSG] = generate_PSAG_angles(meas_num);

%find maximum and minimum transmission
averages= zeros(1,samp);

i = 1;
while i <samp
    disp(strcat('on step',num2str(i)))
    movePSG(xps,PSG(i));
    movePSA(xps,PSA(i));
    start(vid)
    pause(.2)
    im = getdata(vid);
    averages(i)=mean(max(im));
    i = i +1;
    stop(vid)
end


%plot results
x = 1:1:samp;
plot(x, averages,'-o'), xlabel('measurement #'), ylabel('ADU')
grid on

%find maximum position and move PSG/PSA there
m = max(averages);
i = find(averages==m);

movePSG(xps,PSG(i));
movePSA(xps,PSA(i));

PSAmax = PSA(i);
PSGmax = PSG(i);

end

