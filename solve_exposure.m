function [ms,exposures] = solve_exposure(wavelength,COMmono,vid)
%find_exposure this function helps find the PSA/PSG orientation with the
%assumes PSA/PSG are in optimal position

nominalexposure = 0.05;

%set monochromator wavelength

changeWavelength(COMmono,wavelength)

%configure camera settings
src = getselectedsource(vid);
triggerconfig(vid, 'immediate');
vid.FramesPerTrigger = 1;
src.ExposureTimeControl = 'normal';
src.ExposureTime = nominalexposure;

preview(vid)
start(vid)
im=getdata(vid);
m = mean(mean(im));


exposures = .05:.02:.3;
m = zeros(size(exposures));
dim = size(exposures);
dim = dim(2);
for ii = 1:dim
    
    src.ExposureTime = exposures(ii);
    start(vid)
    pause(.4)
    im = getdata(vid);
    m =mean(max(im));
    stop(vid)
    
    ms(ii) = m;

end

plot(exposures,ms,'-o'), xlabel('exposure [s]'), ylabel('average count')
title(strcat(num2str(wavelength),' nm'))
grid on

end

