function refVecs = NormalizedReferenceData(filepath, wavelength, nSteps)

%Kira Hart
%this script reads in and plots the reference detector data


fn      = extractAfter(filepath,"Air_Calibrations\");
wv      = wavelength;

m       = zeros(1,64);
s       = zeros(1,64);

for i = 1:64
    group_name   = strcat('/images/wave',num2str(wv),'/meas',num2str(i));
   
    ref          = h5read(fn,strcat(group_name,'/refdata'));
    
    m(i)         = mean(ref);
    s(i)         = std(ref);
end

refVecs = m/m(1);

%make plot
% meas = 1:64;
% errorbar(meas,m,s,'ko')
% hold on
% xlabel('measurement number','FontSize',14)
% ylabel('Voltage','FontSize',14)
% title(['Average DMM Voltage at ',num2str(wv),'nm'],'FontSize',14)
