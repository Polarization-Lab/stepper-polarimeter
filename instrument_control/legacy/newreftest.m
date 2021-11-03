% waves = 350:10:800;
% exposure = .05;
% 
% src = getselectedsource(vid);
% triggerconfig(vid, 'immediate');
% vid.FramesPerTrigger = 5;
% src.ExposureTimeControl = 'normal';
% src.ExposureTime = exposure;
% 
% camresp = zeros(1,46);
% ref     = zeros(1,46);
% 
% for i = 1:46
%     changeWavelength(COMmono,waves(i))
%     pause(1)
%     
%     [image,r] = take_snapshot(vid, exposure , 5);
%     ref(i) =r;
%     camresp(i) = mean(image,'all');
% end

yyaxis left
plot(waves,camresp);
yyaxis right
plot(waves,ref);
