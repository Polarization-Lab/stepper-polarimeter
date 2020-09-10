function TransmissionMMPlot(LambdaList,mmVecs)

mmNum = [ "00" "01" "02" "03" '10' '11' '12' '13' '20' '21' '22' '23' '30' '31' '32' '33'];
%M00 transmission
%% Average Transmission
mkdir(['Y:\Measurements\Dichroic_Analysis\MM_img\dichroic-45-small-ROI-' date]);
figure(1)
meanTransmission = mean(squeeze(mmVecs(:,1,:,:)),[2 3]);
plot(LambdaList,meanTransmission,'k-*')
title('Avg Transmission over \lambda')
xlabel('\lambda (nm)')
ylabel('Transmission')
saveas(gcf,fullfile(['Y:\Measurements\Dichroic_Analysis\MM_img\dichroic-45-small-ROI-' date],'Avg_Transmission_Lambda.png'));

% MM transmission curves 
figure(2)
for p = 1:16
    subplot(4,4,p)
    plot(LambdaList,mean(squeeze(mmVecs(:,p,:,:)),[2 3]),'b*-');
    xlabel('\lambda (nm)')
    title(strcat('M',mmNum(p)))
    ylim([-1.5 1.5])
end
sgtitle('MM Transmission curves')

