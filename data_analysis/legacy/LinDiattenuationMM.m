function LinDiattenuationMM(mmVecs,LambdaList,nLambda,Label)

for n = 1:nLambda    
    M_00(n,:,:)= squeeze(mmVecs(n,1,:,:));
    M_01(n,:,:)= squeeze(mmVecs(n,2,:,:));
    M_02(n,:,:)= squeeze(mmVecs(n,3,:,:));
    M_03(n,:,:)= squeeze(mmVecs(n,4,:,:));
    M_10(n,:,:)= squeeze(mmVecs(n,5,:,:));
    M_11(n,:,:)= squeeze(mmVecs(n,6,:,:));
    M_12(n,:,:)= squeeze(mmVecs(n,7,:,:));
    M_13(n,:,:)= squeeze(mmVecs(n,8,:,:));
    M_20(n,:,:)= squeeze(mmVecs(n,9,:,:));
    M_21(n,:,:)= squeeze(mmVecs(n,10,:,:));
    M_22(n,:,:)= squeeze(mmVecs(n,11,:,:));
    M_23(n,:,:)= squeeze(mmVecs(n,12,:,:));
    M_30(n,:,:)= squeeze(mmVecs(n,13,:,:)); 
    M_31(n,:,:)= squeeze(mmVecs(n,14,:,:));
    M_32(n,:,:)= squeeze(mmVecs(n,15,:,:));
    M_33(n,:,:)= squeeze(mmVecs(n,16,:,:));
end

M_00(M_00 == 0) = NaN;

mkdir(['d:\Measurements\Dichroic_Analysis\Linear_Diattenuation\dichroic-45-small-ROI-' date]);

for n = 1:nLambda
    %Equation 6.49 page 179 from Russell's book
    LinDia(n,:,:) = sqrt(M_01(n,:,:).^2 + M_02(n,:,:).^2)./M_00(n,:,:);
    %Plot LD map and histogram
    figure(2)
    subplot(1,2,1)
    imshow(squeeze(LinDia(n,:,:)),[0 1],'colormap',parula);colorbar;
    title(['2D Linear Diatten. per pixel ' Label num2str(LambdaList(n))]);
    subplot(1,2,2)
    histogram(squeeze(LinDia(n,:,:)))
    title(['Linear Diattenuation Histogram (' num2str(LambdaList(n)) ')']);
    saveas(gcf,fullfile(['d:\Measurements\Dichroic_Analysis\Linear_Diattenuation\dichroic-45-small-ROI-' date],['LD_Lambda_' num2str(LambdaList(n)) '_' date '.png']));
    close all;
end

avgLD = mean(LinDia,[2 3]);
plot(LambdaList,avgLD,'r*-');
title('Avg Linear Diattenuation over \lambda')
xlabel('\lambda (nm)')
ylabel('Linear Diattenuation')
saveas(gcf,fullfile(['d:\Measurements\Dichroic_Analysis\Linear_Diattenuation\dichroic-45-small-ROI-' date],['Avg_LinDia_Lambda.png']));
