function Polarizance(mmVecs, LambdaList, nLambda, Label)

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

mkdir(['d:\Measurements\Dichroic_Analysis\Polarizance\dichroic-45-small-ROI-' date]);
for n = 1:nLambda
    %Equation 6.51 page 180 from Russell's book
    % Pol = sqrt(mmVecs(5,:,:).^2 + mmVecs(9,:,:).^2 + mmVecs(13,:,:).^2)./mmVecs(1,:,:);
    Pol(n,:,:) = sqrt(M_10(n,:,:).^2 + M_20(n,:,:).^2 + M_30(n,:,:).^2)./M_00(n,:,:);
    %Plot Polarizance map and histogram
    subplot(1,2,1)
    imshow(squeeze(Pol(n,:,:)),[0 1],'colormap',parula);colorbar;
    title(['2D Polarizance per pixel ' Label num2str(LambdaList(n))]);

    subplot(1,2,2)
    histogram(squeeze(Pol(n,:,:)))
%     imhist(squeeze(Pol))
    title(['Polarizance Histogram (' num2str(LambdaList(n)) ')']);
    saveas(gcf,fullfile(['d:\Measurements\Dichroic_Analysis\Polarizance\dichroic-45-small-ROI-' date],['Pol_Lambda_' num2str(LambdaList(n)) '_' date '.png']));
    close all;
end
%%
avgPol = mean(Pol,[2 3]);
plot(LambdaList,avgPol,'g*-');
title('Avg Polarizance over \lambda')
xlabel('\lambda (nm)')
ylabel('Polarizance')
saveas(gcf,fullfile(['d:\Measurements\Dichroic_Analysis\Polarizance\dichroic-45-small-ROI-' date],['Avg_Pol_Lambda.png']));
