function DepolIndex(mmVecs, LambdaList, nLambda, Label)

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

close all;
mkdir(['d:\Measurements\Dichroic_Analysis\DepolIndex\dichroic-45-small-ROI-' date]);
%Sum of squares variable
M_ij = zeros(nLambda,10,10);

for n = 1:nLambda
%create sum of squares and calculate DI. Equation 6.79 page 195 from Russell's book
    for ii = 1:16
        M_ij(n,:,:) = squeeze(mmVecs(n,ii,:,:)).^2 + squeeze(M_ij(n,:,:));
    end
end

for n= 1:nLambda
    % DepolIndex = sqrt(M_ij-squeeze(mmVecs(1,:,:)).^2)./(sqrt(3)*squeeze(mmVecs(1,:,:)));
    for ii = 1:16
        DepolIndex(n,:,:) = sqrt(squeeze(M_ij(n,:,:))-squeeze(M_00(n,:,:)).^2) ./ (sqrt(3)*squeeze(M_00(n,:,:)));
    end
    
    figure(1)
    %Plot DI map
    subplot(1,2,1)
    imshow(squeeze(DepolIndex(n,:,:)),[0 1],'colormap',parula);colorbar;
    title(['2D Depolarization Index per pixel ' Label num2str(LambdaList(n))]);

    subplot(1,2,2)
    histogram(squeeze(DepolIndex(n,:,:)))
    title(['Depolarization Index Histogram (' num2str(LambdaList(n)) ')']);
    saveas(gcf,fullfile(['d:\Measurements\Dichroic_Analysis\DepolIndex\dichroic-45-small-ROI-' date],['DI_Lambda_' num2str(LambdaList(n)) '_' date '.png']));
    close all;
end

avgDI = mean(DepolIndex,[2 3]);
plot(LambdaList,avgDI,'m*-');
title('Depolarization Index over \lambda')
xlabel('\lambda (nm)')
ylabel('Depol. Index')
saveas(gcf,fullfile(['d:\Measurements\Dichroic_Analysis\DepolIndex\dichroic-45-small-ROI-' date],'Avg_DI_Lambda.png'));

