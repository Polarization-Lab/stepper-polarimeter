function mmVecs = CreateMuellerMatrix(Amp,W,nLambda,ROI,ROIimg,nSteps)
AmpArray = reshape(Amp(:),30,5,5);
M = zeros(nLambda,16,ROI*ROI);

for n = 1:nLambda
    Winv = pinv(squeeze(W(n,:,:)));
    M(n,:,:) = Winv*(reshape(squeeze(ROIimg(n,:,:,:)),nSteps,ROI*ROI)./repmat(AmpArray(n,:),[nSteps 1]));
end

for n = 1:nLambda
    mmVecs(n,:,:,:) = reshape(squeeze(M(n,:,:)),16,ROI,ROI);
end

end

