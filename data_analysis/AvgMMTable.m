function TotalMM2 = AvgMMTable(mmVecs,Lambda,LambdaList)
%Only works with single lambda for now
%% Create averaged scalar MM (Table form)
tic
n = find(LambdaList == Lambda);
clear avgMM TotalMM
for p = 1:16
    avgMM(p) = mean(squeeze(mmVecs(n,p,:)))/mean(squeeze(mmVecs(n,1,:)),'all');
end
NormAvgMM = reshape(avgMM,[4 4])'

for p=1:16
    TotalMM(p) = mean(squeeze(mmVecs(n,p,:,:)),'all');
end
TotalMM2 = reshape(TotalMM,[4 4])'
toc