function airMeasurement = LoadAirMeasurements(calibrationFile,LambdaList,nSteps)

disp('Loading measurements')

for ii = 1:length(LambdaList) %for loading all lambda(s)
%Get h5 reference data
    Lambda = LambdaList(ii);
    [~,~,refVecs(ii,:)] = load_refdata(calibrationFilePath,Lambda,nSteps);

%Grab measured data from h5 file
    airMeasurement(ii,:,:,:) = ReadCalImages(calibrationFilePath,Lambda,refVecs(ii,:),nSteps);
end

disp('Finished loading')