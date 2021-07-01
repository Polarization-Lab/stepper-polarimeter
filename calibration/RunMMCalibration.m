calibrationFilePath = 'D:/Measurements/Air_Calibrations/BVO-UV_PSA-LP/BVO-UV_PSA-LPjune302021-30-Jun-2021.h5'; 
xBegROI = 751;
yBegROI = 751;
xEndROI = 760;
yEndROI = 760;
LambdaList = [400:10:800];
nLambda = 41;
nSteps = 64;

[caldata,W,NRMSD,PixelCount] = MM_Calibration(calibrationFilePath,xBegROI,xEndROI,yBegROI,yEndROI,LambdaList,nSteps)

%save calibration data
%savePath = 'D:\Measurements\Air_Calibrations\Summer21\Calibration_Data'

