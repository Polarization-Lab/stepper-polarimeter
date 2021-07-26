function Irrad = GenAirMeasFromStepperParsV2(amp,angles,nSteps,LambdaList)
%This function accepts parameters of the Stepper instrument and returns the air measurement defined by those parameters

%INPUT VARIABLES (Stepper parameters, all angles in radians)
%amp: nLambdaxROI_widthxROI_width
%a_PSG: 1x1 slope of retardance magnitude dispersion
%b_PSG: 1x1 intercept of retardance magnitude dispersion
%a_PSA: 1x1 slope of retardance magnitude dispersion
%b_PSA: 1x1 intercept of retardance magnitude dispersion
%PSG_theta: 1x1
%PSA_theta: 1x1
%PSA_LP: 1x1
%OUTPUT VARIABLES (measurement defined by stepper parameters)
%Iradd: nLambdaxnStepsxROI_widthxROI_width

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nLambda=length(LambdaList);%Number of wavelengths measured
start=1;
a_PSG = angles(start);               %PSG retardance magnitude slope and intercept
start=start+1;
b_PSG = angles(start);               
start=start+1;
a_PSA = angles(start);               %PSA retardance magnitude slope and intercept
start=start+1;
b_PSA = angles(start);               
start=start+1;
PSG_theta = angles(start);                               %PSG retarder orientaiton NOT pixel or wavelength dependent
start=start+1;
PSA_theta = angles(start);                               %PSA retarder orientaiton NOT pixel or wavelength dependent
start=start+1;
PSA_LP = angles(start);                                  %PSA linear Polarizer theta

ROI_width = sqrt(length(amp(:))/nLambda);
amp = reshape(amp,nLambda,ROI_width,ROI_width);
ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 
ThetaMotorAna = 4.9*ThetaMotorGen;

%OUTPUT VARIABLE
Irrad=zeros(nLambda,nSteps,ROI_width,ROI_width);

for n=1:nLambda
    
    PSG_delta=a_PSG*LambdaList(n)+b_PSG;
    PSA_delta=a_PSA*LambdaList(n)+b_PSA;
    
    for kk = 1:nSteps
    
        psgMM = LR(PSG_delta, ThetaMotorGen(kk) + PSG_theta)*LP(0);
        psaMM = LP(PSA_LP)*LR(PSA_delta, ThetaMotorAna(kk) + PSA_theta);
        a = psaMM(1,:);%analyzer vector 1x4
        g = psgMM(:,1);%polarizance vector 4x1
        w = kron(a,g');%1x16 vector 
        Irrad(n,kk,:,:) = squeeze(amp(nLambda,:,:))*(w*reshape(eye(4,4),16,1));%amp is factored out of W due to linearity  
        
    end
end

return

