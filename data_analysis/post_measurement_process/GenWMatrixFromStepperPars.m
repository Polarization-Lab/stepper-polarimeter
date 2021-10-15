function [W,amp] = GenWMatrixFromStepperPars(fn_air,LambdaList)
%This function accepts an air measurement filepath, list of wavelengths,
%and number of retarder steps to return W matrix and amplitude fit per
%wavelength. 

%INPUT VARIABLES 
%fn: filepath to air measurement (which has been fit to system parameters)
%LambdaList: wavelengths of interest
%nSteps: usually 64 - number of retarder positions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nLambda=length(LambdaList);
PSG=h5readatt(fn_air,'/images','PSG_positions');
nSteps = length(PSG);


ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 
ThetaMotorAna = 4.9*ThetaMotorGen;

%OUTPUT VARIABLE
W=zeros(nLambda,nSteps,16);
%Variables that do not depend on wavelength
    a_PSG=h5readatt(fn_air,'/cal/','PSG_retardance_slope');
    b_PSG=h5readatt(fn_air,'/cal/','PSG_retardance_intercept');
    a_PSA=h5readatt(fn_air,'/cal/','PSA_retardance_slope');
    b_PSA=h5readatt(fn_air,'/cal/','PSA_retardance_intercept');
    PSG_theta=h5readatt(fn_air,'/cal/','PSG_theta');
    PSA_theta=h5readatt(fn_air,'/cal/','PSA_theta');
    PSA_LP=h5readatt(fn_air,'/cal/','PSA_LP');
    amp=h5readatt(fn_air,'/cal/','Amp');


for n=1:nLambda%Variables that DO depend on wavelength
    PSG_delta=a_PSG*LambdaList(n)+b_PSG;
    PSA_delta=a_PSA*LambdaList(n)+b_PSA;
    for kk = 1:nSteps
        psgMM = LR(PSG_delta, ThetaMotorGen(kk) + PSG_theta)*LP(0);
        psaMM = LP(PSA_LP)*LR(PSA_delta, ThetaMotorAna(kk) + PSA_theta);
        a = psaMM(1,:);%analyzer vector 1x4
        g = psgMM(:,1);%polarizance vector 4x1
        w = kron(a,g');%1x16 vector 
        W(n,kk,:) = w;        
    end
end

return

