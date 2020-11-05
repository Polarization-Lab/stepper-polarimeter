function [Irrad,W] = AirCalFit(amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,nSteps)

[~,sizeX]=size(amp);
amp = reshape(amp,[1,sqrt(sizeX),sqrt(sizeX)]);
ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 
ThetaMotorAna = 4.9*ThetaMotorGen;

%OUTPUT VARIABLE
Irrad=zeros(nSteps,sqrt(sizeX),sqrt(sizeX));
W=zeros(nSteps,16);

%INPUT VARIABLE
%amp nLambdaxsizeXxsizeY
%PSG_delta 1xlambda
%PSA_delta 1xlambda
%PSG_theta 1x1
%PSA_theta 1x1
%PSA_LP 1x1

for kk = 1:nSteps
    psgMM = LR(PSG_delta, ThetaMotorGen(kk) + PSG_theta)*LP(0);
    psaMM = LP(PSA_LP)*LR(PSA_delta, ThetaMotorAna(kk) + PSA_theta);
    a = psaMM(1,:);%analyzer vector 1x4
    g = psgMM(:,1);%polarizance vector 4x1
    w = kron(a,g');%1x16 vector 
    Irrad(kk,:,:) = squeeze(amp(1,:,:))*(w*reshape(eye(4,4),16,1));%amp is factored out of W due to linearity  
    W(kk,:) = w;

end


return

