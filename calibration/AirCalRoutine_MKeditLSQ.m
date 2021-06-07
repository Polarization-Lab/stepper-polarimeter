function Irrad = AirCalRoutine_MKeditLSQ(amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,nSteps)
nLambda=length(PSG_delta);
sizeX = length(amp)/nLambda;
amp = reshape(amp,nLambda,10,10);
ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 
ThetaMotorAna = 4.9*ThetaMotorGen;
%OUTPUT VARIABLE
Irrad=zeros(nLambda,nSteps,10,10);

%INPUT VARIABLE
%amp nLambdaxsizeXxsizeY
%PSG_delta 1xlambda
%PSA_delta 1xlambda
%PSG_theta 1x1
%PSA_theta 1x1
%PSA_LP 1x1

for n=1:nLambda
    for kk = 1:nSteps
    
        psgMM = LR(PSG_delta(n), ThetaMotorGen(kk) + PSG_theta)*LP(0);
        psaMM = LP(PSA_LP)*LR(PSA_delta(n), ThetaMotorAna(kk) + PSA_theta);
        a = psaMM(1,:);%analyzer vector 1x4
        g = psgMM(:,1);%polarizance vector 4x1
        w = kron(a,g');%1x16 vector 
        Irrad(n,kk,:,:) = squeeze(amp(nLambda,:,:))*(w*reshape(eye(4,4),16,1));%amp is factored out of W due to linearity  
        
    end
end
return

