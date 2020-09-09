function Irrad = AirCalRoutine_MKeditLSQ(amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,nSteps)
nLambda=length(PSG_delta);
sizeX = length(amp)/nLambda;
amp = reshape(amp,nLambda,sqrt(sizeX),sqrt(sizeX));
ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 
ThetaMotorAna = 4.9*ThetaMotorGen;
%OUTPUT VARIABLE
Irrad=zeros(nLambda,nSteps,sqrt(sizeX),sqrt(sizeX));

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
%amp is not in W to save space. To use Winv
%measurement=rand(nLambda,nSteps,sizeX,sizeY);
%M=zero(nLambda,15,sizeX*sizeY);
%[~,W]=AirCalRoutine_MKedit2(blah);
%for n=1:nLambda
%Winv=pinv(squeeze(W(n,:,:));
%M(n,:,:)=Winv*(reshape(measurements(k,:,:,:),nSteps,sizeX*sizeY)./repmat(amp(:)',[nSteps 1]));
%end
%ln38 the estimated amplitude of the source is in the demoninator of the measurements (pixel-by-pixel) before applying Winv b/c "amp" can be factored out of W due to linearity


return

