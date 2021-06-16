function [Irrad,W] = AirCalRoutine_MKedit2(amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,nSteps)
%%%%%%%%%%%%%%%%%%%%%%%%% SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%
% Returns the Irradiance values (although ln 48 currently commented) and unscaled W matrix parameterized by the properties of optical components given as input. The W matrix is not scaled by the amplitude to save space. The unscaled W matrix is the same for all pixels and only needs to be scaled by the amplitude.

%%%%%%%%%%%%%%%%%%%%%%%%% INPUT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%% 
% amp: nLambdaxsizeXxsizeX tensor of amplitude at each wavelength and each pixel

% PSG_delta: nLambdax1 vector retardance magnitude on generator side at each wavelength

% PSG_theta: scalar-valued orientation of retarder on generator side - not wavelength dependent

% PSA_delta: nLambdax1 vector retardance magnitude on analyzer side at each wavelength

% PSA_theta: scalar-valued orientation of retarder on analyzer side - not wavelength dependent

% PSA_LP: scalar-valued orientation of linear polarizer on analyzer side - not wavelength dependent

% nSteps: scalar-valued number of polarimetric measurements

%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
% Irrad: nLambdaxnStepsxnStepsxsizeXxsizeX tensor of polarimetric measurements which are images at each wavelength and PSA/PSG pair

% W: nLambdaxnStepsx16 tensor of polarimetric measurement matrices for each wavelength. NTB that the returned value of W does not depend on amplitude and the usage of the unscaled W is shown in ln54-64.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nLambda=length(PSG_delta);
% [~,sizeX]=size(amp);
sizeX = length(squeeze(amp(1,1,:)));
% amp = reshape(amp,nLambda,sqrt(sizeX),sqrt(sizeX));
ThetaMotorGen = (0:nSteps-1)*2*pi/nSteps; 
ThetaMotorAna = 4.9*ThetaMotorGen;
%OUTPUT VARIABLE
Irrad=zeros(nLambda,nSteps,sizeX,sizeX);
W=zeros(nLambda,nSteps,16);

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
%         Irrad(n,kk,:,:) = squeeze(amp(nLambda,:,:))*(w*reshape(eye(4,4),16,1));%amp is factored out of W due to linearity  
        W(n,kk,:) = w;
        
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

