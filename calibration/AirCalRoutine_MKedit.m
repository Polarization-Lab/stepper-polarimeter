function [Irrad,W] = AirCalRoutine_MKedit(amp,PSG_delta,PSG_theta,PSA_delta,PSA_theta,PSA_LP,PSGangles,PSAangles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% for kk = 1:length(PSGangles)
%     temp = LP(PSA_LP)*LR(PSA_delta, PSAangles(kk) + PSA_theta)*LR(PSG_delta, PSGangles(kk) + PSG_theta)*LP(0)*[amp;0;0;0];
%     %Irrad(kk) = temp(1);
% end


for kk = 1:length(PSGangles)
        psgMM = LR(PSG_delta, -PSGangles(kk) + PSG_theta)*LP(0);
        psaMM = LP(PSA_LP)*LR(PSA_delta, -PSAangles(kk) + PSA_theta);
        a = psaMM(1,:);%analyzer vector 1x4
        g = psgMM(:,1);%polarizance vector 4x1
        w = kron(a,g');%1x16 vector 
        Irrad(kk) = w*amp*reshape(eye(4,4),16,1);%air measurements for given W matrix 
        W(kk,:) = w;
end



return

