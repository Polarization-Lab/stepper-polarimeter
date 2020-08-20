function LR = LR(delta,Theta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Theta = Theta*pi/180;
% LR = [1,0,0,0;
%     0,(cos(2*Theta)^2)+cos(delta)*sin(2*Theta)^2,(1-cos(delta))*cos(2*Theta)*sin(2*Theta),-sin(delta)*sin(2*Theta);
%     0,(1-cos(delta))*cos(2*Theta)*sin(2*Theta),cos(delta)*(cos(2*Theta)^2)+(sin(2*Theta)^2),cos(2*Theta)*sin(delta);
%     0,sin(delta)*sin(2*Theta),-cos(2*Theta)*sin(delta),cos(delta)];
LR = [1,0,0,0;
    0,(cos(2*Theta)^2)+cos(delta)*sin(2*Theta)^2,(1-cos(delta))*cos(2*Theta)*sin(2*Theta),-sin(delta)*sin(2*Theta);
    0,(1-cos(delta))*cos(2*Theta)*sin(2*Theta),cos(delta)*(cos(2*Theta)^2)+(sin(2*Theta)^2),cos(2*Theta)*sin(delta);
    0,sin(delta)*sin(2*Theta),-cos(2*Theta)*sin(delta),cos(delta)];
end
