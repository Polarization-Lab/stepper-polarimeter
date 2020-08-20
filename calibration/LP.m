function LP = LP(Theta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Theta = Theta*pi/180;
LP = [0.5,0.5*cos(2*Theta),0.5*sin(2*Theta),0;
    0.5*cos(2*Theta), 0.5*(cos(2*Theta)^2), 0.5*cos(2*Theta)*sin(2*Theta),0;
    0.5*sin(2*Theta), 0.5*cos(2*Theta)*sin(2*Theta),0.5*(sin(2*Theta)^2),0;
    0,0,0,0];
end

