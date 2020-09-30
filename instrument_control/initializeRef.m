function [dq,ch] = initializeRef()
%initializeRef opens NI daq channel
%   Detailed explanation goes here

dq = daq("ni");
ch = addinput(dq,"Dev2", "ai0","Voltage");

dq.Rate = 100;


end

 