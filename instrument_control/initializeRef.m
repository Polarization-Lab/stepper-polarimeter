function [dq,ch] = initializeRef()
%initializeRef opens NI daq channel
%   Detailed explanation goes here

dq = daq("ni");
ch = addinput(dq,"Dev2", "ai0","Voltage");


ch0 = addoutput(dq,"Dev2","port1/line0:3","Digital");
ch1 = addinput(dq,"Dev2","port0/line0:3","Digital");

end

 