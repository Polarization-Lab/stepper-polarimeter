function [] = NI_shutter(position)
%NI_SHUTTER this function will open and close the external shutter
%   Note this only works in R2020 and beyond
%   Until Hamamatsu driver works in R2020, this will need to be executed
%   externally to the general imaging chain
%   position - bool, if true shutter opens, close if F

%create a DataAquisition and Add Analog input channels 
s = daq.createSession('ni');
addDigitalChannel(s,'dev2','Port1/Line0:3','OutputOnly');


% Write values to output channels (turns on, waits for a return, turns off)
if position
    output = [1 1 1 1];
else 
    output = [0 0 0 0];
end

outputSingleScan(s, output)
pause(5)
release(s)

end

