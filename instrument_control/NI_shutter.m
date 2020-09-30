function [] = NI_shutter(position)
%NI_SHUTTER this function will open and close the external shutter
%   Note this only works in R2020 and beyond
%   position - bool, if true shutter opens, close if F

%create a DataAquisition and Add Analog input channels 
dq = daq("ni");
addoutput(dq,"Dev2","port1/line0:3","Digital")

% Write values to output channels (turns on, waits for a return, turns off)
if position
    output = ones(1, 4);
else
    output = zeros(1, 4);
end

write(dq, output)
pause(.5)

end

