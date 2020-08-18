function [COMdmm] = initializeDMM()
disp('Initializing DMM')
COMdmm=serialport('COM8' , 9600 ,...
    'DataBits' , 7 ,...
    'Parity' , 'even' ,...
    'StopBits' , 2 ,...
    'Timeout', 1);

configureTerminator(COMdmm,"LF")

% ^dmm RS-232 comm. settings must match with computer
%place DMM in known state 
writeline(COMdmm,'SYST:REM'); %Enables remote/computer commands

%input settings
writeline(COMdmm, 'CONF:VOLT:DC 10,0.0001') ;%Configures measurement for dc Volts, 10 volts max, and how many digits min

%set up triggering conditions
writeline(COMdmm, 'TRIG:SOUR BUS') ;    %Current trigger source is matlab software cmd,
writeline(COMdmm, 'SAMP:COUN 10')  ;    %number of samples/measurements taken when a trigger is sent
writeline(COMdmm, 'TRIG:COUN 1');       %How many triggers are to be sent/recieved by DMM

end

