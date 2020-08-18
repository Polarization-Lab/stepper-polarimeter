function [] = prepare_DMM(COMdmm,exposure)
%PREPARE_DMM Summary of this function goes here
%   this must be executed BEFORE trigger_DMM


%match number of measurements with exposure time
num_meas = exposure * 100;

sample_command = ['SAMP:COUN ', num2str(num_meas)];

%set up triggering conditions
writeline(COMdmm, sample_command)  ;    %number of samples/measurements taken when a trigger is sent

%clear buffer and initialize
writeline(COMdmm, '*CLS');  
writeline(COMdmm, 'INIT');   
 

end
