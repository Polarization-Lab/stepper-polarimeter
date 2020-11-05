function [COMmono] = initializeMono()
%INITIALIZEMONO Summary of this function goes here
%Kira Hart khart@optics.arizona.edu
%July 28 2020
%This script is designed to initialize the monochromator
%and display the current wavelength

disp('Initializing Monochromator')

COMmono=serial('COM10');  %init serial comm w/ Monochromator
% set( COMmono , 'BaudRate' , 9600 ,...
%      'DataBits' , 8 ,...
%      'Parity' , 'none' ,...
%      'StopBits' , 1);
%     'Terminator' , 'CR/LF' );
% set( COMmono, 'Timeout', 1);
%%
fopen(COMmono);           %opens comm/VISA w/ Mono, can now send commands
fprintf(COMmono,'WAVE?'); %Sends current wavelength cmd, which is echoed
Mono_Status=fscanf(COMmono);  %reads echo
Mono_Status=fscanf(COMmono);  %reads wavelength (str)
disp('Wavelength is set at [nm]')
disp(Mono_Status)
Mono_Status=str2double(Mono_Status);    %converts to num

pause(2)
end

