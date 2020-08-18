function [num] = read_DMM(COMdmm,exposure)
%READ_DMM Summary of this function goes here
%   Detailed explanation goes here
writeline(COMdmm, 'FETC?') %send internal memory to BUS
pause(exposure*2)
num = readline(COMdmm);
format long
num=str2num(num);

end

