function [num] = trigger_DMM(COMdmm)
%trigger DMM this function returns a reading from the reference detector
%   prepare_DMM must be executed before this 

writeline(COMdmm, '*TRG')
end



    
        
        
   