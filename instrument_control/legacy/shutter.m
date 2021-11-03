function [successQ] = shutter(COMmono,openQ)
%SHUTTER opens or closes the shutter
%   COMmono is the com port to the monochromator
%   openQ is a boolean, true = open, false = close
%   successQ tells if operation was successful 

%assign serial command
if openQ
    cmd = 'OPEN';
    key = 'O';
else
    cmd = 'CLOSE';
    key = 'C';
end
    
%check if shutter is already in correct position
fprintf(COMmono,'SHUTTER?'); %check shutter position
status=fscanf(COMmono);  %reads echo
status=fscanf(COMmono);  %reads response


if(isempty(regexp(status,key,'once'))) %search response for desired key
    disp('moving shutter')
    fprintf(COMmono,join('SHUTTER ',cmd)); %Sends cmd to change shutter
    status=fscanf(COMmono);  %reads echo
    status=fscanf(COMmono);  %reads response
    pause(5)
    if(isempty(regexp(status,key,'once')))
        disp('shutter move unsucessful')
        successQ = false;
    else
        disp('shutter complete')
        successQ = true;
    end
        
else
    disp('shutter already in requested position') %if already in position, do nothing
    successQ = true;
end
    

end

