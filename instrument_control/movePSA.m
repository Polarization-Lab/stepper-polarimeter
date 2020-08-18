function [] = movePSA(xps,newPos)
%MOVEPSA this function changes the position of the PSA
% xps is the initialized XPS instance  
% newPos is the new angular position


%round value
pos = round(newPos,4);

%define the positioners
PSG     = 'Group1' ;
PSGpos  =  'Group1.Pos';
PSA     = 'Group2' ;
PSApos  =  'Group2.Pos';


disp('Moving PSA')
%moves PSA and waits til finished
G2=xps.GroupMoveAbsolute(PSA,pos,1);
[G2,CurrPos2]=xps.GroupPositionCurrentGet(PSA,1);
while double(CurrPos2) ~= pos
    pause(.1)
    [G2,CurrPos2]=xps.GroupPositionCurrentGet(PSA,1);
    CurrPos2=double(CurrPos2);
end

end

