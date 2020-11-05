function J = M2J(mmVecs)
%%%%%%%%%%%%%%%%%%%%%%%%% SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%
% Returns Mueller and Coherency matrix generated for user. Also returns
% Mueller and coherency matrix of closest approximation to non-depolarizing
% Jones-Mueller matrix and assocoated Jones matrix
%%%%%%%%%%%%%%%%%%%%%%%%% INPUT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%% 
% Nterms: positive variable, number of terms in convex sum of Jones
% matrices to create Mueller matrix

%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
% M: Mueller matrix
% C: Coherency matrix
% lambda: elements of convex sum
% Jtot: Jones matrices in convex sum
% C1: rank 1 coherency matrix, truncated to create Mueller-Jones matrix
% MJ: Mueller-Jones matrix which is MNLS solution for non-depolarizing Mueller approximation to M
% J: Jones matrix associated with MJ

M=zeros(4,4);
C=zeros(4,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = [1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 sqrt(-1) -sqrt(-1) 0];
PI = zeros(4,4,4,4);
for n=1:4
	for m=1:4
		PI(n,m,:,:) = 0.5*U*(kron(reshape(U(n,:),2,2),conj(reshape(U(m,:),2,2))))*U';
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% lambda=rand(Nterms,1);
% lambda=lambda/sum(lambda);
% for n=1:Nterms
% 	%Transmittance condition - start:(EQ12 in
% 	%https://doi.org/10.1364/JOSAA.17.000328)- UNITARY equivalent? +
% 	%convex sum
% 	c=10;
%     
% 	while(c>1)
% 		J=.3*randn(2,2)+.3*1i*randn(2,2);
% 		tsq=abs(J).^2;
% 		c=0.5*(sum(tsq(:))+sqrt((tsq(1,1)-tsq(1,2)+tsq(2,1)-tsq(2,2))^2+(4*abs(J(1,1)*conj(J(1,2)+J(2,1)*conj(J(2,2))))^2)));
%     end
%     
%     Jtot(n,:,:)=J;
% 	%Transmittance condition - end
% 	term=lambda(n)*0.5*U*(kron(J,conj(J)))*U';
% 	M=M+term;
% end

M = reshape(mmVecs,31,4,4,10,10);

for ii = 1:10
    for jj = 1:10
        for n=1:4
            for m=1:4
                C(n,m,ii,jj) = 0.25*trace(reshape(PI(n,m,:,:),4,4)*squeeze(M(1,:,:,ii,jj)));
            end
        end
    end
end

for ii = 1:10
    for jj = 1:10
        [UC,S] = svd(C(:,:,ii,jj));
        UC_temp(ii*jj,:,:) = UC;
        S_temp(ii*jj,:,:) = S;
    end
end

J = zeros(100,2,2);

for ii=1:100
    J(ii,:,:)=sqrt(S(1,1))*reshape(squeeze(UC_temp(ii,:,1))*U,2,2)';
end

Jones = reshape(J,2,2,10,10);
A_xx2 = real(squeeze(Jones(1,1,:,:)));
A_xx2 = real(squeeze(Jones(1,2,:,:)));
A_xx2 = real(squeeze(Jones(2,1,:,:)));
A_xx2 = real(squeeze(Jones(2,2,:,:)));

A_xx2 = real(squeeze(Jones(1,1,:,:)));
A_xx2 = real(squeeze(Jones(1,1,:,:)));
A_xx2 = real(squeeze(Jones(1,1,:,:)));
A_xx2 = real(squeeze(Jones(1,1,:,:)));

JonesMatrix(1,:,:) = real(squeeze(Jones(1,1,:,:)));
JonesMatrix(2,:,:) = real(squeeze(Jones(1,2,:,:)));
JonesMatrix(3,:,:) = imag(squeeze(Jones(1,1,:,:)));
JonesMatrix(4,:,:) = imag(squeeze(Jones(1,2,:,:)));
JonesMatrix(5,:,:) = real(squeeze(Jones(2,1,:,:)));
JonesMatrix(6,:,:) = real(squeeze(Jones(2,2,:,:)));
JonesMatrix(7,:,:) = imag(squeeze(Jones(2,1,:,:)));
JonesMatrix(8,:,:) = imag(squeeze(Jones(2,1,:,:)));

Titles = ["A_{xx}","A_{xy}","\phi_{xx}","\phi_{xy}","A_{yx}","A_{yy}","\phi_{yx}","\phi_{yy}"] ;

for p = 1:8
    
subplot(2,4,p)
imshow(squeeze(JonesMatrix(p,:,:)),[],'colormap',GWP)
colorbar;
title(Titles(p));

end

C1=UC*diag([S(1,1), 0, 0, 0])*UC';%Rank 1 coherency matrix is the closest approximation of Mueller matrix to Mueller-Jones matrix
% 
MJ=zeros(4,4);
for n=1:4
	for m=1:4
		MJ(n,m) = trace(reshape(PI(n,m,:,:),4,4)*C1);
	end
end

for ii = 1:100
    MJ2(ii,:,:)=0.5*U*(kron(squeeze(J(ii,:,:)),conj(squeeze(J(ii,:,:)))))*U';%equivalent to MJ
end
% U=U';
% W=[U(1,n)+U(2,n) U(3,n)-1i*U(4,n); U(3,n)+1i*U(4,n) U(1,n)-U(2,n)];
% W*W'-W'*W;

return
