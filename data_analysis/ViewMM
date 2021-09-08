function ViewNormMM( I )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~,~,X,Y]=size(I);
I=reshape(I,16,X,Y);
figure;
I=I./(I(1,:,:));

% for n=1:16
% subplot(4,4,n);
% imagesc(squeeze(I(n,:,:)));colormap gray;axis equal; axis off;
% end

    tempImg=[];
    for n=1:4
        tempRow=[];
        for zz=1:4
            tempSubImg=squeeze(I(zz+(n-1)*4,:,:));
            tempSubImg=[min(I(:))*zeros(X,1)';tempSubImg;min(I(:))*zeros(X,1)'];%Put borders on image
            tempSubImg=[min(I(:))*zeros(Y+2,1),tempSubImg,min(I(:))*zeros(Y+2,1)];
            tempRow=[tempRow;tempSubImg];
        end
        tempImg=[tempImg,tempRow];
    end
    I=reshape(I,16,X*Y);
    t=max(max(abs(I(2:end,:))));
    imagesc(tempImg,[-t,t]);axis square; axis off; colormap(GWP);colorbar;
    
    

end

