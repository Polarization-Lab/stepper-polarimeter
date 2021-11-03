function [] = GenImgfromMM(fn_sample,Si,So)

%Si: 4x1 of  incident Stokes vector, input polarization state
%fn_sample: path to MM sample file 
%function outputs image generated from MM from a given stokes vector 
%NOTE: change title of figure in line 51 for different input vectors

PSG=h5readatt(fn_sample,'/images','PSG_positions');
nSteps = length(PSG);

waves=h5info(fn_sample,'/images/');
waves={waves.Groups.Name};
LambdaList=zeros(1,length(waves));
nLambda=length(LambdaList);
for l=1:nLambda
    extractAfter(waves(l),'wave');
    LambdaList(l)=str2double(char(extractAfter(waves(l),'wave')));
end

info = h5info(fn_sample);
index=length({info.Groups.Name})-5; % because there are 4 components in each index and we want the first one

rect=h5readatt(fn_sample,strcat('/Mueller',num2str(index),'/'),'ROI');
nPixels=(rect(4)+1)^2; 
%% 


%read MM from file, reshapes to get 4x4 needed to do stokes analysis
%code might need to be adjusted for multiple wavelengths in one h5 file?
for n=1:nLambda
    I(n,:,:)=h5read(fn_sample,strcat('/Mueller',num2str(index),'/wave',num2str(LambdaList(n))));
    MM = reshape(I,[4 4 sqrt(nPixels) sqrt(nPixels)]);
end

%analyzes 4x1 incident and 4x1 exiting stokes vector on 4x4 mueller matrix for
%each pixel and compiles into image
for n = 1:sqrt(nPixels)
    for m = 1:sqrt(nPixels)
       %So = MM(:,:,n,m)*Si;%4x1
       image(n,m) = Si.'*MM(:,:,n,m)*So; %dimensions of ROI. Si is transposed
    end
end

%plotting image, change title depending on stokes vector

for n = 1:nLambda
    image = imrotate(image,90); %waveoptics has a different orientation than us
    imagesc(image,[0 4e5]);colormap gray;colorbar;axis off
    figure(n);
    imagesc(image);colormap gray;colorbar
    figure(n);
    title(sprintf('Verically Polarized Light: \\lambda=%d[nm]',LambdaList(n)),'FontSize',20);
end
