function imagevideo = RawImage_movie(Directory,file,lambdas,num_meas)
%Directoy - to h5 file where measurement data is held
%file  - name of .h5 file to be read, input without extension
%lambdas - single value or vector of wavelengths to show images for  

%Directory =  'D:\Measurements\Waveguides\option1\'; %change 
%file = 'SampleMeas-460_530_625-05-Jul-2021'
%lambdas = 530;

ThetaMotorGen = rad2deg((0:num_meas-1)*2*pi/num_meas); 
ThetaMotorAna = rem(4.9*ThetaMotorGen,360);

imagevideo = VideoWriter(strcat(Directory,file,'-Rawimagemovie.mp4'),'MPEG-4');
imagevideo.FrameRate = 16/15;
open(imagevideo);
set(0,'DefaultFigureVisible', 'off');
file = strcat(file,'.h5');

for jj = 1:length(lambdas)
for ii = 1:num_meas 
    group_name   = strcat('/images/wave',num2str(lambdas(jj)),'/meas',num2str(ii));
    images = h5read(strcat(Directory,file),strcat(group_name,'/imagedata'));
	imagesc(images,[0 2^16]);
    colormap(gray);
    colorbar;
	title(sprintf('$$\\overline{\\lambda} =%d[nm],\\theta_{PSG}=%d^\\circ,\\theta_{PSA}=%d^\\circ$$',lambdas(jj),round(ThetaMotorGen(ii)),round(ThetaMotorAna(ii))),'interpreter','latex','FontSize',20);
    axis equal; %set aspect ratio 1:1
    axis off; %get rid of pixel count
    %xlabel(strcat('Exposure Time =',num2str(Exps(jj)),'[s]')) ;
    frame = getframe(gcf);
    writeVideo(imagevideo,frame);
end
hold off
end

set(0,'DefaultFigureVisible', 'on');


close(imagevideo);
