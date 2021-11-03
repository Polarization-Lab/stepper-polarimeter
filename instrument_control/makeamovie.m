function video = makeamovie(Directory,file,lambdas,num_meas,Exps)
%Directoy - to h5 file where measurement data is held
%file  - name of .h5 file to be read, input without extension
%lambdas - single value or vector of wavelengths to show images for  

disp('Generating intensity plot Movie')
IntensityPlot_movie(Directory,file,lambdas,num_meas,Exps);

disp('Generating Reference detector plot Movie')
ReferencePlot_movie(Directory,file,lambdas,num_meas)

disp('Generating raw image movie')
RawImage_movie(Directory,file,lambdas,num_meas);

end


