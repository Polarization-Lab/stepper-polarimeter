function [image_data] = load_imagedata(fn,wavelength,meas_num,xstart,xcount,ystart,ycount)
%LOAD_IMAGEDATA this function loads image data from the hdf5 file for a
%single wavelength. Will average over apporipriate number of frames and
%display associated metadata
%   fn is the filename
%   wavelength is the wavelength in nm
%   meas_num is the number of measurements in the file
%   ROI is the roi data to save [xstart,xstop,ystart,ystop]


%build the group name
group_name = strcat('/images/wave',num2str(wavelength));

%check if data is present and load number of frames, exposure
info = h5info(fn,'/images');
trace_groups = info.Groups;
groups_names = {trace_groups.Name};

if ismember(group_name,groups_names)
    %exposure
    exp=h5readatt(fn,group_name,'exposure_time');
    disp('the exposure time is [sec]:')
    disp(exp)
    
    %fames taken 
    fm=h5readatt(fn,group_name,'frames_per_trigger');
    disp('frames per measurement:')
    disp(fm)
    
else
    disp('this wavelength is not in the .h5 file')
end


%allocate space 
image_data = zeros(xcount,ycount,meas_num);
start = [xstart,ystart,1,1];
count = [xcount,ycount,1,fm];

for i = 1:meas_num
    group_name   = strcat('/images/wave',num2str(wavelength),'/meas',num2str(i));
    images = h5read(fn,strcat(group_name,'/imagedata'),start,count);

    
    %take average of frames
    mean_img = images(:,:,1); 
    for j = 2:fm
        mean_img = mean_img + images(:,:,j);
    end
    mean_img = mean_img./fm;
    
    c=imagesc(mean_img);
    colorbar
    
    image_data(:,:,i) = mean_img ;
    
    clear images;clear mean_img
end

