function [] = read_wavelength_h5(filename,wavelength)
%READ_WAVELENGTHH5  outputs data sotred in wavelength file
%  show measurement parameters for wavelength

%build the group name
group_name = strcat('/images/wave',num2str(wavelength));

%check if data is present
info = h5info(filename,'/images');
trace_groups = info.Groups;
groups_names = {trace_groups.Name};

if ismember(group_name,groups_names)
    %exposure
    exp=h5readatt(filename,group_name,'exposure_time');
    disp('the exposure time is [sec]:')
    disp(exp)
    
    %fames taken 
    fm=h5readatt(filename,group_name,'frames_per_trigger');
    disp('frames per measurement:')
    disp(fm)
    
else
    disp('this wavelength is not in the .h5 file')

end

end

