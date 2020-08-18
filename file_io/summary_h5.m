function [] = summary_h5(filename)
%summary_h5 this function will return the contents of a stepper HDF5 files
%inputs:
%   filename string of filename and path


%collect and display wavelengths in file 
info = h5info(filename,'/images');
trace_groups = info.Groups;
groups_names = {trace_groups.Name};

disp('the wavelength groups present are')
disp(groups_names)

%display info about PSG/PSA
psg=h5readatt(filename,'/images','PSG_positions');
num = size(psg);
disp('number of PSA/PSG orientation measurements')
disp(num(1))

%display user notes
notes=h5readatt(filename,'/images','user_notes');
disp('user notes :')
disp(notes)
end

