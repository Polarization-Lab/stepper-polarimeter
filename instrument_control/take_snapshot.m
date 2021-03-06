function [image, ref] = take_snapshot(vid, exposure , framesPerTrigger)
%take_snapshot takes image and reference data
%   vid video


% Create the object.
ai = daq('ni') ;

% Add one channel for recording the reference
addinput(ai,'Dev2','ai1','Voltage');

% Set the sample rate to 1000 Hz.
ai.Rate = 1000;
% determine how many measurements to aquire 
ref_size = 1 * 1000;

%calculate container sizes
start(vid)
im = getdata(vid);
imdim = size(im);
image = zeros(imdim);
ref = zeros(1,ref_size);
stop(vid)

%aquire image simultaneously
start(vid)
im = getdata(vid);
ref = read(ai, ref_size, "OutputFormat", "Matrix");
stop(vid);

im = double(im);

%calculate average image
image = im(:,:,1,1);
for i = 2:framesPerTrigger
    image = image + im(:,:,1,i);
end
image = image/framesPerTrigger;
% imagesc(image)
%calculate reference amplitude
ref = mean(ref);

end

