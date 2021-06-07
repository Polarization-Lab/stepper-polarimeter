%BWR Blue-White-Red colormap, mostly used for highly oscillatory data
%  bwr(M), a variant of  JET(M), is an M-by-3 matrix containing
%    the default colormap used by CONTOUR, SURF, PCOLOR and IMAGE.
%    The colors begin with dark blue, range through white and end with 
%    dark red.

% by Rodrigo S. Portugal (rosoport.matlab@gmail.com)
% last revision: 19/10/2012

function bkr = bkr(n)

if nargin < 1
   n = size(get(gcf,'colormap'),1);
end

k0 = round([1, n/4, n/2, 3*n/4, n]);

r0 = [255, 128,   0,   0,   25] / 255;
g0 = [ 25,  76,   0,  55,  25] / 255;
b0 = [  25,   0,   0, 128, 255] / 255;

k = 1:n;
r = interp1(k0, r0, k);
g = interp1(k0, g0, k);
b = interp1(k0, b0, k);

bkr = [r',g',b'];