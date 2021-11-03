function h = GWP(m)

if nargin < 1, m = size(get(gcf,'colormap'),1); end
x0=[0,1,2,3,4,5,6,7,8]/8;
ColorCore= [0.5 0 1; 0.625 0.27 1; 0.75 0.65 1; 0.875 1 1; 1 1 1; 0.875 1 0.875;0.65 1 0.65; 0.25 1 0.55; 0 1 0.5];
if mod(m,2)==1,
    h=interp1(x0,ColorCore,[0:1/(m-2):1]);
else
    h1=interp1(x0,ColorCore,[0:1/(m-2):0.5]);
    h2=interp1(x0,ColorCore,0.5+[0:1/(m-2):0.5]);
    h = [h1;h2];
end
