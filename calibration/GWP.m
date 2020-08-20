function h = GWP(m)

if nargin < 1, m = size(get(gcf,'colormap'),1); end
x0=[0,1,2,3,4]/4;
ColorCore= [0.5 0 0.5; 1 0 1; 1 1 1; 0 1 0.5; 0 0.5 0];
if mod(m,2)==1,
    h=interp1(x0,ColorCore,[0:1/(m-2):1]);
else
    h1=interp1(x0,ColorCore,[0:1/(m-2):0.5]);
    h2=interp1(x0,ColorCore,0.5+[0:1/(m-2):0.5]);
    h = [h1;h2];
end