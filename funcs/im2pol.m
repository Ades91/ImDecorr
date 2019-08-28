function imP = im2pol(imR)

rMin=0; 
rMax=1;

[Mr Nr] = size(imR); % size of rectangular image 
xRc = (Mr+1)/2; % co-ordinates of the center of the image 
yRc = (Nr+1)/2; 
sx = (Mr-1)/2; % scale factors 
sy = (Nr-1)/2;

M = 2*Mr; N = 2*Nr;
% imP = zeros(M, N);

dr = (rMax - rMin)/(M-1); 
dth = 2*pi/N;

% loop in radius and 
r=rMin:dr:rMin+(M-1)*dr; 
th=(0:dth:(N-1)*dth)'; 
[r,th]=meshgrid(r,th); 
x=r.*cos(th); 
y=r.*sin(th); 
xR = x*sx + xRc; 
yR = y*sy + yRc; 
imP = interp2(imR, xR, yR,'cubic'); %interpolate (imR, xR, yR);

imP(isnan(imP)) = 0;

end