% apodize the input image in with a cosine edge mask over a length of N
% pixels
function [out,mask] = apodImRect(in,N)

Nx = min(size(in,1),size(in,2));

x = abs(linspace(-Nx/2,Nx/2,Nx));
map = x > Nx/2 - N;

val = mean(mean(in(:)));

d = (-abs(x)- mean(-abs(x(map)))).*map;
d = linmap(d,-pi/2,pi/2);
d(not(map)) = pi/2;
mask = (sin(d)+1)/2;

% make it 2D
if size(in,1) > size(in,2)
    mask = mask.*imresize(mask',[size(in,1) 1],'bilinear');
elseif size(in,1) < size(in,2)
    mask = imresize(mask,[1 size(in,2)],'bilinear').*mask';
else
    mask = mask.*mask';
end
	
out = (in-val).*mask + val;