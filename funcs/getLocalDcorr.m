% Author : Adrien Descloux
% Image quality estimation based on a decorrelation analysis
%
% The analysis is parameter free and allows independent estimation of 
% image resolution and SNR.
%
% [kcMap,A0Map] = getLocalDcorr(im,tileSize,tileOverlap,r,Ng,figID)
% 
% optional input :
%       r : specify the range and sampling of dcorr analysis
%       N : Number of high-pass image
%       figID : if > 0, display results in figure(figID)

function [kcMap,A0Map] = getLocalDcorr(im,tileSize,tileOverlap,r,Ng,figID)

if nargin < 7; figID = 0; end
if nargin < 6; Ng = 10; end
if nargin < 5; r = linspace(0,1,50);end

px = round(linspace(1,size(im,2),ceil(size(im,2)/(tileSize-tileOverlap))));
py = round(linspace(1,size(im,1),ceil(size(im,1)/(tileSize-tileOverlap))));
kcMap = zeros(length(py)-1,length(px)-1);
A0Map = kcMap;
for xx = 1:length(px)-1
	for yy = 1:length(py)-1
        subIm = im(py(yy):py(yy+1),px(xx):px(xx+1),1);
        [kc,A0] = getDcorr(subIm,r,Ng);
        kcMap(yy,xx) = kc;
        A0Map(yy,xx) = A0;
	end
end

if figID
    figure(figID)
    subplot(121)
        imagesc(kcMap); colorbar; title('kcMap')
    subplot(122)
        imagesc(A0Map); colorbar; title('A0Map')
end
    
    