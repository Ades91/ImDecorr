% [kcMax,A0,kcGM,d0,d] = getDcorr(im,r,Ng,figID)
% ---------------------------------------
%
% Estimate the image cut-off frequency based on decorrelation analysis
%
% Inputs:
%  im        	2D image to be analyzed
%  r           	Fourier space sampling of the analysis (default: r = linspace(0,1,50)
%  Ng			Number of high-pass filtering (default: Ng = 10)
%  figID		If figID > 1, curves will be plotted in figure(figID)
%
% Outputs:
%  kcMax        Estimated cut-off frequency of the image in normalized frequency
%  A0			Amplitude of the local maxima of d0
%  kcGM			Estimated cut-off frequency using Geometric-Mean metric
%  d0 			Decorrelation function before high-pass filtering
%  d			All decorrelation functions
%
% ---------------------------------------
%
% A detailled description of the method can be found in : 
% "Descloux, A., K. S. Grußmayer, and A. Radenovic. "Parameter-free image 
% resolution estimation based on decorrelation analysis."
% Nature methods (2019): 1-7."
% 
%   Copyright © 2018 Adrien Descloux - adrien.descloux@epfl.ch,
%   Ecole Polytechnique Federale de Lausanne, LBEN,
%   BM 5.134, Station 17, 1015 Lausanne, Switzerland.
%
%  	This program is free software: you can redistribute it and/or modify
%  	it under the terms of the GNU General Public License as published by
% 	the Free Software Foundation, either version 3 of the License, or
%  	(at your option) any later version.
%
%  	This program is distributed in the hope that it will be useful,
%  	but WITHOUT ANY WARRANTY; without even the implied warranty of
%  	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  	GNU General Public License for more details.
%
% 	You should have received a copy of the GNU General Public License
%  	along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [kcMax,A0,kcGM,d0,d] = getDcorr(im,r,Ng,figID)

if nargin < 4; figID = 0; end
if nargin < 3; Ng = 10;  end
if nargin < 2; r = linspace(0,1,50); end

% input check
if length(r) < 30
    r = linspace(min(r),max(r),30);
end
if Ng < 5
    Ng = 5;
end
%%
im = single(im);
im = im(1:end-not(mod(size(im,1),2)),1:end-not(mod(size(im,2),2))); % odd number of pixels

[X,Y] = meshgrid(linspace(-1,1,size(im,2)),linspace(-1,1,size(im,1)));
R = sqrt(X.^2 + Y.^2);
Nr = length(r);
if isa(im,'gpuArray')
    r = gpuArray(r); 
    R = gpuArray(R);
end
% In : Fourier normalized image
In = fftshift(fftn(fftshift(im))); In = In./abs(In); In(isinf(In)) = 0; In(isnan(In)) = 0;
mask0 = R.^2 < 1^2;
In = mask0.*In; % restric all the analysis to the region r < 1

if figID
    fprintf('Computing dcorr: ');
end
% Ik : Fourier transform of im
Ik = mask0.*fftshift(fftn(fftshift(im)));

c = sqrt(sum(sum(Ik.*conj(Ik))));

r0 = linspace(r(1),r(end),Nr);

for k = length(r0):-1:1
    cc = getCorrcoef(Ik,(R.^2 < r0(k)^2).*In,c);
    if isnan(cc); cc = 0; end
    d0(k) = gather(cc); % gather if input image is gpuArray 
end

[ind,snr0] = getDcorrMax(d0);
res0 = gather(r(ind));

gMax = 2/r0(ind);
if isinf(gMax); gMax = max(size(im,1),size(im,2))/2;end

% search of highest frequency peak
g = [size(im,1)/4, exp(linspace(log(gMax),log(0.15),Ng))];

d = zeros(Nr,2*Ng); kc = res0; SNR = snr0; gm = res0;

for refin = 1:2 % two step refinement

    for h = 1:length(g)
        Ir = Ik.*(1 - exp(-2*g(h)*g(h)*R.^2)); % Fourier Gaussian filtering
        c = sqrt(sum(sum(abs(Ir).^2)));
        
        for k = length(r):-1:1

            if isa(im,'gpuArray')
                cc = getCorrcoef(Ir,In.*(R.^2 < r(k)^2),c);
                if isnan(cc); cc = 0; end
                d(k,h + Ng*(refin-1)) = gather(cc);
            else
                mask = (R.^2 < r(k)^2);
                cc = getCorrcoef(Ir(mask),In(mask),c);
                if isnan(cc); cc = 0; end
                d(k,h + Ng*(refin-1)) = cc;
            end

        end
        
        [ind,snr] = getDcorrMax(d(:,h + Ng*(refin-1)));
        kc(end+1) = gather(r(ind));
        SNR(end+1) = snr;
        gm(end+1) = sqrt(snr*kc(end));
        if figID
        	fprintf('-');
        end
    end

% refining the high-pass threshold and the radius sampling
    if refin == 1

    % high-pass filtering refinement
        [~,indgm] = max(gm);
        [~,indmax] = max(kc);
        ind1 = min(indgm,indmax);
        ind2 = max(indgm,indmax);
        if ind1 == 1
            ind1 = 1;
            ind2 = 1;
        elseif ind2 == length(g)
            ind2 = length(g)-1;
        end
        g1 = g(ind1); g2 = g(ind2+1);
        g = exp(linspace(log(g1),log(g2),Ng+2));
        g = g(2:end-1);
        
        % radius sampling refinement
        [~,ind] = max(kc);
        r1 = kc(ind)-0.05; r2 = kc(ind)+0.3;
        if r1 < 0 ; r1 = 0; end
        if r2 > 1; r2 = 1; end
        r = linspace(r1,min(r2,r(end)),Nr);
        r2 = r;
    end
end
if figID
    fprintf(' -- Computation done -- \n');
end
radAv = getRadAvg(gather(log(abs(Ik)+1)));

% add d0 results to the analysis (usefull for high noise images)
kc(end+1) = gather(res0);
SNR(end+1) = snr0;

% % need at least 0.05 of SNR to be even considered
kc(SNR < 0.05) = 0;
SNR(SNR < 0.05) = 0;

snr = SNR;

% output results computation
if ~isempty(kc)
    % highest resolution found 
    [kcMax,ind] = max(kc);
    AMax = SNR(ind);

    % compute the geometric mean to determine the best res/SNR curve
    gm = sqrt(kc.*SNR);
    [~,ind] = max(gm);

    kcGM = kc(ind);
    A0 = snr0; % average image contrast has to be estimated from original image
else
    kcMax = r(2);
    Amax = 0;
    res = r(2);
    A0 = 0;
end

% results display if figID specified
if figID
    lnwd = 1.5;
    figure(figID);
    plot(r0,d(:,1:Ng),'color',[0.2 0.2 0.2 0.5]);
    hold on
    radAv(1) = radAv(2); %for plot 
    radAv(end) = radAv(end-1);
    plot(linspace(0,1,length(radAv)),linmap(radAv,0,1),'linewidth',lnwd,'color',[1 0 1])
    for n = 1:Ng
        plot(r2,d(:,n+Ng),'color',[0 0 (n-1)/Ng])
    end
    plot(r0,d0,'linewidth',lnwd,'color','g')
    plot([kcMax kcMax],[0 1],'k')
    for k = 1:length(kc)
        plot(kc(k),snr(k),'bx','linewidth',1)
    end
    hold off
    title(['Dcor analysis : res ~ ',num2str(kcMax,4),', SNR ~ ',num2str(A0,4)])
    xlim([0 1]); ylim([0 1])
    xlabel('Normalized spatial frequencies')
    ylabel('C.c. coefficients')
end