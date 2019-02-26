% Author : Adrien Descloux
% Image quality estimation based on a decorrelation analysis
%
% The analysis is parameter free and allows independent estimation of 
% image resolution and SNR.
%
% [kcMax,A0,kcGM,d0,d] = getDcorr(im,r,Ng,figID)
% 
% optional input :
%       r : specify the range and sampling of dcorr analysis
%       Ng : Number of high-pass image
%       figID : if > 0, display results in figure(figID)

function [kcMax,A0,kcGM,d0,d] = getDcorr(im,r,Ng,figID)

if nargin < 4; figID = 0; end
if nargin < 3; Ng = 10;  end
if nargin < 2; r = linspace(0,1,50); end

% input check
if length(r) < 30
    r = linspace(min(r),max(r),30);
end
if Ng < 10
    Ng = 10;
end

im = double(im);
if mod(size(im,1),2) == 0
    im = im(1:end-1,1:end-1); % odd number of pixels
end

if figID
    hwait = waitbar(0,'Computing dcorr');
end

[X,Y] = meshgrid(linspace(-1,1,size(im,2)),linspace(-1,1,size(im,1)));
Nr = length(r);

I = fftshift(fftn(fftshift(im))); I = I./abs(I); I(isinf(I)) = 0; I(isnan(I)) = 0;
mask0 = X.^2 + Y.^2 < 1^2;
I = mask0.*I; % restric all the analysis to the region r == 1

Ik = mask0.*fftshift(fftn(fftshift((im-mean(im(:)))/std(im(:)))));
imr = real(ifftshift(ifftn(ifftshift(Ik))));

% compute dcorr 0 and find its maxima
imt = imr;
Ir = mask0.*fftshift(fftn(fftshift(imt)));
c = sqrt(sum(sum(abs(Ir).^2)));

% t0 = tic;
count = 0;
r0 = linspace(0,1,Nr);
for k = length(r0):-1:1
    rt = r0(k);
    mask  = X.^2 + Y.^2 < rt^2;
        
    temp = mask.*I;
    temp((end+1)/2,(end+1)/2) = 0; % remove the mean
%         cc = sum(sum(abs(Ir.*temp)))/(c.*sqrt(sum(sum(abs(temp).^2)))); % compute the correlation coeficient
    cc = getCorrcoef(Ir,temp,c);
    if isnan(cc); cc = 0; end
        d0(k) = gather(cc); % gather if input image is gpuArray 
        count = count +1;
        if figID
            waitbar(0.1*count/Nr,hwait);
        end
end

[ind,snr0] = getDcorrMax(d0);
res0 = r(ind);

% t1 = toc(t0);disp(['D0 computation time : ',num2str(t1)])

gMax = 2/r0(ind);
if isinf(gMax); gMax = max(size(im,1),size(im,2))/2;end

% automatic search of best geometric mean
g = exp(linspace(log(gMax),log(0.15),Ng));
d = []; kc = []; SNR = []; gm = []; dc = 1;

for refin = 1:2 % two step refinement
for h = 1:length(g)
    imt = imr - imgaussfilt(imr,g(h));
    Ir = mask0.*fftshift(fftn(fftshift(imt)));
    c = sqrt(sum(sum(abs(Ir).^2)));
    for k = length(r):-1:1
        rt = r(k);
        mask  = X.^2 + Y.^2 < rt^2;
        temp = mask.*I;
        temp((end+1)/2,(end+1)/2) = 0; % remove the mean
%         cc = sum(sum(abs(Ir.*temp)))/(c.*sqrt(sum(sum(abs(temp).^2)))); % compute the correlation coeficient
        cc = getCorrcoef(Ir,temp,c);
        if isnan(cc); cc = 0; end
        d(k,dc) = gather(cc); % gather if input image is gpuArray
        count = count+1;
        if figID
            waitbar(0.1 + 0.9*count/(Nr*N*2),hwait);
        end
    end
    dc = dc + 1;
    [ind,snr] = getDcorrMax(d(:,end));
    kc(end+1) = r(ind);
    SNR(end+1) = snr;
    gm(end+1) = sqrt(snr*r(ind)); % 
end

% refining the high-pass threshold and the radius sampling
if refin == 1
    
    % remove unphysical peaks from refinement anaysis
    kc(SNR < 0.05) = 0;
    SNR(SNR < 0.05) = 0;

    gm(SNR == 0) = 0;
    gm(kc == 0) = 0;
    
    % high-pass filtering refinement
    [~,indgm] = max(gm);
    [~,indmax] = max(kc);
    ind1 = min(indgm,indmax);
    ind2 = max(indgm,indmax);
    if ind1 == 1
        ind1 = 2;
    elseif ind2 == length(g)
        ind2 = length(g)-1;
    end
    g1 = (g(ind1)); g2 = (g(ind2+1));
    g = exp(linspace(log(g1),log(g2),Ng));
    
    % radius sampling refinement
    [~,ind] = max(gm);
    r1 = kc(ind)-0.05; r2 = kc(ind)+0.2;
    if r1 < 0 ; r1 = 0; end
    if r2 > 1; r2 = 1; end
    r = linspace(r1,r2,Nr);
    r2 = r;
end
end
% release GPU memory

radAv = getRadAvg(gather(log(abs(Ik)+1)));
    
% gpuDevice(1);

% % keep only refined data
% kc = kc(N+1:end);
% SNR = SNR(N+1:end);

% add d0 results to the analysis (usefull for high noise images)
kc(end+1) = res0;
SNR(end+1) = snr0;

% need at least 0.05 of SNR to be even considered
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

if figID
    waitbar(1,hwait);
    delete(hwait)
end

% results display if figID specified
if figID
    lnwd = 1.5;
    figure(figID);
    r0 = linspace(0,1,size(d,1));
    plot(r0,d(:,1:Ng),'color',[0.2 0.2 0.2 0.5]);
    hold on
    radAv(1) = radAv(2); %for plot 
    radAv(end) = radAv(end-1);
    plot(linspace(0,1,length(radAv)),linmap(radAv,0,1),'linewidth',lnwd,'color',[1 0 1])
    for n = 1:Ng
        plot(r2,d(:,n+Ng:end),'color',[0 0 (n-1)/Ng])
    end
    plot(r0,d0,'linewidth',lnwd,'color','g')
    plot([kcMax kcMax],[0 1],'k')
    for k = 1:length(kc)
        plot(kc(k),snr(k),'bx','linewidth',lnwd)
    end
    hold off
    title(['Dcor analysis : res ~ ',num2str(kcMax,4),', SNR ~ ',num2str(A0,4)])
    xlim([0 1])
    ylim([0 1])
    xlabel('Normalized spatial frequencies')
    ylabel('C.c. coefficients')
end