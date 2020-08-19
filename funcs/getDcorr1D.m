function [kcMax,A0,d0,d] = getDcorr1D(sig,r,Ng,figID)

if nargin < 4; figID = 0; end
if ischar(figID)
    figID = 0;
    fastMode = 1;
else
    fastMode = 0;
end
if nargin < 3; Ng = 10;  end
if nargin < 2; r = linspace(0,1,50); end

% input check
if length(r) < 30
    r = linspace(min(r),max(r),min(30,(numel(sig)+1)/2));
elseif length(r) > (numel(sig)+1)/2
    r = linspace(min(r),max(r),(numel(sig)+1)/2);
end
if Ng < 5
    Ng = 5;
end
if size(sig,1) > 1 && size(sig,2) == 1
    sig = sig';
end
sig = single(sig);
sig = sig(1:end-not(mod(numel(sig),2))); % odd number of pixels

R = abs(linspace(-1,1,numel(sig)));
Nr = length(r);
if isa(sig,'gpuArray')
    r = gpuArray(r); 
    R = gpuArray(R);
end
% Sn : Fourier normalized signal
Sn = fftshift(fft(fftshift(sig))); Sn = Sn./abs(Sn); Sn(isinf(Sn)) = 0; Sn(isnan(Sn)) = 0;
mask0 = R <= 1;
Sn = mask0.*Sn; % restric all the analysis to the region r < 1

if figID
    fprintf('Computing dcorr: ');
end

Sk = mask0.*fftshift(fft(fftshift(sig)));

c = sqrt(sum(sum(Sk.*conj(Sk))));

r0 = linspace(r(1),r(end),Nr);

for k = length(r0):-1:1
    cc = getCorrcoef(Sk,(R < r0(k)).*Sn,c);
    if isnan(cc); cc = 0; end
    d0(k) = gather(cc); % gather if input image is gpuArray 

    if fastMode == 1
        [ind0,snr0] = getDcorrLocalMax(d0(k:end));
        ind0 = ind0 + k-1;
        if ind0 > k % found a local maxima, skip the calculation
            break;
        end
    end
end
if fastMode == 0
%     d0(k:end) = d0(k:end).*sqrt(r0(k:end));
    ind0 = getDcorrLocalMax(d0(k:end));
    snr0 = d0(ind0);
end
k0 = gather(r(ind0));

gMax = 2/r0(ind0);
if isinf(gMax); gMax = max(numel(sig))/2;end

% search of highest frequency peak
g = [numel(sig/4)/4, exp(linspace(log(gMax),log(0.15),Ng))];
d = zeros(Nr,2*Ng); kc = k0; SNR = snr0;
if fastMode == 0
    ind0 = 1;
else
    if ind0 > 1
        ind0 = ind0 -1;
    end
end
for refin = 1:2 % two step refinement

    for h = 1:length(g)
        Sr = Sk.*(1 - exp(-2*g(h)*g(h)*R.^2)); % Fourier Gaussian filtering
        
        c = sqrt(sum(sum(abs(Sr).^2)));
        
        for k = length(r):-1:ind0

            if isa(sig,'gpuArray')
                cc = getCorrcoef(Sr,Sn.*(R < r(k)),c);
                if isnan(cc); cc = 0; end
                d(k,h + Ng*(refin-1)) = gather(cc);
            else
                mask = (R < r(k));
                cc = getCorrcoef(Sr(mask),Sn(mask),c);
                if isnan(cc); cc = 0; end
                d(k,h + Ng*(refin-1)) = cc;
            end
            if fastMode
                [ind,snr] = getDcorrLocalMax(d(k:end,h + Ng*(refin-1)));
                ind = ind +k-1;
                if ind > k % found a local maxima, skip the calculation
                    break;
                end
            end

        end
        if fastMode == 0
%             d(k:end,h + Ng*(refin-1)) = d(k:end,h + Ng*(refin-1)).*sqrt(r(k:end)');
            ind = getDcorrLocalMax(d(k:end,h + Ng*(refin-1)));
            snr = d(ind,h + Ng*(refin-1));
            ind = ind +k-1;
        end
        kc(h + Ng*(refin-1)+1) = gather(r(ind));
        SNR(h + Ng*(refin-1)+1) = snr;
        if figID
        	fprintf('-');
        end
    end

% refining the high-pass threshold and the radius sampling
    if refin == 1

    % high-pass filtering refinement
        indmax = find(kc == max(kc));
        ind1 = indmax(end);
        if ind1 == 1 % peak only without highpass 
            ind1 = 1;
            ind2 = 2;
            g1 = size(sig,1);
            g2 = g(1);
        elseif ind1 >= numel(g)
            ind2 = ind1-1;
            ind1 = ind1-2;
            g1 = g(ind1); g2 = g(ind2);
        else
            ind2 = ind1;
            ind1 = ind1-1;
            g1 = g(ind1); g2 = g(ind2);
        end
        g = exp(linspace(log(g1),log(g2),Ng));
        
        % radius sampling refinement

        r1 = kc(indmax(end))-(r(2)-r(1)); r2 = kc(indmax(end))+0.3;
        if r1 < 0 ; r1 = 0; end
        if r2 > 1; r2 = 1; end
        r = linspace(r1,min(r2,r(end)),Nr);
        ind0 = 1;
        r2 = r;
    end
end
if figID
    fprintf(' -- Computation done -- \n');
end

% add d0 results to the analysis (usefull for high noise images)
kc(end+1) = gather(k0);
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
    A0 = snr0; % average image contrast has to be estimated from original image
else
    kcMax = r(2);
    Amax = 0;
    res = r(2);
    A0 = 0;
end

% results display if figID specified
if figID 
    radAv = log(abs(gather(Sk((end+1)/2:end)))+1);
    lnwd = 1.5;
    figure(figID);
    plot(r0,d(:,1:Ng),'color',[0.2 0.2 0.2 0.5]);
    hold on

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