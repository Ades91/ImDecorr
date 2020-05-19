addpath('funcs')
%% Data import
% import your localization data here

%% data rendering
pps0 = 5;       % projected pixel size of 10 nanometers
fov = 12000;    % FOV of 12um

% render the localization data matrix as a bilinear histogram
% data should have the size [N,2] and be in the range [0,fov]
im = smlmHist(data,pps0,fov);
figure(100);
imagesc(im);colormap(hot);colorbar
title('Bilinear histogram')

% estimate the object area from the rendering.
% changing pps0 can significantly change the estimated sample area
% for comparison with the manuscript (add ref), set pps0 at 5nm
areanm = sum(im(:)>0.5)*pps0*pps0;
areaum = areanm*1e-6;
%% processing loop
% parameters 
pps = 7:2:35;  % pixel size
Nfr = 15;      % number of frames, if Nfr == 1 all the localizations are used
nloc = round(linspace(50000,nLoc,Nfr));

% compute the localization density normalized to the sample area
density = nloc/areaum;

kc = zeros(Nfr,numel(pps)); res = kc;
for nn = 1:Nfr
    for pp = numel(pps):-1:1
        temp = smlmHist(data(1:nloc(nn),:),pps(pp),fov);
%         figure(101);imagesc(temp)
        % 'fast' option 
        kc(nn,pp) = getDcorr(gpuArray(apodImRect(single(temp),20)),linspace(0,1,50),10,'fast');
        % if no local max found, going to smaller pixel size will not help
        % => break from the loop
        if kc(nn,pp) == 0
            break;
        end
        pause(0.05) % for matlab to update figure
    end
    % convert the cut off frequency in nm
    res(nn,:) = 2*pps./kc(nn,:);
end

%% figures
% resolution as a function of the pixel size
tmp = min(res,[],1);
tmp(tmp > 205) = 205; % clip inf and large res value for plots
figure(102)
plot(pps,tmp,'-x','linewidth',1.5)
hold on;plot(pps,2*pps,'k--','linewidth',1.5); hold off
xlabel('Pixel size (nm)'); ylabel('Resolution (nm)')
% resolution as a function of the localization density
tmp = min(res,[],2);
tmp(tmp > 205) = 205; % clip inf and large res value for plots
figure(103)
plot(density,tmp,'-x','linewidth',1.5)
xlabel('Loc. density (#loc./\mu m^2)'); ylabel('Resolution (nm)')