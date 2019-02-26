function [im,path] = loadData(path)
if nargin < 1
    [fname,pname] = uigetfile('*.*');
    path = [pname,filesep,fname];
end
    

keepReading = 1; k = 1;
warning('off')
while keepReading 
    try
        im(:,:,k) = imread(path,k);
        k = k+1;
    catch
        keepReading = 0;
        disp('Finished reading... ')
        disp(['Stack size : ',num2str(size(im))])
    end
end
warning('on')