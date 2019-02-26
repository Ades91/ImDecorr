function cc = getCorrcoef(I1,I2,c1,c2)

if nargin < 4
    c2 = sqrt(sum(sum(abs(I2).^2)));
end
if nargin < 3
	c1 = sqrt(sum(sum(abs(I1).^2)));
end

% N = (numel(I1)-1);
cc = sum(sum(real(I1.*conj(I2))))./((c1*c2));