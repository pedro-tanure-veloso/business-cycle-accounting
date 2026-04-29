function R = xcorr(x, y, maxlag, opt)
% XCORR  Normalized cross-correlation at lags -maxlag .. +maxlag.
%
%   R = XCORR(X, Y, MAXLAG, 'Coef')
%
%   Matches the signature used in gwedges2.m.
%   'Coef' divides by sqrt(Rxx(0)*Ryy(0)), i.e. normalizes to [-1,1].
%   Result is a (2*maxlag+1) x 1 column vector, lag order -maxlag..+maxlag.

if nargin < 3; maxlag = length(x) - 1; end
if nargin < 4; opt = 'none'; end

x = x(:);
y = y(:);
T = length(x);
R = zeros(2*maxlag+1, 1);

for k = -maxlag:maxlag
    idx = k + maxlag + 1;
    if k >= 0
        n = T - k;
        R(idx) = sum(x(1:n) .* y(k+1:k+n));
    else
        n = T + k;
        R(idx) = sum(x(-k+1:-k+n) .* y(1:n));
    end
end

if strncmpi(opt, 'coef', 4)
    Rxx0 = sum(x .* x);
    Ryy0 = sum(y .* y);
    if Rxx0 > 0 && Ryy0 > 0
        R = R / sqrt(Rxx0 * Ryy0);
    end
end

end
