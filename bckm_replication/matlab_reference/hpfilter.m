function [trend, cycle] = hpfilter(y, lambda)
% HPFILTER  Hodrick-Prescott filter (inline implementation).
%
%   [TREND, CYCLE] = HPFILTER(Y, LAMBDA)
%
%   Matches the call signature of MATLAB's Econometrics Toolbox hpfilter.
%   Y      : T x 1 time series
%   LAMBDA : smoothing parameter (1600 for quarterly data)
%   TREND  : HP trend component
%   CYCLE  : HP cyclical component  (= Y - TREND)
%
%   Solution: trend = (I + lambda * D2' * D2) \ y
%   where D2 is the (T-2) x T second-difference matrix.

T = length(y);
if T < 3
    trend = y;
    cycle = zeros(size(y));
    return
end

e  = ones(T, 1);
D2 = spdiags([e, -2*e, e], 0:2, T-2, T);
A  = speye(T) + lambda * (D2' * D2);
trend = A \ y;
cycle = y - trend;

end
