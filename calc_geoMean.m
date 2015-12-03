function M = calc_geoMean(X, dim)
% CALC_GEOMEAN
% This function will calculate the geomean (MLE estimator of mean for
% lognormal distribution) ignoring all NaN values.  Note that only up to 2
% dimensions are supported.
% 
% Inputs:
% X   -  Matrix or vector of numbers to take geomean of
% DIM -  Dimension along which to take geomean (default = 1)
% 
% Outputs:
% M   - Output matrix or vector of geomean values

% Default dimension of 1
if nargin == 1
    dim = 1;
end

% Find all lines which have NaN and removing them
if dim == 1
    realX = X(~any(isnan(X),2),:);
else
    realX = X(:,~any(isnan(X),1));
end

% Calculating geomean of result above
M = geomean(realX, dim);

end
