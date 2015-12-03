function [thetaPoint, betaPoint] = calc_edpParam(stripeDat, imval, imPoint, nEDP)
% CALC_EDPPARAM
% Calculates the geomean and dispersion values (lognormal distribution)
% given the results from a stripe analysis
% 
% Outputs:
% thetaPoint   -  geomean value at IM of interest
% betaPoint    -  dispersion value at IM of interest

% Processing inputs
nstripe = numel(stripeDat);

% Finding the geomean and dispersion of each EDP value, for each stripe
thetaStripe = NaN(nEDP, nstripe);
betaStripe  = NaN(nEDP, nstripe);
for i = 1:nstripe
    thetaStripe(:, i) = calc_geoMean(stripeDat{i},2);
    betaStripe (:, i) = nanstd(log(stripeDat{i}),[],2) ;
end

% Figuring out whether needs to interpolate or extrapolate
interpPoint = imPoint <= max(imval);
extrapPoint = imPoint >  max(imval);

% Initialization
thetaPoint = NaN(nEDP, numel(imPoint));
betaPoint  = NaN(nEDP, numel(imPoint));

% Setting all extrapolated points to be a plateau at the maximum stripe
% result and all 
thetaPoint(:, extrapPoint) = repmat(max(thetaStripe, [], 2),1,sum(extrapPoint));
betaPoint (:, extrapPoint) = repmat(max(betaStripe , [], 2),1,sum(extrapPoint));

% Setting all interpolated points to be linear between given stripe values
for j = 1:nEDP
    thetaPoint(j, interpPoint) = interp1([0; imval], [0; thetaStripe(j,:)'], imPoint(interpPoint))';
    betaPoint (j, interpPoint) = interp1([0; imval], [0;  betaStripe(j,:)'], imPoint(interpPoint))';
end