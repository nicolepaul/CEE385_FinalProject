function [lambda_edp_c, lambda_edp_nc, midedp_im] = calc_edpHazard(imval, stripeDat, nEDP, pcoeff, typeEDP)

% Discretization
% nIM = 100;
% imPoint = linspace(1e-6,max(imval)*2, nIM)';
n = 1000;
pval = linspace(1e-8, 0.99, n)';

% Calcuating collapse
% Determining probability of collapse fragility
fittype = 'MLE';
ngms = size(stripeDat{1},2);
[IM, ~, ~, fittedFC, ~, ~] = get_collapseRisk(pcoeff, imval, stripeDat, fittype, ngms);

nIM = numel(IM);
imPoint = IM';
dIM = IM(2) - IM(1);

% Getting interpolated or extrapolated values
[thetaPoint, betaPoint] = calc_edpParam(stripeDat, imval, imPoint, nEDP);

% Initialization
edp_im    = NaN(n  , nEDP);
pIMdiff   = NaN(nIM, 1);
midedp_im = NaN(n-1, nEDP);
binedp_im = NaN(n-1, nEDP);
pEDP_im = NaN(n-1, nEDP, nIM);
lambda_edp_im = NaN(n-1, nEDP, nIM);
lambda_edp_im_c = NaN(n-1, nEDP, nIM);
lambda_edp_nc    = NaN(n-1, nEDP);
lambda_edp_c    = NaN(n-1, nEDP);

% Calculating deaggregation curve for EDP
uniqEDP = unique(typeEDP);
for i = 1:numel(uniqEDP)
    indEDP = find(typeEDP == uniqEDP(i));
    edp_im(:,indEDP) = repmat(logninv(pval, log(max(thetaPoint(indEDP,round(nIM/2)),[],1)), max(betaPoint(indEDP,round(nIM/2)),[],1)),1,numel(indEDP));
    binedp_im(:,indEDP) = edp_im(2:end,indEDP) - edp_im(1:end-1,indEDP);
    midedp_im(:,indEDP) = edp_im(1:end-1,indEDP) + binedp_im(:,indEDP)/2;
end


for i = 1:nEDP
    for j = 1:nIM
        % Finding probability exceedance of that edp value
        pEDP_im(:,i,j) = 1 - logncdf(midedp_im(:,i), log(thetaPoint(i,j)), betaPoint(i,j));
    end
end

for i = 1:nIM
    % Interpolating collapse fragility curve
    % Seismic hazard curve differentiation
    lnx = log(imPoint(i));
    p = fliplr(pcoeff);
    pIMdiff(i) = abs(((p(2)+2*p(3)*lnx+3*p(4)*lnx^2+4*p(5)*lnx^3)/imPoint(i))*exp(p(1)+p(2)*lnx+p(3)*lnx^2+p(4)*lnx^3+p(5)*lnx^4));
    % EDP hazard curve
    for j = 1:nEDP
        lambda_edp_im(:,j,i)   = pEDP_im(:,j,i)*pIMdiff(i)*dIM;
        lambda_edp_im_c(:,j,i) = (pEDP_im(:,j,i).*(1 - fittedFC(i)) + fittedFC(i))*pIMdiff(i)*dIM;
    end
end

% Summing for all edp values
for i = 1:nEDP
    lambda_edp_nc(:,i) = sum(squeeze(lambda_edp_im(:,i,:)),2);
    lambda_edp_c(:,i)  = sum(squeeze(lambda_edp_im_c(:,i,:)),2);
end

end
