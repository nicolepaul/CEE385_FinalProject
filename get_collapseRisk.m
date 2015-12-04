%% Main function
function [IM, pCIM, deagg, fittedFC, MAF, pc50] = get_collapseRisk(pcoeff, imval, stripeDat, fittype, nGM)

% Calculation will be done with below range of IM values
m = 250;
im = linspace(1e-5, max(imval)*1.75, m)'; % Full range
IM = (im(2:end)+im(1:end-1))/2; % At midpoints

% Determining probability of collapse given IM, p(C|IM)
% Note: Defining "NaN" as collapse
nIM = numel(imval);
pCIM = NaN(size(imval));
for i = 1:nIM
    pCIM(i) = sum(isnan(stripeDat{i}(1,:)))/nGM;
end

% Performing calculations using fit type requested
fc = get_collapseFit(fittype, imval, pCIM, nGM);

% Using fragility curve function on range of IM values
fittedFC = fc(IM);

% Performing numerical integration using closed form solution of seismic 
% hazard curve
pcim = NaN(m-1,1);
pIMdiff = NaN(m-1,1);
deagg = NaN(m-1,1);
for i = 1:m-1
    pcim(i) = fc(IM(i));
    lnx = log(IM(i));
    p = fliplr(pcoeff);
    pIMdiff(i) = abs(((p(2)+2*p(3)*lnx+3*p(4)*lnx^2+4*p(5)*lnx^3)/IM(i))*exp(p(1)+p(2)*lnx+p(3)*lnx^2+p(4)*lnx^3+p(5)*lnx^4));
    deagg(i) = pcim(i)*pIMdiff(i);
end

% Mean annual frequency of collapse
MAF = trapz(IM, deagg);
% Probability of collapse in 50 years
pc50 = 1 - exp(-MAF*50);

end

%% Subfunction
function [fc, FC] = get_collapseFit(fittype, imval, pCIM, nGM)

% Fitting MLE parameters using lognormal distribution
FC.mle = lognfit(imval, 0.5, [], pCIM*nGM);
% Fitting least square paramterers
FC.lsq = lsqnonneg(imval, pCIM);

% Writing fragility curve equation based off of fit type
if strcmpi(fittype, 'Least Squares')
    fc = @(x) FC.lsq.*x;
else
    fc = @(x) logncdf(x, FC.mle(1), FC.mle(2));
end

end