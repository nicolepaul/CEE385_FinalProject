function [IM, EDP, pEDPIM, deagg, fittedFC, MAF, pc50] = get_edpHazard(pcoeff, imval, edpDat, nGM)

% Calculation will be done with below range of IM values
m = 500;
im = linspace(0, 5, m)'; % Full range
IM = (im(2:end)+im(1:end-1))/2; % At midpoints

% Determining probability of collapse given IM, p(C|IM)
% Note: Defining "NaN" as collapse
nIM = numel(imval);
pCIM = NaN(size(imval));
for i = 1:nIM
    pCIM(i) = sum(isnan(edpDat{i}(1,:)))/nGM;
end

% Determining probability of exceeding certain EDP given IM % NOTE: Check
% this calc
nvals = 500;
EDP = linspace(0, max(edpDat{end}(:))*0.75, nvals);
pEDPIM = NaN(nvals, nIM);
for i = 1:nIM
    for j = 1:nvals
        pEDPIM(j,i) = pCIM(i) + sum(max(abs(edpDat{i}),[],1)>EDP(j))/nGM;
    end
end

% Fitting MLE parameters using lognormal distribution
fittedFC = NaN(m-1, nvals);
for j = 1:nvals
    mleparam = lognfit(imval', 0.5, [], pEDPIM(j,:)*nGM');
    fc = @(x) logncdf(x, mleparam(1), mleparam(2));
    fittedFC(:,j) = fc(IM);
end

% Performing numerical integration using closed form solution of seismic 
% hazard curve
pIMdiff = NaN(m-1, nvals);
deagg = NaN(m-1, nvals);
MAF = NaN(nvals,1);
pc50 = NaN(nvals,1);
for j = 1:nvals
    for i = 1:m-1
        lnx = log(IM(i));
        p = fliplr(pcoeff);
        pIMdiff(i) = abs(((p(2)+2*p(3)*lnx+3*p(4)*lnx^2+4*p(5)*lnx^3)/IM(i))*exp(p(1)+p(2)*lnx+p(3)*lnx^2+p(4)*lnx^3+p(5)*lnx^4));
        deagg(i,j) = fittedFC(i,j)*pIMdiff(i);
    end
    MAF(j) = trapz(IM, deagg(:,j));
    pc50(j) = 1 - exp(-MAF(j)*50);
end


end