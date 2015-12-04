function [EL, P, bldgEL, floorEL, compEL, im_ofint, IM, AAL] = calc_lossCurve(frag, imval, stripeDat, nEDP, nfloors, pcoeff, EC, cqty, dqty, dfrag)

%% Collapse Case

% Getting collapse fragility
ngms = size(stripeDat{1},2);
[im_ofint, ~, ~, fittedFC, ~, ~] = get_collapseRisk(pcoeff, imval, stripeDat, fittype, ngms);
n = numel(im_ofint);

% Calculating expected loss for each IM
P.C = fittedFC;
EL.C = P.C*EC.C;
P.NC = 1 - P.C;

%% Demolition Case

% Initialization
dbldgEL = NaN(n, 1);
dfloorEL = NaN(nEDP,n);
pd = NaN(n,1);

% Deaggregation curve of p(D)
pIMdiff = NaN(n,1);
for i = 1:n-1
    lnx = log(im_ofint(i));
    p = fliplr(pcoeff);
    pIMdiff(i) = abs(((p(2)+2*p(3)*lnx+3*p(4)*lnx^2+4*p(5)*lnx^3)/im_ofint(i))*exp(p(1)+p(2)*lnx+p(3)*lnx^2+p(4)*lnx^3+p(5)*lnx^4));
end

% Getting median and dispersion values of EDPs
[thetaPoint, betaPoint] = calc_edpParam(stripeDat, imval, im_ofint, nEDP);


% Looping over all im requested
for i = 1:n
    edp = [thetaPoint(:,i)'; betaPoint(:,i)'];
    [bldgEL_im, floorEL_im, ~, ~, ~, p_fordemo] = calc_lossIntensity(edp, dfrag, dqty, n);
    dbldgEL(i) = bldgEL_im;
    dfloorEL(:,i) = floorEL_im;
    % Calculating total probability of demolition (given no collapse)
    pd(i) = nanmax(p_fordemo,[],2); % Assume max pd across all floors is p(D|IM,NC)
end

% Calculating expected loss for each IM given NC for R
P.D = (pd).*(1-P.C);
EL.D = P.D*EC.D;

%% Non-collapse Case

% Initialization
bldgEL = NaN(n, 1);
floorEL = NaN(nEDP,n);
compEL = NaN(n-1,nEDP,numel(frag),n);

% Looping over all im requested
for i = 1:n
    edp = [thetaPoint(:,i)'; betaPoint(:,i)'];
    [bldgEL_im, floorEL_im, compEL_im] = calc_lossIntensity(edp, frag, cqty, n);
    bldgEL(i) = bldgEL_im;
    floorEL(:,i) = floorEL_im;
    compEL(:,:,:,i) = compEL_im;
end

% Calculating expected loss for each IM given NC for R
EL.R = bldgEL;
P.R = P.NC - P.D;


%% Combination

EL.T = EL.D + EL.C + EL.R;
P.T = P.D + P.C + P.R;
%%
figure;
subplot(2,2,1);
plot(im_ofint, [EL.R EL.D EL.C EL.T], 'LineWidth', 1.2);
grid on;
legend('R','D','C', 'T','Location','best');
xlim([0 max(im_ofint)]);
ylabel('E[L_i|IM]'); xlabel('IM');
title('Expected Loss Due to Each Case');

subplot(2,2,2);
plot(im_ofint, [P.R P.D P.C P.T], 'LineWidth', 1.2);
grid on;
legend('R','D','C', 'T','Location','best');
xlim([0 max(im_ofint)]);
ylabel('p[i|IM]'); xlabel('IM');
title('Probability of Each Case');

subplot(2,2,3);
plot(im_ofint, [EL.R EL.D EL.C  EL.T]./(repmat(EL.T,1,4)), 'LineWidth', 1.2);
grid on;
legend('R','D','C','T','Location','best');
xlim([0 max(im_ofint)]);
ylabel('E[L_i|IM]/E[L_T|IM]'); xlabel('IM');
title('Contribution of Each Case');

subplot(2,2,4);
plot(im_ofint, [EL.R EL.D EL.C  EL.T]./EC.RCN, 'LineWidth', 1.2);
grid on;
legend('R','D','C','T','Location','best');
xlim([0 max(im_ofint)]);
ylabel('E[L_i|IM]/RCN'); xlabel('IM');
title('Expected Loss as Percentage of RCN');

figure;
plot(im_ofint, [EL.R EL.D EL.C  EL.T].*repmat(pIMdiff,1,4), 'Linewidth', 1.2); grid on;
legend('R','D','C','T','Location','best');
xlabel('IM'); ylabel('Deaggregation of Expected Loss');


%%
% Display AAL results

xval = im_ofint;
yval = [EL.R EL.D EL.C  EL.T].*repmat(pIMdiff,1,4);
for i = 1
    indsnan = isnan(xval(:,i)) | isnan(yval(:,i));
    xval(indsnan,:) = [];
    yval(indsnan,:) = [];
end


disp('AAL for each Case (R, D, C, T)');
format bank;
AAL.val = trapz(xval, yval);
disp(AAL.val);
disp('AAL for each Case % of Total');
format short;
AAL.perc = trapz(xval,yval)./trapz(xval(:,end),yval(:,end));
disp(AAL.perc);

IM = im_ofint;
end