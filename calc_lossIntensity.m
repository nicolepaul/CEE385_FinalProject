function [bldgEL_im, floorEL_im, compqEL_im, pDS_edp_im, midedp_im, p_fordemo] = calc_lossIntensity(edp, frag, cqty, n)
% GET_EXPLOSSNC
% Calculates expected loss in building for a given IM level
%
% Inputs:
%  -  edp: theta and beta values for each floor given IM of interest (logn)
%  -  frag: component fragility curves
%  -  cqty: quantity of each component on each floor
%
% Outputs:
%  -  bldgEL_im: total building expected loss given non-collapse, given IM
%  -  floorEL_im: expected loss for each floor given non-collapse, given IM
%  -  comp_EL_im: expected loss for each component given non-collapse, given IM
%  -  pDS_edp_im: probably of each DS given EDP, given IM, given  non-collapse
%  -  midedp_im: midpoint of each EDP bin for each story

% Counting inputs
nstory = size(edp, 2);
ncomp = numel(frag);

% Choosing number of points for EDP discretization
p = linspace(1e-5,0.99999,n)';

% Initialization
edp_im = NaN(n, nstory);
midedp_im = NaN(n-1, nstory);
binedp_im = NaN(n-1, nstory);
pDS_edp_im = cell(ncomp,1);
pEDP_im = NaN(n-1,nstory);
cEDP_im = NaN(n-1,nstory);
qtyEL_im = NaN(n-1,nstory);
cEL_im = NaN(n-1,nstory);
compEL_im = NaN(n-1, nstory, ncomp);
compqEL_im = NaN(nstory, ncomp);
floorEL_im = NaN(nstory,1);
p_fordemo = NaN(1, nstory);


% Calculating loss for each component
for i = 1:nstory
    % Finding range of edp's of interest
    edp_im(:,i) = logninv(p, log(edp(1,i)), edp(2,i));
    binedp_im(:,i) = edp_im(2:end,i) - edp_im(1:end-1,i);
    midedp_im(:,i) = edp_im(1:end-1,i) + binedp_im(:,i);
    % Finding probability of that edp value for that story
    pEDP_im(:,i) = lognpdf(midedp_im(:,i), log(edp(1,i)), edp(2,i)).*binedp_im(:,i);
    cEDP_im(:,i) = logncdf(midedp_im(:,i), log(edp(1,i)), edp(2,i)).*binedp_im(:,i);
    % Determining losses for each component
    for j = 1:ncomp
        pDS_edp_im{j} = NaN(n-1, frag{j}.nds, nstory);
        % Determining probability of each damage state
        for k = frag{j}.nds:-1:1
            if k == frag{j}.nds
                pDS_edp_im{j}(:,k,i) = logncdf(midedp_im(:,i), log(frag{j}.theta(k)), frag{j}.beta(k));
            else
                pDS_edp_im{j}(:,k,i) = logncdf(midedp_im(:,i), log(frag{j}.theta(k)), frag{j}.beta(k)) - pDS_edp_im{j}(:,k+1,i);
            end
        end
        % Multiplying probability of each damage state by expected loss
        compEL_im(:,i,j) = sum(repmat(frag{j}.ctheta,n-1,1).*pDS_edp_im{j}(:,:,i),2);
        compqEL_im(i,j) = sum(compEL_im(:,i,j)*cqty(i,j),1);
    end
    % Multiplying component loss by quantity at that floor
    qtyEL_im(:,i) = sum(squeeze(compEL_im(:,i,:)).*repmat(cqty(i,:),n-1,1),2);
    % Multiplying loss by probability of seeing that edp
    cEL_im(:,i) = qtyEL_im(:,i).*pEDP_im(:,i);
    % Calculating probability for demolition
    p_fordemo(i) = sum(repmat(cqty(i,1),n-1,1).*pDS_edp_im{1}(:,1,i).*pEDP_im(:,i),1);
    % Summing expected losses for all floors
    floorEL_im(i) = sum(cEL_im(:,i));
end

% Summing expected losses across entire building
bldgEL_im = sum(floorEL_im);


end