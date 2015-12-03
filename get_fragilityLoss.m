function [compnames, compedp, fragDat] = get_fragilityLoss(nfragility, nds, fragfile)

% Import fragility csv
nparam = 2; % central value and uncertainty value
ninfo = 2; % damage and loss
textdat = textscan(fopen(fragfile), '%s', 'Delimiter', ',');
numdat = dlmread(fragfile, ',', [1 2 nfragility 1+nds*nparam*ninfo]);
fclose('all');

% Looping through fragility info
ncol = nds*nparam*ninfo + 2;
fragDat = cell(nfragility,1);
compnames = cell(nfragility,1);
compedp = cell(nfragility,1);
inds_theta = 1:nparam:nds*nparam;
inds_beta = 2:nparam:nds*nparam;
inds_ctheta = nds*nparam+1:nparam:nds*nparam+nds*nparam;
inds_cbeta = nds*nparam+2:nparam:nds*nparam+nds*nparam;
m = 1000;
p = linspace(0,1,m);
for i = 1:nfragility
    % Text data
    compnames{i} = textdat{1}{ncol*i+1};
    compedp{i}   = textdat{1}{ncol*i+2};
    % Numeric data
    fragDat{i}.theta  = numdat(i, inds_theta);
    fragDat{i}.beta   = numdat(i, inds_beta);
    fragDat{i}.ctheta = numdat(i, inds_ctheta);
    fragDat{i}.cbeta  = numdat(i, inds_cbeta);
    fragDat{i}.nds = nds;
    % Solving at multiple points for fragility curve plottings
    fragDat{i}.x = NaN(m, nds);
    fragDat{i}.c = NaN(m, nds);
    fragDat{i}.p = p;
    for j = 1:nds
        fragDat{i}.x(:,j) = logninv(p, log(fragDat{i}.theta(j)), fragDat{i}.beta(j));
        fragDat{i}.c(:,j) = logninv(p, log(fragDat{i}.ctheta(j)), fragDat{i}.cbeta(j));
    end
end

end