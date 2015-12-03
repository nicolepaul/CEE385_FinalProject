%% Main function
function shcdat = get_usgsHazard(lat, long, siteclass, period)

% Processing inputs
[site, per] = get_usgsFormat(siteclass, period);


% Forming url for data access
str1 = 'http://geohazards.usgs.gov/hazardtool/curves.php?format=2&lat=';
str2 = '&lon=';
str3 = '&site=';
str4 = '&period=';
url = strcat(str1,lat,str2,long,str3,site,str4,per);

% Saving data locally
filename = strcat('LAT',lat,'LON',long,'SC',site,'T',per,'.csv');
shcFile = fullfile('SHC',filename);
urlwrite(url,shcFile);

% Importing data
try dlmread(shcFile, ',', 1, 5);
    shcdat = dlmread(shcFile, ',', 1, 5)';
catch
    shcdat = dlmread(shcFile, ',', 1, 6)';
end

end

%% Subfunction
function [site, per] = get_usgsFormat(siteclass, period)

% Switching site type to match usgs expected string
switch siteclass
    case 1
        site = '2000';
    case 2
        site = '1150';
    case 3
        site = '760';
    case 4
        site = '537';
    case 5
        site = '360';
    case 6
        site = '259';
    case 7
        site = '180';
end

% Switching period to match usgs expected string
switch period
    case 1
        per = '0p00';
    case 2
        per = '0p10';
    case 3
        per = '0p20';
    case 4
        per = '0p30';
    case 5
        per = '0p50';
    case 6
        per = '0p75';
    case 7
        per = '1p00';
    case 8
        per = '2p00';
    case 9
        per = '3p00';
    case 10
        per = '4p00';
    case 11
        per = '5p00';
end

end