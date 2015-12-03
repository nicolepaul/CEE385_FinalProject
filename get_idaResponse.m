function [imval, stripeDat, imchoice, unitstr, edptext] = get_idaResponse(dirpath, nstripe, nedps, ngms)

% Opening directory with IDA results
FileList = dir(dirpath);
FileList = FileList(arrayfun(@(x) x.name(1), FileList) ~= '.');
N = size(FileList, 1);
% Display warning if nstripe does not match number of stripes found
if nstripe ~= N
    warning('Number of stripes input does not match number of stripes found');
end

% Initialization of variables
imchoice = cell(nstripe, 1);
unitstr = cell(nstripe, 1);
imval = NaN(nstripe, 1);
edptext = cell(nstripe, 1);
stripeDat = cell(nstripe, 1);

% Reading each data file
for i = 1:nstripe
    
    % Reading IM case
    fid = fullfile(dirpath,FileList(i).name);
    imstr = textscan(fopen(fid), '%s', 'Delimiter', ',');
    [tok, ~] = strsplit(imstr{1}{1}, '=');
    imchoice{i} = tok{1};
    unitstr{i} = tok{2}(end); % Assuming last character is units
    imval(i) = str2double(tok{2}(1:end-1));
    
    % Reading each EDP case
    headers = transpose(reshape(imstr{1}, ngms+1, nedps+1));
    edptext{i} = headers(2:end, 1);
    
    % Reading data
    stripeDat{i} = dlmread(fid, ',', 1, 1);
    
end
fclose('all');

end