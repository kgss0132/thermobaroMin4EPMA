function [mask, info] = selectDataset(T)
% liquid.selectDataset
% Build a unique list of datasets and let user pick one.
% Default grouping is by "Experiment" (can be changed below).

% --- choose grouping key:
% Option A (recommended): group by Experiment (often stable across many rows)
groupKey = "Experiment";
% Option B: group by Index (uncomment if Index uniquely identifies dataset)
% groupKey = "Index";

idxCol = T.('Index');
expCol = T.('Experiment');
citCol = T.('Citation');

% Convert to string for robust display/compare
idxS = local_toString(idxCol);
expS = local_toString(expCol);
citS = local_toString(citCol);

% Build unique groups by chosen key
switch groupKey
    case "Experiment"
        keyS = expS;
    case "Index"
        keyS = idxS;
    otherwise
        error('Unknown groupKey: %s', groupKey);
end

[uniqKey, ~, g] = unique(keyS, 'stable');

% For each group, pick representative Index/Experiment/Citation (first row)
repIdx = strings(numel(uniqKey),1);
repExp = strings(numel(uniqKey),1);
repCit = strings(numel(uniqKey),1);
nRows  = zeros(numel(uniqKey),1);

for i = 1:numel(uniqKey)
    rows = (g == i);
    first = find(rows, 1, 'first');
    repIdx(i) = idxS(first);
    repExp(i) = expS(first);
    repCit(i) = citS(first);
    nRows(i) = sum(rows);
end

% Build display list
items = strings(numel(uniqKey),1);
for i = 1:numel(uniqKey)
    items(i) = sprintf('%s | %s | %s  (n=%d)', repIdx(i), repExp(i), repCit(i), nRows(i));
end

% UI select
[sel, ok] = listdlg( ...
    'PromptString', 'Select one dataset (Index | Experiment | Citation):', ...
    'SelectionMode', 'single', ...
    'ListString', cellstr(items), ...
    'ListSize', [900 420]);

if ~ok || isempty(sel)
    error('User cancelled dataset selection.');
end

% mask rows that belong to selected group
mask = (g == sel);

info = struct();
info.groupKey = char(groupKey);
info.keyValue = char(uniqKey(sel));
info.representativeIndex = char(repIdx(sel));
info.representativeExperiment = char(repExp(sel));
info.representativeCitation = char(repCit(sel));
info.nRows = nRows(sel);

end

function s = local_toString(x)
if isstring(x)
    s = x;
elseif ischar(x)
    s = string(cellstr(x));
elseif iscellstr(x)
    s = string(x);
elseif iscategorical(x)
    s = string(x);
elseif isnumeric(x)
    s = string(x);
else
    % fallback: stringify each element
    s = strings(numel(x),1);
    for k = 1:numel(x)
        s(k) = string(x(k));
    end
end
s = s(:);
end
