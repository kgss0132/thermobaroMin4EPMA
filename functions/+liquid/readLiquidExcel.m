function [liq, meta] = readLiquidExcel()
% liquid.readLiquidExcel
% Read ExPetDB-style liquid composition Excel and let user select ONE dataset
% via (Index | Experiment | Citation) list (using liquid.selectDataset).
%
% NOTE
%   Column names may be 'SiO2value' etc. No required-column check here.

% default folder: functions/+liquid/+liquidComp
thisFileDir = fileparts(mfilename('fullpath'));     % .../functions/+liquid
defaultDir  = fullfile(thisFileDir, '+liquidComp'); % .../functions/+liquid/+liquidComp
if ~exist(defaultDir, 'dir')
    defaultDir = thisFileDir;
end

[file, path] = uigetfile({'*.xlsx;*.xls','Excel Files (*.xlsx, *.xls)'}, ...
    'Select liquid composition Excel (ExPetDB format)', defaultDir);

if isequal(file, 0)
    error('User cancelled file selection.');
end

fullpath = fullfile(path, file);

opts = detectImportOptions(fullpath, 'NumHeaderLines', 0);
T = readtable(fullpath, opts);

% Preserve original names across MATLAB versions
T.Properties.VariableNames = opts.VariableNames;

% Select dataset (requires functions/+liquid/selectDataset.m)
[selMask, selInfo] = liquid.selectDataset(T);

liq = T(selMask, :);

meta = struct();
meta.file = fullpath;
meta.selected = selInfo;

end
