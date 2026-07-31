function thermobaroMin4EPMA()
%THERMOBAROMIN4EPMA
% Entry point for thermometry/barometry calculations.
%
% Available modes:
%   1) Thermometer with fixed pressure
%   2) Barometer with fixed temperature
%   3) Thermometer over a pressure range
%   4) Barometer over a temperature range
%
% Range-mode result tables must contain:
%   T_degreeC
%   P_kbar
%
% For range-mode plotting, result rows are separated into individual
% calculation vectors using dataCode variables and range-value resets.

%% Ensure the functions folder is on the MATLAB path
projectRoot = fileparts(mfilename('fullpath'));
functionsDir = fullfile(projectRoot, 'functions');

if isfolder(functionsDir)
    addpath(functionsDir);
end

%% Select the calculation mode
choice = menu( ...
    'Select calculation mode:', ...
    'Thermometer: fixed pressure', ...
    'Barometer: fixed temperature', ...
    'Thermometer: pressure range', ...
    'Barometer: temperature range', ...
    'Cancel');

if ~ismember(choice, 1:4)
    return;
end

mode = getModeConfiguration_(choice);

%% Import the cation dataset
[rawdata_struct, info] = importCationDataset_();

%% Check the selected launcher
if exist(mode.launcher, 'file') ~= 2
    error( ...
        'thermobaroMin4EPMA:LauncherNotFound', ...
        ['%s.m was not found on the MATLAB path.\n' ...
         'Expected launcher for this mode: %s.m'], ...
        mode.launcher, ...
        mode.launcher);
end

%% Run the selected workflow
[results, funcName, groupName] = ...
    feval(mode.launcher, rawdata_struct);

%% Validate returned results
if isempty(results) || ~istable(results) || height(results) == 0
    disp('No results were returned because the calculation was canceled or empty.');
    return;
end

%% Construct the condition tag for the output filename
conditionTag = buildConditionTag_(results, mode);

%% Clean function and group names for the output filename
funcNameClean  = cleanFileToken_(funcName);
groupNameClean = cleanFileToken_(groupName);

%% Construct the output filename
fileParts = { ...
    info.fileName, ...
    mode.calcLabel, ...
    groupNameClean, ...
    funcNameClean, ...
    conditionTag};

fileParts = fileParts(~cellfun(@isempty, fileParts));

outputStem = strjoin(fileParts, '_');

outputDir  = fileparts(info.fileFullPath);
outputFile = fullfile(outputDir, [outputStem '.xlsx']);

%% Export the calculation results
writetable(results, outputFile, 'Sheet', 'Results');

disp(['Data have been exported to: ' outputFile]);

%% Draw and export the P-T plot for range modes
if mode.makeTPPlot
    try
        fig = plotTPResults_( ...
            results, ...
            mode, ...
            funcName, ...
            groupName);

        plotFile = fullfile(outputDir, [outputStem '_TPplot.png']);

        try
            exportgraphics(fig, plotFile, 'Resolution', 300);
        catch
            saveas(fig, plotFile);
        end

        disp(['P-T plot has been exported to: ' plotFile]);

    catch ME
        fprintf(2, ...
            ['WARNING: Results were exported, but the P-T plot could not ' ...
             'be created.\n' ...
             '         %s\n'], ...
            ME.message);
    end
end

end


%% ====================================================================== %
function mode = getModeConfiguration_(choice)
%GETMODECONFIGURATION_ Return settings for the selected calculation mode.

mode = struct();

switch choice

    case 1
        mode.launcher          = 'startThermoCalc_fixedP';
        mode.calcLabel         = 'Geotherm_fixedP';
        mode.displayName       = 'Thermometer: fixed pressure';

        mode.conditionVariable = 'P_kbar';
        mode.conditionType     = 'single';
        mode.conditionPrefix   = '';
        mode.conditionUnit     = 'kbar';

        mode.makeTPPlot        = false;
        mode.sortVariable      = '';

    case 2
        mode.launcher          = 'startBaroCalc_fixedT';
        mode.calcLabel         = 'Geobaro_fixedT';
        mode.displayName       = 'Barometer: fixed temperature';

        mode.conditionVariable = 'T_degreeC';
        mode.conditionType     = 'single';
        mode.conditionPrefix   = '';
        mode.conditionUnit     = 'degC';

        mode.makeTPPlot        = false;
        mode.sortVariable      = '';

    case 3
        mode.launcher          = 'startThermoCalc_rangeP';
        mode.calcLabel         = 'Geotherm_rangeP';
        mode.displayName       = 'Thermometer: pressure range';

        mode.conditionVariable = 'P_kbar';
        mode.conditionType     = 'range';
        mode.conditionPrefix   = 'P';
        mode.conditionUnit     = 'kbar';

        mode.makeTPPlot        = true;
        mode.sortVariable      = 'P_kbar';

    case 4
        mode.launcher          = 'startBaroCalc_rangeT';
        mode.calcLabel         = 'Geobaro_rangeT';
        mode.displayName       = 'Barometer: temperature range';

        mode.conditionVariable = 'T_degreeC';
        mode.conditionType     = 'range';
        mode.conditionPrefix   = 'T';
        mode.conditionUnit     = 'degC';

        mode.makeTPPlot        = true;
        mode.sortVariable      = 'T_degreeC';

    otherwise
        error('Invalid calculation mode.');
end

end


%% ====================================================================== %
function conditionTag = buildConditionTag_(results, mode)
%BUILDCONDITIONTAG_ Create a filename tag from a fixed value or value range.

conditionTag = '';
variableName = mode.conditionVariable;

if ~ismember(variableName, results.Properties.VariableNames)
    fprintf(2, ...
        ['WARNING: The results table does not contain "%s". ' ...
         'The condition tag will be omitted from the filename.\n'], ...
        variableName);
    return;
end

values = results.(variableName);

if ~isnumeric(values)
    fprintf(2, ...
        ['WARNING: The results variable "%s" is not numeric. ' ...
         'The condition tag will be omitted from the filename.\n'], ...
        variableName);
    return;
end

values = values(:);
values = values(isfinite(values));

if isempty(values)
    return;
end

switch mode.conditionType

    case 'single'
        valueString = formatNumberForFilename_(values(1));

        conditionTag = [ ...
            mode.conditionPrefix, ...
            valueString, ...
            mode.conditionUnit];

    case 'range'
        minimumValue = min(values);
        maximumValue = max(values);

        minimumString = formatNumberForFilename_(minimumValue);
        maximumString = formatNumberForFilename_(maximumValue);

        if abs(maximumValue - minimumValue) < 1e-10
            conditionTag = [ ...
                mode.conditionPrefix, ...
                minimumString, ...
                mode.conditionUnit];
        else
            conditionTag = [ ...
                mode.conditionPrefix, ...
                minimumString, ...
                'to', ...
                maximumString, ...
                mode.conditionUnit];
        end

    otherwise
        error('Unknown condition type: %s', mode.conditionType);
end

end


%% ====================================================================== %
function numberString = formatNumberForFilename_(value)
%FORMATNUMBERFORFILENAME_ Convert a number into a filename-safe string.

if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    numberString = '';
    return;
end

if abs(value - round(value)) < 1e-8
    numberString = sprintf('%d', round(value));
else
    numberString = sprintf('%.4f', value);
    numberString = regexprep(numberString, '0+$', '');
    numberString = regexprep(numberString, '\.$', '');
    numberString = strrep(numberString, '.', 'p');
end

end


%% ====================================================================== %
function token = cleanFileToken_(inputText)
%CLEANFILETOKEN_ Remove characters unsuitable for an output filename.

if isempty(inputText)
    token = '';
    return;
end

token = char(string(inputText));
token = regexprep(token, '[^a-zA-Z0-9_-]', '');

end


%% ====================================================================== %
function fig = plotTPResults_(results, mode, funcName, groupName)
%PLOTTPRESULTS_ Draw temperature-pressure results.
%
% Horizontal axis:
%   Temperature in degrees Celsius
%
% Vertical axis:
%   Pressure in kbar
%
% Each calculation vector is plotted separately using a different line
% color and marker symbol.

requiredVariables = {'T_degreeC', 'P_kbar'};

missingVariables = requiredVariables( ...
    ~ismember(requiredVariables, results.Properties.VariableNames));

if ~isempty(missingVariables)
    error( ...
        'The results table must contain the following variables: %s', ...
        strjoin(missingVariables, ', '));
end

temperature = results.T_degreeC;
pressure    = results.P_kbar;

if ~isnumeric(temperature) || ~isnumeric(pressure)
    error('T_degreeC and P_kbar must both be numeric variables.');
end

temperature = temperature(:);
pressure    = pressure(:);

if numel(temperature) ~= height(results) || ...
        numel(pressure) ~= height(results)
    error( ...
        'T_degreeC and P_kbar must contain one value per result-table row.');
end

%% Identify individual calculation vectors
[curveIndex, curveLabels] = identifyDatasetCurves_( ...
    results, ...
    mode.sortVariable);

%% Remove rows with invalid temperature or pressure values
validRows = isfinite(temperature) & isfinite(pressure);

if ~any(validRows)
    error('No finite P-T results were available for plotting.');
end

temperature = temperature(validRows);
pressure    = pressure(validRows);
curveIndex  = curveIndex(validRows);

%% Create the figure
fig = figure( ...
    'Name', mode.displayName, ...
    'Color', 'w');

ax = axes(fig);

hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

%% Prepare line colors and marker symbols
curveNumbers = unique(curveIndex, 'stable');
numberOfCurves = numel(curveNumbers);

lineColors = lines(max(numberOfCurves, 1));

markerList = { ...
    'o', ...
    's', ...
    '^', ...
    'v', ...
    'd', ...
    '>', ...
    '<', ...
    'p', ...
    'h', ...
    'x', ...
    '+'};

%% Prepare the sorting variable
sortValuesAll = [];

if ~isempty(mode.sortVariable) && ...
        ismember(mode.sortVariable, results.Properties.VariableNames)

    candidateSortValues = results.(mode.sortVariable);

    if isnumeric(candidateSortValues) && ...
            numel(candidateSortValues) == height(results)

        sortValuesAll = candidateSortValues(:);
        sortValuesAll = sortValuesAll(validRows);
    end
end

%% Plot each calculation vector separately
for plotIndex = 1:numberOfCurves

    curveNumber = curveNumbers(plotIndex);
    rowsInCurve = curveIndex == curveNumber;

    curveTemperature = temperature(rowsInCurve);
    curvePressure    = pressure(rowsInCurve);

    if ~isempty(sortValuesAll)
        curveSortValues = sortValuesAll(rowsInCurve);
        [~, sortOrder] = sort(curveSortValues);

    else
        [~, sortOrder] = sort(curvePressure);
    end

    curveTemperature = curveTemperature(sortOrder);
    curvePressure    = curvePressure(sortOrder);

    markerIndex = mod(plotIndex - 1, numel(markerList)) + 1;

    plot( ...
        ax, ...
        curveTemperature, ...
        curvePressure, ...
        '-', ...
        'Color', lineColors(plotIndex, :), ...
        'Marker', markerList{markerIndex}, ...
        'LineWidth', 1.4, ...
        'MarkerSize', 5, ...
        'DisplayName', char(curveLabels(curveNumber)));
end

%% Format the axes
xlabel(ax, 'Temperature (°C)');
ylabel(ax, 'Pressure (kbar)');

%% Construct the plot title
titleParts = { ...
    mode.displayName, ...
    char(string(groupName)), ...
    char(string(funcName))};

titleParts = titleParts(~cellfun(@isempty, titleParts));

title( ...
    ax, ...
    strjoin(titleParts, ' | '), ...
    'Interpreter', 'none');

%% Add a legend
if numberOfCurves > 0
    legend( ...
        ax, ...
        'show', ...
        'Location', 'best', ...
        'Interpreter', 'none');
end

hold(ax, 'off');

end


%% ====================================================================== %
function [curveIndex, curveLabels] = identifyDatasetCurves_( ...
        results, ...
        rangeVariable)
%IDENTIFYDATASETCURVES_ Identify individual calculation vectors.
%
% A new curve begins when:
%   1) One or more dataCode values change, or
%   2) The range variable resets from a higher value to a lower value.
%
% The second criterion keeps repeated calculations with the same dataset
% combination as separate curves.

numberOfRows = height(results);

curveIndex  = ones(numberOfRows, 1);
curveLabels = strings(0, 1);

if numberOfRows == 0
    return;
end

variableNames = results.Properties.VariableNames;

%% Find dataset-identifier variables
dataCodeVariables = variableNames( ...
    contains(variableNames, 'dataCode', 'IgnoreCase', true));

newCurve = false(numberOfRows, 1);
newCurve(1) = true;

%% Detect changes in data-code combinations
for variableIndex = 1:numel(dataCodeVariables)

    variableName = dataCodeVariables{variableIndex};

    values = string(results.(variableName));
    values = values(:);
    values(ismissing(values)) = "";

    newCurve(2:end) = ...
        newCurve(2:end) | ...
        values(2:end) ~= values(1:end-1);
end

%% Detect the reset of the pressure or temperature range
if ~isempty(rangeVariable) && ...
        ismember(rangeVariable, variableNames)

    rangeValues = results.(rangeVariable);

    if isnumeric(rangeValues) && ...
            numel(rangeValues) == numberOfRows

        rangeValues = rangeValues(:);

        previousValues = rangeValues(1:end-1);
        currentValues  = rangeValues(2:end);

        validPairs = ...
            isfinite(previousValues) & ...
            isfinite(currentValues);

        resetDetected = false(numberOfRows - 1, 1);

        tolerance = 1e-12;

        resetDetected(validPairs) = ...
            currentValues(validPairs) < ...
            previousValues(validPairs) - tolerance;

        newCurve(2:end) = ...
            newCurve(2:end) | resetDetected;
    end
end

%% Assign sequential curve numbers
curveIndex = cumsum(newCurve);

numberOfCurves = max(curveIndex);
curveLabels = strings(numberOfCurves, 1);
baseLabels  = strings(numberOfCurves, 1);

%% Construct a legend label for each calculation vector
for curveNumber = 1:numberOfCurves

    firstRow = find(curveIndex == curveNumber, 1, 'first');

    if isempty(dataCodeVariables)
        baseLabel = "Dataset " + string(curveNumber);

    else
        labelParts = strings(1, numel(dataCodeVariables));

        for variableIndex = 1:numel(dataCodeVariables)

            variableName = dataCodeVariables{variableIndex};
            value = string(results.(variableName)(firstRow));

            if ismissing(value) || strlength(value) == 0
                value = "missing";
            end

            labelParts(variableIndex) = value;
        end

        baseLabel = strjoin(labelParts, " + ");
    end

    baseLabels(curveNumber) = baseLabel;

    repeatedCount = sum( ...
        baseLabels(1:curveNumber-1) == baseLabel);

    if repeatedCount == 0
        curveLabels(curveNumber) = baseLabel;
    else
        curveLabels(curveNumber) = ...
            baseLabel + " [" + string(repeatedCount + 1) + "]";
    end
end

end


%% ====================================================================== %
function [rawdata_struct, info] = importCationDataset_()
%IMPORTCATIONDATASET_ Read an Excel or CSV file into a struct.
%
% Outputs:
%   rawdata_struct:
%       For an Excel file, fields are sheet names and values are tables.
%       For a CSV file, the field "data" contains the imported table.
%
%   info:
%       fileFullPath
%       fileName
%       ext
%       sheetNames

[fileName, filePath] = uigetfile( ...
    {'*.xlsx;*.csv', 'Excel or CSV files (*.xlsx, *.csv)'}, ...
    'Select a CSV or Excel file', ...
    'MultiSelect', 'off');

if isequal(fileName, 0)
    error('File selection was canceled.');
end

fileFullPath = fullfile(filePath, fileName);

disp([fileFullPath ' has been imported as raw data.']);

[~, baseName, ext] = fileparts(fileFullPath);
ext = lower(ext);

rawdata_struct = struct();
sheetNames = {};

switch ext

    case '.xlsx'
        [~, sheetNames] = xlsfinfo(fileFullPath);

        if isempty(sheetNames)
            error( ...
                'No sheets were found in the selected Excel file: %s', ...
                fileFullPath);
        end

        for sheetIndex = 1:numel(sheetNames)
            tableData = readtable( ...
                fileFullPath, ...
                'Sheet', sheetNames{sheetIndex});

            rawdata_struct.(sheetNames{sheetIndex}) = tableData;
        end

    case '.csv'
        tableData = readtable(fileFullPath);
        rawdata_struct.data = tableData;

    otherwise
        error( ...
            'Unsupported file extension: %s. Select an Excel or CSV file.', ...
            ext);
end

info = struct();
info.fileFullPath = fileFullPath;
info.fileName     = baseName;
info.ext          = ext;
info.sheetNames   = sheetNames;

end
