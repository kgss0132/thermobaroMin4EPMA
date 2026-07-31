function [results, baroFunc, groupName] = startBaroCalc_rangeT(rawdata_struct)
%STARTBAROCALC_RANGET
% Select and run a geobarometer over a specified temperature range.
%
% Workflow:
%   1) Select a barometer group from Geobarometer_list.xlsx.
%   2) Select a barometer by Reference and Name. Type is not displayed.
%   3) Select the temperature unit.
%   4) Enter the minimum and maximum temperatures.
%   5) Generate a temperature vector.
%   6) Run functions/+baro/+<group>/<mfile>.m once using the vector.
%   7) Standardize the pressure variables used for downstream plotting.
%   8) Attach temperature-range metadata to the returned table.
%
% Expected barometer-module call:
%   results = baro.<group>.<mfile>(rawdata_struct, T_degreeC_range)
%
% Outputs:
%   results   : Result table. An empty table is returned if canceled.
%   baroFunc  : Barometer function name without ".m".
%   groupName : Selected sheet/group name.
%
% Required variables for downstream P-T plotting:
%   T_degreeC
%   P_kbar

disp('--------------------------------------------------------------');
disp('=== Geobarometer selection and temperature-range calculation ===');

results   = table();
baroFunc  = '';
groupName = '';

%% Range settings
nRangePoints = 101;

%% Input validation
if nargin < 1
    error('startBaroCalc_rangeT requires rawdata_struct as input.');
end

if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end

%% Locate the barometer-list workbook
excelName = 'Geobarometer_list.xlsx';
excelFile = locateBaroExcel_(excelName);

try
    sh = sheetnames(excelFile);
    sheetNames = cellstr(sh);
catch
    [~, sheetNames] = xlsfinfo(excelFile);
end

if isempty(sheetNames)
    error('No sheets were found in: %s', excelFile);
end

%% 1) Select a barometer group
[idGroup, ok] = listdlg( ...
    'PromptString', 'Select a barometer mineral group:', ...
    'SelectionMode', 'single', ...
    'ListString', sheetNames, ...
    'ListSize', [360 320]);

if ~ok || isempty(idGroup)
    disp('Barometer-group selection was canceled.');
    return;
end

groupName = sheetNames{idGroup};

disp(['Selected barometer group: ' groupName]);

%% Read the selected sheet
moduli = readcell(excelFile, 'Sheet', groupName);

if size(moduli, 1) < 4
    error( ...
        'Sheet "%s" has fewer than four rows. A four-row format is required.', ...
        groupName);
end

sheetWidth = size(moduli, 2);

if sheetWidth < 2
    error( ...
        'Sheet "%s" does not contain any barometer columns.', ...
        groupName);
end

%% 2) Build and display the barometer list
% Display format: Reference (Barometer). Row 3 (Type) is not displayed.
listBaro = cell(1, sheetWidth - 1);

for columnIndex = 2:sheetWidth
    barometerName = toStr_(moduli{1, columnIndex});
    reference     = toStr_(moduli{2, columnIndex});

    listBaro{columnIndex - 1} = ...
        [reference ' (' barometerName ')'];
end

[idBaro, ok] = listdlg( ...
    'PromptString', 'Select a barometer:', ...
    'SelectionMode', 'single', ...
    'ListString', listBaro, ...
    'ListSize', [560 360]);

if ~ok || isempty(idBaro)
    disp('Barometer selection was canceled.');

    results   = table();
    baroFunc  = '';
    groupName = '';
    return;
end

selectedColumn = idBaro + 1;
mfileBase      = toStr_(moduli{4, selectedColumn});

if isempty(mfileBase)
    error( ...
        ['The selected barometer has an empty function name ' ...
         'in sheet "%s".'], ...
        groupName);
end

baroFunc = mfileBase;

disp(['Selected barometer function: ' baroFunc]);

%% 3) Select the temperature unit
unitList = {'degreeC', 'K'};

[idUnit, ok] = listdlg( ...
    'PromptString', 'Select the temperature unit:', ...
    'SelectionMode', 'single', ...
    'ListString', unitList, ...
    'ListSize', [220 160]);

if ~ok || isempty(idUnit)
    disp('Temperature-unit selection was canceled.');

    results   = table();
    baroFunc  = '';
    groupName = '';
    return;
end

TUnit = unitList{idUnit};

%% 4) Enter the minimum and maximum temperatures
prompt = { ...
    sprintf('Enter the minimum temperature [%s]:', TUnit), ...
    sprintf('Enter the maximum temperature [%s]:', TUnit)};

answer = inputdlg( ...
    prompt, ...
    'Temperature range', ...
    [1 45; 1 45], ...
    {'', ''});

if isempty(answer)
    disp('Temperature-range input was canceled.');

    results   = table();
    baroFunc  = '';
    groupName = '';
    return;
end

TInputMin = str2double(strtrim(answer{1}));
TInputMax = str2double(strtrim(answer{2}));

if ~isfinite(TInputMin) || ~isfinite(TInputMax)
    errordlg( ...
        'Both temperature limits must be finite numeric values.', ...
        'Invalid temperature range');

    error('The temperature-range values are invalid.');
end

if TInputMax <= TInputMin
    errordlg( ...
        'The maximum temperature must be greater than the minimum temperature.', ...
        'Invalid temperature range');

    error('The maximum temperature must exceed the minimum temperature.');
end

%% Convert the temperature limits to degrees Celsius
switch TUnit
    case 'degreeC'
        TMinDegreeC = TInputMin;
        TMaxDegreeC = TInputMax;

    case 'K'
        if TInputMin < 0 || TInputMax < 0
            errordlg( ...
                'Temperatures expressed in K must be greater than or equal to zero.', ...
                'Invalid temperature range');

            error('Temperature values in K cannot be negative.');
        end

        TMinDegreeC = TInputMin - 273.15;
        TMaxDegreeC = TInputMax - 273.15;

    otherwise
        error('Unsupported temperature unit: %s', TUnit);
end

TMinK = TMinDegreeC + 273.15;
TMaxK = TMaxDegreeC + 273.15;

if TMinK < 0 || TMaxK < 0
    errordlg( ...
        'The specified temperature range is below absolute zero.', ...
        'Invalid temperature range');

    error('The specified temperature range is below absolute zero.');
end

%% 5) Generate the temperature vector
T_degreeC_range = ...
    linspace(TMinDegreeC, TMaxDegreeC, nRangePoints).';

disp('--------------------------------------------------------------');
disp([ ...
    'Input temperature range: ' ...
    num2str(TInputMin, '%.10g') ' to ' ...
    num2str(TInputMax, '%.10g') ' ' TUnit]);

disp([ ...
    'Converted temperature range: ' ...
    num2str(TMinDegreeC, '%.10g') ' to ' ...
    num2str(TMaxDegreeC, '%.10g') ' degreeC']);

disp([ ...
    'Equivalent Kelvin range: ' ...
    num2str(TMinK, '%.10g') ' to ' ...
    num2str(TMaxK, '%.10g') ' K']);

disp(['Number of temperature points: ' num2str(nRangePoints)]);

%% 6) Run the barometer module
functionsDir = fileparts(mfilename('fullpath'));

if ~contains(path, functionsDir)
    addpath(functionsDir);
end

funcQualified = ['baro.' groupName '.' mfileBase];

disp('--------------------------------------------------------------');
disp(['Running barometer: ' funcQualified]);
disp('Passing the complete temperature vector to the barometer module.');

try
    results = feval( ...
        funcQualified, ...
        rawdata_struct, ...
        T_degreeC_range);

catch ME
    disp('=== Barometer execution failed ===');
    disp(['Function: ' funcQualified]);
    disp(getReport(ME, 'extended'));
    rethrow(ME);
end

%% 7) Standardize pressure variables for downstream plotting
if istable(results) && height(results) > 0

    if ~ismember('P_kbar', results.Properties.VariableNames)

        if ismember('P_kb', results.Properties.VariableNames)
            results.P_kbar = results.P_kb;

        elseif ismember('Pressure_kbar', results.Properties.VariableNames)
            results.P_kbar = results.Pressure_kbar;

        elseif ismember('P_MPa', results.Properties.VariableNames)
            results.P_kbar = results.P_MPa / 100;

        elseif ismember('P_GPa', results.Properties.VariableNames)
            results.P_kbar = results.P_GPa * 10;

        elseif ismember('P_bar', results.Properties.VariableNames)
            results.P_kbar = results.P_bar / 1000;

        else
            warning( ...
                ['The barometer output does not contain a recognizable ' ...
                 'pressure variable. P_kbar could not be created.']);
        end
    end

    if ~ismember('P_MPa', results.Properties.VariableNames) && ...
            ismember('P_kbar', results.Properties.VariableNames)

        results.P_MPa = results.P_kbar * 100;
    end
end

%% 8) Attach temperature-range metadata
if istable(results) && height(results) > 0
    numberOfRows = height(results);

    results.T_range_input_min = repmat(TInputMin, numberOfRows, 1);
    results.T_range_input_max = repmat(TInputMax, numberOfRows, 1);
    results.T_range_input_unit = repmat(string(TUnit), numberOfRows, 1);

    results.T_range_min_degreeC = repmat(TMinDegreeC, numberOfRows, 1);
    results.T_range_max_degreeC = repmat(TMaxDegreeC, numberOfRows, 1);
    results.T_range_min_K       = repmat(TMinK, numberOfRows, 1);
    results.T_range_max_K       = repmat(TMaxK, numberOfRows, 1);
    results.T_range_npoints     = repmat(nRangePoints, numberOfRows, 1);

    metadataVariables = { ...
        'T_range_input_min', ...
        'T_range_input_max', ...
        'T_range_input_unit', ...
        'T_range_min_degreeC', ...
        'T_range_max_degreeC', ...
        'T_range_min_K', ...
        'T_range_max_K', ...
        'T_range_npoints'};

    results = moveMetadataVars_(results, metadataVariables);

    %% Attach the row-by-row temperature values when possible
    if ~ismember('T_degreeC', results.Properties.VariableNames)

        if mod(numberOfRows, nRangePoints) == 0
            numberOfDatasets = numberOfRows / nRangePoints;

            results.T_degreeC = repmat( ...
                T_degreeC_range, ...
                numberOfDatasets, ...
                1);

            results.T_K = results.T_degreeC + 273.15;

            results = moveMetadataVars_( ...
                results, ...
                {'T_degreeC', 'T_K'});
        else
            warning( ...
                ['The result table has %d rows, which is not an integer ' ...
                 'multiple of the %d temperature values. T_degreeC was ' ...
                 'not added because the row correspondence cannot be ' ...
                 'determined.'], ...
                numberOfRows, ...
                nRangePoints);
        end

    else
        validateRangeOutput_( ...
            results.T_degreeC, ...
            TMinDegreeC, ...
            TMaxDegreeC, ...
            'T_degreeC');

        if ~ismember('T_K', results.Properties.VariableNames)
            results.T_K = results.T_degreeC + 273.15;

            results = moveMetadataVars_( ...
                results, ...
                {'T_K'});
        end
    end

    %% Move standardized pressure variables to a convenient position
    pressureVariables = { ...
        'P_kbar', ...
        'P_MPa'};

    results = moveMetadataVars_( ...
        results, ...
        pressureVariables);
end

end


%% ====================================================================== %
function excelFile = locateBaroExcel_(excelName)
%LOCATEBAROEXCEL Locate the geobarometer-list workbook.

thisDir = fileparts(mfilename('fullpath'));

candidate1 = fullfile(thisDir, '+baro', excelName);
candidate2 = fullfile(fileparts(thisDir), '+baro', excelName);
candidate3 = fullfile(pwd, '+baro', excelName);

if exist(candidate1, 'file') == 2
    excelFile = candidate1;
    return;
end

if exist(candidate2, 'file') == 2
    excelFile = candidate2;
    return;
end

if exist(candidate3, 'file') == 2
    excelFile = candidate3;
    return;
end

locatedFile = which(excelName);

if ~isempty(locatedFile)
    excelFile = locatedFile;
    return;
end

error( ...
    'Cannot find %s under +baro or on the MATLAB path.', ...
    excelName);

end


%% ====================================================================== %
function results = moveMetadataVars_(results, variableNames)
%MOVEMETADATAVARS Move metadata columns after the final dataCode column.

existingVariables = variableNames( ...
    ismember(variableNames, results.Properties.VariableNames));

if isempty(existingVariables)
    return;
end

allVariableNames = results.Properties.VariableNames;
dataCodeIndices  = find(contains(allVariableNames, 'dataCode'));

if ~isempty(dataCodeIndices)
    lastDataCodeVariable = allVariableNames{max(dataCodeIndices)};

    results = movevars( ...
        results, ...
        existingVariables, ...
        'After', ...
        lastDataCodeVariable);
else
    results = movevars( ...
        results, ...
        existingVariables, ...
        'Before', ...
        1);
end

end


%% ====================================================================== %
function validateRangeOutput_(values, minimumValue, maximumValue, variableName)
%VALIDATERANGEOUTPUT Check whether returned values lie in the input range.

if ~isnumeric(values)
    warning( ...
        'The returned variable %s is not numeric.', ...
        variableName);
    return;
end

validValues = values(isfinite(values));

if isempty(validValues)
    warning( ...
        'The returned variable %s does not contain finite values.', ...
        variableName);
    return;
end

tolerance = 1e-9;

if any(validValues < minimumValue - tolerance) || ...
        any(validValues > maximumValue + tolerance)

    warning( ...
        ['Some returned %s values are outside the specified range ' ...
         'of %.10g to %.10g.'], ...
        variableName, ...
        minimumValue, ...
        maximumValue);
end

end


%% ====================================================================== %
function outputString = toStr_(inputValue)
%TOSTR Convert readcell output to a scalar character vector safely.

if isempty(inputValue)
    outputString = '';
    return;
end

try
    if any(ismissing(inputValue(:)))
        outputString = '';
        return;
    end
catch
end

if isstring(inputValue)
    outputString = char(inputValue(1));

elseif ischar(inputValue)
    outputString = inputValue;

elseif isnumeric(inputValue) || islogical(inputValue)
    outputString = num2str(inputValue(1));

else
    try
        outputString = char(string(inputValue));
    catch
        outputString = '';
    end
end

outputString = strtrim(outputString);

end
