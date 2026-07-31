function [results, thermoFunc, groupName] = startThermoCalc_rangeP(rawdata_struct)
%STARTTHERMOCALC_RANGEP
% Select and run a geothermometer over a specified pressure range.
%
% Workflow:
%   1) Select a thermometer group from Geothermometer_list.xlsx.
%   2) Select a thermometer by Reference and Name. Type is not displayed.
%   3) Select the pressure unit.
%   4) Enter the minimum and maximum pressures.
%   5) Generate a pressure vector.
%   6) Run functions/+thermo/+<group>/<mfile>.m once using the vector.
%   7) Standardize the temperature variables used for downstream plotting.
%   8) Attach pressure-range metadata to the returned table.
%
% Expected thermometer-module call:
%   results = thermo.<group>.<mfile>(rawdata_struct, P_kbar_range)
%
% Outputs:
%   results    : Result table. An empty table is returned if canceled.
%   thermoFunc : Thermometer function name without ".m".
%   groupName  : Selected sheet/group name.
%
% Required variables for downstream P-T plotting:
%   T_degreeC
%   P_kbar

disp('--------------------------------------------------------------');
disp('=== Geothermometer selection and pressure-range calculation ===');

results    = table();
thermoFunc = '';
groupName  = '';

%% Range settings
nRangePoints = 101;

%% Input validation
if nargin < 1
    error('startThermoCalc_rangeP requires rawdata_struct as input.');
end

if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end

%% Locate the thermometer-list workbook
excelName = 'Geothermometer_list.xlsx';
excelFile = locateThermoExcel_(excelName);

try
    sh = sheetnames(excelFile);
    sheetNames = cellstr(sh);
catch
    [~, sheetNames] = xlsfinfo(excelFile);
end

if isempty(sheetNames)
    error('No sheets were found in: %s', excelFile);
end

%% 1) Select a thermometer group
[idGroup, ok] = listdlg( ...
    'PromptString', 'Select a thermometer mineral group:', ...
    'SelectionMode', 'single', ...
    'ListString', sheetNames, ...
    'ListSize', [360 320]);

if ~ok || isempty(idGroup)
    disp('Thermometer-group selection was canceled.');
    return;
end

groupName = sheetNames{idGroup};

disp(['Selected thermometer group: ' groupName]);

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
        'Sheet "%s" does not contain any thermometer columns.', ...
        groupName);
end

%% 2) Build and display the thermometer list
% Display format: Reference (Thermometer). Row 3 (Type) is not displayed.
listThermo = cell(1, sheetWidth - 1);

for columnIndex = 2:sheetWidth
    thermometerName = toStr_(moduli{1, columnIndex});
    reference       = toStr_(moduli{2, columnIndex});

    listThermo{columnIndex - 1} = ...
        [reference ' (' thermometerName ')'];
end

[idThermo, ok] = listdlg( ...
    'PromptString', 'Select a thermometer:', ...
    'SelectionMode', 'single', ...
    'ListString', listThermo, ...
    'ListSize', [560 360]);

if ~ok || isempty(idThermo)
    disp('Thermometer selection was canceled.');

    results    = table();
    thermoFunc = '';
    groupName  = '';
    return;
end

selectedColumn = idThermo + 1;
mfileBase      = toStr_(moduli{4, selectedColumn});

if isempty(mfileBase)
    error( ...
        ['The selected thermometer has an empty function name ' ...
         'in sheet "%s".'], ...
        groupName);
end

thermoFunc = mfileBase;

disp(['Selected thermometer function: ' thermoFunc]);

%% 3) Select the pressure unit
unitList = {'kbar', 'MPa'};

[idUnit, ok] = listdlg( ...
    'PromptString', 'Select the pressure unit:', ...
    'SelectionMode', 'single', ...
    'ListString', unitList, ...
    'ListSize', [220 160]);

if ~ok || isempty(idUnit)
    disp('Pressure-unit selection was canceled.');

    results    = table();
    thermoFunc = '';
    groupName  = '';
    return;
end

PUnit = unitList{idUnit};

%% 4) Enter the minimum and maximum pressures
prompt = { ...
    sprintf('Enter the minimum pressure [%s]:', PUnit), ...
    sprintf('Enter the maximum pressure [%s]:', PUnit)};

answer = inputdlg( ...
    prompt, ...
    'Pressure range', ...
    [1 45; 1 45], ...
    {'', ''});

if isempty(answer)
    disp('Pressure-range input was canceled.');

    results    = table();
    thermoFunc = '';
    groupName  = '';
    return;
end

PInputMin = str2double(strtrim(answer{1}));
PInputMax = str2double(strtrim(answer{2}));

if ~isfinite(PInputMin) || ~isfinite(PInputMax)
    errordlg( ...
        'Both pressure limits must be finite numeric values.', ...
        'Invalid pressure range');

    error('The pressure-range values are invalid.');
end

if PInputMin < 0 || PInputMax < 0
    errordlg( ...
        'Pressure values must be greater than or equal to zero.', ...
        'Invalid pressure range');

    error('Pressure values cannot be negative.');
end

if PInputMax <= PInputMin
    errordlg( ...
        'The maximum pressure must be greater than the minimum pressure.', ...
        'Invalid pressure range');

    error('The maximum pressure must exceed the minimum pressure.');
end

%% Convert the pressure limits to kbar
switch PUnit
    case 'kbar'
        PMinKbar = PInputMin;
        PMaxKbar = PInputMax;

    case 'MPa'
        PMinKbar = PInputMin / 100;
        PMaxKbar = PInputMax / 100;

    otherwise
        error('Unsupported pressure unit: %s', PUnit);
end

%% 5) Generate the pressure vector
P_kbar_range = linspace(PMinKbar, PMaxKbar, nRangePoints).';

disp('--------------------------------------------------------------');
disp([ ...
    'Input pressure range: ' ...
    num2str(PInputMin, '%.10g') ' to ' ...
    num2str(PInputMax, '%.10g') ' ' PUnit]);

disp([ ...
    'Converted pressure range: ' ...
    num2str(PMinKbar, '%.10g') ' to ' ...
    num2str(PMaxKbar, '%.10g') ' kbar']);

disp(['Number of pressure points: ' num2str(nRangePoints)]);

%% 6) Run the thermometer module
functionsDir = fileparts(mfilename('fullpath'));

if ~contains(path, functionsDir)
    addpath(functionsDir);
end

funcQualified = ['thermo.' groupName '.' mfileBase];

disp('--------------------------------------------------------------');
disp(['Running thermometer: ' funcQualified]);
disp('Passing the complete pressure vector to the thermometer module.');

try
    results = feval( ...
        funcQualified, ...
        rawdata_struct, ...
        P_kbar_range);

catch ME
    disp('=== Thermometer execution failed ===');
    disp(['Function: ' funcQualified]);
    disp(getReport(ME, 'extended'));
    rethrow(ME);
end

%% 7) Standardize temperature variables for downstream plotting
if istable(results) && height(results) > 0

    if ~ismember('T_degreeC', results.Properties.VariableNames)

        if ismember('T_deg', results.Properties.VariableNames)
            results.T_degreeC = results.T_deg;

        elseif ismember('T_C', results.Properties.VariableNames)
            results.T_degreeC = results.T_C;

        elseif ismember( ...
                'Temperature_degreeC', ...
                results.Properties.VariableNames)
            results.T_degreeC = results.Temperature_degreeC;

        elseif ismember('T_K', results.Properties.VariableNames)
            results.T_degreeC = results.T_K - 273.15;

        else
            warning( ...
                ['The thermometer output does not contain a recognizable ' ...
                 'temperature variable. T_degreeC could not be created.']);
        end
    end

    if ~ismember('T_K', results.Properties.VariableNames) && ...
            ismember('T_degreeC', results.Properties.VariableNames)

        results.T_K = results.T_degreeC + 273.15;
    end
end

%% 8) Attach pressure-range metadata
if istable(results) && height(results) > 0
    numberOfRows = height(results);

    results.P_range_input_min = repmat(PInputMin, numberOfRows, 1);
    results.P_range_input_max = repmat(PInputMax, numberOfRows, 1);
    results.P_range_input_unit = repmat(string(PUnit), numberOfRows, 1);

    results.P_range_min_kbar = repmat(PMinKbar, numberOfRows, 1);
    results.P_range_max_kbar = repmat(PMaxKbar, numberOfRows, 1);
    results.P_range_npoints  = repmat(nRangePoints, numberOfRows, 1);

    metadataVariables = { ...
        'P_range_input_min', ...
        'P_range_input_max', ...
        'P_range_input_unit', ...
        'P_range_min_kbar', ...
        'P_range_max_kbar', ...
        'P_range_npoints'};

    results = moveMetadataVars_(results, metadataVariables);

    %% Attach the row-by-row pressure values when possible
    if ~ismember('P_kbar', results.Properties.VariableNames)

        if numberOfRows == nRangePoints
            results.P_kbar = P_kbar_range;
            results.P_MPa  = P_kbar_range * 100;

            results = moveMetadataVars_( ...
                results, ...
                {'P_kbar', 'P_MPa'});
        else
            warning( ...
                ['The result table has %d rows, whereas the pressure ' ...
                 'vector has %d values. P_kbar was not added because ' ...
                 'the row correspondence cannot be determined.'], ...
                numberOfRows, ...
                nRangePoints);
        end

    else
        validateRangeOutput_( ...
            results.P_kbar, ...
            PMinKbar, ...
            PMaxKbar, ...
            'P_kbar');

        if ~ismember('P_MPa', results.Properties.VariableNames)
            results.P_MPa = results.P_kbar * 100;

            results = moveMetadataVars_( ...
                results, ...
                {'P_MPa'});
        end
    end

    %% Move standardized temperature variables to a convenient position
    temperatureVariables = { ...
        'T_degreeC', ...
        'T_K'};

    results = moveMetadataVars_( ...
        results, ...
        temperatureVariables);
end

end


%% ====================================================================== %
function excelFile = locateThermoExcel_(excelName)
%LOCATETHERMOEXCEL Locate the geothermometer-list workbook.

thisDir = fileparts(mfilename('fullpath'));

candidate1 = fullfile(thisDir, '+thermo', excelName);
candidate2 = fullfile(fileparts(thisDir), '+thermo', excelName);
candidate3 = fullfile(pwd, '+thermo', excelName);

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
    'Cannot find %s under +thermo or on the MATLAB path.', ...
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
