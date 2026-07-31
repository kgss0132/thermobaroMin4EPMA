function [results, baroFunc, groupName] = startBaroCalc_fixedT(rawdata_struct)
%STARTBAROCALC_FIXEDT Select and run a geobarometer at a fixed temperature.
%
% Workflow:
%   1) Select a geobarometer group or target mineral from the sheet names
%      in Geobarometer_list.xlsx.
%   2) Select a geobarometer using the Geobarometer and Reference
%      information in the selected sheet. Type is not displayed.
%   3) Select the temperature unit and enter a fixed temperature.
%   4) Convert the input temperature to degrees Celsius and pass it to
%      functions/+baro/+<group>/<mfile>.m as:
%          baro.<group>.<mfile>
%   5) Return the calculation results to the calling function.
%
% Outputs:
%   results   : Result table. An empty table is returned if canceled.
%   baroFunc  : Base name of the selected barometer function.
%               An empty character vector is returned if canceled.
%   groupName : Selected sheet or group name.
%               An empty character vector is returned if canceled.

disp('--------------------------------------------------------------');
disp('=== Geobarometer selection and fixed-temperature calculation ===');

results   = table();
baroFunc  = '';
groupName = '';

%% Validate input
if nargin < 1
    error('startBaroCalc_fixedT requires rawdata_struct as input.');
end

if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end

%% Locate the geobarometer-list workbook
excelName = 'Geobarometer_list.xlsx';
excelFile = locateBaroExcel_(excelName);

% Read sheet names, which represent barometer groups or target minerals.
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
[id_group, ok] = listdlg( ...
    'PromptString', 'Select a barometer mineral group:', ...
    'SelectionMode', 'single', ...
    'ListString', sheetNames, ...
    'ListSize', [360 320]);

if ~ok || isempty(id_group)
    disp('Barometer-group selection was canceled.');
    return;
end

groupName = sheetNames{id_group};

disp(['Selected group: ' groupName]);

%% Read the selected sheet
moduli = readcell(excelFile, 'Sheet', groupName);

if size(moduli, 1) < 4
    error( ...
        'Sheet "%s" has fewer than four rows. A four-row format is required.', ...
        groupName);
end

w = size(moduli, 2);

if w < 2
    error( ...
        'Sheet "%s" does not contain any barometer columns.', ...
        groupName);
end

%% 2) Build the barometer list from rows 1-2 and columns 2-w
% Display format:
%   Reference (Geobarometer)
% Row 3 (Type) remains available in the workbook but is not displayed.
listBaro = cell(1, w - 1);

for c = 2:w
    name = toStr_(moduli{1, c});
    ref  = toStr_(moduli{2, c});

    listBaro{c - 1} = [ref ' (' name ')'];
end

[id_baro, ok] = listdlg( ...
    'PromptString', 'Select a barometer:', ...
    'SelectionMode', 'single', ...
    'ListString', listBaro, ...
    'ListSize', [560 360]);

if ~ok || isempty(id_baro)
    disp('Barometer selection was canceled.');
    return;
end

col = id_baro + 1;

% Read the function name without the ".m" extension.
mfileBase = toStr_(moduli{4, col});

if isempty(mfileBase)
    error( ...
        ['The selected barometer has an empty function name ' ...
         'in sheet "%s".'], ...
        groupName);
end

baroFunc = mfileBase;

disp(['Selected barometer function: ' baroFunc]);

%% 3) Enter the fixed temperature
unitList = {'degreeC', 'K'};

[id_unit, ok] = listdlg( ...
    'PromptString', 'Select the temperature unit:', ...
    'SelectionMode', 'single', ...
    'ListString', unitList, ...
    'ListSize', [220 160]);

if ~ok || isempty(id_unit)
    disp('Temperature-unit selection was canceled.');

    results   = table();
    baroFunc  = '';
    groupName = '';
    return;
end

T_unit = unitList{id_unit};

answer = inputdlg( ...
    {sprintf('Enter the temperature [%s]:', T_unit)}, ...
    'Fixed temperature', ...
    [1 45], ...
    {''});

if isempty(answer)
    disp('Temperature input was canceled.');

    results   = table();
    baroFunc  = '';
    groupName = '';
    return;
end

T_value = str2double(strtrim(answer{1}));

if ~isfinite(T_value)
    errordlg( ...
        'The temperature must be entered as a finite numeric value.', ...
        'Invalid temperature');

    error('Invalid temperature value.');
end

%% Convert the temperature to degrees Celsius
switch T_unit
    case 'degreeC'
        T_degreeC = T_value;

    case 'K'
        if T_value < 0
            errordlg( ...
                'A temperature expressed in K must be greater than or equal to zero.', ...
                'Invalid temperature');

            error('Invalid temperature value in K.');
        end

        T_degreeC = T_value - 273.15;

    otherwise
        error('Unsupported temperature unit: %s', T_unit);
end

T_K = T_degreeC + 273.15;

disp([ ...
    'Input temperature: ' num2str(T_value) ' ' T_unit ...
    '  (= ' num2str(T_degreeC) ' degreeC, ' ...
    num2str(T_K) ' K)']);

%% 4) Run the barometer module
% The module is called as:
%   baro.<groupName>.<mfileBase>

functionsDir = fileparts(mfilename('fullpath'));

if ~contains(path, functionsDir)
    addpath(functionsDir);
end

funcQualified = ['baro.' groupName '.' mfileBase];

disp('--------------------------------------------------------------');
disp(['Running barometer: ' funcQualified]);

try
    results = feval( ...
        funcQualified, ...
        rawdata_struct, ...
        T_degreeC);

catch ME
    disp('=== Barometer execution failed ===');
    disp(['Function: ' funcQualified]);
    disp(getReport(ME, 'extended'));

    rethrow(ME);
end

%% 5) Attach temperature metadata and validate propagation
if istable(results) && height(results) > 0
    n = height(results);

    % Preserve the original user input.
    results.T_input_value = repmat(T_value, n, 1);
    results.T_input_unit  = repmat(string(T_unit), n, 1);

    % Add standardized temperature values.
    results.T_degreeC = repmat(T_degreeC, n, 1);
    results.T_K       = repmat(T_K, n, 1);

    % Move the temperature metadata to a convenient position.
    insertVars = { ...
        'T_input_value', ...
        'T_input_unit', ...
        'T_degreeC', ...
        'T_K'};

    varNames = results.Properties.VariableNames;

    % Find all columns whose names contain "dataCode".
    idxDataCode = find(contains(varNames, 'dataCode'));

    if ~isempty(idxDataCode)
        % Place the metadata after the rightmost dataCode column.
        lastDataCodeCol = varNames{max(idxDataCode)};

        results = movevars( ...
            results, ...
            insertVars, ...
            'After', ...
            lastDataCodeCol);
    else
        % Place the metadata at the beginning if no dataCode column exists.
        results = movevars( ...
            results, ...
            insertVars, ...
            'Before', ...
            1);
    end

    % Confirm that the barometer output uses the specified temperature.
    if ismember('T_degreeC', results.Properties.VariableNames)
        Tmod = results.T_degreeC;
        valid = isfinite(Tmod);

        if any(valid)
            if any(abs(Tmod(valid) - T_degreeC) > 1e-9)
                warning( ...
                    ['Barometer output T_degreeC does not match the ' ...
                     'input temperature. Input = %.6g; ' ...
                     'output unique values = %s'], ...
                    T_degreeC, ...
                    mat2str(unique(Tmod(valid))'));
            end
        end
    end
end

end


%% ====================================================================== %
function excelFile = locateBaroExcel_(excelName)
%LOCATEBAROEXCEL Locate the geobarometer-list workbook.

thisDir = fileparts(mfilename('fullpath'));

% Search priority:
%   1) <functions>/+baro/Geobarometer_list.xlsx
%   2) <project_root>/+baro/Geobarometer_list.xlsx
%   3) <current_folder>/+baro/Geobarometer_list.xlsx
%   4) MATLAB path

candidate1 = fullfile(thisDir, '+baro', excelName);
candidate2 = fullfile(fileparts(thisDir), '+baro', excelName);
candidate3 = fullfile(pwd, '+baro', excelName);

if exist(candidate1, 'file') == 2
    excelFile = candidate1;
    return;

elseif exist(candidate2, 'file') == 2
    excelFile = candidate2;
    return;

elseif exist(candidate3, 'file') == 2
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
function outputString = toStr_(inputValue)
%TOSTR Convert readcell output to a scalar character vector safely.

if isempty(inputValue)
    outputString = '';
    return;
end

% Return an empty character vector for missing values.
try
    if any(ismissing(inputValue(:)))
        outputString = '';
        return;
    end
catch
    % Continue if ismissing is not applicable to the input type.
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
