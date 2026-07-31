function [results, thermoFunc, groupName] = startThermoCalc_fixedP(rawdata_struct)
%STARTTHERMOCALC_FIXEDP Select and run a geothermometer at a fixed pressure.
%
% Workflow:
%   1) Select a geothermometer group or target mineral from the sheet names
%      in Geothermometer_list.xlsx.
%   2) Select a geothermometer using the Geothermometer and Reference
%      information in the selected sheet. Type is not displayed.
%   3) Select the pressure unit and enter a fixed pressure.
%   4) Convert the input pressure to kbar and pass it to
%      functions/+thermo/+<group>/<mfile>.m as:
%          thermo.<group>.<mfile>
%   5) Return the calculation results to the calling function.
%
% Outputs:
%   results    : Result table. An empty table is returned if canceled.
%   thermoFunc : Base name of the selected thermometer function.
%                An empty character vector is returned if canceled.
%   groupName  : Selected sheet or group name.
%                An empty character vector is returned if canceled.

disp('--------------------------------------------------------------');
disp('=== Geothermometer selection and fixed-pressure calculation ===');

results    = table();
thermoFunc = '';
groupName  = '';

%% Validate input
if nargin < 1
    error('startThermoCalc_fixedP requires rawdata_struct as input.');
end

if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end

%% Locate the geothermometer-list workbook
excelName = 'Geothermometer_list.xlsx';
excelFile = locateThermoExcel_(excelName);

% Read sheet names, which represent thermometer groups or target minerals.
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
[id_group, ok] = listdlg( ...
    'PromptString', 'Select a thermometer mineral group:', ...
    'SelectionMode', 'single', ...
    'ListString', sheetNames, ...
    'ListSize', [360 320]);

if ~ok || isempty(id_group)
    disp('Thermometer-group selection was canceled.');
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
        'Sheet "%s" does not contain any thermometer columns.', ...
        groupName);
end

%% 2) Build the thermometer list from rows 1-2 and columns 2-w
% Display format:
%   Reference (Geothermometer)
% Row 3 (Type) remains available in the workbook but is not displayed.
listThermo = cell(1, w - 1);

for c = 2:w
    name = toStr_(moduli{1, c});
    ref  = toStr_(moduli{2, c});

    listThermo{c - 1} = [ref ' (' name ')'];
end

[id_th, ok] = listdlg( ...
    'PromptString', 'Select a thermometer:', ...
    'SelectionMode', 'single', ...
    'ListString', listThermo, ...
    'ListSize', [560 360]);

if ~ok || isempty(id_th)
    disp('Thermometer selection was canceled.');
    return;
end

col = id_th + 1;

% Read the function name without the ".m" extension.
mfileBase = toStr_(moduli{4, col});

if isempty(mfileBase)
    error( ...
        ['The selected thermometer has an empty function name ' ...
         'in sheet "%s".'], ...
        groupName);
end

thermoFunc = mfileBase;

disp(['Selected thermometer function: ' thermoFunc]);

%% 3) Enter the fixed pressure
unitList = {'kbar', 'MPa'};

[id_unit, ok] = listdlg( ...
    'PromptString', 'Select the pressure unit:', ...
    'SelectionMode', 'single', ...
    'ListString', unitList, ...
    'ListSize', [220 160]);

if ~ok || isempty(id_unit)
    disp('Pressure-unit selection was canceled.');

    results    = table();
    thermoFunc = '';
    groupName  = '';
    return;
end

P_unit = unitList{id_unit};

answer = inputdlg( ...
    {sprintf('Enter the pressure [%s]:', P_unit)}, ...
    'Fixed pressure', ...
    [1 45], ...
    {''});

if isempty(answer)
    disp('Pressure input was canceled.');

    results    = table();
    thermoFunc = '';
    groupName  = '';
    return;
end

P_value = str2double(strtrim(answer{1}));

if ~isfinite(P_value) || P_value < 0
    errordlg( ...
        ['The pressure must be entered as a finite numeric value ' ...
         'greater than or equal to zero.'], ...
        'Invalid pressure');

    error('Invalid pressure value.');
end

%% Convert the pressure to kbar
switch P_unit
    case 'kbar'
        P_kbar = P_value;

    case 'MPa'
        % 1 kbar = 100 MPa
        P_kbar = P_value / 100;

    otherwise
        error('Unsupported pressure unit: %s', P_unit);
end

disp([ ...
    'Input pressure: ' num2str(P_value) ' ' P_unit ...
    '  (= ' num2str(P_kbar) ' kbar)']);

%% 4) Run the thermometer module
% The module is called as:
%   thermo.<groupName>.<mfileBase>

functionsDir = fileparts(mfilename('fullpath'));

if ~contains(path, functionsDir)
    addpath(functionsDir);
end

funcQualified = ['thermo.' groupName '.' mfileBase];

disp('--------------------------------------------------------------');
disp(['Running thermometer: ' funcQualified]);

try
    results = feval( ...
        funcQualified, ...
        rawdata_struct, ...
        P_kbar);

catch ME
    disp('=== Thermometer execution failed ===');
    disp(['Function: ' funcQualified]);
    disp(getReport(ME, 'extended'));

    rethrow(ME);
end

%% 5) Attach pressure metadata and validate propagation
if istable(results) && height(results) > 0
    n = height(results);

    % Preserve the original user input.
    results.P_input_value = repmat(P_value, n, 1);
    results.P_input_unit  = repmat(string(P_unit), n, 1);

    % Add standardized pressure values.
    results.P_kbar = repmat(P_kbar, n, 1);
    results.P_MPa  = repmat(P_kbar * 100, n, 1);

    % Move the pressure metadata to a convenient position.
    insertVars = { ...
        'P_input_value', ...
        'P_input_unit', ...
        'P_MPa', ...
        'P_kbar'};

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

    % Confirm that the thermometer output uses the specified pressure.
    if ismember('P_kbar', results.Properties.VariableNames)
        Pmod = results.P_kbar;
        valid = isfinite(Pmod);

        if any(valid)
            if any(abs(Pmod(valid) - P_kbar) > 1e-9)
                warning( ...
                    ['Thermometer output P_kbar does not match the ' ...
                     'input pressure. Input = %.6g; ' ...
                     'output unique values = %s'], ...
                    P_kbar, ...
                    mat2str(unique(Pmod(valid))'));
            end
        end
    end
end

end


%% ====================================================================== %
function excelFile = locateThermoExcel_(excelName)
%LOCATETHERMOEXCEL Locate the geothermometer-list workbook.

thisDir = fileparts(mfilename('fullpath'));

% Search priority:
%   1) <functions>/+thermo/Geothermometer_list.xlsx
%   2) <project_root>/+thermo/Geothermometer_list.xlsx
%   3) <current_folder>/+thermo/Geothermometer_list.xlsx
%   4) MATLAB path

candidate1 = fullfile(thisDir, '+thermo', excelName);
candidate2 = fullfile(fileparts(thisDir), '+thermo', excelName);
candidate3 = fullfile(pwd, '+thermo', excelName);

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
    'Cannot find %s under +thermo or on the MATLAB path.', ...
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
