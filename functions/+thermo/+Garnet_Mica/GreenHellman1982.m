function results = GreenHellman1982(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/GreenHellman1982.m
% Tested with MATLAB R2024b
%
% Garnet-phengite Fe-Mg exchange thermometer
% Green, T.H. and Hellman, P.L. (1982)
% Lithos, 15, 253-266
% DOI: https://doi.org/10.1016/0024-4937(82)90017-2
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% analysis (selected by the user from tables) and calculates temperature
% using the Green & Hellman (1982) garnet-phengite Fe-Mg exchange
% thermometer.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Garnet-Mica pair and appends the result
% to a single output table. A scalar pressure returns one row per selected
% pair. A pressure vector returns one row per pressure value for every pair,
% allowing use from both startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Green & Hellman (1982) studied the garnet-phengite Fe-Mg exchange over an
% overall experimental range of 20-35 kbar and 800-1000 degreeC (abstract,
% p. 253; Composition A results, p. 256). However, the thermodynamic
% parameters used to derive the general equations were assumed constant over
% 20-35 kbar and 800-950 degreeC (p. 263). This implementation therefore uses
% the following conservative warning limits for the published equations:
%
%   Temperature : 800-950 degreeC
%   Pressure    : 20-35 kbar
%
% Experiments at 950-1000 degreeC frequently crystallized trioctahedral or
% approximately 2.5-octahedral mica rather than dioctahedral phengite
% (pp. 254, 256-257). The thermometer must be applied to garnet-phengite
% pairs, not indiscriminately to garnet paired with another mica type.
%
% The exchange relation is non-ideal and depends strongly on bulk-rock CaO
% and mg number. The authors state that it should not be used outside the
% investigated compositional ranges (abstract, p. 253; discussion,
% pp. 262-263):
%
%   Equation (a): low-Ca systems (< 2 wt.% CaO), bulk mg approximately 67;
%                 observed calibration.
%   Equation (b): low-Ca systems (< 2 wt.% CaO), bulk mg approximately
%                 20-30; interpolated calibration. Experiments with mg 23
%                 and 46 indicated little change across mg 23-46 (p. 257).
%   Equation (c): basaltic systems (approximately 10 wt.% CaO), bulk mg
%                 approximately 67; preliminary, interpolated calibration.
%
% The pressure effect was measured directly only for the low-Ca, mg
% approximately 67 composition. Its use for the Fe-rich and Ca-bearing
% calibrations is an assumption that may introduce error (pp. 260, 263).
% The basaltic calibration was described as preliminary and was estimated to
% reproduce temperature to within approximately 50 degreeC (p. 260).
%
% Severe experimental equilibrium problems occur below 800 degreeC, and
% equilibrium must be evaluated carefully even at higher temperature
% (p. 255). Natural garnet and phengite compositions must represent the same
% equilibrium assemblage and must not have been modified during later
% metamorphism (p. 265).
%
% The exchange coefficient uses Fe2+/Mg. If natural phengite Fe3+ is unknown
% and total Fe is treated as Fe2+, calculated temperatures are maximum
% estimates (abstract, p. 253; application, p. 264). In this implementation,
% Fe_cation_apfu is therefore used as Fe2+ and Fe3_cation_apfu is retained
% only as a reported diagnostic value; Fe3+ is not added to Fe2+ in K_D.
% Input data must consequently provide Fe_cation_apfu on an Fe2+ basis.
%
% The effect of Mn in garnet was not calibrated. The authors state that the
% thermometer should not be applied when garnet contains more than
% approximately 5 mol.% spessartine (p. 265). This function reports Mn when
% available but cannot calculate spessartine mol.% from Mn alone reliably.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 20-35 kbar;
%   2) a finite calculated temperature is outside 800-950 degreeC;
%   3) a required thermometer input contains NaN;
%   4) a required input is zero and a valid Fe-Mg ratio cannot be formed; or
%   5) a calculated temperature is NaN or Inf.
%
% The three equations represent different bulk compositions. Their three
% temperatures are not an uncertainty interval and are not interchangeable.
% Select only the result appropriate for the sample bulk composition.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Mica   : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns contain normalized
% cation data.
%
% Required variables in both Garnet and Mica tables:
%   Fe_cation_apfu         % Fe2+ used in K_D
%   Mg_cation_apfu
%
% Optional variables:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Finite mineral-composition values must be non-negative. Negative or Inf
% values stop the calculation with an error. NaN values are retained as
% missing values and propagated through the calculation. Zero is accepted as
% a non-negative analytical value, but a zero Fe2+ or Mg input cannot define
% the logarithmic exchange thermometer; the corresponding results are stored
% as NaN and reported by non-stopping warnings.
%
% Missing optional variables are stored as NaN rather than being interpreted
% as measured zero values.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Green & Hellman, 1982; pp. 263-264)
%
%   K_D = (Fe2+/Mg)_garnet / (Fe2+/Mg)_phengite
%
% (a) Low-Ca systems (< 2 wt.% CaO), bulk mg approximately 67
%     T(K) = (5560 + 0.036*P(bar)) / (ln(K_D) + 4.65)
%
% (b) Low-Ca systems (< 2 wt.% CaO), bulk mg approximately 20-30
%     T(K) = (5680 + 0.036*P(bar)) / (ln(K_D) + 4.48)
%
% (c) Basaltic systems, bulk mg approximately 67
%     T(K) = (5170 + 0.036*P(bar)) / (ln(K_D) + 4.17)
%
% Pressure is supplied in kbar and converted internally to bar.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = GreenHellman1982(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Mica tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Mica pair. Intermediate variables and all
%             three composition-specific temperatures are retained.

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('GreenHellman1982 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

% Standardize scalar, row-vector, and column-vector pressure inputs to one
% column-vector representation for stable table construction.
P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The tables are read but
% are not modified by this function.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_grt = rawdata_struct.Garnet;
dataset_mica = rawdata_struct.Mica;

requiredVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
validateRequiredVariables(dataset_grt, requiredVariables, 'Garnet');
validateRequiredVariables(dataset_mica, requiredVariables, 'Mica');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation is stored temporarily as one table block. Repeated table
% concatenation inside the interactive loop is avoided because it causes
% repeated reallocation and copying of the complete results table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Conservative equation limits documented above.
calibrationT_min_degC = 800;
calibrationT_max_degC = 950;
calibrationP_min_kbar = 20;
calibrationP_max_kbar = 35;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or selects
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
    % The first table column is treated as the display identifier.
    dataCodes_grt = dataset_grt{:, 1};

    [selectedIdx_grt, ok] = listdlg( ...
        'PromptString', 'Please select the Garnet data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_grt)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_grt)
        disp('Selection canceled');
        break;
    end

    selectedCode_grt = dataCodes_grt(selectedIdx_grt);
    disp(['Garnet selected: ' char(string(selectedCode_grt))]);

    % ----- Mica selection -----
    disp('=== Step 4: Selecting a data code from the list (Mica) ===');

    dataCodes_mica = dataset_mica{:, 1};

    [selectedIdx_mica, ok] = listdlg( ...
        'PromptString', 'Please select the Mica data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_mica)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_mica)
        disp('Selection canceled');
        break;
    end

    selectedCode_mica = dataCodes_mica(selectedIdx_mica);
    disp(['Mica selected: ' char(string(selectedCode_mica))]);

    % ----- Calculation -----
    % Garnet and mica are selected independently; matching row indices are not
    % assumed between the two datasets.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_mica = dataset_mica(selectedIdx_mica, :);

    % Validate all present cation variables before calculation. NaN is allowed
    % to pass through unchanged; finite negative values and Inf are rejected.
    validateNonNegativeInputs(selectedData_grt, selectedData_mica);

    % Check the values actually used in K_D. These checks prepare warning
    % messages only and do not replace NaN or zero with another value.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_mica);
    zeroInputNames = findZeroInputs(selectedData_grt, selectedData_mica);

    row = calcTemp(selectedData_grt, selectedData_mica, P_kbar);

    % Repeat the selected identifiers once per pressure value.
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_mica = repmat(string(selectedCode_mica), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_mica'}, 'Before', 1);

    % Store this calculation in a buffered cell. Capacity is doubled only
    % when needed, avoiding results-table growth on every loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': (a) = ' num2str(row.T_a_C) ' degreeC, (b) = ' ...
            num2str(row.T_b_C) ' degreeC, (c) = ' ...
            num2str(row.T_c_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': (a) = ' num2str(row.T_a_C(1)) ' to ' ...
            num2str(row.T_a_C(end)) ' degreeC, (b) = ' ...
            num2str(row.T_b_C(1)) ' to ' num2str(row.T_b_C(end)) ...
            ' degreeC, (c) = ' num2str(row.T_c_C(1)) ' to ' ...
            num2str(row.T_c_C(end)) ' degreeC']);
    end

    % Pressure is common to all selected pairs, so print this warning only
    % once per function call. The calculation is not stopped.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the conservative calibration ' ...
             'range of Green & Hellman (1982): 20-35 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar. See p. 263.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Each equation is checked separately because it represents a different
    % bulk-composition calibration.
    printTemperatureCalibrationWarning(row.T_a_C, '(a)', ...
        calibrationT_min_degC, calibrationT_max_degC, ...
        selectedCode_grt, selectedCode_mica);
    printTemperatureCalibrationWarning(row.T_b_C, '(b)', ...
        calibrationT_min_degC, calibrationT_max_degC, ...
        selectedCode_grt, selectedCode_mica);
    printTemperatureCalibrationWarning(row.T_c_C, '(c)', ...
        calibrationT_min_degC, calibrationT_max_degC, ...
        selectedCode_grt, selectedCode_mica);

    % Report missing required values without stopping the calculation.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; calculated temperatures may be NaN.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Zero is non-negative but cannot define a positive Fe/Mg exchange ratio.
    % The invalid result is intentionally NaN rather than 0 K.
    if ~isempty(zeroInputNames)
        fprintf(2, ...
            ['WARNING: Zero was found in the Fe-Mg thermometer input(s) for %s & %s: %s.\n' ...
             '         A positive K_D could not be defined; affected temperatures remain NaN.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            char(strjoin(zeroInputNames, ', ')));
    end

    % Retain and report any NaN/Inf results caused by invalid logarithmic or
    % denominator conditions. No result is replaced with zero.
    printNonFiniteResultWarning(row, selectedCode_grt, selectedCode_mica);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'GreenHellman1982', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once. Return an empty table if the
% user canceled before performing a calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(dataTable, requiredVariables, mineralLabel)
% validateRequiredVariables
% Confirm that each mineral table contains the variables required by K_D.

missingMask = ~ismember(requiredVariables, dataTable.Properties.VariableNames);
if any(missingMask)
    missingVariables = string(requiredVariables(missingMask));
    error('%s table must contain variable(s): %s', ...
        mineralLabel, char(strjoin(missingVariables, ', ')));
end

end

function nanInputNames = findNaNInputs(data_garnet, data_mica)
% findNaNInputs
% Return required K_D input names containing NaN. No values are modified.

inputNames = ["Garnet.Fe_cation_apfu"; "Garnet.Mg_cation_apfu"; ...
    "Mica.Fe_cation_apfu"; "Mica.Mg_cation_apfu"];
inputValues = [data_garnet.Fe_cation_apfu; data_garnet.Mg_cation_apfu; ...
    data_mica.Fe_cation_apfu; data_mica.Mg_cation_apfu];

nanInputNames = inputNames(isnan(inputValues));

end

function zeroInputNames = findZeroInputs(data_garnet, data_mica)
% findZeroInputs
% Return required K_D input names equal to zero. Zero is allowed by the
% non-negative validation but cannot produce a valid logarithmic K_D.

inputNames = ["Garnet.Fe_cation_apfu"; "Garnet.Mg_cation_apfu"; ...
    "Mica.Fe_cation_apfu"; "Mica.Mg_cation_apfu"];
inputValues = [data_garnet.Fe_cation_apfu; data_garnet.Mg_cation_apfu; ...
    data_mica.Fe_cation_apfu; data_mica.Mg_cation_apfu];

zeroInputNames = inputNames(isfinite(inputValues) & inputValues == 0);

end

function validateNonNegativeInputs(data_garnet, data_mica)
% validateNonNegativeInputs
% Stop when a present finite cation value is negative or when a value is Inf.
% NaN and zero are allowed here and handled by non-stopping warnings.

variableNames = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu', 'Mn_cation_apfu', 'Ca_cation_apfu', ...
    'Ti_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu', ...
    'K_cation_apfu', 'Na_cation_apfu'};
mineralTables = {data_garnet, data_mica};
mineralLabels = {'Garnet', 'Mica'};

maximumInvalidNames = numel(variableNames) * numel(mineralTables);
invalidInputNames = strings(maximumInvalidNames, 1);
nInvalidInputs = 0;

for mineralIndex = 1:numel(mineralTables)
    currentTable = mineralTables{mineralIndex};

    for variableIndex = 1:numel(variableNames)
        variableName = variableNames{variableIndex};

        if ~ismember(variableName, currentTable.Properties.VariableNames)
            continue;
        end

        variableValue = currentTable.(variableName);

        if ~isnumeric(variableValue) || ~isscalar(variableValue)
            error('%s.%s must be a numeric scalar in the selected row.', ...
                mineralLabels{mineralIndex}, variableName);
        end

        if isinf(variableValue) || (isfinite(variableValue) && variableValue < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputNames(nInvalidInputs) = ...
                string(mineralLabels{mineralIndex}) + "." + string(variableName);
        end
    end
end

invalidInputNames = invalidInputNames(1:nInvalidInputs);

if ~isempty(invalidInputNames)
    error(['GreenHellman1982: mineral-composition values must be ' ...
           'non-negative or NaN. Negative or Inf value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_garnet, data_mica, P_kbar)
% calcTemp
% Calculate all three Green & Hellman (1982) estimates for one garnet-mica
% pair at every supplied pressure value.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_bar = P_bar;

% Extract one-row mineral data. Optional missing values remain NaN.
grt = prepareMineralRow(data_garnet, 'Garnet');
mica = prepareMineralRow(data_mica, 'Mica');

% Repeat scalar mineral values once per pressure for a stable output schema.
Fe2_grt = repmat(grt.Fe2, nP, 1);
Fe3_grt = repmat(grt.Fe3, nP, 1);
FeUsed_grt = Fe2_grt;
Mg_grt = repmat(grt.Mg, nP, 1);
Mn_grt = repmat(grt.Mn, nP, 1);
Ca_grt = repmat(grt.Ca, nP, 1);
Al_grt = repmat(grt.Al, nP, 1);
Si_grt = repmat(grt.Si, nP, 1);

Fe2_mica = repmat(mica.Fe2, nP, 1);
Fe3_mica = repmat(mica.Fe3, nP, 1);
FeUsed_mica = Fe2_mica;
Mg_mica = repmat(mica.Mg, nP, 1);
Mn_mica = repmat(mica.Mn, nP, 1);
Ca_mica = repmat(mica.Ca, nP, 1);
Ti_mica = repmat(mica.Ti, nP, 1);
Al_mica = repmat(mica.Al, nP, 1);
Si_mica = repmat(mica.Si, nP, 1);
K_mica = repmat(mica.K, nP, 1);
Na_mica = repmat(mica.Na, nP, 1);

% Green & Hellman (1982):
% K_D = (Fe2+/Mg)_garnet / (Fe2+/Mg)_phengite.
% Arithmetic is allowed to propagate NaN. Zero-derived non-positive or
% non-finite K_D values are explicitly invalidated below so that they cannot
% produce a misleading 0 K (-273.15 degreeC) result through log(0).
FeMg_grt = FeUsed_grt ./ Mg_grt;
FeMg_mica = FeUsed_mica ./ Mg_mica;
KD = FeMg_grt ./ FeMg_mica;

validKD = isfinite(KD) & KD > 0;
lnKD = NaN(nP, 1);
lnKD(validKD) = log(KD(validKD));

% Equation numerators vary with pressure; denominators depend on K_D.
num_a = 5560 + 0.036 .* P_bar;
num_b = 5680 + 0.036 .* P_bar;
num_c = 5170 + 0.036 .* P_bar;

den_a = lnKD + 4.65;
den_b = lnKD + 4.48;
den_c = lnKD + 4.17;

T_a_K = solveTemperature(num_a, den_a);
T_b_K = solveTemperature(num_b, den_b);
T_c_K = solveTemperature(num_c, den_c);

T_a_C = T_a_K - 273.15;
T_b_C = T_b_K - 273.15;
T_c_C = T_c_K - 273.15;

% Pack inputs, intermediate values, and final results.
row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.FeUsed_grt = FeUsed_grt;
row.Mg_grt = Mg_grt;
row.Mn_grt = Mn_grt;
row.Ca_grt = Ca_grt;
row.Al_grt = Al_grt;
row.Si_grt = Si_grt;

row.Fe2_mica = Fe2_mica;
row.Fe3_mica = Fe3_mica;
row.FeUsed_mica = FeUsed_mica;
row.Mg_mica = Mg_mica;
row.Mn_mica = Mn_mica;
row.Ca_mica = Ca_mica;
row.Ti_mica = Ti_mica;
row.Al_mica = Al_mica;
row.Si_mica = Si_mica;
row.K_mica = K_mica;
row.Na_mica = Na_mica;

row.FeMg_grt = FeMg_grt;
row.FeMg_mica = FeMg_mica;
row.KD = KD;
row.lnKD = lnKD;

row.T_a_K = T_a_K;
row.T_a_C = T_a_C;
row.T_b_K = T_b_K;
row.T_b_C = T_b_C;
row.T_c_K = T_c_K;
row.T_c_C = T_c_C;

end

function temperature_K = solveTemperature(numerator, denominator)
% solveTemperature
% Return NaN unless both terms are finite and the denominator is positive.
% This prevents invalid logarithmic inputs from appearing as 0 K.

temperature_K = NaN(size(numerator));
validCalculation = isfinite(numerator) & isfinite(denominator) & ...
    denominator > 0;
temperature_K(validCalculation) = ...
    numerator(validCalculation) ./ denominator(validCalculation);

end

function mineral = prepareMineralRow(dataTable, mineralLabel)
% prepareMineralRow
% Extract one-row mineral cation data. Missing optional variables and present
% NaN values remain NaN; no missing value is converted to zero.

if height(dataTable) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();
mineral.Fe2 = getVarOrError(dataTable, 'Fe_cation_apfu', mineralLabel);
mineral.Mg = getVarOrError(dataTable, 'Mg_cation_apfu', mineralLabel);

mineral.Fe3 = getVarOrNaN(dataTable, 'Fe3_cation_apfu', mineralLabel);
mineral.Mn = getVarOrNaN(dataTable, 'Mn_cation_apfu', mineralLabel);
mineral.Ca = getVarOrNaN(dataTable, 'Ca_cation_apfu', mineralLabel);
mineral.Ti = getVarOrNaN(dataTable, 'Ti_cation_apfu', mineralLabel);
mineral.Al = getVarOrNaN(dataTable, 'Al_cation_apfu', mineralLabel);
mineral.Si = getVarOrNaN(dataTable, 'Si_cation_apfu', mineralLabel);
mineral.K = getVarOrNaN(dataTable, 'K_cation_apfu', mineralLabel);
mineral.Na = getVarOrNaN(dataTable, 'Na_cation_apfu', mineralLabel);

end

function value = getVarOrError(dataTable, variableName, mineralLabel)
% getVarOrError
% Read a required scalar variable while retaining NaN unchanged.

if ~ismember(variableName, dataTable.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = dataTable.(variableName);
validateScalarCationValue(value, variableName, mineralLabel);

end

function value = getVarOrNaN(dataTable, variableName, mineralLabel)
% getVarOrNaN
% Read an optional scalar variable. Missing variables are represented by NaN.

if ismember(variableName, dataTable.Properties.VariableNames)
    value = dataTable.(variableName);
    validateScalarCationValue(value, variableName, mineralLabel);
else
    value = NaN;
end

end

function validateScalarCationValue(value, variableName, mineralLabel)
% validateScalarCationValue
% Accept non-negative finite values and NaN; reject negative values and Inf.

if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end
if isinf(value) || (isfinite(value) && value < 0)
    error('%s.%s must be non-negative or NaN.', mineralLabel, variableName);
end

end

function printTemperatureCalibrationWarning(temperature_degC, equationLabel, ...
        calibrationMinimum, calibrationMaximum, selectedCode_grt, selectedCode_mica)
% printTemperatureCalibrationWarning
% Print a non-stopping range warning for one composition-specific equation.

finiteTemperature = isfinite(temperature_degC);
outsideCalibration = finiteTemperature & ...
    (temperature_degC < calibrationMinimum | ...
     temperature_degC > calibrationMaximum);

if any(outsideCalibration)
    finiteValues = temperature_degC(finiteTemperature);
    fprintf(2, ...
        ['WARNING: Calculated temperature from equation %s is outside the ' ...
         'conservative Green & Hellman (1982) range of 800-950 degreeC. ' ...
         '%d of %d finite point(s) are outside the range; calculated finite ' ...
         'range = %.4g-%.4g degreeC for %s & %s. See p. 263.\n'], ...
        equationLabel, ...
        sum(outsideCalibration), ...
        sum(finiteTemperature), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

end

function printNonFiniteResultWarning(row, selectedCode_grt, selectedCode_mica)
% printNonFiniteResultWarning
% Report NaN/Inf temperature values while retaining them in the output table.

temperatureNames = {'T_a_C', 'T_b_C', 'T_c_C'};
equationLabels = {'(a)', '(b)', '(c)'};

for equationIndex = 1:numel(temperatureNames)
    temperatureValues = row.(temperatureNames{equationIndex});
    invalidTemperature = ~isfinite(temperatureValues);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated with equation %s ' ...
             'for %s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table; the calculation was not stopped.\n'], ...
            equationLabels{equationIndex}, ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            sum(invalidTemperature), ...
            numel(temperatureValues), ...
            sum(isnan(temperatureValues)), ...
            sum(isinf(temperatureValues)));
    end
end

end
