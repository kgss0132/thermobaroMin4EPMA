function results = BreyKohler1990Na(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/BreyKohler1990Na.m
% Tested with MATLAB R2024b
%
% Na-partition two-pyroxene thermometer for garnet lherzolites
% Brey, G.P., Kohler, T. (1990)
% Journal of Petrology, 31(6), 1353–1378
% DOI: https://doi.org/10.1093/petrology/31.6.1353
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis and one
% Clinopyroxene analysis and calculates temperature using the Na-partition
% thermometer of Brey & Kohler (1990).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Opx–Cpx pair and stores the results in
% a single output table.
%
% Both scalar and vector pressure inputs are accepted. For each selected
% Opx–Cpx pair, one output row is returned for every pressure value. This
% allows the function to be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP.
%
% -------------------------------------------------------------------------
% CALIBRATION / APPLICATION RANGE AND APPLICATION NOTES
%
% Brey & Kohler (1990) developed this thermometer from published mineral
% compositions of natural garnet lherzolites rather than from an independent
% set of reversed high-pressure experiments. Temperatures and pressures used
% for the regression were calculated with the authors' TBKN two-pyroxene
% thermometer and PBKN barometer. The Na-partition equation is presented as
% equation (11) on p. 1368.
%
% The natural-rock data discussed for this calibration include garnet
% lherzolites from Vitim, The Thumb, Namibia, and the Kaapvaal craton
% (pp. 1367–1368; Fig. 8). The Namibia data were excluded from the fitting
% because Na in Opx showed large scatter, probably related to strong
% serpentinization (p. 1368).
%
% Brey & Kohler (1990) did not state a strict rectangular experimental
% calibration range for equation (11). The following approximate limits are
% therefore used only for non-stopping range warnings in this implementation:
%
%   Temperature : approximately 700–1500 degreeC
%   Pressure    : approximately 10–70 kbar
%   Composition : natural garnet-lherzolite Opx–Cpx pairs
%
% These approximate limits represent the P–T domain occupied by the natural
% samples shown in Fig. 8c–d on p. 1367; they are not formal experimental
% calibration boundaries stated by the authors.
%
% Temperatures calculated with equation (11) agreed with TBKN temperatures
% to approximately +/-56 degreeC (1 sigma), with no systematic dependence on
% pressure or temperature (p. 1368; Fig. 8d). This agreement is relative to
% calculated TBKN temperatures, not to independently known experimental
% temperatures.
%
% The authors considered the thermometer especially promising at relatively
% low temperatures because Na partitioning provides greater temperature
% resolution than their equations (9) and (10) in that region (p. 1368).
%
% Accurate Na analyses are essential because Na in Opx may be low and close
% to analytical detection limits. Alteration, serpentinization, inclusions,
% poor background correction, or pairing Opx and Cpx from different
% generations may produce misleading temperatures.
%
% The calibration data shown in Fig. 8c on p. 1367 mainly occupy
% approximately -4 <= ln(DNa) <= -1, where DNa = Na_opx / Na_cpx and
% Na_opx < Na_cpx. This approximate compositional interval is used for an
% additional non-stopping warning. It is not a formally stated calibration
% boundary.
%
% This implementation therefore issues non-stopping warnings when:
%   1) input pressure is outside approximately 10–70 kbar,
%   2) a finite calculated temperature is outside approximately
%      700–1500 degreeC,
%   3) finite ln(DNa) is outside approximately -4 to -1,
%   4) an input cation value contains NaN,
%   5) Na_opx or Na_cpx is zero, or
%   6) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%   rawdata_struct.Cpx : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns should contain
% normalized pyroxene cations, preferably on a 6-oxygen basis.
%
% Required variables in both Opx and Cpx tables:
%   Si_cation_apfu         % stored for traceability
%   Al_cation_apfu         % stored for traceability
%   Fe_cation_apfu         % stored as total Fe
%   Mg_cation_apfu         % stored for traceability
%   Ca_cation_apfu         % stored for traceability
%   Na_cation_apfu         % used directly in the thermometer
%
% Optional variables in both tables:
%   Fe3_cation_apfu        % if the COLUMN is absent, Fe3+ is assumed to be 0
%   Mn_cation_apfu         % stored only; absent column is returned as NaN
%   Ti_cation_apfu         % stored only; absent column is returned as NaN
%   Cr_cation_apfu         % stored only; absent column is returned as NaN
%
% Important NaN and sign handling:
% - A NaN value present in an input column is retained as NaN. It is never
%   replaced by zero.
% - If Na_cation_apfu is NaN, DNa, ln(DNa), T_K, and T_deg remain NaN.
% - If an optional stored-only column exists but its selected value is NaN,
%   that NaN remains in the output. Only an entirely absent
%   Fe3_cation_apfu column invokes the explicit assumption Fe3+ = 0.
% - Finite mineral-composition values must be >= 0. Negative finite values
%   stop the calculation with an error.
% - Na_cation_apfu = 0 is allowed to remain in the output, but division and
%   logarithmic terms are undefined. DNa, ln(DNa), T_K, and T_deg are
%   therefore retained as NaN and reported by fprintf warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
%   T(K) = (35000 + 61.5 * P_kbar) / ((ln DNa)^2 + 19.8)
%
% where
%   DNa = Na_opx / Na_cpx
%
% Notes:
% - Pressure is supplied and used directly in kbar.
% - Temperature is calculated in Kelvin and converted to degree Celsius.
% - Na values are read directly from Na_cation_apfu.
% - The squared logarithmic term is mathematically symmetric for DNa and
%   1/DNa, but the calibration data have Na_opx < Na_cpx. Results for
%   DNa >= 1 should therefore be regarded as outside the observed
%   calibration behavior.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = BreyKohler1990Na(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Opx and Cpx tables (see above)
%   P_kbar         : pressure in kbar; finite non-negative numeric scalar or
%                    vector
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Opx–Cpx pair
%

%% Input validation
% Basic argument checks prevent silent failures caused by missing inputs or
% invalid pressure values.
if nargin < 2
    error('BreyKohler1990Na requires (rawdata_struct, P_kbar).');
end

if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end

if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

% Convert pressure to a column vector so all downstream table variables have
% a consistent orientation in both fixed-P and range-P modes.
P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The source tables are
% not modified.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end

if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

dataset_opx = rawdata_struct.Opx;
dataset_cpx = rawdata_struct.Cpx;

% Preserve the original output-variable interface by requiring the same
% basic cation columns as the source implementation. Only Na enters equation
% (11) directly; the other required variables are stored for traceability.
requiredVariables = {'Si_cation_apfu', 'Al_cation_apfu', ...
    'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu'};

validateRequiredVariables(dataset_opx, requiredVariables, 'Opx');
validateRequiredVariables(dataset_cpx, requiredVariables, 'Cpx');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each selected mineral pair produces one table block. Blocks are buffered
% in a preallocated cell array and concatenated only once after the
% interactive loop, avoiding repeated reallocation of the full results table
% on every iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Approximate natural-data P–T domain represented in Fig. 8. These are used
% only for warning messages and are not formal experimental limits.
applicationT_min_degC = 700;
applicationT_max_degC = 1500;
applicationP_min_kbar = 10;
applicationP_max_kbar = 70;

% Approximate compositional domain visible in Fig. 8c.
applicationLnDNa_min = -4;
applicationLnDNa_max = -1;

% Pressure is common to all mineral pairs selected during this function
% call, so its range warning is printed only once.
pressureOutsideApplication = ...
    P_kbar < applicationP_min_kbar | ...
    P_kbar > applicationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Opx) ===');

while true
    % ----- Opx selection -----
    % The first table column is treated as the identifier shown to the user.
    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Opx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Opx selected: ' char(string(selectedCode_opx))]);

    % ----- Cpx selection -----
    disp('=== Step 4: Selecting a data code from the list (Cpx) ===');

    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Cpx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIdx_cpx);
    disp(['Cpx selected: ' char(string(selectedCode_cpx))]);

    % ----- Calculation -----
    % Opx and Cpx are selected independently; row numbers are not assumed to
    % correspond between the two source tables.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_opx = dataset_opx(selectedIdx_opx, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    % Report all measured NaN cation values without interrupting the
    % calculation. Negative finite cation values are prohibited.
    validateNonNegativeInputs(selectedData_opx, selectedData_cpx);
    nanInputNames = findNaNInputs(selectedData_opx, selectedData_cpx);

    row = calcTemp(selectedData_opx, selectedData_cpx, P_kbar);

    % Repeat the selected identifiers for every pressure row.
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_opx', 'dataCode_cpx'}, 'Before', 1);

    % Store this pair as one buffered table block. Capacity is doubled only
    % when necessary, rather than changing the results table size on every
    % loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_opx)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_opx)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ...
            ' degreeC']);
    end

    % Warn once when pressure lies outside the approximate natural-data
    % domain shown in Fig. 8. The calculation and output are retained.
    if any(pressureOutsideApplication) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate natural-data ' ...
             'application domain represented in Fig. 8 of Brey & Kohler ' ...
             '(1990): approximately 10–70 kbar. %d of %d pressure point(s) ' ...
             'are outside the range; input range = %.4g–%.4g kbar. ' ...
             'The paper does not state a strict rectangular experimental ' ...
             'calibration range for equation (11).\n'], ...
            sum(pressureOutsideApplication), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the
    % approximate natural-data domain shown in Fig. 8d. NaN and Inf are
    % reported separately below.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideApplication = finiteTemperature & ...
        (row.T_deg < applicationT_min_degC | ...
         row.T_deg > applicationT_max_degC);

    if any(temperatureOutsideApplication)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximate ' ...
             'natural-data application domain represented in Fig. 8d of ' ...
             'Brey & Kohler (1990): approximately 700–1500 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s & %s. ' ...
             'These are approximate figure-based limits, not formal ' ...
             'experimental calibration boundaries.\n'], ...
            sum(temperatureOutsideApplication), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)));
    end

    % Report finite ln(DNa) values outside the approximate compositional
    % domain shown in Fig. 8c.
    finiteLnDNa = isfinite(row.lnDNa_BKN);
    lnDNaOutsideApplication = finiteLnDNa & ...
        (row.lnDNa_BKN < applicationLnDNa_min | ...
         row.lnDNa_BKN > applicationLnDNa_max);

    if any(lnDNaOutsideApplication)
        finiteValues = row.lnDNa_BKN(finiteLnDNa);
        fprintf(2, ...
            ['WARNING: ln(DNa) is outside the approximate range occupied by ' ...
             'the natural calibration/application data in Fig. 8c of ' ...
             'Brey & Kohler (1990): approximately -4 to -1. ' ...
             'Calculated finite ln(DNa) range = %.4g–%.4g for %s & %s. ' ...
             'Results outside this compositional domain should be treated ' ...
             'with caution.\n'], ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)));
    end

    % Print the exact input cation names containing NaN. fprintf is used
    % instead of warning so the message remains visible even when MATLAB
    % warning display has been disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the input cation value(s) for %s & %s: %s.\n' ...
             '         The NaN value(s) were retained; they were not replaced by zero.\n' ...
             '         NaN in Na_cation_apfu propagates to DNa and temperature, whereas NaN in a stored-only variable remains in that output column.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Zero Na is not negative, but it lies outside the division/logarithmic
    % domain. Report the exact zero input names without stopping.
    zeroNaInputNames = findZeroNaInputs(selectedData_opx, selectedData_cpx);

    if ~isempty(zeroNaInputNames)
        fprintf(2, ...
            ['WARNING: Zero was found in the Na thermometer input(s) for %s & %s: %s.\n' ...
             '         DNa, ln(DNa), T_K, and T_deg were retained as NaN because the ratio/logarithm is undefined.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(zeroNaInputNames, ', ')));
    end

    % Retain and report all NaN/Inf temperature values, whether caused by a
    % NaN input, zero Na, or another invalid mathematical domain.
    invalidTemperature = ~isfinite(row.T_deg);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    % Ask whether another Opx–Cpx pair should be calculated.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'BreyKohler1990Na', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered blocks once. Return an empty table if the user
% canceled before completing a calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(dataset, requiredVariables, mineralLabel)
% validateRequiredVariables
% Confirm that every required table variable exists before interactive
% selection begins.

missingVariables = strings(numel(requiredVariables), 1);
nMissing = 0;

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};

    if ~ismember(variableName, dataset.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingVariables(nMissing) = string(variableName);
    end
end

if nMissing > 0
    missingVariables = missingVariables(1:nMissing);
    error('%s table is missing required variable(s): %s', ...
        mineralLabel, char(strjoin(missingVariables, ', ')));
end

end

function nanInputNames = findNaNInputs(data_opx, data_cpx)
% findNaNInputs
% Return the names of input cation variables containing measured NaN values.
% This function never stops the calculation.

requiredVariables = {'Si_cation_apfu', 'Al_cation_apfu', ...
    'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu'};
optionalVariables = {'Mn_cation_apfu', 'Ti_cation_apfu', ...
    'Cr_cation_apfu', 'Fe3_cation_apfu'};

% Maximum possible count is ten variables for each of the two minerals.
nanInputNamesBuffer = strings(20, 1);
nFound = 0;

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    variableValue = data_opx.(variableName);

    if any(isnan(variableValue(:)))
        nFound = nFound + 1;
        nanInputNamesBuffer(nFound) = "Opx." + string(variableName);
    end
end

for i = 1:numel(optionalVariables)
    variableName = optionalVariables{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        variableValue = data_opx.(variableName);

        if any(isnan(variableValue(:)))
            nFound = nFound + 1;
            nanInputNamesBuffer(nFound) = "Opx." + string(variableName);
        end
    end
end

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    variableValue = data_cpx.(variableName);

    if any(isnan(variableValue(:)))
        nFound = nFound + 1;
        nanInputNamesBuffer(nFound) = "Cpx." + string(variableName);
    end
end

for i = 1:numel(optionalVariables)
    variableName = optionalVariables{i};

    if ismember(variableName, data_cpx.Properties.VariableNames)
        variableValue = data_cpx.(variableName);

        if any(isnan(variableValue(:)))
            nFound = nFound + 1;
            nanInputNamesBuffer(nFound) = "Cpx." + string(variableName);
        end
    end
end

nanInputNames = nanInputNamesBuffer(1:nFound);

end

function zeroNaInputNames = findZeroNaInputs(data_opx, data_cpx)
% findZeroNaInputs
% Return the names of finite Na thermometer inputs equal to zero.

zeroNaInputNamesBuffer = strings(2, 1);
nFound = 0;

Na_opx = data_opx.Na_cation_apfu;
Na_cpx = data_cpx.Na_cation_apfu;

if isfinite(Na_opx) && Na_opx == 0
    nFound = nFound + 1;
    zeroNaInputNamesBuffer(nFound) = "Opx.Na_cation_apfu";
end

if isfinite(Na_cpx) && Na_cpx == 0
    nFound = nFound + 1;
    zeroNaInputNamesBuffer(nFound) = "Cpx.Na_cation_apfu";
end

zeroNaInputNames = zeroNaInputNamesBuffer(1:nFound);

end

function validateNonNegativeInputs(data_opx, data_cpx)
% validateNonNegativeInputs
% Stop when any finite input cation value is negative. Zero and NaN are
% intentionally allowed: zero Na is handled as an invalid mathematical
% domain, and NaN is retained and reported by non-stopping warnings.

requiredVariables = {'Si_cation_apfu', 'Al_cation_apfu', ...
    'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu'};
optionalVariables = {'Mn_cation_apfu', 'Ti_cation_apfu', ...
    'Cr_cation_apfu', 'Fe3_cation_apfu'};

% Maximum possible count is ten variables for each of the two minerals.
invalidInputNamesBuffer = strings(20, 1);
nInvalid = 0;

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    variableValue = data_opx.(variableName);
    validateScalarNumericValue(variableValue, variableName, 'Opx');

    if isfinite(variableValue) && variableValue < 0
        nInvalid = nInvalid + 1;
        invalidInputNamesBuffer(nInvalid) = "Opx." + string(variableName);
    end
end

for i = 1:numel(optionalVariables)
    variableName = optionalVariables{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        variableValue = data_opx.(variableName);
        validateScalarNumericValue(variableValue, variableName, 'Opx');

        if isfinite(variableValue) && variableValue < 0
            nInvalid = nInvalid + 1;
            invalidInputNamesBuffer(nInvalid) = ...
                "Opx." + string(variableName);
        end
    end
end

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    variableValue = data_cpx.(variableName);
    validateScalarNumericValue(variableValue, variableName, 'Cpx');

    if isfinite(variableValue) && variableValue < 0
        nInvalid = nInvalid + 1;
        invalidInputNamesBuffer(nInvalid) = "Cpx." + string(variableName);
    end
end

for i = 1:numel(optionalVariables)
    variableName = optionalVariables{i};

    if ismember(variableName, data_cpx.Properties.VariableNames)
        variableValue = data_cpx.(variableName);
        validateScalarNumericValue(variableValue, variableName, 'Cpx');

        if isfinite(variableValue) && variableValue < 0
            nInvalid = nInvalid + 1;
            invalidInputNamesBuffer(nInvalid) = ...
                "Cpx." + string(variableName);
        end
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNamesBuffer(1:nInvalid);
    error(['BreyKohler1990Na: input cation values must be >= 0. ' ...
           'Negative finite value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_cpx, P_kbar)
% calcTemp
% Calculate the Na-partition temperature for one Opx row, one Cpx row, and
% one or more pressures. One output row is returned for every pressure value.
%
% Inputs:
%   data_opx : 1-row Opx table with normalized cations
%   data_cpx : 1-row Cpx table with normalized cations
%   P_kbar   : pressure in kbar; scalar or vector converted to a column
%
% Output:
%   row : table containing pressure, input compositions, intermediate terms,
%         and calculated temperatures

P_kbar = P_kbar(:);
nP = numel(P_kbar);

% Extract one-row mineral compositions while preserving measured NaN values.
opx = preparePyroxeneRow(data_opx, 'Opx');
cpx = preparePyroxeneRow(data_cpx, 'Cpx');

% Replicate each scalar mineral-composition value once for every pressure.
Si_opx = repmat(opx.Si, nP, 1);
Al_opx = repmat(opx.Al, nP, 1);
Fe_total_opx = repmat(opx.Fe_total, nP, 1);
Fe2_opx = repmat(opx.Fe2, nP, 1);
Fe3_opx = repmat(opx.Fe3, nP, 1);
Mg_opx = repmat(opx.Mg, nP, 1);
Ca_opx = repmat(opx.Ca, nP, 1);
Na_opx = repmat(opx.Na, nP, 1);
Mn_opx = repmat(opx.Mn, nP, 1);
Ti_opx = repmat(opx.Ti, nP, 1);
Cr_opx = repmat(opx.Cr, nP, 1);

Si_cpx = repmat(cpx.Si, nP, 1);
Al_cpx = repmat(cpx.Al, nP, 1);
Fe_total_cpx = repmat(cpx.Fe_total, nP, 1);
Fe2_cpx = repmat(cpx.Fe2, nP, 1);
Fe3_cpx = repmat(cpx.Fe3, nP, 1);
Mg_cpx = repmat(cpx.Mg, nP, 1);
Ca_cpx = repmat(cpx.Ca, nP, 1);
Na_cpx = repmat(cpx.Na, nP, 1);
Mn_cpx = repmat(cpx.Mn, nP, 1);
Ti_cpx = repmat(cpx.Ti, nP, 1);
Cr_cpx = repmat(cpx.Cr, nP, 1);

% Initialize all ratio/logarithmic outputs as NaN. Values are calculated only
% where both measured Na inputs are finite and strictly positive.
DNa_BKN = NaN(nP, 1);
lnDNa_BKN = NaN(nP, 1);

validNaInputs = isfinite(Na_opx) & isfinite(Na_cpx) & ...
    Na_opx > 0 & Na_cpx > 0;

DNa_BKN(validNaInputs) = ...
    Na_opx(validNaInputs) ./ Na_cpx(validNaInputs);

validDNa = isfinite(DNa_BKN) & DNa_BKN > 0;
lnDNa_BKN(validDNa) = log(DNa_BKN(validDNa));

% Published equation (11). The numerator remains finite because P_kbar was
% validated before calculation. Invalid ratio/logarithmic inputs remain NaN.
numerator_Na = 35000 + 61.5 .* P_kbar;
denominator_Na = (lnDNa_BKN .^ 2) + 19.8;

T_K = NaN(nP, 1);
T_deg = NaN(nP, 1);

validTemperatureCalculation = isfinite(numerator_Na) & ...
    isfinite(denominator_Na) & denominator_Na > 0;

T_K(validTemperatureCalculation) = ...
    numerator_Na(validTemperatureCalculation) ./ ...
    denominator_Na(validTemperatureCalculation);
T_deg(validTemperatureCalculation) = ...
    T_K(validTemperatureCalculation) - 273.15;

% Construct the complete output table in one operation so every variable has
% a fixed and consistent size.
row = table( ...
    P_kbar, ...
    Si_opx, Al_opx, Fe_total_opx, Fe2_opx, Fe3_opx, Mg_opx, Ca_opx, ...
    Na_opx, Mn_opx, Ti_opx, Cr_opx, ...
    Si_cpx, Al_cpx, Fe_total_cpx, Fe2_cpx, Fe3_cpx, Mg_cpx, Ca_cpx, ...
    Na_cpx, Mn_cpx, Ti_cpx, Cr_cpx, ...
    DNa_BKN, lnDNa_BKN, numerator_Na, denominator_Na, T_K, T_deg, ...
    'VariableNames', { ...
    'P_kbar', ...
    'Si_opx', 'Al_opx', 'Fe_total_opx', 'Fe2_opx', 'Fe3_opx', ...
    'Mg_opx', 'Ca_opx', 'Na_opx', 'Mn_opx', 'Ti_opx', 'Cr_opx', ...
    'Si_cpx', 'Al_cpx', 'Fe_total_cpx', 'Fe2_cpx', 'Fe3_cpx', ...
    'Mg_cpx', 'Ca_cpx', 'Na_cpx', 'Mn_cpx', 'Ti_cpx', 'Cr_cpx', ...
    'DNa_BKN', 'lnDNa_BKN', 'numerator_Na', 'denominator_Na', ...
    'T_K', 'T_deg'});

end

function px = preparePyroxeneRow(data_px, mineralLabel)
% preparePyroxeneRow
% Extract one pyroxene analysis while preserving measured NaN values.
% Missing stored-only variables are returned as NaN. If the entire Fe3
% column is absent, Fe3+ is explicitly assumed to be zero.

if height(data_px) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

px = struct();

px.Si = getRequiredValue(data_px, 'Si_cation_apfu', mineralLabel);
px.Al = getRequiredValue(data_px, 'Al_cation_apfu', mineralLabel);
px.Fe_total = getRequiredValue(data_px, 'Fe_cation_apfu', mineralLabel);
px.Mg = getRequiredValue(data_px, 'Mg_cation_apfu', mineralLabel);
px.Ca = getRequiredValue(data_px, 'Ca_cation_apfu', mineralLabel);
px.Na = getRequiredValue(data_px, 'Na_cation_apfu', mineralLabel);

px.Mn = getOptionalValue(data_px, 'Mn_cation_apfu', NaN, mineralLabel);
px.Ti = getOptionalValue(data_px, 'Ti_cation_apfu', NaN, mineralLabel);
px.Cr = getOptionalValue(data_px, 'Cr_cation_apfu', NaN, mineralLabel);
px.Fe3 = getOptionalValue(data_px, 'Fe3_cation_apfu', 0, mineralLabel);

% Fe3+ greater than total Fe is physically inconsistent. The check is made
% only when both values are finite so NaN continues to propagate.
if isfinite(px.Fe3) && isfinite(px.Fe_total) && px.Fe3 > px.Fe_total
    error('%s has Fe3_cation_apfu > Fe_cation_apfu.', mineralLabel);
end

px.Fe2 = px.Fe_total - px.Fe3;

end

function value = getRequiredValue(tbl, variableName, mineralLabel)
% getRequiredValue
% Read one required scalar numeric table value without modifying NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);
validateScalarNumericValue(value, variableName, mineralLabel);

end

function value = getOptionalValue(tbl, variableName, missingDefault, mineralLabel)
% getOptionalValue
% Read one optional scalar numeric value. The specified default is used only
% when the entire table column is absent; an existing NaN value is retained.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);
    validateScalarNumericValue(value, variableName, mineralLabel);
else
    value = missingDefault;
end

end

function validateScalarNumericValue(value, variableName, mineralLabel)
% validateScalarNumericValue
% Validate shape and type without rejecting NaN. Negative-value screening is
% handled separately by validateNonNegativeInputs.

if ~isnumeric(value) || ~isscalar(value) || ~isreal(value)
    error('%s variable %s must be a real numeric scalar in a 1-row table.', ...
        mineralLabel, variableName);
end

if isinf(value)
    error('%s variable %s must not be Inf.', mineralLabel, variableName);
end

end
