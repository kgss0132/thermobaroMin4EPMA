function results = BreyKohler1990BKN(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/BreyKohler1990BKN.m
% Tested with MATLAB R2024b
%
% Two-pyroxene thermometer for natural peridotitic systems
% Brey, G.P., Kohler, T. (1990)
% Journal of Petrology, 31(6), 1353–1378
% DOI: https://doi.org/10.1093/petrology/31.6.1353
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis and one
% Clinopyroxene analysis and calculates temperature using the BKN
% two-pyroxene thermometer of Brey & Kohler (1990).
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
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Brey & Kohler (1990) calibrated the final BKN two-pyroxene thermometer
% using experiments in simple systems and in a Cr-bearing natural
% peridotitic composition. The natural-system formulation is presented as
% equation (9) on p. 1363.
%
% The following limits are used for non-stopping range warnings in this
% implementation:
%
%   Temperature : approximately 900–1400 degreeC
%   Pressure    : 10–60 kbar
%   Composition : natural, Mg-rich peridotitic Opx–Cpx assemblages
%
% The natural-system pressure interval of 10–60 kbar is stated with the
% natural-system calibration on p. 1362. The experimental temperature span
% of approximately 900–1400 degreeC is shown in Fig. 6 on p. 1364.
%
% The final BKN thermometer reproduced the natural-system experiments to
% approximately +/-15 degreeC (1 sigma), without a systematic dependence on
% pressure, temperature, or the tested compositional parameters (pp.
% 1363–1364; Fig. 6).
%
% The Fe correction was tested for bulk Mg# values from 1.00 down to at
% least 0.89. The authors suggested that application may extend to bulk Mg#
% values of approximately 0.80–0.85 with minimal error, but this lower range
% was not directly covered by the natural-system calibration experiments
% (p. 1364).
%
% The Na correction used in Ca* is described by equations (7) and (8) on
% p. 1363. Na must therefore be measured and supplied; it must not be
% silently replaced by zero when the measured value is NaN.
%
% The thermometer assumes that the selected Opx and Cpx were coexisting and
% mutually equilibrated. Pairing analyses from different textural or
% chemical generations, exsolved host and lamella compositions, altered
% grains, or grains that re-equilibrated at different times may produce a
% numerically valid but geologically meaningless temperature.
%
% Brey & Kohler (1990) found that the combination of this thermometer with
% their Al-in-Opx barometer gave the most accurate simultaneous P–T
% estimates for natural peridotitic samples, reproducing the experiments to
% approximately +/-20 degreeC and +/-3 kbar (pp. 1374–1375; Fig. 11a).
%
% This implementation therefore issues non-stopping warnings when:
%   1) input pressure is outside 10–60 kbar,
%   2) a finite calculated temperature is outside 900–1400 degreeC,
%   3) a required thermometer input contains NaN, or
%   4) a calculated temperature is NaN or Inf.
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
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu         % total Fe used to calculate Fe2+
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu         % required for the published Na correction
%
% Optional variables in both tables:
%   Fe3_cation_apfu        % if the COLUMN is absent, Fe3+ is assumed to be 0
%   Mn_cation_apfu         % stored only; absent column is returned as NaN
%   Ti_cation_apfu         % stored only; absent column is returned as NaN
%   Cr_cation_apfu         % stored only; absent column is returned as NaN
%
% Important NaN and sign handling:
% - A NaN value present in any required thermometer input is retained as
%   NaN and propagated through the calculation. It is never replaced by 0.
% - If Fe3_cation_apfu exists but its selected value is NaN, that NaN is
%   retained and propagated. Only an entirely absent Fe3_cation_apfu column
%   invokes the explicit assumption Fe3+ = 0.
% - Finite values used in the thermometer must be >= 0. Negative finite
%   values stop the calculation with an error. Zero is allowed, although it
%   may lead to a NaN result when a denominator becomes zero.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
%   T(K) = [23664 + (24.9 + 126.3*XFe_cpx)*P_kbar]
%          / [13.38 + (ln KD)^2 + 11.59*XFe_opx]
%
% where
%   XFe_opx = Fe2_opx / (Fe2_opx + Mg_opx)
%   XFe_cpx = Fe2_cpx / (Fe2_cpx + Mg_cpx)
%
%   KD = (1 - CaStar_cpx) / (1 - CaStar_opx)
%
%   CaStar = Ca_M2 / (1 - Na_M2)
%
% For pyroxene structural formulae normalized to 6 oxygens, this
% implementation uses:
%   Ca_M2 = Ca_cation_apfu
%   Na_M2 = Na_cation_apfu
%
% Fe2+ is calculated as:
%   Fe2 = Fe_total - Fe3
%
% Notes:
% - Pressure is supplied and used directly in kbar.
% - Temperature is calculated in Kelvin and converted to degree Celsius.
% - If CaStar or KD falls outside the mathematical domain of the published
%   logarithmic expression, the affected result is set to NaN and retained.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = BreyKohler1990BKN(rawdata_struct, P_kbar)
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
    error('BreyKohler1990BKN requires (rawdata_struct, P_kbar).');
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

% Fail early when a required table variable is missing. Na is required
% because it enters the published Ca* correction directly.
requiredVariables = {'Si_cation_apfu', 'Al_cation_apfu', ...
    'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu'};

validateRequiredVariables(dataset_opx, requiredVariables, 'Opx');
validateRequiredVariables(dataset_cpx, requiredVariables, 'Cpx');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each selected mineral pair produces one table block. Blocks are buffered
% in a preallocated cell array and concatenated only once after the
% interactive loop, avoiding repeated reallocation of the full results
% table on every iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Natural-system calibration limits used for warning messages.
calibrationT_min_degC = 900;
calibrationT_max_degC = 1400;
calibrationP_min_kbar = 10;
calibrationP_max_kbar = 60;

% Pressure is common to all mineral pairs selected during this function
% call, so its range warning is printed only once.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
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

    % Identify NaN values only in variables that affect the thermometer.
    % Calculation continues and the exact input names are printed later.
    nanInputNames = findNaNInputs(selectedData_opx, selectedData_cpx);

    % Negative finite values are prohibited. NaN is deliberately allowed so
    % it remains NaN and propagates through the calculation.
    validateNonNegativeInputs(selectedData_opx, selectedData_cpx);

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

    % Warn once when any pressure lies outside 10–60 kbar. The calculation
    % and output are retained.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the natural-system calibration ' ...
             'range of Brey & Kohler (1990): 10–60 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the
    % approximate 900–1400 degreeC experimental range. NaN and Inf are
    % reported separately below.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximate ' ...
             'natural-system calibration range of Brey & Kohler (1990): ' ...
             '900–1400 degreeC. %d of %d finite temperature point(s) are ' ...
             'outside the range; calculated finite range = %.4g–%.4g ' ...
             'degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)));
    end

    % Print the exact thermometer input names containing NaN. fprintf is
    % used instead of warning so the message remains visible even when
    % MATLAB warning display has been disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         The NaN value(s) were retained and propagated; they were not replaced by zero.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all NaN/Inf temperature values, whether caused by a
    % NaN input, a zero denominator, or another invalid mathematical domain.
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
        'BreyKohler1990BKN', ...
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
% Return the names of thermometer input variables containing NaN. This
% function never stops the calculation.

requiredThermometerVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu'};
optionalThermometerVariable = 'Fe3_cation_apfu';

% Maximum possible count is five variables for each of the two minerals.
nanInputNamesBuffer = strings(10, 1);
nFound = 0;

for i = 1:numel(requiredThermometerVariables)
    variableName = requiredThermometerVariables{i};
    variableValue = data_opx.(variableName);
    if any(isnan(variableValue(:)))
        nFound = nFound + 1;
        nanInputNamesBuffer(nFound) = "Opx." + string(variableName);
    end
end

if ismember(optionalThermometerVariable, data_opx.Properties.VariableNames)
    variableValue = data_opx.(optionalThermometerVariable);
    if any(isnan(variableValue(:)))
        nFound = nFound + 1;
        nanInputNamesBuffer(nFound) = ...
            "Opx." + string(optionalThermometerVariable);
    end
end

for i = 1:numel(requiredThermometerVariables)
    variableName = requiredThermometerVariables{i};
    variableValue = data_cpx.(variableName);
    if any(isnan(variableValue(:)))
        nFound = nFound + 1;
        nanInputNamesBuffer(nFound) = "Cpx." + string(variableName);
    end
end

if ismember(optionalThermometerVariable, data_cpx.Properties.VariableNames)
    variableValue = data_cpx.(optionalThermometerVariable);
    if any(isnan(variableValue(:)))
        nFound = nFound + 1;
        nanInputNamesBuffer(nFound) = ...
            "Cpx." + string(optionalThermometerVariable);
    end
end

nanInputNames = nanInputNamesBuffer(1:nFound);

end

function validateNonNegativeInputs(data_opx, data_cpx)
% validateNonNegativeInputs
% Stop when a finite thermometer input is negative. Zero and NaN are
% intentionally allowed: zero may generate a NaN result, and NaN is retained
% and reported by non-stopping warnings.

requiredThermometerVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu'};
optionalThermometerVariable = 'Fe3_cation_apfu';

% Maximum possible count is five variables for each of the two minerals.
invalidInputNamesBuffer = strings(10, 1);
nInvalid = 0;

for i = 1:numel(requiredThermometerVariables)
    variableName = requiredThermometerVariables{i};
    variableValue = data_opx.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidInputNamesBuffer(nInvalid) = "Opx." + string(variableName);
    end
end

if ismember(optionalThermometerVariable, data_opx.Properties.VariableNames)
    variableValue = data_opx.(optionalThermometerVariable);
    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidInputNamesBuffer(nInvalid) = ...
            "Opx." + string(optionalThermometerVariable);
    end
end

for i = 1:numel(requiredThermometerVariables)
    variableName = requiredThermometerVariables{i};
    variableValue = data_cpx.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidInputNamesBuffer(nInvalid) = "Cpx." + string(variableName);
    end
end

if ismember(optionalThermometerVariable, data_cpx.Properties.VariableNames)
    variableValue = data_cpx.(optionalThermometerVariable);
    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidInputNamesBuffer(nInvalid) = ...
            "Cpx." + string(optionalThermometerVariable);
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNamesBuffer(1:nInvalid);
    error(['BreyKohler1990BKN: thermometer input values must be >= 0. ' ...
           'Negative finite value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_cpx, P_kbar)
% calcTemp
% Calculate the BKN temperature for one Opx row, one Cpx row, and one or
% more pressures. One output row is returned for every pressure value.
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

% Fe fractions. Zero denominators are allowed to produce NaN, which is
% retained and reported by the caller.
XFe_opx = Fe2_opx ./ (Fe2_opx + Mg_opx);
XFe_cpx = Fe2_cpx ./ (Fe2_cpx + Mg_cpx);

% Na-corrected Ca occupancies on the M2 site.
Ca_M2_opx = Ca_opx;
Na_M2_opx = Na_opx;
Ca_M2_cpx = Ca_cpx;
Na_M2_cpx = Na_cpx;

CaStar_opx = Ca_M2_opx ./ (1 - Na_M2_opx);
CaStar_cpx = Ca_M2_cpx ./ (1 - Na_M2_cpx);

% The logarithmic expression requires finite 0 <= CaStar < 1 values and a
% finite positive KD. Invalid points are converted to NaN rather than
% stopping the calculation or returning a misleading finite temperature.
validCaStar = isfinite(CaStar_opx) & isfinite(CaStar_cpx) & ...
    CaStar_opx >= 0 & CaStar_opx < 1 & ...
    CaStar_cpx >= 0 & CaStar_cpx < 1;

KD_BKN = NaN(nP, 1);
KD_BKN(validCaStar) = ...
    (1 - CaStar_cpx(validCaStar)) ./ ...
    (1 - CaStar_opx(validCaStar));

validKD = isfinite(KD_BKN) & KD_BKN > 0;
lnKD_BKN = NaN(nP, 1);
lnKD_BKN(validKD) = log(KD_BKN(validKD));

% Published BKN temperature equation. NaN values propagate naturally.
numerator_BKN = 23664 + ...
    (24.9 + 126.3 .* XFe_cpx) .* P_kbar;
denominator_BKN = 13.38 + ...
    (lnKD_BKN .^ 2) + 11.59 .* XFe_opx;

T_K = numerator_BKN ./ denominator_BKN;
T_deg = T_K - 273.15;

% Construct the complete output table in one operation so every variable has
% a fixed and consistent size.
row = table( ...
    P_kbar, ...
    Si_opx, Al_opx, Fe_total_opx, Fe2_opx, Fe3_opx, Mg_opx, Ca_opx, ...
    Na_opx, Mn_opx, Ti_opx, Cr_opx, ...
    Si_cpx, Al_cpx, Fe_total_cpx, Fe2_cpx, Fe3_cpx, Mg_cpx, Ca_cpx, ...
    Na_cpx, Mn_cpx, Ti_cpx, Cr_cpx, ...
    Ca_M2_opx, Na_M2_opx, CaStar_opx, ...
    Ca_M2_cpx, Na_M2_cpx, CaStar_cpx, ...
    XFe_opx, XFe_cpx, KD_BKN, lnKD_BKN, ...
    numerator_BKN, denominator_BKN, T_K, T_deg, ...
    'VariableNames', { ...
    'P_kbar', ...
    'Si_opx', 'Al_opx', 'Fe_total_opx', 'Fe2_opx', 'Fe3_opx', ...
    'Mg_opx', 'Ca_opx', 'Na_opx', 'Mn_opx', 'Ti_opx', 'Cr_opx', ...
    'Si_cpx', 'Al_cpx', 'Fe_total_cpx', 'Fe2_cpx', 'Fe3_cpx', ...
    'Mg_cpx', 'Ca_cpx', 'Na_cpx', 'Mn_cpx', 'Ti_cpx', 'Cr_cpx', ...
    'Ca_M2_opx', 'Na_M2_opx', 'CaStar_opx', ...
    'Ca_M2_cpx', 'Na_M2_cpx', 'CaStar_cpx', ...
    'XFe_opx', 'XFe_cpx', 'KD_BKN', 'lnKD_BKN', ...
    'numerator_BKN', 'denominator_BKN', 'T_K', 'T_deg'});

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

% Ca and Na are treated as M2-site occupancies in a 6-oxygen formula. Values
% with Ca + Na > 1 are structurally inconsistent and are rejected when both
% inputs are finite. Equality is allowed to propagate to NaN through the
% mathematical-domain check in calcTemp.
if isfinite(px.Ca) && isfinite(px.Na) && (px.Ca + px.Na) > 1 + 1e-8
    error('%s has Ca_cation_apfu + Na_cation_apfu > 1. Check pyroxene normalization.', ...
        mineralLabel);
end

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
