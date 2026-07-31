function results = Wang2008(rawdata_struct, P_kbar)
% functions/+thermo/+Olivine_Spinel/Wang2008.m
% Tested with MATLAB R2024b
%
% Al partitioning geothermometer between Olivine and Spinel
% Wan, Z., Coogan, L.A., Canil, D. (2008)
% American Mineralogist, 93, 1142–1147
% DOI: https://doi.org/10.2138/am.2008.2758
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis and one Spinel
% analysis (selected by the user from tables) and calculates temperature
% using the Al2O3 partitioning thermometer of Wan et al. (2008).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Ol–Sp pair, and appends results into a
% single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Wan et al. (2008) experimentally calibrated Al2O3 partitioning between
% forsterite-rich olivine and Cr-rich spinel over the following conditions:
%
%   Temperature : 1250–1450 degreeC
%   Pressure    : 100 kPa (approximately 1 bar or 0.001 kbar)
%   Olivine     : Fo87–Fo93, centered approximately on Fo90
%   Spinel YCr  : 0.07–0.69, where YCr = Cr / (Cr + Al)
%   Spinel Fe3+ : < 0.1 atoms per 4 O
%   Spinel Ti   : < 0.025 atoms per 4 O
%   Redox range : approximately FMQ−1.5 to FMQ−1.8
%
% The calibration conditions and equation are summarized in the abstract
% on p. 1142. Experimental methods and analytical limitations are described
% on pp. 1143–1145, and Equation 3 is presented on p. 1145.
%
% Wan et al. (2008) report that Equation 3 reproduces the calibration
% experiments within approximately ±22 degreeC (pp. 1145–1146).
%
% APPLICATION NOTES:
% - The thermometer was directly calibrated only at 100 kPa. Published
%   comparison with experiments at 1 GPa showed no systematic pressure
%   offset, suggesting that pressure has only a minor effect (p. 1146).
%   This implementation therefore uses 0–1 GPa (0–10 kbar) only as an
%   empirical pressure-screening interval, not as a formal calibration
%   range. The equation itself contains no pressure term.
% - Natural mantle xenoliths at approximately 700–1100 degreeC showed no
%   systematic offset from a two-pyroxene thermometer, but with substantial
%   scatter. This is an empirical low-temperature test and does not extend
%   the direct experimental calibration below 1250 degreeC (pp. 1146–1147).
% - Preliminary experiments indicate that the thermometer cannot be used
%   for Fe3+-rich spinel when Fe3+/(Fe3+ + Cr + Al) > 0.15 (p. 1146).
% - Olivine commonly contains only trace Al. Wan et al. (2008) used
%   LA-ICP-MS for low-Al natural olivine because the concentrations were too
%   low for reliable routine electron-microprobe analysis (p. 1144).
% - Olivine and spinel must represent an equilibrium pair from the same
%   textural and compositional domain. Zoned, altered, or oxidized spinel
%   should be screened carefully.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside the empirical 0–1 GPa screening interval,
%   2) a finite calculated temperature is outside 1250–1450 degreeC,
%   3) spinel composition is outside the recommended YCr, Fe3+, or Ti range,
%   4) Fe3+/(Fe3+ + Cr + Al) exceeds 0.15,
%   5) a required input is NaN or zero, or
%   6) a calculated intermediate or temperature value is non-finite.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%   rawdata_struct.Spinel  : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include:
%
%   Olivine table variables:
%     Al2O3
%
%   Spinel table variables:
%     Cr_cation_apfu
%     Al_cation_apfu
%     Al2O3
%     Fe3_cation_apfu
%     Ti_cation_apfu
%
% Al2O3 values are used directly in wt%:
%   KD = Al2O3_ol / Al2O3_sp
%
% Spinel YCr is calculated from atomic proportions:
%   YCr = Cr / (Cr + Al)
%
% Negative finite values in required mineral-composition variables are not
% permitted. Zero values are retained but reported by non-stopping warnings
% because they may produce undefined logarithms or ratios. NaN values are
% retained as missing values, propagated through the calculation, and
% reported by non-stopping warnings. Inf values are not permitted.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Wan et al., 2008; Eq. 3, p. 1145)
%
%   YCr = Cr_sp / (Cr_sp + Al_sp)
%   KD  = Al2O3_ol / Al2O3_sp
%
%   T (degreeC) = 10000 / (0.512 + 0.873*YCr - 0.910*ln(KD)) - 273
%
% Notes:
% - ln is the natural logarithm.
% - P_kbar is accepted for interface compatibility and output traceability,
%   but it is not used in the Wan et al. (2008) thermometer equation.
% - When P_kbar is a vector, one output row is returned per pressure value.
%   Because the equation is pressure independent, the calculated temperature
%   is repeated for each pressure value.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Wang2008(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Olivine and Spinel tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Ol–Sp pair. The output variable set is intended
%             to remain stable for downstream processing.
%

%% Input validation
% Basic argument checks to prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('Wang2008 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. We do not modify the
% tables here; we only read the relevant columns during calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end
if ~isfield(rawdata_struct, 'Spinel') || ~istable(rawdata_struct.Spinel)
    error('rawdata_struct must contain table: rawdata_struct.Spinel');
end

dataset_ol = rawdata_struct.Olivine;
dataset_sp = rawdata_struct.Spinel;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated table concatenation inside the loop is avoided because it forces
% MATLAB to repeatedly reallocate and copy the entire results table.
%
% The cell buffer is preallocated and doubled only when necessary. After the
% interactive loop finishes, all blocks are concatenated once with vertcat.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct experimental calibration limits reported by Wan et al. (2008).
calibrationT_min_degC = 1250;
calibrationT_max_degC = 1450;

% The equation was calibrated at 100 kPa and tested experimentally at 1 GPa.
% The interval below is used only for a non-stopping pressure warning.
pressureScreen_min_GPa = 0;
pressureScreen_max_GPa = 1;

% Pressure is common to all selected mineral pairs in this function call.
% The warning is therefore printed only once, after the first calculation.
P_GPa_input = P_kbar ./ 10;
pressureOutsideScreen = ...
    P_GPa_input < pressureScreen_min_GPa | ...
    P_GPa_input > pressureScreen_max_GPa;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or selects
% "Finish" after a calculation.
disp('=== Step 3: Selecting a data code from the list (Olivine) ===');

while true
    % ----- Olivine selection -----
    % Assumption: the first column stores an identifier to display to the user.
    dataCodes_ol = dataset_ol{:, 1};

    [selectedIdx_ol, ok] = listdlg( ...
        'PromptString', 'Please select the Olivine data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_ol, ...
        'ListSize', [320 320]);

    % If the user cancels, exit the loop gracefully.
    if ~ok || isempty(selectedIdx_ol)
        disp('Selection canceled');
        break;
    end

    selectedCode_ol = dataCodes_ol(selectedIdx_ol);
    disp(['Olivine selected: ' char(string(selectedCode_ol))]);

    % ----- Spinel selection -----
    disp('=== Step 4: Selecting a data code from the list (Spinel) ===');

    dataCodes_sp = dataset_sp{:, 1};

    [selectedIdx_sp, ok] = listdlg( ...
        'PromptString', 'Please select the Spinel data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_sp, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_sp)
        disp('Selection canceled');
        break;
    end

    selectedCode_sp = dataCodes_sp(selectedIdx_sp);
    disp(['Spinel selected: ' char(string(selectedCode_sp))]);

    % ----- Calculation -----
    % IMPORTANT:
    % Olivine and spinel are selected independently; do not assume row indices
    % correspond between the two datasets.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_sp = dataset_sp(selectedIdx_sp, :);

    % Confirm that all variables used by the thermometer are present before
    % inspecting or calculating their values.
    validateRequiredVariables(selectedData_ol, selectedData_sp);

    % NaN values are retained and allowed to propagate. Zero values are also
    % retained, but both are collected for non-stopping warnings.
    nanInputNames = findNaNInputs(selectedData_ol, selectedData_sp);
    zeroInputNames = findZeroInputs(selectedData_ol, selectedData_sp);

    % Negative finite inputs and Inf values are prohibited. NaN and zero are
    % intentionally excluded from the error condition.
    validateNonNegativeInputs(selectedData_ol, selectedData_sp);

    row = calcTemp(selectedData_ol, selectedData_sp, P_kbar);

    % Store the user-selected identifiers for traceability.
    row.dataCode_ol = repmat(string(selectedCode_ol), height(row), 1);
    row.dataCode_sp = repmat(string(selectedCode_sp), height(row), 1);

    % Move identifiers to the front for readability.
    row = movevars(row, {'dataCode_ol','dataCode_sp'}, 'Before', 1);

    % Store this result as one table block. If the preallocated buffer is
    % full, double its capacity rather than resizing on every iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    finiteTemperature = isfinite(row.T_deg);
    if height(row) == 1
        disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    elseif any(finiteTemperature)
        finiteValues = row.T_deg(finiteTemperature);
        disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
            ': ' num2str(min(finiteValues)) ' to ' ...
            num2str(max(finiteValues)) ' degreeC']);
    else
        disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
            ': NaN to NaN degreeC']);
    end

    % Warn once when any input pressure lies outside the empirical pressure
    % screening interval. The calculation is not stopped.
    if any(pressureOutsideScreen) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the empirical pressure-screening ' ...
             'interval used for Wan et al. (2008): 0–1 GPa (0–10 kbar). ' ...
             'The thermometer was directly calibrated at 100 kPa and compared ' ...
             'with published experiments at 1 GPa; the equation contains no ' ...
             'pressure term. %d of %d pressure point(s) are outside the interval; ' ...
             'input range = %.4g–%.4g GPa.\n'], ...
            sum(pressureOutsideScreen), ...
            numel(P_GPa_input), ...
            min(P_GPa_input), ...
            max(P_GPa_input));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the direct
    % experimental calibration range of 1250–1450 degreeC.
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the direct experimental ' ...
             'calibration range of Wan et al. (2008): 1250–1450 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)));
    end

    % Composition-range warnings based on Wan et al. (2008).
    if any(~row.is_YCr_in_range)
        fprintf(2, ...
            ['WARNING: Spinel YCr is outside the recommended range of Wan et al. ' ...
             '(2008): 0.07–0.69 for %s & %s. Calculated YCr = %.6g.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            row.YCr(1));
    end

    if any(~row.is_Fe3_low)
        fprintf(2, ...
            ['WARNING: Spinel Fe3+ is outside the recommended range of Wan et al. ' ...
             '(2008): Fe3+ < 0.1 atoms per 4 O for %s & %s. ' ...
             'Input Fe3+ = %.6g atoms per 4 O.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            row.Fe3_sp(1));
    end

    if any(~row.is_Ti_low)
        fprintf(2, ...
            ['WARNING: Spinel Ti is outside the recommended range of Wan et al. ' ...
             '(2008): Ti < 0.025 atoms per 4 O for %s & %s. ' ...
             'Input Ti = %.6g atoms per 4 O.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            row.Ti_sp(1));
    end

    if any(~row.is_Fe3_fraction_valid)
        fprintf(2, ...
            ['WARNING: Fe3+/(Fe3+ + Cr + Al) exceeds 0.15 for %s & %s ' ...
             '(calculated value = %.6g). Wan et al. (2008, p. 1146) state ' ...
             'that the thermometer cannot be used for such Fe3+-rich spinel.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            row.Fe3_fraction(1));
    end

    % Print a non-stopping warning immediately after the temperature result
    % when any required thermometer input contains NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained as a missing value and the calculation ' ...
             'was continued; the calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Zero is not prohibited by the requested input rule, but it may make a
    % ratio or logarithm undefined. Retain the value and warn without stopping.
    if ~isempty(zeroInputNames)
        fprintf(2, ...
            ['WARNING: Zero was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         Zero was retained, but it may produce an undefined ratio, ' ...
             'logarithm, or physically invalid temperature.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            char(strjoin(zeroInputNames, ', ')));
    end

    % Warn about non-finite or mathematically invalid intermediate values.
    invalidIntermediate = ...
        ~isfinite(row.YCr) | ...
        ~isfinite(row.KD_Al_ol_sp) | ...
        row.KD_Al_ol_sp <= 0 | ...
        ~isfinite(row.lnKD) | ...
        ~isfinite(row.denom) | ...
        row.denom <= 0;

    if any(invalidIntermediate)
        fprintf(2, ...
            ['WARNING: One or more intermediate values are non-finite or outside ' ...
             'the mathematical domain for %s & %s (%d of %d pressure point(s)).\n' ...
             '         Values were retained in the output table and the ' ...
             'calculation was not stopped.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            sum(invalidIntermediate), ...
            numel(invalidIntermediate));
    end

    % Retain a result-based check for NaN/Inf temperatures.
    invalidTemperature = ~isfinite(row.T_deg);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Wang2008', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate the buffered table blocks only once after all selections have
% been completed. Return an empty table when no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(data_olivine, data_spinel)
% validateRequiredVariables
% Confirm that all variables required by the Wan et al. (2008) thermometer
% exist in the selected olivine and spinel tables.

olivineVariables = {'Al2O3'};
spinelVariables = {'Cr_cation_apfu', 'Al_cation_apfu', ...
    'Al2O3', 'Fe3_cation_apfu', 'Ti_cation_apfu'};

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    if ~ismember(variableName, data_olivine.Properties.VariableNames)
        error('Olivine table must contain variable: %s', variableName);
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    if ~ismember(variableName, data_spinel.Properties.VariableNames)
        error('Spinel table must contain variable: %s', variableName);
    end
end

end

function nanInputNames = findNaNInputs(data_olivine, data_spinel)
% findNaNInputs
% Return the names of required thermometer input variables that contain NaN.
% This function does not throw an error for NaN values.

olivineVariables = {'Al2O3'};
spinelVariables = {'Cr_cation_apfu', 'Al_cation_apfu', ...
    'Al2O3', 'Fe3_cation_apfu', 'Ti_cation_apfu'};

maxNames = numel(olivineVariables) + numel(spinelVariables);
nameBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isnan(variableValue(:)))
        nNames = nNames + 1;
        nameBuffer(nNames) = "Olivine." + string(variableName);
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isnan(variableValue(:)))
        nNames = nNames + 1;
        nameBuffer(nNames) = "Spinel." + string(variableName);
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function zeroInputNames = findZeroInputs(data_olivine, data_spinel)
% findZeroInputs
% Return the names of required thermometer input variables containing finite
% zero values. Zero is retained and reported without stopping calculation.

olivineVariables = {'Al2O3'};
spinelVariables = {'Cr_cation_apfu', 'Al_cation_apfu', ...
    'Al2O3', 'Fe3_cation_apfu', 'Ti_cation_apfu'};

maxNames = numel(olivineVariables) + numel(spinelVariables);
nameBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) == 0)
        nNames = nNames + 1;
        nameBuffer(nNames) = "Olivine." + string(variableName);
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) == 0)
        nNames = nNames + 1;
        nameBuffer(nNames) = "Spinel." + string(variableName);
    end
end

zeroInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_olivine, data_spinel)
% validateNonNegativeInputs
% Stop calculation when a required mineral-composition value is negative or
% infinite. NaN and zero are intentionally allowed so they can be retained
% and reported by non-stopping warning messages.

olivineVariables = {'Al2O3'};
spinelVariables = {'Cr_cation_apfu', 'Al_cation_apfu', ...
    'Al2O3', 'Fe3_cation_apfu', 'Ti_cation_apfu'};

maxNames = numel(olivineVariables) + numel(spinelVariables);
negativeBuffer = strings(maxNames, 1);
infiniteBuffer = strings(maxNames, 1);
nNegative = 0;
nInfinite = 0;

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    variableValue = data_olivine.(variableName);

    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nNegative = nNegative + 1;
        negativeBuffer(nNegative) = "Olivine." + string(variableName);
    end

    if any(isinf(variableValue(:)))
        nInfinite = nInfinite + 1;
        infiniteBuffer(nInfinite) = "Olivine." + string(variableName);
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    variableValue = data_spinel.(variableName);

    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nNegative = nNegative + 1;
        negativeBuffer(nNegative) = "Spinel." + string(variableName);
    end

    if any(isinf(variableValue(:)))
        nInfinite = nInfinite + 1;
        infiniteBuffer(nInfinite) = "Spinel." + string(variableName);
    end
end

negativeInputNames = negativeBuffer(1:nNegative);
infiniteInputNames = infiniteBuffer(1:nInfinite);

if ~isempty(negativeInputNames)
    error(['Wang2008: required mineral-composition values must be >= 0. ' ...
           'Negative finite value(s) were found in: ' ...
           char(strjoin(negativeInputNames, ', ')) '.']);
end

if ~isempty(infiniteInputNames)
    error(['Wang2008: required mineral-composition values must not be Inf. ' ...
           'Infinite value(s) were found in: ' ...
           char(strjoin(infiniteInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_olivine, data_spinel, P_kbar)
% calcTemp
% Compute temperature for one olivine row and one spinel row, returning one
% output row per supplied pressure value. Pressure is retained for interface
% compatibility and traceability but is not used in Equation 3.
%
% Inputs:
%   data_olivine : 1-row table containing olivine Al2O3 in wt%
%   data_spinel  : 1-row table containing spinel composition
%   P_kbar       : pressure in kbar; scalar or vector
%
% Output:
%   row : table containing pressure, intermediate variables, applicability
%         flags, and calculated temperatures.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();

% Store pressure for traceability and compatibility with fixed-pressure and
% pressure-range calculation modes.
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% --- Extract and replicate input values ---
% The selected mineral analyses are scalar values. Replication produces one
% output row per pressure value while preserving the pressure-independent
% temperature result.
Al2O3_ol = repmat(data_olivine.Al2O3, nP, 1);
Al2O3_sp = repmat(data_spinel.Al2O3, nP, 1);

Cr_sp  = repmat(data_spinel.Cr_cation_apfu, nP, 1);
Al_sp  = repmat(data_spinel.Al_cation_apfu, nP, 1);
Fe3_sp = repmat(data_spinel.Fe3_cation_apfu, nP, 1);
Ti_sp  = repmat(data_spinel.Ti_cation_apfu, nP, 1);

% --- Wan et al. (2008) thermometer variables ---
YCr = Cr_sp ./ (Cr_sp + Al_sp);
KD = Al2O3_ol ./ Al2O3_sp;
lnKD = log(KD);

% Eq. (3) of Wan et al. (2008).
denom = 0.512 + 0.873 .* YCr - 0.910 .* lnKD;

T_deg = 10000 ./ denom - 273;
T_K = T_deg + 273.15;

% Additional ferric-iron screening criterion discussed on p. 1146.
Fe3_fraction = Fe3_sp ./ (Fe3_sp + Cr_sp + Al_sp);

% --- Applicability flags based on the original paper ---
is_YCr_in_range = (YCr >= 0.07) & (YCr <= 0.69);
is_Fe3_low = Fe3_sp < 0.1;
is_Ti_low = Ti_sp < 0.025;
is_Fe3_fraction_valid = Fe3_fraction <= 0.15;
is_recommended = is_YCr_in_range & is_Fe3_low & is_Ti_low & ...
    is_Fe3_fraction_valid;

is_pressure_screened = P_GPa >= 0 & P_GPa <= 1;
is_temperature_calibrated = ...
    isfinite(T_deg) & T_deg >= 1250 & T_deg <= 1450;

% --- Pack outputs ---
row.Al2O3_ol = Al2O3_ol;
row.Al2O3_sp = Al2O3_sp;

row.Cr_sp = Cr_sp;
row.Al_sp = Al_sp;
row.Fe3_sp = Fe3_sp;
row.Ti_sp = Ti_sp;

row.YCr = YCr;
row.KD_Al_ol_sp = KD;
row.lnKD = lnKD;
row.denom = denom;
row.Fe3_fraction = Fe3_fraction;

row.is_YCr_in_range = is_YCr_in_range;
row.is_Fe3_low = is_Fe3_low;
row.is_Ti_low = is_Ti_low;
row.is_Fe3_fraction_valid = is_Fe3_fraction_valid;
row.is_recommended = is_recommended;
row.is_pressure_screened = is_pressure_screened;
row.is_temperature_calibrated = is_temperature_calibrated;

row.T_K = T_K;
row.T_deg = T_deg;

end
