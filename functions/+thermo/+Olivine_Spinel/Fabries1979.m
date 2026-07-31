function results = Fabries1979(rawdata_struct, P_kbar)
% functions/+thermo/+Olivine_Spinel/Fabries1979.m
% Tested with MATLAB R2024b
%
% Empirical Fe–Mg exchange thermometer between Olivine and Spinel
% Fabriès, J. (1979)
% Contributions to Mineralogy and Petrology, 69, 329–336
% DOI: https://doi.org/10.1007/BF00372258
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis and one Spinel
% analysis (selected by the user from tables) and calculates temperature
% using the empirical Fe–Mg exchange thermometer of Fabriès (1979).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Ol–Sp pair, and stores the results in a
% single output table.
%
% P_kbar may be supplied as either a scalar or a vector. Although pressure is
% not included in the Fabriès (1979) thermometer equation, one output row is
% returned for every pressure value so that this function can be used by both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Fabriès (1979) constructed the empirical calibration using two reference
% isotherms and interpolated intermediate isotherms in proportion to 1/T:
%
%   Temperature : approximately 700–1200 degreeC
%   Pressure    : no formal pressure calibration range was defined
%   Composition : principally spinel peridotites and lherzolitic rocks
%
% The 700 degreeC reference isotherm was obtained from 32 metaperidotite data
% points covering YCr_sp = 0.0–0.7 (p. 330). The 1200 degreeC reference
% isotherm was constrained using volcanic and experimentally equilibrated
% olivine–spinel pairs at approximately 1200–1230 degreeC and 0.5–1.6 GPa
% (pp. 330–331). Additional comparison experiments at approximately
% 850–920 degreeC and 1.0–2.0 GPa are reported in Table 2 on p. 332.
%
% Because Fabriès (1979) did not define a formal pressure applicability
% interval, this implementation uses 0.5–2.0 GPa only as an approximate
% pressure range represented by the calibration and comparison experiments.
% A pressure warning therefore indicates extrapolation beyond the published
% experimental pressure coverage, not violation of a pressure term in the
% thermometer equation.
%
% Important application notes from Fabriès (1979):
%   1) The calibration assumes spinel compositions typical of spinel
%      peridotites and generally low ferric-iron contents, commonly
%      YFe3_sp < 0.05 (p. 330).
%   2) Complete and accurate spinel analyses are required because Fe3+/Fe2+
%      is commonly estimated by stoichiometry. Analytical uncertainty may
%      introduce an uncertainty of approximately +/-50 degreeC (p. 332).
%   3) A calibration developed for a particular compositional range may give
%      erroneous results when applied indiscriminately to other mineral
%      compositions (p. 333).
%   4) Fe–Mg exchange between olivine and spinel may continue during cooling
%      to relatively low temperature. Calculated values may therefore record
%      late-stage exchange-blocking or closure temperatures rather than the
%      original formation or peak-equilibration temperature (pp. 333–335).
%
% This implementation issues non-stopping warnings using fprintf when:
%   1) input pressure is outside the approximate experimental coverage of
%      0.5–2.0 GPa, or
%   2) a finite calculated temperature is outside approximately
%      700–1200 degreeC.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%   rawdata_struct.Spinel  : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following variable names (normalized cations):
%
%   Olivine table variables:
%     Mg_cation_apfu
%     Fe_cation_apfu         % Fe2+ in olivine (assumed)
%
%   Spinel table variables:
%     Mg_cation_apfu
%     Fe_cation_apfu         % Fe2+ in spinel
%     Fe3_cation_apfu        % Fe3+ in spinel
%     Cr_cation_apfu
%     Al_cation_apfu
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Fabriès, 1979; implemented here)
%
% 1) Mole fractions / indices
%   Olivine:
%     XMg_ol = Mg_ol / (Mg_ol + Fe2_ol)
%     XFe_ol = Fe2_ol / (Mg_ol + Fe2_ol)
%
%   Spinel:
%     XMg_sp  = Mg_sp  / (Mg_sp + Fe2_sp)
%     XFe2_sp = Fe2_sp / (Mg_sp + Fe2_sp)
%
%     YCr_sp = Cr_sp / (Cr_sp + Al_sp + Fe3_sp)
%
% 2) Fe–Mg exchange coefficient (Ol–Sp)
%     KD = (XMg_ol * XFe2_sp) / (XFe_ol * XMg_sp)
%
% 3) Temperature solution
%     T(K) = (4250 * YCr_sp + 1343) ...
%            / (ln(KD) + 1.825 * YCr_sp + 0.571)
%
%     T(degreeC) = T(K) - 273.15
%
% Notes:
% - This is the empirical calibration given as equation (17) on p. 331 of
%   Fabriès (1979).
% - P_kbar is not used in the temperature equation. It is retained in the
%   function interface and output for compatibility with fixed-pressure and
%   pressure-range workflows.
% - Accurate Fe3+ estimation in spinel is important because YCr_sp includes
%   Fe3+ in its denominator.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Fabries1979(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Olivine and Spinel tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Ol–Sp pair. Because the equation is pressure
%             independent, temperatures are repeated for all supplied
%             pressure values.
%

%% Input validation
% Basic argument checks to prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('Fabries1979 requires (rawdata_struct, P_kbar).');
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

% Approximate calibration/application limits represented in Fabriès (1979).
calibrationT_min_degC = 700;
calibrationT_max_degC = 1200;
calibrationP_min_GPa = 0.5;
calibrationP_max_GPa = 2.0;

% Pressure is common to all selected mineral pairs in this function call.
% The warning is therefore printed only once, after the first calculation.
P_GPa_input = P_kbar ./ 10;
pressureOutsideCalibration = ...
    P_GPa_input < calibrationP_min_GPa | ...
    P_GPa_input > calibrationP_max_GPa;
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

    % Store the selected input rows so that the thermometer inputs can be
    % checked explicitly for NaN values without interrupting calculation.
    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_sp = dataset_sp(selectedIdx_sp, :);

    % Check only the variables that are actually used in the thermometer.
    % The calculation is intentionally allowed to continue even when NaN is
    % present; a warning message is printed after the temperature result.
    nanInputNames = findNaNInputs(selectedData_ol, selectedData_sp);

    row = calcTemp(selectedData_ol, selectedData_sp, P_kbar);

    % Store the user-selected identifiers for traceability.
    row.dataCode_ol = repmat(string(selectedCode_ol), height(row), 1);
    row.dataCode_sp = repmat(string(selectedCode_sp), height(row), 1);

    % Move identifiers to the front for readability.
    row = movevars(row, {'dataCode_ol','dataCode_sp'}, 'Before', 1);

    % Store this result as one table block. If the preallocated buffer is
    % full, double its capacity. This avoids changing the output-table size
    % during every iteration of the interactive loop.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        oldCapacity = numel(resultBlocks);
        resultBlocks(oldCapacity + 1 : 2 * oldCapacity, 1) = cell(oldCapacity, 1);
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
            ': ' num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside the approximate
    % experimental pressure coverage of 0.5–2.0 GPa. Pressure is not part of
    % the thermometer equation, so this warning does not stop calculation.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate experimental ' ...
             'pressure coverage represented in Fabriès (1979): 0.5–2.0 GPa ' ...
             '(5–20 kbar). Fabriès (1979) did not define a formal pressure ' ...
             'calibration range, and pressure is not used in the equation. ' ...
             '%d of %d pressure point(s) are outside this approximate range; ' ...
             'input range = %.4g–%.4g GPa.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_GPa_input), ...
            min(P_GPa_input), ...
            max(P_GPa_input));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the empirical
    % calibration interval of approximately 700–1200 degreeC. NaN and Inf are
    % handled separately by the non-finite-result warning below.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximate empirical ' ...
             'calibration range of Fabriès (1979): 700–1200 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)));
    end

    % Print a non-stopping warning immediately after the temperature result
    % when any required thermometer input contains NaN. fprintf is used
    % instead of warning so that the message remains visible even when MATLAB
    % warnings have been disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, but the calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.T_deg);

    % Also retain a result-based check for NaN/Inf values caused by reasons
    % other than an explicitly NaN input, such as division by zero or a
    % non-positive exchange coefficient.
    if any(invalidTemperature)
        fprintf(2, ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
                    '(%d of %d points; NaN: %d, Inf: %d).\n' ...
                    '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
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
        'Fabries1979', ...
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
function nanInputNames = findNaNInputs(data_olivine, data_spinel)
% findNaNInputs
% Return the names of required thermometer input variables that contain NaN.
% This function does not throw an error for NaN values; it only prepares a
% warning message for the calling function.

olivineVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
spinelVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', ...
    'Fe3_cation_apfu', 'Cr_cation_apfu', 'Al_cation_apfu'};

% The maximum number of possible entries is known in advance. Preallocation
% avoids increasing the string-array size during each loop iteration.
maxInputCount = numel(olivineVariables) + numel(spinelVariables);
nanInputBuffer = strings(maxInputCount, 1);
nNanInputs = 0;

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Olivine." + string(variableName);
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Spinel." + string(variableName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function row = calcTemp(data_olivine, data_spinel, P_kbar)
% calcTemp
% Compute one Fabriès (1979) temperature estimate for a selected olivine–
% spinel pair and repeat that estimate for every supplied pressure value.
%
% Inputs:
%   data_olivine : 1-row table (olivine apfu cations)
%   data_spinel  : 1-row table (spinel apfu cations)
%   P_kbar       : pressure in kbar; retained for workflow compatibility
%
% Output:
%   row : table containing one row per pressure value, including composition
%         indices, exchange coefficient, equation terms, and temperature.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();

% Store pressure in both kbar and GPa for traceability. Pressure does not
% enter the Fabriès (1979) temperature equation.
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% --- Extract cation data (expected variable names) ---
% The selected analyses are scalar values. Replicate each scalar to the
% pressure-vector length so every table variable has a consistent row count.
Mg_ol  = repmat(data_olivine.Mg_cation_apfu, nP, 1);
Fe2_ol = repmat(data_olivine.Fe_cation_apfu, nP, 1);

Mg_sp  = repmat(data_spinel.Mg_cation_apfu, nP, 1);
Fe2_sp = repmat(data_spinel.Fe_cation_apfu, nP, 1);
Fe3_sp = repmat(data_spinel.Fe3_cation_apfu, nP, 1);
Cr_sp  = repmat(data_spinel.Cr_cation_apfu, nP, 1);
Al_sp  = repmat(data_spinel.Al_cation_apfu, nP, 1);

% --- Compute composition indices used in Fabriès (1979) ---
% Olivine Mg and Fe2+ fractions on the (Mg + Fe2+) basis.
XMg_ol = Mg_ol ./ (Mg_ol + Fe2_ol);
XFe_ol = Fe2_ol ./ (Mg_ol + Fe2_ol);

% Spinel Mg and Fe2+ fractions on the (Mg + Fe2+) basis.
XMg_sp  = Mg_sp ./ (Mg_sp + Fe2_sp);
XFe2_sp = Fe2_sp ./ (Mg_sp + Fe2_sp);

% Chromium fraction among trivalent spinel cations.
YCr_sp = Cr_sp ./ (Cr_sp + Al_sp + Fe3_sp);

% --- Fe–Mg exchange coefficient between olivine and spinel ---
KD_ol_sp = (XMg_ol .* XFe2_sp) ./ (XFe_ol .* XMg_sp);

% --- Solve temperature (Fabriès, 1979, equation 17, p. 331) ---
numerator = 4250 .* YCr_sp + 1343;
denominator = log(KD_ol_sp) + 1.825 .* YCr_sp + 0.571;

T_K   = numerator ./ denominator;
T_deg = T_K - 273.15;

% --- Pack outputs ---
% Store intermediate indices and equation terms for transparency and
% downstream checking.
row.XMg_ol = XMg_ol;
row.XFe_ol = XFe_ol;

row.XMg_sp = XMg_sp;
row.XFe2_sp = XFe2_sp;
row.Fe3_sp_cation = Fe3_sp;

row.YCr_sp = YCr_sp;
row.KD_ol_sp = KD_ol_sp;

row.Numerator = numerator;
row.Denominator = denominator;

row.T_K = T_K;
row.T_deg = T_deg;

end
