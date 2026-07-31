function results = Roeder1979(rawdata_struct, P_kbar)
% functions/+thermo/+Olivine_Spinel/Roeder1979.m
% Tested with MATLAB R2024b
%
% Thermodynamic Fe–Mg exchange thermometer between Olivine and Spinel
% Roeder, P.L., Campbell, I.H., Jamieson, H.E. (1979)
% Contributions to Mineralogy and Petrology, 68, 325–334
% DOI: https://doi.org/10.1007/BF00371554
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis and one Spinel
% analysis (selected by the user from tables) and calculates temperature
% using the Fe2+–Mg exchange equilibrium between olivine and spinel,
% following Roeder et al. (1979), Eq. (3).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Ol–Sp pair, and appends results into a
% single output table.
%
% -------------------------------------------------------------------------
% APPLICATION RANGE AND IMPORTANT LIMITATIONS
%
% Roeder et al. (1979) did not establish a formal experimental calibration
% range in temperature and pressure comparable to later experimentally
% calibrated thermometers. The revised equation was constructed mainly from
% selected thermodynamic data and was evaluated using natural olivine–spinel
% pairs and heating experiments on uncrushed natural samples.
%
% The following calculated temperature ranges are reported in the abstract
% on p. 325:
%
%   Volcanic samples from Kilauea      : 1100–1300 degreeC
%   Basic rocks from layered intrusions:  500–800 degreeC
%
% The 500–800 degreeC values were interpreted as subsolidus Fe2+–Mg
% re-equilibration or closure temperatures in slowly cooled intrusions, not
% as primary magmatic crystallization temperatures (abstract, p. 325).
% Heating experiments on natural samples were conducted at magmatic
% temperatures for periods of two days to four weeks (abstract, p. 325).
%
% Equation (3) is presented on p. 333. It returns temperature in Kelvin and
% contains no pressure term. Roeder et al. (1979) did not define an explicit
% pressure calibration or application range. P_kbar is therefore retained
% only for interface compatibility and output traceability; it does not
% affect the calculated temperature.
%
% Roeder et al. (1979) also note that published free-energy data for spinel
% end-members vary substantially and that temperature formulations based on
% those data remain uncertain until the thermodynamic data are improved
% (p. 333).
%
% Engi and Evans (1980, Contributions to Mineralogy and Petrology, 73,
% 201–203; DOI: https://doi.org/10.1007/BF00371395) subsequently questioned
% the experimental reversibility and validity of Eq. (3), and argued that
% its calculated temperatures may not be meaningful. The thermometer should
% therefore be treated as a historical or comparative calibration rather
% than as a stand-alone temperature constraint.
%
% Because no formal calibration interval was defined, this implementation
% uses the broad 500–1300 degreeC envelope of the natural-sample results on
% p. 325 only as an APPLICATION-RESULT screening range. This is not a direct
% experimental calibration range.
%
% This implementation therefore issues non-stopping fprintf messages when:
%   1) pressure is supplied, because no published pressure range exists and
%      pressure is not used in Eq. (3),
%   2) a finite calculated temperature is outside 500–1300 degreeC,
%   3) a required input contains NaN, or
%   4) a calculated temperature is NaN or Inf.
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
% All finite values in the required mineral-composition variables above must
% be strictly greater than zero. NaN values are retained as missing values,
% propagated through the calculation, and reported by non-stopping fprintf
% messages.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Roeder et al., 1979, Eq. 3, p. 333)
%
% 1) Mole fractions / ratios
%   Olivine:
%     XMg_ol = Mg_ol / (Mg_ol + Fe2_ol)
%     XFe_ol = Fe2_ol / (Mg_ol + Fe2_ol)
%
%   Spinel:
%     XMg_sp  = Mg_sp  / (Mg_sp + Fe2_sp)
%     XFe2_sp = Fe2_sp / (Mg_sp + Fe2_sp)
%
% 2) Distribution coefficient
%     Kd = (XMg_ol / XFe_ol) * (XFe2_sp / XMg_sp)
%
% 3) Trivalent cation fractions in spinel
%     alpha = Cr  / (Cr + Al + Fe3)
%     beta  = Al  / (Cr + Al + Fe3)
%     gamma = Fe3 / (Cr + Al + Fe3)
%
% 4) Temperature equation
%     T(K) =
%       (3480*alpha + 1018*beta - 1720*gamma + 2400) / ...
%       (2.23*alpha + 2.56*beta - 3.08*gamma - 1.47 ...
%        + 1.987*ln(Kd))
%
%     T(degreeC) = T(K) - 273.15
%
% Notes:
% - Equation (3) is printed as t(K) in Roeder et al. (1979, p. 333).
% - Pressure does not occur in the equation. When P_kbar is a vector, the
%   same pressure-independent temperature is repeated for each pressure
%   value so that both startThermoCalc_fixedP and startThermoCalc_rangeP can
%   use a consistent output structure.
% - The result is sensitive to the Fe2+/Fe3+ allocation in spinel and should
%   only be interpreted for texturally and chemically equilibrated Ol–Sp
%   pairs representing the same equilibration stage.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Roeder1979(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Olivine and Spinel tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector;
%                    retained for interface compatibility but not used in
%                    the Roeder et al. temperature equation)
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
    error('Roeder1979 requires (rawdata_struct, P_kbar).');
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

% Roeder et al. (1979) did not define a formal calibration range. These
% limits represent only the broad envelope of the natural-sample temperature
% results reported in the abstract on p. 325.
applicationT_min_degC = 500;
applicationT_max_degC = 1300;

% Pressure is common to all selected mineral pairs in this function call.
% Because Eq. (3) has no pressure term and no pressure range was defined,
% this limitation message is printed only once after the first calculation.
pressureLimitationMessageIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or selects
% 'Finish' after a calculation.
disp('=== Step 3: Selecting a data code from the list (Olivine) ===');

while true
    % ----- Olivine selection -----
    % Assumption: the first column stores an identifier to display to the user.
    dataCodes_ol = dataset_ol{:, 1};

    [selectedIdx_ol, ok] = listdlg( ...
        'PromptString', 'Please select the Olivine data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_ol)), ...
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
        'ListString', cellstr(string(dataCodes_sp)), ...
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

    % Store the selected input rows so that thermometer inputs can be checked
    % explicitly without interrupting calculation for NaN values.
    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_sp = dataset_sp(selectedIdx_sp, :);

    % Check only the variables that are actually used in the thermometer.
    % Calculation is intentionally allowed to continue when NaN is present.
    nanInputNames = findNaNInputs(selectedData_ol, selectedData_sp);

    % All finite mineral-composition inputs used by this thermometer must be
    % strictly positive. NaN values are deliberately excluded from this check
    % so that they remain NaN and propagate through the calculation.
    validatePositiveInputs(selectedData_ol, selectedData_sp);

    row = calcTemp(selectedData_ol, selectedData_sp, P_kbar);

    % Store the user-selected identifiers for traceability.
    row.dataCode_ol = repmat(string(selectedCode_ol), height(row), 1);
    row.dataCode_sp = repmat(string(selectedCode_sp), height(row), 1);

    % Move identifiers to the front for readability.
    row = movevars(row, {'dataCode_ol','dataCode_sp'}, 'Before', 1);

    % Store this result as one table block. If the preallocated buffer is
    % full, double its capacity. This avoids resizing the results table on
    % every loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
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

    % Roeder et al. (1979) did not define a pressure calibration/application
    % range, and pressure is absent from Eq. (3). Therefore a scientifically
    % defensible pressure-outside-range test cannot be performed. Print this
    % limitation once without stopping the calculation.
    if ~pressureLimitationMessageIssued
        fprintf(2, ...
            ['WARNING: Roeder et al. (1979) did not define an explicit pressure ' ...
             'calibration or application range, and pressure is not used in Eq. (3). ' ...
             'Pressure-range screening cannot therefore be performed. ' ...
             'The supplied pressure value(s) are retained only for interface ' ...
             'compatibility and output traceability; input range = %.4g–%.4g kbar.\n'], ...
            min(P_kbar), ...
            max(P_kbar));
        pressureLimitationMessageIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the broad
    % 500–1300 degreeC envelope of natural-sample results reported on p. 325.
    % NaN and Inf are handled separately below.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideApplication = finiteTemperature & ...
        (row.T_deg < applicationT_min_degC | ...
         row.T_deg > applicationT_max_degC);

    if any(temperatureOutsideApplication)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the broad published ' ...
             'application-result envelope of Roeder et al. (1979): ' ...
             '500–1300 degreeC (abstract, p. 325). This interval is not a ' ...
             'formal experimental calibration range. %d of %d finite ' ...
             'temperature point(s) are outside the envelope; calculated ' ...
             'finite range = %.4g–%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideApplication), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)));
    end

    % Print a non-stopping message immediately after the temperature result
    % when any required thermometer input contains NaN. fprintf is used so
    % the message remains visible even when MATLAB warnings are disabled.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, NaN was not replaced by zero, and the calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.T_deg);

    % Retain a result-based check for NaN/Inf values caused by NaN inputs or
    % mathematical singularities. The values remain in the output table.
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
        'Roeder1979', ...
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

maxInputCount = numel(olivineVariables) + numel(spinelVariables);
nanInputNames = strings(maxInputCount, 1);
nNaNInputs = 0;

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isnan(variableValue(:)))
        nNaNInputs = nNaNInputs + 1;
        nanInputNames(nNaNInputs) = "Olivine." + string(variableName);
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isnan(variableValue(:)))
        nNaNInputs = nNaNInputs + 1;
        nanInputNames(nNaNInputs) = "Spinel." + string(variableName);
    end
end

nanInputNames = nanInputNames(1:nNaNInputs);

end

function validatePositiveInputs(data_olivine, data_spinel)
% validatePositiveInputs
% Stop the calculation when a finite required mineral-composition value is
% zero or negative. NaN is intentionally allowed so that it propagates
% through the thermometer calculation and is reported by fprintf messages.

olivineVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
spinelVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', ...
    'Fe3_cation_apfu', 'Cr_cation_apfu', 'Al_cation_apfu'};

maxInputCount = numel(olivineVariables) + numel(spinelVariables);
invalidInputNames = strings(maxInputCount, 1);
nInvalidInputs = 0;

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) <= 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputNames(nInvalidInputs) = ...
            "Olivine." + string(variableName);
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) <= 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputNames(nInvalidInputs) = ...
            "Spinel." + string(variableName);
    end
end

invalidInputNames = invalidInputNames(1:nInvalidInputs);

if ~isempty(invalidInputNames)
    error(['Roeder1979: required mineral-composition values must be > 0. ' ...
           'Zero or negative finite value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_olivine, data_spinel, P_kbar)
% calcTemp
% Compute a pressure-independent temperature estimate for one olivine row
% and one spinel row. When P_kbar is a vector, the same temperature and
% composition indices are repeated for every supplied pressure value.
%
% Inputs:
%   data_olivine : 1-row table (olivine apfu cations)
%   data_spinel  : 1-row table (spinel apfu cations)
%   P_kbar       : pressure in kbar; retained but not used in Eq. (3)
%
% Output:
%   row : table containing one row per supplied pressure value, including
%         intermediate variables and final temperatures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();

% Store pressure only for interface transparency and output traceability.
row.P_kbar = P_kbar;

% --- Extract cation data (expected variable names) ---
% Replicate the selected composition for each supplied pressure value so the
% output shape is compatible with both fixed-pressure and pressure-range
% calculation modes.
Mg_ol  = repmat(data_olivine.Mg_cation_apfu, nP, 1);
Fe2_ol = repmat(data_olivine.Fe_cation_apfu, nP, 1);

Mg_sp  = repmat(data_spinel.Mg_cation_apfu, nP, 1);
Fe2_sp = repmat(data_spinel.Fe_cation_apfu, nP, 1);
Fe3_sp = repmat(data_spinel.Fe3_cation_apfu, nP, 1);
Cr_sp  = repmat(data_spinel.Cr_cation_apfu, nP, 1);
Al_sp  = repmat(data_spinel.Al_cation_apfu, nP, 1);

% --- Compute composition indices used in Roeder et al. (1979) ---
% Olivine Mg and Fe2+ fractions on the (Mg + Fe2+) basis.
den_ol = Mg_ol + Fe2_ol;
XMg_ol = Mg_ol ./ den_ol;
XFe_ol = Fe2_ol ./ den_ol;

% Spinel Mg and Fe2+ fractions on the (Mg + Fe2+) basis.
den_sp_divalent = Mg_sp + Fe2_sp;
XMg_sp  = Mg_sp ./ den_sp_divalent;
XFe2_sp = Fe2_sp ./ den_sp_divalent;

% Spinel Cr, Al, and Fe3+ fractions among trivalent cations.
den_sp_trivalent = Cr_sp + Al_sp + Fe3_sp;
alpha = Cr_sp ./ den_sp_trivalent;
beta  = Al_sp ./ den_sp_trivalent;
gamma = Fe3_sp ./ den_sp_trivalent;

% --- Fe2+–Mg exchange coefficient between olivine and spinel ---
Kd = (XMg_ol ./ XFe_ol) .* (XFe2_sp ./ XMg_sp);
lnKd = log(Kd);

% --- Solve temperature (Roeder et al., 1979, Eq. 3, p. 333) ---
numerator = 3480 .* alpha + 1018 .* beta - 1720 .* gamma + 2400;
denominator = 2.23 .* alpha + 2.56 .* beta - 3.08 .* gamma ...
    - 1.47 + 1.987 .* lnKd;

% Equation (3) returns Kelvin. NaN or Inf values caused by missing inputs or
% mathematical singularities are intentionally retained for later reporting.
T_K = numerator ./ denominator;
T_deg = T_K - 273.15;

% --- Pack outputs ---
% Store intermediate indices and final results so the user can inspect the
% calculation and propagate uncertainties if needed.
row.XMg_ol = XMg_ol;
row.XFe_ol = XFe_ol;

row.XMg_sp  = XMg_sp;
row.XFe2_sp = XFe2_sp;

row.alpha = alpha;
row.beta  = beta;
row.gamma = gamma;

row.Kd = Kd;
row.lnKd = lnKd;

row.numerator = numerator;
row.denominator = denominator;

row.T_K = T_K;
row.T_deg = T_deg;

end
