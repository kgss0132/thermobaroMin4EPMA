function results = Coogan2014(rawdata_struct, P_kbar)
% functions/+thermo/+Olivine_Spinel/Coogan2014.m
% Tested with MATLAB R2024b
%
% Aluminum-in-olivine thermometer between Olivine and Spinel
% Coogan, L.A., Saunders, A.D., Wilson, R.N. (2014)
% Chemical Geology, 368, 1–10
% DOI: https://doi.org/10.1016/j.chemgeo.2014.01.004
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis and one Spinel
% analysis (selected by the user from tables) and calculates temperature
% using the Al-in-olivine thermometer of Coogan et al. (2014).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Ol–Sp pair, and appends results into a
% single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Coogan et al. (2014) updated the experimental Al-in-olivine thermometer
% originally calibrated by Wan et al. (2008). The calibration and tests
% reported by Coogan et al. (2014) cover approximately:
%
%   Temperature : 1250–1450 degreeC
%   Pressure    : 0.1 MPa (approximately 1 bar = 0.001 kbar)
%   Spinel Cr#  : 0–0.69, where Cr# = Cr / (Cr + Al)
%   Spinel Fe3+ : 0–0.11 atoms per 4 oxygens
%
% The original Wan et al. (2008) calibration range and the motivation for
% the additional experiments are summarized on p. 2. The new experiments
% were performed at 0.1 MPa, as described on p. 3. The updated equation,
% expanded compositional calibration range, stated precision, and warning
% against substantial extrapolation are given on p. 4.
%
% The 45 calibration experiments are reproduced within +/-25 degreeC, with
% all but three reproduced within +/-20 degreeC (p. 4). Coogan et al. (2014)
% explicitly state that use substantially outside the calibration range
% should be avoided (p. 4).
%
% IMPORTANT APPLICATION NOTES:
% - The thermometer records the peak temperature of olivine–Cr-spinel
%   co-saturation. It does not directly determine the liquidus temperature,
%   primary melt temperature, or mantle potential temperature (p. 2).
% - In primitive basalts, olivine and spinel generally crystallize early and
%   co-saturation may occur within a few tens of degrees of the liquidus.
%   Archean komatiites can be an important exception because spinel may
%   saturate only after substantial olivine crystallization (p. 2).
% - Use texturally and chemically appropriate equilibrium olivine–spinel
%   pairs. Zoning and mismatched crystal generations can change the inferred
%   temperature; natural-sample zoning is discussed on pp. 5–6.
% - If olivine is P-rich and Al correlates with P, coupled P–Al substitution
%   can raise olivine Al and overestimate temperature. Coogan et al. (2014)
%   recommend extrapolating olivine Al to zero P before calculation (pp. 3–4).
% - Accurate low-level Al analysis in olivine is essential. The authors used
%   long counting times, high beam current, repeated analyses, and careful
%   Z-focus adjustment; poor focus can lower measured Al (p. 2).
% - The effects of silica activity and spinel Fe3+ were small within the
%   tested ranges, generally comparable to the thermometer uncertainty
%   (pp. 3–4). Extrapolation to markedly different compositions or redox
%   states should nevertheless be treated cautiously.
% - This calibration was developed for olivine–Cr-spinel pairs from
%   primitive basaltic to ultramafic magmas. Direct application to
%   high-pressure mantle peridotites or strongly re-equilibrated assemblages
%   lies outside the experimental pressure calibration.
%
% This implementation therefore issues non-stopping warnings when:
%   1) input pressure differs from the 0.1 MPa calibration pressure, or
%   2) a finite calculated temperature is outside 1250–1450 degreeC.
%
% Pressure is not included in the published thermometer equation. P_kbar is
% retained only for compatibility with the common thermometer interface and
% for checking whether the requested pressure matches the experimental
% calibration condition.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%   rawdata_struct.Spinel  : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following variable names:
%
%   Olivine table variables:
%     Al2O3                   % wt.% Al2O3 in olivine
%
%   Spinel table variables:
%     Al2O3                   % wt.% Al2O3 in spinel
%     Cr_cation_apfu
%     Al_cation_apfu
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Coogan et al., 2014; p. 4)
%
%   T(K) = 10000 / (0.575 + 0.884 * Cr# - 0.897 * ln(kd))
%
% where
%
%   kd  = Al2O3_ol / Al2O3_sp
%   Cr# = Cr_sp / (Cr_sp + Al_sp)
%
% Notes:
% - Al2O3 values are used in wt.% as written in the paper.
% - Cr# is calculated from apfu spinel cations.
% - This thermometer does not contain a pressure term.
% - P_kbar is accepted only to retain function-call compatibility with other
%   thermometer functions in the same framework and to issue a calibration
%   pressure warning when appropriate.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Coogan2014(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Olivine and Spinel tables (see above)
%   P_kbar         : pressure in kbar; finite, non-negative scalar or vector
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
    error('Coogan2014 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve datasets
% Extract the required tables from the input struct. We do not modify the
% tables here; we only read the relevant columns during calculation.
disp('=== Step 1: Preparing dataset ===');

if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end
if ~isfield(rawdata_struct, 'Spinel') || ~istable(rawdata_struct.Spinel)
    error('rawdata_struct must contain table: rawdata_struct.Spinel');
end

dataset_ol = rawdata_struct.Olivine;
dataset_sp = rawdata_struct.Spinel;

disp('=== Preparing dataset has been finished ===');

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

% Experimental calibration limits reported by Coogan et al. (2014).
calibrationT_min_degC = 1250;
calibrationT_max_degC = 1450;
calibrationP_kbar = 0.001;  % 0.1 MPa = 1 bar = 0.001 kbar

% Pressure is common to all selected mineral pairs in this function call.
% Because the thermometer equation has no pressure term, this check only
% identifies extrapolation away from the experimental calibration pressure.
pressureTolerance_kbar = max(1e-12, 100 * eps(calibrationP_kbar));
pressureOutsideCalibration = ...
    abs(P_kbar - calibrationP_kbar) > pressureTolerance_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
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

    % Store this result as one table block.
    % If the preallocated buffer is full, double its capacity. This occurs
    % only occasionally and avoids reallocating the full results table after
    % every calculation.
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

    % Warn once when the input pressure differs from the experimental
    % calibration pressure of 0.1 MPa. The calculation is not stopped.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the direct experimental calibration ' ...
             'condition of Coogan et al. (2014): 0.1 MPa ' ...
             '(1 bar; 0.001 kbar). %d of %d pressure point(s) differ from ' ...
             'the calibration pressure; input range = %.6g–%.6g kbar.\n' ...
             '         Pressure is not used in the thermometer equation, so the ' ...
             'calculated temperature is pressure-independent and should be treated ' ...
             'as an extrapolation in pressure.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when the finite calculated temperature lies outside the
    % experimental calibration range of 1250–1450 degreeC. The calculation
    % is not stopped and the value remains in the output table.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the experimental calibration ' ...
             'range of Coogan et al. (2014): 1250–1450 degreeC. ' ...
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

    % Retain a result-based check for NaN/Inf values caused by NaN inputs or
    % other numerical reasons, such as division by zero.
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
        'Coogan2014', ...
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

olivineVariables = {'Al2O3'};
spinelVariables = {'Al2O3', 'Cr_cation_apfu', 'Al_cation_apfu'};

nanInputNames = strings(0, 1);

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isnan(variableValue(:)))
        nanInputNames(end + 1, 1) = "Olivine." + string(variableName); %#ok<AGROW>
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isnan(variableValue(:)))
        nanInputNames(end + 1, 1) = "Spinel." + string(variableName); %#ok<AGROW>
    end
end

end

function row = calcTemp(data_olivine, data_spinel, P_kbar)
% calcTemp
% Compute temperature estimates for one olivine row and one spinel row
% using Coogan et al. (2014). One output row is returned for each supplied
% pressure value, although pressure is not used in the thermometer equation.
%
% Inputs:
%   data_olivine : 1-row table containing olivine Al2O3 in wt.%
%   data_spinel  : 1-row table containing spinel Al2O3 and Cr, Al apfu
%   P_kbar       : pressure in kbar; scalar or vector, stored but not used
%
% Output:
%   row : table with one row per pressure value, containing kd, Cr#,
%         denominator, T(K), and T(degreeC)

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();

% Store input pressure for traceability, although it is not used in the
% published thermometer equation.
row.P_kbar = P_kbar;

% --- Validate required variable names ---
if ~ismember('Al2O3', data_olivine.Properties.VariableNames)
    error('Olivine table must contain variable: Al2O3');
end

requiredSpinelVars = {'Al2O3', 'Cr_cation_apfu', 'Al_cation_apfu'};
for i = 1:numel(requiredSpinelVars)
    if ~ismember(requiredSpinelVars{i}, data_spinel.Properties.VariableNames)
        error(['Spinel table must contain variable: ' requiredSpinelVars{i}]);
    end
end

% --- Extract required data ---
Al2O3_ol = repmat(data_olivine.Al2O3, nP, 1);
Al2O3_sp = repmat(data_spinel.Al2O3, nP, 1);
Cr_sp    = repmat(data_spinel.Cr_cation_apfu, nP, 1);
Al_sp    = repmat(data_spinel.Al_cation_apfu, nP, 1);

% --- Basic validity checks ---
% NaN values are intentionally not replaced by zero and do not stop the
% calculation. MATLAB arithmetic propagates NaN through kd, ln(kd), the
% denominator, and temperature, preserving the missing-value state.
% Finite values that violate the mathematical requirements still stop the
% calculation as before.
if any(Al2O3_ol <= 0)
    error('Olivine Al2O3 must be > 0 for Coogan2014 thermometer.');
end
if any(Al2O3_sp <= 0)
    error('Spinel Al2O3 must be > 0 for Coogan2014 thermometer.');
end
if any((Cr_sp + Al_sp) <= 0)
    error('Spinel Cr + Al must be > 0 to calculate Cr#.');
end

% --- Calculate thermometer parameters ---
% Al partition coefficient between olivine and spinel.
kd = Al2O3_ol ./ Al2O3_sp;
lnkd = log(kd);

% Spinel Cr# on the Cr + Al basis used by Coogan et al. (2014).
Cr_plus_Al_sp = Cr_sp + Al_sp;
Cr_num_sp = Cr_sp ./ Cr_plus_Al_sp;

% Denominator of the empirical temperature expression.
denominator = 0.575 + 0.884 .* Cr_num_sp - 0.897 .* lnkd;

if any(denominator <= 0)
    error('Calculated denominator is <= 0. Check input compositions.');
end

% Temperature in Kelvin and Celsius.
T_K = 10000 ./ denominator;
T_deg = T_K - 273.15;

% --- Pack outputs ---
% Store thermometer inputs and intermediate values for transparency and
% reproducibility.
row.Al2O3_ol = Al2O3_ol;
row.Al2O3_sp = Al2O3_sp;

row.Cr_sp = Cr_sp;
row.Al_sp = Al_sp;
row.Cr_num_sp = Cr_num_sp;

row.kd = kd;
row.lnkd = lnkd;

row.denominator = denominator;
row.T_K = T_K;
row.T_deg = T_deg;

end
