function results = Ballhaus1991(rawdata_struct, P_kbar)
% functions/+thermo/+Olivine_Spinel/Ballhaus1991.m
% Tested with MATLAB R2024b
%
% Thermodynamic Fe–Mg exchange thermometer between Olivine and Spinel
% Ballhaus, C., Berry, R.F., Green, D.H. (1991)
% Contribution to Mineralogy and Petrology, 107, 27–40
% DOI: https://doi.org/10.1007/BF00311183
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis and one Spinel
% analysis (selected by the user from tables) and calculates temperature
% using the Fe–Mg exchange equilibrium between olivine and spinel.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Ol–Sp pair, and appends results into a
% single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Ballhaus et al. (1991) calibrated the underlying olivine–spinel Fe–Mg
% exchange relations using synthetic spinel harzburgite and lherzolite
% assemblages over the following experimental range:
%
%   Temperature : 1040–1300 degreeC
%   Pressure    : 0.3–2.7 GPa
%   Composition : spinel harzburgite to spinel lherzolite assemblages
%   Redox range : a broad, controlled oxygen-fugacity range using the
%                 IW, WCO, NNO, and MH buffers
%
% These experimental conditions are summarized in the abstract on p. 27
% and described in the Experimental technique section on pp. 28–29.
% The olivine–spinel Fe–Mg exchange thermometer equation is presented on
% p. 35.
%
% Ballhaus et al. (1991) also state that their OXYGEN BAROMETER gives
% reasonable results down to approximately 800 degreeC (abstract, p. 27;
% discussion on p. 35). This approximately 800 degreeC statement refers to
% the oxygen barometer and does not extend the direct experimental
% calibration range of the olivine–spinel thermometer below 1040 degreeC.
%
% This implementation therefore issues non-stopping warnings when:
%   1) input pressure is outside 0.3–2.7 GPa, or
%   2) a finite calculated temperature is outside 1040–1300 degreeC.
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
%     Ti_cation_apfu
%     Cr_cation_apfu
%     Al_cation_apfu
%
% All finite values in the required mineral-composition variables above must
% be strictly greater than zero. NaN values are retained as missing values,
% propagated through the calculation, and reported by non-stopping warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% 1) Mole fractions / indices
%   Olivine:
%     XMg_ol = Mg_ol / (Mg_ol + Fe2_ol)
%     XFe_ol = Fe2_ol / (Mg_ol + Fe2_ol)
%
%   Spinel:
%     Mgnum_sp = Mg_sp / (Mg_sp + Fe2_sp)
%     XFe2_sp  = Fe2_sp / (Mg_sp + Fe2_sp)
%     XFe3_sp  = Fe3_sp / (Fe2_sp + Fe3_sp)
%     XTi_sp   = Ti_sp                       (as stored in apfu cations)
%     XCr_sp   = Cr_sp / (Cr_sp + Al_sp + Fe3_sp)
%
% 2) Fe–Mg exchange coefficient (Ol–Sp)
%     K_MgFe(ol-sp) = (XMg_ol * XFe2_sp) / (XFe_ol * Mgnum_sp)
%
% 3) Temperature solution
%     TA = R * ln(K_MgFe) + 4.705
%
%     P_GPa = P_kbar / 10
%
%     TB = (6530 + 280*P_GPa + 7000 + 108*P_GPa) * (1 - 2*XFe_ol)
%          - 1960  * (Mgnum_sp - XFe2_sp)
%          + 16150 * XCr_sp
%          + 25150 * (XFe3_sp + XTi_sp)
%
%     T(K) = TB / TA
%     T(°C) = T(K) - 273
%
% Notes:
% - R is the molar gas constant (J/mol/K): 8.314462618
% - Pressure is supplied in kbar and converted internally to GPa because
%   the pressure terms in Ballhaus et al. (1991) use P in GPa.
% - This implementation follows the algebra exactly as coded for
%   reproducibility; consult Ballhaus et al. (1991) for the original
%   derivation, assumptions, and recommended compositional screening.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Ballhaus1991(rawdata_struct, P_kbar)
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
    error('Ballhaus1991 requires (rawdata_struct, P_kbar).');
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

% Experimental calibration limits reported by Ballhaus et al. (1991).
calibrationT_min_degC = 1040;
calibrationT_max_degC = 1300;
calibrationP_min_GPa = 0.3;
calibrationP_max_GPa = 2.7;

% Pressure is common to all selected mineral pairs in this function call.
% The warning is therefore printed only once, after the first calculation.
P_GPa_input = P_kbar ./ 10;
pressureOutsideCalibration = ...
    P_GPa_input < calibrationP_min_GPa | ...
    P_GPa_input > calibrationP_max_GPa;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or types 'q'
% at the prompt after a calculation.
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

    % All finite mineral-composition inputs used by this thermometer must be
    % strictly positive. NaN values are deliberately excluded from this check
    % so that they remain NaN, propagate through the calculation, and trigger
    % the existing non-stopping warnings after the result is displayed.
    validatePositiveInputs(selectedData_ol, selectedData_sp);

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

    % Warn once when any input pressure lies outside the experimental
    % calibration range of 0.3–2.7 GPa. The calculation is not stopped.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental calibration range ' ...
             'of Ballhaus et al. (1991): 0.3–2.7 GPa (3–27 kbar). ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g GPa.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_GPa_input), ...
            min(P_GPa_input), ...
            max(P_GPa_input));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the
    % experimental calibration range of 1040–1300 degreeC. NaN and Inf are
    % handled separately by the non-finite-result warning below.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the experimental calibration ' ...
             'range of Ballhaus et al. (1991): 1040–1300 degreeC. ' ...
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
    % other than an explicitly NaN input (for example, division by zero).
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
        'Ballhaus1991', ...
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
    'Fe3_cation_apfu', 'Ti_cation_apfu', ...
    'Cr_cation_apfu', 'Al_cation_apfu'};

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

function validatePositiveInputs(data_olivine, data_spinel)
% validatePositiveInputs
% Stop the calculation when a finite required mineral-composition value is
% zero or negative. NaN is intentionally allowed here so that it propagates
% through the thermometer calculation and is reported by the existing
% warning messages.

olivineVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
spinelVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', ...
    'Fe3_cation_apfu', 'Ti_cation_apfu', ...
    'Cr_cation_apfu', 'Al_cation_apfu'};

invalidInputNames = strings(0, 1);

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) <= 0)
        invalidInputNames(end + 1, 1) = ...
            "Olivine." + string(variableName); %#ok<AGROW>
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) <= 0)
        invalidInputNames(end + 1, 1) = ...
            "Spinel." + string(variableName); %#ok<AGROW>
    end
end

if ~isempty(invalidInputNames)
    error(['Ballhaus1991: required mineral-composition values must be > 0. ' ...
           'Zero or negative finite value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_olivine, data_spinel, P_kbar)
% calcTemp
% Compute a single temperature estimate for one olivine row and one spinel
% row, returning a 1-row table containing both intermediate variables and
% final temperatures.
%
% Inputs:
%   data_olivine : 1-row table (olivine apfu cations)
%   data_spinel  : 1-row table (spinel apfu cations)
%   P_kbar       : pressure in kbar; converted internally to GPa
%
% Output:
%   row : 1-row table containing indices, exchange coefficient, TA/TB, and T.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();

% --- Physical constant ---
% Molar gas constant in SI units (J/mol/K).
R = 8.314462618;

% Store both the user-facing pressure unit and the pressure unit used in
% the published equation for transparency and reproducibility.
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.R = repmat(R, nP, 1);

% --- Extract cation data (expected variable names) ---
% Olivine: Mg and Fe are treated as Mg and Fe2+ on the olivine cation basis.
Mg_ol  = repmat(data_olivine.Mg_cation_apfu, nP, 1);
Fe2_ol = repmat(data_olivine.Fe_cation_apfu, nP, 1);

% Spinel: requires both Fe2+ and Fe3+ plus Ti, Cr, Al for correction terms.
Mg_sp  = repmat(data_spinel.Mg_cation_apfu, nP, 1);
Fe2_sp = repmat(data_spinel.Fe_cation_apfu, nP, 1);
Fe3_sp = repmat(data_spinel.Fe3_cation_apfu, nP, 1);
Ti_sp  = repmat(data_spinel.Ti_cation_apfu, nP, 1);
Cr_sp  = repmat(data_spinel.Cr_cation_apfu, nP, 1);
Al_sp  = repmat(data_spinel.Al_cation_apfu, nP, 1);

% --- Compute composition indices used in Ballhaus et al. (1991) ---
% Olivine Mg# expressed as XMg and XFe on the (Mg + Fe2+) basis.
XMg_ol = Mg_ol ./ (Mg_ol + Fe2_ol);
XFe_ol = Fe2_ol ./ (Mg_ol + Fe2_ol);

% Spinel Mg# and ferric ratio terms (definitions follow the implemented algebra).
Mgnum_sp = Mg_sp ./ (Mg_sp + Fe2_sp);
XFe2_sp  = Fe2_sp ./ (Mg_sp + Fe2_sp);
XFe3_sp  = Fe3_sp ./ (Fe2_sp + Fe3_sp);

% In this implementation, Ti is used directly as the normalized cation value.
% (If you want Ti as a true mole fraction on a specific site basis, you must
% pre-normalize accordingly before calling this function.)
XTi_sp   = Ti_sp;

% Chromium fraction among trivalent spinel cations, following the
% definition used in the Ballhaus et al. (1991) thermometer.
XCr_sp   = Cr_sp ./ (Cr_sp + Al_sp + Fe3_sp);

% --- Fe–Mg exchange coefficient between olivine and spinel ---
% K_MgFe(ol-sp) is defined here exactly as in the code header.
KMgFe_ol_sp = (XMg_ol .* XFe2_sp) ./ (XFe_ol .* Mgnum_sp);

% --- Solve temperature (Ballhaus et al., 1991; implemented form) ---
% TA collects the ln(K) term and a constant offset.
TA = R .* log(KMgFe_ol_sp) + 4.705;

% TB collects pressure- and composition-dependent energetic terms.
TB = (6530 + 280*P_GPa + 7000 + 108*P_GPa) .* (1 - 2*XFe_ol) ...
     - 1960 .* (Mgnum_sp - XFe2_sp) ...
     + 16150 .* XCr_sp ...
     + 25150 .* (XFe3_sp + XTi_sp);

% Temperature in Kelvin and Celsius.
T_K   = TB ./ TA;
T_deg = T_K - 273.15;

% --- Pack outputs ---
% Store both intermediate indices and final results, so the user can debug
% inputs and/or propagate uncertainties if needed.
row.XMg_ol = XMg_ol;
row.XFe_ol = XFe_ol;

row.Mgnum_sp = Mgnum_sp;
row.XFe2_sp  = XFe2_sp;
row.XFe3_sp  = XFe3_sp;
row.XTi_sp   = XTi_sp;
row.XCr_sp   = XCr_sp;

row.KMgFe_ol_sp = KMgFe_ol_sp;

row.TA = TA;
row.TB = TB;
row.T_K = T_K;
row.T_deg = T_deg;

end
