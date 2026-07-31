function results = BlundyHolland1990(rawdata_struct, P_kbar)
% functions/+thermo/+Amphibole/BlundyHolland1990.m
% Tested with MATLAB R2024b
%
% Amphibole-Plagioclase thermometer
% Blundy, J.D. and Holland, T.J.B. (1990)
% Contributions to Mineralogy and Petrology, 104, 208–224
% DOI: https://doi.org/10.1007/BF00306444
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Amphibole analysis and one
% Plagioclase analysis (selected by the user from tables) and calculates
% temperature using the Blundy & Holland (1990) amphibole-plagioclase
% thermometer for silica-saturated assemblages.
%
% The function accepts pressure as either a scalar or a vector. Therefore,
% it can be used from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. For every selected Amphibole-Plagioclase pair,
% one output row is produced for each supplied pressure value.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Amphibole-Plagioclase pair and stores
% all result blocks in a single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Blundy & Holland (1990) calibrated the thermometer using experimental
% amphibole-plagioclase data for silica-saturated assemblages. The principal
% application limits recommended by the authors are:
%
%   Temperature       : 500–1100 degreeC
%   Amphibole Si      : Si < 7.8 apfu
%   Plagioclase       : less calcic than An92 (Xan < 0.92)
%   Assemblage        : equilibrium Amphibole + Plagioclase in a
%                       silica-saturated / quartz-saturated rock
%   Typical uncertainty: approximately +/-75 degreeC (2 sigma)
%
% These limits are stated in the abstract on p. 208, discussed in the
% Errors and uncertainties section on pp. 218–219, and summarized in the
% Conclusions on pp. 222–223. The final thermometer equation is given as
% Eq. (4b) on pp. 217–218.
%
% The complete experimental dataset summarized in Table 1 on p. 209 spans
% approximately 1–23 kbar and 400–1150 degreeC. However, the final
% recommended temperature range is 500–1100 degreeC, and the authors state
% that the pressure dependence is poorly constrained. The equation requires
% independently estimated pressure and should not be used as a barometer
% (abstract, p. 208; barometer discussion on p. 218).
%
% Additional cautions:
%   - Below approximately 600 degreeC, ordering, exsolution, and incomplete
%     equilibration in amphibole and/or plagioclase may increase uncertainty
%     (Conclusions, p. 222).
%   - Calculated temperatures may record near-solidus or subsolidus
%     re-equilibration / closure rather than peak crystallization or peak
%     metamorphic temperature (Applications, pp. 219–222).
%   - Amphibole and plagioclase should be texturally and chemically
%     consistent with equilibrium. Zoned, relict, or reaction-related pairs
%     require particular caution (Applications, pp. 219–222).
%   - Silica saturation cannot be determined automatically from the input
%     tables and must be evaluated independently by the user.
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) input pressure is outside the broad experimental dataset range of
%      1–23 kbar,
%   2) a finite calculated temperature is outside 500–1100 degreeC,
%   3) amphibole Si is not < 7.8 apfu,
%   4) plagioclase is not less calcic than An92,
%   5) a calculation input contains NaN, or
%   6) a non-finite temperature is calculated.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole   : table
%   rawdata_struct.Plagioclase : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following normalized cation variables.
%
% Required Amphibole variable used in the thermometer:
%   Si_cation_apfu
%
% Required Plagioclase variables used in the thermometer:
%   Ca_cation_apfu
%   Na_cation_apfu
%
% Optional Plagioclase variable used when present:
%   K_cation_apfu
%
% If K_cation_apfu is absent, K is treated as zero. If the column is present
% but its selected value is NaN, NaN is retained and propagated through the
% calculation; it is not replaced by zero.
%
% Optional Amphibole variables stored in the output when present:
%   Ti_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%   Mn_cation_apfu
%   Cr_cation_apfu
%   Fe3_cation_apfu
%
% Optional Plagioclase variables stored in the output when present:
%   Si_cation_apfu
%   Al_cation_apfu
%
% All finite mineral-composition values used directly in the thermometer
% must be greater than or equal to zero. Finite negative values and infinite
% values stop the calculation with an error. NaN values are retained as
% missing values, propagated through the calculation, and reported by
% non-stopping fprintf warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Blundy & Holland, 1990; Eq. 4b)
%
%   T(K) = (0.677*P - 48.98 + Y) ./ (-0.0429 - R*ln(K))
%
% where
%
%   K = ((Si_amp - 4) / (8 - Si_amp)) * Xab_plag
%
%   Xab_plag = Na / (Ca + Na + K)
%   Xan_plag = Ca / (Ca + Na + K)
%   Xor_plag = K  / (Ca + Na + K)
%
% and
%
%   Y = 0                                 for Xab_plag > 0.5
%     = -8.06 + 25.5*(1 - Xab_plag)^2    for Xab_plag <= 0.5
%
% Units:
%   T in K
%   P in kbar
%   R = 0.0083144 kJ K^-1 mol^-1
%
% -------------------------------------------------------------------------
% Syntax:
%   results = BlundyHolland1990(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Amphibole and Plagioclase tables
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Amphibole-Plagioclase pair. The output variable
%             set is intended to remain stable for downstream processing.
%

%% Input validation
% Basic argument checks prevent silent failures caused by missing inputs or
% invalid pressure values. Pressure may be supplied as a scalar or vector.
if nargin < 2
    error('BlundyHolland1990 requires (rawdata_struct, P_kbar).');
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
% Extract the required tables from the input struct. The tables are not
% modified; only selected rows and relevant columns are read.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Amphibole') || ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end
if ~isfield(rawdata_struct, 'Plagioclase') || ~istable(rawdata_struct.Plagioclase)
    error('rawdata_struct must contain table: rawdata_struct.Plagioclase');
end

dataset_amp = rawdata_struct.Amphibole;
dataset_plag = rawdata_struct.Plagioclase;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated concatenation of the full results table inside the interactive
% loop is avoided because it repeatedly reallocates and copies the table.
%
% The cell buffer is preallocated and doubled only when necessary. All table
% blocks are concatenated once after the interactive loop finishes.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Recommended temperature range and broad experimental dataset pressure
% range reported by Blundy & Holland (1990).
applicationT_min_degC = 500;
applicationT_max_degC = 1100;
experimentalP_min_kbar = 1;
experimentalP_max_kbar = 23;
applicationSi_max_apfu = 7.8;
applicationXan_max = 0.92;

% Pressure is common to every selected mineral pair in this function call.
% The pressure warning is therefore printed only once.
pressureOutsideExperimentalRange = ...
    P_kbar < experimentalP_min_kbar | ...
    P_kbar > experimentalP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or selects
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % ----- Amphibole selection -----
    % Assumption: the first column stores an identifier displayed to the user.
    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, ok] = listdlg( ...
        'PromptString', 'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_amp)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

    % ----- Plagioclase selection -----
    disp('=== Step 4: Selecting a data code from the list (Plagioclase) ===');

    dataCodes_plag = dataset_plag{:, 1};

    [selectedIdx_plag, ok] = listdlg( ...
        'PromptString', 'Please select the Plagioclase data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_plag)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_plag)
        disp('Selection canceled');
        break;
    end

    selectedCode_plag = dataCodes_plag(selectedIdx_plag);
    disp(['Plagioclase selected: ' char(string(selectedCode_plag))]);

    % ----- Calculation -----
    % Amphibole and plagioclase are selected independently; row indices do
    % not need to correspond between the two datasets.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);
    selectedData_plag = dataset_plag(selectedIdx_plag, :);

    % Check only the mineral-composition values used directly in the
    % thermometer. NaN is allowed and propagated; finite negative or infinite
    % values stop the calculation.
    validateNonNegativeInputs(selectedData_amp, selectedData_plag);
    nanInputNames = findNaNInputs(selectedData_amp, selectedData_plag);

    row = calcTemp(selectedData_amp, selectedData_plag, P_kbar, ...
        experimentalP_min_kbar, experimentalP_max_kbar, ...
        applicationT_min_degC, applicationT_max_degC, ...
        applicationSi_max_apfu, applicationXan_max);

    % Store selected identifiers for traceability. A vector pressure input
    % produces one row per pressure, so identifiers are repeated accordingly.
    row.dataCode_amphibole = ...
        repmat(string(selectedCode_amp), height(row), 1);
    row.dataCode_plagioclase = ...
        repmat(string(selectedCode_plag), height(row), 1);
    row = movevars(row, ...
        {'dataCode_amphibole','dataCode_plagioclase'}, 'Before', 1);

    % Store the result as one table block. If the preallocated buffer becomes
    % full, double its capacity. This avoids changing the size of the complete
    % results table on every loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_amp)) ' & ' ...
            char(string(selectedCode_plag)) ': ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_amp)) ' & ' ...
            char(string(selectedCode_plag)) ': ' ...
            num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside the broad experimental
    % dataset range of 1–23 kbar. This range is not presented as a tightly
    % constrained barometric calibration; calculation is not stopped.
    if any(pressureOutsideExperimentalRange) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the broad experimental dataset ' ...
             'range summarized by Blundy & Holland (1990): 1–23 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g kbar. The pressure dependence is poorly ' ...
             'constrained, and this equation should not be used as a barometer.\n'], ...
            sum(pressureOutsideExperimentalRange), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the authors'
    % recommended application range of 500–1100 degreeC.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideApplication = finiteTemperature & ...
        (row.T_deg < applicationT_min_degC | ...
         row.T_deg > applicationT_max_degC);

    if any(temperatureOutsideApplication)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the recommended ' ...
             'application range of Blundy & Holland (1990): 500–1100 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideApplication), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    % Warn when the selected amphibole exceeds the recommended Si limit.
    if isfinite(row.Si_amp(1)) && row.Si_amp(1) >= applicationSi_max_apfu
        fprintf(2, ...
            ['WARNING: Amphibole Si is outside the recommended compositional ' ...
             'range of Blundy & Holland (1990): Si must be < 7.8 apfu. ' ...
             'Selected value = %.4g apfu for %s.\n'], ...
            row.Si_amp(1), ...
            char(string(selectedCode_amp)));
    end

    % Warn when the selected plagioclase is An92 or more calcic.
    if isfinite(row.Xan_plag(1)) && row.Xan_plag(1) >= applicationXan_max
        fprintf(2, ...
            ['WARNING: Plagioclase is outside the recommended compositional ' ...
             'range of Blundy & Holland (1990): Xan must be < 0.92 (An92). ' ...
             'Selected Xan = %.4g for %s.\n'], ...
            row.Xan_plag(1), ...
            char(string(selectedCode_plag)));
    end

    % Print a non-stopping warning immediately after the temperature result
    % when a thermometer input contains NaN. fprintf is used so that the
    % message remains visible even if MATLAB warnings are disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; the calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain a result-based check for NaN/Inf values caused by a NaN input or
    % by a mathematically invalid composition, such as zero component sum.
    invalidTemperature = ~isfinite(row.T_deg);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat the selection and calculation procedure.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'BlundyHolland1990', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after all selections have been
% completed. Return an empty table when no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_amphibole, data_plagioclase)
% findNaNInputs
% Return the names of thermometer input variables containing NaN. This
% function does not throw an error for NaN values; it prepares a warning
% message for the calling function.

amphiboleVariables = {'Si_cation_apfu'};
plagioclaseVariables = {'Ca_cation_apfu', 'Na_cation_apfu', ...
    'K_cation_apfu'};

maxNames = numel(amphiboleVariables) + numel(plagioclaseVariables);
nanInputNamesBuffer = strings(maxNames, 1);
nNaNInputs = 0;

for i = 1:numel(amphiboleVariables)
    variableName = amphiboleVariables{i};
    if ismember(variableName, data_amphibole.Properties.VariableNames)
        variableValue = data_amphibole.(variableName);
        if any(isnan(variableValue(:)))
            nNaNInputs = nNaNInputs + 1;
            nanInputNamesBuffer(nNaNInputs) = ...
                "Amphibole." + string(variableName);
        end
    end
end

for i = 1:numel(plagioclaseVariables)
    variableName = plagioclaseVariables{i};
    if ismember(variableName, data_plagioclase.Properties.VariableNames)
        variableValue = data_plagioclase.(variableName);
        if any(isnan(variableValue(:)))
            nNaNInputs = nNaNInputs + 1;
            nanInputNamesBuffer(nNaNInputs) = ...
                "Plagioclase." + string(variableName);
        end
    end
end

nanInputNames = nanInputNamesBuffer(1:nNaNInputs);

end

function validateNonNegativeInputs(data_amphibole, data_plagioclase)
% validateNonNegativeInputs
% Stop the calculation when a thermometer input contains a finite negative
% value or an infinite value. Zero is allowed. NaN is intentionally allowed
% so that it remains missing, propagates through the calculation, and is
% reported by non-stopping fprintf warnings.

amphiboleVariables = {'Si_cation_apfu'};
plagioclaseVariables = {'Ca_cation_apfu', 'Na_cation_apfu', ...
    'K_cation_apfu'};

maxNames = numel(amphiboleVariables) + numel(plagioclaseVariables);
invalidInputNamesBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

for i = 1:numel(amphiboleVariables)
    variableName = amphiboleVariables{i};
    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        error('Amphibole table must contain variable: %s', variableName);
    end

    variableValue = data_amphibole.(variableName);
    validateNumericScalar(variableValue, variableName);

    if any(isinf(variableValue(:)) | ...
            (isfinite(variableValue(:)) & variableValue(:) < 0))
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputNamesBuffer(nInvalidInputs) = ...
            "Amphibole." + string(variableName);
    end
end

for i = 1:numel(plagioclaseVariables)
    variableName = plagioclaseVariables{i};

    % K_cation_apfu is optional. If absent, it is treated as zero by
    % preparePlagioclaseRow and therefore does not require validation here.
    if ~ismember(variableName, data_plagioclase.Properties.VariableNames)
        if strcmp(variableName, 'K_cation_apfu')
            continue;
        end
        error('Plagioclase table must contain variable: %s', variableName);
    end

    variableValue = data_plagioclase.(variableName);
    validateNumericScalar(variableValue, variableName);

    if any(isinf(variableValue(:)) | ...
            (isfinite(variableValue(:)) & variableValue(:) < 0))
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputNamesBuffer(nInvalidInputs) = ...
            "Plagioclase." + string(variableName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputNamesBuffer(1:nInvalidInputs);
    error(['BlundyHolland1990: thermometer input values must be ' ...
           'non-negative or NaN. Finite negative or infinite value(s) ' ...
           'were found in: ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_amphibole, data_plagioclase, P_kbar, ...
        experimentalP_min_kbar, experimentalP_max_kbar, ...
        applicationT_min_degC, applicationT_max_degC, ...
        applicationSi_max_apfu, applicationXan_max)
% calcTemp
% Compute Blundy & Holland (1990) temperatures for one selected Amphibole-
% Plagioclase pair at one or more pressure values.
%
% Inputs:
%   data_amphibole  : 1-row table containing amphibole cations
%   data_plagioclase: 1-row table containing plagioclase cations
%   P_kbar          : pressure in kbar, supplied as a scalar or vector
%
% Output:
%   row : table containing one row per pressure value, intermediate
%         variables, applicability flags, and final temperatures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

% --- Physical constant ---
R_kJ = 0.0083144;

% --- Prepare selected mineral rows ---
amp = prepareAmphiboleRow(data_amphibole);
plag = preparePlagioclaseRow(data_plagioclase);

% --- Calculate derived mineral variables ---
site_amp = calcAmphiboleVariables(amp);
site_plag = calcPlagioclaseVariables(plag);

Si_amp_scalar = site_amp.Si;
Al_IV_amp_scalar = site_amp.Al_IV;
Al_VI_amp_scalar = site_amp.Al_VI;

Xab_plag_scalar = site_plag.Xab;
Xan_plag_scalar = site_plag.Xan;
Xor_plag_scalar = site_plag.Xor;

% --- Plagioclase non-ideality term ---
% A NaN Xab value remains NaN through either expression.
if Xab_plag_scalar > 0.5
    Y_kJ_scalar = 0.0;
else
    Y_kJ_scalar = -8.06 + 25.5 .* (1 - Xab_plag_scalar).^2;
end

% --- Equilibrium constant term ---
% The calculation is performed only when K is finite and positive. Invalid
% or missing compositions retain NaN rather than being replaced by zero.
K_scalar = ((Si_amp_scalar - 4) ./ (8 - Si_amp_scalar)) .* ...
    Xab_plag_scalar;

if isfinite(K_scalar) && K_scalar > 0
    lnK_scalar = log(K_scalar);
else
    lnK_scalar = NaN;
end

% --- Temperature calculation for scalar or vector pressure ---
num = 0.677 .* P_kbar - 48.98 + Y_kJ_scalar;
denom = repmat(-0.0429 - R_kJ .* lnK_scalar, nP, 1);

T_K = NaN(nP, 1);
validDenominator = isfinite(denom) & abs(denom) > 1e-12;
T_K(validDenominator) = num(validDenominator) ./ denom(validDenominator);
T_deg = T_K - 273.15;

% --- Replicate scalar mineral variables to match vector pressure output ---
Si_amp = repmat(amp.Si, nP, 1);
Ti_amp = repmat(amp.Ti, nP, 1);
Al_amp = repmat(amp.Al, nP, 1);
Fe_amp = repmat(amp.FeT, nP, 1);
Fe3_amp = repmat(amp.Fe3, nP, 1);
Fe2_amp = repmat(amp.Fe2, nP, 1);
Mg_amp = repmat(amp.Mg, nP, 1);
Ca_amp = repmat(amp.Ca, nP, 1);
Na_amp = repmat(amp.Na, nP, 1);
K_amp = repmat(amp.K, nP, 1);
Mn_amp = repmat(amp.Mn, nP, 1);
Cr_amp = repmat(amp.Cr, nP, 1);

Al_IV_amp = repmat(Al_IV_amp_scalar, nP, 1);
Al_VI_amp = repmat(Al_VI_amp_scalar, nP, 1);

Si_plag = repmat(plag.Si, nP, 1);
Al_plag = repmat(plag.Al, nP, 1);
Ca_plag = repmat(plag.Ca, nP, 1);
Na_plag = repmat(plag.Na, nP, 1);
K_plag = repmat(plag.K, nP, 1);

Xab_plag = repmat(Xab_plag_scalar, nP, 1);
Xan_plag = repmat(Xan_plag_scalar, nP, 1);
Xor_plag = repmat(Xor_plag_scalar, nP, 1);
Y_kJ = repmat(Y_kJ_scalar, nP, 1);
K = repmat(K_scalar, nP, 1);
lnK = repmat(lnK_scalar, nP, 1);
R_kJ_output = repmat(R_kJ, nP, 1);

% --- Applicability flags ---
pressureWithinExperimentalRange = ...
    P_kbar >= experimentalP_min_kbar & ...
    P_kbar <= experimentalP_max_kbar;

temperatureWithinApplicationRange = ...
    isfinite(T_deg) & ...
    T_deg >= applicationT_min_degC & ...
    T_deg <= applicationT_max_degC;

amphiboleSiWithinApplicationRange = ...
    isfinite(Si_amp) & Si_amp < applicationSi_max_apfu;

plagioclaseWithinApplicationRange = ...
    isfinite(Xan_plag) & Xan_plag < applicationXan_max;

isApplicable = ...
    pressureWithinExperimentalRange & ...
    temperatureWithinApplicationRange & ...
    amphiboleSiWithinApplicationRange & ...
    plagioclaseWithinApplicationRange;

% --- Pack outputs ---
row = table();
row.P_kbar = P_kbar;
row.R_kJ = R_kJ_output;

row.Si_amp = Si_amp;
row.Ti_amp = Ti_amp;
row.Al_amp = Al_amp;
row.Fe_amp = Fe_amp;
row.Fe3_amp = Fe3_amp;
row.Fe2_amp = Fe2_amp;
row.Mg_amp = Mg_amp;
row.Ca_amp = Ca_amp;
row.Na_amp = Na_amp;
row.K_amp = K_amp;
row.Mn_amp = Mn_amp;
row.Cr_amp = Cr_amp;

row.Al_IV_amp = Al_IV_amp;
row.Al_VI_amp = Al_VI_amp;

row.Si_plag = Si_plag;
row.Al_plag = Al_plag;
row.Ca_plag = Ca_plag;
row.Na_plag = Na_plag;
row.K_plag = K_plag;

row.Xab_plag = Xab_plag;
row.Xan_plag = Xan_plag;
row.Xor_plag = Xor_plag;

row.Y_kJ = Y_kJ;
row.K = K;
row.lnK = lnK;
row.num = num;
row.denom = denom;

row.T_K = T_K;
row.T_deg = T_deg;

row.pressureWithinExperimentalRange = pressureWithinExperimentalRange;
row.temperatureWithinApplicationRange = temperatureWithinApplicationRange;
row.amphiboleSiWithinApplicationRange = ...
    amphiboleSiWithinApplicationRange;
row.plagioclaseWithinApplicationRange = ...
    plagioclaseWithinApplicationRange;
row.isApplicable = isApplicable;

end

function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Extract one selected amphibole analysis from a 1-row table. Missing
% optional columns default to zero to preserve the existing output schema.
% A NaN value in an existing column remains NaN and is never converted to
% zero.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();

amp.Si = getVarOrError(data_amphibole, ...
    'Si_cation_apfu', 'Amphibole');

amp.Ti = getVarOrDefaultByNameList(data_amphibole, ...
    {'Ti_cation_apfu'}, 0);
amp.Al = getVarOrDefaultByNameList(data_amphibole, ...
    {'Al_cation_apfu'}, 0);
amp.FeT = getVarOrDefaultByNameList(data_amphibole, ...
    {'Fe_cation_apfu'}, 0);
amp.Mn = getVarOrDefaultByNameList(data_amphibole, ...
    {'Mn_cation_apfu'}, 0);
amp.Mg = getVarOrDefaultByNameList(data_amphibole, ...
    {'Mg_cation_apfu'}, 0);
amp.Ca = getVarOrDefaultByNameList(data_amphibole, ...
    {'Ca_cation_apfu'}, 0);
amp.Na = getVarOrDefaultByNameList(data_amphibole, ...
    {'Na_cation_apfu'}, 0);
amp.K = getVarOrDefaultByNameList(data_amphibole, ...
    {'K_cation_apfu'}, 0);
amp.Cr = getVarOrDefaultByNameList(data_amphibole, ...
    {'Cr_cation_apfu'}, 0);
amp.Fe3 = getVarOrDefaultByNameList(data_amphibole, ...
    {'Fe3_cation_apfu'}, 0);

% Do not use max(0, ...) here because it could hide invalid values or alter
% missing-value behavior. NaN is retained naturally by subtraction.
amp.Fe2 = amp.FeT - amp.Fe3;

if isfinite(amp.Fe3) && isfinite(amp.FeT) && amp.Fe3 > amp.FeT + 1e-12
    error('Amphibole contains Fe3_cation_apfu > Fe_cation_apfu.');
end

end

function plag = preparePlagioclaseRow(data_plagioclase)
% preparePlagioclaseRow
% Extract one selected plagioclase analysis from a 1-row table. K is treated
% as zero only when K_cation_apfu is absent. A NaN value in an existing
% K_cation_apfu column remains NaN and propagates through the calculation.

if height(data_plagioclase) ~= 1
    error('Plagioclase input must be a 1-row table.');
end

plag = struct();

plag.Ca = getVarOrError(data_plagioclase, ...
    'Ca_cation_apfu', 'Plagioclase');
plag.Na = getVarOrError(data_plagioclase, ...
    'Na_cation_apfu', 'Plagioclase');

plag.Si = getVarOrDefaultByNameList(data_plagioclase, ...
    {'Si_cation_apfu'}, 0);
plag.Al = getVarOrDefaultByNameList(data_plagioclase, ...
    {'Al_cation_apfu'}, 0);
plag.K = getVarOrDefaultByNameList(data_plagioclase, ...
    {'K_cation_apfu'}, 0);

end

function site = calcAmphiboleVariables(amp)
% calcAmphiboleVariables
% Calculate the amphibole variables stored or used by the thermometer.

site = struct();
site.Si = amp.Si;

if isnan(amp.Si)
    site.Al_IV = NaN;
else
    site.Al_IV = max(0, 8 - amp.Si);
end

if isnan(amp.Al) || isnan(site.Al_IV)
    site.Al_VI = NaN;
elseif amp.Al > 0
    site.Al_VI = max(0, amp.Al - site.Al_IV);
else
    site.Al_VI = NaN;
end

end

function site = calcPlagioclaseVariables(plag)
% calcPlagioclaseVariables
% Calculate albite, anorthite, and orthoclase fractions. A zero component
% sum or a NaN component produces NaN fractions; calculation is not stopped.

site = struct();
componentSum = plag.Ca + plag.Na + plag.K;

if isfinite(componentSum) && componentSum > 0
    site.Xab = plag.Na ./ componentSum;
    site.Xan = plag.Ca ./ componentSum;
    site.Xor = plag.K ./ componentSum;
else
    site.Xab = NaN;
    site.Xan = NaN;
    site.Xor = NaN;
end

end

function value = getVarOrError(tbl, varName, mineralLabel)
% getVarOrError
% Retrieve a required numeric scalar variable from a 1-row table. NaN is
% allowed and retained; value-range checks are handled separately.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);
validateNumericScalar(value, varName);

end

function value = getVarOrDefaultByNameList(tbl, varNameList, defaultValue)
% getVarOrDefaultByNameList
% Retrieve the first available numeric scalar variable from a list. Return
% defaultValue only when none of the listed columns exists. NaN in an
% existing column is retained and is not replaced by defaultValue.

value = defaultValue;

for i = 1:numel(varNameList)
    varName = varNameList{i};
    if ismember(varName, tbl.Properties.VariableNames)
        value = tbl.(varName);
        validateNumericScalar(value, varName);
        return;
    end
end

end

function validateNumericScalar(value, varName)
% validateNumericScalar
% Require a numeric scalar in every selected 1-row mineral table column.

if ~isnumeric(value) || ~isscalar(value) || ~isreal(value)
    error('Variable %s must be a real numeric scalar in a 1-row table.', ...
        varName);
end

end
