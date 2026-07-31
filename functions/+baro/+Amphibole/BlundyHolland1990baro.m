function results = BlundyHolland1990baro(rawdata_struct, T_degreeC)
% functions/+baro/+Amphibole/BlundyHolland1990baro.m
% Tested with MATLAB R2024b
%
% Pressure-form calculation based on Amphibole-Plagioclase equilibrium
% Blundy, J.D. and Holland, T.J.B. (1990)
% Contributions to Mineralogy and Petrology, 104, 208-224
% DOI: https://doi.org/10.1007/BF00306444
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Amphibole analysis and one
% Plagioclase analysis (selected by the user from tables) and calculates
% pressure using the pressure-form regression of the edenite-tremolite-
% albite-quartz equilibrium presented by Blundy and Holland (1990).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Amphibole-Plagioclase pair, one
% output row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% IMPORTANT: Blundy and Holland (1990) developed and recommended this
% equilibrium primarily as an AMPHIBOLE-PLAGIOCLASE THERMOMETER, not as a
% reliable barometer. The abstract explicitly states that the pressure
% dependence is poorly constrained and that the equilibria are not suitable
% for barometry (p. 208). The pressure-form regression is presented as
% Eq. (5) on p. 218 and has a residual standard deviation of 1.95 kbar.
% The authors state that pressures cannot be expected to be estimated to
% better than approximately +/-4 kbar (p. 218).
%
% The experimental amphibole dataset summarized in Table 1 spans broadly:
%
%   Temperature : approximately 400-1150 degreeC
%   Pressure    : approximately 1-23 kbar
%   Materials   : basaltic to felsic and peridotitic starting compositions
%
% These dataset limits are summarized in the Introduction and Table 1 on
% pp. 208-209. They are NOT a validated pressure-calibration range for the
% barometer. In this implementation, 1-23 kbar is used only as an
% experimental-data envelope for a non-stopping warning.
%
% For the recommended thermometer formulation, Blundy and Holland (1990)
% give the following practical limits in the Errors and uncertainties
% section on p. 218 and in the Conclusions on pp. 222-223:
%
%   Temperature            : 500-1100 degreeC
%   Amphibole Si            : Si < 7.8 apfu
%   Plagioclase composition : An < 92 mol% (XAn < 0.92)
%   Silica activity         : quartz/silica-saturated assemblages preferred
%
% Errors increase strongly as amphibole Si approaches 8 apfu and as
% plagioclase approaches very calcic compositions. At temperatures below
% approximately 600 degreeC, ordering, unmixing, or incomplete equilibration
% in amphibole and plagioclase may produce additional uncertainty
% (pp. 218-222).
%
% Amphibole and plagioclase must represent an equilibrium pair. Patchy or
% concentric zoning, inherited crystals, multiple amphibole generations,
% subsolidus re-equilibration, and non-adjacent mineral analyses may produce
% misleading pressures (applications discussed on pp. 219-222).
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 500-1100 degreeC,
%   2) finite calculated pressure is outside the 1-23 kbar experimental-data
%      envelope,
%   3) amphibole Si is outside the equation/recommended domain,
%   4) plagioclase is An92 or more calcic,
%   5) a required calculation input contains NaN, or
%   6) a calculated pressure is NaN or Inf.
%
% Even when all checks are passed, calculated pressure must be treated only
% as an approximate reference value because the pressure dependence is
% poorly constrained and the practical uncertainty is approximately
% +/-4 kbar.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole   : table
%   rawdata_struct.Plagioclase : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following normalized-cation variables:
%
%   Amphibole table variables:
%     Si_cation_apfu           % required by the pressure equation
%     Ti_cation_apfu
%     Al_cation_apfu
%     Fe_cation_apfu           % total Fe in the input cation table
%     Mg_cation_apfu
%     Ca_cation_apfu
%     Na_cation_apfu
%
%   Optional Amphibole variables retained in the output when present:
%     Fe3_cation_apfu
%     Mn_cation_apfu
%     K_cation_apfu
%     Cr_cation_apfu
%
%   Plagioclase table variables:
%     Si_cation_apfu
%     Al_cation_apfu
%     Ca_cation_apfu           % required for XAn and XAb
%     Na_cation_apfu           % required for XAb
%     K_cation_apfu            % required for feldspar normalization
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if it makes a logarithm or ratio
% undefined, the resulting NaN/Inf is retained and reported.
%
% No liquid composition is used by this barometer. Therefore, treatment of
% Liq F and Cl and cationTotal_liq is not applicable to this function.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
% Pressure form of the equilibrium:
%
%   Edenite + 4 Quartz = Tremolite + Albite
%
%   P(kbar) =
%     [52.16 - 0.0617*T(K) - R*T(K)*ln(K1) - Z] / 1.218
%
% where:
%
%   K1 = (27/256) * [(Si_amp - 4) / (8 - Si_amp)] * XAb_plag
%
%   XAb_plag = Na / (Ca + Na + K)
%   XAn_plag = Ca / (Ca + Na + K)
%   XOr_plag = K  / (Ca + Na + K)
%
%   Z = 0                                           for XAb > 0.5
%   Z = 24.9*(1 - XAb)^2 - 0.25*24.9               for XAb <= 0.5
%
%   R = 0.008314 kJ K^-1 mol^-1
%
% K1, including the normalization factor 27/256, is given for Model 2 in
% Table 4 (p. 216) and Eq. (4a) (p. 217). The pressure equation and its
% regression parameters are given in Eq. (5) and Table 6b (p. 218).
% The simplified K of Eq. (4b), which omits 27/256, must not be combined
% directly with the Eq. (5)/Table 6b regression constants.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = BlundyHolland1990baro(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing Amphibole and Plagioclase tables
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Amphibole-Plagioclase pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('BlundyHolland1990baro requires (rawdata_struct, T_degreeC).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(T_degreeC) || isempty(T_degreeC) || ~isvector(T_degreeC) || ...
        any(isinf(T_degreeC(:))) || ...
        any(isfinite(T_degreeC(:)) & T_degreeC(:) < 0)
    error(['T_degreeC must be a non-negative numeric scalar or vector. ' ...
           'NaN is allowed and retained, but Inf and finite negative values are prohibited.']);
end

T_degreeC = T_degreeC(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The source tables are
% not modified; selected rows are read during each calculation.
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
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Practical temperature limits recommended for the associated thermometer.
recommendedT_min_degreeC = 500;
recommendedT_max_degreeC = 1100;

% Overall pressure range represented by the experimental dataset. This is
% used only as a warning envelope and is not a validated barometer range.
experimentalP_min_kbar = 1;
experimentalP_max_kbar = 23;

temperatureOutsideRecommended = isfinite(T_degreeC) & ...
    (T_degreeC < recommendedT_min_degreeC | ...
     T_degreeC > recommendedT_max_degreeC);
temperatureWarningIssued = false;
barometerCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % ----- Amphibole selection -----
    % The first table column is used only as the displayed data identifier.
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
    % Amphibole and plagioclase are selected independently; row numbers are
    % not assumed to correspond between the two tables.
    disp('=== Step 5: Calculating the pressure ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);
    selectedData_plag = dataset_plag(selectedIdx_plag, :);

    % Check NaN only in variables directly used by the pressure equation.
    % NaN values are not changed and do not stop the calculation.
    nanInputNames = findNaNInputs(selectedData_amp, selectedData_plag, T_degreeC);

    % Reject only Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(selectedData_amp, selectedData_plag);

    row = calcPressure(selectedData_amp, selectedData_plag, T_degreeC);

    % Repeat identifiers for all temperatures in the current calculation.
    row.dataCode_amphibole = repmat(string(selectedCode_amp), height(row), 1);
    row.dataCode_plagioclase = repmat(string(selectedCode_plag), height(row), 1);
    row = movevars(row, ...
        {'dataCode_amphibole','dataCode_plagioclase'}, 'Before', 1);

    % Store one block per selected mineral pair. Expand the cell buffer only
    % when its current capacity has been exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_amp)) ' & ' ...
            char(string(selectedCode_plag)) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_amp)) ' & ' ...
            char(string(selectedCode_plag)) ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Print the barometer-specific limitation once per function call.
    if ~barometerCautionIssued
        fprintf(2, ...
            ['CAUTION: Blundy and Holland (1990) state that the pressure dependence ' ...
             'of this equilibrium is poorly constrained and is not suitable for ' ...
             'reliable barometry (p. 208). The pressure regression has a residual ' ...
             'standard deviation of 1.95 kbar, corresponding to an approximate ' ...
             'practical uncertainty of +/-4 kbar (p. 218).\n']);
        barometerCautionIssued = true;
    end

    % Input temperature is common to all selected pairs, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideRecommended) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the recommended range associated ' ...
             'with Blundy and Holland (1990): 500-1100 degreeC (p. 218; pp. 222-223). ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideRecommended), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures fall outside the broad pressure
    % envelope represented by the experiments. This is not presented as a
    % validated calibration range.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideExperimentalEnvelope = finitePressure & ...
        (row.P_kbar < experimentalP_min_kbar | ...
         row.P_kbar > experimentalP_max_kbar);

    if any(pressureOutsideExperimentalEnvelope)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the approximately 1-23 kbar ' ...
             'experimental-data envelope summarized by Blundy and Holland (1990; ' ...
             'pp. 208-209). This envelope is not a validated barometer calibration ' ...
             'range. %d of %d finite pressure point(s) are outside the envelope; ' ...
             'calculated finite range = %.4g-%.4g kbar for %s & %s.\n'], ...
            sum(pressureOutsideExperimentalEnvelope), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    % Warn when the selected mineral compositions lie outside the practical
    % limits discussed for the associated thermometer formulation.
    finiteSi = isfinite(row.Si_amp(1));
    if finiteSi && row.Si_amp(1) >= 7.8
        fprintf(2, ...
            ['WARNING: Amphibole Si = %.4g apfu is outside the recommended domain ' ...
             'Si < 7.8 apfu of Blundy and Holland (1990; p. 218) for %s & %s.\n'], ...
            row.Si_amp(1), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    if finiteSi && (row.Si_amp(1) <= 4 || row.Si_amp(1) >= 8)
        fprintf(2, ...
            ['WARNING: Amphibole Si = %.4g apfu is outside the mathematical domain ' ...
             '4 < Si < 8 required for K1 = (27/256)*((Si-4)/(8-Si))*XAb for %s & %s. ' ...
             'The pressure result may be NaN or Inf.\n'], ...
            row.Si_amp(1), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    finiteXan = isfinite(row.Xan_plag(1));
    if finiteXan && row.Xan_plag(1) >= 0.92
        fprintf(2, ...
            ['WARNING: Plagioclase XAn = %.4g is outside the recommended domain ' ...
             'XAn < 0.92 (An < 92 mol%%) of Blundy and Holland (1990; p. 218) ' ...
             'for %s & %s.\n'], ...
            row.Xan_plag(1), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    % List the exact required input names containing NaN. For vector
    % temperature input, the NaN element indices are included.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s & %s: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'pressure may remain NaN.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % A negative finite pressure is mathematically retained but is outside the
    % physical/useful domain and has already triggered the range warning.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & %s ' ...
             '(%d of %d points). The values were retained for diagnostic purposes.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'BlundyHolland1990baro', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered result blocks once. If the user canceled before
% any calculation, return an empty table.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_amphibole, data_plagioclase, T_degreeC)
% findNaNInputs
% Return the names of pressure-equation inputs containing NaN. NaN values do
% not cause an error and are not replaced by zero.

maxNames = 5;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = "T_degreeC(indices=" + indexText + ")";
end

if any(isnan(data_amphibole.Si_cation_apfu(:)))
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Amphibole.Si_cation_apfu";
end

plagioclaseVariables = {'Ca_cation_apfu', 'Na_cation_apfu', 'K_cation_apfu'};
for i = 1:numel(plagioclaseVariables)
    variableName = plagioclaseVariables{i};
    variableValue = data_plagioclase.(variableName);
    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Plagioclase." + string(variableName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_amphibole, data_plagioclase)
% validateNonNegativeInputs
% Reject Inf and finite negative values in variables directly used by the
% pressure equation. Zero and NaN are intentionally allowed and retained.

maxNames = 4;
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

amphiboleValue = data_amphibole.Si_cation_apfu;
if any(isinf(amphiboleValue(:)) | ...
        (isfinite(amphiboleValue(:)) & amphiboleValue(:) < 0))
    nInvalidInputs = nInvalidInputs + 1;
    invalidInputBuffer(nInvalidInputs) = "Amphibole.Si_cation_apfu";
end

plagioclaseVariables = {'Ca_cation_apfu', 'Na_cation_apfu', 'K_cation_apfu'};
for i = 1:numel(plagioclaseVariables)
    variableName = plagioclaseVariables{i};
    variableValue = data_plagioclase.(variableName);
    if any(isinf(variableValue(:)) | ...
            (isfinite(variableValue(:)) & variableValue(:) < 0))
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Plagioclase." + string(variableName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['BlundyHolland1990baro: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_amphibole, data_plagioclase, T_degreeC)
% calcPressure
% Compute pressure for one amphibole row and one plagioclase row at one or
% more input temperatures. NaN values propagate through the calculation.
%
% Inputs:
%   data_amphibole   : 1-row Amphibole table
%   data_plagioclase : 1-row Plagioclase table
%   T_degreeC        : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Physical constant and regression constants.
R_kJ = 0.008314;
WAb_kJ = 24.9;

% Extract one-row cation data. Required columns are checked explicitly;
% optional columns are returned as NaN when absent and are never set to zero.
amp = prepareAmphiboleRow(data_amphibole);
plag = preparePlagioclaseRow(data_plagioclase);

% Plagioclase component fractions. NaN values propagate naturally through
% the sum and ratios. A zero total produces NaN through 0/0.
alkaliEarthSum = plag.Ca + plag.Na + plag.K;
Xab_plag_scalar = plag.Na ./ alkaliEarthSum;
Xan_plag_scalar = plag.Ca ./ alkaliEarthSum;
Xor_plag_scalar = plag.K ./ alkaliEarthSum;

% Equilibrium constant K1 for Eq. (5) / Table 6b. The normalization factor
% 27/256 is required by Table 4 and Eq. (4a). Set ln(K1) to NaN when K1 is
% not finite and strictly positive to avoid complex logarithms.
K_BH90_scalar = (27 / 256) .* ...
    ((amp.Si - 4) ./ (8 - amp.Si)) .* Xab_plag_scalar;
lnK_BH90_scalar = NaN;
if isfinite(K_BH90_scalar) && K_BH90_scalar > 0
    lnK_BH90_scalar = log(K_BH90_scalar);
end

% Plagioclase non-ideality correction. NaN XAb produces NaN Z.
if isnan(Xab_plag_scalar)
    Z_kJ_scalar = NaN;
elseif Xab_plag_scalar > 0.5
    Z_kJ_scalar = 0.0;
else
    Z_kJ_scalar = WAb_kJ .* (1 - Xab_plag_scalar).^2 - ...
        0.25 .* WAb_kJ;
end

% Expand composition-dependent scalars to the temperature-vector length.
Si_amp = repmat(amp.Si, nT, 1);
Ti_amp = repmat(amp.Ti, nT, 1);
Al_amp = repmat(amp.Al, nT, 1);
FeT_amp = repmat(amp.FeT, nT, 1);
Fe2_amp = repmat(amp.Fe2, nT, 1);
Fe3_amp = repmat(amp.Fe3, nT, 1);
Mg_amp = repmat(amp.Mg, nT, 1);
Ca_amp = repmat(amp.Ca, nT, 1);
Na_amp = repmat(amp.Na, nT, 1);
K_amp = repmat(amp.K, nT, 1);
Mn_amp = repmat(amp.Mn, nT, 1);
Cr_amp = repmat(amp.Cr, nT, 1);

Si_plag = repmat(plag.Si, nT, 1);
Al_plag = repmat(plag.Al, nT, 1);
Ca_plag = repmat(plag.Ca, nT, 1);
Na_plag = repmat(plag.Na, nT, 1);
K_plag = repmat(plag.K, nT, 1);

Xab_plag = repmat(Xab_plag_scalar, nT, 1);
Xan_plag = repmat(Xan_plag_scalar, nT, 1);
Xor_plag = repmat(Xor_plag_scalar, nT, 1);
K_BH90 = repmat(K_BH90_scalar, nT, 1);
lnK_BH90 = repmat(lnK_BH90_scalar, nT, 1);
WAb = repmat(WAb_kJ, nT, 1);
Z_kJ = repmat(Z_kJ_scalar, nT, 1);
R_output = repmat(R_kJ, nT, 1);

% Practical site-allocation outputs retained for compatibility and
% diagnostics. NaN inputs remain NaN in these outputs.
site_amp = calcAmphiboleSites(amp);
Al_4_amp = repmat(site_amp.Al_4, nT, 1);
Al_6_amp = repmat(site_amp.Al_6, nT, 1);
Ca_B_amp = repmat(site_amp.Ca_B, nT, 1);
Na_B_amp = repmat(site_amp.Na_B, nT, 1);
Na_A_amp = repmat(site_amp.Na_A, nT, 1);
K_A_amp = repmat(site_amp.K_A, nT, 1);
vA_amp = repmat(site_amp.vA, nT, 1);

% Pressure calculation. No finite-value guard is used here: NaN inputs
% remain NaN and propagate to P_kbar as requested.
P_kbar = (52.16 - 0.0617 .* T_K - ...
    R_kJ .* T_K .* lnK_BH90 - Z_kJ) ./ 1.218;

% Applicability flags. These are diagnostic flags only and do not imply that
% this pressure form is a reliable barometer.
isWithinEquationDomain = ...
    isfinite(Si_amp) & Si_amp > 4 & Si_amp < 8 & ...
    isfinite(Xab_plag) & Xab_plag > 0;

isWithinCompositionDomain = ...
    isWithinEquationDomain & ...
    Si_amp < 7.8 & ...
    isfinite(Xan_plag) & Xan_plag < 0.92;

isWithinRecommendedTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 500 & T_degreeC <= 1100;

isWithinExperimentalPEnvelope = ...
    isfinite(P_kbar) & P_kbar >= 1 & P_kbar <= 23;

isBarometerReliable = false(nT, 1);

% Pack outputs using pre-sized vectors of equal height.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.R_kJ = R_output;

row.Si_amp = Si_amp;
row.Ti_amp = Ti_amp;
row.Al_amp = Al_amp;
row.FeT_amp = FeT_amp;
row.Fe2_amp = Fe2_amp;
row.Fe3_amp = Fe3_amp;
row.Mg_amp = Mg_amp;
row.Ca_amp = Ca_amp;
row.Na_amp = Na_amp;
row.K_amp = K_amp;
row.Mn_amp = Mn_amp;
row.Cr_amp = Cr_amp;

row.Al_4_amp = Al_4_amp;
row.Al_6_amp = Al_6_amp;
row.Ca_B_amp = Ca_B_amp;
row.Na_B_amp = Na_B_amp;
row.Na_A_amp = Na_A_amp;
row.K_A_amp = K_A_amp;
row.vA_amp = vA_amp;

row.Si_plag = Si_plag;
row.Al_plag = Al_plag;
row.Ca_plag = Ca_plag;
row.Na_plag = Na_plag;
row.K_plag = K_plag;

row.Xab_plag = Xab_plag;
row.Xan_plag = Xan_plag;
row.Xor_plag = Xor_plag;

row.K_BH90 = K_BH90;
row.lnK_BH90 = lnK_BH90;
row.WAb_kJ = WAb;
row.Z_kJ = Z_kJ;

row.P_kbar = P_kbar;
row.P_uncertainty_1sigma_kbar = repmat(1.95, nT, 1);
row.P_uncertainty_approx_2sigma_kbar = repmat(4.0, nT, 1);

row.isWithinEquationDomain = isWithinEquationDomain;
row.isWithinCompositionDomain = isWithinCompositionDomain;
row.isWithinRecommendedTRange = isWithinRecommendedTRange;
row.isWithinExperimentalPEnvelope = isWithinExperimentalPEnvelope;

% Backward-compatible aliases retained from the original implementation.
row.isApplicable_composition = isWithinCompositionDomain;
row.isRecommended_T_range = isWithinRecommendedTRange;
row.isBarometerReliable = isBarometerReliable;

end

function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Extract one-row amphibole cation data. Required variables must exist;
% optional variables are retained as NaN when absent.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();
amp.Si = getVarOrError(data_amphibole, 'Si_cation_apfu', 'Amphibole');
amp.Ti = getVarOrError(data_amphibole, 'Ti_cation_apfu', 'Amphibole');
amp.Al = getVarOrError(data_amphibole, 'Al_cation_apfu', 'Amphibole');
amp.FeT = getVarOrError(data_amphibole, 'Fe_cation_apfu', 'Amphibole');
amp.Mg = getVarOrError(data_amphibole, 'Mg_cation_apfu', 'Amphibole');
amp.Ca = getVarOrError(data_amphibole, 'Ca_cation_apfu', 'Amphibole');
amp.Na = getVarOrError(data_amphibole, 'Na_cation_apfu', 'Amphibole');

amp.Fe3 = getVarOrNaN(data_amphibole, 'Fe3_cation_apfu');
amp.Mn = getVarOrNaN(data_amphibole, 'Mn_cation_apfu');
amp.K = getVarOrNaN(data_amphibole, 'K_cation_apfu');
amp.Cr = getVarOrNaN(data_amphibole, 'Cr_cation_apfu');
amp.Fe2 = amp.FeT - amp.Fe3;

% Reject Inf and finite negative values in all extracted cation variables.
fieldNames = fieldnames(amp);
for i = 1:numel(fieldNames)
    value = amp.(fieldNames{i});
    if ~isscalar(value)
        error('Amphibole variable %s must be scalar in a 1-row table.', fieldNames{i});
    end
    if isinf(value) || (isfinite(value) && value < 0)
        error('Amphibole contains an invalid negative or Inf value for %s.', fieldNames{i});
    end
end

% Check Fe3 <= total Fe only when both values are finite. NaN is retained.
if isfinite(amp.Fe3) && isfinite(amp.FeT) && amp.Fe3 > amp.FeT
    error('Amphibole contains Fe3_cation_apfu > Fe_cation_apfu.');
end

end

function plag = preparePlagioclaseRow(data_plagioclase)
% preparePlagioclaseRow
% Extract one-row plagioclase cation data without replacing NaN by zero.

if height(data_plagioclase) ~= 1
    error('Plagioclase input must be a 1-row table.');
end

plag = struct();
plag.Si = getVarOrError(data_plagioclase, 'Si_cation_apfu', 'Plagioclase');
plag.Al = getVarOrError(data_plagioclase, 'Al_cation_apfu', 'Plagioclase');
plag.Ca = getVarOrError(data_plagioclase, 'Ca_cation_apfu', 'Plagioclase');
plag.Na = getVarOrError(data_plagioclase, 'Na_cation_apfu', 'Plagioclase');
plag.K = getVarOrError(data_plagioclase, 'K_cation_apfu', 'Plagioclase');

fieldNames = fieldnames(plag);
for i = 1:numel(fieldNames)
    value = plag.(fieldNames{i});
    if ~isscalar(value)
        error('Plagioclase variable %s must be scalar in a 1-row table.', fieldNames{i});
    end
    if isinf(value) || (isfinite(value) && value < 0)
        error('Plagioclase contains an invalid negative or Inf value for %s.', fieldNames{i});
    end
end

end

function site = calcAmphiboleSites(amp)
% calcAmphiboleSites
% Calculate practical amphibole site-allocation outputs. NaN inputs are
% propagated rather than replaced by zero.

site = struct();

if isnan(amp.Si) || isnan(amp.Al)
    site.Al_4 = NaN;
    site.Al_6 = NaN;
else
    T_deficit = max(0, 8 - amp.Si);
    site.Al_4 = min(amp.Al, T_deficit);
    site.Al_6 = amp.Al - site.Al_4;
end

if isnan(amp.Ca)
    site.Ca_B = NaN;
else
    site.Ca_B = min(amp.Ca, 2);
end

if isnan(amp.Na) || isnan(site.Ca_B)
    site.Na_B = NaN;
else
    site.Na_B = min(amp.Na, max(0, 2 - site.Ca_B));
end

if isnan(amp.Na) || isnan(site.Na_B)
    site.Na_A = NaN;
else
    site.Na_A = max(0, amp.Na - site.Na_B);
end

site.K_A = amp.K;
site.A_occ = site.Na_A + site.K_A;
site.vA = 1 - site.A_occ;

end

function value = getVarOrError(tbl, variableName, mineralLabel)
% getVarOrError
% Retrieve a required scalar table variable without altering NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);
if ~isscalar(value)
    error('Variable %s must be scalar in a 1-row table.', variableName);
end

end

function value = getVarOrNaN(tbl, variableName)
% getVarOrNaN
% Retrieve an optional scalar variable. Missing optional variables are
% represented by NaN, never by zero.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);
    if ~isscalar(value)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end
else
    value = NaN;
end

end
