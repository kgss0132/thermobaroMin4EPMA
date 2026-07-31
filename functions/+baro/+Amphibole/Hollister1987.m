function results = Hollister1987(rawdata_struct, T_degreeC)
% functions/+baro/+Amphibole/Hollister1987.m
% Tested with MATLAB R2024b
%
% Aluminum-in-hornblende empirical igneous geobarometer
% Hollister, L.S., Grissom, G.C., Peters, E.K., Stowell, H.H.,
% and Sisson, V.B. (1987)
% American Mineralogist, 72, 231-239
% DOI: XXXXXXX
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Amphibole analysis and calculates
% pressure using the total-Al-in-hornblende empirical igneous geobarometer
% of Hollister et al. (1987).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Amphibole analysis, one output row
% is returned for every input temperature value.
%
% Temperature does not enter the published pressure equation. T_degreeC is
% carried into the output and is used only for application-range screening.
%
% The function is designed for repeated calculations. After each run, it
% asks whether another Amphibole analysis should be calculated and stores
% all result blocks in one output table.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Hollister et al. (1987) tested and refined the Hammarstrom and Zen (1986)
% total-Al-in-hornblende pressure relation. Their preferred regression is:
%
%   P(kbar) = -4.76 + 5.64 * AlT
%
% where AlT is total Al in hornblende. The regression is given in the
% Discussion on p. 238.
%
% Pressure calibration and uncertainty:
%   Pressure range : 2-8 kbar
%   Uncertainty    : approximately +/-1 kbar, provided the sampling and
%                    petrographic conditions described below are followed
%
% The 2-8 kbar range and +/-1 kbar estimate are stated in the abstract on
% p. 231 and reiterated in the Discussion on p. 238. The study particularly
% strengthens the intermediate-pressure interval of approximately 4-6 kbar.
%
% Temperature context:
%   Amphibole in calc-alkaline compositions may crystallize broadly between
%   approximately 650 and 950 degreeC (p. 232). This interval is NOT a
%   direct experimental calibration range for the pressure equation.
%   Hollister et al. (1987) argue that, above approximately 2 kbar, the
%   final solidification temperature of relevant calc-alkaline plutons
%   varies by no more than about 100 degreeC. Below approximately 2 kbar,
%   temperature effects on hornblende Al may become too large to separate
%   reliably from pressure effects (pp. 232 and 238).
%
% Required petrographic and mineralogical conditions:
%   1) Quartz, plagioclase, hornblende, biotite, orthoclase (K-feldspar),
%      titanite, and magnetite must have crystallized together from a melt
%      (pp. 232-233).
%   2) Only unaltered hornblende RIM compositions should be used because
%      these are the portions most likely to have equilibrated with the
%      final residual melt and the complete mineral assemblage (p. 233;
%      Discussion, p. 238).
%   3) The pressure should be above approximately 2 kbar (pp. 232-233).
%   4) Coexisting plagioclase rim compositions in the calibration were
%      approximately An25-An35 (p. 233).
%   5) Multiple rim analyses adjacent to different minerals should be
%      obtained and averaged. The calibration used averaged rim data rather
%      than a single spot analysis (Analytical methods, pp. 233-234;
%      Discussion, p. 238).
%
% Important cautions:
%   - Orthoclase-deficient rocks may contain hornblende with anomalously high
%     AlT and may therefore yield pressure estimates more than approximately
%     1 kbar too high (discussion of the Carlson Creek pluton, pp. 236-237).
%   - Using hornblende cores instead of rims may introduce an additional
%     error of approximately +/-0.5 kbar (p. 238).
%   - A single named pluton may contain multiple intrusive pulses emplaced
%     over a significant pressure range. Field relations and intrusive
%     history must therefore be considered (Ponder pluton discussion,
%     p. 236).
%   - The numerical flags produced by this function cannot evaluate rim
%     position, alteration, mineral coexistence, plagioclase An content,
%     orthoclase presence, or whether multiple rim analyses were averaged.
%     These conditions must be assessed independently by the user.
%
% Supplementary amphibole screening:
%   The original implementation retained Si < 7.5 and Ca > 1.6 apfu as a
%   simple igneous calcic-amphibole screen inherited from Hammarstrom and
%   Zen (1986). These two numerical checks are retained as diagnostic output
%   for backward compatibility, but they are not sufficient to establish
%   applicability of the Hollister et al. (1987) calibration.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) a finite input temperature is outside the broad 650-950 degreeC
%      crystallization interval discussed on p. 232,
%   2) a finite calculated pressure is outside the 2-8 kbar calibration
%      range,
%   3) the supplementary Si-Ca chemical screen is not satisfied,
%   4) a required calculation or screening input contains NaN,
%   5) the calculated pressure is NaN or Inf, or
%   6) a negative finite pressure is calculated.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include:
%
%   Required Amphibole variables:
%     Si_cation_apfu
%     Ti_cation_apfu
%     Al_cation_apfu          % AlT used by the pressure equation
%     Fe_cation_apfu          % total Fe in the input cation table
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
% Finite cation values must be non-negative. NaN is allowed and retained as
% missing; it is never replaced by zero. Inf and finite negative values are
% prohibited. A NaN Al_cation_apfu value propagates directly to P_kbar.
%
% This barometer does not use a liquid (Liq) composition. Therefore, rules
% concerning exclusion of F and Cl from cationTotal_liq and their exclusion
% from Liq NaN warnings are not applicable to this function.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
% 1) Amphibole Al-site quantities retained for output:
%     AlT  = total Al cations
%     AlIV = max(0, 8 - Si)
%     AlVI = max(0, AlT - AlIV)
%
% 2) Pressure:
%     P(kbar) = 5.64 * AlT - 4.76
%
% Notes:
% - AlT is the only mineral-composition variable in the pressure equation.
% - T_degreeC does not enter the pressure equation.
% - Negative calculated pressure is retained for diagnostic purposes and
%   reported by a non-stopping warning.
% - NaN is retained and propagated rather than interpreted as zero.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Hollister1987(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing an Amphibole table
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per input temperature value for every
%             user-selected Amphibole analysis
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('Hollister1987 requires (rawdata_struct, T_degreeC).');
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

%% 1) Retrieve cation dataset
% Extract the required Amphibole table. The source table is not modified.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Amphibole') || ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end

dataset_amp = rawdata_struct.Amphibole;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-analysis result in a cell buffer and concatenate once
% after the interactive loop. This avoids repeated growth of the full output
% table during every calculation.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Published pressure calibration range.
calibrationP_min_kbar = 2;
calibrationP_max_kbar = 8;

% Broad amphibole crystallization interval discussed by Hollister et al.
% (1987). This is contextual and not a direct pressure calibration range.
contextT_min_degreeC = 650;
contextT_max_degreeC = 950;

temperatureOutsideContext = isfinite(T_degreeC) & ...
    (T_degreeC < contextT_min_degreeC | ...
     T_degreeC > contextT_max_degreeC);

temperatureWarningIssued = false;
applicationCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % ----- Amphibole selection -----
    % The first table column is used only as the displayed identifier.
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

    % ----- Calculation -----
    disp('=== Step 4: Calculating the pressure ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);

    % Reject Inf and finite negative cation values. NaN is allowed and is
    % retained for propagation through the calculation.
    validateNonNegativeInputs(selectedData_amp);

    % Identify NaN values in the pressure equation, temperature input, and
    % supplementary numerical applicability screen.
    nanInputNames = findNaNInputs(selectedData_amp, T_degreeC);

    row = calcPressure(selectedData_amp, T_degreeC);

    % Repeat the selected identifier for every temperature row.
    row.dataCode_amphibole = ...
        repmat(string(selectedCode_amp), height(row), 1);
    row = movevars(row, 'dataCode_amphibole', 'Before', 1);

    % Store one result block per selected analysis. Expand the buffer only
    % when its current capacity is exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated result for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    finitePressure = isfinite(row.P_kbar);
    finiteTemperature = isfinite(row.T_degreeC);

    if height(row) == 1
        if finitePressure
            fprintf('%s: P = %.6g kbar at T = %.6g degreeC\n', ...
                char(string(selectedCode_amp)), ...
                row.P_kbar, row.T_degreeC);
        else
            fprintf('%s: P = %s kbar at T = %.6g degreeC\n', ...
                char(string(selectedCode_amp)), ...
                char(string(row.P_kbar)), row.T_degreeC);
        end
    else
        if any(finitePressure)
            finitePressureValues = row.P_kbar(finitePressure);
            if any(finiteTemperature)
                finiteTemperatureValues = row.T_degreeC(finiteTemperature);
                fprintf('%s: P = %.6g to %.6g kbar for finite T = %.6g to %.6g degreeC\n', ...
                    char(string(selectedCode_amp)), ...
                    min(finitePressureValues), max(finitePressureValues), ...
                    min(finiteTemperatureValues), max(finiteTemperatureValues));
            else
                fprintf('%s: P = %.6g to %.6g kbar; all %d input temperature values are NaN\n', ...
                    char(string(selectedCode_amp)), ...
                    min(finitePressureValues), max(finitePressureValues), ...
                    height(row));
            end
        else
            if any(finiteTemperature)
                finiteTemperatureValues = row.T_degreeC(finiteTemperature);
                fprintf('%s: P = NaN/Inf for all %d temperature point(s), finite T = %.6g to %.6g degreeC\n', ...
                    char(string(selectedCode_amp)), height(row), ...
                    min(finiteTemperatureValues), max(finiteTemperatureValues));
            else
                fprintf('%s: P = NaN/Inf for all %d point(s); all input temperature values are NaN\n', ...
                    char(string(selectedCode_amp)), height(row));
            end
        end
    end

    % Print the major non-numerical application limitations once per call.
    if ~applicationCautionIssued
        fprintf(2, ...
            ['CAUTION: Hollister et al. (1987) calibrated this equation using ' ...
             'averaged, unaltered hornblende rim analyses from calc-alkaline ' ...
             'plutons containing quartz + plagioclase + hornblende + biotite + ' ...
             'orthoclase + titanite + magnetite (pp. 232-233, 238). The program ' ...
             'cannot verify rim position, alteration, phase coexistence, ' ...
             'plagioclase An25-An35, orthoclase presence, or averaging of ' ...
             'multiple rim analyses.\n']);
        applicationCautionIssued = true;
    end

    % Input temperature is common to all selected analyses in the current
    % function call, so this warning is printed only once.
    if any(temperatureOutsideContext) && ~temperatureWarningIssued
        finiteTemperatureValues = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the broad approximately ' ...
             '650-950 degreeC amphibole crystallization interval discussed by ' ...
             'Hollister et al. (1987; p. 232). This interval is not a direct ' ...
             'experimental calibration range for the pressure equation. ' ...
             '%d of %d finite temperature point(s) are outside the interval; ' ...
             'input finite range = %.6g-%.6g degreeC.\n'], ...
            sum(temperatureOutsideContext), ...
            numel(finiteTemperatureValues), ...
            min(finiteTemperatureValues), ...
            max(finiteTemperatureValues));
        temperatureWarningIssued = true;
    end

    % Warn when finite pressure results fall outside the published 2-8 kbar
    % calibration range.
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the Hollister et al. ' ...
             '(1987) calibration range of 2-8 kbar (pp. 231 and 238). ' ...
             '%d of %d finite pressure point(s) are outside the range; ' ...
             'calculated finite range = %.6g-%.6g kbar for %s. ' ...
             'The values are retained in the output table.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_amp)));
    end

    % Retain the original Si-Ca screen as a supplementary diagnostic. This
    % screen alone does not establish applicability of the 1987 calibration.
    if ~row.isIgneousChemicalScreen(1)
        fprintf(2, ...
            ['WARNING: Amphibole %s does not satisfy the supplementary ' ...
             'igneous calcic-amphibole screen retained from the original ' ...
             'implementation: Si < 7.5 and Ca > 1.6 apfu. ' ...
             'Si = %.6g, Ca = %.6g. The calculation was continued. ' ...
             'This screen is not a substitute for the petrographic conditions ' ...
             'specified by Hollister et al. (1987).\n'], ...
            char(string(selectedCode_amp)), ...
            row.Si_amp(1), row.Ca_amp(1));
    end

    % List exact required input names containing NaN. Temperature-vector NaN
    % positions are reported by index.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s: %s.\n' ...
             '         NaN values were retained and were not replaced by zero; ' ...
             'the pressure and/or applicability fields may remain NaN or false.\n'], ...
            char(string(selectedCode_amp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all NaN/Inf pressure results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_amp)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressure is outside the calibration but is retained for
    % diagnosis rather than silently replaced by NaN or zero.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s ' ...
             '(%d of %d points). The value was retained for diagnostic ' ...
             'purposes and is outside the 2-8 kbar calibration range.\n'], ...
            char(string(selectedCode_amp)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another amphibole analysis.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Hollister1987', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered blocks once. Return an empty table when the user
% canceled before completing any calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_amphibole, T_degreeC)
% findNaNInputs
% Return names of NaN inputs used by the pressure equation, temperature
% screening, or supplementary numerical applicability checks. NaN values do
% not cause an error and are never replaced by zero.

maxNames = 4;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

variableNames = { ...
    'Al_cation_apfu', ...
    'Si_cation_apfu', ...
    'Ca_cation_apfu'};

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    variableValue = data_amphibole.(variableName);

    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = ...
            "Amphibole." + string(variableName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_amphibole)
% validateNonNegativeInputs
% Reject Inf and finite negative values in all amphibole cation variables
% read by this implementation. NaN and zero are allowed and retained.

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu'};

optionalVariables = { ...
    'Fe3_cation_apfu', ...
    'Mn_cation_apfu', ...
    'K_cation_apfu', ...
    'Cr_cation_apfu'};

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        error('Amphibole table must contain variable: %s', variableName);
    end
end

allVariables = [requiredVariables, optionalVariables];
invalidInputBuffer = strings(numel(allVariables), 1);
nInvalidInputs = 0;

for i = 1:numel(allVariables)
    variableName = allVariables{i};

    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        continue;
    end

    variableValue = data_amphibole.(variableName);

    if ~isscalar(variableValue)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end

    if any(isinf(variableValue(:)) | ...
            (isfinite(variableValue(:)) & variableValue(:) < 0))
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Amphibole." + string(variableName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['Hollister1987: amphibole cation values must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_amphibole, T_degreeC)
% calcPressure
% Compute pressure for one amphibole row at one or more input temperatures.
% Temperature does not enter the published pressure equation. NaN values
% propagate through all calculations in which they are used.
%
% Inputs:
%   data_amphibole : 1-row Amphibole table
%   T_degreeC      : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

amp = prepareAmphiboleRow(data_amphibole);

% Preserve NaN explicitly in derived Al-site quantities.
AlT_scalar = amp.Al;

AlIV_scalar = 8 - amp.Si;
if isfinite(AlIV_scalar) && AlIV_scalar < 0
    AlIV_scalar = 0;
end

AlVI_scalar = AlT_scalar - AlIV_scalar;
if isfinite(AlVI_scalar) && AlVI_scalar < 0
    AlVI_scalar = 0;
end

% Hollister et al. (1987) regression.
P_kbar_scalar = 5.64 .* AlT_scalar - 4.76;

% Numerical screening and range flags.
isIgneousChemicalScreen_scalar = ...
    isfinite(amp.Si) && isfinite(amp.Ca) && ...
    amp.Si < 7.5 && amp.Ca > 1.6;

isWithinContextTRange = ...
    isfinite(T_degreeC) & ...
    T_degreeC >= 650 & T_degreeC <= 950;

isWithinCalibrationP_scalar = ...
    isfinite(P_kbar_scalar) && ...
    P_kbar_scalar >= 2 && P_kbar_scalar <= 8;

% This is only a numeric diagnostic flag. Full application conditions must
% be assessed from petrography, field relations, and plagioclase data.
isApplicable_numeric = ...
    isWithinContextTRange & ...
    repmat(isWithinCalibrationP_scalar, nT, 1) & ...
    repmat(isIgneousChemicalScreen_scalar, nT, 1);

% Expand scalar composition and pressure outputs to the temperature-vector
% length. All table variables are created at their final size.
row = table();

row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.Si_amp = repmat(amp.Si, nT, 1);
row.Ti_amp = repmat(amp.Ti, nT, 1);
row.AlT_amp = repmat(AlT_scalar, nT, 1);
row.Al_tot_amp = row.AlT_amp;
row.AlIV_amp = repmat(AlIV_scalar, nT, 1);
row.Al_IV_amp = row.AlIV_amp;
row.AlVI_amp = repmat(AlVI_scalar, nT, 1);
row.Al_VI_amp = row.AlVI_amp;

row.FeT_amp = repmat(amp.FeT, nT, 1);
row.Fe2_amp = repmat(amp.Fe2, nT, 1);
row.Fe3_amp = repmat(amp.Fe3, nT, 1);
row.Mg_amp = repmat(amp.Mg, nT, 1);
row.Ca_amp = repmat(amp.Ca, nT, 1);
row.Na_amp = repmat(amp.Na, nT, 1);
row.K_amp = repmat(amp.K, nT, 1);
row.Mn_amp = repmat(amp.Mn, nT, 1);
row.Cr_amp = repmat(amp.Cr, nT, 1);

row.P_kbar = repmat(P_kbar_scalar, nT, 1);
row.P_uncertainty_kbar = ones(nT, 1);
row.AlT_standard_deviation_apfu = repmat(0.13, nT, 1);

row.isWithinContextTRange = isWithinContextTRange;
row.isWithinCalibrationP = ...
    repmat(isWithinCalibrationP_scalar, nT, 1);
row.isIgneousChemicalScreen = ...
    repmat(isIgneousChemicalScreen_scalar, nT, 1);
row.isApplicable_numeric = isApplicable_numeric;

% Backward-compatible output name. This remains a numerical diagnostic and
% does not verify the full petrographic application conditions.
row.isApplicable = isApplicable_numeric;

end

function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Extract one-row amphibole cation data. Required variables must exist.
% Missing optional variables are represented as NaN, never as zero.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();

amp.Si = getVarOrError(data_amphibole, ...
    'Si_cation_apfu', 'Amphibole');
amp.Ti = getVarOrError(data_amphibole, ...
    'Ti_cation_apfu', 'Amphibole');
amp.Al = getVarOrError(data_amphibole, ...
    'Al_cation_apfu', 'Amphibole');
amp.FeT = getVarOrError(data_amphibole, ...
    'Fe_cation_apfu', 'Amphibole');
amp.Mg = getVarOrError(data_amphibole, ...
    'Mg_cation_apfu', 'Amphibole');
amp.Ca = getVarOrError(data_amphibole, ...
    'Ca_cation_apfu', 'Amphibole');
amp.Na = getVarOrError(data_amphibole, ...
    'Na_cation_apfu', 'Amphibole');

amp.Fe3 = getVarOrNaN(data_amphibole, 'Fe3_cation_apfu');
amp.Mn = getVarOrNaN(data_amphibole, 'Mn_cation_apfu');
amp.K = getVarOrNaN(data_amphibole, 'K_cation_apfu');
amp.Cr = getVarOrNaN(data_amphibole, 'Cr_cation_apfu');

% Fe2 remains NaN when total Fe or Fe3 is NaN. Missing Fe3 is not interpreted
% as zero.
if isfinite(amp.Fe3) && isfinite(amp.FeT) && amp.Fe3 > amp.FeT
    error('Amphibole contains Fe3_cation_apfu > Fe_cation_apfu.');
end

amp.Fe2 = amp.FeT - amp.Fe3;

% A negative derived Fe2 value is physically invalid. This can only occur
% when both total Fe and Fe3 are finite and inconsistent.
if isfinite(amp.Fe2) && amp.Fe2 < 0
    error('Amphibole contains a negative derived Fe2 value.');
end

end

function value = getVarOrError(tbl, variableName, mineralLabel)
% getVarOrError
% Retrieve one required scalar variable without altering NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);

if ~isscalar(value)
    error('Variable %s must be scalar in a 1-row table.', variableName);
end

if isinf(value) || (isfinite(value) && value < 0)
    error('Variable %s must be non-negative or NaN.', variableName);
end

end

function value = getVarOrNaN(tbl, variableName)
% getVarOrNaN
% Retrieve one optional scalar variable. Missing variables and explicit NaN
% values remain NaN and are never interpreted as zero.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);

    if ~isscalar(value)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end

    if isinf(value) || (isfinite(value) && value < 0)
        error('Variable %s must be non-negative or NaN.', variableName);
    end
else
    value = NaN;
end

end
