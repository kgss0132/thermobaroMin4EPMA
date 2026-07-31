function results = Putirka2016(rawdata_struct, P_kbar)
% functions/+thermo/+Amphibole/Putirka2016.m
% Tested with MATLAB R2024b
%
% Amphibole-only thermometers for igneous systems
% Putirka, K.D. (2016)
% American Mineralogist, 101, 841–858
% DOI: https://doi.org/10.2138/am-2016-5506
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Amphibole analysis and calculates
% temperature using the amphibole-only thermometers of Putirka (2016):
%
%   Equation (5): pressure-independent amphibole thermometer
%   Equation (6): pressure-dependent amphibole thermometer
%   Equation (8): global pressure-dependent amphibole thermometer
%
% For compatibility with the common P-T plotting routine, Equation (8) is
% copied to the standard output variables T_degreeC and T_K.
%
% The function accepts either a scalar pressure or a pressure vector. It is
% therefore compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is generated for each pressure value
% for every user-selected Amphibole analysis.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka (2016) developed these thermometers primarily for CALCIC IGNEOUS
% AMPHIBOLES. The experimental data and amphibole selection are described in
% the Method section on pp. 844–846. The amphibole-only Equations (5) and (6)
% and their calibration statistics are given on p. 850. Equation (8), which
% uses the combined global experimental data set, is given on p. 854.
%
% Important application notes:
%
%   1) Amphibole compositions must be calculated on a 23-oxygen basis.
%      Total Fe is treated as FeO and used as total Fe cations. Putirka (2016)
%      specifically does not use stoichiometrically estimated Fe3+/Fe2+ for
%      these thermometers (p. 846).
%
%   2) The calibration is intended mainly for calcic igneous amphiboles.
%      Sodic-calcic amphiboles and many Mg-Fe-Mn-Li amphiboles were excluded
%      from calibration or testing because they depart from the principal
%      calcic-amphibole compositional trends (pp. 844–845). Equation (8)
%      explicitly excludes Mg-Fe-Mn-Li amphiboles (p. 854).
%
%   3) Putirka (2016) does not define strict equation-specific minimum and
%      maximum validity limits for Equations (5), (6), and (8). The global
%      experimental database spans approximately 650–1175 degreeC and about
%      0–25 kbar (see Fig. 2 on p. 844, Fig. 5 on p. 849, and discussion of
%      the calibration data on pp. 844–850). These limits are therefore used
%      here only as APPROXIMATE DATA-COVERAGE LIMITS for non-stopping warning
%      messages, not as strict validity boundaries stated by the author.
%
%   4) Reported calibration precision is approximately:
%        Eq. (5): +/-30 degreeC for calibration data and about +/-53 degreeC
%                 for independent test data (Fig. 5 and p. 850).
%        Eq. (6): +/-28 degreeC for calibration data and about +/-52 degreeC
%                 for independent test data (Fig. 5 and p. 850).
%        Eq. (8): +/-47 degreeC for the combined global data set (p. 854).
%      Equation (8) was regressed using the combined data set and therefore
%      does not have the same independent-test structure as Equations (5)
%      and (6).
%
%   5) The restrictions T < 800 degreeC and
%      Fe#Amp = FeAmp/(FeAmp+MgAmp) < 0.65 discussed by Putirka (2016) apply
%      to amphibole BAROMETRY under restricted conditions, not to the
%      amphibole thermometers implemented here (pp. 841–842 and 852–854).
%      They are therefore not used as thermometer validity criteria.
%
%   6) A simple CaAmp >= 1.5 apfu test is retained only as an approximate
%      calcic-amphibole screening flag. Formal amphibole classification
%      requires site allocation following an accepted classification scheme.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside the approximate 0–25 kbar data coverage,
%   2) a finite calculated temperature is outside approximately
%      650–1175 degreeC,
%   3) the selected analysis fails the simple CaAmp >= 1.5 apfu screen,
%   4) a required mineral-composition input is NaN, or
%   5) a calculated temperature is NaN or Inf.
%
% NaN mineral-composition values are retained as missing values and propagate
% through the calculation. They are never replaced by zero. Finite negative
% mineral-composition values are prohibited; zero values are allowed.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The table must include the following
% variables, calculated as cations per formula unit on a 23-oxygen basis:
%
%   Si_cation_apfu
%   Ti_cation_apfu
%   Fe_cation_apfu    % total Fe cations, with total Fe treated as FeO
%   Mg_cation_apfu
%   Na_cation_apfu
%   Ca_cation_apfu
%
% Required mineral-composition values must be numeric. NaN is allowed and is
% propagated. Inf and finite negative values are prohibited. Zero is allowed.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Equation (5), pressure independent:
%   T5(degreeC) = 1781
%                 - 132.74*SiAmp
%                 + 116.6*TiAmp
%                 - 69.41*FeTAmp
%                 + 101.62*NaAmp
%
% Equation (6), pressure dependent:
%   T6(degreeC) = 1687
%                 - 118.7*SiAmp
%                 + 131.56*TiAmp
%                 - 71.41*FeTAmp
%                 + 86.13*NaAmp
%                 + 22.44*P(GPa)
%
% Equation (8), global pressure-dependent model:
%   T8(degreeC) = 1201.4
%                 - 97.93*SiAmp
%                 + 201.82*TiAmp
%                 + 72.85*MgAmp
%                 + 88.9*NaAmp
%                 + 40.65*P(GPa)
%
% Pressure is supplied in kbar and converted internally using:
%   P(GPa) = P(kbar) / 10
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2016_eq8(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing an Amphibole table (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Amphibole analysis
%
% Downstream plotting compatibility:
%   T_degreeC and T_K are standardized aliases of Equation (8), the
%   global pressure-dependent amphibole thermometer. Equation-specific
%   outputs T5_degC/T5_K, T6_degC/T6_K, and T8_degC/T8_K are retained
%   unchanged.
%

%% Input validation
if nargin < 2
    error('Putirka2016_eq8 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation dataset
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Amphibole') || ...
        ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end

dataset_amp = rawdata_struct.Amphibole;

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ca_cation_apfu'};

missingVariables = requiredVariables( ...
    ~ismember(requiredVariables, dataset_amp.Properties.VariableNames));

if ~isempty(missingVariables)
    error('Missing required Amphibole column(s): %s', ...
        strjoin(missingVariables, ', '));
end

if width(dataset_amp) < 1
    error('rawdata_struct.Amphibole must contain at least one column.');
end

if height(dataset_amp) < 1
    error('rawdata_struct.Amphibole contains no analyses.');
end

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Repeated table concatenation inside the interactive loop is avoided.
% Results are stored as table blocks in a preallocated cell buffer and are
% concatenated once after the loop finishes.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max(16, height(dataset_amp));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Approximate coverage of the experimental database used by Putirka (2016).
% These are warning limits only, not strict author-defined validity limits.
approxT_min_degC = 650;
approxT_max_degC = 1175;
approxP_min_kbar = 0;
approxP_max_kbar = 25;

pressureOutsideApproxRange = ...
    P_kbar < approxP_min_kbar | P_kbar > approxP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % Assumption: the first column stores an identifier displayed to the user.
    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, ok] = listdlg( ...
        'PromptString', ...
        'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', string(dataCodes_amp), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

    disp('=== Step 4: Checking the selected Amphibole data ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);

    % NaN values are allowed and retained. Their variable names are collected
    % here so that a non-stopping message can be printed after calculation.
    nanInputNames = findNaNInputs(selectedData_amp);

    % Required variables must be numeric. Inf and finite negative values are
    % prohibited. NaN and zero are intentionally allowed.
    validateCompositionInputs(selectedData_amp);

    disp('=== Step 5: Calculating the temperature ===');

    row = calcTemp(selectedData_amp, P_kbar, ...
        approxT_min_degC, approxT_max_degC, ...
        approxP_min_kbar, approxP_max_kbar);

    row.dataCode_amp = repmat(string(selectedCode_amp), height(row), 1);
    row = movevars(row, 'dataCode_amp', 'Before', 1);

    % Store the result as one table block. The buffer is expanded only when
    % its capacity is exhausted, rather than on every loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo calculated temperature ranges for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    disp([char(string(selectedCode_amp)) ': Eq.(5) = ' ...
        formatTemperatureRange(row.T5_degC) ' degreeC']);
    disp([char(string(selectedCode_amp)) ': Eq.(6) = ' ...
        formatTemperatureRange(row.T6_degC) ' degreeC']);
    disp([char(string(selectedCode_amp)) ': Eq.(8) = ' ...
        formatTemperatureRange(row.T8_degC) ' degreeC']);

    % Pressure warning is common to all selected amphiboles in this function
    % call, so print it only once.
    if any(pressureOutsideApproxRange) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate pressure coverage ' ...
             'of the experimental database used by Putirka (2016): about 0–25 kbar. ' ...
             '%d of %d pressure point(s) are outside this approximate range; ' ...
             'input range = %.4g–%.4g kbar. This is a data-coverage warning, ' ...
             'not a strict validity limit stated for Equations (5), (6), or (8).\n'], ...
            sum(pressureOutsideApproxRange), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % The Ca threshold is only an approximate calcic-amphibole screen.
    finiteCa = isfinite(row.CaAmp);
    if any(finiteCa & ~row.isCalcicApprox)
        fprintf(2, ...
            ['WARNING: The selected analysis %s has CaAmp < 1.5 apfu and fails ' ...
             'the simple calcic-amphibole screening used in this implementation. ' ...
             'Putirka (2016) primarily calibrated these thermometers for calcic ' ...
             'igneous amphiboles. Formal classification requires site allocation.\n'], ...
            char(string(selectedCode_amp)));
    end

    % Warn independently for each thermometer because their calculated
    % temperatures can differ.
    printTemperatureRangeWarning( ...
        row.T5_degC, 'Eq.(5)', selectedCode_amp, ...
        approxT_min_degC, approxT_max_degC);
    printTemperatureRangeWarning( ...
        row.T6_degC, 'Eq.(6)', selectedCode_amp, ...
        approxT_min_degC, approxT_max_degC);
    printTemperatureRangeWarning( ...
        row.T8_degC, 'Eq.(8)', selectedCode_amp, ...
        approxT_min_degC, approxT_max_degC);

    % NaN input warning. Calculation is intentionally not stopped.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         NaN values were retained and propagated; they were not ' ...
             'replaced by zero. Calculated temperatures may therefore be NaN.\n'], ...
            char(string(selectedCode_amp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Result-based warnings also capture Inf or NaN produced by mathematical
    % operations even when no input variable was explicitly NaN.
    printNonfiniteTemperatureWarning( ...
        row.T5_degC, 'Eq.(5)', selectedCode_amp);
    printNonfiniteTemperatureWarning( ...
        row.T6_degC, 'Eq.(6)', selectedCode_amp);
    printNonfiniteTemperatureWarning( ...
        row.T8_degC, 'Eq.(8)', selectedCode_amp);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Putirka2016_eq8', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate table blocks only once. Return an empty table when the user
% finishes without performing a calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_amphibole)
% findNaNInputs
% Return required thermometer variable names that contain NaN. The function
% does not throw an error because NaN is intentionally propagated.

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ca_cation_apfu'};

nVariables = numel(requiredVariables);
nanFlags = false(nVariables, 1);
allNames = strings(nVariables, 1);

for i = 1:nVariables
    variableName = requiredVariables{i};
    variableValue = data_amphibole.(variableName);
    allNames(i) = "Amphibole." + string(variableName);

    if isnumeric(variableValue) && any(isnan(variableValue(:)))
        nanFlags(i) = true;
    end
end

nanInputNames = allNames(nanFlags);

end

function validateCompositionInputs(data_amphibole)
% validateCompositionInputs
% Verify that all required inputs are numeric. NaN and zero are allowed.
% Inf and finite negative values stop the calculation.

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ca_cation_apfu'};

nVariables = numel(requiredVariables);
nonNumericFlags = false(nVariables, 1);
infFlags = false(nVariables, 1);
negativeFlags = false(nVariables, 1);
allNames = strings(nVariables, 1);

for i = 1:nVariables
    variableName = requiredVariables{i};
    variableValue = data_amphibole.(variableName);
    allNames(i) = "Amphibole." + string(variableName);

    if ~isnumeric(variableValue)
        nonNumericFlags(i) = true;
        continue;
    end

    if any(isinf(variableValue(:)))
        infFlags(i) = true;
    end

    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        negativeFlags(i) = true;
    end
end

if any(nonNumericFlags)
    error(['Putirka2016_eq8: required mineral-composition inputs must be numeric. ' ...
           'Non-numeric variable(s): ' ...
           char(strjoin(allNames(nonNumericFlags), ', ')) '.']);
end

if any(infFlags)
    error(['Putirka2016_eq8: Inf is not allowed in required mineral-composition ' ...
           'inputs. Inf was found in: ' ...
           char(strjoin(allNames(infFlags), ', ')) '.']);
end

if any(negativeFlags)
    error(['Putirka2016_eq8: finite negative mineral-composition values are not ' ...
           'allowed. Negative value(s) were found in: ' ...
           char(strjoin(allNames(negativeFlags), ', ')) '. ' ...
           'Zero and NaN values are allowed.']);
end

end

function row = calcTemp(data_amphibole, P_kbar, ...
        approxT_min_degC, approxT_max_degC, ...
        approxP_min_kbar, approxP_max_kbar)
% calcTemp
% Calculate Putirka (2016) amphibole-only temperatures for one selected
% amphibole analysis over a scalar or vector of pressures.
%
% Inputs:
%   data_amphibole : one-row Amphibole table
%   P_kbar         : finite non-negative pressure scalar or vector in kbar
%
% Output:
%   row : table with one row per pressure value

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

% Replicate the selected amphibole composition to match the number of
% pressure points. NaN values are intentionally retained by repmat.
SiAmp = repmat(data_amphibole.Si_cation_apfu, nP, 1);
TiAmp = repmat(data_amphibole.Ti_cation_apfu, nP, 1);
FeTAmp = repmat(data_amphibole.Fe_cation_apfu, nP, 1);
MgAmp = repmat(data_amphibole.Mg_cation_apfu, nP, 1);
NaAmp = repmat(data_amphibole.Na_cation_apfu, nP, 1);
CaAmp = repmat(data_amphibole.Ca_cation_apfu, nP, 1);

% Fe number is retained as a diagnostic output only. The Fe# < 0.65
% criterion discussed by Putirka (2016) is a restricted barometer criterion,
% not a validity criterion for Equations (5), (6), or (8).
FeNumAmp = FeTAmp ./ (FeTAmp + MgAmp);

% Putirka (2016), Equation (5): pressure-independent thermometer.
T5_degC = 1781 ...
    - 132.74 .* SiAmp ...
    + 116.6 .* TiAmp ...
    - 69.41 .* FeTAmp ...
    + 101.62 .* NaAmp;
T5_K = T5_degC + 273.15;

% Putirka (2016), Equation (6): pressure-dependent thermometer.
T6_degC = 1687 ...
    - 118.7 .* SiAmp ...
    + 131.56 .* TiAmp ...
    - 71.41 .* FeTAmp ...
    + 86.13 .* NaAmp ...
    + 22.44 .* P_GPa;
T6_K = T6_degC + 273.15;

% Putirka (2016), Equation (8): global pressure-dependent thermometer.
T8_degC = 1201.4 ...
    - 97.93 .* SiAmp ...
    + 201.82 .* TiAmp ...
    + 72.85 .* MgAmp ...
    + 88.9 .* NaAmp ...
    + 40.65 .* P_GPa;
T8_K = T8_degC + 273.15;

% Standardized downstream temperature variables. The common P-T plotting
% routines require T_degreeC and optionally T_K. Equation (8) is used as
% the representative plotted thermometer; the equation-specific variables
% T5_degC/T5_K, T6_degC/T6_K, and T8_degC/T8_K remain unchanged.
T_degreeC = T8_degC;
T_K = T8_K;

% Approximate screening and data-coverage flags. NaN values evaluate to
% false because they cannot be confirmed to lie inside a range.
isCalcicApprox = isfinite(CaAmp) & CaAmp >= 1.5;
isPWithinApproxRange = isfinite(P_kbar) & ...
    P_kbar >= approxP_min_kbar & P_kbar <= approxP_max_kbar;
isT5WithinApproxRange = isfinite(T5_degC) & ...
    T5_degC >= approxT_min_degC & T5_degC <= approxT_max_degC;
isT6WithinApproxRange = isfinite(T6_degC) & ...
    T6_degC >= approxT_min_degC & T6_degC <= approxT_max_degC;
isT8WithinApproxRange = isfinite(T8_degC) & ...
    T8_degC >= approxT_min_degC & T8_degC <= approxT_max_degC;

row = table( ...
    P_kbar, P_GPa, T_degreeC, T_K, ...
    SiAmp, TiAmp, FeTAmp, MgAmp, NaAmp, CaAmp, ...
    FeNumAmp, ...
    T5_degC, T5_K, ...
    T6_degC, T6_K, ...
    T8_degC, T8_K, ...
    isCalcicApprox, ...
    isPWithinApproxRange, ...
    isT5WithinApproxRange, ...
    isT6WithinApproxRange, ...
    isT8WithinApproxRange, ...
    'VariableNames', { ...
    'P_kbar', 'P_GPa', 'T_degreeC', 'T_K', ...
    'SiAmp', 'TiAmp', 'FeTAmp', 'MgAmp', 'NaAmp', 'CaAmp', ...
    'FeNumAmp', ...
    'T5_degC', 'T5_K', ...
    'T6_degC', 'T6_K', ...
    'T8_degC', 'T8_K', ...
    'isCalcicApprox', ...
    'isPWithinApproxRange', ...
    'isT5WithinApproxRange', ...
    'isT6WithinApproxRange', ...
    'isT8WithinApproxRange'});

end

function printTemperatureRangeWarning(temperatureValues, equationName, ...
        selectedCode, approxT_min_degC, approxT_max_degC)
% printTemperatureRangeWarning
% Print a non-stopping warning when a finite calculated temperature lies
% outside the approximate experimental data coverage.

finiteTemperature = isfinite(temperatureValues);
outsideApproxRange = finiteTemperature & ...
    (temperatureValues < approxT_min_degC | ...
     temperatureValues > approxT_max_degC);

if any(outsideApproxRange)
    finiteValues = temperatureValues(finiteTemperature);
    fprintf(2, ...
        ['WARNING: %s calculated temperature for %s is outside the approximate ' ...
         'temperature coverage of the experimental database used by Putirka ' ...
         '(2016): about %.0f–%.0f degreeC. %d of %d finite temperature point(s) ' ...
         'are outside this approximate range; calculated finite range = ' ...
         '%.4g–%.4g degreeC. This is a data-coverage warning, not a strict ' ...
         'validity limit stated for this equation.\n'], ...
        equationName, ...
        char(string(selectedCode)), ...
        approxT_min_degC, ...
        approxT_max_degC, ...
        sum(outsideApproxRange), ...
        sum(finiteTemperature), ...
        min(finiteValues), ...
        max(finiteValues));
end

end

function printNonfiniteTemperatureWarning(temperatureValues, equationName, ...
        selectedCode)
% printNonfiniteTemperatureWarning
% Print a non-stopping warning for NaN or Inf calculated temperatures.

invalidTemperature = ~isfinite(temperatureValues);

if any(invalidTemperature)
    fprintf(2, ...
        ['WARNING: Non-finite temperature value(s) were calculated by %s for %s ' ...
         '(%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the calculation ' ...
         'has not been stopped.\n'], ...
        equationName, ...
        char(string(selectedCode)), ...
        sum(invalidTemperature), ...
        numel(temperatureValues), ...
        sum(isnan(temperatureValues)), ...
        sum(isinf(temperatureValues)));
end

end

function rangeText = formatTemperatureRange(temperatureValues)
% formatTemperatureRange
% Format one or more temperature values for display without changing data.

if isscalar(temperatureValues)
    rangeText = num2str(temperatureValues, '%.1f');
    return;
end

finiteValues = temperatureValues(isfinite(temperatureValues));

if isempty(finiteValues)
    rangeText = 'NaN to NaN';
else
    rangeText = [num2str(min(finiteValues), '%.1f') ' to ' ...
        num2str(max(finiteValues), '%.1f')];

    if any(~isfinite(temperatureValues))
        rangeText = [rangeText ' (including non-finite value(s))'];
    end
end

end
