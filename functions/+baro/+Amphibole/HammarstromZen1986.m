function results = HammarstromZen1986(rawdata_struct, T_degreeC)
% functions/+baro/+Amphibole/HammarstromZen1986.m
% Tested with MATLAB R2024b
%
% Aluminum-in-hornblende empirical igneous geobarometer
% Hammarstrom, J.M. and Zen, E-an. (1986)
% American Mineralogist, 71, 1297-1313
% DOI: None
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Amphibole analysis and calculates
% pressure using the total-Al-in-hornblende empirical igneous geobarometer
% of Hammarstrom and Zen (1986).
%
% The pressure equation does not contain a temperature term. T_degreeC is
% nevertheless accepted as a finite, non-negative scalar or vector so that
% the function can be called consistently from both startBaroCalc_fixedT
% and startBaroCalc_rangeT. One output row is returned for each input
% temperature value, with the same pressure repeated for the selected
% amphibole analysis.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Amphibole analysis, and appends all
% result blocks into a single output table.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% The empirical regression is:
%
%   P(kbar) = -3.92 + 5.03 * AlT
%
% where AlT is total Al in hornblende, expressed as cations per 23 oxygens.
% The preferred linear regression, its residual error, and limitations are
% discussed on pp. 1306-1308; the equation is listed in Table 6 on p. 1308.
%
% Recommended application:
%   - Primary magmatic calcic hornblende in quartz-saturated, calc-alkaline
%     plutonic rocks.
%   - The recommended buffering assemblage is plagioclase (andesine to
%     oligoclase) + K-feldspar + quartz + hornblende + biotite + sphene +
%     magnetite or ilmenite (application discussion, p. 1310).
%   - Hornblende used for pressure estimation should coexist with quartz.
%     Lack of quartz equilibrium may produce anomalously high pressure
%     estimates because silica activity also affects hornblende Al content
%     (pp. 1303-1304 and p. 1310).
%   - Texturally late, altered, actinolitized, oxidized, or clearly
%     subsolidus-reequilibrated amphibole should be excluded. Late-stage
%     oxidation can lower hornblende Al and yield anomalously low pressure
%     estimates (p. 1310).
%   - Multiple analyses should be used to evaluate core-rim variation,
%     grain-to-grain variation, adjacent-mineral effects, and alteration
%     (p. 1310). A single analysis should not automatically be interpreted
%     as the emplacement pressure of an entire pluton.
%
% Chemical screening used by Hammarstrom and Zen (1986):
%   Si < 7.5 cations per 23 oxygens
%   Ca > 1.6 cations per 23 oxygens
% These filters were used to screen the igneous calcic-amphibole dataset
% and are described on p. 1298. They are necessary but are not, by
% themselves, sufficient proof that the geobarometer is applicable.
%
% Nominal range used in this implementation:
%   Temperature : approximately 700-900 degreeC
%                 This is the normal hydrous tonalite-granodiorite magmatic
%                 temperature interval discussed on p. 1299. Temperature is
%                 used here only for applicability warnings and is not part
%                 of the pressure equation.
%   Pressure    : approximately 0-8 kbar
%                 The natural plutonic calibration is dominated by shallow
%                 complexes near 1-2 kbar and deeper complexes near 8 kbar
%                 (pp. 1299-1300). The authors prefer the linear relation for
%                 low- to moderate-pressure regimes (pp. 1306-1308).
%   High P note : results above 8 kbar require caution. The paper states that
%                 the regression may fail above approximately 10 kbar or for
%                 inappropriate mineral assemblages (pp. 1308-1310).
%
% The linear regression predicts negative pressures when AlT is below about
% 0.779 apfu. Such values are retained in the output and reported by a
% non-stopping warning rather than silently replaced by zero or NaN.
%
% Suggested uncertainty:
%   approximately +/- 3 kbar (abstract, p. 1297).
%   The residual error of the preferred linear regression is 1.86 kbar
%   (Table 6, p. 1308), but +/- 3 kbar is retained as the practical
%   uncertainty recommended in the paper.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The table must include the following
% normalized-cation variables:
%
%   Required for pressure calculation and/or applicability screening:
%     Si_cation_apfu
%     Al_cation_apfu
%     Ca_cation_apfu
%
%   Required for the stable output variable set:
%     Ti_cation_apfu
%     Fe_cation_apfu
%     Mg_cation_apfu
%     Na_cation_apfu
%
%   Optional output variables; missing values are retained as NaN, not zero:
%     Fe3_cation_apfu
%     Mn_cation_apfu
%     K_cation_apfu
%     Cr_cation_apfu
%
% All finite amphibole-cation values read by this function must be
% non-negative. Negative finite values stop the calculation. NaN values are
% retained, propagated through calculations that use them, and reported by
% non-stopping fprintf warnings. Inf values are rejected.
%
% This barometer does not use a liquid (Liq) composition. Therefore, rules
% concerning F and Cl exclusion from cationTotal_liq are not applicable to
% this implementation.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
% 1) Amphibole Al partitioning for output and screening
%     AlT  = total Al cations per 23 oxygens
%     AlIV = max(0, 8 - Si)
%     AlVI = max(0, AlT - AlIV)
%
% 2) Pressure
%     P(kbar) = 5.03 * AlT - 3.92
%
% Notes:
% - Pressure depends only on AlT in this empirical equation.
% - T_degreeC is carried into the output and used for range warnings only.
% - NaN is never replaced by zero.
% - Negative calculated pressure is retained and flagged as outside the
%   meaningful lower range of the calibration.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = HammarstromZen1986(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing an Amphibole table (see above)
%   T_degreeC      : temperature in degreeC; finite, non-negative numeric
%                    scalar or vector
%
% Output:
%   results : table containing one row per input temperature value for every
%             user-selected Amphibole analysis
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or an
% invalid temperature vector.
if nargin < 2
    error('HammarstromZen1986 requires (rawdata_struct, T_degreeC).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(T_degreeC) || isempty(T_degreeC) || ~isvector(T_degreeC) || ...
        any(~isfinite(T_degreeC(:))) || any(T_degreeC(:) < 0)
    error('T_degreeC must be a finite non-negative numeric scalar or vector.');
end

T_degreeC = T_degreeC(:);

%% 1) Retrieve cation dataset
% Extract the required Amphibole table from the input struct.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Amphibole') || ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end

dataset_amp = rawdata_struct.Amphibole;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each selected-analysis result is stored temporarily as one table block.
% Repeated concatenation of the complete results table inside the interactive
% loop is avoided. The buffer is expanded only when its capacity is reached,
% and all result blocks are concatenated once after the loop.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Nominal application limits summarized from Hammarstrom and Zen (1986).
calibrationT_min_degC = 700;
calibrationT_max_degC = 900;
calibrationP_min_kbar = 0;
calibrationP_max_kbar = 8;
highPressureCaution_kbar = 10;

% Temperature is common to all selected amphibole analyses in this call, so
% an out-of-range input-temperature warning is printed only once.
temperatureOutsideCalibration = ...
    T_degreeC < calibrationT_min_degC | ...
    T_degreeC > calibrationT_max_degC;
temperatureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
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

    % ----- Calculation -----
    disp('=== Step 4: Calculating the pressure ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);

    % Check NaN only in variables that directly affect pressure,
    % Al-site calculations, or the published chemical screening.
    nanInputNames = findNaNInputs(selectedData_amp);

    % Reject negative finite cation values and Inf values. NaN is allowed so
    % that it remains missing and propagates through relevant calculations.
    validateNonNegativeInputs(selectedData_amp);

    row = calcPressure(selectedData_amp, T_degreeC);

    % Store the selected identifier for every temperature row.
    row.dataCode_amphibole = ...
        repmat(string(selectedCode_amp), height(row), 1);
    row = movevars(row, 'dataCode_amphibole', 'Before', 1);

    % Store this result as one buffered table block.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the calculated pressure for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    finitePressure = isfinite(row.P_kbar);
    if height(row) == 1
        if finitePressure
            fprintf('%s: P = %.6g kbar at T = %.6g degreeC\n', ...
                char(string(selectedCode_amp)), row.P_kbar, row.T_degreeC);
        else
            fprintf('%s: P = %s kbar at T = %.6g degreeC\n', ...
                char(string(selectedCode_amp)), ...
                char(string(row.P_kbar)), row.T_degreeC);
        end
    else
        if any(finitePressure)
            finiteValues = row.P_kbar(finitePressure);
            fprintf('%s: P = %.6g to %.6g kbar for T = %.6g to %.6g degreeC\n', ...
                char(string(selectedCode_amp)), ...
                min(finiteValues), max(finiteValues), ...
                row.T_degreeC(1), row.T_degreeC(end));
        else
            fprintf(['%s: P = NaN/Inf for all %d temperature point(s), ' ...
                    'T = %.6g to %.6g degreeC\n'], ...
                char(string(selectedCode_amp)), height(row), ...
                row.T_degreeC(1), row.T_degreeC(end));
        end
    end

    % Warn once when any input temperature lies outside the nominal hydrous
    % magmatic range. The equation itself does not use temperature.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        fprintf(2, ...
            ['WARNING: Input temperature is outside the nominal application range ' ...
             'summarized from Hammarstrom and Zen (1986): 700-900 degreeC. ' ...
             '%d of %d temperature point(s) are outside the range; ' ...
             'input range = %.6g-%.6g degreeC. Temperature is used only for ' ...
             'applicability screening and does not enter the pressure equation.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(T_degreeC), ...
            min(T_degreeC), ...
            max(T_degreeC));
        temperatureWarningIssued = true;
    end

    % Warn when the calculated pressure is outside the nominal 0-8 kbar
    % natural-pluton calibration range. Non-finite values are handled below.
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finiteValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the nominal calibration ' ...
             'range of Hammarstrom and Zen (1986): approximately 0-8 kbar. ' ...
             '%d of %d finite pressure point(s) are outside the range; ' ...
             'calculated finite range = %.6g-%.6g kbar for %s. ' ...
             'The values are retained in the output table.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_amp)));
    end

    % Give a stronger warning for the high-pressure region where the paper
    % explicitly cautions that the regression may fail.
    pressureAboveHighCaution = finitePressure & ...
        row.P_kbar > highPressureCaution_kbar;

    if any(pressureAboveHighCaution)
        fprintf(2, ...
            ['WARNING: Calculated pressure exceeds 10 kbar for %s ' ...
             '(%d of %d finite point(s)). Hammarstrom and Zen (1986) state ' ...
             'that the regression may fail above approximately 10 kbar or ' ...
             'for inappropriate mineral assemblages.\n'], ...
            char(string(selectedCode_amp)), ...
            sum(pressureAboveHighCaution), ...
            sum(finitePressure));
    end

    % Report failure of the two chemical screening criteria used to select
    % igneous calcic amphiboles in the original calibration dataset.
    if ~row.isIgneousChemicalScreen(1)
        fprintf(2, ...
            ['WARNING: Amphibole %s does not satisfy the complete chemical ' ...
             'screening criteria used by Hammarstrom and Zen (1986): ' ...
             'Si < 7.5 and Ca > 1.6 cations per 23 oxygens. ' ...
             'Si = %.6g, Ca = %.6g. The calculation was continued.\n'], ...
            char(string(selectedCode_amp)), ...
            row.Si_amp(1), ...
            row.Ca_amp(1));
    end

    % Print a non-stopping warning that identifies every NaN input used by
    % the pressure calculation, Al partitioning, or applicability screening.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s: %s.\n' ...
             '         NaN was retained rather than replaced by zero; ' ...
             'the calculated pressure and/or applicability fields may be NaN or false.\n'], ...
            char(string(selectedCode_amp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report NaN/Inf pressure results without stopping.
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

    disp('--------------------------------------------------');

    % Ask whether to repeat with another amphibole analysis.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'HammarstromZen1986', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered blocks once. Return an empty table when no
% calculation was completed.
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
% Return names of NaN variables that directly affect pressure, Al-site
% calculations, or the published igneous-amphibole screening criteria.
% NaN values do not cause an error in this function.

variableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Ca_cation_apfu'};

nanFlags = false(numel(variableNames), 1);

for i = 1:numel(variableNames)
    variableName = variableNames{i};

    if ismember(variableName, data_amphibole.Properties.VariableNames)
        variableValue = data_amphibole.(variableName);
        nanFlags(i) = any(isnan(variableValue(:)));
    end
end

nanInputNames = "Amphibole." + string(variableNames(nanFlags)).';

end

function validateNonNegativeInputs(data_amphibole)
% validateNonNegativeInputs
% Reject negative finite values and Inf values in all amphibole-cation
% variables read by this implementation. NaN is intentionally allowed and
% is retained for propagation and warning messages.

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
invalidNames = strings(numel(allVariables), 1);
nInvalid = 0;

for i = 1:numel(allVariables)
    variableName = allVariables{i};

    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        continue;
    end

    variableValue = data_amphibole.(variableName);

    if ~isscalar(variableValue)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end

    if any(isinf(variableValue(:))) || ...
            any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidNames(nInvalid) = "Amphibole." + string(variableName);
    end
end

if nInvalid > 0
    invalidNames = invalidNames(1:nInvalid);
    error(['HammarstromZen1986: amphibole-cation values must be ' ...
           'non-negative or NaN. Negative finite or Inf value(s) were ' ...
           'found in: ' char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcPressure(data_amphibole, T_degreeC)
% calcPressure
% Compute pressure for one amphibole analysis and return one output row per
% input temperature. Temperature does not enter the pressure equation.
%
% Inputs:
%   data_amphibole : 1-row Amphibole table
%   T_degreeC      : finite, non-negative scalar or column vector
%
% Output:
%   row : table with one row per temperature value

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);

amp = prepareAmphiboleRow(data_amphibole);

% Preserve NaN explicitly. The conditional clamping is applied only to
% finite derived Al-site values.
AlT_scalar = amp.Al;
AlIV_scalar = 8 - amp.Si;
if isfinite(AlIV_scalar) && AlIV_scalar < 0
    AlIV_scalar = 0;
end

AlVI_scalar = AlT_scalar - AlIV_scalar;
if isfinite(AlVI_scalar) && AlVI_scalar < 0
    AlVI_scalar = 0;
end

% Preferred linear regression of Hammarstrom and Zen (1986).
P_kbar_scalar = 5.03 .* AlT_scalar - 3.92;

% Chemical screening criteria used in the original calibration dataset.
isIgneousChemicalScreen_scalar = ...
    isfinite(amp.Si) && isfinite(amp.Ca) && ...
    amp.Si < 7.5 && amp.Ca > 1.6;

withinCalibrationT = T_degreeC >= 700 & T_degreeC <= 900;
withinCalibrationP_scalar = ...
    isfinite(P_kbar_scalar) && ...
    P_kbar_scalar >= 0 && P_kbar_scalar <= 8;

% Programmatic applicability combines only criteria that can be evaluated
% from the supplied numeric inputs. Quartz equilibrium, magmatic texture,
% alteration, and the full buffering assemblage must be assessed separately
% by the user from petrography and mineral context.
isApplicable = ...
    withinCalibrationT & ...
    repmat(withinCalibrationP_scalar, nT, 1) & ...
    repmat(isIgneousChemicalScreen_scalar, nT, 1);

row = table();

% Temperature is carried through for compatibility with both fixed-T and
% range-T launcher modes.
row.T_degreeC = T_degreeC;
row.T_K = T_degreeC + 273.15;

% Amphibole composition and derived Al-site quantities.
row.Si_amp = repmat(amp.Si, nT, 1);
row.Ti_amp = repmat(amp.Ti, nT, 1);
row.AlT_amp = repmat(AlT_scalar, nT, 1);
row.AlIV_amp = repmat(AlIV_scalar, nT, 1);
row.AlVI_amp = repmat(AlVI_scalar, nT, 1);
row.FeT_amp = repmat(amp.FeT, nT, 1);
row.Fe2_amp = repmat(amp.Fe2, nT, 1);
row.Fe3_amp = repmat(amp.Fe3, nT, 1);
row.Mg_amp = repmat(amp.Mg, nT, 1);
row.Ca_amp = repmat(amp.Ca, nT, 1);
row.Na_amp = repmat(amp.Na, nT, 1);
row.K_amp = repmat(amp.K, nT, 1);
row.Mn_amp = repmat(amp.Mn, nT, 1);
row.Cr_amp = repmat(amp.Cr, nT, 1);

% Pressure and applicability fields.
row.P_kbar = repmat(P_kbar_scalar, nT, 1);
row.P_uncertainty_kbar = repmat(3.0, nT, 1);
row.withinCalibrationT = withinCalibrationT;
row.withinCalibrationP = repmat(withinCalibrationP_scalar, nT, 1);
row.isIgneousChemicalScreen = ...
    repmat(isIgneousChemicalScreen_scalar, nT, 1);
row.isApplicable = isApplicable;

end

function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Extract one-row amphibole cation data. Required variables must exist.
% Missing optional variables are represented as NaN and are never converted
% to zero.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();

amp.Si  = getRequiredValue(data_amphibole, 'Si_cation_apfu', 'Amphibole');
amp.Ti  = getRequiredValue(data_amphibole, 'Ti_cation_apfu', 'Amphibole');
amp.Al  = getRequiredValue(data_amphibole, 'Al_cation_apfu', 'Amphibole');
amp.FeT = getRequiredValue(data_amphibole, 'Fe_cation_apfu', 'Amphibole');
amp.Mg  = getRequiredValue(data_amphibole, 'Mg_cation_apfu', 'Amphibole');
amp.Ca  = getRequiredValue(data_amphibole, 'Ca_cation_apfu', 'Amphibole');
amp.Na  = getRequiredValue(data_amphibole, 'Na_cation_apfu', 'Amphibole');

amp.Mn  = getOptionalValue(data_amphibole, 'Mn_cation_apfu');
amp.K   = getOptionalValue(data_amphibole, 'K_cation_apfu');
amp.Cr  = getOptionalValue(data_amphibole, 'Cr_cation_apfu');
amp.Fe3 = getOptionalValue(data_amphibole, 'Fe3_cation_apfu');

% Fe2+ is retained as NaN when FeT or Fe3 is NaN. No missing value is
% replaced by zero. The consistency check is applied only when both values
% are finite.
if isfinite(amp.Fe3) && isfinite(amp.FeT) && amp.Fe3 > amp.FeT
    error('Amphibole contains Fe3_cation_apfu > Fe_cation_apfu.');
end

amp.Fe2 = amp.FeT - amp.Fe3;

end

function value = getRequiredValue(tbl, variableName, mineralLabel)
% getRequiredValue
% Read one required scalar table variable. NaN is allowed and retained.

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

function value = getOptionalValue(tbl, variableName)
% getOptionalValue
% Read one optional scalar variable. Missing variables and explicit NaN
% values remain NaN; they are not interpreted as zero.

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
