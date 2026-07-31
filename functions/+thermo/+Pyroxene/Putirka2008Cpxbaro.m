function results = Putirka2008Cpxbaro(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/Putirka2008Cpx.m
% Tested with MATLAB R2024b
%
% Clinopyroxene-only thermometer
% Putirka, K.D. (2008), Equation (32d)
% Reviews in Mineralogy and Geochemistry, 69, 61–120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Clinopyroxene analysis and
% calculates temperature using the clinopyroxene-only thermometer of
% Putirka (2008), Equation (32d).
%
% The function accepts a scalar or vector pressure input. It is therefore
% compatible with both startThermoCalc_fixedP and startThermoCalc_rangeP.
% For each selected Clinopyroxene analysis, one output row is returned for
% every pressure value supplied in P_kbar.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Clinopyroxene analysis and stores each
% result block in a preallocated cell buffer. The result blocks are
% concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka (2008) recalibrated the clinopyroxene-only thermometer of Nimis
% and Taylor (2000) using the Nimis and Taylor enstatite-activity model.
% Equation (32d) and the activity definition are given on p. 94.
%
% The clinopyroxene component calculation is based on cations normalized to
% 6 oxygens. Putirka (2008) states that a good clinopyroxene analysis should
% have a cation sum close to 4; the calculation procedure and this quality
% check are described on pp. 88–90 and in Table 3 on p. 89.
%
% Figure 9f on pp. 92–93 reports the following test statistics:
%   Anhydrous experiments : SEE = ±58 degreeC, R^2 = 0.82, n = 910
%   Hydrous experiments   : SEE = ±87 degreeC, R^2 = 0.36, n = 314
%
% The much poorer performance for hydrous experiments indicates that this
% clinopyroxene-only thermometer should be applied cautiously to hydrous
% magmatic systems. A liquid composition is not used by Equation (32d), so
% clinopyroxene–liquid equilibrium cannot be tested within this function.
%
% Putirka (2008) does not provide a separate, strict numerical temperature
% calibration interval for Equation (32d). Figure 9 states that its data are
% the same experimental dataset used in Figure 8, which is restricted to
% pressures below 40 kbar (Figure 8 caption, p. 90). Figure 9f displays the
% model comparison over approximately 600–2400 degreeC.
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) input pressure is greater than 40 kbar, or
%   2) a finite calculated temperature is outside 600–2400 degreeC.
%
% The 600–2400 degreeC interval is an approximate graphical comparison
% envelope from Figure 9f, not a strict calibration boundary stated in the
% text. Results outside this envelope are retained in the output table.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% Supported input styles for Clinopyroxene composition:
%
% (A) Preferred: normalized cation data on a 6-oxygen basis
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   Mn_cation_apfu
%   Ti_cation_apfu
%   Cr_cation_apfu
%   K_cation_apfu          % optional; if absent, assumed zero
%
% (B) Alternative: oxide wt% data
%   SiO2
%   TiO2                   % optional; if absent, assumed zero
%   Al2O3                  % optional; if absent, assumed zero
%   FeO or FeOt
%   MnO                    % optional; if absent, assumed zero
%   MgO
%   CaO
%   Na2O                   % optional; if absent, assumed zero
%   K2O                    % optional; if absent, assumed zero
%   Cr2O3                  % optional; if absent, assumed zero
%
% If the complete required apfu cation-column set exists, those values are
% used directly. Otherwise, cations are calculated from oxide wt% on the
% basis of 6 oxygens. FeO is preferred when both FeO and FeOt columns exist.
%
% NaN values in columns that are used by the thermometer are retained and
% propagated through the calculation; they are never replaced with zero.
% Missing optional columns are assigned zero because no analytical value was
% supplied. Finite negative mineral-composition values are prohibited.
%
% This thermometer does not use a liquid composition. Therefore the policy
% concerning exclusion of F and Cl from cationTotal_liq and from liquid NaN
% warnings is not applicable to this function.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Putirka (2008), Equation (32d), p. 94:
%
%   T(K) = (93100 + 544 * P_kbar) / ...
%          (61.1 + 36.6*XTi_cpx + 10.9*XFe_cpx ...
%           - 0.95*(XAl_cpx + XCr_cpx - XNa_cpx - XK_cpx) ...
%           + 0.395*[ln(aEn_cpx)]^2)
%
% where all X values are cations on a 6-oxygen basis and:
%
%   aEn_cpx = (1 - XCa_cpx - XNa_cpx - XK_cpx) * ...
%             (1 - 0.5*(XAl_cpx + XCr_cpx + XNa_cpx + XK_cpx))
%
% Temperature is reported in Kelvin and degree Celsius.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008Cpx(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing a Clinopyroxene table
%   P_kbar         : pressure in kbar; finite, non-negative numeric scalar
%                    or vector
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Clinopyroxene analysis
%

%% Input validation
% Accept either a single pressure or a pressure vector so that this function
% can be called by both fixed-pressure and pressure-range launchers.
if nargin < 2
    error('Putirka2008Cpx requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve Clinopyroxene dataset
% Extract the required table without modifying its original contents.
disp('=== Step 1: Preparing cpx dataset ===');

if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

if isempty(rawdata_struct.Cpx)
    error('rawdata_struct.Cpx is empty.');
end

dataset_cpx = rawdata_struct.Cpx;

disp('=== Preparing cpx dataset has been finished ===');

%% 2) Initialize output container and applicability limits
% Store each calculation as one table block. Repeated concatenation of the
% complete output table inside the interactive loop is avoided.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Pressure range of the experimental dataset used in Figures 8 and 9.
calibrationP_min_kbar = 0;
calibrationP_max_kbar = 40;

% Approximate graphical comparison envelope shown in Figure 9f. Putirka
% (2008) does not state a strict numerical T calibration interval for Eq.32d.
screeningT_min_degC = 600;
screeningT_max_degC = 2400;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–4) Interactive selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Cpx) ===');

while true
    % ----- Clinopyroxene selection -----
    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Cpx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIdx_cpx);
    disp(['Cpx selected: ' char(string(selectedCode_cpx))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    % Identify NaN values only in variables actually used for the selected
    % input mode. The calculation continues and the names are printed later.
    nanInputNames = findNaNInputs(selectedData_cpx);

    % Stop only for finite negative values or infinite values. NaN values are
    % allowed to propagate through the calculation.
    validateNonNegativeInputs(selectedData_cpx);

    row = calcTemp(selectedData_cpx, P_kbar);

    % Repeat the selected identifier once for every pressure point.
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_cpx'}, 'Before', 1);

    % Store this result as one block. The buffer is enlarged geometrically,
    % not after every iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_cpx)) ': ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_cpx)) ': ' ...
            num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    % Print the pressure-range warning only once because the same pressure
    % vector is used for all selected Clinopyroxene analyses in this call.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental dataset range ' ...
             'used for Putirka (2008) Equation (32d): 0–40 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when a finite temperature lies outside the approximate graphical
    % comparison envelope shown in Figure 9f. The result remains unchanged.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideScreening = finiteTemperature & ...
        (row.T_deg < screeningT_min_degC | ...
         row.T_deg > screeningT_max_degC);

    if any(temperatureOutsideScreening)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximate ' ...
             '600–2400 degreeC comparison envelope shown for Putirka (2008) ' ...
             'Equation (32d), Figure 9f. Putirka (2008) does not state a strict ' ...
             'numerical temperature calibration interval for this equation. ' ...
             '%d of %d finite temperature point(s) are outside the envelope; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s.\n'], ...
            sum(temperatureOutsideScreening), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_cpx)));
    end

    % Print the exact names of input variables that contained NaN. fprintf is
    % used so that the message remains visible even if MATLAB warnings are
    % disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'temperature may be NaN.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Report invalid enstatite activities without stopping the calculation.
    invalidAEn = ~isfinite(row.aEn_cpx) | row.aEn_cpx <= 0;
    if any(invalidAEn)
        fprintf(2, ...
            ['WARNING: aEn_cpx was non-finite or <= 0 for %s at %d of %d ' ...
             'pressure point(s). ln(aEn_cpx) and the corresponding temperature ' ...
             'were retained as NaN.\n'], ...
            char(string(selectedCode_cpx)), ...
            sum(invalidAEn), ...
            numel(row.aEn_cpx));
    end

    % Retain and report all non-finite temperature results.
    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Putirka2008Cpx', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all table blocks once after the interactive loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_cpx)
% findNaNInputs
% Return names of the selected input variables that contain NaN. Missing
% optional variables are not reported because they are intentionally assigned
% zero. F and Cl are irrelevant because this thermometer has no liquid input.

[inputVariableNames, inputMode] = getActiveInputVariables(data_cpx);
nCandidates = numel(inputVariableNames);
nanBuffer = strings(nCandidates, 1);
nNaN = 0;

for i = 1:nCandidates
    variableName = inputVariableNames{i};
    variableValue = data_cpx.(variableName);

    if isnumeric(variableValue) && any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanBuffer(nNaN) = "Cpx." + string(variableName) + ...
            " [" + string(inputMode) + "]";
    end
end

nanInputNames = nanBuffer(1:nNaN);

end

function validateNonNegativeInputs(data_cpx)
% validateNonNegativeInputs
% Stop when an active input contains a finite negative value or Inf. Zero is
% allowed. NaN is intentionally allowed and is propagated by the calculation.

[inputVariableNames, ~] = getActiveInputVariables(data_cpx);
nCandidates = numel(inputVariableNames);
invalidBuffer = strings(nCandidates, 1);
nInvalid = 0;

for i = 1:nCandidates
    variableName = inputVariableNames{i};
    variableValue = data_cpx.(variableName);

    if ~isnumeric(variableValue) || ~isscalar(variableValue)
        error('Cpx variable %s must be a numeric scalar in the selected row.', ...
            variableName);
    end

    if any(isinf(variableValue(:))) || ...
            any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Cpx." + string(variableName);
    end
end

if nInvalid > 0
    invalidInputNames = invalidBuffer(1:nInvalid);
    error(['Putirka2008Cpx: active mineral-composition inputs must not ' ...
           'contain negative finite values or Inf. Invalid value(s) were ' ...
           'found in: ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function [inputVariableNames, inputMode] = getActiveInputVariables(data_cpx)
% getActiveInputVariables
% Identify the exact table columns used for the selected input mode. Existing
% optional variables are included; absent optional variables are omitted.

if hasNormalizedCationColumns(data_cpx)
    inputMode = "apfu";

    requiredNames = { ...
        'Si_cation_apfu', ...
        'Al_cation_apfu', ...
        'Fe_cation_apfu', ...
        'Mg_cation_apfu', ...
        'Ca_cation_apfu', ...
        'Na_cation_apfu', ...
        'Mn_cation_apfu', ...
        'Ti_cation_apfu', ...
        'Cr_cation_apfu'};

    inputVariableNames = requiredNames;
    if ismember('K_cation_apfu', data_cpx.Properties.VariableNames)
        inputVariableNames{end + 1} = 'K_cation_apfu';
    end
else
    inputMode = "oxide";

    requiredNames = {'SiO2', 'MgO', 'CaO'};
    for i = 1:numel(requiredNames)
        if ~ismember(requiredNames{i}, data_cpx.Properties.VariableNames)
            error('Cpx table must contain variable: %s', requiredNames{i});
        end
    end

    if ismember('FeO', data_cpx.Properties.VariableNames)
        ironName = 'FeO';
    elseif ismember('FeOt', data_cpx.Properties.VariableNames)
        ironName = 'FeOt';
    else
        error('Cpx table must contain FeO or FeOt when apfu cations are absent.');
    end

    optionalNames = {'TiO2', 'Al2O3', 'MnO', 'Na2O', 'K2O', 'Cr2O3'};
    inputVariableNames = [requiredNames, {ironName}];

    optionalPresent = ismember(optionalNames, data_cpx.Properties.VariableNames);
    inputVariableNames = [inputVariableNames, optionalNames(optionalPresent)];
end

end

function row = calcTemp(data_cpx, P_kbar)
% calcTemp
% Compute temperatures for one selected Clinopyroxene analysis and a scalar
% or vector of pressures. NaN inputs are retained and propagated.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_kbar ./ 10;

% Prepare one scalar Clinopyroxene composition on a 6-oxygen basis.
cpx = prepareCpxRow(data_cpx);

% Replicate composition values to match the number of pressure points.
Si_cpx = repmat(cpx.Si, nP, 1);
Ti_cpx = repmat(cpx.Ti, nP, 1);
Al_cpx = repmat(cpx.Al, nP, 1);
Fe_cpx = repmat(cpx.Fe, nP, 1);
Mn_cpx = repmat(cpx.Mn, nP, 1);
Mg_cpx = repmat(cpx.Mg, nP, 1);
Ca_cpx = repmat(cpx.Ca, nP, 1);
Na_cpx = repmat(cpx.Na, nP, 1);
K_cpx  = repmat(cpx.K, nP, 1);
Cr_cpx = repmat(cpx.Cr, nP, 1);
cationSum_cpx = repmat(cpx.cationSum, nP, 1);

% Enstatite activity in clinopyroxene after Nimis and Taylor (2000), as
% reproduced by Putirka (2008), Equation (32d), p. 94.
aEn_cpx = (1 - Ca_cpx - Na_cpx - K_cpx) .* ...
          (1 - 0.5 .* (Al_cpx + Cr_cpx + Na_cpx + K_cpx));

% Avoid complex values from log of zero or a negative activity. Such cases
% are represented by NaN and reported by the main function.
ln_aEn_cpx = nan(nP, 1);
validAEn = isfinite(aEn_cpx) & aEn_cpx > 0;
ln_aEn_cpx(validAEn) = log(aEn_cpx(validAEn));

% Putirka (2008), Equation (32d). Note the coefficient 544*P_kbar and the
% squared natural-logarithm term [ln(aEn_cpx)]^2.
denom_T = 61.1 ...
        + 36.6 .* Ti_cpx ...
        + 10.9 .* Fe_cpx ...
        - 0.95 .* (Al_cpx + Cr_cpx - Na_cpx - K_cpx) ...
        + 0.395 .* (ln_aEn_cpx .^ 2);

T_K = (93100 + 544 .* P_kbar) ./ denom_T;
T_deg = T_K - 273.15;

% Pack composition, intermediate variables, and results into a stable output
% table with one row per pressure point.
row.Si_cpx = Si_cpx;
row.Ti_cpx = Ti_cpx;
row.Al_cpx = Al_cpx;
row.Fe_cpx = Fe_cpx;
row.Mn_cpx = Mn_cpx;
row.Mg_cpx = Mg_cpx;
row.Ca_cpx = Ca_cpx;
row.Na_cpx = Na_cpx;
row.K_cpx = K_cpx;
row.Cr_cpx = Cr_cpx;

row.cationSum_cpx = cationSum_cpx;
row.aEn_cpx = aEn_cpx;
row.ln_aEn_cpx = ln_aEn_cpx;
row.ln_aEn_cpx_squared = ln_aEn_cpx .^ 2;
row.denom_T = denom_T;

row.T_K = T_K;
row.T_deg = T_deg;

end

function cpx = prepareCpxRow(data_cpx)
% prepareCpxRow
% Extract one Clinopyroxene row. Prefer apfu cations; otherwise calculate
% cations from oxide wt% on a 6-oxygen basis.

if height(data_cpx) ~= 1
    error('Cpx input must be a 1-row table.');
end

if hasNormalizedCationColumns(data_cpx)
    cpx = getCpxFromNormalizedCations(data_cpx);
else
    cpx = getCpxFromOxides(data_cpx);
end

% All fields must be numeric scalars. NaN is allowed, Inf and finite negative
% values are prohibited. Zero is allowed by policy.
fieldsToCheck = {'Si','Ti','Al','Fe','Mn','Mg','Ca','Na','K','Cr','cationSum'};
for i = 1:numel(fieldsToCheck)
    fieldName = fieldsToCheck{i};
    value = cpx.(fieldName);

    if ~isnumeric(value) || ~isscalar(value)
        error('Cpx contains a non-scalar or non-numeric value for %s.', fieldName);
    end
    if isinf(value) || (isfinite(value) && value < 0)
        error('Cpx contains an invalid negative or infinite value for %s.', fieldName);
    end
end

end

function tf = hasNormalizedCationColumns(tbl)
% hasNormalizedCationColumns
% Return true only when the full required apfu cation-column set is present.

requiredNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

tf = all(ismember(requiredNames, tbl.Properties.VariableNames));

end

function cpx = getCpxFromNormalizedCations(tbl)
% getCpxFromNormalizedCations
% Read apfu cations directly. Existing NaN values are retained. An absent
% optional K_cation_apfu column is assigned zero.

cpx = struct();

cpx.Si = getVarOrError(tbl, 'Si_cation_apfu', 'Cpx');
cpx.Al = getVarOrError(tbl, 'Al_cation_apfu', 'Cpx');
cpx.Fe = getVarOrError(tbl, 'Fe_cation_apfu', 'Cpx');
cpx.Mg = getVarOrError(tbl, 'Mg_cation_apfu', 'Cpx');
cpx.Ca = getVarOrError(tbl, 'Ca_cation_apfu', 'Cpx');
cpx.Na = getVarOrError(tbl, 'Na_cation_apfu', 'Cpx');
cpx.Mn = getVarOrError(tbl, 'Mn_cation_apfu', 'Cpx');
cpx.Ti = getVarOrError(tbl, 'Ti_cation_apfu', 'Cpx');
cpx.Cr = getVarOrError(tbl, 'Cr_cation_apfu', 'Cpx');
cpx.K  = getVarOrZero(tbl, 'K_cation_apfu');

cpx.cationSum = cpx.Si + cpx.Ti + cpx.Al + cpx.Fe + cpx.Mn + ...
                cpx.Mg + cpx.Ca + cpx.Na + cpx.K + cpx.Cr;

end

function cpx = getCpxFromOxides(tbl)
% getCpxFromOxides
% Calculate cations on a 6-oxygen basis from oxide wt%. Existing NaN values
% are retained and propagate through oxygen normalization.

SiO2  = getOxideOrError(tbl, 'SiO2', 'Cpx');
TiO2  = getOxideOrZero(tbl, 'TiO2');
Al2O3 = getOxideOrZero(tbl, 'Al2O3');

if ismember('FeO', tbl.Properties.VariableNames)
    FeO = getOxideOrError(tbl, 'FeO', 'Cpx');
elseif ismember('FeOt', tbl.Properties.VariableNames)
    FeO = getOxideOrError(tbl, 'FeOt', 'Cpx');
else
    error('Cpx table must contain FeO or FeOt when apfu cations are absent.');
end

MnO   = getOxideOrZero(tbl, 'MnO');
MgO   = getOxideOrError(tbl, 'MgO', 'Cpx');
CaO   = getOxideOrError(tbl, 'CaO', 'Cpx');
Na2O  = getOxideOrZero(tbl, 'Na2O');
K2O   = getOxideOrZero(tbl, 'K2O');
Cr2O3 = getOxideOrZero(tbl, 'Cr2O3');

MW = struct();
MW.SiO2  = 60.083;
MW.TiO2  = 79.865;
MW.Al2O3 = 101.961;
MW.FeO   = 71.844;
MW.MnO   = 70.937;
MW.MgO   = 40.304;
MW.CaO   = 56.077;
MW.Na2O  = 61.979;
MW.K2O   = 94.196;
MW.Cr2O3 = 151.990;

mol = struct();
mol.SiO2  = SiO2  ./ MW.SiO2;
mol.TiO2  = TiO2  ./ MW.TiO2;
mol.Al2O3 = Al2O3 ./ MW.Al2O3;
mol.FeO   = FeO   ./ MW.FeO;
mol.MnO   = MnO   ./ MW.MnO;
mol.MgO   = MgO   ./ MW.MgO;
mol.CaO   = CaO   ./ MW.CaO;
mol.Na2O  = Na2O  ./ MW.Na2O;
mol.K2O   = K2O   ./ MW.K2O;
mol.Cr2O3 = Cr2O3 ./ MW.Cr2O3;

oxySum = ...
    2 .* mol.SiO2  + ...
    2 .* mol.TiO2  + ...
    3 .* mol.Al2O3 + ...
    1 .* mol.FeO   + ...
    1 .* mol.MnO   + ...
    1 .* mol.MgO   + ...
    1 .* mol.CaO   + ...
    1 .* mol.Na2O  + ...
    1 .* mol.K2O   + ...
    3 .* mol.Cr2O3;

% If oxySum is NaN, zero, or otherwise invalid, retain the failure as NaN
% instead of stopping or replacing it with a numerical value.
if isnan(oxySum) || ~isfinite(oxySum) || oxySum <= 0
    ORF = NaN;
else
    ORF = 6 ./ oxySum;
end

cpx = struct();
cpx.Si = 1 .* mol.SiO2  .* ORF;
cpx.Ti = 1 .* mol.TiO2  .* ORF;
cpx.Al = 2 .* mol.Al2O3 .* ORF;
cpx.Fe = 1 .* mol.FeO   .* ORF;
cpx.Mn = 1 .* mol.MnO   .* ORF;
cpx.Mg = 1 .* mol.MgO   .* ORF;
cpx.Ca = 1 .* mol.CaO   .* ORF;
cpx.Na = 2 .* mol.Na2O  .* ORF;
cpx.K  = 2 .* mol.K2O   .* ORF;
cpx.Cr = 2 .* mol.Cr2O3 .* ORF;

cpx.cationSum = cpx.Si + cpx.Ti + cpx.Al + cpx.Fe + cpx.Mn + ...
                cpx.Mg + cpx.Ca + cpx.Na + cpx.K + cpx.Cr;

end

function value = getVarOrError(tbl, varName, mineralLabel)
% getVarOrError
% Read a required apfu cation value while retaining NaN.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);
validateScalarCompositionValue(value, mineralLabel, varName);

end

function value = getVarOrZero(tbl, varName)
% getVarOrZero
% Read an optional apfu cation value. An absent column is zero; an existing
% NaN value remains NaN.

if ismember(varName, tbl.Properties.VariableNames)
    value = tbl.(varName);
    validateScalarCompositionValue(value, 'Cpx', varName);
else
    value = 0;
end

end

function value = getOxideOrError(tbl, varName, mineralLabel)
% getOxideOrError
% Read a required oxide value while retaining NaN.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);
validateScalarCompositionValue(value, mineralLabel, varName);

end

function value = getOxideOrZero(tbl, varName)
% getOxideOrZero
% Read an optional oxide value. An absent column is zero; an existing NaN
% value remains NaN.

if ismember(varName, tbl.Properties.VariableNames)
    value = tbl.(varName);
    validateScalarCompositionValue(value, 'Cpx', varName);
else
    value = 0;
end

end

function validateScalarCompositionValue(value, mineralLabel, varName)
% validateScalarCompositionValue
% Require numeric scalar input, allow NaN and zero, and reject finite
% negative values and Inf.

if ~isnumeric(value) || ~isscalar(value)
    error('%s variable %s must be a numeric scalar in the selected row.', ...
        mineralLabel, varName);
end

if isinf(value) || (isfinite(value) && value < 0)
    error('%s contains a negative or infinite value for %s.', ...
        mineralLabel, varName);
end

end
