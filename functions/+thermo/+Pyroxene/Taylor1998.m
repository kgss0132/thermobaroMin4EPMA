function results = Taylor1998(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/Taylor1998.m
% Tested with MATLAB R2024b
%
% Two-pyroxene thermometer
% Taylor, W.R. (1998)
% Neues Jahrbuch fuer Mineralogie - Abhandlungen, 172, 381-408
% DOI: https://doi.org/10.1127/njma/172/1998/381
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis and one
% Clinopyroxene analysis and calculates temperature using the Taylor (1998)
% two-pyroxene thermometer (TA97 formulation).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Opx-Cpx pair, one output row is returned for every pressure
% value supplied in P_kbar. It is therefore compatible with both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Taylor (1998) reports new high-pressure experiments over:
%
%   Temperature : 1050-1260 degreeC
%   Pressure    : 1.0-3.5 GPa (10-35 kbar)
%   Composition : Na2O-, TiO2-, and pyroxene-rich fertile peridotite
%
% These conditions are stated in the abstract on p. 381. The TA97
% two-pyroxene formulation is also presented in the abstract on p. 381.
%
% Taylor (1998) designed TA97 to extend the pressure, temperature, and
% compositional applicability of earlier two-pyroxene thermometers by using
% the new fertile-peridotite experiments together with additional published
% experimental data. However, the paper does not define one simple,
% universal rectangular P-T limit for every composition represented in the
% combined database. This implementation therefore uses the directly stated
% 1.0-3.5 GPa and 1050-1260 degreeC experimental interval as a conservative
% non-stopping screening range rather than as an absolute universal limit.
%
% Important application cautions:
%
%   1) The thermometer was developed for upper-mantle peridotitic and
%      pyroxenitic systems, especially fertile lherzolite and garnet
%      websterite. Application to strongly depleted, highly unusual, or
%      crustal pyroxene compositions is an extrapolation.
%
%   2) The selected Opx and Cpx analyses must represent mutually
%      equilibrated domains. Do not arbitrarily combine porphyroclast cores,
%      rims, neoblasts, exsolution lamellae, or minerals modified during
%      different metamorphic or melt-rock-reaction stages.
%
%   3) Taylor's experiments used graphite inner capsules to minimize Fe
%      loss, longer run durations and fluid or melt to improve equilibration,
%      and low oxygen fugacity so that Fe3+ in silicates was negligible
%      (abstract, p. 381).
%
%   4) Fe_cpx in the equation is treated here as Fe2+. If an
%      Fe3_cation_apfu column is supplied, Fe2 is calculated as total Fe
%      minus Fe3. If the Fe3 column is absent, Fe3 = 0 is assumed. This
%      assumption is most appropriate for reduced mantle pyroxenes similar
%      to those used in the calibration.
%
%   5) Pressure is an explicit term in the thermometer. An independently
%      constrained pressure, a pressure range, or an internally consistent
%      iterative thermobarometric solution should be used.
%
%   6) Pyroxene cations must be normalized consistently, preferably on a
%      6-oxygen basis. The equation requires distinct tetrahedral and
%      octahedral Al terms:
%
%        AlIV = max(0, 2 - Si)
%        AlVI = Al_total - AlIV
%
%      The enstatite activity uses AlVI in its second factor and AlIV in its
%      third factor.
%
%   7) Both calculated enstatite activities must be finite and strictly
%      positive because the thermometer contains ln(aEn_cpx) and
%      ln(aEn_opx).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 1.0-3.5 GPa,
%   2) a finite calculated temperature is outside 1050-1260 degreeC,
%   3) an explicitly stored required or active input is NaN,
%   4) Al-site allocation, enstatite activity, X_ts, or the final
%      thermometer denominator is invalid, or
%   5) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%   rawdata_struct.Cpx : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% Required variables for both Opx and Cpx:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   Ti_cation_apfu
%   Cr_cation_apfu
%
% Optional variables:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%
% Fe_cation_apfu is treated as total Fe. When Fe3_cation_apfu is present:
%
%   Fe2_cation_apfu = Fe_cation_apfu - Fe3_cation_apfu
%
% Existing NaN values are retained and propagated; they are never replaced
% by zero. Zero is assigned only when an optional column itself is absent.
% All finite mineral-composition inputs must be greater than or equal to
% zero. Negative finite values and Inf stop the calculation. A derived
% Fe2 value below zero also stops the calculation.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and liquid NaN warnings is not applicable.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Taylor (1998), TA97 formulation, abstract, p. 381:
%
%   T(K) = [24787 + 678 * P(GPa)] / ...
%          [15.67 + 14.37 * Ti_cpx + 3.69 * Fe_cpx ...
%           - 3.25 * X_ts + (lnKd)^2]
%
% where:
%
%   lnKd = ln(aEn_cpx) - ln(aEn_opx)
%
%   X_ts = (Al + Cr - Na)_cpx
%
%   aEn = (1 - Ca - Na) ...
%       * (1 - AlVI - Cr - Ti) ...
%       * (1 - AlIV/2)^2
%
%   AlIV = max(0, 2 - Si)
%   AlVI = Al_total - AlIV
%
% Pressure is entered as kbar and converted internally:
%
%   P_GPa = P_kbar / 10
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Taylor1998(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Opx and Cpx tables
%   P_kbar         : finite non-negative numeric scalar or vector
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Opx-Cpx pair. T_K, T_degreeC, and T_deg are provided for
%             downstream launcher and plotting compatibility.
%

%% Input validation
if nargin < 2
    error('Taylor1998 requires (rawdata_struct, P_kbar).');
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
disp('=== Step 1: Preparing pyroxene datasets ===');

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end
if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end
if isempty(rawdata_struct.Opx)
    error('rawdata_struct.Opx is empty.');
end
if isempty(rawdata_struct.Cpx)
    error('rawdata_struct.Cpx is empty.');
end

dataset_opx = rawdata_struct.Opx;
dataset_cpx = rawdata_struct.Cpx;

validateRequiredColumns(dataset_opx, 'Opx');
validateRequiredColumns(dataset_cpx, 'Cpx');

disp('=== Preparing pyroxene datasets has been finished ===');

%% 2) Initialize output container and calibration limits
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationP_min_GPa = 1.0;
calibrationP_max_GPa = 3.5;
calibrationT_min_degreeC = 1050;
calibrationT_max_degreeC = 1260;

P_GPa_input = P_kbar ./ 10;
pressureOutsideCalibration = ...
    P_GPa_input < calibrationP_min_GPa | ...
    P_GPa_input > calibrationP_max_GPa;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop and calculation
dataCodes_opx = dataset_opx{:, 1};
displayCodes_opx = cellstr(string(dataCodes_opx));

dataCodes_cpx = dataset_cpx{:, 1};
displayCodes_cpx = cellstr(string(dataCodes_cpx));

disp('=== Step 3: Selecting a data code from the list (Opx) ===');

while true
    % ----- Orthopyroxene selection -----
    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Opx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_opx, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = string(dataCodes_opx(selectedIdx_opx));
    selectedData_opx = dataset_opx(selectedIdx_opx, :);
    disp(['Opx selected: ' char(selectedCode_opx)]);

    % ----- Clinopyroxene selection -----
    disp('=== Step 4: Selecting a data code from the list (Cpx) ===');

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Cpx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_cpx, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = string(dataCodes_cpx(selectedIdx_cpx));
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);
    disp(['Cpx selected: ' char(selectedCode_cpx)]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    nanInputNames = findNaNInputs(selectedData_opx, selectedData_cpx);

    validateNonNegativeInputs(selectedData_opx, selectedData_cpx);

    row = calcTemp(selectedData_opx, selectedData_cpx, P_kbar);

    nRows = height(row);
    row.dataCode_opx = repmat(selectedCode_opx, nRows, 1);
    row.dataCode_cpx = repmat(selectedCode_cpx, nRows, 1);
    row = movevars(row, {'dataCode_opx', 'dataCode_cpx'}, 'Before', 1);

    % Store one completed result block. The complete result table is not
    % enlarged on every interactive iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    printTemperatureSummary( ...
        selectedCode_opx, selectedCode_cpx, row.T_degreeC);

    % Pressure is common to all selected pairs in this function call.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the directly stated ' ...
             'Taylor (1998) experimental range of 1.0-3.5 GPa ' ...
             '(10-35 kbar; abstract, p. 381). %d of %d pressure point(s) ' ...
             'are outside the range; input range = %.4g-%.4g GPa. ' ...
             'The calculation was retained because TA97 also incorporated ' ...
             'additional published experiments.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_GPa_input), ...
            min(P_GPa_input), ...
            max(P_GPa_input));
        pressureWarningIssued = true;
    end

    printTemperatureRangeWarning( ...
        row.T_degreeC, ...
        calibrationT_min_degreeC, ...
        calibrationT_max_degreeC, ...
        selectedCode_opx, selectedCode_cpx);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Taylor (1998) input(s) for ' ...
             '%s & %s: %s.\n' ...
             '         Existing NaN values were retained and propagated. ' ...
             'The temperature may be NaN; variables not entering the final ' ...
             'equation directly may leave the temperature finite.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_cpx), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTerms = findInvalidEquationTerms(row);
    if ~isempty(invalidTerms)
        fprintf(2, ...
            ['WARNING: Invalid Taylor (1998) activity, component, or ' ...
             'equation term(s) were found for %s & %s: %s.\n' ...
             '         Required activity factors, aEn values, X_ts, and the ' ...
             'temperature denominator must be finite and physically valid. ' ...
             'Affected temperatures were retained as NaN.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_cpx), ...
            char(strjoin(invalidTerms, ', ')));
    end

    printNonFiniteTemperatureWarning( ...
        row.T_degreeC, selectedCode_opx, selectedCode_cpx);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Taylor1998', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after the selection loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'primaryTemperatureEquation', 'Taylor (1998) TA97', ...
    'directExperimentalPressure_GPa', [1.0, 3.5], ...
    'directExperimentalTemperature_degreeC', [1050, 1260], ...
    'activityAlTerm', 'AlVI in second factor; AlIV in third factor');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_opx, data_cpx)
% findNaNInputs
% Return names of explicitly stored required or active variables containing
% NaN. Missing optional columns are assigned zero and are not reported.

variables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Fe3_cation_apfu'};

maxNames = 2 .* numel(variables);
nameBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(variables)
    variableName = variables{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        value = data_opx.(variableName);
        if isnumeric(value) && any(isnan(value(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Opx." + string(variableName);
        end
    end

    if ismember(variableName, data_cpx.Properties.VariableNames)
        value = data_cpx.(variableName);
        if isnumeric(value) && any(isnan(value(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Cpx." + string(variableName);
        end
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_opx, data_cpx)
% validateNonNegativeInputs
% Stop when a stored composition value is negative or infinite. Zero is
% allowed. NaN is deliberately allowed and propagated.

variables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Fe3_cation_apfu', ...
    'Mn_cation_apfu'};

maxNames = 2 .* numel(variables);
nameBuffer = strings(maxNames, 1);
nInvalid = 0;

for i = 1:numel(variables)
    variableName = variables{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        value = data_opx.(variableName);
        validateScalarVariable(value, 'Opx', variableName);
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            nameBuffer(nInvalid) = "Opx." + string(variableName);
        end
    end

    if ismember(variableName, data_cpx.Properties.VariableNames)
        value = data_cpx.(variableName);
        validateScalarVariable(value, 'Cpx', variableName);
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            nameBuffer(nInvalid) = "Cpx." + string(variableName);
        end
    end
end

if nInvalid > 0
    invalidNames = nameBuffer(1:nInvalid);
    error(['Taylor1998: pyroxene-composition inputs must not be negative ' ...
           'or infinite. Invalid value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_cpx, P_kbar)
% calcTemp
% Calculate Taylor (1998) TA97 for one selected Opx-Cpx pair and a scalar or
% vector of pressures. NaN inputs and invalid derived activity terms are
% retained as NaN.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

opx = preparePyroxeneRow(data_opx, 'Opx');
cpx = preparePyroxeneRow(data_cpx, 'Cpx');

[AlIV_opx, AlVI_opx] = calcAlSiteTerms(opx.Si, opx.Al);
[AlIV_cpx, AlVI_cpx] = calcAlSiteTerms(cpx.Si, cpx.Al);

% Taylor (1998) enstatite activity factors. The second factor uses AlVI,
% while the third factor uses AlIV.
factorCaNa_opx = 1 - opx.Ca - opx.Na;
factorM1_opx = 1 - AlVI_opx - opx.Cr - opx.Ti;
factorAlIV_opx = 1 - AlIV_opx ./ 2;

factorCaNa_cpx = 1 - cpx.Ca - cpx.Na;
factorM1_cpx = 1 - AlVI_cpx - cpx.Cr - cpx.Ti;
factorAlIV_cpx = 1 - AlIV_cpx ./ 2;

if allPositiveFinite( ...
        factorCaNa_opx, factorM1_opx, factorAlIV_opx)
    aEn_opx = factorCaNa_opx .* factorM1_opx .* ...
        (factorAlIV_opx .^ 2);
else
    aEn_opx = NaN;
end

if allPositiveFinite( ...
        factorCaNa_cpx, factorM1_cpx, factorAlIV_cpx)
    aEn_cpx = factorCaNa_cpx .* factorM1_cpx .* ...
        (factorAlIV_cpx .^ 2);
else
    aEn_cpx = NaN;
end

ln_aEn_opx = safeLogPositive(aEn_opx);
ln_aEn_cpx = safeLogPositive(aEn_cpx);

if isfinite(ln_aEn_opx) && isfinite(ln_aEn_cpx)
    lnKd = ln_aEn_cpx - ln_aEn_opx;
else
    lnKd = NaN;
end

% Tschermak-type Cpx component used in TA97.
X_ts_raw = cpx.Al + cpx.Cr - cpx.Na;
if isfinite(X_ts_raw) && X_ts_raw >= 0
    X_ts = X_ts_raw;
else
    X_ts = NaN;
end

Fe_cpx_term = cpx.Fe2;
Ti_cpx_term = cpx.Ti;

denom_scalar = ...
    15.67 ...
    + 14.37 .* Ti_cpx_term ...
    + 3.69 .* Fe_cpx_term ...
    - 3.25 .* X_ts ...
    + (lnKd .^ 2);

denom_T = repmat(denom_scalar, nP, 1);
numerator_T = 24787 + 678 .* P_GPa;

T_K = nan(nP, 1);
validTemperature = ...
    isfinite(numerator_T) & numerator_T > 0 & ...
    isfinite(denom_T) & denom_T > 0;
T_K(validTemperature) = numerator_T(validTemperature) ./ ...
    denom_T(validTemperature);

T_degreeC = T_K - 273.15;

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PrimaryEquation = repmat("Taylor1998_TA97", nP, 1);

% Orthopyroxene composition.
row.Si_opx = repmat(opx.Si, nP, 1);
row.Al_opx = repmat(opx.Al, nP, 1);
row.FeTot_opx = repmat(opx.FeTot, nP, 1);
row.Fe2_opx = repmat(opx.Fe2, nP, 1);
row.Fe3_opx = repmat(opx.Fe3, nP, 1);
row.Mg_opx = repmat(opx.Mg, nP, 1);
row.Ca_opx = repmat(opx.Ca, nP, 1);
row.Na_opx = repmat(opx.Na, nP, 1);
row.Mn_opx = repmat(opx.Mn, nP, 1);
row.Ti_opx = repmat(opx.Ti, nP, 1);
row.Cr_opx = repmat(opx.Cr, nP, 1);
row.cationSum_opx = repmat(opx.cationSum, nP, 1);

% Clinopyroxene composition.
row.Si_cpx = repmat(cpx.Si, nP, 1);
row.Al_cpx = repmat(cpx.Al, nP, 1);
row.FeTot_cpx = repmat(cpx.FeTot, nP, 1);
row.Fe2_cpx = repmat(cpx.Fe2, nP, 1);
row.Fe3_cpx = repmat(cpx.Fe3, nP, 1);
row.Mg_cpx = repmat(cpx.Mg, nP, 1);
row.Ca_cpx = repmat(cpx.Ca, nP, 1);
row.Na_cpx = repmat(cpx.Na, nP, 1);
row.Mn_cpx = repmat(cpx.Mn, nP, 1);
row.Ti_cpx = repmat(cpx.Ti, nP, 1);
row.Cr_cpx = repmat(cpx.Cr, nP, 1);
row.cationSum_cpx = repmat(cpx.cationSum, nP, 1);

% Al site allocation and enstatite activity terms.
row.AlIV_opx = repmat(AlIV_opx, nP, 1);
row.AlVI_opx = repmat(AlVI_opx, nP, 1);
row.AlIV_cpx = repmat(AlIV_cpx, nP, 1);
row.AlVI_cpx = repmat(AlVI_cpx, nP, 1);

row.factorCaNa_opx = repmat(factorCaNa_opx, nP, 1);
row.factorM1_opx = repmat(factorM1_opx, nP, 1);
row.factorAlIV_opx = repmat(factorAlIV_opx, nP, 1);

row.factorCaNa_cpx = repmat(factorCaNa_cpx, nP, 1);
row.factorM1_cpx = repmat(factorM1_cpx, nP, 1);
row.factorAlIV_cpx = repmat(factorAlIV_cpx, nP, 1);

row.aEn_opx = repmat(aEn_opx, nP, 1);
row.aEn_cpx = repmat(aEn_cpx, nP, 1);
row.ln_aEn_opx = repmat(ln_aEn_opx, nP, 1);
row.ln_aEn_cpx = repmat(ln_aEn_cpx, nP, 1);
row.lnKd = repmat(lnKd, nP, 1);

% TA97 terms and temperature.
row.X_ts_raw = repmat(X_ts_raw, nP, 1);
row.X_ts = repmat(X_ts, nP, 1);
row.Fe_cpx_term = repmat(Fe_cpx_term, nP, 1);
row.Ti_cpx_term = repmat(Ti_cpx_term, nP, 1);
row.numerator_T = numerator_T;
row.denom_T = denom_T;

row.T_K = T_K;
row.T_degreeC = T_degreeC;
row.T_deg = T_degreeC;

end

function px = preparePyroxeneRow(data_px, mineralLabel)
% preparePyroxeneRow
% Extract one pyroxene composition. Existing NaN values remain NaN. Missing
% optional Fe3 and Mn columns are assigned zero.

if height(data_px) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

px = struct();

px.Si = getVarRequired(data_px, 'Si_cation_apfu', mineralLabel);
px.Al = getVarRequired(data_px, 'Al_cation_apfu', mineralLabel);
px.FeTot = getVarRequired(data_px, 'Fe_cation_apfu', mineralLabel);
px.Mg = getVarRequired(data_px, 'Mg_cation_apfu', mineralLabel);
px.Ca = getVarRequired(data_px, 'Ca_cation_apfu', mineralLabel);
px.Na = getVarRequired(data_px, 'Na_cation_apfu', mineralLabel);
px.Ti = getVarRequired(data_px, 'Ti_cation_apfu', mineralLabel);
px.Cr = getVarRequired(data_px, 'Cr_cation_apfu', mineralLabel);

px.Fe3 = getVarOptional( ...
    data_px, 'Fe3_cation_apfu', 0, mineralLabel);
px.Mn = getVarOptional( ...
    data_px, 'Mn_cation_apfu', 0, mineralLabel);

px.Fe2 = px.FeTot - px.Fe3;

if isfinite(px.Fe2) && px.Fe2 < 0
    error(['%s derived Fe2_cation_apfu is negative because ' ...
           'Fe3_cation_apfu exceeds Fe_cation_apfu.'], mineralLabel);
end
if isinf(px.Fe2)
    error('%s derived Fe2_cation_apfu is infinite.', mineralLabel);
end

px.cationSum = ...
    px.Si + px.Al + px.FeTot + px.Mg + px.Ca ...
    + px.Na + px.Mn + px.Ti + px.Cr;

end

function [AlIV, AlVI] = calcAlSiteTerms(Si, AlTotal)
% calcAlSiteTerms
% Calculate tetrahedral and octahedral Al while preserving NaN. A finite
% negative AlVI is treated as an invalid derived site allocation and is
% represented by NaN.

if isfinite(Si)
    AlIV = max(0, 2 - Si);
else
    AlIV = NaN;
end

if isfinite(AlTotal) && isfinite(AlIV)
    AlVI_raw = AlTotal - AlIV;
    if AlVI_raw >= 0
        AlVI = AlVI_raw;
    else
        AlVI = NaN;
    end
else
    AlVI = NaN;
end

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid Al-site, activity-factor, activity, component, and final
% thermometer terms. Composition-derived terms are pressure-independent, so
% only the first output row is needed for those checks.

termBuffer = strings(18, 1);
nTerms = 0;

checks = { ...
    'AlIV_opx', row.AlIV_opx(1), true; ...
    'AlVI_opx', row.AlVI_opx(1), true; ...
    'AlIV_cpx', row.AlIV_cpx(1), true; ...
    'AlVI_cpx', row.AlVI_cpx(1), true; ...
    'factorCaNa_opx', row.factorCaNa_opx(1), false; ...
    'factorM1_opx', row.factorM1_opx(1), false; ...
    'factorAlIV_opx', row.factorAlIV_opx(1), false; ...
    'factorCaNa_cpx', row.factorCaNa_cpx(1), false; ...
    'factorM1_cpx', row.factorM1_cpx(1), false; ...
    'factorAlIV_cpx', row.factorAlIV_cpx(1), false; ...
    'aEn_opx', row.aEn_opx(1), false; ...
    'aEn_cpx', row.aEn_cpx(1), false; ...
    'ln(aEn_opx)', row.ln_aEn_opx(1), true; ...
    'ln(aEn_cpx)', row.ln_aEn_cpx(1), true; ...
    'lnKd', row.lnKd(1), true; ...
    'X_ts', row.X_ts(1), true};

for i = 1:size(checks, 1)
    label = checks{i, 1};
    value = checks{i, 2};
    allowZero = checks{i, 3};

    if allowZero
        invalid = ~isfinite(value) || value < 0;
    else
        invalid = ~isfinite(value) || value <= 0;
    end

    if invalid
        nTerms = nTerms + 1;
        termBuffer(nTerms) = string(label);
    end
end

if any(~isfinite(row.denom_T) | row.denom_T <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "temperature denominator";
end

if any(~isfinite(row.T_K))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "T_K";
end

invalidTerms = termBuffer(1:nTerms);

end

function validateRequiredColumns(tbl, mineralLabel)
% validateRequiredColumns
% Verify the minimum normalized cation columns needed by the thermometer.

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    if ~ismember(variableName, tbl.Properties.VariableNames)
        error('%s table must contain variable: %s', ...
            mineralLabel, variableName);
    end
end

end

function value = getVarRequired(tbl, variableName, mineralLabel)
% getVarRequired
% Read a required numeric scalar while retaining NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);
validateScalarVariable(value, mineralLabel, variableName);

end

function value = getVarOptional( ...
        tbl, variableName, defaultValue, mineralLabel)
% getVarOptional
% Read an optional numeric scalar. An absent column receives defaultValue;
% an explicitly stored NaN remains NaN.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);
    validateScalarVariable(value, mineralLabel, variableName);
else
    value = defaultValue;
end

end

function validateScalarVariable(value, mineralLabel, variableName)
% validateScalarVariable
% Require one numeric scalar. NaN is allowed; negative finite values and Inf
% are handled separately.

if ~isnumeric(value) || ~isscalar(value)
    error('%s variable %s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end

function result = allPositiveFinite(varargin)
% allPositiveFinite
% Return true only when every supplied scalar is finite and greater than zero.

result = true;
for i = 1:nargin
    value = varargin{i};
    if ~isfinite(value) || value <= 0
        result = false;
        return
    end
end

end

function value = safeLogPositive(value)
% safeLogPositive
% Return ln(value) only for a finite, strictly positive argument.

if isfinite(value) && value > 0
    value = log(value);
else
    value = NaN;
end

end

function printTemperatureSummary( ...
        selectedCode_opx, selectedCode_cpx, temperatureValues)
% printTemperatureSummary
% Display one temperature or the first-to-last values for a pressure vector.

if isscalar(temperatureValues)
    disp([char(selectedCode_opx) ' & ' char(selectedCode_cpx) ': ' ...
        num2str(temperatureValues) ' degreeC']);
else
    disp([char(selectedCode_opx) ' & ' char(selectedCode_cpx) ': ' ...
        num2str(temperatureValues(1)) ' to ' ...
        num2str(temperatureValues(end)) ' degreeC']);
end

end

function printTemperatureRangeWarning( ...
        temperatureValues, minimumTemperature, maximumTemperature, ...
        selectedCode_opx, selectedCode_cpx)
% printTemperatureRangeWarning
% Warn when finite temperatures lie outside the directly stated Taylor
% experimental interval. Results are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);
    fprintf(2, ...
        ['WARNING: Calculated Taylor (1998) temperature is outside the ' ...
         'directly stated experimental range of %.4g-%.4g degreeC ' ...
         '(abstract, p. 381). %d of %d finite point(s) are outside the ' ...
         'range; calculated finite range = %.4g-%.4g degreeC for %s & %s. ' ...
         'The result is retained because TA97 also incorporated additional ' ...
         'published experiments.\n'], ...
        minimumTemperature, ...
        maximumTemperature, ...
        sum(outsideMask), ...
        sum(finiteMask), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(selectedCode_opx), ...
        char(selectedCode_cpx));
end

end

function printNonFiniteTemperatureWarning( ...
        temperatureValues, selectedCode_opx, selectedCode_cpx)
% printNonFiniteTemperatureWarning
% Report NaN and Inf temperatures without modifying or deleting them.

invalidMask = ~isfinite(temperatureValues);
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Taylor (1998) temperature values were ' ...
         'calculated for %s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the ' ...
         'calculation has not been stopped.\n'], ...
        char(selectedCode_opx), ...
        char(selectedCode_cpx), ...
        sum(invalidMask), ...
        numel(temperatureValues), ...
        sum(isnan(temperatureValues)), ...
        sum(isinf(temperatureValues)));
end

end
