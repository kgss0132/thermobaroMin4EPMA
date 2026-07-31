function results = Wells1977(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/Wells1977.m
% Intended for MATLAB R2024b
%
% Semi-empirical two-pyroxene thermometer
% Wells, P.R.A. (1977)
% Contributions to Mineralogy and Petrology, 62, 129-139
% DOI: https://doi.org/10.1007/BF00372872
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis and one
% Clinopyroxene analysis and calculates temperature using Equation (5) of
% Wells (1977).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Opx-Cpx pair, one output row is returned for every pressure
% supplied in P_kbar. Equation (5) contains no pressure term, so all pressure
% rows for a given pair contain the same calculated temperature. This
% interface is compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Wells (1977) calibrated the Fe dependence of the multicomponent
% two-pyroxene solvus using 43 experimentally equilibrated Opx-Cpx
% assemblages. The empirical relation and thermometer are given as
% Equations (4) and (5) on p. 135.
%
% DIRECT MULTICOMPONENT CALIBRATION RANGE (p. 135):
%
%   Temperature : 875-1500 degreeC
%   XFe_opx     : 0.0-0.85
%
% where:
%
%   XFe_opx = Fe2_opx / (Fe2_opx + Mg_opx)
%
% BROADER EXPERIMENTAL CONSISTENCY RANGE (pp. 135 and 137-138):
%
%   Temperature : approximately 785-1500 degreeC
%   XFe_opx     : approximately 0.0-1.0
%   Al2O3_cpx   : approximately 0-10 wt%
%
% The wider interval includes Fe-rich hedenbergite-ferrosilite experiments
% that were not all used directly in the regression. This implementation
% therefore uses the more conservative direct multicomponent range of
% 875-1500 degreeC and XFe_opx = 0.0-0.85 for non-stopping screening
% warnings.
%
% CALIBRATION PERFORMANCE (p. 135; Conclusions, pp. 137-138):
%
%   Calibration data:
%     mean absolute deviation approximately 40 degreeC
%     approximately 90% of data within 70 degreeC
%
%   All experimental data considered:
%     mean absolute deviation approximately 64 degreeC
%     approximately 70% of data within 70 degreeC
%
% Wells (1977) concludes that Equation (5) should generally provide
% temperatures accurate to about 70 degreeC only when both composition and
% equilibration conditions fall within the reported calibration ranges.
%
% PRESSURE (pp. 131-132):
%
% Equation (5) contains no pressure term. For the simple diopside-enstatite
% solvus, the experimental database shown by Wells spans approximately
% 1 bar to 40 kbar. Although the diopside limb is moderately pressure
% dependent at high temperature, Wells considered the pressure effect to be
% masked by experimental uncertainty and negligible over that interval.
%
% The 1-bar to 40-kbar interval is not presented as a strict rectangular
% pressure calibration range for every multicomponent composition used in
% Equation (5). It is used here only as a non-stopping pressure-screening
% envelope. P_kbar is retained in the output but is not used in the
% temperature equation.
%
% IMPORTANT APPLICATION CAUTIONS
%
%   1) Opx and Cpx must represent an equilibrated mineral pair. The accuracy
%      statement is conditional on both composition and equilibration
%      conditions being within the calibration domain (pp. 137-138).
%
%   2) Some experiments at temperatures below 1100 degreeC and at Mg-rich
%      compositions differ from Equation (5) by more than 100 degreeC.
%      Wells attributes these discrepancies to possible disequilibrium and
%      analytical difficulties with fine-grained run products (p. 135).
%      Low-temperature Mg-rich natural pairs therefore require especially
%      careful textural and chemical equilibrium screening.
%
%   3) Do not arbitrarily combine Opx and Cpx cores, rims, porphyroclasts,
%      neoblasts, exsolution lamellae, or grains modified during different
%      metamorphic, cooling, or melt-rock-reaction stages.
%
%   4) Na in Cpx affects the adopted activity-composition calculation.
%      Wells shows that arbitrarily assigned Na2O in an incomplete Cpx
%      analysis can produce large temperature errors (p. 135). Measured Na
%      should be used rather than an assumed value.
%
%   5) The calibration includes aluminous Cpx. Wells reports reasonable
%      agreement for experimental Cpx containing substantial Al2O3, and the
%      final stated multicomponent range extends to approximately 10 wt%
%      Al2O3 in Cpx (pp. 133 and 137-138).
%
%   6) The ideal two-site solution model is acknowledged to be physically
%      simplified. Equation (5) is a semi-empirical interpolation of the
%      available phase-equilibrium data, not a complete non-ideal
%      thermodynamic model (pp. 130-131 and 137).
%
%   7) The original Wells calculation adopts the activity-composition
%      relations of Wood and Banno (1973). The supplied thermoCalcMin source
%      code instead estimates the Mg2Si2O6 activity using a simplified
%      Wood-Banno-style site allocation:
%
%        aEn = XMg_M1 * XMg_M2
%
%      This approximation is retained here for continuity with the supplied
%      implementation. Results may therefore differ from implementations
%      using the complete original activity allocation.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside the approximate 1-bar to 40-kbar
%      experimental pressure envelope,
%   2) a finite calculated temperature is outside 875-1500 degreeC,
%   3) XFe_opx is outside the direct calibration range 0.0-0.85,
%   4) an explicitly stored calculation input is NaN,
%   5) a site allocation, activity, logarithm, or Equation (5) denominator
%      is invalid, or
%   6) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%   rawdata_struct.Cpx : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialogs. Pyroxene cations should be normalized
% consistently, preferably on a 6-oxygen basis.
%
% Required variables for both Opx and Cpx:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   Mn_cation_apfu
%   Ti_cation_apfu
%   Cr_cation_apfu
%
% Optional variable:
%   Fe3_cation_apfu
%
% Na, Mn, Ti, and Cr are required here because they enter the simplified
% site-allocation procedure used to estimate the enstatite activities.
% Treating an absent analytical column as a true zero would otherwise mix
% "not measured" with "measured as zero".
%
% Fe_cation_apfu is treated as total Fe. When Fe3_cation_apfu is present:
%
%   Fe2_cation_apfu = Fe_cation_apfu - Fe3_cation_apfu
%
% If the Fe3 column is absent, Fe3 is assumed to be zero. An explicitly
% stored NaN is retained and propagated; it is never replaced by zero.
%
% All finite mineral-composition inputs must be greater than or equal to
% zero. Negative finite values and Inf stop the calculation. A derived Fe2
% value below zero also stops the calculation.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and liquid NaN warnings is not applicable.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Wells (1977), Equation (5), p. 135:
%
%   T(K) = 7341 / [3.355 + 2.44*XFe_opx - lnK]
%
% where:
%
%   XFe_opx = Fe2_opx / (Fe2_opx + Mg_opx)
%
%   lnK = ln(aEn_cpx / aEn_opx)
%
% The supplied implementation estimates:
%
%   aEn = XMg_M1 * XMg_M2
%
% using the following simplified site allocation:
%
%   AlIV = max(0, 2 - Si)
%   AlVI = Al_total - AlIV
%
%   M1 fixed occupants:
%     AlVI + Cr + Ti + Fe3
%
%   M2 fixed occupants:
%     Ca + Na + Mn
%
%   Remaining M1 and M2 capacity is filled by Mg and Fe2 in proportion to
%   Mg/(Mg + Fe2) and Fe2/(Mg + Fe2).
%
% Temperature is calculated in Kelvin and returned in both Kelvin and
% degreeC:
%
%   T_degreeC = T_K - 273.15
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Wells1977(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Opx and Cpx tables
%   P_kbar         : finite non-negative numeric scalar or vector; retained
%                    in output but not used in Equation (5)
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Opx-Cpx pair. T_K, T_degreeC, and T_deg are supplied for
%             downstream launcher and plotting compatibility.
%

%% Input validation
if nargin < 2
    error('Wells1977 requires (rawdata_struct, P_kbar).');
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

%% 2) Initialize output container and screening limits
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct multicomponent calibration range reported on p. 135.
calibrationT_min_degreeC = 875;
calibrationT_max_degreeC = 1500;
calibrationXFeOpx_min = 0.0;
calibrationXFeOpx_max = 0.85;

% Approximate simple-system pressure envelope discussed on pp. 131-132.
% 1 bar = 0.001 kbar.
screeningP_min_kbar = 0.001;
screeningP_max_kbar = 40;

pressureOutsideScreening = ...
    P_kbar < screeningP_min_kbar | P_kbar > screeningP_max_kbar;
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

    % Store one completed result block. The complete output table is not
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
    if any(pressureOutsideScreening) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate 1-bar to ' ...
             '40-kbar simple-system experimental envelope discussed by ' ...
             'Wells (1977, pp. 131-132). %d of %d pressure point(s) are ' ...
             'outside the envelope; input range = %.4g-%.4g kbar. ' ...
             'Equation (5) has no pressure term, and this envelope is not a ' ...
             'strict multicomponent pressure calibration range.\n'], ...
            sum(pressureOutsideScreening), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    printTemperatureRangeWarning( ...
        row.T_degreeC, ...
        calibrationT_min_degreeC, ...
        calibrationT_max_degreeC, ...
        selectedCode_opx, selectedCode_cpx);

    % Direct multicomponent calibration range for Opx XFe is 0.0-0.85.
    if isfinite(row.XFe_opx(1)) && ...
            (row.XFe_opx(1) < calibrationXFeOpx_min || ...
             row.XFe_opx(1) > calibrationXFeOpx_max)
        fprintf(2, ...
            ['WARNING: XFe_opx = %.4g for %s & %s is outside the direct ' ...
             'multicomponent calibration range 0.0-0.85 reported by Wells ' ...
             '(1977, p. 135). Equation (5) is described as broadly ' ...
             'consistent with additional experimental data up to ' ...
             'XFe_opx = 1.0, but the result is an extrapolation beyond the ' ...
             'direct regression range and has been retained.\n'], ...
            row.XFe_opx(1), ...
            char(selectedCode_opx), ...
            char(selectedCode_cpx));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Wells (1977) thermometer ' ...
             'input(s) for %s & %s: %s.\n' ...
             '         Existing NaN values were retained and propagated; ' ...
             'the calculated temperature may be NaN.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_cpx), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTerms = findInvalidEquationTerms(row);
    if ~isempty(invalidTerms)
        fprintf(2, ...
            ['WARNING: Invalid Wells (1977) site-allocation, activity, or ' ...
             'equation term(s) were found for %s & %s: %s.\n' ...
             '         Required site fractions, enstatite activities, ' ...
             'logarithm arguments, and the Equation (5) denominator must ' ...
             'be finite and physically valid. Affected temperatures were ' ...
             'retained as NaN.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_cpx), ...
            char(strjoin(invalidTerms, ', ')));
    end

    printNonFiniteTemperatureWarning( ...
        row.T_degreeC, selectedCode_opx, selectedCode_cpx);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Wells1977', ...
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
    'primaryTemperatureEquation', 'Wells (1977) Eq. (5)', ...
    'directCalibrationTemperature_degreeC', [875, 1500], ...
    'directCalibrationXFeOpx', [0.0, 0.85], ...
    'broaderExperimentalTemperature_degreeC', [785, 1500], ...
    'broaderExperimentalXFeOpx', [0.0, 1.0], ...
    'pressureUsedInEquation', false, ...
    'pressureScreening_kbar', [0.001, 40], ...
    'reportedTypicalAccuracy_degreeC', 70, ...
    'activityImplementation', ...
        'Simplified Wood-Banno-style XMg_M1*XMg_M2 allocation');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_opx, data_cpx)
% findNaNInputs
% Return names of explicitly stored variables used by the implemented site
% allocation or thermometer that contain NaN.

activeVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Fe3_cation_apfu'};

maxNames = 2 .* numel(activeVariables);
nameBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(activeVariables)
    variableName = activeVariables{i};

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
% Stop when a stored mineral-composition input is negative or infinite.
% Zero is allowed. NaN is deliberately allowed and propagated.

variablesToCheck = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Fe3_cation_apfu'};

maxNames = 2 .* numel(variablesToCheck);
nameBuffer = strings(maxNames, 1);
nInvalid = 0;

for i = 1:numel(variablesToCheck)
    variableName = variablesToCheck{i};

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
    error(['Wells1977: pyroxene-composition inputs must not be negative ' ...
           'or infinite. Invalid value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_cpx, P_kbar)
% calcTemp
% Calculate Wells (1977) Equation (5) for one selected Opx-Cpx pair and
% repeat the pressure-independent result for every supplied pressure. NaN
% inputs and invalid derived terms are retained as NaN.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

opx = preparePyroxeneRow(data_opx, 'Opx');
cpx = preparePyroxeneRow(data_cpx, 'Cpx');

site_opx = calcSiteFractions(opx);
site_cpx = calcSiteFractions(cpx);

% Simplified Wood-Banno-style enstatite activity retained from the supplied
% implementation. Both Mg site fractions must be finite and strictly
% positive for a logarithm to be defined.
if isPositiveFinite(site_opx.XMg_M1) && ...
        isPositiveFinite(site_opx.XMg_M2)
    aEn_opx = site_opx.XMg_M1 .* site_opx.XMg_M2;
else
    aEn_opx = NaN;
end

if isPositiveFinite(site_cpx.XMg_M1) && ...
        isPositiveFinite(site_cpx.XMg_M2)
    aEn_cpx = site_cpx.XMg_M1 .* site_cpx.XMg_M2;
else
    aEn_cpx = NaN;
end

ln_aEn_opx = safeLogPositive(aEn_opx);
ln_aEn_cpx = safeLogPositive(aEn_cpx);

if isfinite(ln_aEn_opx) && isfinite(ln_aEn_cpx)
    lnK = ln_aEn_cpx - ln_aEn_opx;
else
    lnK = NaN;
end

% Orthopyroxene Fe number used by Equation (5).
denom_XFe_opx = opx.Fe2 + opx.Mg;
if isfinite(denom_XFe_opx) && denom_XFe_opx > 0 && ...
        isfinite(opx.Fe2)
    XFe_opx = opx.Fe2 ./ denom_XFe_opx;
else
    XFe_opx = NaN;
end

% Wells (1977), Equation (5), p. 135.
denom_scalar = 3.355 + 2.44 .* XFe_opx - lnK;

if isfinite(denom_scalar) && denom_scalar > 0
    T_scalar_K = 7341 ./ denom_scalar;
else
    T_scalar_K = NaN;
end

if ~isfinite(T_scalar_K) || T_scalar_K <= 0
    T_scalar_K = NaN;
end

T_K = repmat(T_scalar_K, nP, 1);
T_degreeC = T_K - 273.15;

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PressureUsedInEquation = false(nP, 1);
row.PrimaryEquation = repmat("Wells1977_Eq5", nP, 1);

% Orthopyroxene composition.
row.Si_opx = repmat(opx.Si, nP, 1);
row.Al_opx = repmat(opx.Al, nP, 1);
row.Fe_total_opx = repmat(opx.Fe_total, nP, 1);
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
row.Fe_total_cpx = repmat(cpx.Fe_total, nP, 1);
row.Fe2_cpx = repmat(cpx.Fe2, nP, 1);
row.Fe3_cpx = repmat(cpx.Fe3, nP, 1);
row.Mg_cpx = repmat(cpx.Mg, nP, 1);
row.Ca_cpx = repmat(cpx.Ca, nP, 1);
row.Na_cpx = repmat(cpx.Na, nP, 1);
row.Mn_cpx = repmat(cpx.Mn, nP, 1);
row.Ti_cpx = repmat(cpx.Ti, nP, 1);
row.Cr_cpx = repmat(cpx.Cr, nP, 1);
row.cationSum_cpx = repmat(cpx.cationSum, nP, 1);

% Orthopyroxene site allocation.
row.AlIV_opx = repmat(site_opx.AlIV, nP, 1);
row.AlVI_opx = repmat(site_opx.AlVI, nP, 1);
row.M1_fixed_opx = repmat(site_opx.M1_fixed, nP, 1);
row.M2_fixed_opx = repmat(site_opx.M2_fixed, nP, 1);
row.M1_remaining_opx = repmat(site_opx.M1_remaining, nP, 1);
row.M2_remaining_opx = repmat(site_opx.M2_remaining, nP, 1);
row.XMg_M1_opx = repmat(site_opx.XMg_M1, nP, 1);
row.XMg_M2_opx = repmat(site_opx.XMg_M2, nP, 1);
row.M1_total_opx = repmat(site_opx.M1_total, nP, 1);
row.M2_total_opx = repmat(site_opx.M2_total, nP, 1);

% Clinopyroxene site allocation.
row.AlIV_cpx = repmat(site_cpx.AlIV, nP, 1);
row.AlVI_cpx = repmat(site_cpx.AlVI, nP, 1);
row.M1_fixed_cpx = repmat(site_cpx.M1_fixed, nP, 1);
row.M2_fixed_cpx = repmat(site_cpx.M2_fixed, nP, 1);
row.M1_remaining_cpx = repmat(site_cpx.M1_remaining, nP, 1);
row.M2_remaining_cpx = repmat(site_cpx.M2_remaining, nP, 1);
row.XMg_M1_cpx = repmat(site_cpx.XMg_M1, nP, 1);
row.XMg_M2_cpx = repmat(site_cpx.XMg_M2, nP, 1);
row.M1_total_cpx = repmat(site_cpx.M1_total, nP, 1);
row.M2_total_cpx = repmat(site_cpx.M2_total, nP, 1);

% Activity, composition, and equation terms.
row.aEn_opx = repmat(aEn_opx, nP, 1);
row.aEn_cpx = repmat(aEn_cpx, nP, 1);
row.ln_aEn_opx = repmat(ln_aEn_opx, nP, 1);
row.ln_aEn_cpx = repmat(ln_aEn_cpx, nP, 1);
row.lnK = repmat(lnK, nP, 1);
row.denom_XFe_opx = repmat(denom_XFe_opx, nP, 1);
row.XFe_opx = repmat(XFe_opx, nP, 1);
row.denom_T = repmat(denom_scalar, nP, 1);

row.T_K = T_K;
row.T_degreeC = T_degreeC;
row.T_deg = T_degreeC;

end

function px = preparePyroxeneRow(data_px, mineralLabel)
% preparePyroxeneRow
% Extract one pyroxene composition. Existing NaN values remain NaN. A
% missing Fe3 column is assigned zero.

if height(data_px) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

px = struct();

px.Si = getVarRequired(data_px, 'Si_cation_apfu', mineralLabel);
px.Al = getVarRequired(data_px, 'Al_cation_apfu', mineralLabel);
px.Fe_total = getVarRequired( ...
    data_px, 'Fe_cation_apfu', mineralLabel);
px.Mg = getVarRequired(data_px, 'Mg_cation_apfu', mineralLabel);
px.Ca = getVarRequired(data_px, 'Ca_cation_apfu', mineralLabel);
px.Na = getVarRequired(data_px, 'Na_cation_apfu', mineralLabel);
px.Mn = getVarRequired(data_px, 'Mn_cation_apfu', mineralLabel);
px.Ti = getVarRequired(data_px, 'Ti_cation_apfu', mineralLabel);
px.Cr = getVarRequired(data_px, 'Cr_cation_apfu', mineralLabel);

px.Fe3 = getVarOptional( ...
    data_px, 'Fe3_cation_apfu', 0, mineralLabel);

px.Fe2 = px.Fe_total - px.Fe3;

if isfinite(px.Fe2) && px.Fe2 < 0
    error(['%s derived Fe2_cation_apfu is negative because ' ...
           'Fe3_cation_apfu exceeds Fe_cation_apfu.'], mineralLabel);
end
if isinf(px.Fe2)
    error('%s derived Fe2_cation_apfu is infinite.', mineralLabel);
end

px.cationSum = ...
    px.Si + px.Al + px.Fe_total + px.Mg + px.Ca ...
    + px.Na + px.Mn + px.Ti + px.Cr;

end

function site = calcSiteFractions(px)
% calcSiteFractions
% Perform the simplified Wood-Banno-style site allocation retained from the
% supplied implementation. Invalid derived terms become NaN so that the
% result can be retained and reported rather than stopping the entire run.

site = struct();

if isfinite(px.Si)
    AlIV = max(0, 2 - px.Si);
else
    AlIV = NaN;
end

if isfinite(px.Al) && isfinite(AlIV)
    AlVI_raw = px.Al - AlIV;
    if AlVI_raw >= 0
        AlVI = AlVI_raw;
    else
        AlVI = NaN;
    end
else
    AlVI = NaN;
end

if all(isfinite([AlVI, px.Cr, px.Ti, px.Fe3]))
    M1_fixed_raw = AlVI + px.Cr + px.Ti + px.Fe3;
    if M1_fixed_raw >= 0 && M1_fixed_raw <= 1
        M1_fixed = M1_fixed_raw;
        M1_remaining = 1 - M1_fixed;
    else
        M1_fixed = NaN;
        M1_remaining = NaN;
    end
else
    M1_fixed = NaN;
    M1_remaining = NaN;
end

if all(isfinite([px.Ca, px.Na, px.Mn]))
    M2_fixed_raw = px.Ca + px.Na + px.Mn;
    if M2_fixed_raw >= 0 && M2_fixed_raw <= 1
        M2_fixed = M2_fixed_raw;
        M2_remaining = 1 - M2_fixed;
    else
        M2_fixed = NaN;
        M2_remaining = NaN;
    end
else
    M2_fixed = NaN;
    M2_remaining = NaN;
end

MgFe_total = px.Mg + px.Fe2;
if isfinite(MgFe_total) && MgFe_total > 0 && ...
        isfinite(px.Mg) && isfinite(px.Fe2) && ...
        isfinite(M1_remaining) && isfinite(M2_remaining)

    Mg_fraction = px.Mg ./ MgFe_total;
    Fe_fraction = px.Fe2 ./ MgFe_total;

    Mg_M1 = M1_remaining .* Mg_fraction;
    Fe2_M1 = M1_remaining .* Fe_fraction;

    Mg_M2 = M2_remaining .* Mg_fraction;
    Fe2_M2 = M2_remaining .* Fe_fraction;

    M1_total = M1_fixed + Mg_M1 + Fe2_M1;
    M2_total = M2_fixed + Mg_M2 + Fe2_M2;

    if isfinite(M1_total) && M1_total > 0
        XMg_M1 = Mg_M1 ./ M1_total;
    else
        XMg_M1 = NaN;
        M1_total = NaN;
    end

    if isfinite(M2_total) && M2_total > 0
        XMg_M2 = Mg_M2 ./ M2_total;
    else
        XMg_M2 = NaN;
        M2_total = NaN;
    end
else
    XMg_M1 = NaN;
    XMg_M2 = NaN;
    M1_total = NaN;
    M2_total = NaN;
end

site.AlIV = AlIV;
site.AlVI = AlVI;
site.M1_fixed = M1_fixed;
site.M2_fixed = M2_fixed;
site.M1_remaining = M1_remaining;
site.M2_remaining = M2_remaining;
site.XMg_M1 = XMg_M1;
site.XMg_M2 = XMg_M2;
site.M1_total = M1_total;
site.M2_total = M2_total;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid site-allocation, activity, logarithm, composition, and
% Equation (5) terms. Logarithms and lnK are allowed to be negative; only
% their finiteness is tested.

termBuffer = strings(24, 1);
nTerms = 0;

positiveChecks = { ...
    'XMg_M1_opx', row.XMg_M1_opx(1); ...
    'XMg_M2_opx', row.XMg_M2_opx(1); ...
    'XMg_M1_cpx', row.XMg_M1_cpx(1); ...
    'XMg_M2_cpx', row.XMg_M2_cpx(1); ...
    'aEn_opx', row.aEn_opx(1); ...
    'aEn_cpx', row.aEn_cpx(1); ...
    'denom_XFe_opx', row.denom_XFe_opx(1)};

for i = 1:size(positiveChecks, 1)
    label = positiveChecks{i, 1};
    value = positiveChecks{i, 2};

    if ~isfinite(value) || value <= 0
        nTerms = nTerms + 1;
        termBuffer(nTerms) = string(label);
    end
end

boundedChecks = { ...
    'AlIV_opx', row.AlIV_opx(1), 0, 2; ...
    'AlVI_opx', row.AlVI_opx(1), 0, 1; ...
    'M1_fixed_opx', row.M1_fixed_opx(1), 0, 1; ...
    'M2_fixed_opx', row.M2_fixed_opx(1), 0, 1; ...
    'AlIV_cpx', row.AlIV_cpx(1), 0, 2; ...
    'AlVI_cpx', row.AlVI_cpx(1), 0, 1; ...
    'M1_fixed_cpx', row.M1_fixed_cpx(1), 0, 1; ...
    'M2_fixed_cpx', row.M2_fixed_cpx(1), 0, 1; ...
    'XFe_opx', row.XFe_opx(1), 0, 1};

for i = 1:size(boundedChecks, 1)
    label = boundedChecks{i, 1};
    value = boundedChecks{i, 2};
    minimumValue = boundedChecks{i, 3};
    maximumValue = boundedChecks{i, 4};

    if ~isfinite(value) || value < minimumValue || value > maximumValue
        nTerms = nTerms + 1;
        termBuffer(nTerms) = string(label);
    end
end

finiteOnlyChecks = { ...
    'ln(aEn_opx)', row.ln_aEn_opx(1); ...
    'ln(aEn_cpx)', row.ln_aEn_cpx(1); ...
    'lnK', row.lnK(1)};

for i = 1:size(finiteOnlyChecks, 1)
    label = finiteOnlyChecks{i, 1};
    value = finiteOnlyChecks{i, 2};

    if ~isfinite(value)
        nTerms = nTerms + 1;
        termBuffer(nTerms) = string(label);
    end
end

if any(~isfinite(row.denom_T) | row.denom_T <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Equation (5) denominator";
end

if any(~isfinite(row.T_K) | row.T_K <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "T_K";
end

invalidTerms = termBuffer(1:nTerms);

end

function validateRequiredColumns(tbl, mineralLabel)
% validateRequiredColumns
% Verify all normalized cation columns used by the simplified activity
% calculation.

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Mn_cation_apfu', ...
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

function result = isPositiveFinite(value)
% isPositiveFinite
% Return true only for a finite scalar greater than zero.

result = isfinite(value) && value > 0;

end

function value = safeLogPositive(value)
% safeLogPositive
% Return the natural logarithm only for a finite, strictly positive value.

if isPositiveFinite(value)
    value = log(value);
else
    value = NaN;
end

end

function printTemperatureSummary( ...
        selectedCode_opx, selectedCode_cpx, temperatureValues)
% printTemperatureSummary
% Display one temperature or the first-to-last values for a pressure vector.

if numel(temperatureValues) == 1
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
% Warn when finite temperatures lie outside the direct multicomponent
% calibration range. Results are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);
    fprintf(2, ...
        ['WARNING: Calculated Wells (1977) temperature is outside the ' ...
         'direct multicomponent calibration range of %.4g-%.4g degreeC ' ...
         '(p. 135). %d of %d finite point(s) are outside the range; ' ...
         'calculated finite range = %.4g-%.4g degreeC for %s & %s. ' ...
         'Wells reports broader consistency with some experimental data ' ...
         'down to approximately 785 degreeC, but the result is outside the ' ...
         'direct regression interval and has been retained.\n'], ...
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
        ['WARNING: Non-finite Wells (1977) temperature values were ' ...
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
