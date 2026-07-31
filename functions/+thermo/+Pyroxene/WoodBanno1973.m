function results = WoodBanno1973(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/WoodBanno1973.m
% Tested with MATLAB R2024b
%
% Semi-empirical Orthopyroxene-Clinopyroxene thermometer
% Wood, B.J. and Banno, S. (1973)
% Contributions to Mineralogy and Petrology, 42, 109-124
% DOI: https://doi.org/10.1007/BF00371501
%
% -------------------------------------------------------------------------
% OVERVIEW
%
% This function interactively pairs one Orthopyroxene analysis and one
% Clinopyroxene analysis and calculates temperature using Equation (27) of
% Wood and Banno (1973).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Opx-Cpx pair, one output row is returned for every pressure
% supplied in P_kbar. Equation (27) contains no explicit pressure term, so
% all pressure rows for one selected pair contain the same calculated
% temperature. This interface is compatible with both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Wood and Banno (1973) developed a semi-empirical two-pyroxene thermometer
% by combining the diopside-enstatite miscibility relation with empirical
% corrections for the effect of Fe2 on the Opx-Cpx miscibility gap.
%
% The idealized-activity model and the pressure discussion are presented on
% p. 117. Equation (27) and the explicit warning against extrapolation are
% presented on p. 119. The experimental and comparison dataset is listed in
% Table 3 on p. 120. Fe3 effects and the authors' accuracy estimate are
% discussed on p. 122.
%
% APPROXIMATE DATASET ENVELOPE USED FOR NON-STOPPING SCREENING
%
%   Temperature : approximately 785-1500 degreeC
%                 Based on the observed/reference temperatures represented
%                 in Table 3 (p. 120).
%
%   Pressure    : approximately 1 bar-30 kbar
%                 Wood and Banno (1973) state that increasing pressure from
%                 1 bar to 30 kbar changes the equilibration temperature by
%                 50 degreeC or less at constant activity ratio, allowing
%                 pressure to be neglected for most purposes (p. 117).
%
%   Opx XFe     : approximately 0.055-0.85 for the direct multicomponent and
%                 Fe-bearing calibration/comparison data in Table 3.
%                 Endmember/simple-system data extend the broader dataset
%                 toward XFe = 0 and XFe = 1 (pp. 119-120).
%
% IMPORTANT:
% The numerical intervals above are practical envelopes represented by the
% source dataset and discussion. Wood and Banno (1973) did not define a
% single formal rectangular temperature-pressure-composition calibration
% field. These intervals are therefore used only for non-stopping screening
% warnings and must not be treated as strict validity boundaries.
%
% IMPORTANT APPLICATION LIMITATIONS
%
%   1) Equation (27) is semi-empirical. Wood and Banno (1973) state that the
%      Fe-correction constants were obtained empirically from the available
%      data and that the approach has little theoretical justification.
%      They explicitly caution that considerable errors may occur outside
%      the temperature-composition region covered by the derivation data
%      (p. 119).
%
%   2) The calculated Mg2Si2O6 component activities are "idealized"
%      activities based on ideal two-site mixing. The authors describe this
%      assumption as crude because a pyroxene solvus is present and state
%      that these quantities cannot be the true component activities
%      (p. 117).
%
%   3) The influence of Al-bearing pyroxene components was not well
%      constrained experimentally. Wood and Banno (1973) therefore assumed
%      ideal behavior for Al-pyroxene components in multicomponent Opx and
%      Cpx (p. 117). Strongly aluminous or otherwise unusual pyroxenes may
%      be outside the reliable behavior of the model.
%
%   4) The selected Opx and Cpx must represent a coexisting equilibrium
%      pair. Do not arbitrarily combine different textural generations,
%      cores and unrelated rims, porphyroclasts and neoblasts, or minerals
%      affected by different exsolution, reaction, alteration, or
%      metasomatic histories.
%
%   5) Inverted pigeonite and exsolved pyroxenes require special care.
%      Wood and Banno (1973) note that composition-dependent inversion
%      relations can cause calculated temperatures to underestimate the
%      actual equilibration temperature (p. 121).
%
%   6) The original natural-sample calculations assumed that all measured
%      Fe was Fe2. The authors estimate that assigning 10 percent of total
%      Fe in both pyroxenes to Fe3 would increase calculated temperatures
%      by approximately 20-30 degreeC (p. 122).
%
%      In this implementation, Fe_cation_apfu is treated as total Fe and:
%
%        Fe2_cation_apfu = Fe_cation_apfu - Fe3_cation_apfu
%
%      If Fe3_cation_apfu is absent, Fe3 is assumed to be zero. An explicitly
%      stored Fe3 NaN is retained and propagated rather than replaced by
%      zero.
%
%   7) Wood and Banno (1973) report that almost all experimental data are
%      reproduced within approximately 60 degreeC (abstract, p. 109;
%      p. 119) and estimate that calculated temperatures should in most
%      cases be accurate to about 70 degreeC (p. 122). This is an empirical
%      accuracy assessment rather than a formal analytical uncertainty.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure lies outside approximately 1 bar-30 kbar,
%   2) a finite calculated temperature lies outside approximately
%      785-1500 degreeC,
%   3) finite Opx XFe lies outside approximately 0.055-0.85,
%   4) an explicitly stored calculation input is NaN,
%   5) a derived site-allocation, activity, logarithm, denominator, or
%      temperature term is invalid, or
%   6) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
%
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%   rawdata_struct.Cpx : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. Remaining columns should contain
% consistently normalized pyroxene cations, preferably on a 6-oxygen basis.
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
% Optional variable for both Opx and Cpx:
%   Fe3_cation_apfu
%
% Existing NaN values are retained and propagated; they are never replaced
% by zero. Zero is assigned only when the optional Fe3 column itself is
% absent. All finite stored mineral-composition inputs must be greater than
% or equal to zero. Negative finite values and Inf stop the calculation. A
% derived Fe2 value below zero also stops the calculation.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and liquid NaN warnings is not applicable.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Wood and Banno (1973), Equation (27):
%
%                       -10202
%   T(K) = --------------------------------------------
%          ln(aEn_cpx/aEn_opx) - 7.65*XFe_opx
%          + 3.88*(XFe_opx)^2 - 4.6
%
% where:
%
%   XFe_opx = Fe2_opx / (Fe2_opx + Mg_opx)
%
% and the idealized Mg2Si2O6 component activity is:
%
%   aEn = XMg_M1 * XMg_M2
%
% SITE-ALLOCATION APPROXIMATION
%
% The blocking ions follow the assignments listed by Wood and Banno (1973)
% on p. 118:
%
%   M1 preferential occupants:
%     AlVI, Cr, Ti, Fe3
%
%   M2 preferential occupants:
%     Ca, Na, Mn
%
% Tetrahedral Al is estimated from Si deficiency:
%
%   AlIV = max(0, 2 - Si)
%   AlVI = Al_total - AlIV
%
% After the blocking ions are assigned, the remaining M1 and M2 capacities
% are filled by Mg and Fe2 using the same bulk Mg/(Mg + Fe2) ratio in both
% sites. This follows the random Fe-Mg distribution approximation described
% by Wood and Banno (1973) on p. 118.
%
% Temperature is returned in both Kelvin and degree Celsius:
%
%   T_degreeC = T_K - 273.15
%
% -------------------------------------------------------------------------
% Syntax:
%   results = WoodBanno1973(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Opx and Cpx tables
%   P_kbar         : finite non-negative numeric scalar or vector; stored in
%                    the output but not used explicitly by Equation (27)
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Opx-Cpx pair. T_K, T_degreeC, and T_deg are supplied for
%             downstream launcher and plotting compatibility.
%

%% Input validation
if nargin < 2
    error('WoodBanno1973 requires (rawdata_struct, P_kbar).');
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

%% 2) Initialize output container and screening settings
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Approximate source-dataset envelope, not a formal rectangular calibration
% field.
screeningT_min_degreeC = 785;
screeningT_max_degreeC = 1500;

% One bar equals 0.001 kbar. Wood and Banno (1973) discuss pressure effects
% from 1 bar to 30 kbar.
screeningP_min_kbar = 0.001;
screeningP_max_kbar = 30;

% Direct multicomponent and Fe-bearing range represented in Table 3.
screeningXFeOpx_min = 0.055;
screeningXFeOpx_max = 0.85;

pressureOutsideScreening = ...
    P_kbar < screeningP_min_kbar | ...
    P_kbar > screeningP_max_kbar;

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
        'PromptString', ...
            'Please select the Opx data you would like to use:', ...
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
        'PromptString', ...
            'Please select the Cpx data you would like to use:', ...
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

    nanInputNames = findNaNInputs( ...
        selectedData_opx, selectedData_cpx);

    validateNonNegativeInputs( ...
        selectedData_opx, selectedData_cpx);

    row = calcTemp( ...
        selectedData_opx, selectedData_cpx, P_kbar);

    nRows = height(row);
    row.dataCode_opx = repmat(selectedCode_opx, nRows, 1);
    row.dataCode_cpx = repmat(selectedCode_cpx, nRows, 1);
    row = movevars( ...
        row, {'dataCode_opx', 'dataCode_cpx'}, 'Before', 1);

    % Store one completed result block. The complete output table is not
    % enlarged on every interactive iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = ...
            [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % ----- Immediate result display -----
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    printTemperatureSummary( ...
        selectedCode_opx, selectedCode_cpx, row.T_degreeC);

    % Pressure is not used in Equation (27), but the paper discusses a
    % pressure-effect interval of approximately 1 bar-30 kbar.
    if any(pressureOutsideScreening) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate ' ...
             '1-bar to 30-kbar interval discussed by Wood and Banno ' ...
             '(1973, p. 117). %d of %d pressure point(s) are outside; ' ...
             'input range = %.6g-%.6g kbar. This is a practical ' ...
             'screening interval, not a formally stated rectangular ' ...
             'calibration range. Equation (27) contains no explicit ' ...
             'pressure term, and the supplied pressures remain in the ' ...
             'output table.\n'], ...
            sum(pressureOutsideScreening), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    printTemperatureScreeningWarning( ...
        row.T_degreeC, ...
        screeningT_min_degreeC, ...
        screeningT_max_degreeC, ...
        selectedCode_opx, selectedCode_cpx);

    printXFeScreeningWarning( ...
        row.XFe_opx, ...
        screeningXFeOpx_min, ...
        screeningXFeOpx_max, ...
        selectedCode_opx, selectedCode_cpx);

    % Existing NaN values are retained and never replaced by zero.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Wood-Banno Equation (27) ' ...
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
            ['WARNING: Invalid Wood-Banno Equation (27) site, activity, ' ...
             'fraction, logarithm, denominator, or temperature term(s) ' ...
             'were found for %s & %s: %s.\n' ...
             '         Affected values were retained as NaN, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_cpx), ...
            char(strjoin(invalidTerms, ', ')));
    end

    printNonFiniteTemperatureWarning( ...
        row.T_degreeC, selectedCode_opx, selectedCode_cpx);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'WoodBanno1973', ...
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
    'primaryTemperatureEquation', ...
        'Wood and Banno (1973) Equation (27)', ...
    'calibrationType', ...
        'Semi-empirical two-pyroxene thermometer', ...
    'pressureUsedInEquation', false, ...
    'approximateDatasetTemperatureRange_degreeC', ...
        [screeningT_min_degreeC, screeningT_max_degreeC], ...
    'approximatePressureDiscussionRange_kbar', ...
        [screeningP_min_kbar, screeningP_max_kbar], ...
    'approximateDirectMulticomponentXFeOpxRange', ...
        [screeningXFeOpx_min, screeningXFeOpx_max], ...
    'rangeStatus', ...
        ['Practical source-dataset screening envelopes; not a formal ' ...
         'rectangular calibration field'], ...
    'reportedTypicalAccuracy_degreeC', 70);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_opx, data_cpx)
% findNaNInputs
% Return names of explicitly stored pyroxene variables used by Equation (27)
% or its site-allocation procedure that contain NaN. Missing optional Fe3
% columns are not reported because they receive the documented default zero.

pyroxeneVariables = { ...
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

nameBuffer = strings(2 .* numel(pyroxeneVariables), 1);
nNames = 0;

for i = 1:numel(pyroxeneVariables)
    variableName = pyroxeneVariables{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        value = data_opx.(variableName);

        if isnumeric(value) && any(isnan(value(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Opx." + string(variableName);
        end
    end
end

for i = 1:numel(pyroxeneVariables)
    variableName = pyroxeneVariables{i};

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
% Stop when a stored pyroxene-composition value is negative or infinite.
% Zero is allowed. NaN is deliberately allowed and propagated.

pyroxeneVariables = { ...
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

nameBuffer = strings(2 .* numel(pyroxeneVariables), 1);
nInvalid = 0;

for i = 1:numel(pyroxeneVariables)
    variableName = pyroxeneVariables{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        value = data_opx.(variableName);
        validateScalarVariable(value, 'Opx', variableName);

        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            nameBuffer(nInvalid) = "Opx." + string(variableName);
        end
    end
end

for i = 1:numel(pyroxeneVariables)
    variableName = pyroxeneVariables{i};

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

    error(['WoodBanno1973: pyroxene-composition inputs must not be ' ...
           'negative or infinite. Invalid value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_cpx, P_kbar)
% calcTemp
% Calculate Wood and Banno (1973) Equation (27) for one selected Opx-Cpx
% pair and repeat the pressure-independent result for every supplied
% pressure. Existing NaN values and invalid derived terms are retained as
% NaN.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

opx = preparePyroxeneRow(data_opx, 'Opx');
cpx = preparePyroxeneRow(data_cpx, 'Cpx');

site_opx = calcSiteFractions(opx);
site_cpx = calcSiteFractions(cpx);

% Idealized enstatite-component activities.
if all(isfinite([site_opx.XMg_M1, site_opx.XMg_M2]))
    aEn_opx = site_opx.XMg_M1 .* site_opx.XMg_M2;
else
    aEn_opx = NaN;
end

if all(isfinite([site_cpx.XMg_M1, site_cpx.XMg_M2]))
    aEn_cpx = site_cpx.XMg_M1 .* site_cpx.XMg_M2;
else
    aEn_cpx = NaN;
end

if ~(isfinite(aEn_opx) && aEn_opx > 0)
    aEn_opx = NaN;
end
if ~(isfinite(aEn_cpx) && aEn_cpx > 0)
    aEn_cpx = NaN;
end

% Opx Fe fraction used by the empirical correction.
denom_XFe_opx = opx.Fe2 + opx.Mg;

if isfinite(denom_XFe_opx) && denom_XFe_opx > 0 && ...
        all(isfinite([opx.Fe2, opx.Mg]))

    XFe_opx = opx.Fe2 ./ denom_XFe_opx;
else
    XFe_opx = NaN;
end

if ~(isfinite(XFe_opx) && XFe_opx >= 0 && XFe_opx <= 1)
    XFe_opx = NaN;
end

% Activity ratio and logarithm.
if all(isfinite([aEn_opx, aEn_cpx])) && ...
        aEn_opx > 0 && aEn_cpx > 0

    aEn_ratio = aEn_cpx ./ aEn_opx;

    if isfinite(aEn_ratio) && aEn_ratio > 0
        ln_aEn_ratio = log(aEn_ratio);
    else
        aEn_ratio = NaN;
        ln_aEn_ratio = NaN;
    end
else
    aEn_ratio = NaN;
    ln_aEn_ratio = NaN;
end

% Wood and Banno (1973), Equation (27).
if all(isfinite([ln_aEn_ratio, XFe_opx]))
    denominator_T = ...
        ln_aEn_ratio ...
        - 7.65 .* XFe_opx ...
        + 3.88 .* (XFe_opx .^ 2) ...
        - 4.6;
else
    denominator_T = NaN;
end

if isfinite(denominator_T) && abs(denominator_T) > 1e-12
    T_scalar_raw_K = -10202 ./ denominator_T;
else
    T_scalar_raw_K = NaN;
end

T_scalar_raw_degreeC = T_scalar_raw_K - 273.15;
T_scalar_K = T_scalar_raw_K;
T_scalar_degreeC = T_scalar_raw_degreeC;

% Non-positive Kelvin is physically invalid and is retained as NaN.
if ~isfinite(T_scalar_K) || T_scalar_K <= 0
    T_scalar_K = NaN;
    T_scalar_degreeC = NaN;
end

row = table();

row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PressureUsedInEquation = false(nP, 1);
row.PrimaryEquation = repmat("WoodBanno1973_Eq27", nP, 1);

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

% Opx site-allocation diagnostics.
row.AlIV_raw_opx = repmat(site_opx.AlIV_raw, nP, 1);
row.AlIV_opx = repmat(site_opx.AlIV, nP, 1);
row.AlVI_raw_opx = repmat(site_opx.AlVI_raw, nP, 1);
row.AlVI_opx = repmat(site_opx.AlVI, nP, 1);
row.M1_fixed_raw_opx = repmat(site_opx.M1_fixed_raw, nP, 1);
row.M1_fixed_opx = repmat(site_opx.M1_fixed, nP, 1);
row.M2_fixed_raw_opx = repmat(site_opx.M2_fixed_raw, nP, 1);
row.M2_fixed_opx = repmat(site_opx.M2_fixed, nP, 1);
row.M1_remaining_opx = repmat(site_opx.M1_remaining, nP, 1);
row.M2_remaining_opx = repmat(site_opx.M2_remaining, nP, 1);
row.MgFe_total_opx = repmat(site_opx.MgFe_total, nP, 1);
row.Mg_fraction_opx = repmat(site_opx.Mg_fraction, nP, 1);
row.Fe2_fraction_opx = repmat(site_opx.Fe2_fraction, nP, 1);
row.Mg_M1_opx = repmat(site_opx.Mg_M1, nP, 1);
row.Mg_M2_opx = repmat(site_opx.Mg_M2, nP, 1);
row.Fe2_M1_opx = repmat(site_opx.Fe2_M1, nP, 1);
row.Fe2_M2_opx = repmat(site_opx.Fe2_M2, nP, 1);
row.M1_total_opx = repmat(site_opx.M1_total, nP, 1);
row.M2_total_opx = repmat(site_opx.M2_total, nP, 1);
row.XMg_M1_opx = repmat(site_opx.XMg_M1, nP, 1);
row.XMg_M2_opx = repmat(site_opx.XMg_M2, nP, 1);

% Cpx site-allocation diagnostics.
row.AlIV_raw_cpx = repmat(site_cpx.AlIV_raw, nP, 1);
row.AlIV_cpx = repmat(site_cpx.AlIV, nP, 1);
row.AlVI_raw_cpx = repmat(site_cpx.AlVI_raw, nP, 1);
row.AlVI_cpx = repmat(site_cpx.AlVI, nP, 1);
row.M1_fixed_raw_cpx = repmat(site_cpx.M1_fixed_raw, nP, 1);
row.M1_fixed_cpx = repmat(site_cpx.M1_fixed, nP, 1);
row.M2_fixed_raw_cpx = repmat(site_cpx.M2_fixed_raw, nP, 1);
row.M2_fixed_cpx = repmat(site_cpx.M2_fixed, nP, 1);
row.M1_remaining_cpx = repmat(site_cpx.M1_remaining, nP, 1);
row.M2_remaining_cpx = repmat(site_cpx.M2_remaining, nP, 1);
row.MgFe_total_cpx = repmat(site_cpx.MgFe_total, nP, 1);
row.Mg_fraction_cpx = repmat(site_cpx.Mg_fraction, nP, 1);
row.Fe2_fraction_cpx = repmat(site_cpx.Fe2_fraction, nP, 1);
row.Mg_M1_cpx = repmat(site_cpx.Mg_M1, nP, 1);
row.Mg_M2_cpx = repmat(site_cpx.Mg_M2, nP, 1);
row.Fe2_M1_cpx = repmat(site_cpx.Fe2_M1, nP, 1);
row.Fe2_M2_cpx = repmat(site_cpx.Fe2_M2, nP, 1);
row.M1_total_cpx = repmat(site_cpx.M1_total, nP, 1);
row.M2_total_cpx = repmat(site_cpx.M2_total, nP, 1);
row.XMg_M1_cpx = repmat(site_cpx.XMg_M1, nP, 1);
row.XMg_M2_cpx = repmat(site_cpx.XMg_M2, nP, 1);

% Activities, empirical correction, and temperature.
row.aEn_opx = repmat(aEn_opx, nP, 1);
row.aEn_cpx = repmat(aEn_cpx, nP, 1);
row.aEn_ratio = repmat(aEn_ratio, nP, 1);
row.ln_aEn_ratio = repmat(ln_aEn_ratio, nP, 1);
row.denom_XFe_opx = repmat(denom_XFe_opx, nP, 1);
row.XFe_opx = repmat(XFe_opx, nP, 1);
row.denominator_T = repmat(denominator_T, nP, 1);

row.T_Eq27_raw_K = repmat(T_scalar_raw_K, nP, 1);
row.T_Eq27_raw_degreeC = repmat(T_scalar_raw_degreeC, nP, 1);
row.T_K = repmat(T_scalar_K, nP, 1);
row.T_degreeC = repmat(T_scalar_degreeC, nP, 1);
row.T_deg = row.T_degreeC;

end

function px = preparePyroxeneRow(data_px, mineralLabel)
% preparePyroxeneRow
% Extract one pyroxene composition. Existing NaN values remain NaN. Missing
% Fe3 is assigned zero. Fe_cation_apfu is treated as total Fe.

if height(data_px) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

px = struct();

px.Si = getVarRequired( ...
    data_px, 'Si_cation_apfu', mineralLabel);
px.Al = getVarRequired( ...
    data_px, 'Al_cation_apfu', mineralLabel);
px.Fe_total = getVarRequired( ...
    data_px, 'Fe_cation_apfu', mineralLabel);
px.Mg = getVarRequired( ...
    data_px, 'Mg_cation_apfu', mineralLabel);
px.Ca = getVarRequired( ...
    data_px, 'Ca_cation_apfu', mineralLabel);
px.Na = getVarRequired( ...
    data_px, 'Na_cation_apfu', mineralLabel);
px.Mn = getVarRequired( ...
    data_px, 'Mn_cation_apfu', mineralLabel);
px.Ti = getVarRequired( ...
    data_px, 'Ti_cation_apfu', mineralLabel);
px.Cr = getVarRequired( ...
    data_px, 'Cr_cation_apfu', mineralLabel);

px.Fe3 = getVarOptional( ...
    data_px, 'Fe3_cation_apfu', 0, mineralLabel);

px.Fe2 = deriveFe2(px.Fe_total, px.Fe3, mineralLabel);

end

function Fe2 = deriveFe2(Fe_total, Fe3, mineralLabel)
% deriveFe2
% Derive Fe2 from total Fe and Fe3. NaN is propagated. A finite negative
% derived Fe2 value is prohibited.

Fe2 = Fe_total - Fe3;

if isfinite(Fe2) && Fe2 < 0
    error(['%s derived Fe2_cation_apfu is negative because ' ...
           'Fe3_cation_apfu exceeds Fe_cation_apfu.'], ...
        mineralLabel);
end
if isinf(Fe2)
    error('%s derived Fe2_cation_apfu is infinite.', mineralLabel);
end

end

function site = calcSiteFractions(px)
% calcSiteFractions
% Calculate the idealized two-site Mg fractions used by Wood and Banno
% (1973). Invalid derived values become NaN so that the output row is
% retained and can be diagnosed.

site = struct();

% Tetrahedral Al from Si deficiency. Si values slightly above 2 produce
% AlIV = 0, following the conventional approximation. NaN remains NaN.
if isfinite(px.Si)
    AlIV_raw = 2 - px.Si;

    if AlIV_raw < 0
        AlIV = 0;
    elseif AlIV_raw <= 2
        AlIV = AlIV_raw;
    else
        AlIV = NaN;
    end
else
    AlIV_raw = NaN;
    AlIV = NaN;
end

% Octahedral Al. Negative or excessive values indicate inconsistent cation
% normalization/site allocation and are retained as NaN.
if all(isfinite([px.Al, AlIV]))
    AlVI_raw = px.Al - AlIV;

    if AlVI_raw >= 0 && AlVI_raw <= 1
        AlVI = AlVI_raw;
    else
        AlVI = NaN;
    end
else
    AlVI_raw = NaN;
    AlVI = NaN;
end

% Fixed occupants of the M1 and M2 sites.
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
    M1_fixed_raw = NaN;
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
    M2_fixed_raw = NaN;
    M2_fixed = NaN;
    M2_remaining = NaN;
end

% Random Mg-Fe2 distribution over the remaining capacities.
MgFe_total = px.Mg + px.Fe2;

if isfinite(MgFe_total) && MgFe_total > 0 && ...
        all(isfinite([px.Mg, px.Fe2, M1_remaining, M2_remaining]))

    Mg_fraction = px.Mg ./ MgFe_total;
    Fe2_fraction = px.Fe2 ./ MgFe_total;

    Mg_M1 = M1_remaining .* Mg_fraction;
    Fe2_M1 = M1_remaining .* Fe2_fraction;

    Mg_M2 = M2_remaining .* Mg_fraction;
    Fe2_M2 = M2_remaining .* Fe2_fraction;

    M1_total = M1_fixed + Mg_M1 + Fe2_M1;
    M2_total = M2_fixed + Mg_M2 + Fe2_M2;

    if isfinite(M1_total) && M1_total > 0
        XMg_M1 = Mg_M1 ./ M1_total;
    else
        M1_total = NaN;
        XMg_M1 = NaN;
    end

    if isfinite(M2_total) && M2_total > 0
        XMg_M2 = Mg_M2 ./ M2_total;
    else
        M2_total = NaN;
        XMg_M2 = NaN;
    end
else
    Mg_fraction = NaN;
    Fe2_fraction = NaN;
    Mg_M1 = NaN;
    Fe2_M1 = NaN;
    Mg_M2 = NaN;
    Fe2_M2 = NaN;
    M1_total = NaN;
    M2_total = NaN;
    XMg_M1 = NaN;
    XMg_M2 = NaN;
end

if ~(isfinite(XMg_M1) && XMg_M1 >= 0 && XMg_M1 <= 1)
    XMg_M1 = NaN;
end
if ~(isfinite(XMg_M2) && XMg_M2 >= 0 && XMg_M2 <= 1)
    XMg_M2 = NaN;
end

site.AlIV_raw = AlIV_raw;
site.AlIV = AlIV;
site.AlVI_raw = AlVI_raw;
site.AlVI = AlVI;

site.M1_fixed_raw = M1_fixed_raw;
site.M1_fixed = M1_fixed;
site.M2_fixed_raw = M2_fixed_raw;
site.M2_fixed = M2_fixed;

site.M1_remaining = M1_remaining;
site.M2_remaining = M2_remaining;

site.MgFe_total = MgFe_total;
site.Mg_fraction = Mg_fraction;
site.Fe2_fraction = Fe2_fraction;

site.Mg_M1 = Mg_M1;
site.Fe2_M1 = Fe2_M1;
site.Mg_M2 = Mg_M2;
site.Fe2_M2 = Fe2_M2;

site.M1_total = M1_total;
site.M2_total = M2_total;
site.XMg_M1 = XMg_M1;
site.XMg_M2 = XMg_M2;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid site-allocation, activity, fraction, logarithm,
% denominator, and temperature terms.

termBuffer = strings(32, 1);
nTerms = 0;

checkFiniteBounded = { ...
    'Opx AlIV', row.AlIV_opx(1), 0, 2; ...
    'Opx AlVI', row.AlVI_opx(1), 0, 1; ...
    'Opx M1 fixed occupancy', row.M1_fixed_opx(1), 0, 1; ...
    'Opx M2 fixed occupancy', row.M2_fixed_opx(1), 0, 1; ...
    'Opx XMg_M1', row.XMg_M1_opx(1), 0, 1; ...
    'Opx XMg_M2', row.XMg_M2_opx(1), 0, 1; ...
    'Cpx AlIV', row.AlIV_cpx(1), 0, 2; ...
    'Cpx AlVI', row.AlVI_cpx(1), 0, 1; ...
    'Cpx M1 fixed occupancy', row.M1_fixed_cpx(1), 0, 1; ...
    'Cpx M2 fixed occupancy', row.M2_fixed_cpx(1), 0, 1; ...
    'Cpx XMg_M1', row.XMg_M1_cpx(1), 0, 1; ...
    'Cpx XMg_M2', row.XMg_M2_cpx(1), 0, 1; ...
    'Opx XFe', row.XFe_opx(1), 0, 1};

for i = 1:size(checkFiniteBounded, 1)
    label = checkFiniteBounded{i, 1};
    value = checkFiniteBounded{i, 2};
    minimumValue = checkFiniteBounded{i, 3};
    maximumValue = checkFiniteBounded{i, 4};

    if ~isfinite(value) || ...
            value < minimumValue || value > maximumValue

        nTerms = nTerms + 1;
        termBuffer(nTerms) = string(label);
    end
end

checkPositive = { ...
    'Opx Mg+Fe2 denominator', row.MgFe_total_opx(1); ...
    'Cpx Mg+Fe2 denominator', row.MgFe_total_cpx(1); ...
    'Opx idealized aEn', row.aEn_opx(1); ...
    'Cpx idealized aEn', row.aEn_cpx(1); ...
    'aEn_cpx/aEn_opx', row.aEn_ratio(1)};

for i = 1:size(checkPositive, 1)
    label = checkPositive{i, 1};
    value = checkPositive{i, 2};

    if ~isfinite(value) || value <= 0
        nTerms = nTerms + 1;
        termBuffer(nTerms) = string(label);
    end
end

if ~isfinite(row.ln_aEn_ratio(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "ln(aEn_cpx/aEn_opx)";
end

if ~isfinite(row.denominator_T(1)) || ...
        abs(row.denominator_T(1)) <= 1e-12

    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Equation (27) denominator";
end

if any(~isfinite(row.T_K) | row.T_K <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "T_K";
end

invalidTerms = unique(termBuffer(1:nTerms), 'stable');

end

function printTemperatureSummary( ...
        selectedCode_opx, selectedCode_cpx, temperatureValues)
% printTemperatureSummary
% Display one temperature or first-to-last values for a pressure vector.

label = [char(selectedCode_opx) ' & ' char(selectedCode_cpx)];

if isscalar(temperatureValues)
    disp([label ': ' num2str(temperatureValues) ' degreeC']);
else
    disp([label ': ' ...
        num2str(temperatureValues(1)) ' to ' ...
        num2str(temperatureValues(end)) ' degreeC']);
end

end

function printTemperatureScreeningWarning( ...
        temperatureValues, minimumTemperature, maximumTemperature, ...
        selectedCode_opx, selectedCode_cpx)
% printTemperatureScreeningWarning
% Warn when finite temperatures lie outside the approximate Table 3 dataset
% envelope. Results are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);

    fprintf(2, ...
        ['WARNING: Calculated Wood-Banno temperature is outside the ' ...
         'approximate %.4g-%.4g degreeC source-dataset envelope ' ...
         'represented in Table 3 of Wood and Banno (1973, p. 120). ' ...
         '%d of %d finite point(s) are outside; calculated finite range ' ...
         '= %.6g-%.6g degreeC for %s & %s. This is a practical ' ...
         'screening envelope, not a formally stated rectangular ' ...
         'calibration range. The result has been retained.\n'], ...
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

function printXFeScreeningWarning( ...
        XFeValues, minimumXFe, maximumXFe, ...
        selectedCode_opx, selectedCode_cpx)
% printXFeScreeningWarning
% Warn when finite Opx XFe lies outside the approximate direct
% multicomponent/Fe-bearing range represented in Table 3. Endmember data
% extend toward zero and one, so the result is retained.

finiteMask = isfinite(XFeValues);
outsideMask = finiteMask & ...
    (XFeValues < minimumXFe | XFeValues > maximumXFe);

if any(outsideMask)
    finiteValues = XFeValues(finiteMask);

    fprintf(2, ...
        ['CAUTION: Opx XFe is outside the approximate %.4g-%.4g range ' ...
         'represented by the direct multicomponent and Fe-bearing ' ...
         'calibration/comparison data in Wood and Banno (1973, Table 3, ' ...
         'p. 120). Calculated finite XFe range = %.6g-%.6g for %s & %s. ' ...
         'Simple-system/endmember data extend toward XFe = 0 and 1, so ' ...
         'this message is a composition-screening caution rather than a ' ...
         'strict rejection. The result has been retained.\n'], ...
        minimumXFe, ...
        maximumXFe, ...
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
        ['WARNING: Non-finite Wood-Banno Equation (27) temperature ' ...
         'values were calculated for %s & %s (%d of %d points; ' ...
         'NaN: %d, Inf: %d).\n' ...
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

function validateRequiredColumns(tbl, mineralLabel)
% validateRequiredColumns
% Verify normalized cation columns required for Wood-Banno site allocation.

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
    error('%s table must contain variable: %s', ...
        mineralLabel, variableName);
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
% are handled by validateNonNegativeInputs.

if ~isnumeric(value) || ~isscalar(value)
    error('%s variable %s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end
