function results = Putirka2008CpxLiq(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Pyroxene/Putirka2008CpxLiq.m
% Tested with MATLAB R2024b
%
% Clinopyroxene-Liquid thermometer
% Putirka, K.D. (2008), Equation (33)
% Reviews in Mineralogy and Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Clinopyroxene analysis from
% rawdata_struct.Cpx and combines it with one row from a liquid dataset
% loaded by liquid.readLiquidExcel(). Temperature is calculated using the
% globally calibrated Clinopyroxene-Liquid thermometer of Putirka (2008),
% Equation (33).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Clinopyroxene-Liquid pair, one output row is returned for
% every pressure value supplied in P_kbar. It is therefore compatible with
% both startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka (2008) presents Equation (33) on p. 94. The equation is a global
% calibration of the Clinopyroxene-Liquid Jd-DiHd exchange thermometer using
% experiments conducted at pressures below 70 kbar (p. 94; Figure 9a on
% p. 92).
%
% Figure 9a reports the following validation statistics for Equation (33):
%   All data, P < 70 kbar : SEE = +/-45 degreeC, R^2 = 0.92, n = 1174
%   Anhydrous experiments : SEE = +/-46 degreeC, R^2 = 0.90, n = 854
%   Hydrous experiments   : SEE = +/-42 degreeC, R^2 = 0.86, n = 320
%
% Putirka (2008) does not list one strict minimum-maximum temperature box for
% Equation (33). Figure 9a displays the calibration/validation comparison
% over approximately 600-2400 degreeC. This graphical interval is used only
% as a non-stopping screening range; it is not presented here as a formal
% universal calibration boundary.
%
% Important application cautions:
%   1) A texturally and compositionally plausible equilibrium Cpx-Liquid
%      pair must be used. Pairing a crystal with a non-equilibrium liquid can
%      produce a formally calculated but geologically meaningless result.
%
%   2) Fe-Mg exchange alone does not prove Na-Al or Ca-Na equilibrium.
%      Nevertheless, Putirka (2008) notes that deviations in Fe-Mg exchange
%      correlate with deviations in DiHd and EnFs components. The experimental
%      mean is KD(Fe-Mg)cpx-liq = 0.28 +/- 0.08, with a total observed range
%      of 0.04-0.68; the temperature-dependent relation is given by Eq. (35)
%      on p. 94.
%
%   3) The liquid-only Cpx-saturation temperature from Eq. (34) can be
%      compared with the Eq. (33) result as an additional equilibrium and
%      saturation check (p. 94).
%
%   4) Clinopyroxene components must be calculated on a 6-oxygen basis using
%      the same allocation procedure as the original calibration. The Cpx
%      component scheme is given in Table 3 on pp. 88-89. A good pyroxene
%      analysis should have a total cation sum close to 4.
%
%   5) Liquid cation fractions must be calculated in the same way as in the
%      original calibration: on an anhydrous basis, without first
%      renormalizing oxide wt% to 100%. H2O is entered separately in wt%
%      (component-calculation discussion on pp. 68-71).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure does not satisfy the published P < 70 kbar condition,
%   2) a finite calculated temperature lies outside the approximate
%      600-2400 degreeC Figure-9a comparison envelope,
%   3) an explicitly stored calculation input is NaN,
%   4) a logarithm, ratio, or denominator required by Eq. (33) is invalid,
%      or
%   5) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST column of the Cpx table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% Clinopyroxene oxide columns:
%   Required for normalization:
%     SiO2, MgO, CaO, and FeO or FeOt
%
%   Optional; an absent column is assigned zero:
%     TiO2, Al2O3, MnO, Na2O, K2O, Cr2O3
%
% The liquid table is loaded by liquid.readLiquidExcel(). The following
% oxide names are recognized, including common variants such as
% "SiO2value", "SiO2 value", and "SiO2_value":
%
%   SiO2, TiO2, Al2O3, FeO or FeOt, MnO, MgO, CaO, Na2O, K2O,
%   V2O3, Cr2O3, NiO, P2O5, SO3, Fe2O3, H2O or H2Ot, F, Cl
%
% Missing optional columns are assigned zero. An explicitly stored NaN is
% retained and propagated; it is never replaced by zero. All finite
% calculation inputs must be greater than or equal to zero. Negative finite
% values and Inf stop the calculation.
%
% F and Cl are deliberately excluded from cationTotal_liq because they are
% not oxide cation components in the Eq. (33) anhydrous normalization used
% here. F and Cl are also excluded from the NaN-input warning and from the
% negative-input validation because they do not enter the calculation.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Putirka (2008), Equation (33), p. 94:
%
%   10^4 / T(K) =
%       7.53
%     - 0.14 * ln[(XJd_cpx * XCaO_liq * XFm_liq) /
%                 (XDiHd_cpx * XNaO0.5_liq * XAlO1.5_liq)]
%     + 0.07 * H2O_liq
%     - 14.9 * XCaO_liq * XSiO2_liq
%     - 0.08 * ln(XTiO2_liq)
%     - 3.62 * (XNaO0.5_liq + XKO0.5_liq)
%     - 1.1 * MgNum_liq
%     - 0.18 * ln(XEnFs_cpx)
%     - 0.027 * P_kbar
%
% where:
%   XFm_liq   = XFeO_liq + XMnO_liq + XMgO_liq
%   MgNum_liq = XMgO_liq / (XMgO_liq + XFeO_liq)
%
% H2O_liq is wt% H2O and is excluded from cationTotal_liq.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008CpxLiq(rawdata_struct, P_kbar)
%   results = Putirka2008CpxLiq(rawdata_struct, P_kbar, ...
%       'LiquidRow', n)
%
% Inputs:
%   rawdata_struct : struct containing a Clinopyroxene table
%   P_kbar         : pressure in kbar; finite, non-negative numeric scalar
%                    or vector
%
% Name-value option:
%   LiquidRow      : positive integer row number in the loaded liquid table.
%                    If empty, row 1 is used.
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Clinopyroxene-Liquid pair. T33_K and T33_C retain the original
%             equation-specific names. T_K, T_degreeC, and T_deg are
%             standardized aliases for launcher and plotting compatibility.
%

%% Input validation
% Accept scalar or vector pressure input from both fixed-P and range-P
% launchers.
if nargin < 2
    error('Putirka2008CpxLiq requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end
if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end
if isempty(rawdata_struct.Cpx)
    error('rawdata_struct.Cpx is empty.');
end

P_kbar = P_kbar(:);

%% Options
ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve Clinopyroxene and liquid datasets
disp('=== Step 1: Preparing Clinopyroxene and liquid datasets ===');

dataset_cpx = rawdata_struct.Cpx;
MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();

if isempty(liqAll) || ~istable(liqAll)
    error('Selected liquid dataset is empty or is not a table.');
end

% These columns are necessary to establish a 6-oxygen Cpx normalization.
validateRequiredColumns(dataset_cpx, ...
    {{'SiO2'}, {'MgO'}, {'CaO'}, {'FeO', 'FeOt'}}, 'Clinopyroxene');

disp('=== Preparing Clinopyroxene and liquid datasets has been finished ===');

%% 2) Select liquid row and initialize output container
disp('=== Step 2: Preparing liquid selection and output container ===');

if isempty(liquidRowOpt)
    idxLiq = 1;
    if height(liqAll) > 1
        fprintf(2, ...
            ['WARNING: The selected liquid dataset contains %d rows. ' ...
             'Liquid row 1 is being used. Specify ''LiquidRow'' to select ' ...
             'a different row.\n'], ...
            height(liqAll));
    end
else
    idxLiq = liquidRowOpt;
    if idxLiq > height(liqAll)
        error(['Requested LiquidRow (%d) exceeds the number of rows in the ' ...
               'selected liquid dataset (%d).'], idxLiq, height(liqAll));
    end
end

selectedData_liq = liqAll(idxLiq, :);

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Equation (33) global calibration uses experiments at P < 70 kbar.
calibrationP_max_kbar = 70;
pressureOutsideCalibration = P_kbar >= calibrationP_max_kbar;
pressureWarningIssued = false;

% Approximate graphical comparison envelope shown in Figure 9a.
screeningT_min_degC = 600;
screeningT_max_degC = 2400;

% Missing H2O is an explicit zero-water assumption, not a measured value.
[~, ~, h2oColumnFound] = getOxideValue( ...
    selectedData_liq, {'H2O', 'H2Ot'}, 0, false, 'Liquid');
if ~h2oColumnFound
    fprintf(2, ...
        ['WARNING: No H2O or H2Ot column was found in the selected liquid ' ...
         'dataset. H2O_liq = 0 wt%% is being assumed.\n']);
end

disp(['Liquid selected: Row ' num2str(idxLiq)]);
disp('=== Preparing liquid selection and output container has been finished ===');

%% 3-4) Interactive Clinopyroxene selection and calculation
disp('=== Step 3: Selecting a data code from the list (Cpx) ===');

dataCodes_cpx = dataset_cpx{:, 1};
displayCodes_cpx = cellstr(string(dataCodes_cpx));

while true
    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Cpx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_cpx, ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = string(dataCodes_cpx(selectedIdx_cpx));
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    disp(['Cpx selected: ' char(selectedCode_cpx)]);
    disp('=== Step 4: Calculating the temperature ===');

    % Identify explicitly stored NaN values in active calculation inputs.
    % F and Cl are intentionally absent from this check.
    nanInputNames = findNaNInputs(selectedData_cpx, selectedData_liq);

    % Reject negative finite values and Inf. NaN and zero are allowed so
    % that NaN can propagate and zero-dependent invalid logs can be reported
    % without silently changing the input.
    validateNonNegativeInputs(selectedData_cpx, selectedData_liq);

    row = calcTemp(selectedData_cpx, selectedData_liq, P_kbar, MWinfo);

    nRows = height(row);
    row.dataCode_cpx = repmat(selectedCode_cpx, nRows, 1);
    row.dataRow_liq = repmat(idxLiq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_cpx', 'dataRow_liq'}, 'Before', 1);

    % Store one completed result block. The full output table is not resized
    % on every interactive iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    printTemperatureSummary('Eq. (33)', row.T33_C);

    % Pressure is common to all selected Cpx rows, so warn only once.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental calibration ' ...
             'condition used for Putirka (2008) Equation (33): P < 70 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; input range ' ...
             '= %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn for finite temperatures outside the Figure-9a display envelope.
    printTemperatureRangeWarning( ...
        row.T33_C, screeningT_min_degC, screeningT_max_degC, ...
        selectedCode_cpx, idxLiq);

    % Report exact names of explicitly NaN inputs. Values remain NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s ' ...
             'and Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; the ' ...
             'calculated temperature may be NaN.\n'], ...
            char(selectedCode_cpx), ...
            idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Give a specific warning for invalid logarithm, ratio, or denominator
    % terms. The result is retained as NaN rather than stopping.
    invalidEquationTerms = findInvalidEquationTerms(row);
    if ~isempty(invalidEquationTerms)
        fprintf(2, ...
            ['WARNING: Invalid Equation (33) term(s) were found for %s and ' ...
             'Liquid row %d: %s.\n' ...
             '         Required logarithm arguments, ratios, and the final ' ...
             '10^4/T denominator must be finite and > 0. Affected ' ...
             'temperatures were retained as NaN.\n'], ...
            char(selectedCode_cpx), ...
            idxLiq, ...
            char(strjoin(invalidEquationTerms, ', ')));
    end

    printNonFiniteTemperatureWarning( ...
        row.T33_C, selectedCode_cpx, idxLiq);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Cpx selection (same Liquid dataset)?', ...
        'Putirka2008CpxLiq', ...
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
    'liquid', metaLiq, ...
    'primaryTemperatureEquation', 'Putirka2008 Eq. (33)', ...
    'liquidRow', idxLiq, ...
    'normalizationNote', ...
    'F and Cl are excluded from cationTotal_liq and from NaN-input warnings.');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_cpx, data_liquid)
% findNaNInputs
% Return names of explicitly stored calculation inputs whose selected value
% is NaN. Missing optional columns are not listed because they are assigned
% zero only when the column itself is absent. F and Cl are excluded.

cpxAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'Cr2O3'}};
cpxLabels = [ ...
    "Cpx.SiO2", "Cpx.TiO2", "Cpx.Al2O3", "Cpx.FeO/FeOt", ...
    "Cpx.MnO", "Cpx.MgO", "Cpx.CaO", "Cpx.Na2O", ...
    "Cpx.K2O", "Cpx.Cr2O3"];

liqAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'V2O3'}, {'Cr2O3'}, ...
    {'NiO'}, {'P2O5'}, {'SO3'}, {'Fe2O3'}, {'H2O', 'H2Ot'}};
liqLabels = [ ...
    "Liquid.SiO2", "Liquid.TiO2", "Liquid.Al2O3", ...
    "Liquid.FeO/FeOt", "Liquid.MnO", "Liquid.MgO", ...
    "Liquid.CaO", "Liquid.Na2O", "Liquid.K2O", ...
    "Liquid.V2O3", "Liquid.Cr2O3", "Liquid.NiO", ...
    "Liquid.P2O5", "Liquid.SO3", "Liquid.Fe2O3", ...
    "Liquid.H2O/H2Ot"];

maxNames = numel(cpxAliases) + numel(liqAliases);
nanNamesBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(cpxAliases)
    [value, ~, columnFound] = getOxideValue( ...
        data_cpx, cpxAliases{i}, 0, false, 'Cpx');
    if columnFound && isnan(value)
        nNames = nNames + 1;
        nanNamesBuffer(nNames) = cpxLabels(i);
    end
end

for i = 1:numel(liqAliases)
    [value, ~, columnFound] = getOxideValue( ...
        data_liquid, liqAliases{i}, 0, false, 'Liquid');
    if columnFound && isnan(value)
        nNames = nNames + 1;
        nanNamesBuffer(nNames) = liqLabels(i);
    end
end

nanInputNames = nanNamesBuffer(1:nNames);

end

function validateNonNegativeInputs(data_cpx, data_liquid)
% validateNonNegativeInputs
% Stop when an explicitly stored calculation input is negative or infinite.
% Zero is allowed. NaN is deliberately allowed and propagated. F and Cl are
% excluded because they do not enter cationTotal_liq or Equation (33).

cpxAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'Cr2O3'}};
cpxLabels = [ ...
    "Cpx.SiO2", "Cpx.TiO2", "Cpx.Al2O3", "Cpx.FeO/FeOt", ...
    "Cpx.MnO", "Cpx.MgO", "Cpx.CaO", "Cpx.Na2O", ...
    "Cpx.K2O", "Cpx.Cr2O3"];

liqAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'V2O3'}, {'Cr2O3'}, ...
    {'NiO'}, {'P2O5'}, {'SO3'}, {'Fe2O3'}, {'H2O', 'H2Ot'}};
liqLabels = [ ...
    "Liquid.SiO2", "Liquid.TiO2", "Liquid.Al2O3", ...
    "Liquid.FeO/FeOt", "Liquid.MnO", "Liquid.MgO", ...
    "Liquid.CaO", "Liquid.Na2O", "Liquid.K2O", ...
    "Liquid.V2O3", "Liquid.Cr2O3", "Liquid.NiO", ...
    "Liquid.P2O5", "Liquid.SO3", "Liquid.Fe2O3", ...
    "Liquid.H2O/H2Ot"];

maxNames = numel(cpxAliases) + numel(liqAliases);
invalidNamesBuffer = strings(maxNames, 1);
nInvalid = 0;

for i = 1:numel(cpxAliases)
    [value, ~, columnFound] = getOxideValue( ...
        data_cpx, cpxAliases{i}, 0, false, 'Cpx');
    if columnFound && (isinf(value) || (isfinite(value) && value < 0))
        nInvalid = nInvalid + 1;
        invalidNamesBuffer(nInvalid) = cpxLabels(i);
    end
end

for i = 1:numel(liqAliases)
    [value, ~, columnFound] = getOxideValue( ...
        data_liquid, liqAliases{i}, 0, false, 'Liquid');
    if columnFound && (isinf(value) || (isfinite(value) && value < 0))
        nInvalid = nInvalid + 1;
        invalidNamesBuffer(nInvalid) = liqLabels(i);
    end
end

if nInvalid > 0
    invalidInputNames = invalidNamesBuffer(1:nInvalid);
    error(['Putirka2008CpxLiq: calculation inputs must not be negative ' ...
           'or infinite. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_cpx, data_liquid, P_kbar, MWinfo)
% calcTemp
% Calculate Putirka (2008) Equation (33) for one selected Cpx-Liquid pair
% and a scalar or vector of pressures. Output rows correspond one-to-one
% with input pressure values. NaN inputs are retained and propagated.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

cpx = prepareCpxRow(data_cpx, MWinfo);
liq = prepareLiquidRow(data_liquid, MWinfo);

% Terms that require positive finite arguments are converted to NaN when
% invalid. This prevents complex values and retains the failed result.
exchangeTerm = (cpx.XJd .* liq.XCaO .* liq.XFm) ./ ...
    (cpx.XDiHd .* liq.XNaO0_5 .* liq.XAlO1_5);
lnExchangeTerm = safeLogPositive(exchangeTerm);
lnXTiO2_liq = safeLogPositive(liq.XTiO2);
lnXEnFs_cpx = safeLogPositive(cpx.XEnFs);

mgDenominator = liq.XMgO + liq.XFeO;
if isfinite(liq.XMgO) && isfinite(liq.XFeO) && ...
        isfinite(mgDenominator) && mgDenominator > 0
    MgNum_liq = liq.XMgO ./ mgDenominator;
else
    MgNum_liq = NaN;
end

% Putirka (2008), Equation (33), p. 94. The H2O coefficient is positive.
invT = ...
    7.53 ...
    - 0.14 .* lnExchangeTerm ...
    + 0.07 .* liq.H2O ...
    - 14.9 .* liq.XCaO .* liq.XSiO2 ...
    - 0.08 .* lnXTiO2_liq ...
    - 3.62 .* (liq.XNaO0_5 + liq.XKO0_5) ...
    - 1.1 .* MgNum_liq ...
    - 0.18 .* lnXEnFs_cpx ...
    - 0.027 .* P_kbar;

% A finite positive value of 10^4/T is required for a physical Kelvin
% temperature. Invalid values are represented by NaN and reported later.
T33_K = nan(nP, 1);
validInvT = isfinite(invT) & invT > 0;
T33_K(validInvT) = 10000 ./ invT(validInvT);
T33_C = T33_K - 273.15;

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PrimaryEquation = repmat("Putirka2008_Eq33", nP, 1);

% Selected Clinopyroxene oxides.
row.SiO2_cpx = repmat(cpx.SiO2, nP, 1);
row.TiO2_cpx = repmat(cpx.TiO2, nP, 1);
row.Al2O3_cpx = repmat(cpx.Al2O3, nP, 1);
row.FeO_cpx = repmat(cpx.FeO, nP, 1);
row.MnO_cpx = repmat(cpx.MnO, nP, 1);
row.MgO_cpx = repmat(cpx.MgO, nP, 1);
row.CaO_cpx = repmat(cpx.CaO, nP, 1);
row.Na2O_cpx = repmat(cpx.Na2O, nP, 1);
row.K2O_cpx = repmat(cpx.K2O, nP, 1);
row.Cr2O3_cpx = repmat(cpx.Cr2O3, nP, 1);

% Clinopyroxene cations and components on a 6-oxygen basis.
row.XSi_cpx = repmat(cpx.XSi, nP, 1);
row.XTi_cpx = repmat(cpx.XTi, nP, 1);
row.XAl_cpx = repmat(cpx.XAl, nP, 1);
row.XFe_cpx = repmat(cpx.XFe, nP, 1);
row.XMn_cpx = repmat(cpx.XMn, nP, 1);
row.XMg_cpx = repmat(cpx.XMg, nP, 1);
row.XCa_cpx = repmat(cpx.XCa, nP, 1);
row.XNa_cpx = repmat(cpx.XNa, nP, 1);
row.XK_cpx = repmat(cpx.XK, nP, 1);
row.XCr_cpx = repmat(cpx.XCr, nP, 1);
row.cationSum_cpx = repmat(cpx.cationSum, nP, 1);

row.XAlIV_cpx = repmat(cpx.XAlIV, nP, 1);
row.XAlVI_cpx = repmat(cpx.XAlVI, nP, 1);
row.XFe3_cpx = repmat(cpx.XFe3, nP, 1);
row.XJd_cpx = repmat(cpx.XJd, nP, 1);
row.XCaTs_cpx = repmat(cpx.XCaTs, nP, 1);
row.XCaTi_cpx = repmat(cpx.XCaTi, nP, 1);
row.XCrCaTs_cpx = repmat(cpx.XCrCaTs, nP, 1);
row.XDiHd_cpx = repmat(cpx.XDiHd, nP, 1);
row.XEnFs_cpx = repmat(cpx.XEnFs, nP, 1);

% Selected liquid oxides. F and Cl are retained only for traceability and
% do not enter cationTotal_liq, validation, or NaN-input warnings.
row.SiO2_liq = repmat(liq.SiO2, nP, 1);
row.TiO2_liq = repmat(liq.TiO2, nP, 1);
row.Al2O3_liq = repmat(liq.Al2O3, nP, 1);
row.FeO_liq = repmat(liq.FeO, nP, 1);
row.MnO_liq = repmat(liq.MnO, nP, 1);
row.MgO_liq = repmat(liq.MgO, nP, 1);
row.CaO_liq = repmat(liq.CaO, nP, 1);
row.Na2O_liq = repmat(liq.Na2O, nP, 1);
row.K2O_liq = repmat(liq.K2O, nP, 1);
row.V2O3_liq = repmat(liq.V2O3, nP, 1);
row.Cr2O3_liq = repmat(liq.Cr2O3, nP, 1);
row.NiO_liq = repmat(liq.NiO, nP, 1);
row.P2O5_liq = repmat(liq.P2O5, nP, 1);
row.SO3_liq = repmat(liq.SO3, nP, 1);
row.Fe2O3_liq = repmat(liq.Fe2O3, nP, 1);
row.F_liq = repmat(liq.F, nP, 1);
row.Cl_liq = repmat(liq.Cl, nP, 1);
row.H2O_liq = repmat(liq.H2O, nP, 1);

% Liquid cation fractions. F, Cl, and H2O are excluded from the denominator.
row.cationTotal_liq = repmat(liq.cationTotal, nP, 1);
row.XSiO2_liq = repmat(liq.XSiO2, nP, 1);
row.XTiO2_liq = repmat(liq.XTiO2, nP, 1);
row.XAlO1_5_liq = repmat(liq.XAlO1_5, nP, 1);
row.XFeO_liq = repmat(liq.XFeO, nP, 1);
row.XMnO_liq = repmat(liq.XMnO, nP, 1);
row.XMgO_liq = repmat(liq.XMgO, nP, 1);
row.XCaO_liq = repmat(liq.XCaO, nP, 1);
row.XNaO0_5_liq = repmat(liq.XNaO0_5, nP, 1);
row.XKO0_5_liq = repmat(liq.XKO0_5, nP, 1);
row.XFm_liq = repmat(liq.XFm, nP, 1);
row.MgNum_liq = repmat(MgNum_liq, nP, 1);

% Equation terms and final temperature.
row.exchangeTerm_Eq33 = repmat(exchangeTerm, nP, 1);
row.lnExchangeTerm_Eq33 = repmat(lnExchangeTerm, nP, 1);
row.lnXTiO2_liq = repmat(lnXTiO2_liq, nP, 1);
row.lnXEnFs_cpx = repmat(lnXEnFs_cpx, nP, 1);
row.invT_Eq33 = invT;
row.T33_K = T33_K;
row.T33_C = T33_C;

% Standardized launcher and plotting aliases.
row.T_K = T33_K;
row.T_degreeC = T33_C;
row.T_deg = T33_C;

end

function cpx = prepareCpxRow(data_cpx, MWinfo)
% prepareCpxRow
% Calculate Clinopyroxene cations and components on a 6-oxygen basis.
% Existing NaN values remain NaN and propagate through the normalization.

[cpx.SiO2, ~, ~] = getOxideValue( ...
    data_cpx, {'SiO2'}, NaN, true, 'Cpx');
[cpx.TiO2, ~, ~] = getOxideValue( ...
    data_cpx, {'TiO2'}, 0, false, 'Cpx');
[cpx.Al2O3, ~, ~] = getOxideValue( ...
    data_cpx, {'Al2O3'}, 0, false, 'Cpx');
[cpx.FeO, ~, ~] = getOxideValue( ...
    data_cpx, {'FeO', 'FeOt'}, NaN, true, 'Cpx');
[cpx.MnO, ~, ~] = getOxideValue( ...
    data_cpx, {'MnO'}, 0, false, 'Cpx');
[cpx.MgO, ~, ~] = getOxideValue( ...
    data_cpx, {'MgO'}, NaN, true, 'Cpx');
[cpx.CaO, ~, ~] = getOxideValue( ...
    data_cpx, {'CaO'}, NaN, true, 'Cpx');
[cpx.Na2O, ~, ~] = getOxideValue( ...
    data_cpx, {'Na2O'}, 0, false, 'Cpx');
[cpx.K2O, ~, ~] = getOxideValue( ...
    data_cpx, {'K2O'}, 0, false, 'Cpx');
[cpx.Cr2O3, ~, ~] = getOxideValue( ...
    data_cpx, {'Cr2O3'}, 0, false, 'Cpx');

mol = struct();
mol.SiO2 = cpx.SiO2 ./ MWinfo.MW.SiO2;
mol.TiO2 = cpx.TiO2 ./ MWinfo.MW.TiO2;
mol.Al2O3 = cpx.Al2O3 ./ MWinfo.MW.Al2O3;
mol.FeO = cpx.FeO ./ MWinfo.MW.FeO;
mol.MnO = cpx.MnO ./ MWinfo.MW.MnO;
mol.MgO = cpx.MgO ./ MWinfo.MW.MgO;
mol.CaO = cpx.CaO ./ MWinfo.MW.CaO;
mol.Na2O = cpx.Na2O ./ MWinfo.MW.Na2O;
mol.K2O = cpx.K2O ./ MWinfo.MW.K2O;
mol.Cr2O3 = cpx.Cr2O3 ./ MWinfo.MW.Cr2O3;

oxySum = ...
    2 .* mol.SiO2 ...
    + 2 .* mol.TiO2 ...
    + 3 .* mol.Al2O3 ...
    + mol.FeO ...
    + mol.MnO ...
    + mol.MgO ...
    + mol.CaO ...
    + mol.Na2O ...
    + mol.K2O ...
    + 3 .* mol.Cr2O3;

if isfinite(oxySum) && oxySum > 0
    ORF = 6 ./ oxySum;
else
    ORF = NaN;
end

cpx.XSi = mol.SiO2 .* ORF;
cpx.XTi = mol.TiO2 .* ORF;
cpx.XAl = 2 .* mol.Al2O3 .* ORF;
cpx.XFe = mol.FeO .* ORF;
cpx.XMn = mol.MnO .* ORF;
cpx.XMg = mol.MgO .* ORF;
cpx.XCa = mol.CaO .* ORF;
cpx.XNa = 2 .* mol.Na2O .* ORF;
cpx.XK = 2 .* mol.K2O .* ORF;
cpx.XCr = 2 .* mol.Cr2O3 .* ORF;

cpx.cationSum = ...
    cpx.XSi + cpx.XTi + cpx.XAl + cpx.XFe + cpx.XMn ...
    + cpx.XMg + cpx.XCa + cpx.XNa + cpx.XK + cpx.XCr;

% Component allocation follows Putirka (2008), Table 3. NaN is preserved;
% finite negative intermediate components are truncated to zero as in the
% original implementation's normative allocation.
cpx.XAlIV = truncateNegativePreserveNaN(2 - cpx.XSi);
cpx.XAlVI = truncateNegativePreserveNaN(cpx.XAl - cpx.XAlIV);
cpx.XFe3 = truncateNegativePreserveNaN( ...
    cpx.XNa + cpx.XAlIV - cpx.XAlVI - 2 .* cpx.XTi - cpx.XCr);
cpx.XJd = truncateNegativePreserveNaN(minPreserveNaN(cpx.XAlVI, cpx.XNa));
cpx.XCaTs = truncateNegativePreserveNaN(cpx.XAlVI - cpx.XJd);

if isnan(cpx.XAlIV) || isnan(cpx.XCaTs)
    cpx.XCaTi = NaN;
elseif cpx.XAlIV > cpx.XCaTs
    cpx.XCaTi = truncateNegativePreserveNaN((cpx.XAlIV - cpx.XCaTs) ./ 2);
else
    cpx.XCaTi = 0;
end

cpx.XCrCaTs = truncateNegativePreserveNaN(cpx.XCr ./ 2);
cpx.XDiHd = truncateNegativePreserveNaN( ...
    cpx.XCa - cpx.XCaTi - cpx.XCaTs - cpx.XCrCaTs);
cpx.XEnFs = truncateNegativePreserveNaN( ...
    (cpx.XFe + cpx.XMg - cpx.XDiHd) ./ 2);

end

function liq = prepareLiquidRow(data_liquid, MWinfo)
% prepareLiquidRow
% Calculate anhydrous liquid cation fractions. F and Cl are read only for
% output traceability and are excluded from cationTotal_liq. H2O is also
% excluded from the denominator and is used separately in wt%.

[liq.SiO2, ~, ~] = getOxideValue( ...
    data_liquid, {'SiO2'}, 0, false, 'Liquid');
[liq.TiO2, ~, ~] = getOxideValue( ...
    data_liquid, {'TiO2'}, 0, false, 'Liquid');
[liq.Al2O3, ~, ~] = getOxideValue( ...
    data_liquid, {'Al2O3'}, 0, false, 'Liquid');
[liq.FeO, ~, ~] = getOxideValue( ...
    data_liquid, {'FeO', 'FeOt'}, 0, false, 'Liquid');
[liq.MnO, ~, ~] = getOxideValue( ...
    data_liquid, {'MnO'}, 0, false, 'Liquid');
[liq.MgO, ~, ~] = getOxideValue( ...
    data_liquid, {'MgO'}, 0, false, 'Liquid');
[liq.CaO, ~, ~] = getOxideValue( ...
    data_liquid, {'CaO'}, 0, false, 'Liquid');
[liq.Na2O, ~, ~] = getOxideValue( ...
    data_liquid, {'Na2O'}, 0, false, 'Liquid');
[liq.K2O, ~, ~] = getOxideValue( ...
    data_liquid, {'K2O'}, 0, false, 'Liquid');
[liq.V2O3, ~, ~] = getOxideValue( ...
    data_liquid, {'V2O3'}, 0, false, 'Liquid');
[liq.Cr2O3, ~, ~] = getOxideValue( ...
    data_liquid, {'Cr2O3'}, 0, false, 'Liquid');
[liq.NiO, ~, ~] = getOxideValue( ...
    data_liquid, {'NiO'}, 0, false, 'Liquid');
[liq.P2O5, ~, ~] = getOxideValue( ...
    data_liquid, {'P2O5'}, 0, false, 'Liquid');
[liq.SO3, ~, ~] = getOxideValue( ...
    data_liquid, {'SO3'}, 0, false, 'Liquid');
[liq.Fe2O3, ~, ~] = getOxideValue( ...
    data_liquid, {'Fe2O3'}, 0, false, 'Liquid');
[liq.H2O, ~, ~] = getOxideValue( ...
    data_liquid, {'H2O', 'H2Ot'}, 0, false, 'Liquid');

% F and Cl are intentionally not used below.
[liq.F, ~, ~] = getOxideValue( ...
    data_liquid, {'F'}, 0, false, 'Liquid');
[liq.Cl, ~, ~] = getOxideValue( ...
    data_liquid, {'Cl'}, 0, false, 'Liquid');

n = struct();
n.SiO2 = liq.SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
n.TiO2 = liq.TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
n.Al2O3 = liq.Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n.FeO = liq.FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
n.MnO = liq.MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
n.MgO = liq.MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
n.CaO = liq.CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
n.Na2O = liq.Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
n.K2O = liq.K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
n.V2O3 = liq.V2O3 .* MWinfo.Cat.V2O3 ./ MWinfo.MW.V2O3;
n.Cr2O3 = liq.Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
n.NiO = liq.NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;
n.P2O5 = liq.P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;
n.SO3 = liq.SO3 .* MWinfo.Cat.SO3 ./ MWinfo.MW.SO3;
n.Fe2O3 = liq.Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

% F and Cl are explicitly excluded from this sum.
liq.cationTotal = ...
    n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + n.MgO ...
    + n.CaO + n.Na2O + n.K2O + n.V2O3 + n.Cr2O3 + n.NiO ...
    + n.P2O5 + n.SO3 + n.Fe2O3;

if isfinite(liq.cationTotal) && liq.cationTotal > 0
    liq.XSiO2 = n.SiO2 ./ liq.cationTotal;
    liq.XTiO2 = n.TiO2 ./ liq.cationTotal;
    liq.XAlO1_5 = n.Al2O3 ./ liq.cationTotal;
    liq.XFeO = n.FeO ./ liq.cationTotal;
    liq.XMnO = n.MnO ./ liq.cationTotal;
    liq.XMgO = n.MgO ./ liq.cationTotal;
    liq.XCaO = n.CaO ./ liq.cationTotal;
    liq.XNaO0_5 = n.Na2O ./ liq.cationTotal;
    liq.XKO0_5 = n.K2O ./ liq.cationTotal;
else
    liq.XSiO2 = NaN;
    liq.XTiO2 = NaN;
    liq.XAlO1_5 = NaN;
    liq.XFeO = NaN;
    liq.XMnO = NaN;
    liq.XMgO = NaN;
    liq.XCaO = NaN;
    liq.XNaO0_5 = NaN;
    liq.XKO0_5 = NaN;
end

liq.XFm = liq.XFeO + liq.XMnO + liq.XMgO;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify Equation (33) terms that prevent a finite positive temperature.
% Only the first row is required for pressure-independent composition terms;
% invT_Eq33 is checked over all pressure rows.

maxTerms = 10;
termBuffer = strings(maxTerms, 1);
nTerms = 0;

if ~isfinite(row.XJd_cpx(1)) || row.XJd_cpx(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XJd_cpx";
end
if ~isfinite(row.XDiHd_cpx(1)) || row.XDiHd_cpx(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XDiHd_cpx";
end
if ~isfinite(row.XEnFs_cpx(1)) || row.XEnFs_cpx(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XEnFs_cpx";
end
if ~isfinite(row.XCaO_liq(1)) || row.XCaO_liq(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XCaO_liq";
end
if ~isfinite(row.XNaO0_5_liq(1)) || row.XNaO0_5_liq(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XNaO0_5_liq";
end
if ~isfinite(row.XAlO1_5_liq(1)) || row.XAlO1_5_liq(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XAlO1_5_liq";
end
if ~isfinite(row.XTiO2_liq(1)) || row.XTiO2_liq(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XTiO2_liq";
end
if ~isfinite(row.XFm_liq(1)) || row.XFm_liq(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XFm_liq";
end
if ~isfinite(row.MgNum_liq(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "MgNum_liq";
end
if any(~isfinite(row.invT_Eq33) | row.invT_Eq33 <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "10^4/T denominator";
end

invalidTerms = termBuffer(1:nTerms);

end

function value = safeLogPositive(value)
% safeLogPositive
% Return ln(value) only when value is finite and strictly positive.

if isfinite(value) && value > 0
    value = log(value);
else
    value = NaN;
end

end

function value = truncateNegativePreserveNaN(value)
% truncateNegativePreserveNaN
% Preserve NaN; truncate only finite negative normative components to zero.

if isnan(value)
    return
end
if isfinite(value) && value < 0
    value = 0;
end

end

function value = minPreserveNaN(a, b)
% minPreserveNaN
% Return NaN if either argument is NaN; otherwise return min(a,b).

if isnan(a) || isnan(b)
    value = NaN;
else
    value = min(a, b);
end

end

function validateRequiredColumns(tbl, aliasGroups, phaseLabel)
% validateRequiredColumns
% Ensure that each required oxide or alias group exists in the table.

for i = 1:numel(aliasGroups)
    aliases = aliasGroups{i};
    name = findOxideColumn(tbl.Properties.VariableNames, aliases);
    if isempty(name)
        error('%s table must contain one of: %s', ...
            phaseLabel, strjoin(aliases, ', '));
    end
end

end

function [value, columnName, columnFound] = getOxideValue( ...
        data_table, aliases, defaultValue, required, phaseLabel)
% getOxideValue
% Retrieve one scalar oxide value using canonicalized column-name matching.
% An absent optional column receives defaultValue. An existing NaN remains
% NaN. Numeric strings are converted; empty or non-numeric stored values
% become NaN rather than silently becoming zero.

columnName = findOxideColumn(data_table.Properties.VariableNames, aliases);
columnFound = ~isempty(columnName);

if ~columnFound
    if required
        error('%s table must contain one of: %s', ...
            phaseLabel, strjoin(aliases, ', '));
    end
    value = defaultValue;
    return
end

raw = data_table.(columnName);
value = scalarToDoublePreserveNaN(raw);

end

function columnName = findOxideColumn(variableNames, aliases)
% findOxideColumn
% Match oxide names after removing spaces, underscores, and hyphens. Both
% bare oxide names and names ending in "value" are recognized.

canonicalVariables = cell(size(variableNames));
for i = 1:numel(variableNames)
    canonicalVariables{i} = canonicalizeName(variableNames{i});
end

columnName = '';
for i = 1:numel(aliases)
    canonicalAlias = canonicalizeName(aliases{i});
    targets = {[canonicalAlias 'value'], canonicalAlias};

    for j = 1:numel(targets)
        idx = find(strcmp(canonicalVariables, targets{j}), 1, 'first');
        if ~isempty(idx)
            columnName = variableNames{idx};
            return
        end
    end
end

end

function name = canonicalizeName(name)
% canonicalizeName
% Convert a table variable name to the canonical form used for matching.

name = lower(char(string(name)));
name = strrep(name, ' ', '');
name = strrep(name, '_', '');
name = strrep(name, '-', '');

end

function value = scalarToDoublePreserveNaN(raw)
% scalarToDoublePreserveNaN
% Convert a selected one-row table value to a numeric scalar while retaining
% missing or explicitly NaN values as NaN.

if isempty(raw)
    value = NaN;
    return
end

if isnumeric(raw) || islogical(raw)
    value = double(raw(1));
    return
end

if isstring(raw)
    if ismissing(raw(1))
        value = NaN;
    else
        value = str2double(raw(1));
    end
    return
end

if ischar(raw)
    value = str2double(string(raw));
    return
end

if iscell(raw)
    if isempty(raw{1})
        value = NaN;
        return
    end
    value = scalarToDoublePreserveNaN(raw{1});
    return
end

if iscategorical(raw)
    value = str2double(string(raw(1)));
    return
end

value = NaN;

end

function row = attachLiquidIDs(row, data_liquid)
% attachLiquidIDs
% Attach common liquid identifiers and repeat them for every pressure row.

nRows = height(row);
variableNames = data_liquid.Properties.VariableNames;

if ismember('Index', variableNames)
    value = data_liquid.('Index');
    if isnumeric(value) || islogical(value)
        row.liq_Index = repmat(double(value(1)), nRows, 1);
    else
        row.liq_Index = repmat(string(value(1)), nRows, 1);
    end
end
if ismember('Experiment', variableNames)
    value = data_liquid.('Experiment');
    row.liq_Experiment = repmat(string(value(1)), nRows, 1);
end
if ismember('Citation', variableNames)
    value = data_liquid.('Citation');
    row.liq_Citation = repmat(string(value(1)), nRows, 1);
end

end

function printTemperatureSummary(label, temperatureValues)
% printTemperatureSummary
% Display either one temperature or the first-to-last range.

if numel(temperatureValues) == 1
    disp([label ': ' num2str(temperatureValues) ' degreeC']);
else
    disp([label ': ' num2str(temperatureValues(1)) ' to ' ...
        num2str(temperatureValues(end)) ' degreeC']);
end

end

function printTemperatureRangeWarning( ...
        temperatureValues, minimumTemperature, maximumTemperature, ...
        selectedCode_cpx, idxLiq)
% printTemperatureRangeWarning
% Warn when finite temperatures lie outside the approximate Figure-9a
% comparison envelope. Results are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);
    fprintf(2, ...
        ['WARNING: Calculated Equation (33) temperature is outside the ' ...
         'approximate %.4g-%.4g degreeC comparison envelope shown in ' ...
         'Putirka (2008), Figure 9a. This interval is a graphical ' ...
         'validation envelope, not a strict universal calibration boundary. ' ...
         '%d of %d finite temperature point(s) are outside the envelope; ' ...
         'calculated finite range = %.4g-%.4g degreeC for %s and Liquid ' ...
         'row %d.\n'], ...
        minimumTemperature, ...
        maximumTemperature, ...
        sum(outsideMask), ...
        sum(finiteMask), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(selectedCode_cpx), ...
        idxLiq);
end

end

function printNonFiniteTemperatureWarning( ...
        temperatureValues, selectedCode_cpx, idxLiq)
% printNonFiniteTemperatureWarning
% Report NaN and Inf temperatures without modifying or deleting them.

invalidMask = ~isfinite(temperatureValues);
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Equation (33) temperature values were ' ...
         'calculated for %s and Liquid row %d (%d of %d points; NaN: %d, ' ...
         'Inf: %d).\n' ...
         '         These values remain in the output table, and the ' ...
         'calculation has not been stopped.\n'], ...
        char(selectedCode_cpx), ...
        idxLiq, ...
        sum(invalidMask), ...
        numel(temperatureValues), ...
        sum(isnan(temperatureValues)), ...
        sum(isinf(temperatureValues)));
end

end
