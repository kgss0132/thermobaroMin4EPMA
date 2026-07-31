function results = Putirka2008CpxSat(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Pyroxene/Putirka2008CpxSat.m
% Tested with MATLAB R2024b
%
% Clinopyroxene saturation thermometer
% Putirka, K.D. (2008), Equation (34)
% Reviews in Mineralogy and Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function calculates the temperature at which a selected silicate
% liquid becomes saturated with Clinopyroxene using Putirka (2008),
% Equation (34). The liquid dataset is loaded by liquid.readLiquidExcel().
%
% Equation (34) uses only liquid composition and pressure; no measured
% Clinopyroxene composition is required. rawdata_struct is retained in the
% function interface for compatibility with the common thermometer launchers.
%
% The function accepts either a scalar pressure or a pressure vector. For
% each calculation, one output row is returned for every pressure value
% supplied in P_kbar. It is therefore compatible with both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka (2008) presents the liquid-only Clinopyroxene saturation
% thermometer as Equation (34) on p. 94. Figure 9b on p. 92 evaluates the
% equation using experiments conducted at pressures below 70 kbar.
%
% Figure 9b reports the following validation statistics:
%   All data, P < 70 kbar : SEE = +/-45 degreeC, R^2 = 0.93, n = 1186
%   Anhydrous experiments : SEE = +/-45 degreeC, R^2 = 0.91, n = 861
%   Hydrous experiments   : SEE = +/-45 degreeC, R^2 = 0.85, n = 325
%
% Putirka (2008) does not list one strict minimum-maximum temperature box for
% Equation (34). Figure 9b displays the model comparison over approximately
% 600-2400 degreeC. This graphical interval is used only as a non-stopping
% screening range; it is not a formal universal calibration boundary.
%
% Important application cautions:
%   1) Equation (34) returns a CLINOPYROXENE SATURATION TEMPERATURE for the
%      selected liquid at the specified pressure. It does not by itself show
%      that an observed Clinopyroxene crystal was in equilibrium with that
%      liquid.
%
%   2) Putirka (2008) recommends comparing the liquid-only saturation
%      temperature from Equation (34) with the Cpx-Liquid temperature from
%      Equation (33) as an additional equilibrium and saturation check
%      (p. 94).
%
%   3) The selected composition must represent a plausible silicate liquid.
%      Whole-rock compositions, glasses, melt inclusions, and reconstructed
%      liquids are not automatically interchangeable. Treatment error caused
%      by choosing an inappropriate liquid composition can exceed the
%      experimental model error.
%
%   4) Liquid cation fractions must be calculated in the same way as in the
%      original calibration: on an anhydrous basis, without first
%      renormalizing oxide wt% to 100%. H2O is entered separately in wt%
%      (component-calculation discussion on pp. 68-71).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure does not satisfy the published P < 70 kbar condition,
%   2) a finite calculated temperature lies outside the approximate
%      600-2400 degreeC Figure-9b comparison envelope,
%   3) an explicitly stored calculation input is NaN,
%   4) a logarithm or final 10^4/T denominator required by Equation (34)
%      is invalid, or
%   5) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
%
% rawdata_struct:
%   A struct is required for compatibility with the common launcher
%   interface, but Equation (34) does not use mineral data from this struct.
%
% P_kbar:
%   Finite, non-negative numeric scalar or vector.
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
% not oxide cation components in the Equation (34) anhydrous normalization
% used here. F and Cl are also excluded from the NaN-input warning and from
% negative-input validation because they do not enter the calculation.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Putirka (2008), Equation (34), p. 94:
%
%   10^4 / T(K) =
%       6.39
%     + 0.076 * H2O_liq
%     - 5.55 * XCaO_liq * XSiO2_liq
%     - 0.386 * ln(XMgO_liq)
%     - 0.046 * P_kbar
%     + 2.2e-4 * P_kbar^2
%
% H2O_liq is wt% H2O and is excluded from cationTotal_liq.
% Liquid X terms are anhydrous cation fractions.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008CpxSat(rawdata_struct, P_kbar)
%   results = Putirka2008CpxSat(rawdata_struct, P_kbar, ...
%       'LiquidRow', n)
%
% Inputs:
%   rawdata_struct : struct retained for launcher-interface compatibility
%   P_kbar         : pressure in kbar; finite, non-negative numeric scalar
%                    or vector
%
% Name-value option:
%   LiquidRow      : positive integer row number in the loaded liquid table.
%                    If empty, row 1 is used.
%
% Output:
%   results : table containing one row per pressure value for every completed
%             calculation. T34_K and T34_C retain the equation-specific
%             names. T_K, T_degreeC, and T_deg are standardized aliases for
%             launcher and plotting compatibility.
%

%% Input validation
% Accept scalar or vector pressure input from both fixed-P and range-P
% launchers.
if nargin < 2
    error('Putirka2008CpxSat requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% Options
ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve liquid dataset
disp('=== Step 1: Preparing liquid dataset ===');

MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();

if isempty(liqAll) || ~istable(liqAll)
    error('Selected liquid dataset is empty or is not a table.');
end

disp('=== Preparing liquid dataset has been finished ===');

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

% Equation (34) validation uses experiments at P < 70 kbar.
calibrationP_max_kbar = 70;
pressureOutsideCalibration = P_kbar >= calibrationP_max_kbar;
pressureWarningIssued = false;

% Approximate graphical comparison envelope shown in Figure 9b.
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

%% 3) Repeated calculation loop
disp('=== Step 3: Calculating the Clinopyroxene saturation temperature ===');

while true
    % Identify explicitly stored NaN values in active calculation inputs.
    % F and Cl are intentionally absent from this check.
    nanInputNames = findNaNInputs(selectedData_liq);

    % Reject negative finite values and Inf. NaN and zero are allowed so
    % that NaN can propagate and zero-dependent invalid logs can be reported
    % without silently changing the input.
    validateNonNegativeInputs(selectedData_liq);

    row = calcTemp(selectedData_liq, P_kbar, MWinfo);

    nRows = height(row);
    row.dataRow_liq = repmat(idxLiq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataRow_liq'}, 'Before', 1);

    % Store one completed result block. The full output table is not resized
    % on every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    printTemperatureSummary('Eq. (34)', row.T34_C);

    % Pressure is common to all repeated calculations, so warn only once.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental calibration ' ...
             'condition used for Putirka (2008) Equation (34): P < 70 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; input range ' ...
             '= %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn for finite temperatures outside the Figure-9b display envelope.
    printTemperatureRangeWarning( ...
        row.T34_C, screeningT_min_degC, screeningT_max_degC, idxLiq);

    % Report exact names of explicitly NaN inputs. Values remain NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for Liquid ' ...
             'row %d: %s.\n' ...
             '         NaN values were retained and propagated; the ' ...
             'calculated temperature may be NaN.\n'], ...
            idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Give a specific warning for invalid logarithm or denominator terms.
    % The result is retained as NaN rather than stopping.
    invalidEquationTerms = findInvalidEquationTerms(row);
    if ~isempty(invalidEquationTerms)
        fprintf(2, ...
            ['WARNING: Invalid Equation (34) term(s) were found for Liquid ' ...
             'row %d: %s.\n' ...
             '         XMgO_liq and the final 10^4/T denominator must be ' ...
             'finite and > 0. Affected temperatures were retained as NaN.\n'], ...
            idxLiq, ...
            char(strjoin(invalidEquationTerms, ', ')));
    end

    printNonFiniteTemperatureWarning(row.T34_C, idxLiq);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another calculation (same Liquid row and pressure input)?', ...
        'Putirka2008CpxSat', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after the calculation loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'liquid', metaLiq, ...
    'primaryTemperatureEquation', 'Putirka2008 Eq. (34)', ...
    'liquidRow', idxLiq, ...
    'normalizationNote', ...
    'F and Cl are excluded from cationTotal_liq and from NaN-input warnings.');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_liquid)
% findNaNInputs
% Return names of explicitly stored calculation inputs whose selected value
% is NaN. Missing optional columns are not listed because they are assigned
% zero only when the column itself is absent. F and Cl are excluded.

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

nanNamesBuffer = strings(numel(liqAliases), 1);
nNames = 0;

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

function validateNonNegativeInputs(data_liquid)
% validateNonNegativeInputs
% Stop when an explicitly stored calculation input is negative or infinite.
% Zero is allowed. NaN is deliberately allowed and propagated. F and Cl are
% excluded because they do not enter cationTotal_liq or Equation (34).

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

invalidNamesBuffer = strings(numel(liqAliases), 1);
nInvalid = 0;

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
    error(['Putirka2008CpxSat: calculation inputs must not be negative ' ...
           'or infinite. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_liquid, P_kbar, MWinfo)
% calcTemp
% Calculate Putirka (2008) Equation (34) for one selected liquid and a
% scalar or vector of pressures. Output rows correspond one-to-one with
% input pressure values. NaN inputs are retained and propagated.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

liq = prepareLiquidRow(data_liquid, MWinfo);

% Equation (34) contains ln(XMgO_liq), which requires a positive finite
% argument. Invalid values become NaN to prevent complex or infinite output.
lnXMgO_liq = safeLogPositive(liq.XMgO);

% Putirka (2008), Equation (34), p. 94. The quadratic pressure coefficient
% is 2.2e-4 when P is entered in kbar.
invT = ...
    6.39 ...
    + 0.076 .* liq.H2O ...
    - 5.55 .* liq.XCaO .* liq.XSiO2 ...
    - 0.386 .* lnXMgO_liq ...
    - 0.046 .* P_kbar ...
    + 2.2e-4 .* (P_kbar .^ 2);

% A finite positive value of 10^4/T is required for a physical Kelvin
% temperature. Invalid values are represented by NaN and reported later.
T34_K = nan(nP, 1);
validInvT = isfinite(invT) & invT > 0;
T34_K(validInvT) = 10000 ./ invT(validInvT);
T34_C = T34_K - 273.15;

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PrimaryEquation = repmat("Putirka2008_Eq34", nP, 1);

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

% Equation terms and final temperature.
row.lnXMgO_liq = repmat(lnXMgO_liq, nP, 1);
row.invT_Eq34 = invT;
row.T34_K = T34_K;
row.T34_C = T34_C;

% Standardized launcher and plotting aliases.
row.T_K = T34_K;
row.T_degreeC = T34_C;
row.T_deg = T34_C;

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

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify Equation (34) terms that prevent a finite positive temperature.
% Only the first row is required for pressure-independent composition terms;
% invT_Eq34 is checked over all pressure rows.

maxTerms = 3;
termBuffer = strings(maxTerms, 1);
nTerms = 0;

if ~isfinite(row.cationTotal_liq(1)) || row.cationTotal_liq(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "cationTotal_liq";
end
if ~isfinite(row.XMgO_liq(1)) || row.XMgO_liq(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XMgO_liq";
end
if any(~isfinite(row.invT_Eq34) | row.invT_Eq34 <= 0)
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
        temperatureValues, minimumTemperature, maximumTemperature, idxLiq)
% printTemperatureRangeWarning
% Warn when finite temperatures lie outside the approximate Figure-9b
% comparison envelope. Results are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);
    fprintf(2, ...
        ['WARNING: Calculated Equation (34) temperature is outside the ' ...
         'approximate %.4g-%.4g degreeC comparison envelope shown in ' ...
         'Putirka (2008), Figure 9b. This interval is a graphical ' ...
         'validation envelope, not a strict universal calibration boundary. ' ...
         '%d of %d finite temperature point(s) are outside the envelope; ' ...
         'calculated finite range = %.4g-%.4g degreeC for Liquid row %d.\n'], ...
        minimumTemperature, ...
        maximumTemperature, ...
        sum(outsideMask), ...
        sum(finiteMask), ...
        min(finiteValues), ...
        max(finiteValues), ...
        idxLiq);
end

end

function printNonFiniteTemperatureWarning(temperatureValues, idxLiq)
% printNonFiniteTemperatureWarning
% Report NaN and Inf temperatures without modifying or deleting them.

invalidMask = ~isfinite(temperatureValues);
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Equation (34) temperature values were ' ...
         'calculated for Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the ' ...
         'calculation has not been stopped.\n'], ...
        idxLiq, ...
        sum(invalidMask), ...
        numel(temperatureValues), ...
        sum(isnan(temperatureValues)), ...
        sum(isinf(temperatureValues)));
end

end
