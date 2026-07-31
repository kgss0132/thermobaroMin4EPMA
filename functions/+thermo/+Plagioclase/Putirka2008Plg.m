function results = Putirka2008Plg(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Feldspar/Putirka2008Plg.m
% Tested with MATLAB R2024b
%
% Plagioclase-Liquid thermometers and plagioclase-saturation thermometer
% Putirka, K.D. (2008)
% Reviews in Mineralogy and Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Plagioclase analysis from
% rawdata_struct.Plagioclase and combines it with one row from a liquid
% dataset loaded by liquid.readLiquidExcel().
%
% The following equations from Putirka (2008) are calculated:
%   Eq. (23)  : Plagioclase-Liquid equilibrium thermometer
%   Eq. (24a) : Globally recalibrated Plagioclase-Liquid thermometer
%   Eq. (26)  : Plagioclase-saturation thermometer using liquid only
%
% Eq. (24c) is not included because it is an ALKALI-FELDSPAR saturation
% thermometer rather than a Plagioclase thermometer (pp. 79-82).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Plagioclase-Liquid pair, one output row is returned for each
% pressure value. It can therefore be called by both startThermoCalc_fixedP
% and startThermoCalc_rangeP.
%
% Eq. (24a), the 2008 global Plagioclase-Liquid calibration, is assigned to
% the standardized output variables T_K, T_degreeC, and T_deg. Individual
% outputs for Eq. (23), Eq. (24a), and Eq. (26) are retained separately.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka (2008) presents the Plagioclase-Liquid equations on p. 79 and
% compares their performance in Figure 5 on pp. 80-81:
%
%   Eq. (23)  : SEE approximately +/-43 degreeC; n = 1186
%   Eq. (24a) : SEE approximately +/-36 degreeC; n = 1190
%   Eq. (26)  : SEE approximately +/-43 degreeC; n = 1203
%
% Eq. (23) is reproduced from Putirka (2005). Its underlying experimental
% database spans approximately:
%
%   Temperature : 850-1430 degreeC
%   Pressure    : 0.001-27 kbar
%   Liquid SiO2 : 42-73 wt%
%   Conditions  : anhydrous and hydrous, Plagioclase-saturated experiments
%
% Putirka (2008) does not provide a separate strict minimum-maximum
% calibration box for Eq. (24a) or Eq. (26). Figure 5 on p. 80 displays
% validation data over approximately 650-1600 degreeC. This implementation
% therefore uses:
%
%   Eq. (23) warning range       : 850-1430 degreeC
%   Eq. (24a)/(26) warning range : 650-1600 degreeC
%   Practical pressure range     : 0.001-27 kbar
%
% The 650-1600 degreeC limits are a practical Figure-5 validation envelope,
% not strict universal calibration boundaries.
%
% Important application cautions:
%   1) Use a texturally and compositionally plausible equilibrium pair of
%      Plagioclase and melt. Pairing a crystal with a non-equilibrium liquid
%      is identified as a major source of treatment error on pp. 104-105.
%
%   2) Compare Eq. (23) or Eq. (24a) with the liquid-only Plagioclase
%      saturation temperature from Eq. (26). Agreement within model error
%      supports an approach to equilibrium (Figure 5 caption, pp. 80-81;
%      Toba application, pp. 100-104).
%
%   3) H2O_liq is entered independently in wt% and is excluded from the
%      anhydrous liquid cation-fraction calculation (pp. 68-71 and p. 79).
%      An absent H2O/H2Ot column is treated as 0 wt%, but an explicitly
%      stored NaN is retained and propagated.
%
%   4) Liquid and mineral components must be calculated in exactly the same
%      way as in the original calibration. Putirka (2008) emphasizes this
%      requirement on p. 68. Liquid cation fractions are calculated on an
%      anhydrous basis without first renormalizing oxide wt% to 100%
%      (pp. 68-71).
%
%   5) Na loss from open-furnace experiments or during hydrous-glass
%      microanalysis can bias Plagioclase-Liquid calculations. The problem
%      is discussed for feldspar models on p. 82.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 0.001-27 kbar,
%   2) a finite Eq. (23) temperature is outside 850-1430 degreeC,
%   3) a finite Eq. (24a) or Eq. (26) temperature is outside the practical
%      Figure-5 validation envelope of 650-1600 degreeC,
%   4) finite liquid SiO2 is outside 42-73 wt%,
%   5) a calculation input contains NaN, or
%   6) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Plagioclase : table
%
% The FIRST column of the Plagioclase table is treated as an identifier
% ("data code") displayed in the selection dialog.
%
% Required Plagioclase oxide columns:
%   SiO2, Al2O3, CaO, Na2O
%
% Optional Plagioclase oxide columns:
%   TiO2, FeO or FeOt, MnO, MgO, K2O, Cr2O3
%
% The liquid table is loaded by liquid.readLiquidExcel().
%
% Required liquid oxide columns:
%   SiO2, Al2O3, CaO, Na2O
%
% Optional liquid oxide columns:
%   TiO2, FeO or FeOt, MnO, MgO, K2O, Cr2O3, H2O or H2Ot
%
% Missing optional columns are assigned 0. An explicitly stored NaN is not
% replaced by 0. All finite calculation inputs must be greater than or equal
% to zero. Negative finite values and Inf stop the calculation. NaN values
% are retained, propagated, and reported by non-stopping fprintf warnings.
%
% -------------------------------------------------------------------------
% COMPONENT CALCULATIONS
%
% Liquid cation fractions are calculated on an anhydrous basis using:
%   SiO2, TiO2, AlO1.5, FeO, MnO, MgO, CaO,
%   NaO0.5, KO0.5, and CrO1.5
%
% H2O is excluded from the cation-fraction denominator and is entered
% separately in wt%. Fe2O3, P2O5, SO3, F, Cl, V2O3, and NiO are not included
% in the denominator used by these equations.
%
% Plagioclase components are calculated as:
%   XAn_pl = Ca / (Ca + Na + K)
%   XAb_pl = Na / (Ca + Na + K)
%   XOr_pl = K  / (Ca + Na + K)
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATIONS
%
% Eq. (23), p. 79:
%   T(K) = 10^4 / D23
%
%   D23 = 6.12
%       + 0.257 * ln[XAn_pl / (XCaO_liq*XAlO1.5_liq^2*XSiO2_liq^2)]
%       - 3.166 * XCaO_liq
%       - 3.137 * XAlO1.5_liq/(XAlO1.5_liq + XSiO2_liq)
%       + 1.216 * XAb_pl^2
%       - 2.475e-2 * P_kbar
%       + 0.2166 * H2O_liq
%
% Eq. (24a), p. 79:
%   T(K) = 10^4 / D24a
%
%   D24a = 6.4706
%        + 0.3128 * ln[XAn_pl/(XCaO_liq*XAlO1.5_liq^2*XSiO2_liq^2)]
%        - 8.103 * XSiO2_liq
%        + 4.872 * XKO0.5_liq
%        + 1.5346 * XAb_pl^2
%        + 8.661 * XSiO2_liq^2
%        - 3.341e-2 * P_kbar
%        + 0.18047 * H2O_liq
%
% Eq. (26), p. 83:
%   T(K) = 10^4 / D26
%
%   D26 = 10.86
%       - 9.7654 * XSiO2_liq
%       + 4.241 * XCaO_liq
%       - 55.56 * XCaO_liq*XAlO1.5_liq
%       + 37.50 * XKO0.5_liq*XAlO1.5_liq
%       + 11.206 * XSiO2_liq^3
%       - 3.151e-2 * P_kbar
%       + 0.1709 * H2O_liq
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008Plg(rawdata_struct, P_kbar)
%   results = Putirka2008Plg(rawdata_struct, P_kbar, 'LiquidRow', n)
%
% Inputs:
%   rawdata_struct : struct containing a Plagioclase table
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Name-value option:
%   LiquidRow      : positive integer row number in the loaded liquid table.
%                    If empty, row 1 is used.
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Plagioclase-Liquid pair. Standard T_K/T_degreeC/T_deg values
%             correspond to Eq. (24a).
%

%% Input validation
% Basic argument checks allow scalar and vector pressure inputs from both
% fixed-P and range-P launchers.
if nargin < 2
    error('Putirka2008Plg requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

if ~isfield(rawdata_struct, 'Plagioclase') || ...
        ~istable(rawdata_struct.Plagioclase)
    error('rawdata_struct must contain table: rawdata_struct.Plagioclase');
end

%% Options
ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve Plagioclase and liquid datasets
disp('=== Step 1: Preparing Plagioclase and liquid datasets ===');

MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();

if isempty(liqAll) || ~istable(liqAll)
    error('Selected liquid dataset is empty or is not a table.');
end

dataset_pl = rawdata_struct.Plagioclase;

% Required columns are checked before the interactive loop. Their values may
% still be NaN because NaN must be retained and propagated.
validateRequiredColumns(dataset_pl, ...
    {{'SiO2'}, {'Al2O3'}, {'CaO'}, {'Na2O'}}, 'Plagioclase');
validateRequiredColumns(liqAll, ...
    {{'SiO2'}, {'Al2O3'}, {'CaO'}, {'Na2O'}}, 'Liquid');

disp('=== Preparing Plagioclase and liquid datasets has been finished ===');

%% 2) Select the liquid row and initialize the output container
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

% Practical warning limits described in the header.
pressureMin_kbar = 0.001;
pressureMax_kbar = 27;
temperatureEq23Min_degC = 850;
temperatureEq23Max_degC = 1430;
temperatureGlobalMin_degC = 650;
temperatureGlobalMax_degC = 1600;
liquidSiO2Min_wt = 42;
liquidSiO2Max_wt = 73;

pressureOutsideRange = ...
    P_kbar < pressureMin_kbar | P_kbar > pressureMax_kbar;
pressureWarningIssued = false;

SiO2_liq_forRange = getLiqOxRequired(selectedData_liq, {'SiO2'});
liquidCompositionWarningIssued = false;

[~, ~, h2oColumnFound] = getOxideValue( ...
    selectedData_liq, {'H2O', 'H2Ot'}, 0, false, 'Liquid');
if ~h2oColumnFound
    fprintf(2, ...
        ['WARNING: No H2O or H2Ot column was found in the selected liquid ' ...
         'dataset. H2O_liq = 0 wt%% is being assumed.\n']);
end

disp(['Liquid selected: Row ' num2str(idxLiq)]);
disp('=== Preparing liquid selection and output container has been finished ===');

%% 3-4) Interactive Plagioclase selection and calculation
disp('=== Step 3: Selecting a data code from the list (Plagioclase) ===');

dataCodes_pl = dataset_pl{:, 1};
displayCodes_pl = cellstr(string(dataCodes_pl));

while true
    [selectedIdx_pl, ok] = listdlg( ...
        'PromptString', 'Please select the Plagioclase data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_pl, ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_pl)
        disp('Selection canceled');
        break;
    end

    selectedCode_pl = string(dataCodes_pl(selectedIdx_pl));
    selectedData_pl = dataset_pl(selectedIdx_pl, :);

    disp(['Plagioclase selected: ' char(selectedCode_pl)]);
    disp('=== Step 4: Calculating the temperature ===');

    % Explicit NaN values are identified without stopping the calculation.
    nanInputNames = findNaNInputs(selectedData_pl, selectedData_liq);

    % Negative finite values and Inf are prohibited. Zero is allowed. NaN is
    % deliberately permitted so that it propagates to the result.
    validateNonNegativeInputs(selectedData_pl, selectedData_liq);

    row = calcTemp(selectedData_pl, selectedData_liq, P_kbar, MWinfo);

    nRows = height(row);
    row.dataCode_pl = repmat(selectedCode_pl, nRows, 1);
    row.dataRow_liq = repmat(idxLiq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_pl', 'dataRow_liq'}, 'Before', 1);

    % Store one completed result block. This avoids repeated expansion and
    % copying of the full results table on every interactive iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo calculated temperatures.
    disp('--------------------------------------------------');
    disp('=== Temperatures were calculated: ===');
    printTemperatureSummary('Eq. (23)', row.T_Eq23_C);
    printTemperatureSummary('Eq. (24a)', row.T_Eq24a_C);
    printTemperatureSummary('Eq. (26)', row.T_Eq26_C);

    % Pressure warning is common to all selected Plagioclase rows and is
    % therefore printed only once.
    if any(pressureOutsideRange) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the practical experimental-' ...
             'database range used for Putirka Plagioclase-Liquid thermometry: ' ...
             '0.001-27 kbar. %d of %d pressure point(s) are outside the ' ...
             'range; input range = %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideRange), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % The selected liquid is common to all Plagioclase rows, so this warning
    % is printed only once.
    if ~liquidCompositionWarningIssued && isfinite(SiO2_liq_forRange) && ...
            (SiO2_liq_forRange < liquidSiO2Min_wt || ...
             SiO2_liq_forRange > liquidSiO2Max_wt)
        fprintf(2, ...
            ['WARNING: Liquid SiO2 is outside the experimental-database ' ...
             'range associated with Eq. (23): 42-73 wt%%. Selected liquid ' ...
             'SiO2 = %.4g wt%% (Liquid row %d).\n'], ...
            SiO2_liq_forRange, idxLiq);
        liquidCompositionWarningIssued = true;
    end

    printTemperatureRangeWarning( ...
        row.T_Eq23_C, temperatureEq23Min_degC, temperatureEq23Max_degC, ...
        'Eq. (23)', 'experimental-database', ...
        selectedCode_pl, idxLiq);

    printTemperatureRangeWarning( ...
        row.T_Eq24a_C, temperatureGlobalMin_degC, temperatureGlobalMax_degC, ...
        'Eq. (24a)', 'Figure-5 validation-envelope', ...
        selectedCode_pl, idxLiq);

    printTemperatureRangeWarning( ...
        row.T_Eq26_C, temperatureGlobalMin_degC, temperatureGlobalMax_degC, ...
        'Eq. (26)', 'Figure-5 validation-envelope', ...
        selectedCode_pl, idxLiq);

    % List all explicitly NaN calculation inputs. Values remain NaN and are
    % not replaced by zero.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s ' ...
             'and Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; one or more ' ...
             'calculated temperatures may be NaN.\n'], ...
            char(selectedCode_pl), ...
            idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    printNonFiniteTemperatureWarning( ...
        row.T_Eq23_C, 'Eq. (23)', selectedCode_pl, idxLiq);
    printNonFiniteTemperatureWarning( ...
        row.T_Eq24a_C, 'Eq. (24a)', selectedCode_pl, idxLiq);
    printNonFiniteTemperatureWarning( ...
        row.T_Eq26_C, 'Eq. (26)', selectedCode_pl, idxLiq);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Plagioclase selection (same Liquid dataset)?', ...
        'Putirka2008Plg', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks once after the selection loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'liquid', metaLiq, ...
    'primaryTemperatureEquation', 'Putirka2008 Eq. (24a)', ...
    'excludedEquation', ...
    'Eq. (24c) excluded because it is an alkali-feldspar saturation thermometer.');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_plagioclase, data_liquid)
% findNaNInputs
% Return names of explicitly stored calculation inputs whose selected value
% is NaN. Missing optional columns are not listed because they are assigned
% zero only when the column itself is absent.

plagAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'Cr2O3'}};
plagLabels = [ ...
    "Plagioclase.SiO2", "Plagioclase.TiO2", ...
    "Plagioclase.Al2O3", "Plagioclase.FeO/FeOt", ...
    "Plagioclase.MnO", "Plagioclase.MgO", ...
    "Plagioclase.CaO", "Plagioclase.Na2O", ...
    "Plagioclase.K2O", "Plagioclase.Cr2O3"];

liqAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'Cr2O3'}, ...
    {'H2O', 'H2Ot'}};
liqLabels = [ ...
    "Liquid.SiO2", "Liquid.TiO2", "Liquid.Al2O3", ...
    "Liquid.FeO/FeOt", "Liquid.MnO", "Liquid.MgO", ...
    "Liquid.CaO", "Liquid.Na2O", "Liquid.K2O", ...
    "Liquid.Cr2O3", "Liquid.H2O/H2Ot"];

maxNames = numel(plagAliases) + numel(liqAliases);
nanNamesBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(plagAliases)
    [value, ~, columnFound] = getOxideValue( ...
        data_plagioclase, plagAliases{i}, 0, false, 'Plagioclase');
    if columnFound && isnan(value)
        nNames = nNames + 1;
        nanNamesBuffer(nNames) = plagLabels(i);
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

function validateNonNegativeInputs(data_plagioclase, data_liquid)
% validateNonNegativeInputs
% Stop when an explicitly stored calculation input is negative or infinite.
% Zero is allowed. NaN is deliberately allowed and propagated.

plagAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'Cr2O3'}};
plagLabels = [ ...
    "Plagioclase.SiO2", "Plagioclase.TiO2", ...
    "Plagioclase.Al2O3", "Plagioclase.FeO/FeOt", ...
    "Plagioclase.MnO", "Plagioclase.MgO", ...
    "Plagioclase.CaO", "Plagioclase.Na2O", ...
    "Plagioclase.K2O", "Plagioclase.Cr2O3"];

liqAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'Cr2O3'}, ...
    {'H2O', 'H2Ot'}};
liqLabels = [ ...
    "Liquid.SiO2", "Liquid.TiO2", "Liquid.Al2O3", ...
    "Liquid.FeO/FeOt", "Liquid.MnO", "Liquid.MgO", ...
    "Liquid.CaO", "Liquid.Na2O", "Liquid.K2O", ...
    "Liquid.Cr2O3", "Liquid.H2O/H2Ot"];

maxNames = numel(plagAliases) + numel(liqAliases);
invalidNamesBuffer = strings(maxNames, 1);
nInvalid = 0;

for i = 1:numel(plagAliases)
    [value, ~, columnFound] = getOxideValue( ...
        data_plagioclase, plagAliases{i}, 0, false, 'Plagioclase');
    if columnFound && (isinf(value) || (isfinite(value) && value < 0))
        nInvalid = nInvalid + 1;
        invalidNamesBuffer(nInvalid) = plagLabels(i);
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
    error(['Putirka2008Plg: calculation inputs must not be negative or ' ...
           'infinite. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_plagioclase, data_liquid, P_kbar, MWinfo)
% calcTemp
% Calculate Putirka (2008) Eq. (23), Eq. (24a), and Eq. (26) for one
% selected Plagioclase-Liquid pair and a scalar or vector of pressures.
%
% Output rows correspond one-to-one with input pressure values.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

pl = preparePlagioclaseRow(data_plagioclase, MWinfo);
liq = prepareLiquidRow(data_liquid, MWinfo);

KAn_pl_liq = pl.XAn ./ ...
    (liq.XCaO .* (liq.XAlO1_5 .^ 2) .* (liq.XSiO2 .^ 2));
lnKAn_pl_liq = log(KAn_pl_liq);

denEq23 = ...
    6.12 ...
    + 0.257 .* lnKAn_pl_liq ...
    - 3.166 .* liq.XCaO ...
    - 3.137 .* (liq.XAlO1_5 ./ (liq.XAlO1_5 + liq.XSiO2)) ...
    + 1.216 .* (pl.XAb .^ 2) ...
    - 2.475e-2 .* P_kbar ...
    + 0.2166 .* liq.H2O;

denEq24a = ...
    6.4706 ...
    + 0.3128 .* lnKAn_pl_liq ...
    - 8.103 .* liq.XSiO2 ...
    + 4.872 .* liq.XKO0_5 ...
    + 1.5346 .* (pl.XAb .^ 2) ...
    + 8.661 .* (liq.XSiO2 .^ 2) ...
    - 3.341e-2 .* P_kbar ...
    + 0.18047 .* liq.H2O;

denEq26 = ...
    10.86 ...
    - 9.7654 .* liq.XSiO2 ...
    + 4.241 .* liq.XCaO ...
    - 55.56 .* (liq.XCaO .* liq.XAlO1_5) ...
    + 37.50 .* (liq.XKO0_5 .* liq.XAlO1_5) ...
    + 11.206 .* (liq.XSiO2 .^ 3) ...
    - 3.151e-2 .* P_kbar ...
    + 0.1709 .* liq.H2O;

% NaN and Inf are not replaced. MATLAB arithmetic propagates them naturally.
T_Eq23_K = 1e4 ./ denEq23;
T_Eq23_C = T_Eq23_K - 273.15;

T_Eq24a_K = 1e4 ./ denEq24a;
T_Eq24a_C = T_Eq24a_K - 273.15;

T_Eq26_K = 1e4 ./ denEq26;
T_Eq26_C = T_Eq26_K - 273.15;

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PrimaryEquation = repmat("Putirka2008_Eq24a", nP, 1);

% Selected Plagioclase oxides.
row.SiO2_pl = repmat(pl.SiO2, nP, 1);
row.TiO2_pl = repmat(pl.TiO2, nP, 1);
row.Al2O3_pl = repmat(pl.Al2O3, nP, 1);
row.FeO_pl = repmat(pl.FeO, nP, 1);
row.MnO_pl = repmat(pl.MnO, nP, 1);
row.MgO_pl = repmat(pl.MgO, nP, 1);
row.CaO_pl = repmat(pl.CaO, nP, 1);
row.Na2O_pl = repmat(pl.Na2O, nP, 1);
row.K2O_pl = repmat(pl.K2O, nP, 1);
row.Cr2O3_pl = repmat(pl.Cr2O3, nP, 1);

% Plagioclase cation fractions and end-member components.
row.cationTotal_pl = repmat(pl.cationTotal, nP, 1);
row.XSiO2_pl = repmat(pl.XSiO2, nP, 1);
row.XTiO2_pl = repmat(pl.XTiO2, nP, 1);
row.XAlO1_5_pl = repmat(pl.XAlO1_5, nP, 1);
row.XFeO_pl = repmat(pl.XFeO, nP, 1);
row.XMnO_pl = repmat(pl.XMnO, nP, 1);
row.XMgO_pl = repmat(pl.XMgO, nP, 1);
row.XCaO_pl = repmat(pl.XCaO, nP, 1);
row.XNaO0_5_pl = repmat(pl.XNaO0_5, nP, 1);
row.XKO0_5_pl = repmat(pl.XKO0_5, nP, 1);
row.XCrO1_5_pl = repmat(pl.XCrO1_5, nP, 1);
row.XAn_pl = repmat(pl.XAn, nP, 1);
row.XAb_pl = repmat(pl.XAb, nP, 1);
row.XOr_pl = repmat(pl.XOr, nP, 1);

% Selected liquid oxides.
row.SiO2_liq = repmat(liq.SiO2, nP, 1);
row.TiO2_liq = repmat(liq.TiO2, nP, 1);
row.Al2O3_liq = repmat(liq.Al2O3, nP, 1);
row.FeO_liq = repmat(liq.FeO, nP, 1);
row.MnO_liq = repmat(liq.MnO, nP, 1);
row.MgO_liq = repmat(liq.MgO, nP, 1);
row.CaO_liq = repmat(liq.CaO, nP, 1);
row.Na2O_liq = repmat(liq.Na2O, nP, 1);
row.K2O_liq = repmat(liq.K2O, nP, 1);
row.Cr2O3_liq = repmat(liq.Cr2O3, nP, 1);
row.H2O_liq_wt = repmat(liq.H2O, nP, 1);

% Liquid cation fractions.
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
row.XCrO1_5_liq = repmat(liq.XCrO1_5, nP, 1);

% Equilibrium and denominator terms.
row.KAn_pl_liq = repmat(KAn_pl_liq, nP, 1);
row.lnKAn_pl_liq = repmat(lnKAn_pl_liq, nP, 1);
row.den_Eq23 = denEq23;
row.den_Eq24a = denEq24a;
row.den_Eq26 = denEq26;

% Individual equation outputs.
row.T_Eq23_K = T_Eq23_K;
row.T_Eq23_C = T_Eq23_C;
row.T_Eq24a_K = T_Eq24a_K;
row.T_Eq24a_C = T_Eq24a_C;
row.T_Eq26_K = T_Eq26_K;
row.T_Eq26_C = T_Eq26_C;

% Differences useful for evaluating Plagioclase-Liquid equilibrium.
row.DeltaT_Eq23_minus_Eq26_C = T_Eq23_C - T_Eq26_C;
row.DeltaT_Eq24a_minus_Eq26_C = T_Eq24a_C - T_Eq26_C;

% Standardized launcher and plotting outputs use Eq. (24a).
row.T_K = T_Eq24a_K;
row.T_degreeC = T_Eq24a_C;
row.T_deg = T_Eq24a_C;

end

function pl = preparePlagioclaseRow(data_plagioclase, MWinfo)
% preparePlagioclaseRow
% Convert one Plagioclase oxide row to anhydrous cation fractions and
% calculate An-Ab-Or components. Explicit NaN values remain NaN.

SiO2  = getMineralOxRequired(data_plagioclase, {'SiO2'});
TiO2  = getMineralOxOptional(data_plagioclase, {'TiO2'}, 0);
Al2O3 = getMineralOxRequired(data_plagioclase, {'Al2O3'});
FeO   = getMineralOxOptional(data_plagioclase, {'FeO', 'FeOt'}, 0);
MnO   = getMineralOxOptional(data_plagioclase, {'MnO'}, 0);
MgO   = getMineralOxOptional(data_plagioclase, {'MgO'}, 0);
CaO   = getMineralOxRequired(data_plagioclase, {'CaO'});
Na2O  = getMineralOxRequired(data_plagioclase, {'Na2O'});
K2O   = getMineralOxOptional(data_plagioclase, {'K2O'}, 0);
Cr2O3 = getMineralOxOptional(data_plagioclase, {'Cr2O3'}, 0);

nSiO2  = SiO2  .* MWinfo.Cat.SiO2  ./ MWinfo.MW.SiO2;
nTiO2  = TiO2  .* MWinfo.Cat.TiO2  ./ MWinfo.MW.TiO2;
nAl2O3 = Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
nFeO   = FeO   .* MWinfo.Cat.FeO   ./ MWinfo.MW.FeO;
nMnO   = MnO   .* MWinfo.Cat.MnO   ./ MWinfo.MW.MnO;
nMgO   = MgO   .* MWinfo.Cat.MgO   ./ MWinfo.MW.MgO;
nCaO   = CaO   .* MWinfo.Cat.CaO   ./ MWinfo.MW.CaO;
nNa2O  = Na2O  .* MWinfo.Cat.Na2O  ./ MWinfo.MW.Na2O;
nK2O   = K2O   .* MWinfo.Cat.K2O   ./ MWinfo.MW.K2O;
nCr2O3 = Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;

cationTotal = nSiO2 + nTiO2 + nAl2O3 + nFeO + nMnO + ...
    nMgO + nCaO + nNa2O + nK2O + nCr2O3;

alkaliSiteTotal = nCaO + nNa2O + nK2O;

pl = struct();
pl.SiO2 = SiO2;
pl.TiO2 = TiO2;
pl.Al2O3 = Al2O3;
pl.FeO = FeO;
pl.MnO = MnO;
pl.MgO = MgO;
pl.CaO = CaO;
pl.Na2O = Na2O;
pl.K2O = K2O;
pl.Cr2O3 = Cr2O3;

pl.cationTotal = cationTotal;
pl.XSiO2 = nSiO2 ./ cationTotal;
pl.XTiO2 = nTiO2 ./ cationTotal;
pl.XAlO1_5 = nAl2O3 ./ cationTotal;
pl.XFeO = nFeO ./ cationTotal;
pl.XMnO = nMnO ./ cationTotal;
pl.XMgO = nMgO ./ cationTotal;
pl.XCaO = nCaO ./ cationTotal;
pl.XNaO0_5 = nNa2O ./ cationTotal;
pl.XKO0_5 = nK2O ./ cationTotal;
pl.XCrO1_5 = nCr2O3 ./ cationTotal;

pl.XAn = nCaO ./ alkaliSiteTotal;
pl.XAb = nNa2O ./ alkaliSiteTotal;
pl.XOr = nK2O ./ alkaliSiteTotal;

end

function liq = prepareLiquidRow(data_liquid, MWinfo)
% prepareLiquidRow
% Convert one liquid oxide row to anhydrous cation fractions. H2O is stored
% separately in wt%. Explicit NaN values remain NaN.

SiO2  = getLiqOxRequired(data_liquid, {'SiO2'});
TiO2  = getLiqOxOptional(data_liquid, {'TiO2'}, 0);
Al2O3 = getLiqOxRequired(data_liquid, {'Al2O3'});
FeO   = getLiqOxOptional(data_liquid, {'FeO', 'FeOt'}, 0);
MnO   = getLiqOxOptional(data_liquid, {'MnO'}, 0);
MgO   = getLiqOxOptional(data_liquid, {'MgO'}, 0);
CaO   = getLiqOxRequired(data_liquid, {'CaO'});
Na2O  = getLiqOxRequired(data_liquid, {'Na2O'});
K2O   = getLiqOxOptional(data_liquid, {'K2O'}, 0);
Cr2O3 = getLiqOxOptional(data_liquid, {'Cr2O3'}, 0);
H2O   = getLiqOxOptional(data_liquid, {'H2O', 'H2Ot'}, 0);

nSiO2  = SiO2  .* MWinfo.Cat.SiO2  ./ MWinfo.MW.SiO2;
nTiO2  = TiO2  .* MWinfo.Cat.TiO2  ./ MWinfo.MW.TiO2;
nAl2O3 = Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
nFeO   = FeO   .* MWinfo.Cat.FeO   ./ MWinfo.MW.FeO;
nMnO   = MnO   .* MWinfo.Cat.MnO   ./ MWinfo.MW.MnO;
nMgO   = MgO   .* MWinfo.Cat.MgO   ./ MWinfo.MW.MgO;
nCaO   = CaO   .* MWinfo.Cat.CaO   ./ MWinfo.MW.CaO;
nNa2O  = Na2O  .* MWinfo.Cat.Na2O  ./ MWinfo.MW.Na2O;
nK2O   = K2O   .* MWinfo.Cat.K2O   ./ MWinfo.MW.K2O;
nCr2O3 = Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;

cationTotal = nSiO2 + nTiO2 + nAl2O3 + nFeO + nMnO + ...
    nMgO + nCaO + nNa2O + nK2O + nCr2O3;

liq = struct();
liq.SiO2 = SiO2;
liq.TiO2 = TiO2;
liq.Al2O3 = Al2O3;
liq.FeO = FeO;
liq.MnO = MnO;
liq.MgO = MgO;
liq.CaO = CaO;
liq.Na2O = Na2O;
liq.K2O = K2O;
liq.Cr2O3 = Cr2O3;
liq.H2O = H2O;

liq.cationTotal = cationTotal;
liq.XSiO2 = nSiO2 ./ cationTotal;
liq.XTiO2 = nTiO2 ./ cationTotal;
liq.XAlO1_5 = nAl2O3 ./ cationTotal;
liq.XFeO = nFeO ./ cationTotal;
liq.XMnO = nMnO ./ cationTotal;
liq.XMgO = nMgO ./ cationTotal;
liq.XCaO = nCaO ./ cationTotal;
liq.XNaO0_5 = nNa2O ./ cationTotal;
liq.XKO0_5 = nK2O ./ cationTotal;
liq.XCrO1_5 = nCr2O3 ./ cationTotal;

end

function printTemperatureSummary(equationLabel, temperature_degC)
% printTemperatureSummary
% Display one value or the first-to-last range for a temperature vector.

if numel(temperature_degC) == 1
    disp([equationLabel ': ' num2str(temperature_degC) ' degreeC']);
else
    disp([equationLabel ': ' num2str(temperature_degC(1)) ' to ' ...
        num2str(temperature_degC(end)) ' degreeC']);
end

end

function printTemperatureRangeWarning(temperature_degC, rangeMin_degC, ...
        rangeMax_degC, equationLabel, rangeDescription, ...
        dataCode_pl, liquidRow)
% printTemperatureRangeWarning
% Print a non-stopping warning for finite temperatures outside a specified
% calibration or validation range.

finiteTemperature = isfinite(temperature_degC);
outsideRange = finiteTemperature & ...
    (temperature_degC < rangeMin_degC | temperature_degC > rangeMax_degC);

if any(outsideRange)
    finiteValues = temperature_degC(finiteTemperature);
    fprintf(2, ...
        ['WARNING: Calculated %s temperature is outside the %s range: ' ...
         '%.4g-%.4g degreeC. %d of %d finite temperature point(s) are ' ...
         'outside the range; calculated finite range = %.4g-%.4g degreeC ' ...
         'for %s and Liquid row %d.\n'], ...
        equationLabel, ...
        rangeDescription, ...
        rangeMin_degC, ...
        rangeMax_degC, ...
        sum(outsideRange), ...
        sum(finiteTemperature), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(dataCode_pl), ...
        liquidRow);
end

end

function printNonFiniteTemperatureWarning(temperature_degC, equationLabel, ...
        dataCode_pl, liquidRow)
% printNonFiniteTemperatureWarning
% Report NaN and Inf temperatures without deleting or replacing them.

invalidTemperature = ~isfinite(temperature_degC);

if any(invalidTemperature)
    fprintf(2, ...
        ['WARNING: Non-finite %s temperature values were calculated for %s ' ...
         'and Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the ' ...
         'calculation has not been stopped.\n'], ...
        equationLabel, ...
        char(dataCode_pl), ...
        liquidRow, ...
        sum(invalidTemperature), ...
        numel(temperature_degC), ...
        sum(isnan(temperature_degC)), ...
        sum(isinf(temperature_degC)));
end

end

function row = attachLiquidIDs(row, data_liquid)
% attachLiquidIDs
% Copy common liquid identifiers to every pressure row.

nRows = height(row);
variableNames = data_liquid.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    row.liq_Index = repmat(data_liquid.('Index'), nRows, 1);
end
if any(strcmp(variableNames, 'Experiment'))
    row.liq_Experiment = repmat(string(data_liquid.('Experiment')), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    row.liq_Citation = repmat(string(data_liquid.('Citation')), nRows, 1);
end

end

function validateRequiredColumns(data_table, aliasGroups, phaseName)
% validateRequiredColumns
% Confirm that every required oxide has at least one accepted column alias.

missingBuffer = strings(numel(aliasGroups), 1);
nMissing = 0;

for i = 1:numel(aliasGroups)
    aliases = aliasGroups{i};
    [~, matchedName] = findOxideColumnAliases( ...
        data_table.Properties.VariableNames, aliases);
    if isempty(matchedName)
        nMissing = nMissing + 1;
        missingBuffer(nMissing) = strjoin(string(aliases), '/');
    end
end

if nMissing > 0
    missingNames = missingBuffer(1:nMissing);
    error('%s table is missing required oxide column(s): %s.', ...
        phaseName, char(strjoin(missingNames, ', ')));
end

end

function value = getMineralOxRequired(data_table, aliases)
% getMineralOxRequired
% Read a required Plagioclase oxide. Existing NaN remains NaN.

[value, ~, ~] = getOxideValue( ...
    data_table, aliases, NaN, true, 'Plagioclase');

end

function value = getMineralOxOptional(data_table, aliases, defaultValue)
% getMineralOxOptional
% Read an optional Plagioclase oxide. The default is used only when no
% accepted column alias exists; existing NaN remains NaN.

[value, ~, ~] = getOxideValue( ...
    data_table, aliases, defaultValue, false, 'Plagioclase');

end

function value = getLiqOxRequired(data_table, aliases)
% getLiqOxRequired
% Read a required liquid oxide. Existing NaN remains NaN.

[value, ~, ~] = getOxideValue( ...
    data_table, aliases, NaN, true, 'Liquid');

end

function value = getLiqOxOptional(data_table, aliases, defaultValue)
% getLiqOxOptional
% Read an optional liquid oxide. The default is used only when no accepted
% column alias exists; existing NaN remains NaN.

[value, ~, ~] = getOxideValue( ...
    data_table, aliases, defaultValue, false, 'Liquid');

end

function [value, matchedName, columnFound] = getOxideValue( ...
        data_table, aliases, defaultValue, isRequired, phaseName)
% getOxideValue
% Resolve accepted oxide aliases and convert a one-row table value to a
% scalar double without replacing explicit NaN.

if ischar(aliases) || isstring(aliases)
    aliases = cellstr(string(aliases));
end

[~, matchedName] = findOxideColumnAliases( ...
    data_table.Properties.VariableNames, aliases);
columnFound = ~isempty(matchedName);

if ~columnFound
    if isRequired
        error('%s table must contain oxide column: %s.', ...
            phaseName, char(strjoin(string(aliases), ' or ')));
    end
    value = defaultValue;
    return
end

value = toScalarDouble(data_table.(matchedName));

end

function [columnIndex, matchedName] = findOxideColumnAliases(varNames, aliases)
% findOxideColumnAliases
% Find the first table variable matching an accepted oxide alias. Matching
% is case-insensitive and ignores spaces, underscores, and hyphens.

canonicalVarNames = strings(numel(varNames), 1);
for i = 1:numel(varNames)
    canonicalVarNames(i) = canonicalizeName(varNames{i});
end

columnIndex = [];
matchedName = '';

for i = 1:numel(aliases)
    canonicalAlias = canonicalizeName(aliases{i});
    targets = [canonicalAlias + "value", canonicalAlias];

    for j = 1:numel(targets)
        idx = find(canonicalVarNames == targets(j), 1, 'first');
        if ~isempty(idx)
            columnIndex = idx;
            matchedName = varNames{idx};
            return
        end
    end
end

end

function canonicalName = canonicalizeName(inputName)
% canonicalizeName
% Convert a table-variable or oxide name to the form used for matching.

canonicalName = lower(string(inputName));
canonicalName = replace(canonicalName, " ", "");
canonicalName = replace(canonicalName, "_", "");
canonicalName = replace(canonicalName, "-", "");

end

function value = toScalarDouble(rawValue)
% toScalarDouble
% Convert the first value from a one-row table variable to double. Missing,
% empty, or non-convertible values become NaN. NaN is never replaced here.

value = NaN;

if isempty(rawValue)
    return
end

if isnumeric(rawValue) || islogical(rawValue)
    value = double(rawValue(1));
    return
end

if iscell(rawValue)
    if isempty(rawValue{1})
        return
    end
    rawValue = rawValue{1};
end

if iscategorical(rawValue)
    rawValue = string(rawValue(1));
elseif isstring(rawValue)
    if ismissing(rawValue(1))
        return
    end
    rawValue = rawValue(1);
elseif ischar(rawValue)
    rawValue = string(rawValue);
end

if isnumeric(rawValue) || islogical(rawValue)
    value = double(rawValue(1));
    return
end

convertedValue = str2double(string(rawValue));
if ~isempty(convertedValue)
    value = convertedValue(1);
end

end
