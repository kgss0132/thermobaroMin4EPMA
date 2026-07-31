function results = Zhang2020TiLiq(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Quartz/Zhang2020TiLiq.m
% Tested with MATLAB R2024b
%
% Quartz-Liquid Ti thermobarometer / geothermometer
% Zhang, C., Li, X., Almeev, R.R., Horn, I., Behrens, H., and Holtz, F. (2020)
% Earth and Planetary Science Letters, 538, 116213
% DOI: https://doi.org/10.1016/j.epsl.2020.116213
%
% -------------------------------------------------------------------------
% OVERVIEW
%
% This function interactively selects one Quartz analysis and combines it
% with one row from an externally selected Liquid dataset. Temperature is
% calculated using the Quartz-Liquid Ti partitioning model of
% Zhang et al. (2020), Equation (8).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Quartz-Liquid pair, one output row is returned for every
% pressure supplied in P_kbar. This interface is compatible with both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% The Liquid dataset is loaded once with liquid.readLiquidExcel(). By
% default, row 1 is used. A different row can be specified with the
% optional 'LiquidRow' name-value input. The selected Liquid row is held
% constant while the user selects additional Quartz analyses.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Zhang et al. (2020) calibrated Ti partitioning between Quartz and hydrous
% high-silica melt using experiments at shallow crustal conditions. The
% safest direct experimental range for Equation (8) is:
%
%   Temperature : 700-900 degreeC
%   Pressure    : 0.5-4 kbar
%   Quartz Ti   : approximately 90-700 ppm elemental Ti
%   Liquid Ti   : approximately 839-3836 ppm elemental Ti
%   Melt type   : hydrous high-silica, rhyolitic, strongly peraluminous
%                 aluminosilicate melt
%   Melt SiO2   : high-silica melt; the associated rutile-solubility model
%                 is recommended for SiO2 >70 wt% after normalization to
%                 an anhydrous total of 100 wt%
%
% The overall experimental range is stated in the abstract on p. 1.
% Starting materials and melt compositions are described on p. 2.
% Experimental conditions and measured Quartz and glass compositions are
% listed in Table 1 on p. 3. Analytical methods are described on pp. 3-4.
% Attainment of near-equilibrium conditions is discussed on pp. 5-7.
% The rutile-saturated Quartz model is Equation (3) on p. 7.
% The TiO2-solubility model for high-silica melt is Equation (4) on
% pp. 8-9. The activity relation is Equation (7), and the Quartz-Liquid
% partitioning thermobarometer implemented here is Equation (8), on p. 9.
% Uncertainty and natural applications are discussed on pp. 9-11.
%
% IMPORTANT APPLICATION LIMITATIONS
%
%   1) Equation (8) was developed for Quartz coexisting with hydrous,
%      high-silica, rhyolitic melt at shallow crustal pressures. It was not
%      calibrated as a universal Quartz-Liquid thermometer for mafic,
%      intermediate, strongly alkaline, carbonatitic, metamorphic-fluid, or
%      hydrothermal systems (abstract and Introduction, pp. 1-2).
%
%   2) The direct paired Quartz-melt experiments span 700-900 degreeC and
%      0.5-4 kbar. Although the associated rutile-solubility model includes
%      compiled data extending to approximately 1000 degreeC and 10 kbar,
%      Equation (8) is most safely applied inside the direct Quartz-melt
%      range. Values outside 700-900 degreeC or 0.5-4 kbar are treated here
%      as extrapolations (Sections 4.2-4.5, pp. 6-10).
%
%   3) The selected Quartz domain and Liquid composition must represent a
%      coexisting equilibrium pair. A Quartz core must not be paired
%      automatically with a melt inclusion or matrix glass representing a
%      later growth stage. Post-entrapment crystallization, diffusive
%      modification, boundary-layer effects, decrepitation, or
%      re-equilibration can invalidate the pairing (Sections 4.1, 4.5, and
%      4.6, pp. 5-11).
%
%   4) Equation (8) eliminates the need to enter aTiO2 explicitly by
%      combining the Quartz Ti model with the melt TiO2-solubility model.
%      Its derivation nevertheless assumes ideal TiO2 activity behavior in
%      the silicate melt, i.e., a liquid TiO2 activity coefficient of
%      approximately 1. The equation is valid only when this assumption is
%      reasonable (Section 4.5, p. 9).
%
%   5) The analyzed Quartz Ti must represent Ti dissolved in the Quartz
%      lattice. Zhang et al. (2020) excluded experiments with high-Ti
%      starting glasses because rutile inclusions in Quartz prevented
%      reliable measurement of dissolved Ti (Starting materials, p. 2).
%      Rutile inclusions, rutile lamellae, Ti-oxide contamination, or glass
%      inclusions in an analytical volume can overestimate Quartz Ti.
%
%   6) The Liquid composition must represent the melt that coexisted with
%      the analyzed Quartz domain. Matrix glass affected by eruption,
%      degassing, crystallization, alteration, alkali loss, or mixing may
%      not reproduce the trapped equilibrium melt composition.
%
%   7) The associated TiO2-solubility model was recommended for silicic
%      melts containing more than 70 wt% SiO2 after normalization to an
%      anhydrous total of 100 wt% (Section 4.4, p. 9). Application to lower
%      SiO2 melts is a compositional extrapolation.
%
%   8) The experimental starting melts were Qz-Ab-Or-rich and strongly
%      peraluminous, with an aluminum saturation index of approximately
%      1.6-1.8 (Starting materials, p. 2). Similar numerical P-T values do
%      not guarantee applicability to substantially different melt
%      compositions.
%
%   9) The FM compositional parameter is:
%
%          Na + K + 2Ca + 2Mg + 2Fe
%     FM = --------------------------
%                    Si * Al
%
%      where the chemical symbols are cation molar fractions. Total Fe must
%      be counted once. This implementation uses FeO + Fe2O3 when a valid
%      split-iron analysis is available; otherwise it uses FeOt as total Fe
%      expressed as FeO. FeOt is never added to FeO + Fe2O3.
%
%  10) F and Cl are not included in cationTotal_liq, are not used in FM,
%      and are excluded from NaN-input diagnostics. Their possible presence
%      in the Liquid table therefore does not change this calculation.
%
%  11) Analytical precision is critical because Equation (8) uses the ratio
%      of Quartz Ti to Liquid Ti. The reported Ti detection limits were
%      approximately 38 ppm for Quartz by EPMA and approximately 10 ppm by
%      fs-LA-ICP-MS (Analytical methods, pp. 3-4). Small relative errors in
%      the Ti ratio can cause substantial P-T shifts. The abstract notes
%      that pressure precision of approximately +/-0.2 kbar for an input
%      temperature uncertainty of +/-25 degreeC requires high-precision Ti
%      measurements (p. 1).
%
%  12) Natural Quartz-Liquid Ti ratios are commonly approximately
%      0.09-0.13 in the examples examined by Zhang et al. (2020). At about
%      800 degreeC, this relatively small variation can correspond to a
%      pressure change from approximately 5 to 1 kbar (abstract, p. 1;
%      Section 4.6 and Fig. 7, p. 10). Pair selection and analytical
%      uncertainty must therefore be evaluated carefully.
%
%  13) Missing optional oxide columns are treated as components not supplied
%      and contribute zero to the cation total. In contrast, an oxide column
%      that is present but contains NaN remains NaN, is reported, and
%      propagates through cationTotal_liq, FM, and temperature. A measured
%      NaN is never silently replaced by zero.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside 0.5-4 kbar,
%   2) a finite calculated temperature is outside 700-900 degreeC,
%   3) anhydrous-normalized Liquid SiO2 is not greater than 70 wt%,
%   4) finite Quartz or Liquid Ti lies outside the approximate experimental
%      envelope,
%   5) an explicitly used Quartz or Liquid input is NaN,
%   6) the Ti conversion, cation total, cation fractions, FM, logarithm,
%      pressure term, denominator, or temperature term is invalid, or
%   7) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
%
% rawdata_struct must contain:
%   rawdata_struct.Quartz : table
%
% The FIRST column of the Quartz table is treated as an identifier
% ("data code") displayed in the selection dialog.
%
% Required Quartz variable:
%   Ti_cation_apfu
%
% Ti_cation_apfu must represent Ti atoms per formula unit normalized to
% O = 2. It is not a wt% value and is not multiplied directly by 1e4.
% It is converted to elemental Ti ppm using a binary TiO2-SiO2 mass
% relation before Equation (8) is evaluated.
%
% Required Liquid variables:
%   SiO2
%   TiO2
%   Al2O3
%   Na2O
%   K2O
%
% Optional Liquid variables:
%   CaO
%   MgO
%   FeO and/or Fe2O3
%   FeOt, FeOT, or equivalent total-Fe-as-FeO column
%   MnO
%   V2O3
%   Cr2O3
%   NiO
%   P2O5
%   SO3
%
% Column-name matching ignores spaces, underscores, and hyphens and accepts
% either the oxide name or the oxide name followed by "value".
%
% F and Cl columns, if present, are deliberately ignored.
%
% Existing NaN values in calculation inputs are retained and propagated;
% they are never replaced by zero. Finite calculation inputs must not be
% negative. Negative finite values and Inf stop the calculation. Zero is
% retained, but expressions requiring a positive Ti concentration, cation
% total, Si fraction, Al fraction, or denominator return NaN with
% non-stopping diagnostics.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Zhang et al. (2020), Equation (8), p. 9:
%
%   log10(C_Ti_Qtz / C_Ti_Liq) =
%       -1.1963
%       + (1058.1 - 520.4*P^0.2)/T
%       - 0.1155*FM
%
% Rearranged to calculate temperature:
%
%                 1058.1 - 520.4*P^0.2
%   T(K) = ------------------------------------------------
%          log10(C_Ti_Qtz/C_Ti_Liq) + 1.1963 + 0.1155*FM
%
%   T(degreeC) = T(K) - 273.15
%
% where:
%   C_Ti_Qtz : elemental Ti concentration in Quartz, ppm by weight
%   C_Ti_Liq : elemental Ti concentration in Liquid, ppm by weight
%   P        : pressure in kbar
%   T        : temperature in K
%
% Liquid TiO2 wt% is converted to elemental Ti ppm by:
%
%   C_Ti_Liq =
%       TiO2(wt%) * [M_Ti/M_TiO2] * 1e4
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Zhang2020TiLiq(rawdata_struct, P_kbar)
%   results = Zhang2020TiLiq(rawdata_struct, P_kbar, ...
%       'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing a Quartz table
%   P_kbar         : finite non-negative numeric scalar or vector
%
% Optional name-value input:
%   'LiquidRow'    : finite positive integer selecting one row from the
%                    Liquid dataset. Default: row 1.
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Quartz-Liquid pair. Generic T_K, T_degreeC, T_degC, and T_deg
%             variables are supplied in addition to the original
%             T_Zhang2020_Eq8_K and T_Zhang2020_Eq8_C names.
%

%% Input validation
if nargin < 2
    error('Zhang2020TiLiq requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

if ~isfield(rawdata_struct, 'Quartz') || ~istable(rawdata_struct.Quartz)
    error('rawdata_struct must contain table: rawdata_struct.Quartz');
end
if isempty(rawdata_struct.Quartz)
    error('rawdata_struct.Quartz is empty.');
end

dataset_qtz = rawdata_struct.Quartz;

if ~ismember('Ti_cation_apfu', dataset_qtz.Properties.VariableNames)
    error('Quartz table must contain variable: Ti_cation_apfu');
end

%% Optional inputs
ip = inputParser;
ip.addParameter('LiquidRow', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && ...
     x >= 1 && mod(x, 1) == 0));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Load molecular weights and Liquid dataset
disp('=== Step 1: Preparing Quartz and Liquid datasets ===');

MWinfo = liquid.getMolarWeights();

[liqAll, metaLiq] = liquid.readLiquidExcel();

if ~istable(liqAll) || isempty(liqAll)
    error('Selected Liquid dataset must be a non-empty table.');
end

if isempty(liquidRowOpt)
    idxLiq = 1;

    if height(liqAll) > 1
        fprintf(2, ...
            ['CAUTION: The selected Liquid dataset contains %d rows. ' ...
             'LiquidRow was not specified, so row 1 will be used for all ' ...
             'Quartz selections.\n'], ...
            height(liqAll));
    end
else
    idxLiq = liquidRowOpt;

    if idxLiq > height(liqAll)
        error(['Requested LiquidRow (%d) exceeds the number of rows in ' ...
               'the selected Liquid dataset (%d).'], ...
            idxLiq, height(liqAll));
    end
end

selectedData_liq = liqAll(idxLiq, :);
liquidInputs = extractLiquidInputs(selectedData_liq);

disp(['Liquid selected: Row ' num2str(idxLiq)]);
disp('=== Preparing Quartz and Liquid datasets has been finished ===');

%% 2) Initialize output container and calibration screening settings
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Safest direct paired Quartz-melt experimental limits.
calibrationT_min_degreeC = 700;
calibrationT_max_degreeC = 900;
calibrationP_min_kbar = 0.5;
calibrationP_max_kbar = 4;

% Approximate experimental composition envelopes from Table 1.
screeningQtzTi_min_ppm = 90;
screeningQtzTi_max_ppm = 700;
screeningLiqTi_min_ppm = 839;
screeningLiqTi_max_ppm = 3836;
screeningFM_min = 0.682;
screeningFM_max = 1.276;
screeningSiO2_min_anhydrous_wtpercent = 70;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

if isscalar(P_kbar)
    disp(['Pressure P = ' num2str(P_kbar) ' kbar']);
else
    disp(['Pressure P = ' num2str(P_kbar(1)) ' to ' ...
        num2str(P_kbar(end)) ' kbar']);
end

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive Quartz selection and calculation
dataCodes_qtz = dataset_qtz{:, 1};
displayCodes_qtz = cellstr(string(dataCodes_qtz));

disp('=== Step 3: Selecting a data code from the list (Quartz) ===');

while true
    [idxQtz, ok] = listdlg( ...
        'PromptString', ...
            'Please select the Quartz data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_qtz, ...
        'ListSize', [420 360]);

    if ~ok || isempty(idxQtz)
        disp('Selection canceled');
        break;
    end

    codeQtz = string(dataCodes_qtz(idxQtz));
    selectedData_qtz = dataset_qtz(idxQtz, :);

    disp(['Quartz selected: ' char(codeQtz)]);
    disp(['Liquid selected: Row ' num2str(idxLiq)]);
    disp('=== Step 4: Checking calculation inputs ===');

    % Report explicitly used NaN inputs without replacing them by zero.
    nanInputNames = findNaNInputs(selectedData_qtz, liquidInputs);

    % Stop only for negative finite values or Inf. NaN and zero are retained.
    validateNonNegativeInputs(selectedData_qtz, liquidInputs);

    disp('=== Step 5: Calculating the temperature ===');

    row = calcTemp( ...
        selectedData_qtz, liquidInputs, P_kbar, MWinfo);

    nRows = height(row);
    row.dataCode_qtz = repmat(codeQtz, nRows, 1);
    row.dataRow_liq = repmat(idxLiq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_qtz', 'dataRow_liq'}, 'Before', 1);

    % Store one completed table block. Repeated enlargement of the complete
    % output table inside the interactive loop is avoided.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = ...
            [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % ----- Immediate result display -----
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    printTemperatureSummary(codeQtz, idxLiq, row.T_degreeC);

    % Pressure is common to every selected pair in this function call, so
    % the pressure warning is printed only once.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the safest direct paired ' ...
             'Quartz-melt experimental range of Zhang et al. (2020): ' ...
             '0.5-4 kbar (abstract, p. 1; Table 1, p. 3). %d of %d ' ...
             'pressure point(s) are outside; input range = %.6g-%.6g ' ...
             'kbar. The supplied pressures and calculated results have ' ...
             'been retained as extrapolations.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    printTemperatureCalibrationWarning( ...
        row.T_degreeC, ...
        calibrationT_min_degreeC, ...
        calibrationT_max_degreeC, ...
        codeQtz, ...
        idxLiq);

    printCompositionScreeningWarnings( ...
        row, ...
        screeningQtzTi_min_ppm, ...
        screeningQtzTi_max_ppm, ...
        screeningLiqTi_min_ppm, ...
        screeningLiqTi_max_ppm, ...
        screeningFM_min, ...
        screeningFM_max, ...
        screeningSiO2_min_anhydrous_wtpercent, ...
        codeQtz, ...
        idxLiq);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Zhang et al. (2020) ' ...
             'Quartz-Liquid calculation input(s) for Quartz %s and ' ...
             'Liquid row %d: %s.\n' ...
             '         Existing NaN values were retained and propagated; ' ...
             'the calculated temperature may be NaN. F and Cl are not ' ...
             'calculation inputs and were excluded from this check.\n'], ...
            char(codeQtz), ...
            idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTerms = findInvalidEquationTerms(row);

    if ~isempty(invalidTerms)
        fprintf(2, ...
            ['WARNING: Invalid Zhang et al. (2020) Quartz-Liquid ' ...
             'conversion, cation-total, FM, logarithm, pressure, ' ...
             'denominator, or temperature term(s) were found for Quartz ' ...
             '%s and Liquid row %d: %s.\n' ...
             '         Affected values were retained as NaN, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(codeQtz), ...
            idxLiq, ...
            char(strjoin(invalidTerms, ', ')));
    end

    printNonFiniteTemperatureWarning( ...
        row.T_degreeC, codeQtz, idxLiq);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        ['Continue with another Quartz selection using Liquid row ' ...
         num2str(idxLiq) ' ?'], ...
        'Zhang2020TiLiq', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered result blocks only once after the selection loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'liquid', ...
        metaLiq, ...
    'liquidRow', ...
        idxLiq, ...
    'primaryTemperatureEquation', ...
        'Zhang et al. (2020) Equation (8)', ...
    'doi', ...
        'https://doi.org/10.1016/j.epsl.2020.116213', ...
    'experimentalCalibrationTemperatureRange_degreeC', ...
        [calibrationT_min_degreeC, calibrationT_max_degreeC], ...
    'experimentalCalibrationPressureRange_kbar', ...
        [calibrationP_min_kbar, calibrationP_max_kbar], ...
    'approximateQuartzTiRange_ppm', ...
        [screeningQtzTi_min_ppm, screeningQtzTi_max_ppm], ...
    'approximateLiquidTiRange_ppm', ...
        [screeningLiqTi_min_ppm, screeningLiqTi_max_ppm], ...
    'approximateExperimentalFMRange', ...
        [screeningFM_min, screeningFM_max], ...
    'recommendedMinimumAnhydrousSiO2_wtpercent', ...
        screeningSiO2_min_anhydrous_wtpercent, ...
    'pressureUsedInEquation', ...
        true, ...
    'liquidCompositionUsed', ...
        true, ...
    'FIncludedInCationTotal', ...
        false, ...
    'ClIncludedInCationTotal', ...
        false, ...
    'calibrationSystem', ...
        'hydrous high-silica rhyolitic peraluminous melt');

disp('=== Zhang2020TiLiq calculation has been finished! ===');

end

%% ---- local functions ----
function liquidInputs = extractLiquidInputs(data_liq)
% extractLiquidInputs
% Read one selected Liquid row while distinguishing:
%   1) an absent optional column, represented by zero contribution, and
%   2) a present column containing NaN, retained as NaN.
%
% F and Cl are intentionally not read because they are excluded from both
% cationTotal_liq and NaN-input diagnostics.

maxUsedInputs = 16;
usedNames = strings(maxUsedInputs, 1);
usedValues = NaN(maxUsedInputs, 1);
nUsed = 0;

[liquidInputs.SiO2, name, present] = ...
    getLiquidValue(data_liq, {'SiO2'}, true, NaN);
[nUsed, usedNames, usedValues] = registerUsedInput( ...
    nUsed, usedNames, usedValues, name, ...
    liquidInputs.SiO2, present, true);

[liquidInputs.TiO2, name, present] = ...
    getLiquidValue(data_liq, {'TiO2'}, true, NaN);
[nUsed, usedNames, usedValues] = registerUsedInput( ...
    nUsed, usedNames, usedValues, name, ...
    liquidInputs.TiO2, present, true);

[liquidInputs.Al2O3, name, present] = ...
    getLiquidValue(data_liq, {'Al2O3'}, true, NaN);
[nUsed, usedNames, usedValues] = registerUsedInput( ...
    nUsed, usedNames, usedValues, name, ...
    liquidInputs.Al2O3, present, true);

[liquidInputs.Na2O, name, present] = ...
    getLiquidValue(data_liq, {'Na2O'}, true, NaN);
[nUsed, usedNames, usedValues] = registerUsedInput( ...
    nUsed, usedNames, usedValues, name, ...
    liquidInputs.Na2O, present, true);

[liquidInputs.K2O, name, present] = ...
    getLiquidValue(data_liq, {'K2O'}, true, NaN);
[nUsed, usedNames, usedValues] = registerUsedInput( ...
    nUsed, usedNames, usedValues, name, ...
    liquidInputs.K2O, present, true);

optionalOxides = { ...
    'CaO', {'CaO'}; ...
    'MgO', {'MgO'}; ...
    'MnO', {'MnO'}; ...
    'V2O3', {'V2O3'}; ...
    'Cr2O3', {'Cr2O3'}; ...
    'NiO', {'NiO'}; ...
    'P2O5', {'P2O5'}; ...
    'SO3', {'SO3'}};

for i = 1:size(optionalOxides, 1)
    fieldName = optionalOxides{i, 1};
    aliases = optionalOxides{i, 2};

    [value, name, present] = ...
        getLiquidValue(data_liq, aliases, false, 0);

    liquidInputs.(fieldName) = value;

    [nUsed, usedNames, usedValues] = registerUsedInput( ...
        nUsed, usedNames, usedValues, name, value, present, false);
end

% Read alternative Fe representations without registering them yet. Only
% the representation selected below is included in the calculation and
% NaN-input diagnostics.
[FeO_raw, FeO_name, FeO_present] = ...
    getLiquidValue(data_liq, {'FeO'}, false, 0);
[Fe2O3_raw, Fe2O3_name, Fe2O3_present] = ...
    getLiquidValue(data_liq, {'Fe2O3'}, false, 0);
[FeOt_raw, FeOt_name, FeOt_present] = ...
    getLiquidValue(data_liq, ...
        {'FeOt', 'FeOT', 'FeOtot', 'FeOtotal', 'FeOTotal'}, ...
        false, 0);

splitFePresent = FeO_present || Fe2O3_present;
splitFeFinite = ...
    (~FeO_present || isfinite(FeO_raw)) && ...
    (~Fe2O3_present || isfinite(Fe2O3_raw));

if splitFePresent && splitFeFinite
    % Prefer a valid explicitly split FeO + Fe2O3 analysis.
    liquidInputs.FeSource = "FeO+Fe2O3";
    liquidInputs.FeO = FeO_raw;
    liquidInputs.Fe2O3 = Fe2O3_raw;
    liquidInputs.FeOt = FeOt_raw;

    [nUsed, usedNames, usedValues] = registerUsedInput( ...
        nUsed, usedNames, usedValues, FeO_name, ...
        FeO_raw, FeO_present, false);
    [nUsed, usedNames, usedValues] = registerUsedInput( ...
        nUsed, usedNames, usedValues, Fe2O3_name, ...
        Fe2O3_raw, Fe2O3_present, false);

elseif FeOt_present
    % Use total Fe as FeO when split Fe is unavailable or contains NaN.
    liquidInputs.FeSource = "FeOt";
    liquidInputs.FeO = FeO_raw;
    liquidInputs.Fe2O3 = Fe2O3_raw;
    liquidInputs.FeOt = FeOt_raw;

    [nUsed, usedNames, usedValues] = registerUsedInput( ...
        nUsed, usedNames, usedValues, FeOt_name, ...
        FeOt_raw, true, false);

elseif splitFePresent
    % No finite FeOt alternative exists; retain NaN in the split analysis.
    liquidInputs.FeSource = "FeO+Fe2O3";
    liquidInputs.FeO = FeO_raw;
    liquidInputs.Fe2O3 = Fe2O3_raw;
    liquidInputs.FeOt = 0;

    [nUsed, usedNames, usedValues] = registerUsedInput( ...
        nUsed, usedNames, usedValues, FeO_name, ...
        FeO_raw, FeO_present, false);
    [nUsed, usedNames, usedValues] = registerUsedInput( ...
        nUsed, usedNames, usedValues, Fe2O3_name, ...
        Fe2O3_raw, Fe2O3_present, false);

else
    liquidInputs.FeSource = "none";
    liquidInputs.FeO = 0;
    liquidInputs.Fe2O3 = 0;
    liquidInputs.FeOt = 0;
end

liquidInputs.usedNames = usedNames(1:nUsed);
liquidInputs.usedValues = usedValues(1:nUsed);

end

function [nUsed, usedNames, usedValues] = registerUsedInput( ...
        nUsed, usedNames, usedValues, sourceName, value, present, required)
% registerUsedInput
% Add a required input or a present optional input to fixed-size diagnostic
% buffers. Absent optional components contribute zero but are not reported
% as measured inputs.

if required || present
    nUsed = nUsed + 1;

    if nUsed > numel(usedNames)
        error('Internal Liquid-input diagnostic buffer capacity exceeded.');
    end

    if strlength(string(sourceName)) == 0
        sourceName = "unnamed";
    end

    usedNames(nUsed) = "Liquid." + string(sourceName);
    usedValues(nUsed) = value;
end

end

function nanInputNames = findNaNInputs(data_qtz, liquidInputs)
% findNaNInputs
% Return names of calculation inputs that are NaN. F and Cl are deliberately
% absent from liquidInputs.usedNames and are therefore excluded.

qtzValue = data_qtz.Ti_cation_apfu;
validateScalarVariable(qtzValue, 'Quartz', 'Ti_cation_apfu');

maxNames = 1 + numel(liquidInputs.usedNames);
nameBuffer = strings(maxNames, 1);
nNames = 0;

if isnan(qtzValue)
    nNames = nNames + 1;
    nameBuffer(nNames) = "Quartz.Ti_cation_apfu";
end

for i = 1:numel(liquidInputs.usedNames)
    if isnan(liquidInputs.usedValues(i))
        nNames = nNames + 1;
        nameBuffer(nNames) = liquidInputs.usedNames(i);
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_qtz, liquidInputs)
% validateNonNegativeInputs
% Stop for negative finite or infinite inputs used in the calculation.
% NaN and zero are deliberately allowed and handled by non-stopping
% diagnostics after calculation.

qtzValue = data_qtz.Ti_cation_apfu;
validateScalarVariable(qtzValue, 'Quartz', 'Ti_cation_apfu');

maxNames = 1 + numel(liquidInputs.usedNames);
invalidNames = strings(maxNames, 1);
nInvalid = 0;

if isinf(qtzValue) || (isfinite(qtzValue) && qtzValue < 0)
    nInvalid = nInvalid + 1;
    invalidNames(nInvalid) = "Quartz.Ti_cation_apfu";
end

for i = 1:numel(liquidInputs.usedNames)
    value = liquidInputs.usedValues(i);

    if isinf(value) || (isfinite(value) && value < 0)
        nInvalid = nInvalid + 1;
        invalidNames(nInvalid) = liquidInputs.usedNames(i);
    end
end

if nInvalid > 0
    invalidNames = invalidNames(1:nInvalid);

    error([ ...
        'Zhang2020TiLiq: calculation inputs must not be negative or ' ...
        'infinite. Invalid value(s) were found in: ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_qtz, liquidInputs, P_kbar, MWinfo)
% calcTemp
% Calculate Zhang et al. (2020) Equation (8) for one selected Quartz row,
% one selected Liquid row, and every supplied pressure. Existing NaN values
% and invalid derived terms are retained as NaN.
%
% F and Cl are not included in the Liquid cation total.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
% Number of supplied pressure points. Keep this name distinct from
% nP_cation, the P2O5-derived phosphorus cation amount below.
nPressure = numel(P_kbar);

% --- Atomic and molecular weights ---
M_Ti = 47.867;
M_TiO2 = MWinfo.MW.TiO2;
M_SiO2 = MWinfo.MW.SiO2;

% --- Quartz Ti conversion ---
X_TiO2_qtz_scalar = data_qtz.Ti_cation_apfu;
validateScalarVariable( ...
    X_TiO2_qtz_scalar, 'Quartz', 'Ti_cation_apfu');

Ti_Qtz_ppm_scalar = NaN;
Ti_Qtz_conversion_valid_scalar = false;

if isnan(X_TiO2_qtz_scalar)
    Ti_Qtz_conversion_valid_scalar = false;
elseif isfinite(X_TiO2_qtz_scalar) && ...
        X_TiO2_qtz_scalar >= 0 && X_TiO2_qtz_scalar <= 1

    formulaMass_qtz = ...
        X_TiO2_qtz_scalar .* M_TiO2 + ...
        (1 - X_TiO2_qtz_scalar) .* M_SiO2;

    if isfinite(formulaMass_qtz) && formulaMass_qtz > 0
        Ti_Qtz_ppm_scalar = ...
            (X_TiO2_qtz_scalar .* M_Ti ./ formulaMass_qtz) .* 1e6;

        Ti_Qtz_conversion_valid_scalar = ...
            isfinite(Ti_Qtz_ppm_scalar) && ...
            Ti_Qtz_ppm_scalar >= 0;
    end
end

% --- Liquid Ti conversion ---
Ti_Liq_ppm_scalar = NaN;
Ti_Liq_conversion_valid_scalar = false;

if isnan(liquidInputs.TiO2)
    Ti_Liq_conversion_valid_scalar = false;
elseif isfinite(liquidInputs.TiO2) && liquidInputs.TiO2 >= 0
    Ti_Liq_ppm_scalar = ...
        liquidInputs.TiO2 .* (M_Ti ./ M_TiO2) .* 1e4;

    Ti_Liq_conversion_valid_scalar = ...
        isfinite(Ti_Liq_ppm_scalar) && Ti_Liq_ppm_scalar >= 0;
end

% --- Liquid cation molar amounts ---
nSi = liquidInputs.SiO2 ./ MWinfo.MW.SiO2;
nTi = liquidInputs.TiO2 ./ MWinfo.MW.TiO2;
nAl = 2 .* liquidInputs.Al2O3 ./ MWinfo.MW.Al2O3;
nMn = liquidInputs.MnO ./ MWinfo.MW.MnO;
nMg = liquidInputs.MgO ./ MWinfo.MW.MgO;
nCa = liquidInputs.CaO ./ MWinfo.MW.CaO;
nNa = 2 .* liquidInputs.Na2O ./ MWinfo.MW.Na2O;
nK = 2 .* liquidInputs.K2O ./ MWinfo.MW.K2O;
nV = 2 .* liquidInputs.V2O3 ./ MWinfo.MW.V2O3;
nCr = 2 .* liquidInputs.Cr2O3 ./ MWinfo.MW.Cr2O3;
nNi = liquidInputs.NiO ./ MWinfo.MW.NiO;
% Phosphorus cation amount from P2O5. Do not name this nP because
% nPressure stores the number of pressure points used for allocation.
nP_cation = 2 .* liquidInputs.P2O5 ./ MWinfo.MW.P2O5;
nS = liquidInputs.SO3 ./ MWinfo.MW.SO3;

if strcmp(liquidInputs.FeSource, "FeOt")
    nFe = liquidInputs.FeOt ./ MWinfo.MW.FeO;
    FeO_used = liquidInputs.FeOt;
    Fe2O3_used = 0;
elseif strcmp(liquidInputs.FeSource, "FeO+Fe2O3")
    nFe = ...
        liquidInputs.FeO ./ MWinfo.MW.FeO + ...
        2 .* liquidInputs.Fe2O3 ./ MWinfo.MW.Fe2O3;
    FeO_used = liquidInputs.FeO;
    Fe2O3_used = liquidInputs.Fe2O3;
else
    nFe = 0;
    FeO_used = 0;
    Fe2O3_used = 0;
end

% F and Cl are deliberately excluded here.
cationTotal_scalar = ...
    nSi + nTi + nAl + nFe + nMn + nMg + nCa + nNa + nK + ...
    nV + nCr + nNi + nP_cation + nS;

% --- Cation molar fractions and FM ---
XSi_scalar = NaN;
XAl_scalar = NaN;
XNa_scalar = NaN;
XK_scalar = NaN;
XCa_scalar = NaN;
XMg_scalar = NaN;
XFe_scalar = NaN;
FM_scalar = NaN;

cationTotal_valid_scalar = ...
    isfinite(cationTotal_scalar) && cationTotal_scalar > 0;

if cationTotal_valid_scalar
    XSi_scalar = nSi ./ cationTotal_scalar;
    XAl_scalar = nAl ./ cationTotal_scalar;
    XNa_scalar = nNa ./ cationTotal_scalar;
    XK_scalar = nK ./ cationTotal_scalar;
    XCa_scalar = nCa ./ cationTotal_scalar;
    XMg_scalar = nMg ./ cationTotal_scalar;
    XFe_scalar = nFe ./ cationTotal_scalar;

    validFMInputs = ...
        isfinite(XSi_scalar) && XSi_scalar > 0 && ...
        isfinite(XAl_scalar) && XAl_scalar > 0 && ...
        isfinite(XNa_scalar) && ...
        isfinite(XK_scalar) && ...
        isfinite(XCa_scalar) && ...
        isfinite(XMg_scalar) && ...
        isfinite(XFe_scalar);

    if validFMInputs
        FM_scalar = ...
            (XNa_scalar + XK_scalar + ...
             2 .* XCa_scalar + 2 .* XMg_scalar + 2 .* XFe_scalar) ./ ...
            (XSi_scalar .* XAl_scalar);
    end
end

% --- Anhydrous oxide total and normalized SiO2 screening value ---
anhydrousOxideTotal_scalar = ...
    liquidInputs.SiO2 + liquidInputs.TiO2 + liquidInputs.Al2O3 + ...
    liquidInputs.Na2O + liquidInputs.K2O + liquidInputs.CaO + ...
    liquidInputs.MgO + FeO_used + Fe2O3_used + ...
    liquidInputs.MnO + liquidInputs.V2O3 + liquidInputs.Cr2O3 + ...
    liquidInputs.NiO + liquidInputs.P2O5 + liquidInputs.SO3;

SiO2_anhydrous_wtpercent_scalar = NaN;

if isfinite(anhydrousOxideTotal_scalar) && ...
        anhydrousOxideTotal_scalar > 0 && ...
        isfinite(liquidInputs.SiO2)

    SiO2_anhydrous_wtpercent_scalar = ...
        liquidInputs.SiO2 ./ anhydrousOxideTotal_scalar .* 100;
end

% --- Repeat scalar composition values for every pressure ---
X_TiO2_qtz = repmat(X_TiO2_qtz_scalar, nPressure, 1);
Ti_Qtz_ppm = repmat(Ti_Qtz_ppm_scalar, nPressure, 1);
Ti_Qtz_conversion_valid = ...
    repmat(Ti_Qtz_conversion_valid_scalar, nPressure, 1);

Ti_Liq_ppm = repmat(Ti_Liq_ppm_scalar, nPressure, 1);
Ti_Liq_conversion_valid = ...
    repmat(Ti_Liq_conversion_valid_scalar, nPressure, 1);

cationTotal_liq = repmat(cationTotal_scalar, nPressure, 1);
cationTotal_valid = repmat(cationTotal_valid_scalar, nPressure, 1);

XSi = repmat(XSi_scalar, nPressure, 1);
XAl = repmat(XAl_scalar, nPressure, 1);
XNa = repmat(XNa_scalar, nPressure, 1);
XK = repmat(XK_scalar, nPressure, 1);
XCa = repmat(XCa_scalar, nPressure, 1);
XMg = repmat(XMg_scalar, nPressure, 1);
XFe = repmat(XFe_scalar, nPressure, 1);
FM = repmat(FM_scalar, nPressure, 1);

anhydrousOxideTotal = ...
    repmat(anhydrousOxideTotal_scalar, nPressure, 1);
SiO2_anhydrous_wtpercent = ...
    repmat(SiO2_anhydrous_wtpercent_scalar, nPressure, 1);

% --- Equation (8) ---
Ti_Qtz_over_Ti_Liq = NaN(nPressure, 1);
log_ratio = NaN(nPressure, 1);
P_power_0p2 = NaN(nPressure, 1);
numerator = NaN(nPressure, 1);
denominator = NaN(nPressure, 1);
T_raw_K = NaN(nPressure, 1);
T_raw_degreeC = NaN(nPressure, 1);
T_K = NaN(nPressure, 1);
T_degreeC = NaN(nPressure, 1);

P_power_0p2(:) = P_kbar .^ 0.2;
numerator(:) = 1058.1 - 520.4 .* P_power_0p2;

validTiRatioInputs = ...
    isfinite(Ti_Qtz_ppm) & Ti_Qtz_ppm > 0 & ...
    isfinite(Ti_Liq_ppm) & Ti_Liq_ppm > 0;

Ti_Qtz_over_Ti_Liq(validTiRatioInputs) = ...
    Ti_Qtz_ppm(validTiRatioInputs) ./ ...
    Ti_Liq_ppm(validTiRatioInputs);

validTiRatio = ...
    validTiRatioInputs & ...
    isfinite(Ti_Qtz_over_Ti_Liq) & ...
    Ti_Qtz_over_Ti_Liq > 0;

log_ratio(validTiRatio) = ...
    log10(Ti_Qtz_over_Ti_Liq(validTiRatio));

validDenominatorInputs = ...
    validTiRatio & ...
    isfinite(FM) & ...
    isfinite(P_power_0p2) & ...
    isfinite(numerator);

denominator(validDenominatorInputs) = ...
    log_ratio(validDenominatorInputs) + ...
    1.1963 + ...
    0.1155 .* FM(validDenominatorInputs);

validDenominator = ...
    validDenominatorInputs & ...
    isfinite(denominator) & ...
    abs(denominator) > 1e-12;

T_raw_K(validDenominator) = ...
    numerator(validDenominator) ./ denominator(validDenominator);
T_raw_degreeC(validDenominator) = ...
    T_raw_K(validDenominator) - 273.15;

% Non-positive Kelvin is physically invalid. Preserve raw solutions for
% diagnosis, but return NaN as the accepted temperature.
validTemperature = ...
    isfinite(T_raw_K) & T_raw_K > 0;

T_K(validTemperature) = T_raw_K(validTemperature);
T_degreeC(validTemperature) = T_raw_degreeC(validTemperature);

% --- Pack outputs ---
row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.P_power_0p2 = P_power_0p2;
row.PressureUsedInEquation = true(nPressure, 1);

row.X_TiO2_qtz = X_TiO2_qtz;
row.Ti_cation_apfu_qtz = X_TiO2_qtz;
row.Ti_Qtz_ppm = Ti_Qtz_ppm;
row.Ti_Qtz_conversion_valid = Ti_Qtz_conversion_valid;

row.SiO2_liq = repmat(liquidInputs.SiO2, nPressure, 1);
row.TiO2_liq = repmat(liquidInputs.TiO2, nPressure, 1);
row.Al2O3_liq = repmat(liquidInputs.Al2O3, nPressure, 1);
row.Na2O_liq = repmat(liquidInputs.Na2O, nPressure, 1);
row.K2O_liq = repmat(liquidInputs.K2O, nPressure, 1);
row.CaO_liq = repmat(liquidInputs.CaO, nPressure, 1);
row.MgO_liq = repmat(liquidInputs.MgO, nPressure, 1);
row.FeO_liq = repmat(liquidInputs.FeO, nPressure, 1);
row.Fe2O3_liq = repmat(liquidInputs.Fe2O3, nPressure, 1);
row.FeOt_liq = repmat(liquidInputs.FeOt, nPressure, 1);
row.FeSource_liq = repmat(string(liquidInputs.FeSource), nPressure, 1);
row.MnO_liq = repmat(liquidInputs.MnO, nPressure, 1);
row.V2O3_liq = repmat(liquidInputs.V2O3, nPressure, 1);
row.Cr2O3_liq = repmat(liquidInputs.Cr2O3, nPressure, 1);
row.NiO_liq = repmat(liquidInputs.NiO, nPressure, 1);
row.P2O5_liq = repmat(liquidInputs.P2O5, nPressure, 1);
row.SO3_liq = repmat(liquidInputs.SO3, nPressure, 1);

row.Ti_Liq_ppm = Ti_Liq_ppm;
row.Ti_Liq_conversion_valid = Ti_Liq_conversion_valid;

row.cationTotal_liq = cationTotal_liq;
row.cationTotal_liq_valid = cationTotal_valid;
row.F_included_in_cationTotal_liq = false(nPressure, 1);
row.Cl_included_in_cationTotal_liq = false(nPressure, 1);

row.XSi_liq = XSi;
row.XAl_liq = XAl;
row.XNa_liq = XNa;
row.XK_liq = XK;
row.XCa_liq = XCa;
row.XMg_liq = XMg;
row.XFe_liq = XFe;
row.FM = FM;

row.anhydrousOxideTotal_liq = anhydrousOxideTotal;
row.SiO2_anhydrous_wtpercent_liq = ...
    SiO2_anhydrous_wtpercent;

row.Ti_Qtz_over_Ti_Liq = Ti_Qtz_over_Ti_Liq;
row.log_TiQtz_over_TiLiq = log_ratio;
row.numerator_Eq8 = numerator;
row.denominator_Eq8 = denominator;

row.T_raw_K = T_raw_K;
row.T_raw_degreeC = T_raw_degreeC;
row.T_Zhang2020_Eq8_K = T_K;
row.T_Zhang2020_Eq8_C = T_degreeC;
row.T_K = T_K;
row.T_degreeC = T_degreeC;
row.T_degC = T_degreeC;
row.T_deg = T_degreeC;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid conversion, cation-total, FM, logarithm, denominator, and
% temperature terms. Repeated pressure rows are summarized by term name.

termBuffer = strings(15, 1);
nTerms = 0;

if any(~row.Ti_Qtz_conversion_valid)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Quartz Ti_cation_apfu-to-Ti_ppm conversion";
end

if any(~isfinite(row.X_TiO2_qtz) | ...
        row.X_TiO2_qtz <= 0 | row.X_TiO2_qtz > 1)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = ...
        "X_TiO2_qtz (must satisfy 0 < X <= 1)";
end

if any(~isfinite(row.Ti_Qtz_ppm) | row.Ti_Qtz_ppm <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Ti_Qtz_ppm (> 0 required)";
end

if any(~row.Ti_Liq_conversion_valid)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Liquid TiO2-to-elemental-Ti conversion";
end

if any(~isfinite(row.Ti_Liq_ppm) | row.Ti_Liq_ppm <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Ti_Liq_ppm (> 0 required)";
end

if any(~row.cationTotal_liq_valid | ...
        ~isfinite(row.cationTotal_liq) | ...
        row.cationTotal_liq <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "cationTotal_liq";
end

if any(~isfinite(row.XSi_liq) | row.XSi_liq <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XSi_liq (> 0 required)";
end

if any(~isfinite(row.XAl_liq) | row.XAl_liq <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XAl_liq (> 0 required)";
end

if any(~isfinite(row.FM))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "FM";
end

if any(~isfinite(row.Ti_Qtz_over_Ti_Liq) | ...
        row.Ti_Qtz_over_Ti_Liq <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Ti_Qtz_ppm/Ti_Liq_ppm";
end

if any(~isfinite(row.log_TiQtz_over_TiLiq))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "log10(Ti_Qtz_ppm/Ti_Liq_ppm)";
end

if any(~isfinite(row.P_power_0p2))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "P_kbar^0.2";
end

if any(~isfinite(row.numerator_Eq8))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Equation (8) numerator";
end

if any(~isfinite(row.denominator_Eq8) | ...
        abs(row.denominator_Eq8) <= 1e-12)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Equation (8) denominator";
end

if any(~isfinite(row.T_K) | row.T_K <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "accepted T_K";
end

invalidTerms = unique(termBuffer(1:nTerms), 'stable');

end

function printTemperatureSummary(codeQtz, idxLiq, temperatureValues)
% printTemperatureSummary
% Display one temperature or the first-to-last values for a pressure vector.

label = [char(codeQtz) ' & Liquid row ' num2str(idxLiq)];

if isscalar(temperatureValues)
    disp([label ': ' num2str(temperatureValues) ' degreeC']);
else
    disp([label ': ' ...
        num2str(temperatureValues(1)) ' to ' ...
        num2str(temperatureValues(end)) ' degreeC']);
end

end

function printTemperatureCalibrationWarning( ...
        temperatureValues, minimumTemperature, maximumTemperature, ...
        codeQtz, idxLiq)
% printTemperatureCalibrationWarning
% Warn when finite temperatures lie outside the safest direct paired
% Quartz-melt experimental range. Results are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);

    fprintf(2, ...
        ['WARNING: Calculated Zhang et al. (2020) Equation (8) ' ...
         'temperature is outside the safest direct paired Quartz-melt ' ...
         'experimental range of %.4g-%.4g degreeC (abstract, p. 1; ' ...
         'Table 1, p. 3). %d of %d finite point(s) are outside; ' ...
         'calculated finite range = %.6g-%.6g degreeC for Quartz %s ' ...
         'and Liquid row %d. The result has been retained as an ' ...
         'extrapolation.\n'], ...
        minimumTemperature, ...
        maximumTemperature, ...
        sum(outsideMask), ...
        sum(finiteMask), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(codeQtz), ...
        idxLiq);
end

end

function printCompositionScreeningWarnings( ...
        row, qtzTiMin, qtzTiMax, liqTiMin, liqTiMax, ...
        FMmin, FMmax, minimumSiO2, codeQtz, idxLiq)
% printCompositionScreeningWarnings
% Print non-stopping composition cautions. These are empirical dataset
% envelopes, not strict independent rectangular calibration limits.

finiteQtzTi = row.Ti_Qtz_ppm(isfinite(row.Ti_Qtz_ppm) & ...
    row.Ti_Qtz_ppm > 0);

if ~isempty(finiteQtzTi) && ...
        any(finiteQtzTi < qtzTiMin | finiteQtzTi > qtzTiMax)

    fprintf(2, ...
        ['CAUTION: Quartz elemental Ti is outside the approximate ' ...
         '%.6g-%.6g ppm experimental envelope in Zhang et al. (2020, ' ...
         'Table 1, p. 3). Input = %.6g ppm for Quartz %s. The result ' ...
         'has been retained.\n'], ...
        qtzTiMin, qtzTiMax, finiteQtzTi(1), char(codeQtz));
end

finiteLiqTi = row.Ti_Liq_ppm(isfinite(row.Ti_Liq_ppm) & ...
    row.Ti_Liq_ppm > 0);

if ~isempty(finiteLiqTi) && ...
        any(finiteLiqTi < liqTiMin | finiteLiqTi > liqTiMax)

    fprintf(2, ...
        ['CAUTION: Liquid elemental Ti is outside the approximate ' ...
         '%.6g-%.6g ppm experimental envelope in Zhang et al. (2020, ' ...
         'Table 1, p. 3). Input = %.6g ppm for Liquid row %d. The ' ...
         'result has been retained.\n'], ...
        liqTiMin, liqTiMax, finiteLiqTi(1), idxLiq);
end

finiteFM = row.FM(isfinite(row.FM));

if ~isempty(finiteFM) && ...
        any(finiteFM < FMmin | finiteFM > FMmax)

    fprintf(2, ...
        ['CAUTION: Liquid FM is outside the approximate %.6g-%.6g range ' ...
         'represented by the experiments in Zhang et al. (2020, Table 1, ' ...
         'p. 3). Input FM = %.6g for Liquid row %d. The result has been ' ...
         'retained.\n'], ...
        FMmin, FMmax, finiteFM(1), idxLiq);
end

finiteSiO2 = row.SiO2_anhydrous_wtpercent_liq( ...
    isfinite(row.SiO2_anhydrous_wtpercent_liq));

if ~isempty(finiteSiO2) && any(finiteSiO2 <= minimumSiO2)
    fprintf(2, ...
        ['CAUTION: Anhydrous-normalized Liquid SiO2 is %.6g wt%% for ' ...
         'Liquid row %d. Zhang et al. (2020) recommend the associated ' ...
         'high-silica melt model for SiO2 > %.6g wt%% (Section 4.4, ' ...
         'p. 9). The result has been retained as a compositional ' ...
         'extrapolation.\n'], ...
        finiteSiO2(1), idxLiq, minimumSiO2);
end

end

function printNonFiniteTemperatureWarning( ...
        temperatureValues, codeQtz, idxLiq)
% printNonFiniteTemperatureWarning
% Report NaN and Inf temperatures without modifying or deleting them.

invalidMask = ~isfinite(temperatureValues);

if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Zhang et al. (2020) Equation (8) ' ...
         'temperature values were calculated for Quartz %s and Liquid ' ...
         'row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the ' ...
         'calculation has not been stopped.\n'], ...
        char(codeQtz), ...
        idxLiq, ...
        sum(invalidMask), ...
        numel(temperatureValues), ...
        sum(isnan(temperatureValues)), ...
        sum(isinf(temperatureValues)));
end

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Attach common Liquid identifiers and repeat them to match the number of
% pressure rows.

nRows = height(row);
vn = data_liq.Properties.VariableNames;

if any(strcmp(vn, 'Index'))
    value = data_liq.('Index');
    row.liq_Index = repeatTableScalar(value, nRows, false);
end
if any(strcmp(vn, 'Experiment'))
    value = data_liq.('Experiment');
    row.liq_Experiment = repeatTableScalar(value, nRows, true);
end
if any(strcmp(vn, 'Citation'))
    value = data_liq.('Citation');
    row.liq_Citation = repeatTableScalar(value, nRows, true);
end
if any(strcmp(vn, 'Sample'))
    value = data_liq.('Sample');
    row.liq_Sample = repeatTableScalar(value, nRows, true);
end

end

function output = repeatTableScalar(raw, nRows, forceString)
% repeatTableScalar
% Convert one table-row value to a repeated column suitable for assignment
% to an nRows-high table.

if forceString
    output = repmat(string(raw(1)), nRows, 1);
    return;
end

if isnumeric(raw) || islogical(raw)
    output = repmat(raw(1), nRows, 1);
else
    output = repmat(string(raw(1)), nRows, 1);
end

end

function [value, sourceName, present] = getLiquidValue( ...
        data_liq, aliases, required, absentDefault)
% getLiquidValue
% Locate an oxide column using one or more aliases. A present NaN is
% retained. An absent optional column returns absentDefault.

sourceName = findOxideColumn( ...
    data_liq.Properties.VariableNames, aliases);
present = ~isempty(sourceName);

if ~present
    if required
        error('Selected Liquid row must contain variable: %s', aliases{1});
    end

    value = absentDefault;
    sourceName = aliases{1};
    return;
end

value = toScalarDoublePreserveNaN(data_liq.(sourceName));

end

function name = findOxideColumn(varNames, aliases)
% findOxideColumn
% Match aliases after removing spaces, underscores, and hyphens. Either the
% canonical oxide name or the oxide name followed by "value" is accepted.

canon = cell(size(varNames));

for i = 1:numel(varNames)
    canon{i} = canonicalName(varNames{i});
end

name = '';

for a = 1:numel(aliases)
    ox = canonicalName(aliases{a});
    targets = {[ox 'value'], ox};

    for t = 1:numel(targets)
        idx = find(strcmp(canon, targets{t}), 1, 'first');

        if ~isempty(idx)
            name = varNames{idx};
            return;
        end
    end
end

end

function output = canonicalName(inputName)
% canonicalName
% Canonicalize a table-variable name for tolerant oxide matching.

output = lower(char(string(inputName)));
output = strrep(output, ' ', '');
output = strrep(output, '_', '');
output = strrep(output, '-', '');

end

function value = toScalarDoublePreserveNaN(raw)
% toScalarDoublePreserveNaN
% Convert one table cell/value to a numeric scalar. Empty, missing, or
% non-convertible values become NaN and are retained for diagnostics.

value = NaN;

if isempty(raw)
    return;
end

if isnumeric(raw) || islogical(raw)
    value = double(raw(1));
    return;
end

if isstring(raw)
    if ismissing(raw(1))
        return;
    end

    value = str2double(raw(1));
    return;
end

if ischar(raw)
    value = str2double(string(raw));
    return;
end

if iscell(raw)
    if isempty(raw{1})
        return;
    end

    item = raw{1};

    if isnumeric(item) || islogical(item)
        value = double(item(1));
        return;
    end

    if isstring(item)
        if ismissing(item(1))
            return;
        end

        value = str2double(item(1));
        return;
    end

    if ischar(item)
        value = str2double(string(item));
        return;
    end
end

end

function validateScalarVariable(value, mineralLabel, variableName)
% validateScalarVariable
% Require one numeric scalar from the selected table row. NaN is allowed;
% negative finite values and Inf are handled separately.

if ~isnumeric(value) || ~isscalar(value)
    error('%s variable %s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end
