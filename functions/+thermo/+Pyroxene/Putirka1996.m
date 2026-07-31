function results = Putirka1996(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Pyroxene/Putirka1996.m
% Tested with MATLAB R2024b
%
% Clinopyroxene-liquid thermometers T1-T4
% Putirka, K., Johnson, M., Kinzler, R., Longhi, J. and Walker, D. (1996)
% Contributions to Mineralogy and Petrology, 123, 92-108
% DOI: https://doi.org/10.1007/s004100050145
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Clinopyroxene analysis and pairs
% it with one liquid composition loaded by liquid.readLiquidExcel. It then
% calculates the four Cpx-liquid temperatures T1-T4 listed in Table 5 of
% Putirka et al. (1996).
%
% The function accepts either a scalar pressure or a pressure vector and is
% therefore compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP. For every selected Cpx-liquid pair, the output
% contains one row per input pressure value.
%
% T1 and T3 are pressure independent and are repeated for every pressure
% row. T2 and T4 contain pressure terms and are calculated independently at
% every pressure. Because T2 is the pressure-dependent DiHd-Jd thermometer
% used in the natural-sample application and generally performs better than
% T1, the standard output columns T_K and T_deg are assigned from T2. All
% four original outputs are retained as T1_K/T1_C through T4_K/T4_C.
%
% The function is designed for repeated calculations. Each result block is
% stored in a fixed-size preallocated cell buffer and all blocks are
% concatenated once after the interactive loop. No result array is enlarged
% inside the loop.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka et al. (1996) calibrated Cpx-liquid thermometers for mafic igneous
% compositions. Their new experiments covered:
%
%   Pressure    : 8-30 kbar
%   Temperature : 1100-1475 degreeC
%
% The final T1-T4 regressions also include Walter and Presnall (1994) data
% and selected 1-bar experiments. The title and intended principal pressure
% range are 0-30 kbar. Independent test data span approximately 1 bar to
% 50 kbar. Therefore 30-50 kbar should be treated as tested extrapolation,
% whereas pressures above 50 kbar are beyond the reported test range.
%
% The 1100-1475 degreeC interval is the range of the new experiments, not a
% strict rectangular limit for every datum included in the final regression.
% Temperatures outside this interval are reported as outside the principal
% new-experiment range but are not automatically rejected.
%
% Model performance reported in the paper includes:
%
%   T1: regression SEE = 26.8 K; pressure independent; DiHd-Jd exchange
%   T2: regression SEE = 23.8 K; pressure dependent; DiHd-Jd exchange
%   T3: regression SEE = 37.5 K; pressure independent; DiHd-CaTs exchange
%   T4: regression SEE = 28.1 K; pressure dependent; DiHd-CaTs exchange
%
% The independent test-data errors are approximately 40 K for the
% pressure-independent models and approximately 30-34 K for the
% pressure-dependent models. Natural applications should therefore use a
% practical uncertainty of at least approximately 30-40 degreeC rather than
% relying only on the regression SEE.
%
% Relevant locations in the original paper:
%
%   p. 92       : abstract; new-experiment P-T range and model uncertainties
%   pp. 95-97   : experimental methods, quench effects, and equilibrium rims
%   pp. 97-99   : liquid cation fractions and Cpx component calculations
%   p. 100      : Table 5, equations T1-T4
%   pp. 100-104 : regression and independent-test performance
%   pp. 105-106 : natural-sample application and whole-rock/liquid checks
%
% Important application cautions:
%
%   1) The selected Cpx and liquid must represent an equilibrium pair.
%      Putirka et al. (1996) found that disequilibrium Cpx compositions were
%      probably the largest source of experimental model error. Cpx rims in
%      direct contact with homogeneous glass were preferred (pp. 96-97,
%      102). Pairing unrelated cores, rims, antecrysts, xenocrysts, or altered
%      glass may produce geologically meaningless temperatures.
%
%   2) A whole-rock analysis may be used as a liquid only when it reasonably
%      represents the melt from which the selected Cpx crystallized. Crystal
%      accumulation, fractionation, mixing, or open-system behavior can
%      invalidate the pairing. The Mauna Kea application tested and, where
%      justified, corrected whole-rock compositions before thermobarometry
%      (pp. 105-106).
%
%   3) Liquid components must be calculated as cation fractions. Cpx cations
%      and components must follow the 6-oxygen normative scheme described on
%      pp. 97-99. Mixing a different component-allocation scheme with the
%      Table 5 equations will not reproduce the original calibration.
%
%   4) T1 and T2 require positive Jd and DiHd Cpx components and positive
%      Na, Al, Ca, and Fe+Mg liquid cation fractions. T3 and T4 additionally
%      require a positive CaTs component. When CaTs is zero, T3 and T4 are
%      mathematically undefined and remain NaN.
%
%   5) The experimental Cpx contained little Cr and minimal Fe3+. The paper
%      did not assign a Cr-bearing Cpx component, and acmite was considered
%      negligible because the experiments used graphite capsules (p. 99).
%      Strongly oxidized or compositionally unusual Cpx-liquid pairs are
%      therefore extrapolative.
%
%   6) At 1 bar, Na volatilization during long experiments, sluggish
%      equilibration below approximately 1050 degreeC, and sector zoning can
%      obscure equilibrium compositions (pp. 102-103). Similar analytical or
%      petrological effects should be considered in natural samples.
%
% This implementation issues non-stopping fprintf warnings when:
%
%   1) pressure is outside 0-30 kbar;
%   2) pressure is above the approximately 50-kbar independent-test limit;
%   3) a finite T1-T4 result is outside 1100-1475 degreeC;
%   4) any value actually used in the calculation is NaN;
%   5) a logarithm or inverse-temperature term is outside its domain; or
%   6) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST Cpx-table column is treated as an identifier ("data code") for
% the selection dialog. Cpx analyses are expected as oxide wt.% columns.
% Column names may be either oxide names such as SiO2 or standard liquid/
% mineral names such as SiO2Value.
%
% Required Cpx oxides:
%   SiO2, Al2O3, MgO, CaO, Na2O, and FeO or FeOt
%
% Optional Cpx oxides, treated as zero only when the column is absent:
%   TiO2, MnO, K2O, Cr2O3
%
% The liquid dataset is loaded with:
%   [liqAll, metaLiq] = liquid.readLiquidExcel();
%
% Required liquid oxides:
%   SiO2, Al2O3, MgO, CaO, Na2O, and FeO or FeOt
%
% Optional liquid oxides, treated as zero only when the column is absent:
%   TiO2, MnO, K2O, V2O3, Cr2O3, NiO, P2O5, SO3, Fe2O3, F, Cl
%
% If an oxide column is present and its selected value is NaN, the NaN is
% retained and propagated. It is never replaced by zero. If both FeO and
% FeOt columns exist, FeO is used; a present NaN FeO value is retained and
% is not replaced by FeOt. FeOt is used only when the FeO column is absent.
%
% F and Cl are anions. They are retained as raw output values for
% traceability but are excluded from cationTotal_liq, from temperature
% calculation, from NaN-input warnings, and from negative-value validation.
%
% Every finite oxide value used in the Cpx or liquid calculation must be
% greater than or equal to zero. Finite negative values and Inf are rejected.
% Zero and NaN are allowed as raw values; invalid equation domains return
% NaN temperatures and non-stopping warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Liquid oxides are converted to cation fractions. Define:
%
%   Fm_liq    = XFeO_liq + XMgO_liq
%   MgNum_liq = XMgO_liq / (XMgO_liq + XFeO_liq)
%
% Cpx components are calculated on a 6-oxygen basis:
%
%   AlIV  = max(2 - Si, 0)
%   AlVI  = max(Al_total - AlIV, 0)
%   Jd    = min(Na, AlVI)
%   CaTs  = max(AlVI - Jd, 0)
%   CaTi  = max((AlIV - CaTs)/2, 0), when AlIV > CaTs
%   DiHd  = max(Ca - CaTs - CaTi, 0)
%   EnFs  = max((Fe + Mg - DiHd)/2, 0)
%
% Equilibrium terms:
%
%   K_DiHd_Jd = (Jd_cpx * XCa_liq * Fm_liq) /
%                (DiHd_cpx * XNa_liq * XAl_liq)
%
%   K_DiHd_CaTs = (CaTs_cpx * XSi_liq * Fm_liq) /
%                  (DiHd_cpx * XAl_liq^2)
%
% Putirka et al. (1996), Table 5:
%
%   10000/T1 = 6.73 - 0.26*ln(K_DiHd_Jd)
%                    - 0.86*ln(MgNum_liq)
%                    + 0.52*ln(XCa_liq)
%
%   10000/T2 = 6.59 - 0.16*ln(K_DiHd_Jd)
%                    - 0.65*ln(MgNum_liq)
%                    + 0.23*ln(XCa_liq)
%                    - 0.02*P_kbar
%
%   10000/T3 = 6.92 - 0.18*ln(K_DiHd_CaTs)
%                    - 0.84*ln(MgNum_liq)
%                    - 0.29*ln(1/XAl_liq^2)
%
%   10000/T4 = 7.20 - 0.04*ln(K_DiHd_CaTs)
%                    - 0.59*ln(MgNum_liq)
%                    - 0.18*ln(1/XAl_liq^2)
%                    - 0.03*P_kbar
%
% Pressure is in kbar, temperature is in Kelvin, and ln is natural log.
% The positive sign of the T1 Ca term is retained exactly as printed in
% Table 5; this corrects the negative sign in the original script supplied
% for modification.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka1996(rawdata_struct, P_kbar)
%   results = Putirka1996(rawdata_struct, P_kbar, 'LiquidRow', n)
%
% Inputs:
%   rawdata_struct : struct containing the Cpx table described above
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Optional name-value input:
%   LiquidRow      : positive integer liquid-row index. When omitted, row 1
%                    is used and an fprintf message is printed if the loaded
%                    liquid dataset contains multiple rows.
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Cpx-liquid pair. NaN and Inf results are retained. Standard
%             T_K and T_deg columns contain the T2 result.
%

%% Input validation
if nargin < 2
    error('Putirka1996 requires (rawdata_struct, P_kbar).');
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

P_kbar = P_kbar(:);

%% Options
ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve Cpx and liquid datasets
disp('=== Step 1: Preparing Cpx and liquid datasets ===');

dataset_cpx = rawdata_struct.Cpx;
MWinfo = liquid.getMolarWeights();

[liqAll, metaLiq] = liquid.readLiquidExcel();
if isempty(liqAll) || ~istable(liqAll)
    error('Selected liquid dataset is empty or is not a table.');
end

if isempty(liquidRowOpt)
    idxLiq = 1;
    if height(liqAll) > 1
        fprintf(2, ...
            ['WARNING: The liquid dataset contains %d rows. Putirka1996 ' ...
             'uses row 1 because LiquidRow was not specified.\n'], ...
            height(liqAll));
    end
else
    idxLiq = liquidRowOpt;
    if idxLiq > height(liqAll)
        error('Requested LiquidRow (%d) exceeds liquid rows (%d).', ...
            idxLiq, height(liqAll));
    end
end

selectedData_liq = liqAll(idxLiq, :);
liqOxides = prepareLiquidOxides(selectedData_liq);

disp(['Liquid selected: Row ' num2str(idxLiq)]);
disp('=== Preparing Cpx and liquid datasets has been finished ===');

%% 2) Initialize output container
% A fixed-size buffer avoids all loop-dependent result-array resizing. If
% the limit is reached, completed results are returned without enlarging the
% array.
disp('=== Step 2: Preparing output container ===');

maxResultBlocks = max(1024, height(dataset_cpx));
resultBlocks = cell(maxResultBlocks, 1);
nResultBlocks = 0;

% Principal ranges used for non-stopping warnings.
principalP_min_kbar = 0;
principalP_max_kbar = 30;
strongP_max_kbar = 50;
principalT_min_degC = 1100;
principalT_max_degC = 1475;

pressureOutsidePrincipal = P_kbar < principalP_min_kbar | ...
    P_kbar > principalP_max_kbar;
pressureOutsideTested = P_kbar > strongP_max_kbar;
pressureWarningsIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive Cpx selection loop and calculation
disp('=== Step 3: Selecting a data code from the list (Cpx) ===');

while true
    % Prevent result-buffer resizing inside the loop.
    if nResultBlocks >= maxResultBlocks
        fprintf(2, ...
            ['WARNING: The fixed result-buffer limit of %d selections was ' ...
             'reached. Completed calculations will be returned without ' ...
             'enlarging the result array.\n'], maxResultBlocks);
        break;
    end

    dataCodes_cpx = dataset_cpx{:, 1};
    [idxCpx, ok] = listdlg( ...
        'PromptString', 'Please select the Cpx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(idxCpx)
        disp('Selection canceled');
        break;
    end

    codeCpx = dataCodes_cpx(idxCpx);
    disp(['Cpx selected: ' char(string(codeCpx))]);

    disp('=== Step 4: Calculating the temperature ===');

    selectedData_cpx = dataset_cpx(idxCpx, :);
    cpxOxides = prepareCpxOxides(selectedData_cpx);

    % F and Cl are intentionally absent from these checks because they are
    % not used in the liquid cation total or in T1-T4.
    nanInputNames = findNaNInputs(cpxOxides, liqOxides);
    validateInputValues(cpxOxides, liqOxides);

    row = calcTemp(cpxOxides, liqOxides, P_kbar, MWinfo);

    row.dataCode_cpx = repmat(string(codeCpx), height(row), 1);
    row.dataRow_liq = repmat(idxLiq, height(row), 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_cpx', 'dataRow_liq'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    resultBlocks{nResultBlocks} = row;

    % Echo all four model results.
    disp('--------------------------------------------------');
    disp('=== Temperatures were calculated: ===');
    displayTemperatureResult('T1', row.T1_C);
    displayTemperatureResult('T2', row.T2_C);
    displayTemperatureResult('T3', row.T3_C);
    displayTemperatureResult('T4', row.T4_C);

    % Pressure warnings are common to all selected Cpx rows and are printed
    % once per function call.
    if ~pressureWarningsIssued
        if any(pressureOutsidePrincipal)
            fprintf(2, ...
                ['WARNING: Input pressure is outside the principal 0-30 kbar ' ...
                 'range of Putirka et al. (1996). %d of %d pressure ' ...
                 'point(s) are outside; input range = %.6g-%.6g kbar. ' ...
                 'Values from 30-50 kbar are tested extrapolation rather ' ...
                 'than the principal calibration range.\n'], ...
                sum(pressureOutsidePrincipal), numel(P_kbar), ...
                min(P_kbar), max(P_kbar));
        end

        if any(pressureOutsideTested)
            fprintf(2, ...
                ['WARNING: Input pressure exceeds the approximately 50-kbar ' ...
                 'upper limit of the independent test data reported by ' ...
                 'Putirka et al. (1996). %d of %d pressure point(s) exceed ' ...
                 '50 kbar. This is strong pressure extrapolation.\n'], ...
                sum(pressureOutsideTested), numel(P_kbar));
        end
        pressureWarningsIssued = true;
    end

    printTemperatureRangeWarning(row.T1_C, 'T1', codeCpx, ...
        principalT_min_degC, principalT_max_degC);
    printTemperatureRangeWarning(row.T2_C, 'T2', codeCpx, ...
        principalT_min_degC, principalT_max_degC);
    printTemperatureRangeWarning(row.T3_C, 'T3', codeCpx, ...
        principalT_min_degC, principalT_max_degC);
    printTemperatureRangeWarning(row.T4_C, 'T4', codeCpx, ...
        principalT_min_degC, principalT_max_degC);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in Putirka et al. (1996) thermometer ' ...
             'input(s) for Cpx %s and liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; they were ' ...
             'not replaced by zero. F and Cl are excluded because they do ' ...
             'not enter the calculation.\n'], ...
            char(string(codeCpx)), idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    printNonFiniteWarning(row.T1_C, 'T1', codeCpx, idxLiq);
    printNonFiniteWarning(row.T2_C, 'T2', codeCpx, idxLiq);
    printNonFiniteWarning(row.T3_C, 'T3', codeCpx, idxLiq);
    printNonFiniteWarning(row.T4_C, 'T4', codeCpx, idxLiq);

    printDomainWarning(row.T1_domain_valid, 'T1', codeCpx, idxLiq, ...
        'positive K_DiHd_Jd, MgNum_liq, XCa_liq, and 10000/T1');
    printDomainWarning(row.T2_domain_valid, 'T2', codeCpx, idxLiq, ...
        'positive K_DiHd_Jd, MgNum_liq, XCa_liq, and 10000/T2');
    printDomainWarning(row.T3_domain_valid, 'T3', codeCpx, idxLiq, ...
        'positive K_DiHd_CaTs, MgNum_liq, XAl_liq, and 10000/T3');
    printDomainWarning(row.T4_domain_valid, 'T4', codeCpx, idxLiq, ...
        'positive K_DiHd_CaTs, MgNum_liq, XAl_liq, and 10000/T4');

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection (same Liquid dataset)?', ...
        'Putirka1996', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate result blocks once after all selections are complete.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct('liquid', metaLiq);
disp('=== Putirka1996 finished ===');

end

%% ---- local functions ----
function cpxOxides = prepareCpxOxides(data_cpx)
% prepareCpxOxides
% Read one selected Cpx oxide row. Missing optional columns are assigned the
% documented zero defaults. Present NaN values are preserved.

if height(data_cpx) ~= 1
    error('Cpx input must be a 1-row table.');
end

cpxOxides = struct();
cpxOxides.SiO2 = getOxideRequired(data_cpx, 'SiO2', 'Cpx');
cpxOxides.TiO2 = getOxideOptional(data_cpx, 'TiO2', 0);
cpxOxides.Al2O3 = getOxideRequired(data_cpx, 'Al2O3', 'Cpx');
cpxOxides.FeO = getFeORequired(data_cpx, 'Cpx');
cpxOxides.MnO = getOxideOptional(data_cpx, 'MnO', 0);
cpxOxides.MgO = getOxideRequired(data_cpx, 'MgO', 'Cpx');
cpxOxides.CaO = getOxideRequired(data_cpx, 'CaO', 'Cpx');
cpxOxides.Na2O = getOxideRequired(data_cpx, 'Na2O', 'Cpx');
cpxOxides.K2O = getOxideOptional(data_cpx, 'K2O', 0);
cpxOxides.Cr2O3 = getOxideOptional(data_cpx, 'Cr2O3', 0);

end

function liqOxides = prepareLiquidOxides(data_liq)
% prepareLiquidOxides
% Read one selected liquid oxide row. F and Cl are retained only as raw
% metadata and are excluded from calculation and validation.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liqOxides = struct();
liqOxides.SiO2 = getOxideRequired(data_liq, 'SiO2', 'Liquid');
liqOxides.TiO2 = getOxideOptional(data_liq, 'TiO2', 0);
liqOxides.Al2O3 = getOxideRequired(data_liq, 'Al2O3', 'Liquid');
liqOxides.FeO = getFeORequired(data_liq, 'Liquid');
liqOxides.MnO = getOxideOptional(data_liq, 'MnO', 0);
liqOxides.MgO = getOxideRequired(data_liq, 'MgO', 'Liquid');
liqOxides.CaO = getOxideRequired(data_liq, 'CaO', 'Liquid');
liqOxides.Na2O = getOxideRequired(data_liq, 'Na2O', 'Liquid');
liqOxides.K2O = getOxideOptional(data_liq, 'K2O', 0);
liqOxides.V2O3 = getOxideOptional(data_liq, 'V2O3', 0);
liqOxides.Cr2O3 = getOxideOptional(data_liq, 'Cr2O3', 0);
liqOxides.NiO = getOxideOptional(data_liq, 'NiO', 0);
liqOxides.P2O5 = getOxideOptional(data_liq, 'P2O5', 0);
liqOxides.SO3 = getOxideOptional(data_liq, 'SO3', 0);
liqOxides.Fe2O3 = getOxideOptional(data_liq, 'Fe2O3', 0);
liqOxides.F = getOxideOptional(data_liq, 'F', 0);
liqOxides.Cl = getOxideOptional(data_liq, 'Cl', 0);

end

function nanInputNames = findNaNInputs(cpxOxides, liqOxides)
% findNaNInputs
% Return names of NaN values actually used in Cpx normalization, liquid
% cation normalization, or T1-T4. F and Cl are intentionally excluded.

cpxFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};
liqFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'V2O3', 'Cr2O3', ...
    'NiO', 'P2O5', 'SO3', 'Fe2O3'};

nanBuffer = strings(numel(cpxFields) + numel(liqFields), 1);
nNan = 0;

for i = 1:numel(cpxFields)
    fieldName = cpxFields{i};
    if isnan(cpxOxides.(fieldName))
        nNan = nNan + 1;
        nanBuffer(nNan) = "Cpx." + string(fieldName);
    end
end

for i = 1:numel(liqFields)
    fieldName = liqFields{i};
    if isnan(liqOxides.(fieldName))
        nNan = nNan + 1;
        nanBuffer(nNan) = "Liquid." + string(fieldName);
    end
end

nanInputNames = nanBuffer(1:nNan);

end

function validateInputValues(cpxOxides, liqOxides)
% validateInputValues
% Reject finite negative values and Inf in all oxide values used by the
% calculation. Zero and NaN are allowed. F and Cl are excluded.

cpxFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};
liqFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'V2O3', 'Cr2O3', ...
    'NiO', 'P2O5', 'SO3', 'Fe2O3'};

invalidBuffer = strings(numel(cpxFields) + numel(liqFields), 1);
nInvalid = 0;

for i = 1:numel(cpxFields)
    fieldName = cpxFields{i};
    value = cpxOxides.(fieldName);
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Cpx." + string(fieldName);
    end
end

for i = 1:numel(liqFields)
    fieldName = liqFields{i};
    value = liqOxides.(fieldName);
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Liquid." + string(fieldName);
    end
end

if nInvalid > 0
    invalidNames = invalidBuffer(1:nInvalid);
    error(['Putirka1996: finite oxide inputs used by the thermometer ' ...
           'must be greater than or equal to zero, and Inf is not ' ...
           'permitted. Invalid value(s): ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(cpxOxides, liqOxides, P_kbar, MWinfo)
% calcTemp
% Calculate T1-T4 for one selected Cpx-liquid pair over a scalar or vector
% of pressures. One output table row is returned for each pressure value.

P_kbar = P_kbar(:);
nP = numel(P_kbar);
P_GPa = P_kbar ./ 10;

cpx = calcCpxComponents(cpxOxides, MWinfo);
liq = calcLiquidFractions(liqOxides, MWinfo);

% Liquid terms used by Table 5.
XFm_liq = liq.XFeO + liq.XMgO;
MgNum_liq = liq.XMgO ./ (liq.XMgO + liq.XFeO);
XNa_liq = liq.XNaO0_5;
XAl_liq = liq.XAlO1_5;
XCa_liq = liq.XCaO;
XSi_liq = liq.XSiO2;

K_DiHd_Jd = (cpx.XJd .* XCa_liq .* XFm_liq) ./ ...
    (cpx.XDiHd .* XNa_liq .* XAl_liq);
K_DiHd_CaTs = (cpx.XCaTs .* XSi_liq .* XFm_liq) ./ ...
    (cpx.XDiHd .* (XAl_liq .^ 2));

% T1: pressure-independent DiHd-Jd thermometer.
T1_baseValid = isPositiveFinite([K_DiHd_Jd, MgNum_liq, XCa_liq]);
invT1_scalar = NaN;
if T1_baseValid
    invT1_scalar = 6.73 ...
        - 0.26 .* log(K_DiHd_Jd) ...
        - 0.86 .* log(MgNum_liq) ...
        + 0.52 .* log(XCa_liq);
end
T1_domain_scalar = T1_baseValid && isfinite(invT1_scalar) && ...
    invT1_scalar > 0;
invT1 = repmat(invT1_scalar, nP, 1);
T1_domain_valid = repmat(T1_domain_scalar, nP, 1);
[T1_K, T1_C] = invTtoT(invT1, T1_domain_valid);

% T2: pressure-dependent DiHd-Jd thermometer.
T2_baseValid = T1_baseValid;
invT2 = NaN(nP, 1);
if T2_baseValid
    invT2 = 6.59 ...
        - 0.16 .* log(K_DiHd_Jd) ...
        - 0.65 .* log(MgNum_liq) ...
        + 0.23 .* log(XCa_liq) ...
        - 0.02 .* P_kbar;
end
T2_domain_valid = T2_baseValid & isfinite(invT2) & invT2 > 0;
[T2_K, T2_C] = invTtoT(invT2, T2_domain_valid);

% T3: pressure-independent DiHd-CaTs thermometer.
T3_baseValid = isPositiveFinite( ...
    [K_DiHd_CaTs, MgNum_liq, XAl_liq]);
invT3_scalar = NaN;
if T3_baseValid
    invT3_scalar = 6.92 ...
        - 0.18 .* log(K_DiHd_CaTs) ...
        - 0.84 .* log(MgNum_liq) ...
        - 0.29 .* log(1 ./ (XAl_liq .^ 2));
end
T3_domain_scalar = T3_baseValid && isfinite(invT3_scalar) && ...
    invT3_scalar > 0;
invT3 = repmat(invT3_scalar, nP, 1);
T3_domain_valid = repmat(T3_domain_scalar, nP, 1);
[T3_K, T3_C] = invTtoT(invT3, T3_domain_valid);

% T4: pressure-dependent DiHd-CaTs thermometer.
T4_baseValid = T3_baseValid;
invT4 = NaN(nP, 1);
if T4_baseValid
    invT4 = 7.20 ...
        - 0.04 .* log(K_DiHd_CaTs) ...
        - 0.59 .* log(MgNum_liq) ...
        - 0.18 .* log(1 ./ (XAl_liq .^ 2)) ...
        - 0.03 .* P_kbar;
end
T4_domain_valid = T4_baseValid & isfinite(invT4) & invT4 > 0;
[T4_K, T4_C] = invTtoT(invT4, T4_domain_valid);

% Pack outputs.
row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% Raw Cpx oxide values.
row.SiO2_cpx = repmat(cpxOxides.SiO2, nP, 1);
row.TiO2_cpx = repmat(cpxOxides.TiO2, nP, 1);
row.Al2O3_cpx = repmat(cpxOxides.Al2O3, nP, 1);
row.FeO_cpx = repmat(cpxOxides.FeO, nP, 1);
row.MnO_cpx = repmat(cpxOxides.MnO, nP, 1);
row.MgO_cpx = repmat(cpxOxides.MgO, nP, 1);
row.CaO_cpx = repmat(cpxOxides.CaO, nP, 1);
row.Na2O_cpx = repmat(cpxOxides.Na2O, nP, 1);
row.K2O_cpx = repmat(cpxOxides.K2O, nP, 1);
row.Cr2O3_cpx = repmat(cpxOxides.Cr2O3, nP, 1);

% Cpx cations and components.
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
row.oxygenSum_cpx = repmat(cpx.oxygenSum, nP, 1);
row.XAlIV_cpx = repmat(cpx.XAlIV, nP, 1);
row.XAlVI_cpx = repmat(cpx.XAlVI, nP, 1);
row.XJd_cpx = repmat(cpx.XJd, nP, 1);
row.XCaTs_cpx = repmat(cpx.XCaTs, nP, 1);
row.XCaTi_cpx = repmat(cpx.XCaTi, nP, 1);
row.XDiHd_cpx = repmat(cpx.XDiHd, nP, 1);
row.XEnFs_cpx = repmat(cpx.XEnFs, nP, 1);

% Raw liquid oxide values. F and Cl are retained only for traceability.
row.SiO2_liq = repmat(liqOxides.SiO2, nP, 1);
row.TiO2_liq = repmat(liqOxides.TiO2, nP, 1);
row.Al2O3_liq = repmat(liqOxides.Al2O3, nP, 1);
row.FeO_liq = repmat(liqOxides.FeO, nP, 1);
row.MnO_liq = repmat(liqOxides.MnO, nP, 1);
row.MgO_liq = repmat(liqOxides.MgO, nP, 1);
row.CaO_liq = repmat(liqOxides.CaO, nP, 1);
row.Na2O_liq = repmat(liqOxides.Na2O, nP, 1);
row.K2O_liq = repmat(liqOxides.K2O, nP, 1);
row.V2O3_liq = repmat(liqOxides.V2O3, nP, 1);
row.Cr2O3_liq = repmat(liqOxides.Cr2O3, nP, 1);
row.NiO_liq = repmat(liqOxides.NiO, nP, 1);
row.P2O5_liq = repmat(liqOxides.P2O5, nP, 1);
row.SO3_liq = repmat(liqOxides.SO3, nP, 1);
row.Fe2O3_liq = repmat(liqOxides.Fe2O3, nP, 1);
row.F_liq = repmat(liqOxides.F, nP, 1);
row.Cl_liq = repmat(liqOxides.Cl, nP, 1);

% Liquid cation totals and fractions.
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
row.XFm_liq = repmat(XFm_liq, nP, 1);
row.MgNum_liq = repmat(MgNum_liq, nP, 1);

% Equilibrium terms and inverse-temperature expressions.
row.K_DiHd_Jd = repmat(K_DiHd_Jd, nP, 1);
row.K_DiHd_CaTs = repmat(K_DiHd_CaTs, nP, 1);
row.invT1_10000_over_K = invT1;
row.invT2_10000_over_K = invT2;
row.invT3_10000_over_K = invT3;
row.invT4_10000_over_K = invT4;
row.T1_domain_valid = T1_domain_valid;
row.T2_domain_valid = T2_domain_valid;
row.T3_domain_valid = T3_domain_valid;
row.T4_domain_valid = T4_domain_valid;

% Original model outputs.
row.T1_K = T1_K;
row.T1_C = T1_C;
row.T2_K = T2_K;
row.T2_C = T2_C;
row.T3_K = T3_K;
row.T3_C = T3_C;
row.T4_K = T4_K;
row.T4_C = T4_C;

% Standard thermoCalcMin output uses the preferred pressure-dependent T2.
row.T_K = T2_K;
row.T_deg = T2_C;

end

function cpx = calcCpxComponents(ox, MWinfo)
% calcCpxComponents
% Convert Cpx oxide wt.% to 6-oxygen cations and allocate the components
% used by Putirka et al. (1996), preserving NaN values.

mol = struct();
mol.SiO2 = ox.SiO2 ./ MWinfo.MW.SiO2;
mol.TiO2 = ox.TiO2 ./ MWinfo.MW.TiO2;
mol.Al2O3 = ox.Al2O3 ./ MWinfo.MW.Al2O3;
mol.FeO = ox.FeO ./ MWinfo.MW.FeO;
mol.MnO = ox.MnO ./ MWinfo.MW.MnO;
mol.MgO = ox.MgO ./ MWinfo.MW.MgO;
mol.CaO = ox.CaO ./ MWinfo.MW.CaO;
mol.Na2O = ox.Na2O ./ MWinfo.MW.Na2O;
mol.K2O = ox.K2O ./ MWinfo.MW.K2O;
mol.Cr2O3 = ox.Cr2O3 ./ MWinfo.MW.Cr2O3;

oxygenSum = 2 .* mol.SiO2 + 2 .* mol.TiO2 + 3 .* mol.Al2O3 + ...
    mol.FeO + mol.MnO + mol.MgO + mol.CaO + mol.Na2O + ...
    mol.K2O + 3 .* mol.Cr2O3;

if isfinite(oxygenSum) && oxygenSum > 0
    ORF = 6 ./ oxygenSum;
else
    ORF = NaN;
end

cpx = struct();
cpx.oxygenSum = oxygenSum;
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

cpx.cationSum = cpx.XSi + cpx.XTi + cpx.XAl + cpx.XFe + ...
    cpx.XMn + cpx.XMg + cpx.XCa + cpx.XNa + cpx.XK + cpx.XCr;

if isnan(cpx.XSi)
    cpx.XAlIV = NaN;
else
    cpx.XAlIV = max(2 - cpx.XSi, 0);
end

if isnan(cpx.XAl) || isnan(cpx.XAlIV)
    cpx.XAlVI = NaN;
else
    cpx.XAlVI = max(cpx.XAl - cpx.XAlIV, 0);
end

if isnan(cpx.XNa) || isnan(cpx.XAlVI)
    cpx.XJd = NaN;
else
    cpx.XJd = max(min(cpx.XNa, cpx.XAlVI), 0);
end

if isnan(cpx.XAlVI) || isnan(cpx.XJd)
    cpx.XCaTs = NaN;
else
    cpx.XCaTs = max(cpx.XAlVI - cpx.XJd, 0);
end

if isnan(cpx.XAlIV) || isnan(cpx.XCaTs)
    cpx.XCaTi = NaN;
elseif cpx.XAlIV > cpx.XCaTs
    cpx.XCaTi = max((cpx.XAlIV - cpx.XCaTs) ./ 2, 0);
else
    cpx.XCaTi = 0;
end

if isnan(cpx.XCa) || isnan(cpx.XCaTs) || isnan(cpx.XCaTi)
    cpx.XDiHd = NaN;
else
    cpx.XDiHd = max(cpx.XCa - cpx.XCaTs - cpx.XCaTi, 0);
end

if isnan(cpx.XFe) || isnan(cpx.XMg) || isnan(cpx.XDiHd)
    cpx.XEnFs = NaN;
else
    cpx.XEnFs = max((cpx.XFe + cpx.XMg - cpx.XDiHd) ./ 2, 0);
end

end

function liq = calcLiquidFractions(ox, MWinfo)
% calcLiquidFractions
% Convert liquid oxide wt.% to cation fractions. F and Cl are deliberately
% omitted because they are anions and do not belong in cationTotal_liq.

n = struct();
n.SiO2 = ox.SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
n.TiO2 = ox.TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
n.Al2O3 = ox.Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n.FeO = ox.FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
n.MnO = ox.MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
n.MgO = ox.MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
n.CaO = ox.CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
n.Na2O = ox.Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
n.K2O = ox.K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
n.V2O3 = ox.V2O3 .* MWinfo.Cat.V2O3 ./ MWinfo.MW.V2O3;
n.Cr2O3 = ox.Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
n.NiO = ox.NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;
n.P2O5 = ox.P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;
n.SO3 = ox.SO3 .* MWinfo.Cat.SO3 ./ MWinfo.MW.SO3;
n.Fe2O3 = ox.Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

cationTotal = n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + ...
    n.MgO + n.CaO + n.Na2O + n.K2O + n.V2O3 + n.Cr2O3 + ...
    n.NiO + n.P2O5 + n.SO3 + n.Fe2O3;

if isfinite(cationTotal) && cationTotal > 0
    divisor = cationTotal;
else
    divisor = NaN;
end

liq = struct();
liq.cationTotal = cationTotal;
liq.XSiO2 = n.SiO2 ./ divisor;
liq.XTiO2 = n.TiO2 ./ divisor;
liq.XAlO1_5 = n.Al2O3 ./ divisor;
liq.XFeO = n.FeO ./ divisor;
liq.XMnO = n.MnO ./ divisor;
liq.XMgO = n.MgO ./ divisor;
liq.XCaO = n.CaO ./ divisor;
liq.XNaO0_5 = n.Na2O ./ divisor;
liq.XKO0_5 = n.K2O ./ divisor;

end

function [T_K, T_C] = invTtoT(invT, domainValid)
% invTtoT
% Convert 10000/T values to Kelvin and degreeC only where the equation
% domain is valid. Invalid values remain NaN.

T_K = NaN(size(invT));
T_C = NaN(size(invT));
valid = domainValid & isfinite(invT) & invT > 0;
T_K(valid) = 10000 ./ invT(valid);
T_C(valid) = T_K(valid) - 273.15;

end

function tf = isPositiveFinite(values)
% isPositiveFinite
% True only when every supplied scalar is finite and strictly positive.

tf = all(isfinite(values)) && all(values > 0);

end

function value = getOxideRequired(tbl, oxide, phaseLabel)
% getOxideRequired
% Read one required oxide scalar. A present NaN is preserved.

columnName = findOxideColumn(tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    error('%s table must contain oxide variable: %s', phaseLabel, oxide);
end

value = toScalarDoublePreserveNaN(tbl.(columnName), columnName);

end

function value = getOxideOptional(tbl, oxide, defaultValue)
% getOxideOptional
% Return defaultValue only when the optional column is absent. Present NaN
% is preserved and never replaced by zero.

columnName = findOxideColumn(tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    value = defaultValue;
else
    value = toScalarDoublePreserveNaN(tbl.(columnName), columnName);
end

end

function value = getFeORequired(tbl, phaseLabel)
% getFeORequired
% Use FeO when its column exists, preserving a present NaN. Use FeOt only
% when the FeO column is absent. This prevents a NaN FeO value from being
% silently replaced by FeOt.

feOColumn = findOxideColumn(tbl.Properties.VariableNames, 'FeO');
if ~isempty(feOColumn)
    value = toScalarDoublePreserveNaN(tbl.(feOColumn), feOColumn);
    return;
end

feOtColumn = findOxideColumn(tbl.Properties.VariableNames, 'FeOt');
if ~isempty(feOtColumn)
    value = toScalarDoublePreserveNaN(tbl.(feOtColumn), feOtColumn);
    return;
end

error('%s table must contain FeO or FeOt.', phaseLabel);

end

function columnName = findOxideColumn(varNames, oxide)
% findOxideColumn
% Match oxide names such as SiO2 or SiO2Value while ignoring spaces,
% underscores, and hyphens.

canonicalNames = cell(size(varNames));
for i = 1:numel(varNames)
    canonicalNames{i} = canonicalizeName(varNames{i});
end

oxideCanonical = canonicalizeName(oxide);
targets = {[oxideCanonical 'value'], oxideCanonical};

columnName = '';
for i = 1:numel(targets)
    idx = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(idx)
        columnName = varNames{idx};
        return;
    end
end

end

function canonical = canonicalizeName(name)
% canonicalizeName
% Produce a case-insensitive oxide-column key.

canonical = lower(char(string(name)));
canonical = strrep(canonical, ' ', '');
canonical = strrep(canonical, '_', '');
canonical = strrep(canonical, '-', '');

end

function value = toScalarDoublePreserveNaN(rawValue, variableName)
% toScalarDoublePreserveNaN
% Convert one table value to a real numeric scalar. Missing or unparsable
% values become NaN and are not converted to zero.

if isempty(rawValue)
    value = NaN;
    return;
end

if isnumeric(rawValue)
    if numel(rawValue) ~= 1 || ~isreal(rawValue)
        error('Variable %s must contain one real numeric value.', variableName);
    end
    value = double(rawValue);
    return;
end

if islogical(rawValue)
    if numel(rawValue) ~= 1
        error('Variable %s must contain one scalar value.', variableName);
    end
    value = double(rawValue);
    return;
end

if iscell(rawValue)
    if numel(rawValue) ~= 1
        error('Variable %s must contain one scalar value.', variableName);
    end
    value = toScalarDoublePreserveNaN(rawValue{1}, variableName);
    return;
end

if isstring(rawValue)
    if numel(rawValue) ~= 1 || ismissing(rawValue)
        value = NaN;
        return;
    end
    value = str2double(rawValue);
    return;
end

if ischar(rawValue)
    value = str2double(string(rawValue));
    return;
end

if iscategorical(rawValue)
    if numel(rawValue) ~= 1 || isundefined(rawValue)
        value = NaN;
        return;
    end
    value = str2double(string(rawValue));
    return;
end

error('Variable %s must be numeric or convertible to a numeric scalar.', ...
    variableName);

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Replicate common liquid identifiers to match the pressure-vector height.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    row.liq_Index = repeatTableValue(data_liq.('Index'), nRows);
end
if any(strcmp(variableNames, 'Experiment'))
    row.liq_Experiment = repmat(string(data_liq.('Experiment')), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    row.liq_Citation = repmat(string(data_liq.('Citation')), nRows, 1);
end

end

function repeated = repeatTableValue(value, nRows)
% repeatTableValue
% Repeat one scalar table value without changing its basic MATLAB type.

if isnumeric(value) || islogical(value) || isstring(value) || ...
        iscategorical(value) || isdatetime(value) || isduration(value)
    repeated = repmat(value, nRows, 1);
elseif iscell(value)
    repeated = repmat(value, nRows, 1);
elseif ischar(value)
    repeated = repmat(string(value), nRows, 1);
else
    repeated = repmat(string(value), nRows, 1);
end

end

function displayTemperatureResult(modelName, values)
% displayTemperatureResult
% Print one value or the first-to-last range for immediate feedback.

if isscalar(values)
    disp([modelName ': ' num2str(values) ' degreeC']);
else
    disp([modelName ': ' num2str(values(1)) ' to ' ...
        num2str(values(end)) ' degreeC']);
end

end

function printTemperatureRangeWarning(values, modelName, codeCpx, ...
        minimumTemperature, maximumTemperature)
% printTemperatureRangeWarning
% Warn for finite temperatures outside the principal new-experiment range.

finiteMask = isfinite(values);
outsideMask = finiteMask & ...
    (values < minimumTemperature | values > maximumTemperature);

if any(outsideMask)
    finiteValues = values(finiteMask);
    fprintf(2, ...
        ['WARNING: %s temperature is outside the principal new-experiment ' ...
         'range of Putirka et al. (1996): 1100-1475 degreeC. %d of %d ' ...
         'finite point(s) are outside; calculated finite range = ' ...
         '%.6g-%.6g degreeC for Cpx %s. This interval is not a strict ' ...
         'rectangular limit for every datum in the final regression.\n'], ...
        modelName, sum(outsideMask), sum(finiteMask), ...
        min(finiteValues), max(finiteValues), char(string(codeCpx)));
end

end

function printNonFiniteWarning(values, modelName, codeCpx, idxLiq)
% printNonFiniteWarning
% Retain and report NaN or Inf temperature outputs.

invalidMask = ~isfinite(values);
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite %s temperature values were calculated for ' ...
         'Cpx %s and liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table and the ' ...
         'calculation has not been stopped.\n'], ...
        modelName, char(string(codeCpx)), idxLiq, ...
        sum(invalidMask), numel(values), sum(isnan(values)), ...
        sum(isinf(values)));
end

end

function printDomainWarning(domainValid, modelName, codeCpx, idxLiq, ...
        domainDescription)
% printDomainWarning
% Report pressure rows for which a model is mathematically undefined.

invalidMask = ~domainValid;
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Putirka et al. (1996) %s is outside its mathematical ' ...
         'domain for Cpx %s and liquid row %d at %d of %d pressure ' ...
         'point(s). The model requires %s. Corresponding temperatures ' ...
         'remain NaN.\n'], ...
        modelName, char(string(codeCpx)), idxLiq, ...
        sum(invalidMask), numel(invalidMask), domainDescription);
end

end
