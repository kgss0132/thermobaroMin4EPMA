function results = Putirka2003(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Pyroxene/Putirka2003.m
% Tested with MATLAB R2024b
%
% Clinopyroxene-liquid thermometer, Table 4 Model B
% Putirka, K.D., Mikaelian, H., Ryerson, F. and Shaw, H. (2003)
% American Mineralogist, 88, 1542-1554
% DOI: https://doi.org/10.2138/am-2003-1017
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Clinopyroxene (Cpx) analysis and
% pairs it with one liquid composition loaded by liquid.readLiquidExcel. It
% calculates temperature using Model B in Table 4 of Putirka et al. (2003).
%
% The function accepts either a scalar pressure or a pressure vector and is
% therefore compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP. For every selected Cpx-liquid pair, the output
% contains one row per input pressure value.
%
% The function is designed for repeated calculations. Each result block is
% stored in a fixed-size preallocated cell buffer and all blocks are
% concatenated once after the interactive loop. No result array is enlarged
% inside the loop.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka et al. (2003) developed Model B to extend Cpx-liquid
% thermobarometry from mainly basaltic liquids to mafic, evolved,
% high-SiO2, and nominally volatile-bearing compositions.
%
% The new experiments reported in the paper cover approximately:
%
%   Pressure    : 10-35 kbar
%   Temperature : 850-1300 degreeC
%   Liquid SiO2 : up to 71.3 wt.%
%
% These values describe the direct new-experiment envelope in Tables 1-3,
% not strict rectangular limits for every datum used in the final Model B
% regression. The regression also includes selected data from earlier
% studies and therefore covers a broader P-T-composition space.
%
% Model B performance reported in the paper:
%
%   Calibration data:
%     n    = 94
%     R2   = approximately 0.96-0.97
%     SEE  = approximately 33-34 K
%
%   Independent test data:
%     n    = 346
%     R2   = 0.84
%     SEE  = 63 K
%
% Relevant locations in the original paper:
%
%   p. 1542      : abstract, model accuracy, evolved and volatile-bearing use
%   pp. 1544-1546: liquid fractions, Cpx components, and Table 4 Model B
%   pp. 1547-1549: calibration/test data and model performance
%   pp. 1549-1552: Fe-Mg equilibrium checks and natural-sample application
%
% Important application cautions:
%
%   1) The selected Cpx and liquid must represent an equilibrium pair.
%      Experimental Cpx rim compositions in direct contact with homogeneous
%      glass were used preferentially. Unrelated cores, rims, xenocrysts,
%      antecrysts, mixed magmas, altered glass, and reaction zones may return
%      geologically meaningless temperatures (pp. 1547-1548).
%
%   2) A whole-rock composition is not automatically an equilibrium liquid.
%      It is appropriate only when it reasonably represents the melt from
%      which the selected Cpx crystallized. Crystal accumulation,
%      fractionation, assimilation, and magma mixing can invalidate a
%      Cpx-whole-rock pair (pp. 1549-1552).
%
%   3) Putirka et al. (2003) recommend independent equilibrium screening.
%      Their Cpx-liquid Fe-Mg exchange data have mean KD = 0.275 and standard
%      deviation = 0.067; the approximately 3-sigma-filtered experimental
%      range is 0.105-0.488 (p. 1550 and Table 4, Models C-D). This function
%      reports KD and prints a non-stopping warning outside 0.105-0.488.
%
%   4) Liquid components are cation fractions of SiO2, TiO2, AlO1.5, FeO,
%      MnO, MgO, CaO, NaO0.5, KO0.5, and CrO1.5 (pp. 1545-1546). H2O is not
%      included, and oxide wt.% values are not first renormalized to 100.
%      V, Ni, P, S, Fe2O3, F, and Cl are retained as raw output values when
%      available but are not included in the published liquid cation total.
%
%   5) Cpx cations are normalized to 6 oxygens and should sum to about 4.
%      Jd is the lesser of Na and AlVI. DiHd is the Ca remaining after
%      allocating CaTs, CaTi, and CrCaTs. The 2003 scheme differs from the
%      1996 scheme by subtracting Ca used in CrCaTs from DiHd (p. 1546).
%
%   6) Acmite is not explicitly allocated because calculated Fe3+ is
%      uncertain. The charge-balance Fe3+ estimate is retained only as a
%      rough stoichiometric diagnostic and is not used directly in Model B
%      (p. 1547).
%
%   7) The logarithmic equation requires positive Jd, DiHd, liquid Na, Al,
%      Ca, Si, Fe+Mg, and Mg'. Zero or negative values make the model
%      mathematically undefined. In this implementation, such cases produce
%      NaN temperatures and non-stopping fprintf warnings.
%
% This implementation issues non-stopping fprintf warnings when:
%
%   1) input pressure is outside the direct 10-35 kbar experiment envelope;
%   2) input pressure exceeds the approximately 110-kbar outer limit of the
%      full independent experimental compilation;
%   3) a finite calculated temperature is outside 850-1300 degreeC;
%   4) any value actually used in the calculation is NaN;
%   5) the Model B equation is outside its mathematical domain;
%   6) the calculated temperature is NaN or Inf; or
%   7) finite Cpx-liquid Fe-Mg KD lies outside 0.105-0.488.
%
% The P-T warnings identify extrapolation relative to the direct new
% experiments; they do not imply that every value outside this rectangle is
% automatically invalid because the final regression includes additional
% experiments.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST Cpx-table column is treated as an identifier ("data code") for
% the selection dialog. Cpx analyses are expected as oxide wt.% columns.
% Column names may be oxide names such as SiO2 or names such as SiO2Value.
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
%   TiO2, MnO, K2O, Cr2O3
%
% Additional liquid values retained only as raw output metadata:
%   V2O3, NiO, P2O5, SO3, Fe2O3, F, Cl, H2O
%
% If an oxide column is present and its selected value is NaN, the NaN is
% retained. It is never replaced by zero. If both FeO and FeOt columns exist,
% FeO is used; a present NaN FeO value is retained and is not replaced by
% FeOt. FeOt is used only when the FeO column is absent.
%
% F and Cl are anions. They are retained as raw output values for
% traceability but are excluded from cationTotal_liq, temperature calculation,
% NaN-input warnings, and negative-value validation. H2O is likewise excluded
% from anhydrous cation-fraction normalization. V, Ni, P, S, and Fe2O3 are
% also excluded because they are not part of the liquid-component list used
% by Putirka et al. (2003) Model B.
%
% Every finite oxide value used in the Cpx or liquid calculation must be
% greater than or equal to zero. Finite negative values and Inf are rejected.
% Zero and NaN are allowed as raw values; invalid equation domains return
% NaN temperatures and non-stopping warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Liquid oxides are converted to the published anhydrous cation fractions.
% Define:
%
%   Fm_liq       = XFeO_liq + XMgO_liq
%   MgPrime_liq  = XMgO_liq / (XMgO_liq + XFeO_liq)
%
% Cpx components are calculated on a 6-oxygen basis:
%
%   AlIV    = max(2 - Si, 0)
%   AlVI    = max(Al_total - AlIV, 0)
%   Jd      = min(Na, AlVI)
%   CaTs    = max(AlVI - Jd, 0)
%   CaTi    = max((AlIV - CaTs)/2, 0), when AlIV > CaTs
%   CrCaTs  = Cr/2
%   DiHd    = max(Ca - CaTi - CaTs - CrCaTs, 0)
%
% Exchange term:
%
%   K_DiHd_Jd = (Jd_cpx * XCa_liq * XFm_liq) /
%                (DiHd_cpx * XNa_liq * XAl_liq)
%
% Putirka et al. (2003), Table 4 Model B:
%
%   10000/T = 4.60
%             - 0.437*ln(K_DiHd_Jd)
%             - 0.654*ln(MgPrime_liq)
%             - 0.326*ln(XNa_liq)
%             - 0.00632*P_kbar
%             - 0.920*ln(XSi_liq)
%             + 0.274*ln(Jd_cpx)
%
% Pressure is in kbar, temperature is in Kelvin, and ln is natural log.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2003(rawdata_struct, P_kbar)
%   results = Putirka2003(rawdata_struct, P_kbar, 'LiquidRow', n)
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
%             Cpx-liquid pair. NaN and Inf results are retained. T2003_K and
%             T2003_C are retained, and standard T_K and T_deg columns contain
%             the same Model B result.
%

%% Input validation
if nargin < 2
    error('Putirka2003 requires (rawdata_struct, P_kbar).');
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
            ['WARNING: The liquid dataset contains %d rows. Putirka2003 ' ...
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

% Direct new-experiment ranges used for non-stopping warnings.
directP_min_kbar = 10;
directP_max_kbar = 35;
outerTestP_max_kbar = 110;
directT_min_degC = 850;
directT_max_degC = 1300;

pressureOutsideDirect = P_kbar < directP_min_kbar | ...
    P_kbar > directP_max_kbar;
pressureOutsideOuterTest = P_kbar > outerTestP_max_kbar;
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

    % Only variables that enter Cpx normalization, published liquid
    % normalization, or Model B are included. F and Cl are excluded.
    nanInputNames = findNaNInputs(cpxOxides, liqOxides);
    validateInputValues(cpxOxides, liqOxides);

    row = calcTemp(cpxOxides, liqOxides, P_kbar, MWinfo);

    row.dataCode_cpx = repmat(string(codeCpx), height(row), 1);
    row.dataRow_liq = repmat(idxLiq, height(row), 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_cpx', 'dataRow_liq'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    resultBlocks{nResultBlocks} = row;

    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    displayTemperatureResult(row.T2003_C);

    % Pressure warnings are common to all selected Cpx rows and are printed
    % once per function call.
    if ~pressureWarningsIssued
        if any(pressureOutsideDirect)
            fprintf(2, ...
                ['WARNING: Input pressure is outside the direct new-experiment ' ...
                 'range of Putirka et al. (2003): 10-35 kbar. %d of %d ' ...
                 'pressure point(s) are outside; input range = %.6g-%.6g ' ...
                 'kbar. These bounds are not strict rectangular limits for ' ...
                 'the full Model B regression.\n'], ...
                sum(pressureOutsideDirect), numel(P_kbar), ...
                min(P_kbar), max(P_kbar));
        end

        if any(pressureOutsideOuterTest)
            fprintf(2, ...
                ['WARNING: Input pressure exceeds the approximately 110-kbar ' ...
                 'outer limit of the experimental compilation tested by ' ...
                 'Putirka et al. (2003). %d of %d pressure point(s) exceed ' ...
                 '110 kbar. This is strong pressure extrapolation.\n'], ...
                sum(pressureOutsideOuterTest), numel(P_kbar));
        end
        pressureWarningsIssued = true;
    end

    printTemperatureRangeWarning(row.T2003_C, codeCpx, ...
        directT_min_degC, directT_max_degC);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in Putirka et al. (2003) thermometer ' ...
             'input(s) for Cpx %s and liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; they were ' ...
             'not replaced by zero. F and Cl are excluded because they do ' ...
             'not enter cationTotal_liq or Model B.\n'], ...
            char(string(codeCpx)), idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    printNonFiniteWarning(row.T2003_C, codeCpx, idxLiq);
    printDomainWarning(row.model_domain_valid, codeCpx, idxLiq);
    printKdWarning(row.KD_FeMg_cpx_liq, codeCpx, idxLiq);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection (same Liquid dataset)?', ...
        'Putirka2003', ...
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
disp('=== Putirka2003 finished ===');

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
% Read one selected liquid oxide row. Only the published Model B component
% list enters normalization. Other values are retained only as raw output.

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
liqOxides.Cr2O3 = getOxideOptional(data_liq, 'Cr2O3', 0);

% Raw metadata excluded from published cation normalization and Model B.
liqOxides.V2O3 = getOxideOptional(data_liq, 'V2O3', 0);
liqOxides.NiO = getOxideOptional(data_liq, 'NiO', 0);
liqOxides.P2O5 = getOxideOptional(data_liq, 'P2O5', 0);
liqOxides.SO3 = getOxideOptional(data_liq, 'SO3', 0);
liqOxides.Fe2O3 = getOxideOptional(data_liq, 'Fe2O3', 0);
liqOxides.F = getOxideOptional(data_liq, 'F', 0);
liqOxides.Cl = getOxideOptional(data_liq, 'Cl', 0);
liqOxides.H2O = getOxideOptional(data_liq, 'H2O', 0);

end

function nanInputNames = findNaNInputs(cpxOxides, liqOxides)
% findNaNInputs
% Return names of NaN values actually used in Cpx normalization, published
% liquid normalization, or Model B. F and Cl are intentionally excluded.

cpxFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};
liqFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};

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
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};

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
    error(['Putirka2003: finite oxide inputs used by the thermometer ' ...
           'must be greater than or equal to zero, and Inf is not ' ...
           'permitted. Invalid value(s): ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(cpxOxides, liqOxides, P_kbar, MWinfo)
% calcTemp
% Calculate Table 4 Model B for one selected Cpx-liquid pair over a scalar
% or vector of pressures. One output table row is returned per pressure.

P_kbar = P_kbar(:);
nP = numel(P_kbar);
P_GPa = P_kbar ./ 10;

cpx = calcCpxComponents(cpxOxides, MWinfo);
liq = calcLiquidFractions(liqOxides, MWinfo);

% Liquid terms used by Model B.
XNa_liq = liq.XNaO0_5;
XAl_liq = liq.XAlO1_5;
XCa_liq = liq.XCaO;
XSi_liq = liq.XSiO2;
XFm_liq = liq.XFeO + liq.XMgO;
MgPrime_liq = liq.XMgO ./ (liq.XMgO + liq.XFeO);

exchangeTerm = (cpx.XJd .* XCa_liq .* XFm_liq) ./ ...
    (cpx.XDiHd .* XNa_liq .* XAl_liq);

% Fe-Mg equilibrium coefficient as defined in Putirka et al. (2003).
KD_FeMg = (liq.XMgO .* cpx.XFe) ./ (cpx.XMg .* liq.XFeO);

% Base logarithmic domain is independent of pressure.
baseDomainValid = isPositiveFinite([cpx.XJd, cpx.XDiHd, ...
    XCa_liq, XNa_liq, XAl_liq, XSi_liq, XFm_liq, ...
    MgPrime_liq, exchangeTerm]);

invT = NaN(nP, 1);
if baseDomainValid
    invT = 4.60 ...
        - 4.37e-1 .* log(exchangeTerm) ...
        - 6.54e-1 .* log(MgPrime_liq) ...
        - 3.26e-1 .* log(XNa_liq) ...
        - 6.32e-3 .* P_kbar ...
        - 9.20e-1 .* log(XSi_liq) ...
        + 2.74e-1 .* log(cpx.XJd);
end

modelDomainValid = baseDomainValid & isfinite(invT) & invT > 0;
[T2003_K, T2003_C] = invTtoT(invT, modelDomainValid);

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
row.XFe3_chargeBalance_cpx = repmat(cpx.XFe3_chargeBalance, nP, 1);
row.XFe3_cpx = repmat(cpx.XFe3, nP, 1);
row.XJd_cpx = repmat(cpx.XJd, nP, 1);
row.XCaTs_cpx = repmat(cpx.XCaTs, nP, 1);
row.XCaTi_cpx = repmat(cpx.XCaTi, nP, 1);
row.XCrCaTs_cpx = repmat(cpx.XCrCaTs, nP, 1);
row.XDiHd_cpx = repmat(cpx.XDiHd, nP, 1);
row.XEnFs_cpx = repmat(cpx.XEnFs, nP, 1);
row.XFmCaTs_cpx = repmat(cpx.XFmCaTs, nP, 1);
row.XFmTi_cpx = repmat(cpx.XFmTi, nP, 1);

% Raw liquid oxide values. Excluded values are retained for traceability.
row.SiO2_liq = repmat(liqOxides.SiO2, nP, 1);
row.TiO2_liq = repmat(liqOxides.TiO2, nP, 1);
row.Al2O3_liq = repmat(liqOxides.Al2O3, nP, 1);
row.FeO_liq = repmat(liqOxides.FeO, nP, 1);
row.MnO_liq = repmat(liqOxides.MnO, nP, 1);
row.MgO_liq = repmat(liqOxides.MgO, nP, 1);
row.CaO_liq = repmat(liqOxides.CaO, nP, 1);
row.Na2O_liq = repmat(liqOxides.Na2O, nP, 1);
row.K2O_liq = repmat(liqOxides.K2O, nP, 1);
row.Cr2O3_liq = repmat(liqOxides.Cr2O3, nP, 1);
row.V2O3_liq = repmat(liqOxides.V2O3, nP, 1);
row.NiO_liq = repmat(liqOxides.NiO, nP, 1);
row.P2O5_liq = repmat(liqOxides.P2O5, nP, 1);
row.SO3_liq = repmat(liqOxides.SO3, nP, 1);
row.Fe2O3_liq = repmat(liqOxides.Fe2O3, nP, 1);
row.F_liq = repmat(liqOxides.F, nP, 1);
row.Cl_liq = repmat(liqOxides.Cl, nP, 1);
row.H2O_liq = repmat(liqOxides.H2O, nP, 1);

% Published liquid cation total and fractions.
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

% Derived terms and diagnostics.
row.XNa_liq = repmat(XNa_liq, nP, 1);
row.XAl_liq = repmat(XAl_liq, nP, 1);
row.XCa_liq = repmat(XCa_liq, nP, 1);
row.XSi_liq = repmat(XSi_liq, nP, 1);
row.XFm_liq = repmat(XFm_liq, nP, 1);
row.MgPrime_liq = repmat(MgPrime_liq, nP, 1);
row.exchangeTerm_DiHd_Jd = repmat(exchangeTerm, nP, 1);
row.KD_FeMg_cpx_liq = repmat(KD_FeMg, nP, 1);
row.invT_10000_over_K = invT;
row.model_domain_valid = modelDomainValid;

% Model B outputs.
row.T2003_K = T2003_K;
row.T2003_C = T2003_C;
row.T_K = T2003_K;
row.T_deg = T2003_C;

end

function cpx = calcCpxComponents(ox, MWinfo)
% calcCpxComponents
% Convert Cpx oxide wt.% to 6-oxygen cations and allocate components using
% Putirka et al. (2003), preserving NaN values.

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

if any(isnan([cpx.XNa, cpx.XAlIV, cpx.XAlVI, cpx.XTi, cpx.XCr]))
    cpx.XFe3_chargeBalance = NaN;
    cpx.XFe3 = NaN;
else
    cpx.XFe3_chargeBalance = cpx.XNa + cpx.XAlIV - cpx.XAlVI ...
        - 2 .* cpx.XTi - cpx.XCr;
    cpx.XFe3 = max(cpx.XFe3_chargeBalance, 0);
end

% Step 1: jadeite.
if isnan(cpx.XAlVI) || isnan(cpx.XNa)
    cpx.XJd = NaN;
else
    cpx.XJd = max(min(cpx.XAlVI, cpx.XNa), 0);
end

% Step 2: Ca-Tschermak.
if isnan(cpx.XAlVI) || isnan(cpx.XJd)
    cpx.XCaTs = NaN;
else
    cpx.XCaTs = max(cpx.XAlVI - cpx.XJd, 0);
end

% Step 3: CaTi component.
if isnan(cpx.XAlIV) || isnan(cpx.XCaTs)
    cpx.XCaTi = NaN;
elseif cpx.XAlIV > cpx.XCaTs
    cpx.XCaTi = max((cpx.XAlIV - cpx.XCaTs) ./ 2, 0);
else
    cpx.XCaTi = 0;
end

% Step 4: CrCaTs component.
if isnan(cpx.XCr)
    cpx.XCrCaTs = NaN;
else
    cpx.XCrCaTs = max(cpx.XCr ./ 2, 0);
end

% Steps 5-6: DiHd and remaining components, including Ca-poor branch.
XFm = cpx.XFe + cpx.XMg;
if any(isnan([cpx.XCa, cpx.XCaTs, cpx.XCaTi, cpx.XCrCaTs, XFm]))
    cpx.XDiHd = NaN;
    cpx.XEnFs = NaN;
    cpx.XFmCaTs = NaN;
    cpx.XFmTi = NaN;
elseif cpx.XCa >= (cpx.XCaTs + cpx.XCaTi)
    cpx.XDiHd = max(cpx.XCa - cpx.XCaTi - cpx.XCaTs ...
        - cpx.XCrCaTs, 0);
    cpx.XFmCaTs = 0;
    cpx.XFmTi = 0;
    cpx.XEnFs = max((XFm - cpx.XDiHd) ./ 2, 0);
else
    cpx.XDiHd = 0;
    cpx.XCaTs = cpx.XCa;
    VIAlex = max(cpx.XAlVI - cpx.XCaTs, 0);
    cpx.XFmCaTs = max(VIAlex - cpx.XCaTs, 0);
    cpx.XFmTi = max((cpx.XAlIV - cpx.XCaTs ...
        - cpx.XFmCaTs) ./ 2, 0);
    cpx.XEnFs = max((XFm - cpx.XFmCaTs - cpx.XFmTi) ./ 2, 0);
end

end

function liq = calcLiquidFractions(ox, MWinfo)
% calcLiquidFractions
% Convert the published Model B liquid oxide list to cation fractions.
% H2O, V, Ni, P, S, Fe2O3, F, and Cl are deliberately omitted.

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
n.Cr2O3 = ox.Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;

cationTotal = n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + ...
    n.MgO + n.CaO + n.Na2O + n.K2O + n.Cr2O3;

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
liq.XCrO1_5 = n.Cr2O3 ./ divisor;

end

function [T_K, T_C] = invTtoT(invT, domainValid)
% invTtoT
% Convert 10000/T values to Kelvin and degreeC only where Model B is valid.

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
% when the FeO column is absent.

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

function displayTemperatureResult(values)
% displayTemperatureResult
% Print one temperature or the first-to-last range for immediate feedback.

if isscalar(values)
    disp(['T2003: ' num2str(values) ' degreeC']);
else
    disp(['T2003: ' num2str(values(1)) ' to ' ...
        num2str(values(end)) ' degreeC']);
end

end

function printTemperatureRangeWarning(values, codeCpx, ...
        minimumTemperature, maximumTemperature)
% printTemperatureRangeWarning
% Warn for finite temperatures outside the direct new-experiment range.

finiteMask = isfinite(values);
outsideMask = finiteMask & ...
    (values < minimumTemperature | values > maximumTemperature);

if any(outsideMask)
    finiteValues = values(finiteMask);
    fprintf(2, ...
        ['WARNING: Putirka et al. (2003) Model B temperature is outside ' ...
         'the direct new-experiment range: 850-1300 degreeC. %d of %d ' ...
         'finite point(s) are outside; calculated finite range = ' ...
         '%.6g-%.6g degreeC for Cpx %s. These bounds are not strict ' ...
         'rectangular limits for the full regression.\n'], ...
        sum(outsideMask), sum(finiteMask), min(finiteValues), ...
        max(finiteValues), char(string(codeCpx)));
end

end

function printNonFiniteWarning(values, codeCpx, idxLiq)
% printNonFiniteWarning
% Retain and report NaN or Inf temperature outputs.

invalidMask = ~isfinite(values);
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Putirka et al. (2003) Model B temperature ' ...
         'values were calculated for Cpx %s and liquid row %d (%d of %d ' ...
         'points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table and the ' ...
         'calculation has not been stopped.\n'], ...
        char(string(codeCpx)), idxLiq, sum(invalidMask), numel(values), ...
        sum(isnan(values)), sum(isinf(values)));
end

end

function printDomainWarning(domainValid, codeCpx, idxLiq)
% printDomainWarning
% Report pressure rows for which Model B is mathematically undefined.

invalidMask = ~domainValid;
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Putirka et al. (2003) Model B is outside its ' ...
         'mathematical domain for Cpx %s and liquid row %d at %d of %d ' ...
         'pressure point(s). Positive finite Jd, DiHd, liquid Na, Al, ' ...
         'Ca, Si, Fe+Mg, MgPrime, exchange term, and 10000/T are ' ...
         'required. Corresponding temperatures remain NaN.\n'], ...
        char(string(codeCpx)), idxLiq, sum(invalidMask), ...
        numel(invalidMask));
end

end

function printKdWarning(values, codeCpx, idxLiq)
% printKdWarning
% Warn when finite Fe-Mg KD falls outside the approximately 3-sigma-filtered
% experimental range reported by Putirka et al. (2003).

finiteMask = isfinite(values);
outsideMask = finiteMask & (values < 0.105 | values > 0.488);
if any(outsideMask)
    finiteValues = values(finiteMask);
    fprintf(2, ...
        ['WARNING: Cpx-liquid Fe-Mg KD is outside the approximately ' ...
         '3-sigma-filtered experimental range of Putirka et al. (2003): ' ...
         '0.105-0.488. Finite KD range = %.6g-%.6g for Cpx %s and ' ...
         'liquid row %d. The selected pair may not represent equilibrium.\n'], ...
        min(finiteValues), max(finiteValues), char(string(codeCpx)), idxLiq);
end

end
