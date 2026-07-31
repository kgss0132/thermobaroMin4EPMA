function results = Putirka2008CpxLiqbaro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Pyroxene/Putirka2008CpxLiqbaro.m
% Compatibility target: MATLAB R2024b
%
% Clinopyroxene-Liquid Al-partitioning barometer, Equation (32c)
% Putirka, K.D. (2008)
% Reviews in Mineralogy and Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Clinopyroxene analysis selected
% from rawdata_struct.Cpx with one Liquid analysis read by
% liquid.readLiquidExcel and calculates pressure using Equation (32c) of
% Putirka (2008; pp. 91-94).
%
% Equation (32c) is based on the partitioning of Al between Clinopyroxene
% and Liquid. It uses temperature, Cpx CaTs and total Al, anhydrous Liquid
% cation fractions, and the H2O content of the equilibrium Liquid in wt%.
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Clinopyroxene-Liquid pair, one
% output row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% CALIBRATION, PERFORMANCE, AND APPLICATION NOTES
%
% Putirka (2008; p. 94) calibrated Equation (32c) using experimental data
% from Putirka et al. (1996), Kinzler and Grove (1992), Sisson and Grove
% (1993a,b), and Walter and Presnall (1994). Reported statistics are:
%
%   Calibration data : R2 = 0.95, SEE = +/-1.5 kbar, n = 99
%   All available data: R2 = 0.82, SEE = +/-5.0 kbar, n = 1303
%
% Important application cautions:
%
%   1) The selected Clinopyroxene and Liquid must represent an equilibrium
%      pair. Incorrect pairing, crystal accumulation, magma mixing,
%      antecrysts, xenocrysts, or altered glass may produce meaningless P.
%
%   2) Equation (32c) is sensitive to Al partitioning. Putirka (2008;
%      pp. 108-110, Fig. 14) showed that rapid cooling can increase the
%      apparent Cpx-Liquid Al partition coefficient and cause pressure
%      overestimation. Cooling-rate disequilibrium must therefore be
%      considered independently.
%
%   3) Slowly diffusing Al may retain an earlier, deeper P-T record. The
%      input temperature should correspond to the temperature represented
%      by the Al partitioning rather than an unrelated late-stage or
%      eruption temperature.
%
%   4) Later testing of H2O-poor tholeiitic experiments found increasingly
%      positive bias at low pressure, although the calibration SEE remains
%      +/-1.5 kbar. Near-surface estimates require particular caution.
%
%   5) Liquid components are anhydrous cation fractions. H2O is excluded
%      from cationTotal_liq and enters Equation (32c) separately in wt%.
%      F and Cl are also excluded from cationTotal_liq.
%
%   6) Clinopyroxene cations and components are calculated on a six-oxygen
%      basis. A good analysis should have a cation sum close to 4.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside the broad approximately
%      679-1650 degreeC Cpx-saturated experimental envelope,
%   2) finite calculated pressure is outside the broad 0-80 kbar global
%      experimental envelope,
%   3) a calculation input contains NaN,
%   4) the Al partitioning ratio is outside its mathematical domain, or
%   5) a calculated pressure is NaN, Inf, or negative.
%
% These broad P-T envelopes are contextual warnings and are not strict
% Equation (32c)-specific rectangular calibration limits.
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained, propagated, and reported. Inf and finite negative
% values are rejected. An absent optional oxide column is treated as zero;
% a present NaN is never silently replaced by zero.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST Cpx-table column is treated as the identifier displayed in the
% selection dialog. Common oxide names and names ending in "_value" are
% recognized. Cpx FeO is preferred; FeOt is used only when FeO is absent.
%
% Required Cpx oxides:
%   SiO2, Al2O3, MgO, CaO, and FeO or FeOt
%
% Optional Cpx oxides, treated as zero only when the column is absent:
%   TiO2, MnO, Na2O, K2O, Cr2O3
%
% Liquid compositions are read with liquid.readLiquidExcel. The following
% Liquid values are required by Equation (32c) or by liquid normalization:
%   SiO2, Al2O3, CaO, H2O, and FeO or FeOt
%
% Other Liquid oxide columns are optional and are included in the anhydrous
% cation total when present. The H2O column is required because Equation
% (32c) contains an explicit H2O term; an unknown H2O value must not be
% silently interpreted as 0 wt%.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION
%
% Putirka (2008), Equation (32c), pp. 91-94:
%
%   P(kbar) =
%       -57.9
%       + 0.0475*T_K
%       - 40.6*XFeO_Liq
%       - 47.7*XCaTs_Cpx
%       + 0.676*H2O_Liq(wt%)
%       - 153*XCaO_Liq*XSiO2_Liq
%       + 6.89*(XAl_Cpx/XAlO1.5_Liq)
%
% where:
%   XAl_Cpx       = total Al atoms in Cpx on a six-oxygen basis
%                 = XAl(IV)_Cpx + XAl(VI)_Cpx
%   XAlO1.5_Liq   = Al cation fraction in the anhydrous Liquid
%   XFeO_Liq      = total FeO-equivalent cation fraction in the Liquid
%   XCaO_Liq      = Ca cation fraction in the anhydrous Liquid
%   XSiO2_Liq     = Si cation fraction in the anhydrous Liquid
%   H2O_Liq       = Liquid H2O in wt%
%
% The Cpx CaTs component is calculated using the Putirka (2008) six-oxygen
% normative scheme:
%
%   XAlIV  = max(2 - XSi, 0)
%   XAlVI  = max(XAl_total - XAlIV, 0)
%   XJd    = min(XNa, XAlVI)
%   XCaTs  = max(XAlVI - XJd, 0)
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008CpxLiqbaro(rawdata_struct, T_degreeC)
%   results = Putirka2008CpxLiqbaro(..., 'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing a Cpx table
%   T_degreeC      : non-negative numeric scalar or vector in degreeC.
%                    NaN is allowed; Inf and finite negative values are not.
%
% Name-value option:
%   'LiquidRow'    : positive integer scalar selecting a Liquid row.
%                    Default [] uses row 1.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Clinopyroxene-Liquid pair.
%

%% Input validation
if nargin < 2
    error('Putirka2008CpxLiqbaro requires (rawdata_struct, T_degreeC).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(T_degreeC) || isempty(T_degreeC) || ~isvector(T_degreeC) || ...
        any(isinf(T_degreeC(:))) || ...
        any(isfinite(T_degreeC(:)) & T_degreeC(:) < 0)
    error(['T_degreeC must be a non-negative numeric scalar or vector. ' ...
           'NaN is allowed and retained, but Inf and finite negative values are prohibited.']);
end
T_degreeC = T_degreeC(:);

%% Options
ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve Clinopyroxene and Liquid datasets
disp('=== Step 1: Preparing Clinopyroxene and Liquid datasets ===');

if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

dataset_cpx = rawdata_struct.Cpx;
if height(dataset_cpx) < 1
    error('rawdata_struct.Cpx is empty.');
end

MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();
if isempty(liqAll) || ~istable(liqAll)
    error('Selected Liquid dataset is empty or is not a table.');
end

if isempty(liquidRowOpt)
    selectedIdx_liq = 1;
    if height(liqAll) > 1
        fprintf(2, ...
            ['WARNING: Liquid dataset has multiple rows (%d). Row 1 will be used. ' ...
             'Specify ''LiquidRow'' to select a different row.\n'], ...
            height(liqAll));
    end
else
    selectedIdx_liq = liquidRowOpt;
    if selectedIdx_liq > height(liqAll)
        error(['Requested LiquidRow (%d) exceeds the number of rows in the ' ...
               'selected Liquid dataset (%d).'], ...
            selectedIdx_liq, height(liqAll));
    end
end

selectedData_liq = liqAll(selectedIdx_liq, :);
checkRequiredLiquidColumns(selectedData_liq);

disp(['Liquid selected: Row ' num2str(selectedIdx_liq)]);
disp('=== Preparing Clinopyroxene and Liquid datasets has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a cell buffer and concatenate once
% after the interactive loop. This avoids repeated growth of the results
% table during normal operation.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max(16, height(dataset_cpx));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

broadT_min_degreeC = 679;
broadT_max_degreeC = 1650;
broadP_min_kbar = 0;
broadP_max_kbar = 80;

temperatureOutsideBroadRange = isfinite(T_degreeC) & ...
    (T_degreeC < broadT_min_degreeC | T_degreeC > broadT_max_degreeC);
temperatureWarningIssued = false;
modelCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Clinopyroxene) ===');

while true
    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Clinopyroxene data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIdx_cpx);
    disp(['Clinopyroxene selected: ' char(string(selectedCode_cpx))]);

    disp('=== Step 4: Checking calculation inputs ===');
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    nanInputNames = findNaNInputs( ...
        selectedData_cpx, selectedData_liq, T_degreeC);
    validateNonNegativeInputs(selectedData_cpx, selectedData_liq);

    disp('=== Step 5: Calculating the pressure ===');
    row = calcPressure( ...
        selectedData_cpx, selectedData_liq, T_degreeC, MWinfo);

    nRows = height(row);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_cpx', 'dataRow_liq'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        expandedBuffer = cell(2 * numel(resultBlocks), 1);
        expandedBuffer(1:numel(resultBlocks)) = resultBlocks;
        resultBlocks = expandedBuffer;
    end
    resultBlocks{nResultBlocks} = row;

    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');
    displayPressureResult(row.P_kbar, selectedCode_cpx, selectedIdx_liq);

    if ~modelCautionIssued
        fprintf(2, ...
            ['CAUTION: Putirka (2008) Equation (32c) requires an equilibrium ' ...
             'Clinopyroxene-Liquid pair and is sensitive to Al disequilibrium. ' ...
             'Rapid cooling can raise DAl(Cpx-Liq) and overestimate pressure ' ...
             '(Putirka 2008, Fig. 14). The H2O term uses Liquid H2O in wt%%.\n']);
        modelCautionIssued = true;
    end

    if any(temperatureOutsideBroadRange) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the broad approximately ' ...
             '679-1650 degreeC Cpx-saturated experimental envelope summarized ' ...
             'by Putirka (2008). %d of %d finite point(s) are outside; finite ' ...
             'input range = %.4g-%.4g degreeC. This is not a strict Equation ' ...
             '(32c)-specific rectangular limit.\n'], ...
            sum(temperatureOutsideBroadRange), numel(finiteTemperature), ...
            min(finiteTemperature), max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    printPressureRangeWarning(row.P_kbar, selectedCode_cpx, ...
        selectedIdx_liq, broadP_min_kbar, broadP_max_kbar);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in Putirka (2008) Equation (32c) input(s) ' ...
             'for %s & Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated. H2O, F, and Cl ' ...
             'are excluded from cationTotal_liq; only H2O enters Equation ' ...
             '(32c), separately in wt%%.\n'], ...
            char(string(selectedCode_cpx)), selectedIdx_liq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    printNonFinitePressureWarning( ...
        row.P_kbar, selectedCode_cpx, selectedIdx_liq);
    printNegativePressureWarning( ...
        row.P_kbar, selectedCode_cpx, selectedIdx_liq);

    invalidDomain = ~row.Eq32c_domain_valid;
    if any(invalidDomain)
        fprintf(2, ...
            ['WARNING: Equation (32c) is outside its mathematical domain for ' ...
             '%s & Liquid row %d at %d of %d point(s). The calculation ' ...
             'requires finite temperature, Cpx and Liquid terms, a positive ' ...
             'Liquid cation total, and XAlO1.5_Liq > 0. Corresponding pressure ' ...
             'values remain NaN or Inf.\n'], ...
            char(string(selectedCode_cpx)), selectedIdx_liq, ...
            sum(invalidDomain), numel(invalidDomain));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Clinopyroxene selection (same Liquid row)?', ...
        'Putirka2008CpxLiqbaro', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct('liquid', metaLiq);
disp('=== Putirka2008CpxLiqbaro finished ===');

end

%% ---- local functions ----
function checkRequiredLiquidColumns(data_liq)
% checkRequiredLiquidColumns
% Stop before the Cpx selection loop when an Equation (32c) Liquid input is
% absent. A present NaN is allowed and is reported after calculation.

requiredOxides = {'SiO2', 'Al2O3', 'CaO', 'H2O'};
missingBuffer = strings(numel(requiredOxides) + 1, 1);
nMissing = 0;

for i = 1:numel(requiredOxides)
    if isempty(findOxideColumn( ...
            data_liq.Properties.VariableNames, requiredOxides{i}))
        nMissing = nMissing + 1;
        missingBuffer(nMissing) = string(requiredOxides{i});
    end
end

feOColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeO');
feOtColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeOt');
if isempty(feOColumn) && isempty(feOtColumn)
    nMissing = nMissing + 1;
    missingBuffer(nMissing) = "FeO or FeOt";
end

if nMissing > 0
    missingNames = missingBuffer(1:nMissing);
    error(['Liquid table must contain Equation (32c) input column(s): ' ...
           char(strjoin(missingNames, ', ')) '.']);
end

end

function nanInputNames = findNaNInputs(data_cpx, data_liq, T_degreeC)
% findNaNInputs
% Return names of Equation (32c) calculation inputs containing NaN. H2O is
% checked because it enters as wt%. F and Cl are intentionally excluded.

maxNames = 48;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

cpxOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO','Na2O','K2O','Cr2O3'};
for i = 1:numel(cpxOxides)
    columnName = findOxideColumn(data_cpx.Properties.VariableNames, cpxOxides{i});
    if ~isempty(columnName) && isnan(toScalarDouble(data_cpx.(columnName)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Cpx." + string(columnName);
    end
end

feColumnName = findOxideColumn(data_cpx.Properties.VariableNames, 'FeO');
if isempty(feColumnName)
    feColumnName = findOxideColumn(data_cpx.Properties.VariableNames, 'FeOt');
end
if ~isempty(feColumnName) && isnan(toScalarDouble(data_cpx.(feColumnName)))
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Cpx." + string(feColumnName);
end

liqOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO','Na2O','K2O', ...
    'V2O3','Cr2O3','NiO','P2O5','SO3','Fe2O3','H2O'};
for i = 1:numel(liqOxides)
    columnName = findOxideColumn(data_liq.Properties.VariableNames, liqOxides{i});
    if ~isempty(columnName) && isnan(toScalarDouble(data_liq.(columnName)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Liquid." + string(columnName);
    end
end

feLiqColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeO');
if isempty(feLiqColumn)
    feLiqColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeOt');
end
if ~isempty(feLiqColumn) && isnan(toScalarDouble(data_liq.(feLiqColumn)))
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Liquid." + string(feLiqColumn);
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_cpx, data_liq)
% validateNonNegativeInputs
% Reject Inf and finite negative values. Zero and NaN are retained. F and
% Cl are excluded because they do not enter the calculation.

maxNames = 48;
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

cpxOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO','Na2O','K2O','Cr2O3'};
for i = 1:numel(cpxOxides)
    columnName = findOxideColumn(data_cpx.Properties.VariableNames, cpxOxides{i});
    if ~isempty(columnName)
        value = toScalarDouble(data_cpx.(columnName));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = "Cpx." + string(columnName);
        end
    end
end

feColumnName = findOxideColumn(data_cpx.Properties.VariableNames, 'FeO');
if isempty(feColumnName)
    feColumnName = findOxideColumn(data_cpx.Properties.VariableNames, 'FeOt');
end
if ~isempty(feColumnName)
    value = toScalarDouble(data_cpx.(feColumnName));
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = "Cpx." + string(feColumnName);
    end
end

liqOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO','Na2O','K2O', ...
    'V2O3','Cr2O3','NiO','P2O5','SO3','Fe2O3','H2O'};
for i = 1:numel(liqOxides)
    columnName = findOxideColumn(data_liq.Properties.VariableNames, liqOxides{i});
    if ~isempty(columnName)
        value = toScalarDouble(data_liq.(columnName));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = "Liquid." + string(columnName);
        end
    end
end

feLiqColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeO');
if isempty(feLiqColumn)
    feLiqColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeOt');
end
if ~isempty(feLiqColumn)
    value = toScalarDouble(data_liq.(feLiqColumn));
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = "Liquid." + string(feLiqColumn);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['Putirka2008CpxLiqbaro: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative values are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_cpx, data_liq, T_degreeC, MWinfo)
% calcPressure
% Calculate Putirka (2008) Equation (32c) for one selected
% Clinopyroxene-Liquid pair over a scalar or vector of temperatures.

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

cpx = prepareCpxRow(data_cpx, MWinfo);
liq = extractLiquidRow(data_liq, MWinfo);

% Expand composition-dependent scalars to the temperature-vector length.
XSi_cpx = repmat(cpx.XSi, nT, 1);
XTi_cpx = repmat(cpx.XTi, nT, 1);
XAl_cpx = repmat(cpx.XAl, nT, 1);
XFe_cpx = repmat(cpx.XFe, nT, 1);
XMn_cpx = repmat(cpx.XMn, nT, 1);
XMg_cpx = repmat(cpx.XMg, nT, 1);
XCa_cpx = repmat(cpx.XCa, nT, 1);
XNa_cpx = repmat(cpx.XNa, nT, 1);
XK_cpx = repmat(cpx.XK, nT, 1);
XCr_cpx = repmat(cpx.XCr, nT, 1);
cationSum_cpx = repmat(cpx.cationSum, nT, 1);
XAlIV_cpx = repmat(cpx.XAlIV, nT, 1);
XAlVI_cpx = repmat(cpx.XAlVI, nT, 1);
XFe3_cpx = repmat(cpx.XFe3, nT, 1);
XJd_cpx = repmat(cpx.XJd, nT, 1);
XCaTs_cpx = repmat(cpx.XCaTs, nT, 1);
XCaTi_cpx = repmat(cpx.XCaTi, nT, 1);
XCrCaTs_cpx = repmat(cpx.XCrCaTs, nT, 1);
XDiHd_cpx = repmat(cpx.XDiHd, nT, 1);
XEnFs_cpx = repmat(cpx.XEnFs, nT, 1);

SiO2_liq = repmat(liq.SiO2, nT, 1);
TiO2_liq = repmat(liq.TiO2, nT, 1);
Al2O3_liq = repmat(liq.Al2O3, nT, 1);
FeO_liq = repmat(liq.FeO, nT, 1);
MnO_liq = repmat(liq.MnO, nT, 1);
MgO_liq = repmat(liq.MgO, nT, 1);
CaO_liq = repmat(liq.CaO, nT, 1);
Na2O_liq = repmat(liq.Na2O, nT, 1);
K2O_liq = repmat(liq.K2O, nT, 1);
V2O3_liq = repmat(liq.V2O3, nT, 1);
Cr2O3_liq = repmat(liq.Cr2O3, nT, 1);
NiO_liq = repmat(liq.NiO, nT, 1);
P2O5_liq = repmat(liq.P2O5, nT, 1);
SO3_liq = repmat(liq.SO3, nT, 1);
Fe2O3_liq = repmat(liq.Fe2O3, nT, 1);
H2O_liq = repmat(liq.H2O, nT, 1);
F_liq = repmat(liq.F, nT, 1);
Cl_liq = repmat(liq.Cl, nT, 1);

cationTotal_liq = repmat(liq.cationTotal, nT, 1);
XSiO2_liq = repmat(liq.XSiO2, nT, 1);
XTiO2_liq = repmat(liq.XTiO2, nT, 1);
XAlO1_5_liq = repmat(liq.XAlO1_5, nT, 1);
XFeO_liq = repmat(liq.XFeO, nT, 1);
XMnO_liq = repmat(liq.XMnO, nT, 1);
XMgO_liq = repmat(liq.XMgO, nT, 1);
XCaO_liq = repmat(liq.XCaO, nT, 1);
XNaO0_5_liq = repmat(liq.XNaO0_5, nT, 1);
XKO0_5_liq = repmat(liq.XKO0_5, nT, 1);

Al_partition_ratio = XAl_cpx ./ XAlO1_5_liq;
CaSi_product_liq = XCaO_liq .* XSiO2_liq;

% Putirka (2008), Equation (32c), pp. 91-94. The printed coefficient
% 0.676 for the H2O term is retained.
P_kbar = ...
    -57.9 ...
    + 0.0475 .* T_K ...
    - 40.6 .* XFeO_liq ...
    - 47.7 .* XCaTs_cpx ...
    + 0.676 .* H2O_liq ...
    - 153 .* CaSi_product_liq ...
    + 6.89 .* Al_partition_ratio;

P_GPa = P_kbar ./ 10;

Eq32c_domain_valid = ...
    isfinite(T_K) & ...
    isfinite(cationTotal_liq) & cationTotal_liq > 0 & ...
    isfinite(XFeO_liq) & ...
    isfinite(XCaTs_cpx) & ...
    isfinite(H2O_liq) & ...
    isfinite(XCaO_liq) & ...
    isfinite(XSiO2_liq) & ...
    isfinite(XAl_cpx) & ...
    isfinite(XAlO1_5_liq) & XAlO1_5_liq > 0 & ...
    isfinite(Al_partition_ratio);

% Pack outputs using equal-length vectors.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.XSi_cpx = XSi_cpx;
row.XTi_cpx = XTi_cpx;
row.XAl_cpx = XAl_cpx;
row.XFe_cpx = XFe_cpx;
row.XMn_cpx = XMn_cpx;
row.XMg_cpx = XMg_cpx;
row.XCa_cpx = XCa_cpx;
row.XNa_cpx = XNa_cpx;
row.XK_cpx = XK_cpx;
row.XCr_cpx = XCr_cpx;
row.cationSum_cpx = cationSum_cpx;
row.XAlIV_cpx = XAlIV_cpx;
row.XAlVI_cpx = XAlVI_cpx;
row.XFe3_cpx = XFe3_cpx;
row.XJd_cpx = XJd_cpx;
row.XCaTs_cpx = XCaTs_cpx;
row.XCaTi_cpx = XCaTi_cpx;
row.XCrCaTs_cpx = XCrCaTs_cpx;
row.XDiHd_cpx = XDiHd_cpx;
row.XEnFs_cpx = XEnFs_cpx;

row.SiO2_liq = SiO2_liq;
row.TiO2_liq = TiO2_liq;
row.Al2O3_liq = Al2O3_liq;
row.FeO_liq = FeO_liq;
row.MnO_liq = MnO_liq;
row.MgO_liq = MgO_liq;
row.CaO_liq = CaO_liq;
row.Na2O_liq = Na2O_liq;
row.K2O_liq = K2O_liq;
row.V2O3_liq = V2O3_liq;
row.Cr2O3_liq = Cr2O3_liq;
row.NiO_liq = NiO_liq;
row.P2O5_liq = P2O5_liq;
row.SO3_liq = SO3_liq;
row.Fe2O3_liq = Fe2O3_liq;
row.H2O_liq = H2O_liq;
row.F_liq = F_liq;
row.Cl_liq = Cl_liq;

row.cationTotal_liq = cationTotal_liq;
row.XSiO2_liq = XSiO2_liq;
row.XTiO2_liq = XTiO2_liq;
row.XAlO1_5_liq = XAlO1_5_liq;
row.XFeO_liq = XFeO_liq;
row.XMnO_liq = XMnO_liq;
row.XMgO_liq = XMgO_liq;
row.XCaO_liq = XCaO_liq;
row.XNaO0_5_liq = XNaO0_5_liq;
row.XKO0_5_liq = XKO0_5_liq;

row.Al_partition_ratio_CpxLiq = Al_partition_ratio;
row.CaSi_product_liq = CaSi_product_liq;
row.Eq32c_domain_valid = Eq32c_domain_valid;

% Primary pressure variables used by the common barometer launchers.
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% Explicit Equation (32c) aliases.
row.PEq32c_kbar = P_kbar;
row.PEq32c_GPa = P_GPa;

row.P_SEE_calibration_kbar = repmat(1.5, nT, 1);
row.P_SEE_allData_kbar = repmat(5.0, nT, 1);
row.P_R2_calibration = repmat(0.95, nT, 1);
row.P_R2_allData = repmat(0.82, nT, 1);
row.isWithinBroadTRange = isfinite(T_degreeC) & ...
    T_degreeC >= 679 & T_degreeC <= 1650;
row.isWithinBroadPRange = isfinite(P_kbar) & ...
    P_kbar >= 0 & P_kbar <= 80;

end

function cpx = prepareCpxRow(data_cpx, MWinfo)
% prepareCpxRow
% Extract one-row Clinopyroxene oxide data and calculate six-oxygen cations
% and normative components. Existing NaN values are retained.

if height(data_cpx) ~= 1
    error('Clinopyroxene input must be a 1-row table.');
end

SiO2 = getMineralOxRequired(data_cpx, 'SiO2', 'Clinopyroxene');
TiO2 = getMineralOxOptional(data_cpx, 'TiO2', 0);
Al2O3 = getMineralOxRequired(data_cpx, 'Al2O3', 'Clinopyroxene');
MnO = getMineralOxOptional(data_cpx, 'MnO', 0);
MgO = getMineralOxRequired(data_cpx, 'MgO', 'Clinopyroxene');
CaO = getMineralOxRequired(data_cpx, 'CaO', 'Clinopyroxene');
Na2O = getMineralOxOptional(data_cpx, 'Na2O', 0);
K2O = getMineralOxOptional(data_cpx, 'K2O', 0);
Cr2O3 = getMineralOxOptional(data_cpx, 'Cr2O3', 0);

% Prefer FeO when the column exists. If FeO exists but is NaN, retain that
% NaN rather than silently replacing it with FeOt. Use FeOt only when the
% FeO column is absent.
feOColumn = findOxideColumn(data_cpx.Properties.VariableNames, 'FeO');
if ~isempty(feOColumn)
    FeO = toScalarDouble(data_cpx.(feOColumn));
else
    feOtColumn = findOxideColumn(data_cpx.Properties.VariableNames, 'FeOt');
    if isempty(feOtColumn)
        error('Clinopyroxene table must contain FeO or FeOt.');
    end
    FeO = toScalarDouble(data_cpx.(feOtColumn));
end

molProp = struct();
molProp.SiO2 = SiO2 ./ MWinfo.MW.SiO2;
molProp.TiO2 = TiO2 ./ MWinfo.MW.TiO2;
molProp.Al2O3 = Al2O3 ./ MWinfo.MW.Al2O3;
molProp.FeO = FeO ./ MWinfo.MW.FeO;
molProp.MnO = MnO ./ MWinfo.MW.MnO;
molProp.MgO = MgO ./ MWinfo.MW.MgO;
molProp.CaO = CaO ./ MWinfo.MW.CaO;
molProp.Na2O = Na2O ./ MWinfo.MW.Na2O;
molProp.K2O = K2O ./ MWinfo.MW.K2O;
molProp.Cr2O3 = Cr2O3 ./ MWinfo.MW.Cr2O3;

oxygenSum = ...
    2 .* molProp.SiO2 + ...
    2 .* molProp.TiO2 + ...
    3 .* molProp.Al2O3 + ...
    molProp.FeO + ...
    molProp.MnO + ...
    molProp.MgO + ...
    molProp.CaO + ...
    molProp.Na2O + ...
    molProp.K2O + ...
    3 .* molProp.Cr2O3;

oxygenRenormalizationFactor = 6 ./ oxygenSum;

XSi = molProp.SiO2 .* oxygenRenormalizationFactor;
XTi = molProp.TiO2 .* oxygenRenormalizationFactor;
XAl = 2 .* molProp.Al2O3 .* oxygenRenormalizationFactor;
XFe = molProp.FeO .* oxygenRenormalizationFactor;
XMn = molProp.MnO .* oxygenRenormalizationFactor;
XMg = molProp.MgO .* oxygenRenormalizationFactor;
XCa = molProp.CaO .* oxygenRenormalizationFactor;
XNa = 2 .* molProp.Na2O .* oxygenRenormalizationFactor;
XK = 2 .* molProp.K2O .* oxygenRenormalizationFactor;
XCr = 2 .* molProp.Cr2O3 .* oxygenRenormalizationFactor;

cationSum = XSi + XTi + XAl + XFe + XMn + XMg + XCa + XNa + XK + XCr;

XAlIV = maxPreserveNaN(2 - XSi, 0);
XAlVI = maxPreserveNaN(XAl - XAlIV, 0);

XFe3 = XNa + XAlIV - XAlVI - 2 .* XTi - XCr;
XFe3 = maxPreserveNaN(XFe3, 0);

XJd = minPreserveNaN(XAlVI, XNa);
XJd = maxPreserveNaN(XJd, 0);

XCaTs = maxPreserveNaN(XAlVI - XJd, 0);
XCaTi = maxPreserveNaN((XAlIV - XCaTs) ./ 2, 0);
XCrCaTs = maxPreserveNaN(XCr ./ 2, 0);
XDiHd = maxPreserveNaN(XCa - XCaTi - XCaTs - XCrCaTs, 0);
XEnFs = maxPreserveNaN((XFe + XMg - XDiHd) ./ 2, 0);

cpx = struct();
cpx.XSi = XSi;
cpx.XTi = XTi;
cpx.XAl = XAl;
cpx.XFe = XFe;
cpx.XMn = XMn;
cpx.XMg = XMg;
cpx.XCa = XCa;
cpx.XNa = XNa;
cpx.XK = XK;
cpx.XCr = XCr;
cpx.cationSum = cationSum;
cpx.XAlIV = XAlIV;
cpx.XAlVI = XAlVI;
cpx.XFe3 = XFe3;
cpx.XJd = XJd;
cpx.XCaTs = XCaTs;
cpx.XCaTi = XCaTi;
cpx.XCrCaTs = XCrCaTs;
cpx.XDiHd = XDiHd;
cpx.XEnFs = XEnFs;

end

function liq = extractLiquidRow(data_liq, MWinfo)
% extractLiquidRow
% Extract one-row Liquid oxide data and calculate anhydrous cation
% fractions. H2O enters Equation (32c) as wt% but is excluded from
% cationTotal_liq. F and Cl are retained only as diagnostics.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liq = struct();
liq.SiO2 = getLiquidOxRequired(data_liq, 'SiO2', 'Liquid');
liq.TiO2 = getLiquidOxOptional(data_liq, 'TiO2', 0);
liq.Al2O3 = getLiquidOxRequired(data_liq, 'Al2O3', 'Liquid');

feOColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeO');
if ~isempty(feOColumn)
    liq.FeO = toScalarDouble(data_liq.(feOColumn));
else
    feOtColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeOt');
    if isempty(feOtColumn)
        error('Liquid table must contain FeO or FeOt.');
    end
    liq.FeO = toScalarDouble(data_liq.(feOtColumn));
end

liq.MnO = getLiquidOxOptional(data_liq, 'MnO', 0);
liq.MgO = getLiquidOxOptional(data_liq, 'MgO', 0);
liq.CaO = getLiquidOxRequired(data_liq, 'CaO', 'Liquid');
liq.Na2O = getLiquidOxOptional(data_liq, 'Na2O', 0);
liq.K2O = getLiquidOxOptional(data_liq, 'K2O', 0);
liq.V2O3 = getLiquidOxOptional(data_liq, 'V2O3', 0);
liq.Cr2O3 = getLiquidOxOptional(data_liq, 'Cr2O3', 0);
liq.NiO = getLiquidOxOptional(data_liq, 'NiO', 0);
liq.P2O5 = getLiquidOxOptional(data_liq, 'P2O5', 0);
liq.SO3 = getLiquidOxOptional(data_liq, 'SO3', 0);
liq.Fe2O3 = getLiquidOxOptional(data_liq, 'Fe2O3', 0);
liq.H2O = getLiquidOxRequired(data_liq, 'H2O', 'Liquid');
liq.F = getLiquidOxOptional(data_liq, 'F', NaN);
liq.Cl = getLiquidOxOptional(data_liq, 'Cl', NaN);

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

liq.cationTotal = ...
    n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + n.MgO + ...
    n.CaO + n.Na2O + n.K2O + n.V2O3 + n.Cr2O3 + n.NiO + ...
    n.P2O5 + n.SO3 + n.Fe2O3;

liq.XSiO2 = n.SiO2 ./ liq.cationTotal;
liq.XTiO2 = n.TiO2 ./ liq.cationTotal;
liq.XAlO1_5 = n.Al2O3 ./ liq.cationTotal;
liq.XFeO = n.FeO ./ liq.cationTotal;
liq.XMnO = n.MnO ./ liq.cationTotal;
liq.XMgO = n.MgO ./ liq.cationTotal;
liq.XCaO = n.CaO ./ liq.cationTotal;
liq.XNaO0_5 = n.Na2O ./ liq.cationTotal;
liq.XKO0_5 = n.K2O ./ liq.cationTotal;

end

function displayPressureResult(values, selectedCode_cpx, selectedIdx_liq)
% displayPressureResult
% Print one value or the first-to-last pressure range.

if isscalar(values)
    disp(['Equation 32c | ' char(string(selectedCode_cpx)) ...
        ' & Liquid row ' num2str(selectedIdx_liq) ': P = ' ...
        num2str(values) ' kbar']);
else
    disp(['Equation 32c | ' char(string(selectedCode_cpx)) ...
        ' & Liquid row ' num2str(selectedIdx_liq) ': P = ' ...
        num2str(values(1)) ' to ' num2str(values(end)) ' kbar']);
end

end

function printPressureRangeWarning(values, selectedCode_cpx, ...
        selectedIdx_liq, minimumPressure, maximumPressure)
% printPressureRangeWarning
% Warn for finite pressures outside the broad 0-80 kbar global envelope.

finiteMask = isfinite(values);
outsideMask = finiteMask & ...
    (values < minimumPressure | values > maximumPressure);

if any(outsideMask)
    finiteValues = values(finiteMask);
    fprintf(2, ...
        ['WARNING: Equation (32c) pressure is outside the broad 0-80 kbar ' ...
         'global experimental envelope summarized by Putirka (2008). %d of ' ...
         '%d finite point(s) are outside; finite calculated range = ' ...
         '%.4g-%.4g kbar for %s & Liquid row %d. This is not a strict ' ...
         'Equation (32c)-specific calibration range.\n'], ...
        sum(outsideMask), sum(finiteMask), ...
        min(finiteValues), max(finiteValues), ...
        char(string(selectedCode_cpx)), selectedIdx_liq);
end

end

function printNonFinitePressureWarning(values, selectedCode_cpx, selectedIdx_liq)
% printNonFinitePressureWarning
% Retain and report NaN or Inf pressure outputs.

invalidMask = ~isfinite(values);
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Equation (32c) pressure values were calculated ' ...
         'for %s & Liquid row %d (%d of %d points; NaN: %d, Inf: %d). ' ...
         'Values remain in the output table.\n'], ...
        char(string(selectedCode_cpx)), selectedIdx_liq, ...
        sum(invalidMask), numel(values), ...
        sum(isnan(values)), sum(isinf(values)));
end

end

function printNegativePressureWarning(values, selectedCode_cpx, selectedIdx_liq)
% printNegativePressureWarning
% Retain negative finite values for diagnostics and report them. Putirka
% (2008) notes that small negative values can fall within error near 1 atm.

negativeMask = isfinite(values) & values < 0;
if any(negativeMask)
    fprintf(2, ...
        ['WARNING: Negative finite Equation (32c) pressure was calculated ' ...
         'for %s & Liquid row %d (%d of %d points). Values were retained ' ...
         'for diagnostic purposes; small negative values may lie within ' ...
         'model error near 1 atm, whereas large negative values require ' ...
         'input and equilibrium checks.\n'], ...
        char(string(selectedCode_cpx)), selectedIdx_liq, ...
        sum(negativeMask), numel(values));
end

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Copy available Liquid identifiers to every temperature row without
% changing the output table height.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    rawValue = data_liq.('Index');
    if isnumeric(rawValue) || islogical(rawValue)
        row.liq_Index = repmat(double(rawValue(1)), nRows, 1);
    else
        row.liq_Index = repmat(string(rawValue(1)), nRows, 1);
    end
end
if any(strcmp(variableNames, 'Experiment'))
    rawValue = data_liq.('Experiment');
    row.liq_Experiment = repmat(string(rawValue(1)), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    rawValue = data_liq.('Citation');
    row.liq_Citation = repmat(string(rawValue(1)), nRows, 1);
end

end

function value = getMineralOxRequired(data_tbl, oxide, mineralLabel)
% getMineralOxRequired
% Retrieve a required scalar oxide value. A present NaN is retained.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    error('%s table must contain variable: %s', mineralLabel, oxide);
end
value = toScalarDouble(data_tbl.(columnName));

end

function value = getMineralOxOptional(data_tbl, oxide, missingDefault)
% getMineralOxOptional
% Retrieve an optional scalar oxide value. An absent column uses the stated
% backward-compatible default; a present NaN remains NaN.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    value = missingDefault;
else
    value = toScalarDouble(data_tbl.(columnName));
end

end

function value = getLiquidOxRequired(data_tbl, oxide, phaseLabel)
% getLiquidOxRequired
% Retrieve a required scalar Liquid oxide value. A present NaN is retained.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    error('%s table must contain variable: %s', phaseLabel, oxide);
end
value = toScalarDouble(data_tbl.(columnName));

end

function value = getLiquidOxOptional(data_tbl, oxide, missingDefault)
% getLiquidOxOptional
% Retrieve an optional scalar Liquid oxide value. An absent column uses the
% stated default; a present NaN remains NaN.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    value = missingDefault;
else
    value = toScalarDouble(data_tbl.(columnName));
end

end

function columnName = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match common oxide variable names after removing spaces, underscores, and
% hyphens. Both "Oxide" and "Oxide_value" forms are recognized.

canonicalNames = cell(size(variableNames));
for i = 1:numel(variableNames)
    nameText = lower(variableNames{i});
    nameText = strrep(nameText, ' ', '');
    nameText = strrep(nameText, '_', '');
    nameText = strrep(nameText, '-', '');
    canonicalNames{i} = nameText;
end

canonicalOxide = lower(oxide);
canonicalOxide = strrep(canonicalOxide, ' ', '');
canonicalOxide = strrep(canonicalOxide, '_', '');
canonicalOxide = strrep(canonicalOxide, '-', '');

targetNames = {[canonicalOxide 'value'], canonicalOxide};
columnName = '';

for i = 1:numel(targetNames)
    matchedIndex = find(strcmp(canonicalNames, targetNames{i}), 1, 'first');
    if ~isempty(matchedIndex)
        columnName = variableNames{matchedIndex};
        return
    end
end

end

function value = toScalarDouble(rawValue)
% toScalarDouble
% Convert a one-row table value to a scalar double. Missing, empty, or
% non-numeric text is represented by NaN and is not converted to zero.

value = NaN;

if isempty(rawValue)
    return
end

if isnumeric(rawValue) || islogical(rawValue)
    value = double(rawValue(1));
    return
end

if isstring(rawValue)
    if ismissing(rawValue(1))
        return
    end
    value = str2double(rawValue(1));
    return
end

if ischar(rawValue)
    value = str2double(string(rawValue));
    return
end

if iscategorical(rawValue)
    value = str2double(string(rawValue(1)));
    return
end

if iscell(rawValue)
    if isempty(rawValue{1})
        return
    end
    cellValue = rawValue{1};
    if isnumeric(cellValue) || islogical(cellValue)
        value = double(cellValue(1));
    elseif isstring(cellValue) || ischar(cellValue) || iscategorical(cellValue)
        value = str2double(string(cellValue));
    end
end

end

function value = maxPreserveNaN(a, b)
% maxPreserveNaN
% Scalar maximum that explicitly preserves NaN.

if isnan(a) || isnan(b)
    value = NaN;
else
    value = max(a, b);
end

end

function value = minPreserveNaN(a, b)
% minPreserveNaN
% Scalar minimum that explicitly preserves NaN.

if isnan(a) || isnan(b)
    value = NaN;
else
    value = min(a, b);
end

end
