function results = Putirka2008CpxJdLiq(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Pyroxene/Putirka2008CpxJdLiq.m
% Compatibility target: MATLAB R2024b
%
% Clinopyroxene-Liquid Jd barometers, Equations (30) and (31)
% Putirka, K. D. (2008)
% Reviews in Mineralogy and Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Clinopyroxene analysis selected
% from rawdata_struct.Cpx with one Liquid analysis read by
% liquid.readLiquidExcel. It calculates pressure using Putirka (2008)
% Equations (30) and (31), both based on Jd-Liquid equilibrium.
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Clinopyroxene-Liquid pair, one
% output row is returned for every input temperature value.
%
% Equation (30) is the calibration based on the selected experimental
% studies listed on p. 91. Equation (31) is the global calibration of the
% clinopyroxene-saturated experimental database. The standard launcher
% outputs P_kbar and P_GPa are assigned from Equation (31), because it has
% the smaller global error. Results from both equations are retained as
% P30_kbar/P30_GPa and P31_kbar/P31_GPa.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE, PERFORMANCE, AND APPLICATION NOTES
%
% Putirka (2008; pp. 88-94, Figs. 8-9) evaluated Cpx-Liquid barometers over
% a broad experimental database. Equation (31) was globally calibrated
% after excluding experiments at 1 atm, except Grove and Juster (1989), and
% experiments above 40 kbar (p. 91). The following broad envelopes are used
% here only for non-stopping warnings:
%
%   Temperature : approximately 679-1650 degreeC
%   Pressure    : 0-40 kbar
%
% These are broad database envelopes and are not strict rectangular limits
% for every observation used in either equation.
%
% Reported model performance (Fig. 9):
%
%   Equation (30), calibration data:
%       R2 = 0.97, SEE = +/-1.6 kbar, n = 117
%   Equation (30), all data at P < 40 kbar:
%       R2 = 0.83, SEE = +/-3.6 kbar, n = 1066
%   Equation (31), all data at P < 40 kbar:
%       R2 = 0.89, SEE = +/-2.9 kbar, n = 1176
%
% Important application cautions:
%
%   1) The selected Clinopyroxene and Liquid must represent an equilibrium
%      pair. Putirka (2008, pp. 94 and 107-110) recommends equilibrium tests
%      before interpreting P-T results. The broad experimental mean is
%      KD(Fe-Mg)Cpx-Liq = 0.28 +/- 0.08, but matching Fe-Mg exchange alone
%      does not guarantee equilibrium of Na-Al-Ca components.
%
%   2) Whole-rock compositions should be used as Liquid only when they
%      reasonably represent the melt from which the selected Cpx
%      crystallized. Crystal accumulation, mixing, fractionation, xenocrysts,
%      and antecrysts can invalidate a nominal Cpx-Liquid pair.
%
%   3) Liquid components are anhydrous cation fractions. H2O enters the
%      equations separately in wt%. Clinopyroxene cations and normative
%      components are calculated on a six-oxygen basis following Table 3.
%
%   4) A good Cpx analysis should have a six-oxygen cation sum close to 4.
%      Large deviations may indicate analytical problems, mixed analyses,
%      alteration, or a non-pyroxene composition.
%
%   5) Equation (30) requires positive FeO, MgO, KO0.5, and DiHd terms,
%      because their natural logarithms occur explicitly. Equation (31)
%      requires positive DiHd, EnFs, and total Cpx Al. Both equations require
%      a positive Jd-Liquid equilibrium term.
%
%   6) The models return systematically positive pressures for many 1-atm
%      experiments. Figure 9 reports means of approximately 2.8 +/-3.0 kbar
%      for Equation (30) and 2.7 +/-3.2 kbar for Equation (31). Near-surface
%      estimates must therefore be interpreted with particular caution.
%
%   7) Model error and error caused by incorrect mineral-melt pairing are
%      separate. Averaging many calculations does not remove systematic
%      error from disequilibrium or inappropriate compositions.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 679-1650 degreeC,
%   2) finite calculated pressure is outside 0-40 kbar,
%   3) a calculation input contains NaN,
%   4) an equation logarithm argument is outside its domain, or
%   5) a calculated pressure is NaN, Inf, or negative.
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained, propagated, and reported. Inf and finite negative
% inputs are rejected. An absent optional oxide column is treated as zero;
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
% Liquid compositions are read with liquid.readLiquidExcel. H2O is read as
% wt% and defaults to zero only when the column is absent. H2O, F, and Cl do
% not enter cationTotal_liq. F and Cl are retained only for traceability.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION
%
% Define the Jd-Liquid equilibrium term:
%
%   K_Jd_Liq = XJd_Cpx / ...
%       (XNaO0.5_Liq * XAlO1.5_Liq * XSiO2_Liq^2)
%
% and:
%
%   MgNum_Liq = XMgO_Liq / (XMgO_Liq + XFeO_Liq)
%
% Putirka (2008), Equation (30), p. 91:
%
%   P30(kbar) =
%       -48.7
%       + 271*(T_K/10000)
%       + 32*(T_K/10000)*ln(K_Jd_Liq)
%       - 8.2*ln(XFeO_Liq)
%       + 4.6*ln(XMgO_Liq)
%       - 0.96*ln(XKO0.5_Liq)
%       - 2.2*ln(XDiHd_Cpx)
%       - 31*MgNum_Liq
%       + 56*(XNaO0.5_Liq + XKO0.5_Liq)
%       + 0.76*H2O_Liq(wt%)
%
% Putirka (2008), Equation (31), p. 91:
%
%   P31(kbar) =
%       -40.73
%       + 358*(T_K/10000)
%       + 21.69*(T_K/10000)*ln(K_Jd_Liq)
%       - 105.7*XCaO_Liq
%       - 165.5*(XNaO0.5_Liq + XKO0.5_Liq)^2
%       - 50.15*XSiO2_Liq*(XFeO_Liq + XMgO_Liq)
%       - 3.178*ln(XDiHd_Cpx)
%       - 2.205*ln(XEnFs_Cpx)
%       + 0.864*ln(XAl_Cpx)
%       + 0.3962*H2O_Liq(wt%)
%
% Pressure is returned in kbar and GPa. Temperature is supplied in degreeC
% and converted internally to Kelvin.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008CpxJdLiq(rawdata_struct, T_degreeC)
%   results = Putirka2008CpxJdLiq(..., 'LiquidRow', rowNumber)
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
    error('Putirka2008CpxJdLiq requires (rawdata_struct, T_degreeC).');
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
disp(['Liquid selected: Row ' num2str(selectedIdx_liq)]);
disp('=== Preparing Clinopyroxene and Liquid datasets has been finished ===');

%% 2) Initialize output container
disp('=== Step 2: Preparing output container ===');
initialBufferCapacity = max(16, height(dataset_cpx));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

broadT_min_degreeC = 679;
broadT_max_degreeC = 1650;
broadP_min_kbar = 0;
broadP_max_kbar = 40;

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
    disp('=== Pressures were calculated: ===');
    displayPressureResult('Equation 30', row.P30_kbar, ...
        selectedCode_cpx, selectedIdx_liq);
    displayPressureResult('Equation 31', row.P31_kbar, ...
        selectedCode_cpx, selectedIdx_liq);

    if ~modelCautionIssued
        fprintf(2, ...
            ['CAUTION: Putirka (2008) Equations (30) and (31) require an ' ...
             'equilibrium Clinopyroxene-Liquid pair. Equation (31) is the ' ...
             'global calibration and is assigned to the standard P_kbar ' ...
             'output. The models commonly return approximately 2-3 kbar for ' ...
             '1-atm experiments (Fig. 9).\n']);
        modelCautionIssued = true;
    end

    if any(temperatureOutsideBroadRange) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the broad approximately ' ...
             '679-1650 degreeC envelope represented by the Cpx-saturated ' ...
             'experimental database summarized by Putirka (2008, Figs. 8-9). ' ...
             '%d of %d finite point(s) are outside; finite input range = ' ...
             '%.4g-%.4g degreeC. This is not a strict rectangular limit.\n'], ...
            sum(temperatureOutsideBroadRange), numel(finiteTemperature), ...
            min(finiteTemperature), max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    printPressureRangeWarning(row.P30_kbar, 'Equation 30', ...
        selectedCode_cpx, selectedIdx_liq, broadP_min_kbar, broadP_max_kbar);
    printPressureRangeWarning(row.P31_kbar, 'Equation 31', ...
        selectedCode_cpx, selectedIdx_liq, broadP_min_kbar, broadP_max_kbar);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in Putirka (2008) barometer input(s) for ' ...
             '%s & Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated. H2O is excluded ' ...
             'from cationTotal_liq, and F and Cl do not enter the equations.\n'], ...
            char(string(selectedCode_cpx)), selectedIdx_liq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    printNonFinitePressureWarning(row.P30_kbar, 'Equation 30', ...
        selectedCode_cpx, selectedIdx_liq);
    printNonFinitePressureWarning(row.P31_kbar, 'Equation 31', ...
        selectedCode_cpx, selectedIdx_liq);
    printNegativePressureWarning(row.P30_kbar, 'Equation 30', ...
        selectedCode_cpx, selectedIdx_liq);
    printNegativePressureWarning(row.P31_kbar, 'Equation 31', ...
        selectedCode_cpx, selectedIdx_liq);

    invalidJdLog = ~isfinite(row.logArg_Jd_Liq) | row.logArg_Jd_Liq <= 0;
    if any(invalidJdLog)
        fprintf(2, ...
            ['WARNING: The Jd-Liquid logarithm argument is non-positive or ' ...
             'non-finite for %s & Liquid row %d (%d of %d points). ' ...
             'Corresponding pressures remain NaN or Inf.\n'], ...
            char(string(selectedCode_cpx)), selectedIdx_liq, ...
            sum(invalidJdLog), numel(invalidJdLog));
    end

    invalidEq30 = ~row.Eq30_domain_valid;
    if any(invalidEq30)
        fprintf(2, ...
            ['WARNING: Equation (30) is outside its mathematical domain for ' ...
             '%s & Liquid row %d at %d of %d point(s). It requires positive ' ...
             'K_Jd_Liq, XFeO_Liq, XMgO_Liq, XKO0.5_Liq, and XDiHd_Cpx.\n'], ...
            char(string(selectedCode_cpx)), selectedIdx_liq, ...
            sum(invalidEq30), numel(invalidEq30));
    end

    invalidEq31 = ~row.Eq31_domain_valid;
    if any(invalidEq31)
        fprintf(2, ...
            ['WARNING: Equation (31) is outside its mathematical domain for ' ...
             '%s & Liquid row %d at %d of %d point(s). It requires positive ' ...
             'K_Jd_Liq, XDiHd_Cpx, XEnFs_Cpx, and total XAl_Cpx.\n'], ...
            char(string(selectedCode_cpx)), selectedIdx_liq, ...
            sum(invalidEq31), numel(invalidEq31));
    end

    disp('--------------------------------------------------');
    userAction = questdlg( ...
        'Continue with another Clinopyroxene selection (same Liquid row)?', ...
        'Putirka2008CpxJdLiq', ...
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
disp('=== Putirka2008CpxJdLiq finished ===');
end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_cpx, data_liq, T_degreeC)
% findNaNInputs
% Return names of calculation inputs containing NaN. H2O is checked because
% it enters both equations, but H2O is excluded from liquid normalization.

maxNames = 40;
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
% Reject Inf and finite negative values. Zero and NaN are retained.

maxNames = 40;
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
    error(['Putirka2008CpxJdLiq: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative values are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end
end

function row = calcPressure(data_cpx, data_liq, T_degreeC, MWinfo)
% calcPressure
% Calculate Putirka (2008) Equations (30) and (31) for one selected
% Clinopyroxene-Liquid pair over a scalar or vector of temperatures.

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

cpx = prepareCpxRow(data_cpx, MWinfo);
liq = extractLiquidRow(data_liq, MWinfo);

logArgScalar = cpx.XJd ./ ...
    (liq.XNaO0_5 .* liq.XAlO1_5 .* (liq.XSiO2 .^ 2));
lnK_Jd_Liq_scalar = log(logArgScalar);
MgNum_Liq_scalar = liq.XMgO ./ (liq.XMgO + liq.XFeO);

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
MgNum_Liq = repmat(MgNum_Liq_scalar, nT, 1);
logArg_Jd_Liq = repmat(logArgScalar, nT, 1);
lnK_Jd_Liq = repmat(lnK_Jd_Liq_scalar, nT, 1);

% Putirka (2008), Equation (30), p. 91.
P30_kbar = ...
    -48.7 ...
    + 271 .* (T_K ./ 10000) ...
    + 32 .* (T_K ./ 10000) .* lnK_Jd_Liq ...
    - 8.2 .* log(XFeO_liq) ...
    + 4.6 .* log(XMgO_liq) ...
    - 0.96 .* log(XKO0_5_liq) ...
    - 2.2 .* log(XDiHd_cpx) ...
    - 31 .* MgNum_Liq ...
    + 56 .* (XNaO0_5_liq + XKO0_5_liq) ...
    + 0.76 .* H2O_liq;

% Putirka (2008), Equation (31), p. 91. Printed coefficients are retained
% rather than the rounded coefficients used by some later implementations.
P31_kbar = ...
    -40.73 ...
    + 358 .* (T_K ./ 10000) ...
    + 21.69 .* (T_K ./ 10000) .* lnK_Jd_Liq ...
    - 105.7 .* XCaO_liq ...
    - 165.5 .* ((XNaO0_5_liq + XKO0_5_liq) .^ 2) ...
    - 50.15 .* XSiO2_liq .* (XFeO_liq + XMgO_liq) ...
    - 3.178 .* log(XDiHd_cpx) ...
    - 2.205 .* log(XEnFs_cpx) ...
    + 0.864 .* log(XAl_cpx) ...
    + 0.3962 .* H2O_liq;

P30_GPa = P30_kbar ./ 10;
P31_GPa = P31_kbar ./ 10;

Eq30_domain_valid = ...
    isfinite(T_K) & isfinite(logArg_Jd_Liq) & logArg_Jd_Liq > 0 & ...
    isfinite(XFeO_liq) & XFeO_liq > 0 & ...
    isfinite(XMgO_liq) & XMgO_liq > 0 & ...
    isfinite(XKO0_5_liq) & XKO0_5_liq > 0 & ...
    isfinite(XDiHd_cpx) & XDiHd_cpx > 0;

Eq31_domain_valid = ...
    isfinite(T_K) & isfinite(logArg_Jd_Liq) & logArg_Jd_Liq > 0 & ...
    isfinite(XDiHd_cpx) & XDiHd_cpx > 0 & ...
    isfinite(XEnFs_cpx) & XEnFs_cpx > 0 & ...
    isfinite(XAl_cpx) & XAl_cpx > 0;

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
row.MgNum_liq = MgNum_Liq;

row.logArg_Jd_Liq = logArg_Jd_Liq;
row.lnK_Jd_Liq = lnK_Jd_Liq;
row.Eq30_domain_valid = Eq30_domain_valid;
row.Eq31_domain_valid = Eq31_domain_valid;

row.P30_kbar = P30_kbar;
row.P30_GPa = P30_GPa;
row.P31_kbar = P31_kbar;
row.P31_GPa = P31_GPa;

% Standard barometer output uses the global Equation (31).
row.P_kbar = P31_kbar;
row.P_GPa = P31_GPa;

% Explicit aliases for compatibility and clarity.
row.PEq30_kbar = P30_kbar;
row.PEq30_GPa = P30_GPa;
row.PEq31_kbar = P31_kbar;
row.PEq31_GPa = P31_GPa;

row.P30_SEE_calibration_kbar = repmat(1.6, nT, 1);
row.P30_SEE_allData_kbar = repmat(3.6, nT, 1);
row.P31_SEE_allData_kbar = repmat(2.9, nT, 1);
row.isWithinBroadTRange = isfinite(T_degreeC) & ...
    T_degreeC >= 679 & T_degreeC <= 1650;
row.P30_isWithinBroadPRange = isfinite(P30_kbar) & ...
    P30_kbar >= 0 & P30_kbar <= 40;
row.P31_isWithinBroadPRange = isfinite(P31_kbar) & ...
    P31_kbar >= 0 & P31_kbar <= 40;
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
Al2O3 = getMineralOxOptional(data_cpx, 'Al2O3', 0);
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
% fractions. H2O enters Equations (30) and (31) as wt% but is excluded from
% cationTotal_liq. F and Cl are retained only as diagnostics.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liq = struct();
liq.SiO2 = getLiquidOxOptional(data_liq, 'SiO2', 0);
liq.TiO2 = getLiquidOxOptional(data_liq, 'TiO2', 0);
liq.Al2O3 = getLiquidOxOptional(data_liq, 'Al2O3', 0);

feOColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeO');
if ~isempty(feOColumn)
    liq.FeO = toScalarDouble(data_liq.(feOColumn));
else
    feOtColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeOt');
    if isempty(feOtColumn)
        liq.FeO = 0;
    else
        liq.FeO = toScalarDouble(data_liq.(feOtColumn));
    end
end

liq.MnO = getLiquidOxOptional(data_liq, 'MnO', 0);
liq.MgO = getLiquidOxOptional(data_liq, 'MgO', 0);
liq.CaO = getLiquidOxOptional(data_liq, 'CaO', 0);
liq.Na2O = getLiquidOxOptional(data_liq, 'Na2O', 0);
liq.K2O = getLiquidOxOptional(data_liq, 'K2O', 0);
liq.V2O3 = getLiquidOxOptional(data_liq, 'V2O3', 0);
liq.Cr2O3 = getLiquidOxOptional(data_liq, 'Cr2O3', 0);
liq.NiO = getLiquidOxOptional(data_liq, 'NiO', 0);
liq.P2O5 = getLiquidOxOptional(data_liq, 'P2O5', 0);
liq.SO3 = getLiquidOxOptional(data_liq, 'SO3', 0);
liq.Fe2O3 = getLiquidOxOptional(data_liq, 'Fe2O3', 0);
liq.H2O = getLiquidOxOptional(data_liq, 'H2O', 0);
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

function displayPressureResult(modelName, values, selectedCode_cpx, selectedIdx_liq)
% displayPressureResult
% Print one value or the first-to-last pressure range.

if isscalar(values)
    disp([modelName ' | ' char(string(selectedCode_cpx)) ' & Liquid row ' ...
        num2str(selectedIdx_liq) ': P = ' num2str(values) ' kbar']);
else
    disp([modelName ' | ' char(string(selectedCode_cpx)) ' & Liquid row ' ...
        num2str(selectedIdx_liq) ': P = ' num2str(values(1)) ' to ' ...
        num2str(values(end)) ' kbar']);
end
end

function printPressureRangeWarning(values, modelName, selectedCode_cpx, ...
        selectedIdx_liq, minimumPressure, maximumPressure)
% printPressureRangeWarning
% Warn for finite pressures outside the broad 0-40 kbar database envelope.

finiteMask = isfinite(values);
outsideMask = finiteMask & ...
    (values < minimumPressure | values > maximumPressure);
if any(outsideMask)
    finiteValues = values(finiteMask);
    fprintf(2, ...
        ['WARNING: %s pressure is outside the broad 0-40 kbar envelope used ' ...
         'for Putirka (2008) Cpx-Liquid evaluation. %d of %d finite point(s) ' ...
         'are outside; finite calculated range = %.4g-%.4g kbar for %s & ' ...
         'Liquid row %d.\n'], ...
        modelName, sum(outsideMask), sum(finiteMask), ...
        min(finiteValues), max(finiteValues), ...
        char(string(selectedCode_cpx)), selectedIdx_liq);
end
end

function printNonFinitePressureWarning(values, modelName, selectedCode_cpx, selectedIdx_liq)
% printNonFinitePressureWarning
% Retain and report NaN or Inf pressure outputs.

invalidMask = ~isfinite(values);
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite %s pressure values were calculated for %s & ' ...
         'Liquid row %d (%d of %d points; NaN: %d, Inf: %d). Values remain ' ...
         'in the output table.\n'], ...
        modelName, char(string(selectedCode_cpx)), selectedIdx_liq, ...
        sum(invalidMask), numel(values), sum(isnan(values)), sum(isinf(values)));
end
end

function printNegativePressureWarning(values, modelName, selectedCode_cpx, selectedIdx_liq)
% printNegativePressureWarning
% Retain negative finite values for diagnostics and report them.

negativeMask = isfinite(values) & values < 0;
if any(negativeMask)
    fprintf(2, ...
        ['WARNING: Negative finite %s pressure was calculated for %s & ' ...
         'Liquid row %d (%d of %d points). Values were retained for ' ...
         'diagnostic purposes.\n'], ...
        modelName, char(string(selectedCode_cpx)), selectedIdx_liq, ...
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

function value = getLiquidOxOptional(data_tbl, oxide, missingDefault)
% getLiquidOxOptional
% Retrieve an optional scalar Liquid oxide value. An absent column uses the
% stated backward-compatible default; a present NaN remains NaN.

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
