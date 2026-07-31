function results = Putirka2008CpxH2Obaro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Pyroxene/Putirka2008CpxH2Obaro.m
% Compatibility target: MATLAB R2024b
%
% H2O-corrected Clinopyroxene barometer, Equation (32b)
% Putirka, K.D. (2008)
% Reviews in Mineralogy and Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Clinopyroxene analysis selected
% from rawdata_struct.Cpx with the H2O content of one equilibrium Liquid
% analysis read by liquid.readLiquidExcel. It calculates pressure using
% Equation (32b) of Putirka (2008; p. 91).
%
% Equation (32b) uses:
%   1) temperature,
%   2) Clinopyroxene composition and calculated M1-M2 site populations, and
%   3) H2O in the equilibrium liquid, expressed in wt%.
%
% The remaining liquid major-element composition does not enter Equation
% (32b). The Liquid row is nevertheless retained in the output so that the
% H2O source and available identifiers remain traceable.
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Clinopyroxene-Liquid pair, one
% output row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Putirka (2008) recalibrated the Nimis (1995) Cpx-only approach using
% experiments spanning approximately 0.001-80 kbar (p. 91). Equation (32b)
% removes the systematic H2O-related error retained by Equation (32a), but
% requires an estimate of the H2O content of the liquid in equilibrium with
% the selected Cpx.
%
% For the global test dataset at P < 40 kbar, Putirka (2008) reports:
%
%   Equation (32b): R2 = 0.91, SEE = +/-2.6 kbar, n = 1173
%
% The broad Cpx experimental database discussed around Equations (30)-(32)
% spans approximately 679-1650 degreeC. This temperature interval and the
% 0.001-80 kbar experimental pressure span are used here only as
% non-stopping warning envelopes; they are not strict rectangular
% Equation (32b)-specific calibration limits.
%
% Important application cautions:
%
%   1) The selected Cpx and the H2O estimate must represent the same
%      equilibrium liquid. H2O from an unrelated whole rock, melt inclusion,
%      or late-stage liquid may yield a geologically meaningless pressure.
%
%   2) Equation (32b) requires an independently estimated temperature.
%      Temperature uncertainty propagates into pressure through both T and
%      ln(T).
%
%   3) H2O is entered in wt%, exactly as used by the published regression.
%      H2O is not converted to a cation fraction.
%
%   4) Cpx cations and components are calculated on a six-oxygen basis.
%      A cation sum near four is an important normalization and analytical
%      quality check.
%
%   5) Equation (32b) requires XJd > 0 because ln(XJd) occurs explicitly.
%      It also requires physically meaningful Fe-Mg distributions between
%      M1 and M2 sites.
%
%   6) The M1-M2 Fe-Mg allocation follows the Nimis (1995) structural model
%      summarized for Equation (32b). The physically relevant quadratic
%      root is used, and invalid site populations are retained as NaN.
%
%   7) Passing the numerical range and domain checks does not prove
%      Cpx-Liquid equilibrium or guarantee applicability to a natural
%      sample.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 679-1650 degreeC,
%   2) finite calculated pressure is outside 0.001-80 kbar,
%   3) a required calculation input contains NaN,
%   4) XJd or an M1-M2 site quantity is outside its domain,
%   5) the six-oxygen Cpx cation sum differs noticeably from four, or
%   6) a calculated pressure is NaN, Inf, or negative.
%
% Finite inputs used in the calculation must be non-negative. NaN is
% allowed, retained, propagated, and reported. Inf and finite negative
% values are rejected.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST Cpx-table column is treated as the identifier displayed in the
% selection dialog. Cpx analyses are read as oxide wt.% using common names
% or names ending in "_value".
%
% Required Cpx oxides:
%   SiO2, MgO, CaO, and FeO or FeOt
%
% Optional Cpx oxides, treated as zero only when the column is absent:
%   TiO2, Al2O3, MnO, Na2O, K2O, Cr2O3
%
% A present NaN value is retained and never silently replaced by zero.
% FeO is preferred when present; FeOt is used only when FeO is absent.
%
% The Liquid dataset is read using:
%   [liqAll, metaLiq] = liquid.readLiquidExcel();
%
% Required Liquid variable:
%   H2O
%
% H2O may be named H2O, H2O_value, or H2OValue. It is interpreted as wt%.
% Unlike optional major-element columns in some other barometers, an absent
% H2O column causes an error rather than defaulting to zero.
%
% -------------------------------------------------------------------------
% CLINOPYROXENE COMPONENT AND SITE CALCULATIONS
%
% Cpx oxide wt.% values are normalized to six oxygens. Define:
%
%   XAlIV = max(2 - XSi, 0)
%   XAlVI = max(XAl - XAlIV, 0)
%   XFe3  = max(XNa + XAlIV - XAlVI - 2*XTi - XCr, 0)
%   XJd   = min(XAlVI, XNa)
%   XCaTs = max(XAlVI - XJd, 0)
%   XCaTi = max((XAlIV - XCaTs)/2, 0)
%   XCrCaTs = max(XCr/2, 0)
%   XDiHd = max(XCa - XCaTi - XCaTs - XCrCaTs, 0)
%   XEnFs = max((XFe + XMg - XDiHd)/2, 0)
%
% For the Nimis (1995) M1-M2 allocation:
%
%   CNM = XCa + XNa + XMn
%   R3+ = XAlVI + XTi + XCr + XFe3
%
%   KD(M1-M2 Fe-Mg) =
%       exp(0.238*R3+ + 0.289*CNM - 2.315)
%
% where:
%
%   KD = [XFe2+(M1)*XMg(M2)] / [XFe2+(M2)*XMg(M1)]
%
% Define:
%
%   Fe2_total = XFe - XFe3
%   M2_capacity = 1 - CNM
%   x = XFe2+(M2)
%
% Substitution of:
%
%   XFe2+(M1) = Fe2_total - x
%   XMg(M2)   = M2_capacity - x
%   XMg(M1)   = XMg - M2_capacity + x
%
% gives:
%
%   a*x^2 + b*x + c = 0
%
%   a = KD - 1
%   b = KD*(XMg - M2_capacity) + Fe2_total + M2_capacity
%   c = -Fe2_total*M2_capacity
%
% The physically admissible root is selected by requiring non-negative M1
% and M2 Fe-Mg populations and conservation of total Fe2+ and Mg. If a is
% numerically zero, the corresponding linear equation b*x+c=0 is used.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION
%
% Putirka (2008), Equation (32b), p. 91:
%
%   P(kbar) =
%       1458
%     + 0.197*T(K)
%     - 241*ln[T(K)]
%     + 0.453*H2O_Liq(wt%)
%     + 55.5*XAlVI_cpx
%     + 8.05*XFe_cpx
%     - 277*XK_cpx
%     + 18*XJd_cpx
%     + 44.1*XDiHd_cpx
%     + 2.2*ln[XJd_cpx]
%     - 27.7*(XAl_cpx)^2
%     + 97.3*(XFe_M2_cpx)^2
%     + 30.7*(XMg_M2_cpx)^2
%     - 27.6*(XDiHd_cpx)^2
%
% Natural logarithms are used. Temperature is in Kelvin, H2O is in wt%,
% and pressure is returned in kbar.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008CpxH2Obaro(rawdata_struct, T_degreeC)
%   results = Putirka2008CpxH2Obaro(..., 'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing the Cpx table described above
%   T_degreeC      : non-negative numeric scalar or vector in degreeC.
%                    NaN is allowed; Inf and finite negative values are not.
%
% Name-value option:
%   'LiquidRow'    : positive integer scalar selecting one Liquid row.
%                    Default [] uses row 1.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Cpx-Liquid(H2O) pair.
%

%% Input validation
if nargin < 2
    error('Putirka2008CpxH2Obaro requires (rawdata_struct, T_degreeC).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(T_degreeC) || isempty(T_degreeC) || ~isvector(T_degreeC) || ...
        any(isinf(T_degreeC(:))) || ...
        any(isfinite(T_degreeC(:)) & T_degreeC(:) < 0)
    error(['T_degreeC must be a non-negative numeric scalar or vector. ' ...
           'NaN is allowed and retained, but Inf and finite negative ' ...
           'values are prohibited.']);
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
            ['WARNING: Liquid dataset has multiple rows (%d). Row 1 will ' ...
             'be used. Specify ''LiquidRow'' to select another row.\n'], ...
            height(liqAll));
    end
else
    selectedIdx_liq = liquidRowOpt;
    if selectedIdx_liq > height(liqAll)
        error(['Requested LiquidRow (%d) exceeds the number of rows in ' ...
               'the selected Liquid dataset (%d).'], ...
            selectedIdx_liq, height(liqAll));
    end
end

selectedData_liq = liqAll(selectedIdx_liq, :);
H2O_liq_wt = getLiquidH2ORequired(selectedData_liq);

disp(['Liquid selected: Row ' num2str(selectedIdx_liq) ...
    ', H2O = ' num2str(H2O_liq_wt) ' wt%']);
disp('=== Preparing Clinopyroxene and Liquid datasets has been finished ===');

%% 2) Initialize output container
disp('=== Step 2: Preparing output container ===');

% Fixed-size buffer avoids changing the cell-array size inside the loop.
maxResultBlocks = max(1024, height(dataset_cpx));
resultBlocks = cell(maxResultBlocks, 1);
nResultBlocks = 0;

% Broad experimental envelopes used only for non-stopping warnings.
experimentalT_min_degreeC = 679;
experimentalT_max_degreeC = 1650;
experimentalP_min_kbar = 0.001;
experimentalP_max_kbar = 80;

temperatureOutsideExperimental = isfinite(T_degreeC) & ...
    (T_degreeC < experimentalT_min_degreeC | ...
     T_degreeC > experimentalT_max_degreeC);
temperatureWarningIssued = false;
modelCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Clinopyroxene) ===');

while true
    if nResultBlocks >= maxResultBlocks
        fprintf(2, ...
            ['WARNING: The fixed result-buffer limit of %d selections was ' ...
             'reached. Completed calculations will be returned without ' ...
             'enlarging the result array.\n'], ...
            maxResultBlocks);
        break;
    end

    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', ...
        'Please select the Clinopyroxene data you would like to use:', ...
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
        selectedData_cpx, H2O_liq_wt, T_degreeC, MWinfo);

    nRows = height(row);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_cpx', 'dataRow_liq'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    resultBlocks{nResultBlocks} = row;

    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_cpx)) ' & Liquid row ' ...
            num2str(selectedIdx_liq) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_cpx)) ' & Liquid row ' ...
            num2str(selectedIdx_liq) ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    if ~modelCautionIssued
        fprintf(2, ...
            ['CAUTION: Putirka (2008) Equation (32b) requires the H2O ' ...
             'content of the liquid in equilibrium with the selected Cpx ' ...
             'and an independently appropriate temperature. The selected ' ...
             'Liquid row contributes H2O only; its other major-element ' ...
             'values do not enter Equation (32b).\n']);
        modelCautionIssued = true;
    end

    if any(temperatureOutsideExperimental) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        if ~isempty(finiteTemperature)
            fprintf(2, ...
                ['WARNING: Input temperature is outside the approximate ' ...
                 '679-1650 degreeC span of the broad Cpx experimental ' ...
                 'database discussed by Putirka (2008). This is not a ' ...
                 'strict Equation (32b)-specific limit. %d of %d finite ' ...
                 'temperature point(s) are outside; input finite range = ' ...
                 '%.4g-%.4g degreeC.\n'], ...
                sum(temperatureOutsideExperimental), ...
                numel(finiteTemperature), ...
                min(finiteTemperature), ...
                max(finiteTemperature));
        end
        temperatureWarningIssued = true;
    end

    finitePressure = isfinite(row.P_kbar);
    pressureOutsideExperimental = finitePressure & ...
        (row.P_kbar < experimentalP_min_kbar | ...
         row.P_kbar > experimentalP_max_kbar);

    if any(pressureOutsideExperimental)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the approximately ' ...
             '0.001-80 kbar experimental span reported around Putirka ' ...
             '(2008) Equations (32a-b; p. 91). %d of %d finite pressure ' ...
             'point(s) are outside; calculated finite range = ' ...
             '%.4g-%.4g kbar for %s & Liquid row %d.\n'], ...
            sum(pressureOutsideExperimental), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq);
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in Equation (32b) input(s) for %s ' ...
             '& Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; the ' ...
             'calculated pressure may remain NaN.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    finiteCationSum = isfinite(row.cationSum_cpx);
    cationSumOutsideDiagnostic = finiteCationSum & ...
        abs(row.cationSum_cpx - 4) > 0.05;

    if any(cationSumOutsideDiagnostic)
        fprintf(2, ...
            ['WARNING: Cpx six-oxygen cation sum differs from 4 by more ' ...
             'than 0.05 for %s (cation sum = %.6g). The threshold is a ' ...
             'diagnostic, not a published calibration limit.\n'], ...
            char(string(selectedCode_cpx)), ...
            row.cationSum_cpx(find(finiteCationSum, 1, 'first')));
    end

    invalidSiteNames = findInvalidSiteQuantities(row);
    if ~isempty(invalidSiteNames)
        fprintf(2, ...
            ['WARNING: Equation (32b) component, logarithm, or M1-M2 site ' ...
             'quantity is outside its valid domain for %s & Liquid row ' ...
             '%d: %s.\n' ...
             '         The affected values remain NaN or are retained for ' ...
             'diagnostic purposes.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            char(strjoin(invalidSiteNames, ', ')));
    end

    negativeDiscriminant = ...
        isfinite(row.M1M2_discriminant) & row.M1M2_discriminant < 0;
    if any(negativeDiscriminant)
        fprintf(2, ...
            ['WARNING: Negative M1-M2 quadratic discriminant was ' ...
             'calculated for %s & Liquid row %d. The corresponding site ' ...
             'populations and pressures remain NaN.\n'], ...
            char(string(selectedCode_cpx)), selectedIdx_liq);
    end

    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s ' ...
             '& Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & ' ...
             'Liquid row %d (%d of %d points). The values were retained ' ...
             'for diagnostic purposes.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Clinopyroxene selection (same Liquid row)?', ...
        'Putirka2008CpxH2Obaro', ...
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

results.Properties.UserData = struct( ...
    'liquid', metaLiq, ...
    'primaryPressureEquation', 'Putirka2008 Eq. (32b)', ...
    'pressureUnit', 'kbar', ...
    'temperatureUnitInEquation', 'K', ...
    'H2OUnit', 'wt%', ...
    'calibrationSEE_kbar', 2.6, ...
    'broadExperimentalPressureRange_kbar', [0.001 80]);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_cpx, data_liq, T_degreeC)
% findNaNInputs
% Return names of Equation (32b) inputs containing NaN. Only Cpx oxides,
% liquid H2O, and temperature are calculation inputs.

maxNames = 16;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

cpxOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO', ...
    'Na2O','K2O','Cr2O3'};

for i = 1:numel(cpxOxides)
    oxide = cpxOxides{i};
    columnName = findOxideColumn(data_cpx.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDouble(data_cpx.(columnName));
        if isnan(value)
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = "Cpx." + string(columnName);
        end
    end
end

feColumnName = findOxideColumn(data_cpx.Properties.VariableNames, 'FeO');
if isempty(feColumnName)
    feColumnName = findOxideColumn(data_cpx.Properties.VariableNames, 'FeOt');
end
if ~isempty(feColumnName)
    value = toScalarDouble(data_cpx.(feColumnName));
    if isnan(value)
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Cpx." + string(feColumnName);
    end
end

h2oColumnName = findOxideColumn(data_liq.Properties.VariableNames, 'H2O');
if ~isempty(h2oColumnName)
    value = toScalarDouble(data_liq.(h2oColumnName));
    if isnan(value)
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Liquid." + string(h2oColumnName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_cpx, data_liq)
% validateNonNegativeInputs
% Reject Inf and finite negative values used by Equation (32b). Zero and
% NaN are retained.

maxNames = 12;
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

cpxOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO', ...
    'Na2O','K2O','Cr2O3'};

for i = 1:numel(cpxOxides)
    oxide = cpxOxides{i};
    columnName = findOxideColumn(data_cpx.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDouble(data_cpx.(columnName));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = ...
                "Cpx." + string(columnName);
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
        invalidInputBuffer(nInvalidInputs) = ...
            "Cpx." + string(feColumnName);
    end
end

h2oColumnName = findOxideColumn(data_liq.Properties.VariableNames, 'H2O');
if ~isempty(h2oColumnName)
    value = toScalarDouble(data_liq.(h2oColumnName));
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Liquid." + string(h2oColumnName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['Putirka2008CpxH2Obaro: calculation inputs must be ' ...
           'non-negative. NaN is allowed, but Inf and finite negative ' ...
           'value(s) are prohibited. Invalid input(s): ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_cpx, H2O_liq_wt, T_degreeC, MWinfo)
% calcPressure
% Compute Putirka (2008) Equation (32b) for one Cpx row and one liquid-H2O
% value over one or more input temperatures.

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

cpx = prepareCpxRow(data_cpx, MWinfo);

% Natural logarithm of Jd. Non-positive values are outside the equation
% domain and remain NaN rather than producing a complex result.
lnXJd_scalar = NaN;
if isfinite(cpx.XJd) && cpx.XJd > 0
    lnXJd_scalar = log(cpx.XJd);
elseif isinf(cpx.XJd) && cpx.XJd > 0
    lnXJd_scalar = Inf;
end

% Nimis (1995) M1-M2 Fe-Mg site allocation.
CNM_scalar = cpx.XCa + cpx.XNa + cpx.XMn;
R3plus_scalar = cpx.XAlVI + cpx.XTi + cpx.XCr + cpx.XFe3;
KD_M1M2_scalar = exp( ...
    0.238 .* R3plus_scalar + 0.289 .* CNM_scalar - 2.315);

Fe2_total_scalar = cpx.XFe - cpx.XFe3;
M2_FeMg_capacity_scalar = 1 - CNM_scalar;

% Nimis (1995) defines:
%
%   KD = [Fe2+(M1)*Mg(M2)]/[Fe2+(M2)*Mg(M1)]
%
% Let x = Fe2+(M2). Conservation of Fe2+ and Mg gives the quadratic
% coefficients below.
a_scalar = KD_M1M2_scalar - 1;
b_scalar = ...
    KD_M1M2_scalar .* ...
        (cpx.XMg - M2_FeMg_capacity_scalar) ...
    + Fe2_total_scalar ...
    + M2_FeMg_capacity_scalar;
c_scalar = -Fe2_total_scalar .* M2_FeMg_capacity_scalar;

discriminant_scalar = b_scalar.^2 - 4 .* a_scalar .* c_scalar;
sqrtDiscriminant_scalar = NaN;
rootPlus_scalar = NaN;
rootMinus_scalar = NaN;
XFeM2_scalar = NaN;

if isfinite(discriminant_scalar) && discriminant_scalar >= 0
    sqrtDiscriminant_scalar = sqrt(discriminant_scalar);

    scale = max([1, abs(a_scalar), abs(b_scalar), abs(c_scalar)]);
    if abs(a_scalar) <= 100 .* eps(scale)
        if isfinite(b_scalar) && b_scalar ~= 0
            rootPlus_scalar = -c_scalar ./ b_scalar;
        end
    else
        rootPlus_scalar = ...
            (-b_scalar + sqrtDiscriminant_scalar) ./ (2 .* a_scalar);
        rootMinus_scalar = ...
            (-b_scalar - sqrtDiscriminant_scalar) ./ (2 .* a_scalar);
    end

    % Select the root that satisfies site capacities and cation
    % conservation. Usually only one root is physically admissible.
    candidateRoots = [rootPlus_scalar, rootMinus_scalar];
    candidateResiduals = [Inf, Inf];

    for iRoot = 1:numel(candidateRoots)
        candidate = candidateRoots(iRoot);
        if ~isfinite(candidate)
            continue;
        end

        candidateFeM2 = candidate;
        candidateMgM2 = M2_FeMg_capacity_scalar - candidateFeM2;
        candidateFeM1 = Fe2_total_scalar - candidateFeM2;
        candidateMgM1 = cpx.XMg - candidateMgM2;

        tolerance = 1e-10;
        isPhysical = ...
            candidateFeM2 >= -tolerance && ...
            candidateMgM2 >= -tolerance && ...
            candidateFeM1 >= -tolerance && ...
            candidateMgM1 >= -tolerance;

        if isPhysical
            % Evaluate the original Nimis exchange relationship. This also
            % provides deterministic selection if rounding makes both roots
            % appear admissible.
            denominator = candidateFeM2 .* candidateMgM1;
            if isfinite(denominator) && denominator > 0
                calculatedKD = ...
                    (candidateFeM1 .* candidateMgM2) ./ denominator;
                candidateResiduals(iRoot) = ...
                    abs(calculatedKD - KD_M1M2_scalar);
            elseif candidateFeM1 == 0 && candidateMgM2 == 0
                candidateResiduals(iRoot) = 0;
            end
        end
    end

    [minimumResidual, selectedRootIndex] = min(candidateResiduals);
    if isfinite(minimumResidual)
        XFeM2_scalar = candidateRoots(selectedRootIndex);
    end
elseif isinf(discriminant_scalar) && discriminant_scalar > 0
    sqrtDiscriminant_scalar = Inf;
end

XMgM2_scalar = M2_FeMg_capacity_scalar - XFeM2_scalar;
XFeM1_scalar = Fe2_total_scalar - XFeM2_scalar;
XMgM1_scalar = cpx.XMg - XMgM2_scalar;

% Expand composition-dependent scalars to the temperature-vector length.
H2O_liq = repmat(H2O_liq_wt, nT, 1);

SiO2_cpx = repmat(cpx.SiO2, nT, 1);
TiO2_cpx = repmat(cpx.TiO2, nT, 1);
Al2O3_cpx = repmat(cpx.Al2O3, nT, 1);
FeO_cpx = repmat(cpx.FeO, nT, 1);
MnO_cpx = repmat(cpx.MnO, nT, 1);
MgO_cpx = repmat(cpx.MgO, nT, 1);
CaO_cpx = repmat(cpx.CaO, nT, 1);
Na2O_cpx = repmat(cpx.Na2O, nT, 1);
K2O_cpx = repmat(cpx.K2O, nT, 1);
Cr2O3_cpx = repmat(cpx.Cr2O3, nT, 1);

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
lnXJd_cpx = repmat(lnXJd_scalar, nT, 1);

CNM_cpx = repmat(CNM_scalar, nT, 1);
R3plus_cpx = repmat(R3plus_scalar, nT, 1);
KD_M1M2_FeMg = repmat(KD_M1M2_scalar, nT, 1);
Fe2_total_cpx = repmat(Fe2_total_scalar, nT, 1);
M2_FeMg_capacity = repmat(M2_FeMg_capacity_scalar, nT, 1);
M1M2_a = repmat(a_scalar, nT, 1);
M1M2_b = repmat(b_scalar, nT, 1);
M1M2_c = repmat(c_scalar, nT, 1);
M1M2_discriminant = repmat(discriminant_scalar, nT, 1);
M1M2_sqrtDiscriminant = repmat(sqrtDiscriminant_scalar, nT, 1);
M1M2_rootPlus = repmat(rootPlus_scalar, nT, 1);
M1M2_rootMinus = repmat(rootMinus_scalar, nT, 1);
XFeM2_cpx = repmat(XFeM2_scalar, nT, 1);
XMgM2_cpx = repmat(XMgM2_scalar, nT, 1);
XFeM1_cpx = repmat(XFeM1_scalar, nT, 1);
XMgM1_cpx = repmat(XMgM1_scalar, nT, 1);

% Putirka (2008), Equation (32b), p. 91.
P_kbar = ...
    1458 ...
    + 0.197 .* T_K ...
    - 241 .* log(T_K) ...
    + 0.453 .* H2O_liq ...
    + 55.5 .* XAlVI_cpx ...
    + 8.05 .* XFe_cpx ...
    - 277 .* XK_cpx ...
    + 18 .* XJd_cpx ...
    + 44.1 .* XDiHd_cpx ...
    + 2.2 .* lnXJd_cpx ...
    - 27.7 .* (XAl_cpx .^ 2) ...
    + 97.3 .* (XFeM2_cpx .^ 2) ...
    + 30.7 .* (XMgM2_cpx .^ 2) ...
    - 27.6 .* (XDiHd_cpx .^ 2);

P_GPa = P_kbar ./ 10;

siteTolerance = 1e-10;
isWithinEquationDomain = ...
    isfinite(T_K) & T_K > 0 & ...
    isfinite(H2O_liq) & H2O_liq >= 0 & ...
    isfinite(XJd_cpx) & XJd_cpx > 0 & ...
    isfinite(lnXJd_cpx) & ...
    isfinite(CNM_cpx) & CNM_cpx >= 0 & CNM_cpx <= 1 + siteTolerance & ...
    isfinite(Fe2_total_cpx) & Fe2_total_cpx >= -siteTolerance & ...
    isfinite(M1M2_discriminant) & M1M2_discriminant >= 0 & ...
    isfinite(XFeM2_cpx) & XFeM2_cpx >= -siteTolerance & ...
    isfinite(XMgM2_cpx) & XMgM2_cpx >= -siteTolerance & ...
    isfinite(XFeM1_cpx) & XFeM1_cpx >= -siteTolerance & ...
    isfinite(XMgM1_cpx) & XMgM1_cpx >= -siteTolerance;

isWithinExperimentalTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 679 & T_degreeC <= 1650;

isWithinExperimentalPRange = ...
    isfinite(P_kbar) & P_kbar >= 0.001 & P_kbar <= 80;

isCationSumNearFour = ...
    isfinite(cationSum_cpx) & abs(cationSum_cpx - 4) <= 0.05;

row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.H2O_liq_wt = H2O_liq;
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.P_SEE_kbar = repmat(2.6, nT, 1);
row.P_SEE_GPa = repmat(0.26, nT, 1);

row.SiO2_cpx = SiO2_cpx;
row.TiO2_cpx = TiO2_cpx;
row.Al2O3_cpx = Al2O3_cpx;
row.FeO_cpx = FeO_cpx;
row.MnO_cpx = MnO_cpx;
row.MgO_cpx = MgO_cpx;
row.CaO_cpx = CaO_cpx;
row.Na2O_cpx = Na2O_cpx;
row.K2O_cpx = K2O_cpx;
row.Cr2O3_cpx = Cr2O3_cpx;

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
row.lnXJd_cpx = lnXJd_cpx;

row.CNM_cpx = CNM_cpx;
row.R3plus_cpx = R3plus_cpx;
row.KD_M1M2_FeMg = KD_M1M2_FeMg;
row.Fe2_total_cpx = Fe2_total_cpx;
row.M2_FeMg_capacity = M2_FeMg_capacity;
row.M1M2_a = M1M2_a;
row.M1M2_b = M1M2_b;
row.M1M2_c = M1M2_c;
row.M1M2_discriminant = M1M2_discriminant;
row.M1M2_sqrtDiscriminant = M1M2_sqrtDiscriminant;
row.M1M2_rootPlus = M1M2_rootPlus;
row.M1M2_rootMinus = M1M2_rootMinus;
row.XFeM1_cpx = XFeM1_cpx;
row.XMgM1_cpx = XMgM1_cpx;
row.XFeM2_cpx = XFeM2_cpx;
row.XMgM2_cpx = XMgM2_cpx;

row.isWithinEquationDomain = isWithinEquationDomain;
row.isWithinExperimentalTRange = isWithinExperimentalTRange;
row.isWithinExperimentalPRange = isWithinExperimentalPRange;
row.isCationSumNearFour = isCationSumNearFour;

end

function cpx = prepareCpxRow(data_cpx, MWinfo)
% prepareCpxRow
% Read one Cpx oxide row, normalize to six oxygens, and calculate the
% components used by Putirka (2008) Equation (32b).

if height(data_cpx) ~= 1
    error('Clinopyroxene input must be a 1-row table.');
end

cpx = struct();

cpx.SiO2 = getMineralOxRequired(data_cpx, 'SiO2', 'Clinopyroxene');
cpx.TiO2 = getMineralOxOptional(data_cpx, 'TiO2', 0);
cpx.Al2O3 = getMineralOxOptional(data_cpx, 'Al2O3', 0);
cpx.FeO = getFeORequired(data_cpx, 'Clinopyroxene');
cpx.MnO = getMineralOxOptional(data_cpx, 'MnO', 0);
cpx.MgO = getMineralOxRequired(data_cpx, 'MgO', 'Clinopyroxene');
cpx.CaO = getMineralOxRequired(data_cpx, 'CaO', 'Clinopyroxene');
cpx.Na2O = getMineralOxOptional(data_cpx, 'Na2O', 0);
cpx.K2O = getMineralOxOptional(data_cpx, 'K2O', 0);
cpx.Cr2O3 = getMineralOxOptional(data_cpx, 'Cr2O3', 0);

molProp = struct();
molProp.SiO2 = cpx.SiO2 ./ MWinfo.MW.SiO2;
molProp.TiO2 = cpx.TiO2 ./ MWinfo.MW.TiO2;
molProp.Al2O3 = cpx.Al2O3 ./ MWinfo.MW.Al2O3;
molProp.FeO = cpx.FeO ./ MWinfo.MW.FeO;
molProp.MnO = cpx.MnO ./ MWinfo.MW.MnO;
molProp.MgO = cpx.MgO ./ MWinfo.MW.MgO;
molProp.CaO = cpx.CaO ./ MWinfo.MW.CaO;
molProp.Na2O = cpx.Na2O ./ MWinfo.MW.Na2O;
molProp.K2O = cpx.K2O ./ MWinfo.MW.K2O;
molProp.Cr2O3 = cpx.Cr2O3 ./ MWinfo.MW.Cr2O3;

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

if isfinite(oxygenSum) && oxygenSum > 0
    oxygenRenormalizationFactor = 6 ./ oxygenSum;
else
    oxygenRenormalizationFactor = NaN;
end

cpx.oxygenSum = oxygenSum;
cpx.oxygenRenormalizationFactor = oxygenRenormalizationFactor;

cpx.XSi = molProp.SiO2 .* oxygenRenormalizationFactor;
cpx.XTi = molProp.TiO2 .* oxygenRenormalizationFactor;
cpx.XAl = 2 .* molProp.Al2O3 .* oxygenRenormalizationFactor;
cpx.XFe = molProp.FeO .* oxygenRenormalizationFactor;
cpx.XMn = molProp.MnO .* oxygenRenormalizationFactor;
cpx.XMg = molProp.MgO .* oxygenRenormalizationFactor;
cpx.XCa = molProp.CaO .* oxygenRenormalizationFactor;
cpx.XNa = 2 .* molProp.Na2O .* oxygenRenormalizationFactor;
cpx.XK = 2 .* molProp.K2O .* oxygenRenormalizationFactor;
cpx.XCr = 2 .* molProp.Cr2O3 .* oxygenRenormalizationFactor;

cpx.cationSum = ...
    cpx.XSi + cpx.XTi + cpx.XAl + cpx.XFe + cpx.XMn + ...
    cpx.XMg + cpx.XCa + cpx.XNa + cpx.XK + cpx.XCr;

cpx.XAlIV = maxPreserveNaN(2 - cpx.XSi, 0);
cpx.XAlVI = maxPreserveNaN(cpx.XAl - cpx.XAlIV, 0);

cpx.XFe3 = ...
    cpx.XNa + cpx.XAlIV - cpx.XAlVI - 2 .* cpx.XTi - cpx.XCr;
cpx.XFe3 = maxPreserveNaN(cpx.XFe3, 0);

cpx.XJd = minPreserveNaN(cpx.XAlVI, cpx.XNa);
cpx.XJd = maxPreserveNaN(cpx.XJd, 0);

cpx.XCaTs = maxPreserveNaN(cpx.XAlVI - cpx.XJd, 0);
cpx.XCaTi = maxPreserveNaN((cpx.XAlIV - cpx.XCaTs) ./ 2, 0);
cpx.XCrCaTs = maxPreserveNaN(cpx.XCr ./ 2, 0);
cpx.XDiHd = maxPreserveNaN( ...
    cpx.XCa - cpx.XCaTi - cpx.XCaTs - cpx.XCrCaTs, 0);
cpx.XEnFs = maxPreserveNaN( ...
    (cpx.XFe + cpx.XMg - cpx.XDiHd) ./ 2, 0);

end

function invalidSiteNames = findInvalidSiteQuantities(row)
% findInvalidSiteQuantities
% Return names of composition-dependent terms outside the Equation (32b)
% domain. The first row is sufficient because these values do not vary with
% the input temperature vector.

maxNames = 15;
invalidBuffer = strings(maxNames, 1);
nInvalid = 0;

if isempty(row)
    invalidSiteNames = strings(0, 1);
    return;
end

if ~isfinite(row.XJd_cpx(1)) || row.XJd_cpx(1) <= 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XJd_cpx";
end
if ~isfinite(row.lnXJd_cpx(1))
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "lnXJd_cpx";
end
if ~isfinite(row.XAlVI_cpx(1)) || row.XAlVI_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XAlVI_cpx";
end
if ~isfinite(row.XDiHd_cpx(1)) || row.XDiHd_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XDiHd_cpx";
end
if ~isfinite(row.CNM_cpx(1)) || row.CNM_cpx(1) < 0 || row.CNM_cpx(1) > 1
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "CNM_cpx";
end
if ~isfinite(row.R3plus_cpx(1)) || row.R3plus_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "R3plus_cpx";
end
if ~isfinite(row.KD_M1M2_FeMg(1)) || row.KD_M1M2_FeMg(1) <= 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "KD_M1M2_FeMg";
end
if ~isfinite(row.Fe2_total_cpx(1)) || row.Fe2_total_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "Fe2_total_cpx";
end
if ~isfinite(row.M2_FeMg_capacity(1)) || row.M2_FeMg_capacity(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "M2_FeMg_capacity";
end
if ~isfinite(row.M1M2_discriminant(1)) || row.M1M2_discriminant(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "M1M2_discriminant";
end
if ~isfinite(row.XFeM2_cpx(1)) || row.XFeM2_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XFeM2_cpx";
end
if ~isfinite(row.XMgM2_cpx(1)) || row.XMgM2_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XMgM2_cpx";
end
if ~isfinite(row.XFeM1_cpx(1)) || row.XFeM1_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XFeM1_cpx";
end
if ~isfinite(row.XMgM1_cpx(1)) || row.XMgM1_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XMgM1_cpx";
end
if ~isfinite(row.H2O_liq_wt(1)) || row.H2O_liq_wt(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "H2O_liq_wt";
end

invalidSiteNames = invalidBuffer(1:nInvalid);

end

function H2O = getLiquidH2ORequired(data_liq)
% getLiquidH2ORequired
% Read H2O in wt%. The H2O column is mandatory for Equation (32b).

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

columnName = findOxideColumn(data_liq.Properties.VariableNames, 'H2O');
if isempty(columnName)
    error(['Liquid table must contain H2O for Putirka (2008) Equation ' ...
           '(32b). Recognized forms include H2O, H2O_value, and H2OValue.']);
end

H2O = toScalarDouble(data_liq.(columnName));

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Copy common Liquid identifiers to every temperature row.

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
% Retrieve one required scalar oxide value. A present NaN is retained.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    error('%s table must contain variable: %s', mineralLabel, oxide);
end
value = toScalarDouble(data_tbl.(columnName));

end

function value = getMineralOxOptional(data_tbl, oxide, missingDefault)
% getMineralOxOptional
% Retrieve one optional scalar oxide value. An absent column uses the
% stated default; a present NaN remains NaN.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    value = missingDefault;
else
    value = toScalarDouble(data_tbl.(columnName));
end

end

function value = getFeORequired(data_tbl, mineralLabel)
% getFeORequired
% Prefer FeO when present. Use FeOt only when FeO is absent.

feOColumn = findOxideColumn(data_tbl.Properties.VariableNames, 'FeO');
if ~isempty(feOColumn)
    value = toScalarDouble(data_tbl.(feOColumn));
    return;
end

feOtColumn = findOxideColumn(data_tbl.Properties.VariableNames, 'FeOt');
if ~isempty(feOtColumn)
    value = toScalarDouble(data_tbl.(feOtColumn));
    return;
end

error('%s table must contain FeO or FeOt.', mineralLabel);

end

function columnName = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match oxide names after removing spaces, underscores, and hyphens.
% Both "Oxide" and "Oxide_value" forms are recognized.

canonicalNames = cell(size(variableNames));
for i = 1:numel(variableNames)
    canonicalNames{i} = canonicalizeName(variableNames{i});
end

canonicalOxide = canonicalizeName(oxide);
targetNames = {[canonicalOxide 'value'], canonicalOxide};

columnName = '';
for i = 1:numel(targetNames)
    matchedIndex = find(strcmp(canonicalNames, targetNames{i}), 1, 'first');
    if ~isempty(matchedIndex)
        columnName = variableNames{matchedIndex};
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

function value = toScalarDouble(rawValue)
% toScalarDouble
% Convert a one-row table value to a scalar double. Missing, empty, or
% non-numeric text becomes NaN and is not converted to zero.

value = NaN;

if isempty(rawValue)
    return;
end

if isnumeric(rawValue) || islogical(rawValue)
    if numel(rawValue) ~= 1
        error('Expected one scalar value.');
    end
    value = double(rawValue);
    return;
end

if isstring(rawValue)
    if numel(rawValue) ~= 1 || ismissing(rawValue)
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
        return;
    end
    value = str2double(string(rawValue));
    return;
end

if iscell(rawValue)
    if numel(rawValue) ~= 1 || isempty(rawValue{1})
        return;
    end
    value = toScalarDouble(rawValue{1});
    return;
end

error('Table value must be numeric or convertible to a numeric scalar.');

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
