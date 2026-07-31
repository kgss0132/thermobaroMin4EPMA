function results = Putirka1996baro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Pyroxene/Putirka1996baro.m
% Tested with MATLAB R2024b
%
% Clinopyroxene-Liquid barometers of Putirka et al. (1996)
% Putirka, K., Johnson, M., Kinzler, R., Longhi, J. and Walker, D. (1996)
% Contributions to Mineralogy and Petrology, 123, 92-108
% DOI: https://doi.org/10.1007/s004100050145
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Clinopyroxene analysis and pairs
% it with one row from a separately loaded Liquid dataset. It calculates
% pressure using the P1 and P2 clinopyroxene-liquid barometers presented in
% Table 5 of Putirka et al. (1996). P1 is the preferred pressure model of
% the original paper (p. 100).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Clinopyroxene-Liquid pair, one
% output row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Putirka et al. (1996) calibrated clinopyroxene-liquid thermobarometers
% using new experiments and previously published equilibrium experiments.
% The new experiments span approximately (abstract, p. 92; Table 2 and
% experimental description, pp. 95-96):
%
%   Temperature : 1100-1475 degreeC
%   Pressure    : 8-30 kbar
%   Materials   : mafic liquids including MORB-like basalt, Woodlark Basin
%                 basalt, ankaramite, ugandite, and synthetic SCAM+Na data
%
% These limits are used in this implementation as non-stopping warning
% ranges. Although the paper title refers to 0-30 kbar, 1-bar experiments
% were not included in the barometer regression and were used only as an
% extrapolation test (p. 102). Results below the main 8-kbar calibration
% range must therefore be treated as extrapolated values.
%
% Model P1 is preferred for pressure-temperature estimates of natural
% samples and has a regression standard error of estimate of approximately
% 1.36 kbar; model P2 has an SEE of approximately 1.51 kbar (pp. 100-101).
% The Summary states that individual equilibrium Clinopyroxene-Liquid pairs
% are accurate to approximately +/-1.4 kbar, with potentially improved
% precision when multiple equilibrium pairs are averaged (p. 107).
%
% Clinopyroxene and liquid must represent an equilibrium pair. The authors
% identify disequilibrium clinopyroxene compositions as the most important
% source of experimental pressure error (pp. 101-102). Zoned crystals,
% inherited crystals, reaction rims, mixed magma populations, or analyses
% that do not represent the same crystallization stage may therefore produce
% misleading pressures.
%
% When whole-rock composition is used as a liquid proxy, it must represent
% the melt from which the selected clinopyroxene crystallized. Putirka et al.
% (1996) explicitly tested and, where necessary, corrected whole-rock
% compositions before applying the barometer to Hawaiian samples
% (pp. 105-106). Crystal accumulation, melt mixing, and open-system
% evolution can invalidate a whole-rock liquid proxy.
%
% Additional cautions discussed in the original paper include:
%   1) At 1 bar, individual estimates show large scatter; Na loss, sluggish
%      equilibration below about 1050 degreeC, and sector zoning can cause
%      erroneous pressures (p. 102).
%   2) Extremely FeO-rich and SiO2-poor liquids (approximately 27-33 wt%%
%      FeO and 34-36 wt%% SiO2) yielded erroneously high pressures in test
%      data and may lie outside the compositional applicability of the model
%      (p. 101).
%   3) A single 50-kbar test experiment was reproduced successfully, but
%      this does not constitute a broadly calibrated >30-kbar range
%      (pp. 100-101).
%   4) The experimental clinopyroxenes were produced in graphite capsules;
%      the acmite component was considered negligible. Strongly oxidized
%      clinopyroxenes with substantial Fe3+ were not explicitly calibrated
%      by the component scheme (p. 99).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 1100-1475 degreeC,
%   2) finite calculated pressure is outside 8-30 kbar,
%   3) the selected liquid is extremely FeO-rich or SiO2-poor relative to
%      the problematic test compositions discussed on p. 101,
%   4) a required calculation input contains NaN, or
%   5) a calculated pressure is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table containing clinopyroxene oxide analyses
%
% The FIRST column is treated as an identifier ("data code") displayed in
% the selection dialog. Oxide columns are located case-insensitively after
% removing spaces, underscores, and hyphens. Either "SiO2" or
% "SiO2_value"-style names are accepted.
%
% Required Clinopyroxene oxide variables:
%   SiO2, MgO, CaO, and either FeO or FeOt
%
% Optional Clinopyroxene oxide variables:
%   TiO2, Al2O3, MnO, Na2O, K2O, Cr2O3
%
% A Liquid table is loaded with liquid.readLiquidExcel(). The barometer uses
% liquid cation fractions calculated from the following oxides when present:
%   SiO2, TiO2, Al2O3, FeO, MnO, MgO, CaO, Na2O, K2O,
%   V2O3, Cr2O3, NiO, P2O5, SO3, Fe2O3
%
% F and Cl are intentionally excluded from cationTotal_liq and are also
% excluded from the NaN-input warning list. Their values are retained in the
% output table for reference when present.
%
% Missing optional oxide columns are interpreted as zero. In contrast, NaN
% values in existing input columns are never replaced by zero. They are
% retained, propagated through the calculation, and listed in a non-stopping
% fprintf warning. Inf values and finite negative values in calculation
% inputs are prohibited. Zero is retained; if it makes a logarithm or ratio
% undefined, the resulting NaN or Inf is retained and reported.
%
% -------------------------------------------------------------------------
% CLINOPYROXENE COMPONENT CALCULATION
%
% Clinopyroxene cations are normalized to 6 oxygens. Components follow the
% normative scheme described by Putirka et al. (1996, pp. 97-99):
%
%   AlIV = 2 - Si
%   AlVI = Al(total) - AlIV
%   Jd   = min(Na, AlVI)
%   CaTs = AlVI - Jd
%   CaTi = (AlIV - CaTs)/2, where AlIV exceeds CaTs
%   DiHd = Ca - CaTs - CaTi
%   EnFs = [Fe + Mg - DiHd]/2
%
% Derived negative component amounts are bounded at zero in accordance with
% the normative allocation scheme. NaN values remain NaN.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATIONS (Table 5, p. 100)
%
% Model P1, preferred:
%
%   P1(kbar) = -54.3
%              + 299*T(K)/1e4
%              + 36.4*T(K)/1e4
%                * ln[Jd_cpx / (Si_liq^2*Na_liq*Al_liq)]
%              + 367*(Na_liq*Al_liq)
%
% Model P2:
%
%   P2(kbar) = -50.7
%              + 394*T(K)/1e4
%              + 36.4*T(K)/1e4
%                * ln[Jd_cpx / (Si_liq^2*Na_liq*Al_liq)]
%              - 20.0*T(K)/1e4*ln[1/(Na_liq*Al_liq)]
%
% Liquid quantities are cation fractions. Natural logarithms are used.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka1996baro(rawdata_struct, T_degreeC)
%   results = Putirka1996baro(..., 'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing a Cpx table
%   T_degreeC      : non-negative numeric scalar or vector in degreeC.
%                    NaN is allowed and retained; Inf is prohibited.
%
% Name-value option:
%   'LiquidRow'    : positive integer row number in the loaded Liquid table.
%                    Default is the first Liquid row.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Clinopyroxene-Liquid pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('Putirka1996baro requires (rawdata_struct, T_degreeC).');
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
if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

T_degreeC = T_degreeC(:);

%% Options
ip = inputParser;
ip.FunctionName = 'Putirka1996baro';
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOption = ip.Results.LiquidRow;

%% 1) Retrieve Clinopyroxene and Liquid datasets
disp('=== Step 1: Preparing Clinopyroxene and Liquid datasets ===');

dataset_cpx = rawdata_struct.Cpx;
MWinfo = liquid.getMolarWeights();

[liquidAll, metaLiquid] = liquid.readLiquidExcel();
if isempty(liquidAll) || ~istable(liquidAll)
    error('Selected Liquid dataset is empty or is not a table.');
end

if isempty(liquidRowOption)
    selectedLiquidRow = 1;
    if height(liquidAll) > 1
        fprintf(2, ...
            ['WARNING: Liquid dataset contains %d rows. The first row is being used. ' ...
             'Specify ''LiquidRow'' to select another row.\n'], ...
            height(liquidAll));
    end
else
    selectedLiquidRow = liquidRowOption;
    if selectedLiquidRow > height(liquidAll)
        error('Requested LiquidRow (%d) exceeds Liquid table height (%d).', ...
            selectedLiquidRow, height(liquidAll));
    end
end

disp('=== Preparing Clinopyroxene and Liquid datasets has been finished ===');

%% 2) Initialize output container and warning ranges
% Store each selected-pair result in a cell buffer and concatenate only once
% after the interactive loop. This avoids repeated growth of the output
% table at every loop iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degreeC = 1100;
calibrationT_max_degreeC = 1475;
calibrationP_min_kbar = 8;
calibrationP_max_kbar = 30;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);
temperatureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive Clinopyroxene selection + calculation
disp('=== Step 3: Selecting a data code from the list (Clinopyroxene) ===');

while true
    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIndex_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Clinopyroxene data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIndex_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIndex_cpx);
    selectedData_cpx = dataset_cpx(selectedIndex_cpx, :);
    selectedData_liq = liquidAll(selectedLiquidRow, :);

    disp(['Clinopyroxene selected: ' char(string(selectedCode_cpx))]);
    disp(['Liquid selected: Row ' num2str(selectedLiquidRow)]);
    disp('=== Step 4: Preparing selected composition data ===');

    % Read selected compositions once. Missing optional columns become zero,
    % whereas NaN values in existing columns remain NaN.
    cpxOxides = prepareCpxOxides(selectedData_cpx);
    liquidOxides = prepareLiquidOxides(selectedData_liq);

    % List NaN values only for variables that contribute to the calculation.
    % Liquid F and Cl are deliberately excluded.
    nanInputNames = findNaNInputs(cpxOxides, liquidOxides, T_degreeC);

    % Reject Inf and finite negative inputs. NaN and zero remain unchanged.
    validateNonNegativeInputs(cpxOxides, liquidOxides);

    disp('=== Step 5: Calculating pressure ===');
    row = calcPressure(cpxOxides, liquidOxides, T_degreeC, MWinfo);

    % Repeat identifiers for every input temperature value.
    nRows = height(row);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), nRows, 1);
    row.dataRow_liq = repmat(selectedLiquidRow, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_cpx', 'dataRow_liq'}, 'Before', 1);

    % Store one block per selected pair. Expand capacity only when needed,
    % rather than changing the result size on every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo calculated values for immediate inspection.
    disp('--------------------------------------------------');
    displayPressureSummary(selectedCode_cpx, selectedLiquidRow, row);

    % Input temperature is common to all selected Clinopyroxene rows, so the
    % warning is issued only once per function call.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the 1100-1475 degreeC range ' ...
             'of the new Putirka et al. (1996) experiments (p. 92; pp. 95-96). ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn for each model when finite pressures lie outside the main
    % 8-30 kbar calibration envelope of the new experiments.
    printPressureRangeWarning(row.P1_kbar, 'P1', ...
        calibrationP_min_kbar, calibrationP_max_kbar, ...
        selectedCode_cpx, selectedLiquidRow);
    printPressureRangeWarning(row.P2_kbar, 'P2', ...
        calibrationP_min_kbar, calibrationP_max_kbar, ...
        selectedCode_cpx, selectedLiquidRow);

    % Problematic liquid-composition tests discussed on p. 101.
    if isfinite(liquidOxides.FeO) && liquidOxides.FeO >= 27
        fprintf(2, ...
            ['WARNING: Liquid FeO = %.4g wt%% is within or above the approximately ' ...
             '27-33 wt%% FeO range for which Putirka et al. (1996) reported ' ...
             'erroneously high pressures in test compositions (p. 101). ' ...
             'Clinopyroxene %s; Liquid row %d.\n'], ...
            liquidOxides.FeO, char(string(selectedCode_cpx)), selectedLiquidRow);
    end

    if isfinite(liquidOxides.SiO2) && liquidOxides.SiO2 <= 36
        fprintf(2, ...
            ['WARNING: Liquid SiO2 = %.4g wt%% is within or below the approximately ' ...
             '34-36 wt%% SiO2 range for which Putirka et al. (1996) reported ' ...
             'erroneously high pressures in test compositions (p. 101). ' ...
             'Clinopyroxene %s; Liquid row %d.\n'], ...
            liquidOxides.SiO2, char(string(selectedCode_cpx)), selectedLiquidRow);
    end

    % Report exact input names containing NaN. Values are retained and
    % propagated rather than replaced by zero.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for ' ...
             'Clinopyroxene %s and Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; calculated ' ...
             'pressure may remain NaN. Liquid F and Cl are excluded from this ' ...
             'warning because they are not used in cationTotal_liq.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedLiquidRow, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report non-finite results from each model.
    printNonFinitePressureWarning(row.P1_kbar, 'P1', ...
        selectedCode_cpx, selectedLiquidRow);
    printNonFinitePressureWarning(row.P2_kbar, 'P2', ...
        selectedCode_cpx, selectedLiquidRow);

    % Negative finite pressure is retained for diagnostics but is outside the
    % physical and calibrated range.
    printNegativePressureWarning(row.P1_kbar, 'P1', ...
        selectedCode_cpx, selectedLiquidRow);
    printNegativePressureWarning(row.P2_kbar, 'P2', ...
        selectedCode_cpx, selectedLiquidRow);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Clinopyroxene selection (same Liquid row)?', ...
        'Putirka1996baro', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered result blocks once. If selection was canceled before
% any calculation, return an empty table.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct('liquid', metaLiquid);
disp('=== Putirka1996baro calculation has been finished! ===');

end

%% ---- local functions ----
function cpxOxides = prepareCpxOxides(data_cpx)
% prepareCpxOxides
% Extract one Clinopyroxene oxide analysis. Missing optional variables are
% represented by zero, while NaN in an existing variable is preserved.

if height(data_cpx) ~= 1
    error('Clinopyroxene input must be a 1-row table.');
end

cpxOxides = struct();
cpxOxides.SiO2 = getMineralOxideRequired(data_cpx, 'SiO2');
cpxOxides.TiO2 = getMineralOxideOptional(data_cpx, 'TiO2', 0);
cpxOxides.Al2O3 = getMineralOxideOptional(data_cpx, 'Al2O3', 0);

% Use FeO when the column exists. Use FeOt only when FeO is absent. A NaN
% FeO value is not silently replaced by FeOt.
feOColumn = findOxideColumn(data_cpx.Properties.VariableNames, 'FeO');
if ~isempty(feOColumn)
    cpxOxides.FeO = toScalarDouble(data_cpx.(feOColumn));
    cpxOxides.FeSource = "FeO";
else
    feOtColumn = findOxideColumn(data_cpx.Properties.VariableNames, 'FeOt');
    if isempty(feOtColumn)
        error('Selected Clinopyroxene row must contain FeO or FeOt.');
    end
    cpxOxides.FeO = toScalarDouble(data_cpx.(feOtColumn));
    cpxOxides.FeSource = "FeOt";
end

cpxOxides.MnO = getMineralOxideOptional(data_cpx, 'MnO', 0);
cpxOxides.MgO = getMineralOxideRequired(data_cpx, 'MgO');
cpxOxides.CaO = getMineralOxideRequired(data_cpx, 'CaO');
cpxOxides.Na2O = getMineralOxideOptional(data_cpx, 'Na2O', 0);
cpxOxides.K2O = getMineralOxideOptional(data_cpx, 'K2O', 0);
cpxOxides.Cr2O3 = getMineralOxideOptional(data_cpx, 'Cr2O3', 0);

end

function liquidOxides = prepareLiquidOxides(data_liq)
% prepareLiquidOxides
% Extract one Liquid oxide analysis. Missing optional variables are zero;
% NaN in an existing variable is retained. F and Cl are retained only as
% reference outputs and are not used in cationTotal_liq.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liquidOxides = struct();
liquidOxides.SiO2 = getLiquidOxideOptional(data_liq, 'SiO2', 0);
liquidOxides.TiO2 = getLiquidOxideOptional(data_liq, 'TiO2', 0);
liquidOxides.Al2O3 = getLiquidOxideOptional(data_liq, 'Al2O3', 0);
liquidOxides.FeO = getLiquidOxideOptional(data_liq, 'FeO', 0);
liquidOxides.MnO = getLiquidOxideOptional(data_liq, 'MnO', 0);
liquidOxides.MgO = getLiquidOxideOptional(data_liq, 'MgO', 0);
liquidOxides.CaO = getLiquidOxideOptional(data_liq, 'CaO', 0);
liquidOxides.Na2O = getLiquidOxideOptional(data_liq, 'Na2O', 0);
liquidOxides.K2O = getLiquidOxideOptional(data_liq, 'K2O', 0);
liquidOxides.V2O3 = getLiquidOxideOptional(data_liq, 'V2O3', 0);
liquidOxides.Cr2O3 = getLiquidOxideOptional(data_liq, 'Cr2O3', 0);
liquidOxides.NiO = getLiquidOxideOptional(data_liq, 'NiO', 0);
liquidOxides.P2O5 = getLiquidOxideOptional(data_liq, 'P2O5', 0);
liquidOxides.SO3 = getLiquidOxideOptional(data_liq, 'SO3', 0);
liquidOxides.Fe2O3 = getLiquidOxideOptional(data_liq, 'Fe2O3', 0);

% F and Cl are intentionally excluded from calculation, validation, and NaN
% warnings, but are retained in output for reference.
liquidOxides.F = getLiquidOxideOptional(data_liq, 'F', 0);
liquidOxides.Cl = getLiquidOxideOptional(data_liq, 'Cl', 0);

end

function nanInputNames = findNaNInputs(cpxOxides, liquidOxides, T_degreeC)
% findNaNInputs
% Return names of calculation inputs containing NaN. Liquid F and Cl are
% intentionally excluded because they do not contribute to cationTotal_liq.

maxNames = 1 + 10 + 15;
nanInputBuffer = strings(maxNames, 1);
nNaNInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNaNInputs = nNaNInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNaNInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

cpxFields = {'SiO2','TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O','Cr2O3'};
for i = 1:numel(cpxFields)
    fieldName = cpxFields{i};
    if isnan(cpxOxides.(fieldName))
        nNaNInputs = nNaNInputs + 1;
        if strcmp(fieldName, 'FeO')
            nanInputBuffer(nNaNInputs) = ...
                "Cpx." + cpxOxides.FeSource;
        else
            nanInputBuffer(nNaNInputs) = "Cpx." + string(fieldName);
        end
    end
end

liquidFields = {'SiO2','TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O', ...
    'V2O3','Cr2O3','NiO','P2O5','SO3','Fe2O3'};
for i = 1:numel(liquidFields)
    fieldName = liquidFields{i};
    if isnan(liquidOxides.(fieldName))
        nNaNInputs = nNaNInputs + 1;
        nanInputBuffer(nNaNInputs) = "Liquid." + string(fieldName);
    end
end

nanInputNames = nanInputBuffer(1:nNaNInputs);

end

function validateNonNegativeInputs(cpxOxides, liquidOxides)
% validateNonNegativeInputs
% Reject Inf and finite negative values in calculation inputs. NaN and zero
% are intentionally allowed. Liquid F and Cl are not calculation inputs and
% are excluded from this validation.

maxNames = 10 + 15;
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

cpxFields = {'SiO2','TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O','Cr2O3'};
for i = 1:numel(cpxFields)
    fieldName = cpxFields{i};
    value = cpxOxides.(fieldName);
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        if strcmp(fieldName, 'FeO')
            invalidInputBuffer(nInvalidInputs) = ...
                "Cpx." + cpxOxides.FeSource;
        else
            invalidInputBuffer(nInvalidInputs) = "Cpx." + string(fieldName);
        end
    end
end

liquidFields = {'SiO2','TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O', ...
    'V2O3','Cr2O3','NiO','P2O5','SO3','Fe2O3'};
for i = 1:numel(liquidFields)
    fieldName = liquidFields{i};
    value = liquidOxides.(fieldName);
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = "Liquid." + string(fieldName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['Putirka1996baro: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(cpxOxides, liquidOxides, T_degreeC, MWinfo)
% calcPressure
% Calculate P1 and P2 for one Clinopyroxene-Liquid pair at one or more input
% temperatures. NaN and zero are not replaced; they propagate naturally.

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% ---- Clinopyroxene cations and components
cpx = calculateCpxComponents(cpxOxides, MWinfo);

% ---- Liquid cation amounts
% F and Cl are deliberately absent from this calculation.
n = struct();
n.SiO2 = liquidOxides.SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
n.TiO2 = liquidOxides.TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
n.Al2O3 = liquidOxides.Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n.FeO = liquidOxides.FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
n.MnO = liquidOxides.MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
n.MgO = liquidOxides.MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
n.CaO = liquidOxides.CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
n.Na2O = liquidOxides.Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
n.K2O = liquidOxides.K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
n.V2O3 = liquidOxides.V2O3 .* MWinfo.Cat.V2O3 ./ MWinfo.MW.V2O3;
n.Cr2O3 = liquidOxides.Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
n.NiO = liquidOxides.NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;
n.P2O5 = liquidOxides.P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;
n.SO3 = liquidOxides.SO3 .* MWinfo.Cat.SO3 ./ MWinfo.MW.SO3;
n.Fe2O3 = liquidOxides.Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

cationTotal_liq_scalar = ...
    n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + n.MgO + ...
    n.CaO + n.Na2O + n.K2O + n.V2O3 + n.Cr2O3 + n.NiO + ...
    n.P2O5 + n.SO3 + n.Fe2O3;

% Direct division is intentional. A zero or NaN cation total produces
% non-finite cation fractions that are retained and reported.
XSiO2_liq_scalar = n.SiO2 ./ cationTotal_liq_scalar;
XTiO2_liq_scalar = n.TiO2 ./ cationTotal_liq_scalar;
XAlO1_5_liq_scalar = n.Al2O3 ./ cationTotal_liq_scalar;
XFeO_liq_scalar = n.FeO ./ cationTotal_liq_scalar;
XMnO_liq_scalar = n.MnO ./ cationTotal_liq_scalar;
XMgO_liq_scalar = n.MgO ./ cationTotal_liq_scalar;
XCaO_liq_scalar = n.CaO ./ cationTotal_liq_scalar;
XNaO0_5_liq_scalar = n.Na2O ./ cationTotal_liq_scalar;
XKO0_5_liq_scalar = n.K2O ./ cationTotal_liq_scalar;

XJd_cpx_scalar = cpx.XJd;
NaAl_liq_scalar = XNaO0_5_liq_scalar .* XAlO1_5_liq_scalar;
K_Jd_liq_scalar = XJd_cpx_scalar ./ ...
    ((XSiO2_liq_scalar.^2) .* XNaO0_5_liq_scalar .* XAlO1_5_liq_scalar);

% Natural logarithms are evaluated directly. Zero can therefore produce
% +/-Inf, and NaN remains NaN; both are retained for diagnostics.
lnK_Jd_liq_scalar = log(K_Jd_liq_scalar);
lnInvNaAl_liq_scalar = log(1 ./ NaAl_liq_scalar);

% Expand all composition-dependent scalars to the temperature-vector length.
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
oxygenSum_cpx = repmat(cpx.oxygenSum, nT, 1);

XAlIV_cpx = repmat(cpx.XAlIV, nT, 1);
XAlVI_cpx = repmat(cpx.XAlVI, nT, 1);
XJd_cpx = repmat(cpx.XJd, nT, 1);
XCaTs_cpx = repmat(cpx.XCaTs, nT, 1);
XCaTi_cpx = repmat(cpx.XCaTi, nT, 1);
XDiHd_cpx = repmat(cpx.XDiHd, nT, 1);
XEnFs_cpx = repmat(cpx.XEnFs, nT, 1);

SiO2_liq = repmat(liquidOxides.SiO2, nT, 1);
TiO2_liq = repmat(liquidOxides.TiO2, nT, 1);
Al2O3_liq = repmat(liquidOxides.Al2O3, nT, 1);
FeO_liq = repmat(liquidOxides.FeO, nT, 1);
MnO_liq = repmat(liquidOxides.MnO, nT, 1);
MgO_liq = repmat(liquidOxides.MgO, nT, 1);
CaO_liq = repmat(liquidOxides.CaO, nT, 1);
Na2O_liq = repmat(liquidOxides.Na2O, nT, 1);
K2O_liq = repmat(liquidOxides.K2O, nT, 1);
V2O3_liq = repmat(liquidOxides.V2O3, nT, 1);
Cr2O3_liq = repmat(liquidOxides.Cr2O3, nT, 1);
NiO_liq = repmat(liquidOxides.NiO, nT, 1);
P2O5_liq = repmat(liquidOxides.P2O5, nT, 1);
SO3_liq = repmat(liquidOxides.SO3, nT, 1);
Fe2O3_liq = repmat(liquidOxides.Fe2O3, nT, 1);
F_liq = repmat(liquidOxides.F, nT, 1);
Cl_liq = repmat(liquidOxides.Cl, nT, 1);

cationTotal_liq = repmat(cationTotal_liq_scalar, nT, 1);
XSiO2_liq = repmat(XSiO2_liq_scalar, nT, 1);
XTiO2_liq = repmat(XTiO2_liq_scalar, nT, 1);
XAlO1_5_liq = repmat(XAlO1_5_liq_scalar, nT, 1);
XFeO_liq = repmat(XFeO_liq_scalar, nT, 1);
XMnO_liq = repmat(XMnO_liq_scalar, nT, 1);
XMgO_liq = repmat(XMgO_liq_scalar, nT, 1);
XCaO_liq = repmat(XCaO_liq_scalar, nT, 1);
XNaO0_5_liq = repmat(XNaO0_5_liq_scalar, nT, 1);
XKO0_5_liq = repmat(XKO0_5_liq_scalar, nT, 1);

K_Jd_liq = repmat(K_Jd_liq_scalar, nT, 1);
lnK_Jd_liq = repmat(lnK_Jd_liq_scalar, nT, 1);
NaAl_liq = repmat(NaAl_liq_scalar, nT, 1);
lnInvNaAl_liq = repmat(lnInvNaAl_liq_scalar, nT, 1);

% ---- Putirka et al. (1996) pressure models
P1_kbar = -54.3 ...
    + 299 .* T_K ./ 1e4 ...
    + 36.4 .* T_K ./ 1e4 .* lnK_Jd_liq ...
    + 367 .* NaAl_liq;

P2_kbar = -50.7 ...
    + 394 .* T_K ./ 1e4 ...
    + 36.4 .* T_K ./ 1e4 .* lnK_Jd_liq ...
    - 20.0 .* T_K ./ 1e4 .* lnInvNaAl_liq;

P1_GPa = P1_kbar ./ 10;
P1_MPa = P1_kbar .* 100;
P2_GPa = P2_kbar ./ 10;
P2_MPa = P2_kbar .* 100;

isWithinCalibrationTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 1100 & T_degreeC <= 1475;
isWithinCalibrationPRange_P1 = ...
    isfinite(P1_kbar) & P1_kbar >= 8 & P1_kbar <= 30;
isWithinCalibrationPRange_P2 = ...
    isfinite(P2_kbar) & P2_kbar >= 8 & P2_kbar <= 30;
isWithinEquationDomain = ...
    isfinite(K_Jd_liq) & K_Jd_liq > 0 & ...
    isfinite(NaAl_liq) & NaAl_liq > 0;

% Pack equal-length vectors into the output table.
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
row.oxygenSum_cpx = oxygenSum_cpx;

row.XAlIV_cpx = XAlIV_cpx;
row.XAlVI_cpx = XAlVI_cpx;
row.XJd_cpx = XJd_cpx;
row.XCaTs_cpx = XCaTs_cpx;
row.XCaTi_cpx = XCaTi_cpx;
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

row.K_Jd_liq = K_Jd_liq;
row.lnK_Jd_liq = lnK_Jd_liq;
row.NaAl_liq = NaAl_liq;
row.lnInvNaAl_liq = lnInvNaAl_liq;

row.P1_kbar = P1_kbar;
row.P1_GPa = P1_GPa;
row.P1_MPa = P1_MPa;
row.P1_SEE_kbar = repmat(1.36, nT, 1);
row.P2_kbar = P2_kbar;
row.P2_GPa = P2_GPa;
row.P2_MPa = P2_MPa;
row.P2_SEE_kbar = repmat(1.51, nT, 1);

row.P_preferred_kbar = P1_kbar;
row.P_preferred_GPa = P1_GPa;
row.P_preferred_MPa = P1_MPa;

% Standard aliases used by the barometer launchers and plotting routines.
row.P_kbar = P1_kbar;
row.P_GPa = P1_GPa;
row.P_MPa = P1_MPa;
row.preferredel = repmat("P1", nT, 1);

row.isWithinEquationDomain = isWithinEquationDomain;
row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinCalibrationPRange_P1 = isWithinCalibrationPRange_P1;
row.isWithinCalibrationPRange_P2 = isWithinCalibrationPRange_P2;

end

function cpx = calculateCpxComponents(cpxOxides, MWinfo)
% calculateCpxComponents
% Normalize Clinopyroxene cations to 6 oxygens and calculate normative
% components. NaN values are propagated and never replaced by zero.

molProp = struct();
molProp.SiO2 = cpxOxides.SiO2 ./ MWinfo.MW.SiO2;
molProp.TiO2 = cpxOxides.TiO2 ./ MWinfo.MW.TiO2;
molProp.Al2O3 = cpxOxides.Al2O3 ./ MWinfo.MW.Al2O3;
molProp.FeO = cpxOxides.FeO ./ MWinfo.MW.FeO;
molProp.MnO = cpxOxides.MnO ./ MWinfo.MW.MnO;
molProp.MgO = cpxOxides.MgO ./ MWinfo.MW.MgO;
molProp.CaO = cpxOxides.CaO ./ MWinfo.MW.CaO;
molProp.Na2O = cpxOxides.Na2O ./ MWinfo.MW.Na2O;
molProp.K2O = cpxOxides.K2O ./ MWinfo.MW.K2O;
molProp.Cr2O3 = cpxOxides.Cr2O3 ./ MWinfo.MW.Cr2O3;

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

XAlIV = preserveNaNMaximum(2 - XSi, 0);
XAlVI = preserveNaNMaximum(XAl - XAlIV, 0);
XJd = preserveNaNMinimum(XNa, XAlVI);
XJd = preserveNaNMaximum(XJd, 0);
XCaTs = preserveNaNMaximum(XAlVI - XJd, 0);

if isnan(XAlIV) || isnan(XCaTs)
    XCaTi = NaN;
elseif XAlIV > XCaTs
    XCaTi = preserveNaNMaximum((XAlIV - XCaTs) ./ 2, 0);
else
    XCaTi = 0;
end

XDiHd = preserveNaNMaximum(XCa - XCaTs - XCaTi, 0);
XEnFs = preserveNaNMaximum((XFe + XMg - XDiHd) ./ 2, 0);

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
cpx.oxygenSum = oxygenSum;
cpx.XAlIV = XAlIV;
cpx.XAlVI = XAlVI;
cpx.XJd = XJd;
cpx.XCaTs = XCaTs;
cpx.XCaTi = XCaTi;
cpx.XDiHd = XDiHd;
cpx.XEnFs = XEnFs;

end

function outputValue = preserveNaNMaximum(inputValue, lowerBound)
% preserveNaNMaximum
% Apply a lower bound without allowing NaN to be interpreted as zero.

if isnan(inputValue)
    outputValue = NaN;
else
    outputValue = max(inputValue, lowerBound);
end

end

function outputValue = preserveNaNMinimum(valueA, valueB)
% preserveNaNMinimum
% Return the smaller scalar while preserving NaN explicitly.

if isnan(valueA) || isnan(valueB)
    outputValue = NaN;
else
    outputValue = min(valueA, valueB);
end

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Attach available Liquid identifiers and repeat them to match table height.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    value = data_liq.('Index');
    row.liq_Index = repeatTableValue(value, nRows);
end
if any(strcmp(variableNames, 'Experiment'))
    value = string(data_liq.('Experiment'));
    row.liq_Experiment = repmat(value, nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    value = string(data_liq.('Citation'));
    row.liq_Citation = repmat(value, nRows, 1);
end

end

function repeatedValue = repeatTableValue(value, nRows)
% repeatTableValue
% Repeat a scalar table value without changing its fundamental type.

if ~isscalar(value)
    error('Liquid identifier value must be scalar in a 1-row table.');
end

if iscell(value)
    repeatedValue = repmat(value, nRows, 1);
elseif isstring(value)
    repeatedValue = repmat(value, nRows, 1);
elseif ischar(value)
    repeatedValue = repmat(string(value), nRows, 1);
else
    repeatedValue = repmat(value, nRows, 1);
end

end

function displayPressureSummary(selectedCode_cpx, selectedLiquidRow, row)
% displayPressureSummary
% Print scalar or vector pressure summaries without assuming finite values.

label = [char(string(selectedCode_cpx)) ' & Liquid row ' ...
    num2str(selectedLiquidRow)];

if height(row) == 1
    disp([label ': P1 = ' num2str(row.P1_kbar) ' kbar; P2 = ' ...
        num2str(row.P2_kbar) ' kbar']);
else
    disp([label ': P1 = ' num2str(row.P1_kbar(1)) ' to ' ...
        num2str(row.P1_kbar(end)) ' kbar']);
    disp([label ': P2 = ' num2str(row.P2_kbar(1)) ' to ' ...
        num2str(row.P2_kbar(end)) ' kbar']);
end

end

function printPressureRangeWarning(pressure_kbar, modelName, ...
    minimumPressure, maximumPressure, selectedCode_cpx, selectedLiquidRow)
% printPressureRangeWarning
% Print a non-stopping warning for finite pressure values outside the main
% calibration range.

finitePressure = isfinite(pressure_kbar);
outsideRange = finitePressure & ...
    (pressure_kbar < minimumPressure | pressure_kbar > maximumPressure);

if any(outsideRange)
    finiteValues = pressure_kbar(finitePressure);
    fprintf(2, ...
        ['WARNING: %s calculated pressure is outside the main 8-30 kbar ' ...
         'range of the new Putirka et al. (1996) experiments (p. 92; ' ...
         'pp. 95-96). %d of %d finite pressure point(s) are outside the ' ...
         'range; calculated finite range = %.4g-%.4g kbar for ' ...
         'Clinopyroxene %s and Liquid row %d. Values are retained as ' ...
         'extrapolations.\n'], ...
        modelName, ...
        sum(outsideRange), ...
        sum(finitePressure), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(string(selectedCode_cpx)), ...
        selectedLiquidRow);
end

end

function printNonFinitePressureWarning(pressure_kbar, modelName, ...
    selectedCode_cpx, selectedLiquidRow)
% printNonFinitePressureWarning
% Retain and report NaN and Inf pressure values.

nonFinitePressure = ~isfinite(pressure_kbar);
if any(nonFinitePressure)
    fprintf(2, ...
        ['WARNING: Non-finite %s pressure values were calculated for ' ...
         'Clinopyroxene %s and Liquid row %d (%d of %d points; NaN: %d, ' ...
         'Inf: %d). These values remain in the output table and the ' ...
         'calculation has not been stopped.\n'], ...
        modelName, ...
        char(string(selectedCode_cpx)), ...
        selectedLiquidRow, ...
        sum(nonFinitePressure), ...
        numel(pressure_kbar), ...
        sum(isnan(pressure_kbar)), ...
        sum(isinf(pressure_kbar)));
end

end

function printNegativePressureWarning(pressure_kbar, modelName, ...
    selectedCode_cpx, selectedLiquidRow)
% printNegativePressureWarning
% Retain and report negative finite pressure values.

negativePressure = isfinite(pressure_kbar) & pressure_kbar < 0;
if any(negativePressure)
    fprintf(2, ...
        ['WARNING: Negative finite %s pressure was calculated for ' ...
         'Clinopyroxene %s and Liquid row %d (%d of %d points). The values ' ...
         'were retained for diagnostic purposes.\n'], ...
        modelName, ...
        char(string(selectedCode_cpx)), ...
        selectedLiquidRow, ...
        sum(negativePressure), ...
        numel(pressure_kbar));
end

end

function value = getMineralOxideRequired(dataTable, oxideName)
% getMineralOxideRequired
% Retrieve a required scalar oxide without altering NaN.

columnName = findOxideColumn(dataTable.Properties.VariableNames, oxideName);
if isempty(columnName)
    error('Selected Clinopyroxene row must contain variable: %s', oxideName);
end

value = toScalarDouble(dataTable.(columnName));

end

function value = getMineralOxideOptional(dataTable, oxideName, defaultValue)
% getMineralOxideOptional
% Return defaultValue only when the column is absent. NaN in an existing
% column is preserved.

columnName = findOxideColumn(dataTable.Properties.VariableNames, oxideName);
if isempty(columnName)
    value = defaultValue;
else
    value = toScalarDouble(dataTable.(columnName));
end

end

function value = getLiquidOxideOptional(dataTable, oxideName, defaultValue)
% getLiquidOxideOptional
% Return defaultValue only when the column is absent. NaN in an existing
% column is preserved.

columnName = findOxideColumn(dataTable.Properties.VariableNames, oxideName);
if isempty(columnName)
    value = defaultValue;
else
    value = toScalarDouble(dataTable.(columnName));
end

end

function columnName = findOxideColumn(variableNames, oxideName)
% findOxideColumn
% Locate oxide or oxide_value columns after canonicalizing variable names.

canonicalNames = cell(size(variableNames));
for i = 1:numel(variableNames)
    nameText = lower(variableNames{i});
    nameText = strrep(nameText, ' ', '');
    nameText = strrep(nameText, '_', '');
    nameText = strrep(nameText, '-', '');
    canonicalNames{i} = nameText;
end

canonicalOxide = lower(oxideName);
canonicalOxide = strrep(canonicalOxide, ' ', '');
canonicalOxide = strrep(canonicalOxide, '_', '');
canonicalOxide = strrep(canonicalOxide, '-', '');

targetNames = {[canonicalOxide 'value'], canonicalOxide};
columnName = '';

for i = 1:numel(targetNames)
    matchingIndex = find(strcmp(canonicalNames, targetNames{i}), 1, 'first');
    if ~isempty(matchingIndex)
        columnName = variableNames{matchingIndex};
        return
    end
end

end

function value = toScalarDouble(rawValue)
% toScalarDouble
% Convert one table value to a scalar double. Missing/blank/non-numeric text
% becomes NaN and is never converted to zero.

if isempty(rawValue)
    value = NaN;
    return
end

if isnumeric(rawValue) || islogical(rawValue)
    if ~isscalar(rawValue)
        error('Selected oxide value must be scalar in a 1-row table.');
    end
    value = double(rawValue);
    return
end

if isstring(rawValue)
    if ~isscalar(rawValue)
        error('Selected oxide value must be scalar in a 1-row table.');
    end
    if ismissing(rawValue)
        value = NaN;
    else
        value = str2double(rawValue);
    end
    return
end

if ischar(rawValue)
    value = str2double(string(rawValue));
    return
end

if iscell(rawValue)
    if numel(rawValue) ~= 1
        error('Selected oxide value must be scalar in a 1-row table.');
    end
    cellValue = rawValue{1};
    if isempty(cellValue)
        value = NaN;
    elseif isnumeric(cellValue) || islogical(cellValue)
        if ~isscalar(cellValue)
            error('Selected oxide cell value must be scalar.');
        end
        value = double(cellValue);
    elseif isstring(cellValue) || ischar(cellValue)
        value = str2double(string(cellValue));
    else
        value = NaN;
    end
    return
end

value = NaN;

end
