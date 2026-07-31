function results = LesselPutirka2015CpxLiqbaro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Pyroxene/LesselPutirka2015CpxLiqbaro.m
% Compatibility target: MATLAB R2024b
%
% Clinopyroxene-Liquid barometer
% Lessel, J. and Putirka, K. (2015)
% American Mineralogist, 100, 2163-2171
% DOI: https://doi.org/10.2138/am-2015-4732
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Clinopyroxene analysis selected
% from rawdata_struct.Cpx with one Liquid analysis read by
% liquid.readLiquidExcel and calculates pressure using Equation (1) of
% Lessel and Putirka (2015; p. 2165).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Clinopyroxene-Liquid pair, one
% output row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% IMPORTANT: This barometer was calibrated specifically for MARTIAN
% igneous compositions and martian meteorite analog experiments. Martian
% liquids are characteristically richer in FeO and poorer in Al2O3 than
% common terrestrial basaltic liquids (pp. 2163-2164). Application to
% terrestrial magmatic systems, or to compositions unlike the experimental
% martian database, is compositional extrapolation.
%
% The combined experimental database used in the paper spans broadly
% (pp. 2164-2165; Table 1):
%
%   Temperature       : 950-1540 degreeC
%   Pressure          : approximately 1 atm-2.3 GPa
%                       (approximately 0.001-23 kbar)
%   Liquid SiO2       : 40.2-66.16 wt%
%   Liquid MgO        : 0.62-24.52 wt%
%   Liquid FeO        : 2.80-30.2 wt%
%   Liquid Al2O3      : 2.97-20.5 wt%
%   Liquid Na2O + K2O : 0.19-6.77 wt%
%
% These are ranges for the COMBINED experimental database, not strict
% equation-specific calibration limits for Equation (1). They are used in
% this implementation only as non-stopping warning envelopes.
%
% Equation (1) was calibrated using clinopyroxene-saturated partial-melting
% experiments on martian bulk compositions (p. 2165). It reproduces the
% calibration pressures with R2 = 0.91 and RMSE = +/-0.17 GPa (n = 42),
% and predicts the independent test data with R2 = 0.94 and RMSE =
% +/-0.22 GPa (n = 9; pp. 2164-2165, Table 2 and Fig. 2).
%
% Clinopyroxene and Liquid must represent an equilibrium pair. Lessel and
% Putirka (2015) recommend checking Fe-Mg exchange using their Equation
% (13) on p. 2168:
%
%   KD(Fe-Mg)Cpx-Liq = 0.32 - 0.02*(Na2O + K2O)_liq(wt%)
%
% with SEE approximately +/-0.03 after exclusion of greater-than-3-sigma
% outliers. Their KD calculations use a liquid Fe3+/Fe2+ correction based
% on Kress and Carmichael (1991, Eq. 7). This function does not perform that
% redox-dependent equilibrium screening automatically; users must verify
% equilibrium independently before interpreting the calculated pressure.
%
% Volatile-bearing experiments were included in the overall calibration and
% test database and were treated in the same manner as other experiments
% (p. 2165). However, no separate H2O, F, or Cl correction is included in
% Equation (1). In this implementation, Liquid F and Cl are excluded from
% cationTotal_liq and are excluded from the NaN-input warning check.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 950-1540 degreeC,
%   2) finite calculated pressure is outside approximately 0.001-23 kbar,
%   3) finite liquid composition is outside the broad combined-database
%      ranges listed above,
%   4) a calculation input contains NaN, or
%   5) a calculated pressure is NaN or Inf.
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if it makes a logarithm or ratio
% undefined, the resulting NaN or Inf is retained and reported.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST column of the Cpx table is treated as an identifier ("data
% code") displayed in the selection dialog. Clinopyroxene oxide columns are
% read using common names such as SiO2, TiO2, Al2O3, FeO or FeOt, MnO, MgO,
% CaO, Na2O, K2O, and Cr2O3. Columns with a "_value" suffix are also
% recognized.
%
% Liquid compositions are obtained using liquid.readLiquidExcel. The
% Equation (1) calculation uses liquid cation fractions normalized to the
% sum of cations contributed by:
%
%   SiO2, TiO2, Al2O3, FeO, MnO, MgO, CaO, Na2O, K2O,
%   V2O3, Cr2O3, NiO, P2O5, SO3, and Fe2O3.
%
% F and Cl are intentionally excluded from cationTotal_liq. Existing NaN
% values in calculation-input columns remain NaN and are never replaced by
% zero. An optional oxide column that is entirely absent from the input
% table is treated as zero for backward compatibility; a present column
% containing NaN remains NaN and propagates through the calculation.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
% Lessel and Putirka (2015), Equation (1), p. 2165:
%
%   P(GPa) =
%     -412.7
%     + 9.667e-6*T(K)*ln[XJd_Cpx / ...
%       (XAlO1.5_Liq*XNaO0.5_Liq*(XSiO2_Liq)^2)]
%     + 722.1*XSiO2_Liq
%     + 5.496*XMgO_Liq
%     - 195.0*ln(XSiO2_Liq)
%     - 334.7*(XSiO2_Liq)^2
%     + 686.3*(XTiO2_Liq)^2
%     + 30.77*(XAlO1.5_Liq)^2
%     - 57.14*(XCaO_Liq)^2
%
% NOTE: The coefficient of XSiO2_Liq is 722.1 in the original paper, not
% 772.1. The value is corrected in this implementation.
%
% Clinopyroxene components are calculated on a six-oxygen basis using the
% normative scheme described by Putirka (2008), as specified by Lessel and
% Putirka (2015; pp. 2164-2165).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = LesselPutirka2015CpxLiqbaro(rawdata_struct, T_degreeC)
%   results = LesselPutirka2015CpxLiqbaro(..., 'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing a Cpx table
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Name-value option:
%   'LiquidRow'    : positive integer scalar selecting a row in the loaded
%                    Liquid table. Default [] uses row 1.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Clinopyroxene-Liquid pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('LesselPutirka2015CpxLiqbaro requires (rawdata_struct, T_degreeC).');
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
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max(16, height(dataset_cpx));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Broad ranges of the combined experimental database. These are used only
% as warning envelopes and are not strict Equation (1)-specific limits.
experimentalT_min_degreeC = 950;
experimentalT_max_degreeC = 1540;
experimentalP_min_kbar = 0.00101325;
experimentalP_max_kbar = 23;

temperatureOutsideExperimental = isfinite(T_degreeC) & ...
    (T_degreeC < experimentalT_min_degreeC | ...
     T_degreeC > experimentalT_max_degreeC);
temperatureWarningIssued = false;
modelCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Clinopyroxene) ===');

while true
    % ----- Clinopyroxene selection -----
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

    % ----- Calculation -----
    disp('=== Step 4: Checking calculation inputs ===');

    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    % Check NaN only in variables that directly enter cation normalization,
    % component calculation, temperature, or Equation (1). Liquid F and Cl
    % are intentionally excluded from this check.
    nanInputNames = findNaNInputs( ...
        selectedData_cpx, selectedData_liq, T_degreeC);

    % Reject Inf and finite negative values in calculation inputs. NaN and
    % zero are retained. Liquid F and Cl are not calculation inputs here.
    validateNonNegativeInputs(selectedData_cpx, selectedData_liq);

    disp('=== Step 5: Calculating the pressure ===');
    row = calcPressure( ...
        selectedData_cpx, selectedData_liq, T_degreeC, MWinfo);

    % Repeat identifiers for all temperatures in the current calculation.
    nRows = height(row);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_cpx', 'dataRow_liq'}, 'Before', 1);

    % Store one block per selected pair. Expand only when the preallocated
    % cell-buffer capacity is exhausted, rather than on every loop pass.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        expandedBuffer = cell(2 * numel(resultBlocks), 1);
        expandedBuffer(1:numel(resultBlocks)) = resultBlocks;
        resultBlocks = expandedBuffer;
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
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

    % Print the model-specific limitation once per function call.
    if ~modelCautionIssued
        fprintf(2, ...
            ['CAUTION: Lessel and Putirka (2015) Equation (1) was calibrated ' ...
             'specifically for martian igneous compositions (pp. 2163-2165). ' ...
             'Application to terrestrial or compositionally dissimilar systems is ' ...
             'extrapolation. Clinopyroxene-Liquid equilibrium must be checked ' ...
             'independently (Eq. 13, p. 2168).\n']);
        modelCautionIssued = true;
    end

    % Input temperature is common to all selected pairs, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideExperimental) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the 950-1540 degreeC range ' ...
             'of the combined experimental database of Lessel and Putirka ' ...
             '(2015; pp. 2164-2165). This is not a strict Equation (1)-specific ' ...
             'range. %d of %d finite temperature point(s) are outside the range; ' ...
             'input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideExperimental), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures fall outside the broad pressure
    % envelope represented by the combined experiments.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideExperimental = finitePressure & ...
        (row.P_kbar < experimentalP_min_kbar | ...
         row.P_kbar > experimentalP_max_kbar);

    if any(pressureOutsideExperimental)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the approximately ' ...
             '1-atm-2.3-GPa (0.001-23 kbar) range of the combined experimental ' ...
             'database of Lessel and Putirka (2015; pp. 2164-2165). This is not ' ...
             'a strict Equation (1)-specific range. %d of %d finite pressure ' ...
             'point(s) are outside the range; calculated finite range = ' ...
             '%.4g-%.4g kbar for %s & Liquid row %d.\n'], ...
            sum(pressureOutsideExperimental), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq);
    end

    % Warn when finite liquid compositions lie outside the broad combined
    % experimental database ranges listed in the paper.
    printLiquidCompositionWarnings( ...
        row, selectedCode_cpx, selectedIdx_liq);

    % List the exact calculation-input names containing NaN. For vector
    % temperature input, the NaN element indices are included. F and Cl are
    % excluded from this warning by design.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s & ' ...
             'Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'pressure may remain NaN. Liquid F and Cl were excluded from this ' ...
             'check and from cationTotal_liq.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s & ' ...
             'Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressure is retained for diagnostics but is outside
    % the physical/useful domain and also triggers the range warning above.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & ' ...
             'Liquid row %d (%d of %d points). The values were retained for ' ...
             'diagnostic purposes.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    % Report a non-positive or non-finite logarithm argument explicitly.
    invalidLogArgument = ~isfinite(row.logArg_Eq1) | row.logArg_Eq1 <= 0;
    if any(invalidLogArgument)
        fprintf(2, ...
            ['WARNING: Equation (1) logarithm argument is non-positive or ' ...
             'non-finite for %s & Liquid row %d (%d of %d points). The resulting ' ...
             'NaN or Inf values were retained.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            sum(invalidLogArgument), ...
            numel(row.logArg_Eq1));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Clinopyroxene selection (same Liquid row)?', ...
        'LesselPutirka2015CpxLiqbaro', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered result blocks once. If the user canceled before
% any calculation, return an empty table.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct('liquid', metaLiq);
disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_cpx, data_liq, T_degreeC)
% findNaNInputs
% Return names of Equation (1) calculation inputs containing NaN. NaN values
% do not cause an error and are not replaced by zero. Liquid F and Cl are
% intentionally excluded.

maxNames = 32;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

% Cpx oxides used for six-oxygen normalization and component calculation.
cpxOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO','Na2O','K2O','Cr2O3'};
for i = 1:numel(cpxOxides)
    oxide = cpxOxides{i};
    columnName = findOxideColumn(data_cpx.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDouble(data_cpx.(columnName));
        if isnan(value)
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = ...
                "Cpx." + string(columnName);
        end
    end
end

% Use FeO when that column exists; otherwise use FeOt.
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

% Liquid oxides included in cationTotal_liq. F and Cl are excluded.
liqOxides = {'SiO2','TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O', ...
    'V2O3','Cr2O3','NiO','P2O5','SO3','Fe2O3'};
for i = 1:numel(liqOxides)
    oxide = liqOxides{i};
    columnName = findOxideColumn(data_liq.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDouble(data_liq.(columnName));
        if isnan(value)
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = ...
                "Liquid." + string(columnName);
        end
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_cpx, data_liq)
% validateNonNegativeInputs
% Reject Inf and finite negative values in raw oxide variables used by the
% calculation. Zero and NaN are intentionally allowed and retained. Liquid
% F and Cl are not included because they are excluded from cationTotal_liq.

maxNames = 32;
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

cpxOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO','Na2O','K2O','Cr2O3'};
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

liqOxides = {'SiO2','TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O', ...
    'V2O3','Cr2O3','NiO','P2O5','SO3','Fe2O3'};
for i = 1:numel(liqOxides)
    oxide = liqOxides{i};
    columnName = findOxideColumn(data_liq.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDouble(data_liq.(columnName));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = ...
                "Liquid." + string(columnName);
        end
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['LesselPutirka2015CpxLiqbaro: calculation inputs must be ' ...
           'non-negative. NaN is allowed, but Inf and finite negative ' ...
           'value(s) are prohibited. Invalid input(s): ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_cpx, data_liq, T_degreeC, MWinfo)
% calcPressure
% Compute pressure for one Clinopyroxene row and one Liquid row at one or
% more input temperatures. NaN values propagate through the calculation.
%
% Inputs:
%   data_cpx    : 1-row Clinopyroxene table
%   data_liq    : 1-row Liquid table
%   T_degreeC   : scalar or vector temperature in degreeC
%   MWinfo      : molar-weight and cation-number structure
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Extract and normalize one-row Cpx and Liquid data. Existing NaN values are
% preserved. Missing optional columns use zero only for backward
% compatibility.
cpx = prepareCpxRow(data_cpx, MWinfo);
liq = extractLiquidRow(data_liq, MWinfo);

% Equation (1) logarithm argument. XJd and all liquid fractions are scalars
% for a selected pair. Zero values are retained: log(0) becomes -Inf, while
% 0/0 becomes NaN, and both are reported without stopping the calculation.
logArgScalar = cpx.XJd ./ ...
    (liq.XAlO1_5 .* liq.XNaO0_5 .* (liq.XSiO2 .^ 2));
logTermScalar = log(logArgScalar);

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

logArg_Eq1 = repmat(logArgScalar, nT, 1);
logTerm_Eq1 = repmat(logTermScalar, nT, 1);

% Lessel and Putirka (2015) Equation (1), p. 2165. No finite-value guard is
% used: NaN and Inf propagate and remain in the output table.
P_GPa = ...
    -412.7 ...
    + 9.667e-6 .* T_K .* logTerm_Eq1 ...
    + 722.1 .* XSiO2_liq ...
    + 5.496 .* XMgO_liq ...
    - 195.0 .* log(XSiO2_liq) ...
    - 334.7 .* (XSiO2_liq .^ 2) ...
    + 686.3 .* (XTiO2_liq .^ 2) ...
    + 30.77 .* (XAlO1_5_liq .^ 2) ...
    - 57.14 .* (XCaO_liq .^ 2);

P_kbar = P_GPa .* 10;

% Applicability and diagnostic flags. Combined-database ranges are warning
% envelopes only, not strict Equation (1)-specific calibration limits.
isWithinEquationDomain = ...
    isfinite(logArg_Eq1) & logArg_Eq1 > 0 & ...
    isfinite(XSiO2_liq) & XSiO2_liq > 0;

isWithinExperimentalTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 950 & T_degreeC <= 1540;

isWithinExperimentalPRange = ...
    isfinite(P_kbar) & P_kbar >= 0.00101325 & P_kbar <= 23;

finiteBroadComposition = ...
    isfinite(SiO2_liq) & isfinite(MgO_liq) & ...
    isfinite(FeO_liq) & isfinite(Al2O3_liq) & ...
    isfinite(Na2O_liq) & isfinite(K2O_liq);

isWithinBroadLiquidRange = finiteBroadComposition & ...
    SiO2_liq >= 40.2 & SiO2_liq <= 66.16 & ...
    MgO_liq >= 0.62 & MgO_liq <= 24.52 & ...
    FeO_liq >= 2.80 & FeO_liq <= 30.2 & ...
    Al2O3_liq >= 2.97 & Al2O3_liq <= 20.5 & ...
    (Na2O_liq + K2O_liq) >= 0.19 & ...
    (Na2O_liq + K2O_liq) <= 6.77;

% Pack outputs using pre-sized vectors of equal height.
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

row.logArg_Eq1 = logArg_Eq1;
row.logTerm_Eq1 = logTerm_Eq1;

% Primary pressure variables used by the common barometer launchers.
row.P_GPa = P_GPa;
row.P_kbar = P_kbar;

% Backward-compatible aliases retained from the original implementation.
row.PEq1_GPa = P_GPa;
row.PEq1_kbar = P_kbar;

row.P_RMSE_calibration_GPa = repmat(0.17, nT, 1);
row.P_RMSE_test_GPa = repmat(0.22, nT, 1);
row.P_RMSE_calibration_kbar = repmat(1.7, nT, 1);
row.P_RMSE_test_kbar = repmat(2.2, nT, 1);

row.isWithinEquationDomain = isWithinEquationDomain;
row.isWithinExperimentalTRange = isWithinExperimentalTRange;
row.isWithinExperimentalPRange = isWithinExperimentalPRange;
row.isWithinBroadLiquidRange = isWithinBroadLiquidRange;

% Backward-compatible/general aliases.
row.isRecommended_T_range = isWithinExperimentalTRange;
row.isApplicable_composition = isWithinBroadLiquidRange;

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
% Extract one-row Liquid oxide data and calculate cation fractions. F and Cl
% are retained only as output diagnostics and are excluded from
% cationTotal_liq. Existing NaN values in included oxides are retained.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liq = struct();
liq.SiO2 = getLiquidOxOptional(data_liq, 'SiO2', 0);
liq.TiO2 = getLiquidOxOptional(data_liq, 'TiO2', 0);
liq.Al2O3 = getLiquidOxOptional(data_liq, 'Al2O3', 0);
liq.FeO = getLiquidOxOptional(data_liq, 'FeO', 0);
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

% F and Cl do not enter cationTotal_liq and do not trigger NaN warnings.
% Missing values are retained as NaN for diagnostic output.
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

% F and Cl are intentionally absent from this sum.
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

function printLiquidCompositionWarnings(row, selectedCode_cpx, selectedIdx_liq)
% printLiquidCompositionWarnings
% Print non-stopping warnings when finite Liquid values are outside the broad
% combined experimental-database ranges reported on p. 2165.

SiO2 = row.SiO2_liq(1);
MgO = row.MgO_liq(1);
FeO = row.FeO_liq(1);
Al2O3 = row.Al2O3_liq(1);
Na2O = row.Na2O_liq(1);
K2O = row.K2O_liq(1);

maxMessages = 5;
messageBuffer = strings(maxMessages, 1);
nMessages = 0;

if isfinite(SiO2) && (SiO2 < 40.2 || SiO2 > 66.16)
    nMessages = nMessages + 1;
    messageBuffer(nMessages) = sprintf('SiO2 = %.4g wt%% (range 40.2-66.16)', SiO2);
end
if isfinite(MgO) && (MgO < 0.62 || MgO > 24.52)
    nMessages = nMessages + 1;
    messageBuffer(nMessages) = sprintf('MgO = %.4g wt%% (range 0.62-24.52)', MgO);
end
if isfinite(FeO) && (FeO < 2.80 || FeO > 30.2)
    nMessages = nMessages + 1;
    messageBuffer(nMessages) = sprintf('FeO = %.4g wt%% (range 2.80-30.2)', FeO);
end
if isfinite(Al2O3) && (Al2O3 < 2.97 || Al2O3 > 20.5)
    nMessages = nMessages + 1;
    messageBuffer(nMessages) = sprintf('Al2O3 = %.4g wt%% (range 2.97-20.5)', Al2O3);
end

alkalis = Na2O + K2O;
if isfinite(alkalis) && (alkalis < 0.19 || alkalis > 6.77)
    nMessages = nMessages + 1;
    messageBuffer(nMessages) = sprintf( ...
        'Na2O + K2O = %.4g wt%% (range 0.19-6.77)', alkalis);
end

if nMessages > 0
    messages = messageBuffer(1:nMessages);
    fprintf(2, ...
        ['WARNING: Liquid composition is outside the broad combined ' ...
         'experimental-database range of Lessel and Putirka (2015; p. 2165) ' ...
         'for %s & Liquid row %d: %s. These are not strict Equation ' ...
         '(1)-specific limits.\n'], ...
        char(string(selectedCode_cpx)), ...
        selectedIdx_liq, ...
        char(strjoin(messages, '; ')));
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
