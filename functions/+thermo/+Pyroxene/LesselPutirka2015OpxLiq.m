function results = LesselPutirka2015OpxLiq(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Pyroxene/LesselPutirka2015OpxLiq.m
% Tested with MATLAB R2024b
%
% Orthopyroxene-Liquid thermometer for martian igneous compositions
% Lessel, J. and Putirka, K. (2015), equation (4)
% American Mineralogist, 100, 2163-2171
% DOI: https://doi.org/10.2138/am-2015-4732
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Orthopyroxene analysis and pairs
% it with one row from a separately loaded Liquid dataset. Temperature is
% calculated using equation (4) of Lessel and Putirka (2015).
%
% The function accepts either a scalar pressure or a pressure vector and is
% therefore compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP. Equation (4) contains no pressure term, so the same
% calculated temperature is repeated once for every input pressure value.
% Pressure is retained in the output for interface compatibility and
% traceability only.
%
% F and Cl may be retained in the output as raw Liquid values, but they are
% anions and are excluded from the Liquid cation-fraction denominator. They
% are also excluded from NaN and validity checks because they do not enter
% the thermometer calculation.
%
% The function is designed for repeated calculations. Each selected-pair
% result is stored temporarily in a preallocated cell buffer, and all result
% blocks are concatenated only once after the interactive loop finishes.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Lessel and Putirka (2015) developed this thermometer specifically from
% experiments on martian meteorites and martian-analog bulk compositions.
% Martian basaltic liquids are generally richer in FeO and poorer in Al2O3
% than common terrestrial basaltic liquids. Application to ordinary
% terrestrial igneous compositions is therefore extrapolative.
%
% The complete experimental compilation used in the paper spans:
%
%   Pressure     : approximately 0.0001-2.3 GPa
%                  (approximately 0.001-23 kbar)
%   Temperature  : approximately 950-1540 degreeC
%   Liquid SiO2  : 40.2-66.16 wt.%
%   Liquid MgO   : 0.62-24.52 wt.%
%   Liquid FeO   : 2.80-30.2 wt.%
%   Liquid Al2O3 : 2.97-20.5 wt.%
%   Total alkalis: 0.19-6.77 wt.%
%
% These overall limits are reported in the Methodology section on
% pp. 2164-2165. They combine several mineral-liquid and mineral-mineral
% models and are not strict rectangular limits for equation (4) alone.
%
% From the Opx-bearing source experiments listed in Table 1, the approximate
% phase-specific source-study envelope is:
%
%   Pressure    : approximately 0.0001-2.3 GPa
%   Temperature : approximately 950-1540 degreeC
%
% This source-study envelope is used below for non-stopping warnings. It is
% not presented by the authors as a strict rectangular calibration boundary
% for the final equation (4) regression subset.
%
% Equation (4) was calibrated with 74 data and reproduced temperature with
% R2 = 0.94 and RMSE = 33 K. An independent test set of 14 data was predicted
% with R2 = 0.96 and RMSE = 27 K (Table 2 on p. 2164; equation (4) and
% Figure 5 on p. 2166). The test-set RMSE of approximately 27 K provides a
% useful precision estimate for compositions similar to the calibration set.
%
% Mineral and liquid components must be calculated in the same way as in
% the paper: liquid components are anhydrous cation fractions, and Opx
% cations are calculated on a 6-oxygen basis (pp. 2164-2166).
%
% The selected Opx and Liquid must represent an equilibrated pair. Lessel
% and Putirka (2015) report that a single mean Opx-Liquid Fe-Mg exchange
% coefficient should not be used because KD varies strongly with Liquid K2O.
% They recommend the composition-dependent equilibrium relation:
%
%   KD(Fe-Mg)Opx-Liq = 0.32 - 0.05*K2O_liq,wt%
%
% with a standard error of estimate of approximately 0.03 (p. 2167,
% equation (12)). This equilibrium test is not calculated automatically in
% this implementation because rigorous use requires a consistent treatment
% of Fe2+/Fe3+ in the Liquid.
%
% Whole-rock compositions should be used as liquids only when they are a
% defensible proxy for the melt that equilibrated with the selected Opx.
% Cumulate whole-rock compositions, evolved mesostasis, altered glass, and
% unrelated crystal-liquid pairs may yield geologically meaningless results.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside the approximate Opx-bearing source-study
%      envelope of 0.0001-2.3 GPa, or
%   2) a finite calculated temperature is outside the approximate
%      Opx-bearing source-study envelope of 950-1540 degreeC.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%
% The FIRST column of the Opx table is treated as an identifier ("data
% code") displayed in the selection dialog. Opx oxide columns may use either
% the oxide name itself or the oxide name followed by "Value". Spaces,
% underscores, and hyphens in column names are ignored when matching.
%
% Opx oxide variables used by the component calculation:
%   SiO2              % required column
%   TiO2              % absent column is treated as 0
%   Al2O3             % absent column is treated as 0
%   FeO or FeOt       % at least one column is required
%   MnO               % absent column is treated as 0
%   MgO               % required column
%   CaO               % absent column is treated as 0
%   Na2O              % absent column is treated as 0
%   K2O               % absent column is treated as 0
%   Cr2O3             % absent column is treated as 0
%
% If FeO is present, its selected value is used directly, including NaN. If
% FeO is absent, FeOt is used. A present NaN is never replaced by FeOt or 0.
%
% The Liquid dataset is loaded with liquid.readLiquidExcel(). The following
% oxide columns are included in the anhydrous cation-fraction normalization:
%
%   SiO2 TiO2 Al2O3 FeO MnO MgO CaO Na2O K2O
%   V2O3 Cr2O3 NiO P2O5 SO3 Fe2O3
%
% F and Cl may be retained as raw output values when present, but they are
% anions and are excluded from cationTotal_liq and from thermometer-input
% NaN/validity checks. Their values do not affect equation (4).
%
% A missing optional oxide column is treated as zero. A column that is
% present but contains NaN remains NaN and propagates through the
% calculation. Every selected calculation input containing NaN is reported
% by fprintf.
%
% All finite oxide values used by the calculation must be greater than or
% equal to zero. Finite negative values and Inf are rejected. Zero is allowed
% as an input, but a zero in a logarithm or denominator may produce a NaN
% temperature; the result is retained and reported without stopping.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Lessel and Putirka (2015), equation (4):
%
%   1 / T(K) =
%       2.234e-4
%     - 3.982e-6 * DMgO_OpxLiq
%     - 9.394e-4 * Na_Opx
%     - 2.058e-4 * ln(XMgO_liq)
%     + 3.564e-5 * ln(XCaO_liq)
%     + 2.186e-4 * ln(Mg_Opx)
%     - 1.958e-5 * ln(AlIV_Opx)
%     + 0.002298 * (AlIV_Opx^2)
%
% where:
%   DMgO_OpxLiq = Mg_Opx / XMgO_liq
%   AlIV_Opx     = 2 - Si_Opx
%
% Natural logarithms are used. The pressure input is not used by equation
% (4), but one output row is returned for each pressure value.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = LesselPutirka2015OpxLiq(rawdata_struct, P_kbar)
%   results = LesselPutirka2015OpxLiq(..., 'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing the Opx table described above
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Name-value option:
%   LiquidRow      : positive integer row index or [] (default []). When [],
%                    row 1 of the selected Liquid dataset is used.
%
% Output:
%   results : table containing one row per input pressure value for every
%             selected Opx-Liquid pair. NaN and Inf results are retained.
%

%% Input validation
if nargin < 2
    error('LesselPutirka2015OpxLiq requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end
if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end

P_kbar = P_kbar(:);

%% Options
ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == floor(x)));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve mineral and liquid datasets
disp('=== Step 1: Preparing mineral and liquid datasets ===');

dataset_opx = rawdata_struct.Opx;
MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();

if ~istable(liqAll) || isempty(liqAll)
    error('Selected Liquid dataset must be a non-empty table.');
end

% Resolve the liquid row once because the same liquid dataset is retained
% throughout the interactive Opx-selection loop.
if isempty(liquidRowOpt)
    idxLiq = 1;
    if height(liqAll) > 1
        fprintf(2, ...
            ['WARNING: The selected Liquid dataset contains %d rows. ' ...
             'LiquidRow was not specified, so row 1 will be used.\n'], ...
            height(liqAll));
    end
else
    idxLiq = liquidRowOpt;
    if idxLiq > height(liqAll)
        error(['Requested LiquidRow (%d) exceeds the number of rows in the ' ...
               'selected Liquid dataset (%d).'], idxLiq, height(liqAll));
    end
end

selectedData_liq = liqAll(idxLiq, :);
disp(['Liquid selected: Row ' num2str(idxLiq)]);
disp('=== Preparing mineral and liquid datasets has been finished ===');

%% 2) Initialize output container
% Each selected-pair result is stored as one table block. The buffer is
% preallocated and is enlarged only when its capacity is exhausted, rather
% than changing the complete results table on every loop iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max(16, height(dataset_opx));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Approximate Opx-bearing source-study envelope inferred from Table 1 of
% Lessel and Putirka (2015). These are warning limits, not strict boundaries
% explicitly defined for the final equation (4) regression subset.
applicationP_min_GPa = 0.0001;
applicationP_max_GPa = 2.3;
applicationT_min_degC = 950;
applicationT_max_degC = 1540;

P_GPa_input = P_kbar ./ 10;
pressureOutsideRange = P_GPa_input < applicationP_min_GPa | ...
    P_GPa_input > applicationP_max_GPa;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive selection loop and calculation
disp('=== Step 3: Selecting a data code from the list (Opx) ===');

while true
    % ----- Opx selection -----
    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Opx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Opx selected: ' char(string(selectedCode_opx))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    selectedData_opx = dataset_opx(selectedIdx_opx, :);

    % Identify NaN values in every present oxide column used by the mineral
    % and Liquid normalization. F and Cl are intentionally excluded.
    nanInputNames = findNaNInputs(selectedData_opx, selectedData_liq);

    % Reject finite negative inputs and Inf. Zero and NaN are deliberately
    % allowed so equation-domain failures remain non-stopping results.
    validateInputValues(selectedData_opx, selectedData_liq);

    row = calcTemp(selectedData_opx, selectedData_liq, P_kbar, MWinfo);

    % Store identifiers once per pressure row for traceability.
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row.dataRow_liq = repmat(idxLiq, height(row), 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_opx', 'dataRow_liq'}, 'Before', 1);

    % Store this result block without repeatedly concatenating the complete
    % output table during the interactive loop.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_opx)) ' & Liquid row ' ...
            num2str(idxLiq) ': ' num2str(row.TEq4_C) ' degreeC']);
    else
        disp([char(string(selectedCode_opx)) ' & Liquid row ' ...
            num2str(idxLiq) ': ' num2str(row.TEq4_C(1)) ' to ' ...
            num2str(row.TEq4_C(end)) ' degreeC']);
    end

    % Warn once when input pressure lies outside the approximate source-study
    % envelope. Pressure is not used in equation (4), but the warning records
    % whether the requested P context lies outside the source experiments.
    if any(pressureOutsideRange) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate Opx-bearing ' ...
             'source-study envelope of Lessel and Putirka (2015): ' ...
             '0.0001-2.3 GPa (0.001-23 kbar). %d of %d pressure point(s) ' ...
             'are outside the range; input range = %.6g-%.6g GPa. ' ...
             'Pressure is not used in equation (4).\n'], ...
            sum(pressureOutsideRange), numel(P_GPa_input), ...
            min(P_GPa_input), max(P_GPa_input));
        pressureWarningIssued = true;
    end

    % Warn when finite calculated temperatures lie outside the approximate
    % Opx-bearing source-study envelope.
    finiteTemperature = isfinite(row.TEq4_C);
    temperatureOutsideRange = finiteTemperature & ...
        (row.TEq4_C < applicationT_min_degC | ...
         row.TEq4_C > applicationT_max_degC);

    if any(temperatureOutsideRange)
        finiteValues = row.TEq4_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximate ' ...
             'Opx-bearing source-study envelope of Lessel and Putirka ' ...
             '(2015): 950-1540 degreeC. %d of %d finite temperature ' ...
             'point(s) are outside the range; calculated finite range = ' ...
             '%.6g-%.6g degreeC for Opx %s and Liquid row %d.\n'], ...
            sum(temperatureOutsideRange), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_opx)), idxLiq);
    end

    % Report all calculation inputs that contained NaN. Values remain NaN
    % and are never replaced by zero.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for Opx %s ' ...
             'and Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; they were ' ...
             'not replaced by zero.\n'], ...
            char(string(selectedCode_opx)), idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report non-finite temperature results.
    invalidTemperature = ~isfinite(row.TEq4_C);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for ' ...
             'Opx %s and Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_opx)), idxLiq, ...
            sum(invalidTemperature), numel(row.TEq4_C), ...
            sum(isnan(row.TEq4_C)), sum(isinf(row.TEq4_C)));
    end

    % Report equation-domain failures separately so NaN causes are visible.
    invalidEquationDomain = ~row.equation_domain_valid;
    if any(invalidEquationDomain)
        fprintf(2, ...
            ['WARNING: Lessel and Putirka (2015) equation (4) was outside ' ...
             'its mathematical domain for Opx %s and Liquid row %d. ' ...
             'XMgO_liq, XCaO_liq, Mg_Opx, AlIV_Opx, and DMgO_OpxLiq ' ...
             'must be finite and > 0; Na_Opx must be finite and >= 0; ' ...
             'and 1/T must be finite and > 0.\n'], ...
            char(string(selectedCode_opx)), idxLiq);
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another Opx analysis and the same liquid row.
    userAction = questdlg( ...
        'Continue with another data selection (same Liquid dataset)?', ...
        'LesselPutirka2015OpxLiq', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all result blocks once after all selections are complete.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct('liquid', metaLiq);
disp('=== LesselPutirka2015OpxLiq finished ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_opx, data_liq)
% findNaNInputs
% Return names of present Opx and Liquid oxide inputs whose selected values
% are NaN or cannot be converted to a numeric scalar. Missing optional
% columns are not listed because they follow the documented zero assumption.
% Liquid F and Cl are excluded because they do not enter the calculation.

opxOxides = {'SiO2', 'TiO2', 'Al2O3', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};
liqOxides = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', ...
    'SO3', 'Fe2O3'};

maxEntries = numel(opxOxides) + numel(liqOxides) + 1;
nanBuffer = strings(maxEntries, 1);
nNan = 0;

for i = 1:numel(opxOxides)
    oxide = opxOxides{i};
    columnName = findOxideColumn(data_opx.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDoublePreserveNaN(data_opx.(columnName));
        if isnan(value)
            nNan = nNan + 1;
            nanBuffer(nNan) = "Opx." + string(columnName);
        end
    end
end

% FeO is used when present; FeOt is used only when FeO is absent.
ironColumnName = findOxideColumn(data_opx.Properties.VariableNames, 'FeO');
if isempty(ironColumnName)
    ironColumnName = findOxideColumn(data_opx.Properties.VariableNames, 'FeOt');
end
if ~isempty(ironColumnName)
    ironValue = toScalarDoublePreserveNaN(data_opx.(ironColumnName));
    if isnan(ironValue)
        nNan = nNan + 1;
        nanBuffer(nNan) = "Opx." + string(ironColumnName);
    end
end

for i = 1:numel(liqOxides)
    oxide = liqOxides{i};
    columnName = findOxideColumn(data_liq.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDoublePreserveNaN(data_liq.(columnName));
        if isnan(value)
            nNan = nNan + 1;
            nanBuffer(nNan) = "Liquid." + string(columnName);
        end
    end
end

nanInputNames = nanBuffer(1:nNan);

end

function validateInputValues(data_opx, data_liq)
% validateInputValues
% Reject finite negative values and Inf in all present oxide inputs used by
% the calculation. Zero and NaN are allowed and handled without stopping.
% Liquid F and Cl are excluded because they do not enter the calculation.

opxOxides = {'SiO2', 'TiO2', 'Al2O3', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};
liqOxides = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', ...
    'SO3', 'Fe2O3'};

maxEntries = numel(opxOxides) + numel(liqOxides) + 1;
invalidBuffer = strings(maxEntries, 1);
nInvalid = 0;

for i = 1:numel(opxOxides)
    oxide = opxOxides{i};
    columnName = findOxideColumn(data_opx.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDoublePreserveNaN(data_opx.(columnName));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            invalidBuffer(nInvalid) = "Opx." + string(columnName);
        end
    end
end

% Validate only the iron column actually used by prepareOpxRow.
ironColumnName = findOxideColumn(data_opx.Properties.VariableNames, 'FeO');
if isempty(ironColumnName)
    ironColumnName = findOxideColumn(data_opx.Properties.VariableNames, 'FeOt');
end
if ~isempty(ironColumnName)
    ironValue = toScalarDoublePreserveNaN(data_opx.(ironColumnName));
    if isinf(ironValue) || (isfinite(ironValue) && ironValue < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Opx." + string(ironColumnName);
    end
end

for i = 1:numel(liqOxides)
    oxide = liqOxides{i};
    columnName = findOxideColumn(data_liq.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDoublePreserveNaN(data_liq.(columnName));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            invalidBuffer(nInvalid) = "Liquid." + string(columnName);
        end
    end
end

if nInvalid > 0
    invalidNames = invalidBuffer(1:nInvalid);
    error(['LesselPutirka2015OpxLiq: finite oxide inputs must be ' ...
           'greater than or equal to zero, and Inf is not permitted. ' ...
           'Invalid value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_liq, P_kbar, MWinfo)
% calcTemp
% Calculate Lessel and Putirka (2015) equation (4) for one selected
% Opx-Liquid pair. The same composition-dependent temperature is repeated
% once for every pressure value because equation (4) is pressure independent.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

% --- Opx cations on a 6-oxygen basis ---
opx = prepareOpxRow(data_opx, MWinfo);

% --- Liquid oxide values in wt.% ---
% Missing optional columns are treated as zero. Present NaN values remain
% NaN and propagate through the anhydrous cation-fraction normalization.
SiO2 = getOxideOptional(data_liq, 'SiO2', 0);
TiO2 = getOxideOptional(data_liq, 'TiO2', 0);
Al2O3 = getOxideOptional(data_liq, 'Al2O3', 0);
FeO = getOxideOptional(data_liq, 'FeO', 0);
MnO = getOxideOptional(data_liq, 'MnO', 0);
MgO = getOxideOptional(data_liq, 'MgO', 0);
CaO = getOxideOptional(data_liq, 'CaO', 0);
Na2O = getOxideOptional(data_liq, 'Na2O', 0);
K2O = getOxideOptional(data_liq, 'K2O', 0);
V2O3 = getOxideOptional(data_liq, 'V2O3', 0);
Cr2O3 = getOxideOptional(data_liq, 'Cr2O3', 0);
NiO = getOxideOptional(data_liq, 'NiO', 0);
P2O5 = getOxideOptional(data_liq, 'P2O5', 0);
SO3 = getOxideOptional(data_liq, 'SO3', 0);
Fe2O3 = getOxideOptional(data_liq, 'Fe2O3', 0);

% F and Cl are retained only as raw output values. They do not enter the
% cation-fraction denominator or any validation/NaN warning checks.
F = getOxideOptional(data_liq, 'F', NaN);
Cl = getOxideOptional(data_liq, 'Cl', NaN);

% --- Liquid cation proportions on an anhydrous basis ---
n = struct();
n.SiO2 = SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
n.TiO2 = TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
n.Al2O3 = Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n.FeO = FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
n.MnO = MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
n.MgO = MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
n.CaO = CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
n.Na2O = Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
n.K2O = K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
n.V2O3 = V2O3 .* MWinfo.Cat.V2O3 ./ MWinfo.MW.V2O3;
n.Cr2O3 = Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
n.NiO = NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;
n.P2O5 = P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;
n.SO3 = SO3 .* MWinfo.Cat.SO3 ./ MWinfo.MW.SO3;
n.Fe2O3 = Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

% F and Cl are anions and are intentionally excluded from cationTotal_liq.
cationTotal_liq = n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + ...
    n.MgO + n.CaO + n.Na2O + n.K2O + n.V2O3 + n.Cr2O3 + ...
    n.NiO + n.P2O5 + n.SO3 + n.Fe2O3;

Xliq = struct();
Xliq.SiO2 = n.SiO2 ./ cationTotal_liq;
Xliq.TiO2 = n.TiO2 ./ cationTotal_liq;
Xliq.AlO1_5 = n.Al2O3 ./ cationTotal_liq;
Xliq.FeO = n.FeO ./ cationTotal_liq;
Xliq.MnO = n.MnO ./ cationTotal_liq;
Xliq.MgO = n.MgO ./ cationTotal_liq;
Xliq.CaO = n.CaO ./ cationTotal_liq;
Xliq.NaO0_5 = n.Na2O ./ cationTotal_liq;
Xliq.KO0_5 = n.K2O ./ cationTotal_liq;

% --- Lessel and Putirka (2015), equation (4) terms ---
DMgO_OpxLiq = opx.XMg ./ Xliq.MgO;
AlIV_Opx = opx.XAlIV;
Na_Opx = opx.XNa;
Mg_Opx = opx.XMg;

% Prevent invalid logarithms and make all domain failures explicit as NaN.
equationDomainValid_scalar = ...
    isfinite(DMgO_OpxLiq) && DMgO_OpxLiq > 0 && ...
    isfinite(Na_Opx) && Na_Opx >= 0 && ...
    isfinite(Mg_Opx) && Mg_Opx > 0 && ...
    isfinite(AlIV_Opx) && AlIV_Opx > 0 && ...
    isfinite(Xliq.MgO) && Xliq.MgO > 0 && ...
    isfinite(Xliq.CaO) && Xliq.CaO > 0;

if equationDomainValid_scalar
    invT_scalar = ...
        2.234e-4 ...
        - 3.982e-6 .* DMgO_OpxLiq ...
        - 9.394e-4 .* Na_Opx ...
        - 2.058e-4 .* log(Xliq.MgO) ...
        + 3.564e-5 .* log(Xliq.CaO) ...
        + 2.186e-4 .* log(Mg_Opx) ...
        - 1.958e-5 .* log(AlIV_Opx) ...
        + 0.002298 .* (AlIV_Opx .^ 2);
else
    invT_scalar = NaN;
end

if isfinite(invT_scalar) && invT_scalar > 0
    TEq4_K_scalar = 1 ./ invT_scalar;
    TEq4_C_scalar = TEq4_K_scalar - 273.15;
else
    TEq4_K_scalar = NaN;
    TEq4_C_scalar = NaN;
    equationDomainValid_scalar = false;
end

% --- Pack outputs ---
% Every composition-dependent scalar is replicated to nP rows so downstream
% fixed-pressure and pressure-range workflows receive the same table shape.
row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.pressure_used_in_equation = false(nP, 1);

row.XSi_opx = repmat(opx.XSi, nP, 1);
row.XTi_opx = repmat(opx.XTi, nP, 1);
row.XAl_opx = repmat(opx.XAl, nP, 1);
row.XFe_opx = repmat(opx.XFe, nP, 1);
row.XMn_opx = repmat(opx.XMn, nP, 1);
row.XMg_opx = repmat(opx.XMg, nP, 1);
row.XCa_opx = repmat(opx.XCa, nP, 1);
row.XNa_opx = repmat(opx.XNa, nP, 1);
row.XK_opx = repmat(opx.XK, nP, 1);
row.XCr_opx = repmat(opx.XCr, nP, 1);
row.cationSum_opx = repmat(opx.cationSum, nP, 1);
row.XAlIV_opx = repmat(opx.XAlIV, nP, 1);
row.XAlVI_opx = repmat(opx.XAlVI, nP, 1);

row.SiO2_liq = repmat(SiO2, nP, 1);
row.TiO2_liq = repmat(TiO2, nP, 1);
row.Al2O3_liq = repmat(Al2O3, nP, 1);
row.FeO_liq = repmat(FeO, nP, 1);
row.MnO_liq = repmat(MnO, nP, 1);
row.MgO_liq = repmat(MgO, nP, 1);
row.CaO_liq = repmat(CaO, nP, 1);
row.Na2O_liq = repmat(Na2O, nP, 1);
row.K2O_liq = repmat(K2O, nP, 1);
row.V2O3_liq = repmat(V2O3, nP, 1);
row.Cr2O3_liq = repmat(Cr2O3, nP, 1);
row.NiO_liq = repmat(NiO, nP, 1);
row.P2O5_liq = repmat(P2O5, nP, 1);
row.SO3_liq = repmat(SO3, nP, 1);
row.F_liq = repmat(F, nP, 1);
row.Cl_liq = repmat(Cl, nP, 1);
row.Fe2O3_liq = repmat(Fe2O3, nP, 1);
row.cationTotal_liq = repmat(cationTotal_liq, nP, 1);

row.XSiO2_liq = repmat(Xliq.SiO2, nP, 1);
row.XTiO2_liq = repmat(Xliq.TiO2, nP, 1);
row.XAlO1_5_liq = repmat(Xliq.AlO1_5, nP, 1);
row.XFeO_liq = repmat(Xliq.FeO, nP, 1);
row.XMnO_liq = repmat(Xliq.MnO, nP, 1);
row.XMgO_liq = repmat(Xliq.MgO, nP, 1);
row.XCaO_liq = repmat(Xliq.CaO, nP, 1);
row.XNaO0_5_liq = repmat(Xliq.NaO0_5, nP, 1);
row.XKO0_5_liq = repmat(Xliq.KO0_5, nP, 1);

row.DMgO_OpxLiq = repmat(DMgO_OpxLiq, nP, 1);
row.AlIV_Opx = repmat(AlIV_Opx, nP, 1);
row.Na_Opx = repmat(Na_Opx, nP, 1);
row.Mg_Opx = repmat(Mg_Opx, nP, 1);
row.invT_Eq4 = repmat(invT_scalar, nP, 1);
row.equation_domain_valid = repmat(equationDomainValid_scalar, nP, 1);
row.TEq4_K = repmat(TEq4_K_scalar, nP, 1);
row.TEq4_C = repmat(TEq4_C_scalar, nP, 1);

% Ballhaus-style standardized temperature aliases for downstream launchers.
row.T_K = row.TEq4_K;
row.T_deg = row.TEq4_C;

end

function opx = prepareOpxRow(data_opx, MWinfo)
% prepareOpxRow
% Convert one selected Opx oxide analysis to cations per 6 oxygens. Present
% NaN values remain NaN. Missing optional columns follow the documented zero
% assumption.

SiO2 = getOxideRequired(data_opx, 'SiO2', 'Opx');
TiO2 = getOxideOptional(data_opx, 'TiO2', 0);
Al2O3 = getOxideOptional(data_opx, 'Al2O3', 0);
MnO = getOxideOptional(data_opx, 'MnO', 0);
MgO = getOxideRequired(data_opx, 'MgO', 'Opx');
CaO = getOxideOptional(data_opx, 'CaO', 0);
Na2O = getOxideOptional(data_opx, 'Na2O', 0);
K2O = getOxideOptional(data_opx, 'K2O', 0);
Cr2O3 = getOxideOptional(data_opx, 'Cr2O3', 0);

% FeO is preferred when the column exists. A present NaN is preserved. FeOt
% is used only when no FeO column exists.
FeO_name = findOxideColumn(data_opx.Properties.VariableNames, 'FeO');
if ~isempty(FeO_name)
    FeO = toScalarDoublePreserveNaN(data_opx.(FeO_name));
else
    FeOt_name = findOxideColumn(data_opx.Properties.VariableNames, 'FeOt');
    if isempty(FeOt_name)
        error('Selected Opx row must contain either FeO or FeOt.');
    end
    FeO = toScalarDoublePreserveNaN(data_opx.(FeOt_name));
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
    molProp.FeO + molProp.MnO + molProp.MgO + molProp.CaO + ...
    molProp.Na2O + molProp.K2O + 3 .* molProp.Cr2O3;

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

XAlIV = clampDerivedNonnegative(2 - XSi);
XAlVI = clampDerivedNonnegative(XAl - XAlIV);

opx = struct();
opx.XSi = XSi;
opx.XTi = XTi;
opx.XAl = XAl;
opx.XFe = XFe;
opx.XMn = XMn;
opx.XMg = XMg;
opx.XCa = XCa;
opx.XNa = XNa;
opx.XK = XK;
opx.XCr = XCr;
opx.cationSum = cationSum;
opx.XAlIV = XAlIV;
opx.XAlVI = XAlVI;

end

function value = clampDerivedNonnegative(value)
% clampDerivedNonnegative
% Clamp a finite derived component to zero when negative while preserving
% NaN. This applies only to derived components, not raw inputs.

if isnan(value)
    return;
end
if isfinite(value) && value < 0
    value = 0;
end

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Add common Liquid metadata fields, replicated to match the pressure rows.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if ismember('Index', variableNames)
    raw = data_liq.('Index');
    if isnumeric(raw) || islogical(raw)
        row.liq_Index = repmat(raw(1), nRows, 1);
    else
        row.liq_Index = repmat(string(raw(1)), nRows, 1);
    end
end
if ismember('Experiment', variableNames)
    rawExperiment = data_liq.('Experiment');
    row.liq_Experiment = repmat(string(rawExperiment(1)), nRows, 1);
end
if ismember('Citation', variableNames)
    rawCitation = data_liq.('Citation');
    row.liq_Citation = repmat(string(rawCitation(1)), nRows, 1);
end

end

function value = getOxideRequired(data_tbl, oxide, materialLabel)
% getOxideRequired
% Read one required oxide scalar. NaN is preserved. Missing columns stop the
% calculation because the required input cannot be reconstructed.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    error('%s table must contain variable: %s', materialLabel, oxide);
end
value = toScalarDoublePreserveNaN(data_tbl.(columnName));

end

function value = getOxideOptional(data_tbl, oxide, defaultValue)
% getOxideOptional
% Return defaultValue only when the oxide column is absent. A present NaN is
% preserved and propagated through the calculation.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    value = defaultValue;
else
    value = toScalarDoublePreserveNaN(data_tbl.(columnName));
end

end

function columnName = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match oxide columns after removing spaces, underscores, and hyphens. Both
% "oxide" and "oxideValue" forms are accepted.

canonicalNames = cell(size(variableNames));
for i = 1:numel(variableNames)
    canonicalNames{i} = canonicalizeName(variableNames{i});
end

canonicalOxide = canonicalizeName(oxide);
targets = {[canonicalOxide 'value'], canonicalOxide};

columnName = '';
for i = 1:numel(targets)
    idx = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(idx)
        columnName = variableNames{idx};
        return;
    end
end

end

function canonical = canonicalizeName(name)
% canonicalizeName
% Convert a variable name to the matching form used by findOxideColumn.

canonical = lower(char(string(name)));
canonical = strrep(canonical, ' ', '');
canonical = strrep(canonical, '_', '');
canonical = strrep(canonical, '-', '');

end

function value = toScalarDoublePreserveNaN(raw)
% toScalarDoublePreserveNaN
% Convert one selected table value to a numeric scalar. Empty, missing, or
% unparseable values become NaN. NaN is never replaced by a default value.

if isempty(raw)
    value = NaN;
    return;
end

if istable(raw)
    error('A selected oxide value cannot be a nested table.');
end

if iscell(raw)
    if numel(raw) ~= 1
        error('Selected oxide input must contain exactly one value.');
    end
    value = toScalarDoublePreserveNaN(raw{1});
    return;
end

if isnumeric(raw) || islogical(raw)
    if numel(raw) ~= 1
        error('Selected oxide input must be scalar.');
    end
    value = double(raw);
    return;
end

if isstring(raw)
    if numel(raw) ~= 1
        error('Selected oxide input must be scalar.');
    end
    if ismissing(raw)
        value = NaN;
    else
        value = str2double(raw);
    end
    return;
end

if ischar(raw)
    value = str2double(string(raw));
    return;
end

if iscategorical(raw)
    if numel(raw) ~= 1 || isundefined(raw)
        value = NaN;
    else
        value = str2double(string(raw));
    end
    return;
end

error('Unsupported data type for a selected oxide input: %s', class(raw));

end
