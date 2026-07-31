function results = LesselPutirka2015MgO(rawdata_struct, P_kbar, varargin)
% functions/+thermo/LesselPutirka2015MgO.m
% Tested with MATLAB R2024b
%
% Empirical liquid-MgO thermometer
% Lessel, J. and Putirka, K.D. (2015), Equation (17)
% American Mineralogist, 100, 2163-2171
% DOI: https://doi.org/10.2138/am-2015-4732
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one mineral analysis for traceability
% and one liquid analysis, then calculates temperature from liquid MgO using
% Equation (17) of Lessel and Putirka (2015):
%
%   T(degreeC) = 1011 + 29.8 * MgO_liq(wt%)
%
% Mineral chemistry is not used in Equation (17). The selected Cpx, Opx, or
% Olivine row is retained only as an identifier for compatibility with the
% existing interactive thermometer workflow.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to select another mineral row and appends the result to a
% single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL DATA RANGE AND APPLICATION NOTES
%
% Lessel and Putirka (2015) developed their thermobarometers specifically
% for martian igneous compositions, which generally have higher FeO and
% lower Al2O3 than terrestrial basalts. Application of Equation (17) to
% terrestrial magmas was not calibrated or validated in the paper
% (abstract and Introduction, pp. 2163-2164).
%
% Equation (17) is described as a simple empirical MgO thermometer on
% p. 2168 and as less precise than the mineral-liquid thermometers on
% pp. 2168-2169. The paper does not report an Equation (17)-specific R2,
% RMSE, SEE, calibration-data count, or formal calibration range.
%
% The complete martian experimental database used in the study spans:
%
%   Temperature : 950-1540 degreeC
%   Pressure    : 1 atm-2.3 GPa (approximately 0.001-23 kbar)
%   SiO2_liq    : 40.2-66.16 wt%
%   MgO_liq     : 0.62-24.52 wt%
%   FeO_liq     : 2.80-30.2 wt%
%   Al2O3_liq   : 2.97-20.5 wt%
%   Total alkalis: 0.19-6.77 wt%
%
% These study-wide ranges are reported on pp. 2164-2165. They are used in
% this implementation only as contextual screening limits and must not be
% interpreted as a formally reported calibration range specific to
% Equation (17).
%
% For the high-MgO Yamato 980459 composition, Equation (17) gives about
% 1550 degreeC. The authors state that this estimate is probably too high
% relative to the more precise thermometer of Putirka et al. (2007)
% (p. 2168). In contrast, Equation (17) reproduced the experimental result
% for a lower-temperature Gusev basalt reasonably well. The authors
% therefore only tentatively suggest that Equation (17) may be especially
% useful at lower temperatures and pressures (p. 2169).
%
% MgO must represent a silicate liquid or a defensible parental-liquid
% composition. A whole-rock MgO value should not be used for a cumulate,
% crystal-rich rock, or strongly fractionated/altered sample unless the
% whole-rock composition can reasonably be treated as a liquid. The paper
% explicitly makes this assumption when applying whole-rock compositions
% as liquids (p. 2168).
%
% Equation (17) contains no pressure term and provides no mineral-liquid
% equilibrium test. P_kbar is accepted only for compatibility with
% startThermoCalc_fixedP and startThermoCalc_rangeP. One output row is
% returned for every supplied pressure value, although the calculated
% temperature is identical at every pressure for a given liquid MgO value.
%
% This implementation issues non-stopping warnings using fprintf when:
%   1) input pressure lies outside the study-wide range of approximately
%      1 atm-2.3 GPa;
%   2) finite liquid MgO lies outside the study-wide range of
%      0.62-24.52 wt%;
%   3) finite calculated temperature lies outside the study-wide range of
%      950-1540 degreeC;
%   4) MgO is NaN; or
%   5) the calculated temperature is NaN or Inf.
%
% NaN MgO is retained and propagated through Equation (17); it is never
% replaced by zero. A finite liquid MgO value below zero is prohibited.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain at least one of the following mineral tables:
%   rawdata_struct.Cpx
%   rawdata_struct.Opx
%   rawdata_struct.Olivine
%
% The first available table is used in the order Cpx, Opx, Olivine. Its
% FIRST column is treated as an identifier ("data code") displayed in the
% selection dialog. The mineral table is used only for traceability.
%
% The liquid dataset selected through liquid.readLiquidExcel must contain:
%   MgO or MgO_value : liquid MgO in wt%
%
% Missing or explicitly NaN MgO remains NaN. Other liquid oxides are not
% used by Equation (17).
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Lessel and Putirka (2015), Equation (17), p. 2168:
%
%   T(degreeC) = 1011 + 29.8 * MgO_liq(wt%)
%   T(K)       = T(degreeC) + 273.15
%
% -------------------------------------------------------------------------
% Syntax:
%   results = LesselPutirka2015MgO(rawdata_struct, P_kbar)
%   results = LesselPutirka2015MgO(rawdata_struct, P_kbar, ...
%       'LiquidRow', liquidRow)
%
% Inputs:
%   rawdata_struct : struct containing at least one Cpx, Opx, or Olivine
%                    table (see above)
%   P_kbar         : pressure in kbar (finite, non-negative numeric scalar
%                    or vector); retained for interface compatibility but
%                    not used in Equation (17)
%
% Name-value option:
%   'LiquidRow'    : positive integer row number or [] (default []). If
%                    empty, row 1 of the selected liquid dataset is used.
%
% Output:
%   results : table containing one row per pressure value for every selected
%             mineral-liquid record. NaN and Inf results are retained.

%% Input validation
% Basic argument checks prevent silent failures caused by missing inputs or
% invalid pressure values. Both fixed-pressure and pressure-range launchers
% are supported by accepting a scalar or vector P_kbar input.
if nargin < 2
    error('LesselPutirka2015MgO requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve mineral dataset for traceability
% Equation (17) does not use mineral chemistry. One available mineral table
% is nevertheless retained to preserve the established selection interface
% and identify the record associated with each calculation.
disp('=== Step 1: Preparing mineral dataset for traceability ===');

[mineralField, dataset_min] = localPickMineralDataset(rawdata_struct);

disp('=== Preparing mineral dataset has been finished ===');

%% 2) Parse options and retrieve liquid dataset
% LiquidRow is optional because the existing liquid importer may return a
% dataset containing more than one row.
disp('=== Step 2: Preparing liquid dataset ===');

ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

[liqAll, metaLiq] = liquid.readLiquidExcel();
if ~istable(liqAll) || isempty(liqAll)
    error('Selected liquid dataset must be a non-empty table.');
end

if isempty(liquidRowOpt)
    idxLiq = 1;
    if height(liqAll) > 1
        fprintf(2, ['WARNING: The selected liquid dataset contains %d rows. ' ...
                    'Liquid row 1 will be used.\n'], height(liqAll));
    end
else
    idxLiq = liquidRowOpt;
    if idxLiq > height(liqAll)
        error(['Requested LiquidRow (%d) exceeds the number of rows in the ' ...
               'selected liquid dataset (%d).'], idxLiq, height(liqAll));
    end
end

selectedData_liq = liqAll(idxLiq, :);

disp(['Liquid selected: Row ' num2str(idxLiq)]);
disp('=== Preparing liquid dataset has been finished ===');

%% 3) Initialize output container and contextual screening limits
% Each calculation is stored temporarily as one table block. Repeated table
% concatenation inside the interactive loop is avoided. The cell buffer is
% preallocated and doubled only when its capacity is exhausted; all result
% blocks are concatenated once after the loop.
disp('=== Step 3: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Study-wide experimental database limits from Lessel and Putirka (2015),
% pp. 2164-2165. These are contextual limits, not an explicitly reported
% Equation (17)-specific calibration range.
studyT_min_degreeC = 950;
studyT_max_degreeC = 1540;
studyP_min_GPa = 0.0001;
studyP_max_GPa = 2.3;
studyMgO_min_wt = 0.62;
studyMgO_max_wt = 24.52;

P_GPa_input = P_kbar ./ 10;
pressureOutsideStudyRange = ...
    P_GPa_input < studyP_min_GPa | P_GPa_input > studyP_max_GPa;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 4-5) Interactive selection loop and calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp(['=== Step 4: Selecting a data code from the list (' mineralField ') ===']);

while true
    % ----- Mineral selection for traceability -----
    dataCodes_min = dataset_min{:, 1};

    [selectedIdx_min, ok] = listdlg( ...
        'PromptString', ['Please select the ' mineralField ...
        ' data you would like to use:'], ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_min, ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_min)
        disp('Selection canceled');
        break;
    end

    selectedCode_min = string(dataCodes_min(selectedIdx_min));
    disp([mineralField ' selected: ' char(selectedCode_min)]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_min = dataset_min(selectedIdx_min, :);

    % Extract only the variable used in Equation (17). Missing, malformed,
    % or explicitly NaN MgO is represented by NaN and is not replaced by 0.
    MgO_liq = localGetLiqOxOptional(selectedData_liq, 'MgO', NaN);

    % A finite negative MgO value is physically invalid and is prohibited.
    % Zero is allowed by the input check but will trigger the contextual
    % composition-range warning below.
    validateNonnegativeMgO(MgO_liq);

    nanInputNames = findNaNInputs(MgO_liq);

    row = calcTemp(selectedData_min, selectedData_liq, P_kbar);

    nRows = height(row);
    row.dataCode_mineral = repmat(selectedCode_min, nRows, 1);
    row.mineralType = repmat(string(mineralField), nRows, 1);
    row.dataRow_liq = repmat(idxLiq, nRows, 1);
    row = localAttachLiquidIDs(row, selectedData_liq);
    row = movevars(row, ...
        {'mineralType', 'dataCode_mineral', 'dataRow_liq'}, 'Before', 1);

    % Store the result in the preallocated block buffer. Capacity is doubled
    % only when required, rather than changing the results-table size during
    % every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(selectedCode_min) ': ' num2str(row.TEq17_C) ' degreeC']);
    else
        disp([char(selectedCode_min) ': ' num2str(row.TEq17_C(1)) ...
            ' to ' num2str(row.TEq17_C(end)) ' degreeC']);
    end

    % Pressure is not used in Equation (17). This warning only identifies
    % extrapolation beyond the study-wide experimental database and is
    % printed once for the current function call.
    if any(pressureOutsideStudyRange) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the study-wide experimental ' ...
             'data range reported by Lessel and Putirka (2015), pp. 2164-2165: ' ...
             'approximately 1 atm-2.3 GPa (0.001-23 kbar). ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g GPa. Pressure is not used in Equation (17), ' ...
             'and this is not an Equation (17)-specific calibration limit.\n'], ...
            sum(pressureOutsideStudyRange), numel(P_GPa_input), ...
            min(P_GPa_input), max(P_GPa_input));
        pressureWarningIssued = true;
    end

    % Warn when the finite MgO input lies outside the composition range of
    % the study-wide experimental database. Calculation is not stopped.
    if isfinite(MgO_liq) && ...
            (MgO_liq < studyMgO_min_wt || MgO_liq > studyMgO_max_wt)
        fprintf(2, ...
            ['WARNING: Liquid MgO for %s is outside the study-wide experimental ' ...
             'data range reported by Lessel and Putirka (2015), pp. 2164-2165: ' ...
             '0.62-24.52 wt%%; input MgO = %.6g wt%%. This is not an ' ...
             'Equation (17)-specific calibration limit.\n'], ...
            char(selectedCode_min), MgO_liq);
    end

    % Warn when finite temperature results lie outside the study-wide
    % experimental temperature range. NaN and Inf are reported separately.
    finiteTemperature = isfinite(row.TEq17_C);
    temperatureOutsideStudyRange = finiteTemperature & ...
        (row.TEq17_C < studyT_min_degreeC | ...
         row.TEq17_C > studyT_max_degreeC);

    if any(temperatureOutsideStudyRange)
        finiteValues = row.TEq17_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the study-wide ' ...
             'experimental data range reported by Lessel and Putirka (2015), ' ...
             'pp. 2164-2165: 950-1540 degreeC. %d of %d finite temperature ' ...
             'point(s) are outside the range; calculated finite range = ' ...
             '%.4g-%.4g degreeC for %s. This is not an Equation (17)-specific ' ...
             'calibration limit.\n'], ...
            sum(temperatureOutsideStudyRange), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), char(selectedCode_min));
    end

    % Report the exact calculation inputs that contain NaN. The NaN remains
    % unchanged and propagates through Equation (17).
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         NaN was retained and propagated; the calculation was ' ...
             'not stopped.\n'], ...
            char(selectedCode_min), char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite temperature results, including those
    % caused by NaN MgO or another unexpected numerical condition.
    invalidTemperature = ~isfinite(row.TEq17_C);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(selectedCode_min), sum(invalidTemperature), ...
            numel(row.TEq17_C), sum(isnan(row.TEq17_C)), ...
            sum(isinf(row.TEq17_C)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection (same liquid dataset)?', ...
        'LesselPutirka2015MgO', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered table blocks once. Return an empty table if the
% user canceled before performing any calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct('liquid', metaLiq);
disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function [mineralField, dataset_min] = localPickMineralDataset(rawdata_struct)
% localPickMineralDataset
% Return the first available mineral table in the order Cpx, Opx, Olivine.
% The selected table is used only for interactive traceability.

candidateFields = {'Cpx', 'Opx', 'Olivine'};
mineralField = '';
dataset_min = table();

for i = 1:numel(candidateFields)
    fieldName = candidateFields{i};
    if isfield(rawdata_struct, fieldName) && ...
            istable(rawdata_struct.(fieldName)) && ...
            ~isempty(rawdata_struct.(fieldName))
        mineralField = fieldName;
        dataset_min = rawdata_struct.(fieldName);
        break;
    end
end

if isempty(mineralField)
    error(['rawdata_struct must contain at least one non-empty mineral table: ' ...
           'rawdata_struct.Cpx, rawdata_struct.Opx, or rawdata_struct.Olivine.']);
end

end

function nanInputNames = findNaNInputs(MgO_liq)
% findNaNInputs
% Return the names of Equation (17) input variables that contain NaN.

if any(isnan(MgO_liq(:)))
    nanInputNames = "Liquid.MgO";
else
    nanInputNames = strings(0, 1);
end

end

function validateNonnegativeMgO(MgO_liq)
% validateNonnegativeMgO
% Stop when liquid MgO is below zero. NaN is intentionally accepted so that
% it remains missing, propagates through the calculation, and is reported by
% the non-stopping warning messages.

if any(MgO_liq(:) < 0)
    error(['LesselPutirka2015MgO: liquid MgO used in Equation (17) ' ...
           'must be greater than or equal to zero. Input MgO = %.6g wt%%.'], ...
        MgO_liq);
end

end

function row = calcTemp(data_min, data_liq, P_kbar)
% calcTemp
% Calculate Equation (17) for one liquid MgO value and every supplied
% pressure. Pressure and the mineral row are retained only for traceability.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

% Read the liquid composition columns retained by the original function.
% Only MgO enters Equation (17). Missing and explicitly NaN values remain
% NaN rather than being silently replaced by zero.
SiO2  = localGetLiqOxOptional(data_liq, 'SiO2', NaN);
TiO2  = localGetLiqOxOptional(data_liq, 'TiO2', NaN);
Al2O3 = localGetLiqOxOptional(data_liq, 'Al2O3', NaN);
FeO   = localGetLiqOxOptional(data_liq, 'FeO', NaN);
MnO   = localGetLiqOxOptional(data_liq, 'MnO', NaN);
MgO   = localGetLiqOxOptional(data_liq, 'MgO', NaN);
CaO   = localGetLiqOxOptional(data_liq, 'CaO', NaN);
Na2O  = localGetLiqOxOptional(data_liq, 'Na2O', NaN);
K2O   = localGetLiqOxOptional(data_liq, 'K2O', NaN);
Cr2O3 = localGetLiqOxOptional(data_liq, 'Cr2O3', NaN);
NiO   = localGetLiqOxOptional(data_liq, 'NiO', NaN);
Fe2O3 = localGetLiqOxOptional(data_liq, 'Fe2O3', NaN);
H2O   = localGetLiqOxOptional(data_liq, 'H2O', NaN);

% Replicate each pressure-independent liquid value once so every output
% variable has nP rows. NaN remains NaN during these operations.
SiO2_liq_vector  = repmat(SiO2, nP, 1);
TiO2_liq_vector  = repmat(TiO2, nP, 1);
Al2O3_liq_vector = repmat(Al2O3, nP, 1);
FeO_liq_vector   = repmat(FeO, nP, 1);
MnO_liq_vector   = repmat(MnO, nP, 1);
MgO_liq_vector   = repmat(MgO, nP, 1);
CaO_liq_vector   = repmat(CaO, nP, 1);
Na2O_liq_vector  = repmat(Na2O, nP, 1);
K2O_liq_vector   = repmat(K2O, nP, 1);
Cr2O3_liq_vector = repmat(Cr2O3, nP, 1);
NiO_liq_vector   = repmat(NiO, nP, 1);
Fe2O3_liq_vector = repmat(Fe2O3, nP, 1);
H2O_liq_vector   = repmat(H2O, nP, 1);

% Lessel and Putirka (2015), Equation (17), p. 2168.
TEq17_C = 1011 + 29.8 .* MgO_liq_vector;
TEq17_K = TEq17_C + 273.15;

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.SiO2_liq = SiO2_liq_vector;
row.TiO2_liq = TiO2_liq_vector;
row.Al2O3_liq = Al2O3_liq_vector;
row.FeO_liq = FeO_liq_vector;
row.MnO_liq = MnO_liq_vector;
row.MgO_liq = MgO_liq_vector;
row.CaO_liq = CaO_liq_vector;
row.Na2O_liq = Na2O_liq_vector;
row.K2O_liq = K2O_liq_vector;
row.Cr2O3_liq = Cr2O3_liq_vector;
row.NiO_liq = NiO_liq_vector;
row.Fe2O3_liq = Fe2O3_liq_vector;
row.H2O_liq = H2O_liq_vector;
row.TEq17_C = TEq17_C;
row.TEq17_K = TEq17_K;

% Mineral information is not used by the equation and is retained only to
% make that behavior explicit in the output.
row.trace_hasMineralRow = repmat(~isempty(data_min), nP, 1);
row.trace_mineralRowHeight = repmat(height(data_min), nP, 1);

end

function row = localAttachLiquidIDs(row, data_liq)
% localAttachLiquidIDs
% Replicate optional liquid identifiers to match the pressure-vector height.

variableNames = data_liq.Properties.VariableNames;
nRows = height(row);

if any(strcmp(variableNames, 'Index'))
    identifierValue = data_liq.('Index');
    row.liq_Index = repmat(identifierValue(1, :), nRows, 1);
end
if any(strcmp(variableNames, 'Experiment'))
    identifierValue = string(data_liq.('Experiment'));
    row.liq_Experiment = repmat(identifierValue(1), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    identifierValue = string(data_liq.('Citation'));
    row.liq_Citation = repmat(identifierValue(1), nRows, 1);
end

end

function value = localGetLiqOxOptional(data_liq, oxide, defaultValue)
% localGetLiqOxOptional
% Read a scalar liquid-oxide value without replacing NaN by zero.

columnName = localFindOxideColumn(data_liq.Properties.VariableNames, oxide);

if isempty(columnName)
    value = defaultValue;
    return;
end

rawValue = data_liq.(columnName);
value = localToScalarDouble(rawValue, defaultValue);

end

function columnName = localFindOxideColumn(variableNames, oxide)
% localFindOxideColumn
% Match either an oxide name or the corresponding oxide_value name while
% ignoring spaces, underscores, hyphens, and letter case.

canonicalNames = cell(size(variableNames));
for i = 1:numel(variableNames)
    textValue = lower(variableNames{i});
    textValue = strrep(textValue, ' ', '');
    textValue = strrep(textValue, '_', '');
    textValue = strrep(textValue, '-', '');
    canonicalNames{i} = textValue;
end

canonicalOxide = lower(oxide);
canonicalOxide = strrep(canonicalOxide, ' ', '');
canonicalOxide = strrep(canonicalOxide, '_', '');
canonicalOxide = strrep(canonicalOxide, '-', '');

targets = {[canonicalOxide 'value'], canonicalOxide};
columnName = '';

for i = 1:numel(targets)
    matchIndex = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(matchIndex)
        columnName = variableNames{matchIndex};
        return;
    end
end

end

function value = localToScalarDouble(rawValue, defaultValue)
% localToScalarDouble
% Convert the first table entry to double. Numeric NaN and textual "NaN"
% are preserved as NaN rather than replaced with the default or with zero.

value = defaultValue;

if isempty(rawValue)
    return;
end

if isnumeric(rawValue) || islogical(rawValue)
    value = double(rawValue(1));
    return;
end

if isstring(rawValue)
    if ismissing(rawValue(1))
        return;
    end
    convertedValue = str2double(rawValue(1));
    if ~isnan(convertedValue) || strcmpi(strtrim(rawValue(1)), "NaN")
        value = convertedValue;
    end
    return;
end

if ischar(rawValue)
    convertedValue = str2double(string(rawValue));
    if ~isnan(convertedValue) || strcmpi(strtrim(string(rawValue)), "NaN")
        value = convertedValue;
    end
    return;
end

if iscell(rawValue)
    if isempty(rawValue{1})
        return;
    end
    firstValue = rawValue{1};
    if isnumeric(firstValue) || islogical(firstValue)
        value = double(firstValue(1));
        return;
    end
    if isstring(firstValue) || ischar(firstValue)
        textValue = string(firstValue);
        if ismissing(textValue(1))
            return;
        end
        convertedValue = str2double(textValue(1));
        if ~isnan(convertedValue) || strcmpi(strtrim(textValue(1)), "NaN")
            value = convertedValue;
        end
    end
end

end
