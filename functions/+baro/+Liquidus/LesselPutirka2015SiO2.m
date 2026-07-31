function results = LesselPutirka2015SiO2(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Pyroxene/LesselPutirka2015SiO2.m
% Tested with MATLAB R2024b
%
% Si-activity barometer for martian igneous compositions
% Lessel, J. and Putirka, K. (2015), Equation (8)
% American Mineralogist, 100, 2163-2171
% DOI: https://doi.org/10.2138/am-2015-4732
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function calculates pressure from liquid composition and independently
% supplied temperature using the Si-activity barometer of Lessel and Putirka
% (2015), Equation (8).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected mineral-context row and liquid row,
% one output row is returned for every input temperature value.
%
% IMPORTANT: Equation (8) itself uses LIQUID COMPOSITION ONLY. The selected
% Cpx, Opx, or Olivine row is retained only as a contextual identifier and is
% not used numerically in the pressure equation. Lessel and Putirka (2015)
% state that the Si-activity calibration uses liquid compositions from
% experiments in which only olivine and orthopyroxene had formed (p. 2166).
% A selected Cpx row therefore does not demonstrate that the phase-assemblage
% condition of the calibration is satisfied.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% IMPORTANT: This barometer was calibrated specifically for MARTIAN igneous
% compositions, which are generally richer in FeO and poorer in Al2O3 than
% common terrestrial basalts (pp. 2163-2164). Application to ordinary
% terrestrial basaltic compositions is an extrapolation beyond the intended
% compositional domain.
%
% The complete experimental database used in the study spans approximately
% (Methodology and Table 1, pp. 2164-2165):
%
%   Temperature       : 950-1540 degreeC
%   Pressure          : 1 atm-2.3 GPa (approximately 0.001-23 kbar)
%   SiO2              : 40.2-66.16 wt%
%   MgO               : 0.62-24.52 wt%
%   FeO               : 2.80-30.2 wt%
%   Al2O3             : 2.97-20.5 wt%
%   Na2O + K2O        : 0.19-6.77 wt%
%
% These are STUDY-WIDE experimental ranges. The exact minimum and maximum
% values of the 29 experiments used specifically for Equation (8) are not
% tabulated separately in the main article. Accordingly, the limits above are
% used only as non-stopping warning envelopes and must not be interpreted as
% rigorously defined Equation (8)-specific calibration limits.
%
% Equation (8) was based on 29 experiments: 22 calibration data and 7
% independent test data (pp. 2164-2165). It reproduces calibration pressures
% with R^2 = 0.91 and RMSE = +/-0.16 GPa, and predicts test pressures with
% R^2 = 0.92 and RMSE = +/-0.18 GPa (p. 2166; Fig. 9 on p. 2167).
%
% The calibration uses liquids from experiments in which only olivine and
% orthopyroxene had formed (p. 2166). The liquid composition should therefore
% represent the appropriate Ol-Opx-saturated magmatic stage. Whole-rock
% composition should be used as a liquid proxy only when petrologically
% justified; cumulate, mixed, altered, or strongly fractionated whole-rock
% compositions may yield misleading pressures.
%
% Lessel and Putirka (2015) applied Equation (8) only to P-T estimates that
% passed mineral-liquid equilibrium tests (p. 2168). For martian Ol-liquid
% pairs, they recommend KD(Fe-Mg)Ol-liq = 0.36 +/- 0.02 (p. 2167). Opx-liquid
% equilibrium is composition dependent and should not be evaluated using a
% single universal mean KD value (p. 2167).
%
% The SiO2 activity expression includes an Al2O3 correction. Lower Al2O3
% decreases the resulting pressure relative to formulations without this
% correction (pp. 2163 and 2168). Accurate SiO2, Al2O3, TiO2, MnO, MgO,
% Na2O, and K2O analyses are therefore particularly important.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 950-1540 degreeC,
%   2) finite calculated pressure is outside approximately 0.001-23 kbar,
%   3) finite liquid composition is outside the study-wide ranges above,
%   4) a Cpx row is selected as the contextual mineral row,
%   5) a calculation input contains NaN, or
%   6) a calculated pressure is NaN or Inf.
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if zero makes a ratio, inverse power,
% or logarithm undefined, the resulting NaN or Inf is retained and reported.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain at least one of the following tables:
%   rawdata_struct.Cpx
%   rawdata_struct.Opx
%   rawdata_struct.Olivine
%
% The FIRST column of the selected mineral table is treated as an identifier
% ("data code") displayed in the selection dialog. The selected mineral row
% is not used by Equation (8); it is retained only for traceability.
%
% Liquid composition is read through liquid.readLiquidExcel(). The liquid
% table should contain wt% oxide columns for the following calculation
% inputs:
%
%   SiO2, TiO2, Al2O3, FeO or FeOt, MnO, MgO, CaO, Na2O, K2O,
%   V2O3, Cr2O3, NiO, P2O5, SO3, and Fe2O3
%
% Missing or non-numeric values are represented by NaN and are never replaced
% by zero. Because cation fractions are normalized to the total cations, NaN
% in any included oxide propagates through cationTotal_liq and pressure.
%
% F and Cl may be retained in the output for reference, but they are excluded
% from cationTotal_liq and excluded from the NaN-input warning, as required by
% this implementation.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (Lessel and Putirka, 2015, Eq. 8; p. 2166)
%
%   P(GPa) =
%       -11.16
%       -184.3*XMnO_liq
%       -5.268*XMgO_liq
%       +21.18*aSiO2_liq
%       +4.961*(XNaO0.5_liq + XKO0.5_liq)
%       -3.577e-3*T(K)*ln(aSiO2_liq)
%
% where the SiO2 activity is calculated following Beattie (1993):
%
%   aSiO2_liq =
%       (3*XSiO2_liq)^(-2)
%       *(1 - XAlO1.5_liq)^(7/2)
%       *(1 - XTiO2_liq)^7
%
% Liquid components are cation fractions, following Putirka (2008). F and Cl
% are not included in the cation-fraction denominator in this implementation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = LesselPutirka2015SiO2(rawdata_struct, T_degreeC)
%   results = LesselPutirka2015SiO2(rawdata_struct, T_degreeC, ...
%       'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing at least one Cpx, Opx, or Olivine table
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Optional name-value input:
%   LiquidRow      : positive integer row number in the selected liquid table.
%                    If omitted, row 1 is used.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected mineral-context row.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('LesselPutirka2015SiO2 requires (rawdata_struct, T_degreeC).');
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

ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOption = ip.Results.LiquidRow;

%% 1) Retrieve mineral-context and liquid datasets
% The selected mineral row is retained only as an identifier. Equation (8)
% uses liquid composition and temperature only.
disp('=== Step 1: Preparing mineral-context and liquid datasets ===');

[mineralField, dataset_min] = pickMineralDataset(rawdata_struct);

MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();
if isempty(liqAll)
    error('Selected liquid dataset is empty.');
end

if isempty(liquidRowOption)
    if height(liqAll) > 1
        fprintf(2, ...
            ['WARNING: Liquid dataset has multiple rows (%d). Row 1 will be used. ' ...
             'Specify ''LiquidRow'' to use another row.\n'], ...
            height(liqAll));
    end
    selectedIdx_liq = 1;
else
    selectedIdx_liq = liquidRowOption;
    if selectedIdx_liq > height(liqAll)
        error(['Requested LiquidRow (%d) exceeds the number of rows in the ' ...
               'selected liquid dataset (%d).'], ...
            selectedIdx_liq, height(liqAll));
    end
end

disp('=== Preparing mineral-context and liquid datasets has been finished ===');

%% 2) Initialize output container and warning envelopes
% Store each calculation result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Study-wide experimental ranges. Exact Equation (8)-specific minima and
% maxima are not tabulated separately in the main article.
studyT_min_degreeC = 950;
studyT_max_degreeC = 1540;
studyP_min_kbar = 0.001;
studyP_max_kbar = 23;

temperatureOutsideStudyRange = isfinite(T_degreeC) & ...
    (T_degreeC < studyT_min_degreeC | T_degreeC > studyT_max_degreeC);
temperatureWarningIssued = false;
phaseAssemblageWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp(['=== Step 3: Selecting a data code from the list (' mineralField ') ===']);

while true
    % ----- Mineral-context selection -----
    % The first table column is used only as the displayed data identifier.
    dataCodes_min = dataset_min{:, 1};

    [selectedIdx_min, ok] = listdlg( ...
        'PromptString', ['Please select the ' mineralField ...
        ' data you would like to retain as context:'], ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_min)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_min)
        disp('Selection canceled');
        break;
    end

    selectedCode_min = dataCodes_min(selectedIdx_min);
    disp([mineralField ' selected: ' char(string(selectedCode_min))]);
    disp(['Liquid selected: Row ' num2str(selectedIdx_liq)]);

    % ----- Calculation -----
    disp('=== Step 4: Checking inputs and calculating pressure ===');

    selectedData_min = dataset_min(selectedIdx_min, :);
    selectedData_liq = liqAll(selectedIdx_liq, :);

    liquidComposition = prepareLiquidComposition(selectedData_liq);

    % Check NaN only in variables used by the pressure calculation and its
    % cation-fraction denominator. F and Cl are intentionally excluded.
    nanInputNames = findNaNInputs(liquidComposition, T_degreeC);

    % Reject Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(liquidComposition);

    row = calcPressure(selectedData_min, liquidComposition, ...
        T_degreeC, MWinfo);

    % Repeat identifiers for all temperatures in the current calculation.
    nRows = height(row);
    row.mineralType = repmat(string(mineralField), nRows, 1);
    row.dataCode_mineral = repmat(string(selectedCode_min), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, ...
        {'mineralType', 'dataCode_mineral', 'dataRow_liq'}, 'Before', 1);

    % Store one block per selected contextual mineral row. Expand the cell
    % buffer only when its current capacity has been exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_min)) ' & Liquid Row ' ...
            num2str(selectedIdx_liq) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_min)) ' & Liquid Row ' ...
            num2str(selectedIdx_liq) ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Warn once when Cpx is retained as the contextual mineral row. Equation
    % (8) was calibrated using liquids for which only Ol and Opx had formed.
    if strcmp(mineralField, 'Cpx') && ~phaseAssemblageWarningIssued
        fprintf(2, ...
            ['WARNING: A Cpx row was selected only as a contextual identifier. ' ...
             'Lessel and Putirka (2015) state that Equation (8) uses liquid ' ...
             'compositions from experiments in which only olivine and ' ...
             'orthopyroxene had formed (p. 2166). The selected Cpx composition ' ...
             'is not used in the pressure calculation.\n']);
        phaseAssemblageWarningIssued = true;
    end

    % Input temperature is common to all selected rows, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideStudyRange) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the study-wide experimental ' ...
             'range of 950-1540 degreeC reported by Lessel and Putirka (2015; ' ...
             'pp. 2164-2165). The exact Equation (8)-specific range is not ' ...
             'tabulated separately. %d of %d finite temperature point(s) are ' ...
             'outside the range; input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideStudyRange), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures fall outside the study-wide
    % pressure envelope. This is not a rigorously defined Eq. (8) range.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideStudyRange = finitePressure & ...
        (row.P_kbar < studyP_min_kbar | row.P_kbar > studyP_max_kbar);

    if any(pressureOutsideStudyRange)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the approximately ' ...
             '0.001-23 kbar study-wide experimental envelope reported by ' ...
             'Lessel and Putirka (2015; pp. 2164-2165). The exact Equation ' ...
             '(8)-specific range is not tabulated separately. %d of %d finite ' ...
             'pressure point(s) are outside the envelope; calculated finite ' ...
             'range = %.4g-%.4g kbar for %s & Liquid Row %d.\n'], ...
            sum(pressureOutsideStudyRange), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_min)), ...
            selectedIdx_liq);
    end

    % Warn when the finite liquid composition lies outside the study-wide
    % compositional ranges. NaN inputs are handled separately below.
    printCompositionRangeWarnings(liquidComposition, ...
        selectedCode_min, selectedIdx_liq);

    % List exact calculation inputs containing NaN. F and Cl are intentionally
    % absent from this list and do not affect cationTotal_liq.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s & ' ...
             'Liquid Row %d: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'pressure may remain NaN. F and Cl are excluded from both ' ...
             'cationTotal_liq and this NaN warning.\n'], ...
            char(string(selectedCode_min)), ...
            selectedIdx_liq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s & ' ...
             'Liquid Row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_min)), ...
            selectedIdx_liq, ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressures are retained for diagnosis but are outside
    % the physical and experimental domain.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & ' ...
             'Liquid Row %d (%d of %d points). The values were retained for ' ...
             'diagnostic purposes.\n'], ...
            char(string(selectedCode_min)), ...
            selectedIdx_liq, ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another mineral-context selection (same Liquid row)?', ...
        'LesselPutirka2015SiO2', ...
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
function [mineralField, dataset_min] = pickMineralDataset(rawdata_struct)
% pickMineralDataset
% Select the first available contextual mineral table using the original
% priority order. Equation (8) does not use mineral composition numerically.

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
           'rawdata_struct.Cpx, rawdata_struct.Opx, or rawdata_struct.Olivine']);
end

end

function liquidComposition = prepareLiquidComposition(data_liq)
% prepareLiquidComposition
% Extract one-row liquid oxide data without replacing NaN by zero. FeO is
% used when available; FeOt is used only when FeO is NaN. F and Cl are read
% for reference but excluded from the cation total and NaN-input checks.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liquidComposition = struct();
liquidComposition.SiO2 = getLiquidOxide(data_liq, 'SiO2');
liquidComposition.TiO2 = getLiquidOxide(data_liq, 'TiO2');
liquidComposition.Al2O3 = getLiquidOxide(data_liq, 'Al2O3');

FeO_direct = getLiquidOxide(data_liq, 'FeO');
FeOt = getLiquidOxide(data_liq, 'FeOt');
if isnan(FeO_direct)
    liquidComposition.FeO = FeOt;
    liquidComposition.FeO_source = "FeOt";
else
    liquidComposition.FeO = FeO_direct;
    liquidComposition.FeO_source = "FeO";
end

liquidComposition.MnO = getLiquidOxide(data_liq, 'MnO');
liquidComposition.MgO = getLiquidOxide(data_liq, 'MgO');
liquidComposition.CaO = getLiquidOxide(data_liq, 'CaO');
liquidComposition.Na2O = getLiquidOxide(data_liq, 'Na2O');
liquidComposition.K2O = getLiquidOxide(data_liq, 'K2O');
liquidComposition.V2O3 = getLiquidOxide(data_liq, 'V2O3');
liquidComposition.Cr2O3 = getLiquidOxide(data_liq, 'Cr2O3');
liquidComposition.NiO = getLiquidOxide(data_liq, 'NiO');
liquidComposition.P2O5 = getLiquidOxide(data_liq, 'P2O5');
liquidComposition.SO3 = getLiquidOxide(data_liq, 'SO3');
liquidComposition.Fe2O3 = getLiquidOxide(data_liq, 'Fe2O3');

% F and Cl are retained only for reference. They are not calculation inputs.
liquidComposition.F = getLiquidOxide(data_liq, 'F');
liquidComposition.Cl = getLiquidOxide(data_liq, 'Cl');

end

function nanInputNames = findNaNInputs(liquidComposition, T_degreeC)
% findNaNInputs
% Return names of pressure-calculation inputs containing NaN. F and Cl are
% intentionally excluded. NaN values are not changed and do not stop the
% calculation.

calculationFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', ...
    'Fe2O3'};

maxNames = numel(calculationFields) + 1;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

for i = 1:numel(calculationFields)
    fieldName = calculationFields{i};
    fieldValue = liquidComposition.(fieldName);
    if isnan(fieldValue)
        nNanInputs = nNanInputs + 1;
        if strcmp(fieldName, 'FeO')
            nanInputBuffer(nNanInputs) = "Liq.FeO/FeOt";
        else
            nanInputBuffer(nNanInputs) = "Liq." + string(fieldName);
        end
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(liquidComposition)
% validateNonNegativeInputs
% Reject Inf and finite negative values in variables used by the pressure
% calculation or cation total. Zero and NaN are intentionally allowed and
% retained. F and Cl are excluded because they are not calculation inputs.

calculationFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', ...
    'Fe2O3'};

invalidInputBuffer = strings(numel(calculationFields), 1);
nInvalidInputs = 0;

for i = 1:numel(calculationFields)
    fieldName = calculationFields{i};
    fieldValue = liquidComposition.(fieldName);
    if isinf(fieldValue) || (isfinite(fieldValue) && fieldValue < 0)
        nInvalidInputs = nInvalidInputs + 1;
        if strcmp(fieldName, 'FeO')
            invalidInputBuffer(nInvalidInputs) = "Liq.FeO/FeOt";
        else
            invalidInputBuffer(nInvalidInputs) = "Liq." + string(fieldName);
        end
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['LesselPutirka2015SiO2: calculation inputs must be ' ...
           'non-negative. NaN is allowed, but Inf and finite negative ' ...
           'value(s) are prohibited. Invalid input(s): ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_min, liquidComposition, T_degreeC, MWinfo)
% calcPressure
% Compute pressure for one liquid row at one or more input temperatures.
% The selected mineral row is used only for traceability. NaN values are not
% replaced and propagate naturally through the cation total and equation.
%
% Inputs:
%   data_min            : 1-row contextual mineral table
%   liquidComposition   : scalar oxide-composition struct
%   T_degreeC           : scalar or vector temperature in degreeC
%   MWinfo              : molar-weight and cation-number information
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Convert wt% oxides to cation amounts. F and Cl are intentionally omitted.
n = struct();
n.SiO2 = liquidComposition.SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
n.TiO2 = liquidComposition.TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
n.Al2O3 = liquidComposition.Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n.FeO = liquidComposition.FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
n.MnO = liquidComposition.MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
n.MgO = liquidComposition.MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
n.CaO = liquidComposition.CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
n.Na2O = liquidComposition.Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
n.K2O = liquidComposition.K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
n.V2O3 = liquidComposition.V2O3 .* MWinfo.Cat.V2O3 ./ MWinfo.MW.V2O3;
n.Cr2O3 = liquidComposition.Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
n.NiO = liquidComposition.NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;
n.P2O5 = liquidComposition.P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;
n.SO3 = liquidComposition.SO3 .* MWinfo.Cat.SO3 ./ MWinfo.MW.SO3;
n.Fe2O3 = liquidComposition.Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

% Do not use sum(...,'omitnan'): any NaN input must propagate to the total.
cationTotal_liq_scalar = ...
    n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + n.MgO + ...
    n.CaO + n.Na2O + n.K2O + n.V2O3 + n.Cr2O3 + n.NiO + ...
    n.P2O5 + n.SO3 + n.Fe2O3;

% Liquid cation fractions. A zero or NaN total naturally produces NaN/Inf.
XSiO2_scalar = n.SiO2 ./ cationTotal_liq_scalar;
XTiO2_scalar = n.TiO2 ./ cationTotal_liq_scalar;
XAlO1_5_scalar = n.Al2O3 ./ cationTotal_liq_scalar;
XMnO_scalar = n.MnO ./ cationTotal_liq_scalar;
XMgO_scalar = n.MgO ./ cationTotal_liq_scalar;
XNaO0_5_scalar = n.Na2O ./ cationTotal_liq_scalar;
XKO0_5_scalar = n.K2O ./ cationTotal_liq_scalar;

% Beattie (1993) SiO2 activity expression as printed in Lessel and Putirka
% (2015), Eq. (8), p. 2166. The inverse-square term is (3*XSiO2)^(-2).
aSiO2_liq_scalar = ...
    ((3 .* XSiO2_scalar) .^ (-2)) .* ...
    ((1 - XAlO1_5_scalar) .^ (7 ./ 2)) .* ...
    ((1 - XTiO2_scalar) .^ 7);

% Expand composition-dependent scalars to the temperature-vector length.
SiO2_liq = repmat(liquidComposition.SiO2, nT, 1);
TiO2_liq = repmat(liquidComposition.TiO2, nT, 1);
Al2O3_liq = repmat(liquidComposition.Al2O3, nT, 1);
FeO_liq = repmat(liquidComposition.FeO, nT, 1);
MnO_liq = repmat(liquidComposition.MnO, nT, 1);
MgO_liq = repmat(liquidComposition.MgO, nT, 1);
CaO_liq = repmat(liquidComposition.CaO, nT, 1);
Na2O_liq = repmat(liquidComposition.Na2O, nT, 1);
K2O_liq = repmat(liquidComposition.K2O, nT, 1);
V2O3_liq = repmat(liquidComposition.V2O3, nT, 1);
Cr2O3_liq = repmat(liquidComposition.Cr2O3, nT, 1);
NiO_liq = repmat(liquidComposition.NiO, nT, 1);
P2O5_liq = repmat(liquidComposition.P2O5, nT, 1);
SO3_liq = repmat(liquidComposition.SO3, nT, 1);
Fe2O3_liq = repmat(liquidComposition.Fe2O3, nT, 1);
F_liq = repmat(liquidComposition.F, nT, 1);
Cl_liq = repmat(liquidComposition.Cl, nT, 1);
FeO_source = repmat(string(liquidComposition.FeO_source), nT, 1);

cationTotal_liq = repmat(cationTotal_liq_scalar, nT, 1);
XSiO2_liq = repmat(XSiO2_scalar, nT, 1);
XTiO2_liq = repmat(XTiO2_scalar, nT, 1);
XAlO1_5_liq = repmat(XAlO1_5_scalar, nT, 1);
XMnO_liq = repmat(XMnO_scalar, nT, 1);
XMgO_liq = repmat(XMgO_scalar, nT, 1);
XNaO0_5_liq = repmat(XNaO0_5_scalar, nT, 1);
XKO0_5_liq = repmat(XKO0_5_scalar, nT, 1);
aSiO2_liq = repmat(aSiO2_liq_scalar, nT, 1);

% Pressure calculation. No finite-value guard is used: NaN and Inf values
% propagate and remain available for diagnosis in the output table.
P_GPa = ...
    -11.16 ...
    - 184.3 .* XMnO_liq ...
    - 5.268 .* XMgO_liq ...
    + 21.18 .* aSiO2_liq ...
    + 4.961 .* (XNaO0_5_liq + XKO0_5_liq) ...
    - 3.577e-3 .* T_K .* log(aSiO2_liq);

P_kbar = P_GPa .* 10;

% Diagnostic applicability flags. Study-wide ranges are warning envelopes,
% not exact Equation (8)-specific calibration limits.
isWithinStudyTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 950 & T_degreeC <= 1540;
isWithinStudyPRange = ...
    isfinite(P_kbar) & P_kbar >= 0.001 & P_kbar <= 23;
isWithinStudyCompositionRange = ...
    isfinite(SiO2_liq) & SiO2_liq >= 40.2 & SiO2_liq <= 66.16 & ...
    isfinite(MgO_liq) & MgO_liq >= 0.62 & MgO_liq <= 24.52 & ...
    isfinite(FeO_liq) & FeO_liq >= 2.80 & FeO_liq <= 30.2 & ...
    isfinite(Al2O3_liq) & Al2O3_liq >= 2.97 & Al2O3_liq <= 20.5 & ...
    isfinite(Na2O_liq + K2O_liq) & ...
    (Na2O_liq + K2O_liq) >= 0.19 & ...
    (Na2O_liq + K2O_liq) <= 6.77;

isWithinEquationDomain = ...
    isfinite(cationTotal_liq) & cationTotal_liq > 0 & ...
    isfinite(aSiO2_liq) & aSiO2_liq > 0;

% Pack outputs using pre-sized vectors of equal height.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.SiO2_liq = SiO2_liq;
row.TiO2_liq = TiO2_liq;
row.Al2O3_liq = Al2O3_liq;
row.FeO_liq = FeO_liq;
row.FeO_source = FeO_source;
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
row.XMnO_liq = XMnO_liq;
row.XMgO_liq = XMgO_liq;
row.XNaO0_5_liq = XNaO0_5_liq;
row.XKO0_5_liq = XKO0_5_liq;
row.aSiO2_liq = aSiO2_liq;

% General launcher-compatible pressure names.
row.P_GPa = P_GPa;
row.P_kbar = P_kbar;

% Backward-compatible Equation (8)-specific aliases.
row.PEq8_GPa = P_GPa;
row.PEq8_kbar = P_kbar;

row.P_RMSE_calibration_GPa = repmat(0.16, nT, 1);
row.P_RMSE_test_GPa = repmat(0.18, nT, 1);
row.isWithinEquationDomain = isWithinEquationDomain;
row.isWithinStudyTRange = isWithinStudyTRange;
row.isWithinStudyPRange = isWithinStudyPRange;
row.isWithinStudyCompositionRange = isWithinStudyCompositionRange;

row.trace_hasMineralRow = true(nT, 1);
row.trace_mineralRowHeight = repmat(height(data_min), nT, 1);
row.trace_mineralUsedInEquation = false(nT, 1);
row.trace_FClExcludedFromCationTotal = true(nT, 1);

end

function printCompositionRangeWarnings(liquidComposition, selectedCode_min, selectedIdx_liq)
% printCompositionRangeWarnings
% Print non-stopping warnings for finite liquid values outside the study-wide
% compositional ranges on p. 2165. These are not exact Eq. (8)-specific limits.

componentNames = {'SiO2', 'MgO', 'FeO', 'Al2O3', 'TotalAlkalis'};
componentValues = [ ...
    liquidComposition.SiO2, ...
    liquidComposition.MgO, ...
    liquidComposition.FeO, ...
    liquidComposition.Al2O3, ...
    liquidComposition.Na2O + liquidComposition.K2O];
minimumValues = [40.2, 0.62, 2.80, 2.97, 0.19];
maximumValues = [66.16, 24.52, 30.2, 20.5, 6.77];

for i = 1:numel(componentNames)
    value = componentValues(i);
    if isfinite(value) && ...
            (value < minimumValues(i) || value > maximumValues(i))
        fprintf(2, ...
            ['WARNING: Liquid %s = %.4g wt%% is outside the study-wide ' ...
             'experimental range %.4g-%.4g wt%% reported by Lessel and ' ...
             'Putirka (2015; p. 2165) for %s & Liquid Row %d. The exact ' ...
             'Equation (8)-specific range is not tabulated separately.\n'], ...
            componentNames{i}, ...
            value, ...
            minimumValues(i), ...
            maximumValues(i), ...
            char(string(selectedCode_min)), ...
            selectedIdx_liq);
    end
end

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Repeat optional liquid identifiers to match the number of output rows.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    rawIndex = data_liq.('Index');
    if isnumeric(rawIndex) || islogical(rawIndex)
        row.liq_Index = repmat(double(rawIndex(1)), nRows, 1);
    else
        row.liq_Index = repmat(string(rawIndex(1)), nRows, 1);
    end
end

if any(strcmp(variableNames, 'Experiment'))
    rawExperiment = data_liq.('Experiment');
    row.liq_Experiment = repmat(string(rawExperiment(1)), nRows, 1);
end

if any(strcmp(variableNames, 'Citation'))
    rawCitation = data_liq.('Citation');
    row.liq_Citation = repmat(string(rawCitation(1)), nRows, 1);
end

end

function value = getLiquidOxide(data_liq, oxide)
% getLiquidOxide
% Retrieve one liquid oxide value. Missing, empty, or non-numeric entries are
% represented by NaN and are never replaced by zero.

columnName = findOxideColumn(data_liq.Properties.VariableNames, oxide);
if isempty(columnName)
    value = NaN;
    return;
end

value = toScalarDouble(data_liq.(columnName));

end

function columnName = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match oxide column names while ignoring spaces, underscores, and hyphens.

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

function value = toScalarDouble(rawValue)
% toScalarDouble
% Convert the first entry of a one-row table variable to double. Missing,
% empty, or non-numeric values return NaN.

value = NaN;

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
    parsedValue = str2double(rawValue(1));
    if ~isnan(parsedValue)
        value = parsedValue;
    end
    return;
end

if ischar(rawValue)
    parsedValue = str2double(string(rawValue));
    if ~isnan(parsedValue)
        value = parsedValue;
    end
    return;
end

if iscategorical(rawValue)
    if isundefined(rawValue(1))
        return;
    end
    parsedValue = str2double(string(rawValue(1)));
    if ~isnan(parsedValue)
        value = parsedValue;
    end
    return;
end

if iscell(rawValue)
    if isempty(rawValue{1})
        return;
    end
    cellValue = rawValue{1};
    if isnumeric(cellValue) || islogical(cellValue)
        value = double(cellValue(1));
        return;
    end
    if isstring(cellValue) || ischar(cellValue) || iscategorical(cellValue)
        parsedValue = str2double(string(cellValue));
        if ~isnan(parsedValue)
            value = parsedValue;
        end
    end
end

end
