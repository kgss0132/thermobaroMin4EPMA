function results = LesselPutirka2015OpxLiqbaro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Pyroxene/LesselPutirka2015OpxLiqbaro.m
% Compatibility target: MATLAB R2024b
%
% Orthopyroxene-Liquid barometer
% Lessel, J. and Putirka, K. (2015)
% American Mineralogist, 100, 2163-2171
% DOI: https://doi.org/10.2138/am-2015-4732
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis selected
% from rawdata_struct.Opx with one Liquid analysis read by
% liquid.readLiquidExcel and calculates pressure using Equation (3) of
% Lessel and Putirka (2015; pp. 2165-2166).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Orthopyroxene-Liquid pair, one
% output row is returned for every input temperature value.
%
% IMPORTANT: Equation (3) contains no temperature term. Therefore, for one
% selected Orthopyroxene-Liquid pair, the same calculated pressure is
% repeated for every input temperature. T_degreeC and T_K are retained in
% the output so that the common fixed-T and range-T launchers can use the
% same function interface.
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
% equation-specific calibration limits for Equation (3). They are used in
% this implementation only as non-stopping warning envelopes. Because
% Equation (3) is temperature-independent, the temperature warning is a
% contextual extrapolation warning and does not alter the calculated P.
%
% Equation (3) was calibrated using Orthopyroxene-Liquid equilibria from
% experiments on martian bulk compositions. It reproduces the calibration
% pressures with R2 = 0.89 and RMSE = +/-0.20 GPa (n = 75), and predicts
% the independent test data with R2 = 0.91 and RMSE = +/-0.21 GPa
% (n = 14; pp. 2164 and 2166, Table 2 and Fig. 4).
%
% Orthopyroxene and Liquid must represent an equilibrium pair. Lessel and
% Putirka (2015) state that the global mean KD(Fe-Mg)Opx-Liq is sensitive to
% liquid composition and should not be used by itself as an equilibrium
% test. They recommend their Equation (12) on p. 2167:
%
%   KD(Fe-Mg)Opx-Liq = 0.32 - 0.05*K2O_Liq(wt%)
%
% with SEE approximately +/-0.03. Their KD calculations use a liquid
% Fe3+/Fe2+ correction based on Kress and Carmichael (1991, Eq. 7). This
% implementation outputs the observed and predicted KD values as diagnostic
% quantities, but users must consider redox treatment and petrographic
% evidence when deciding whether an Opx-Liquid pair is in equilibrium.
%
% Volatile-bearing experiments were included in the overall calibration and
% test database and were treated in the same manner as other experiments
% (p. 2165). However, no separate H2O, F, or Cl correction is included in
% Equation (3). In this implementation, Liquid F and Cl are excluded from
% cationTotal_liq and are excluded from the NaN-input warning check.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 950-1540 degreeC,
%   2) finite calculated pressure is outside approximately 0.001-23 kbar,
%   3) finite liquid composition is outside the broad combined-database
%      ranges listed above,
%   4) the observed Opx-Liquid KD differs from Equation (12) by more than
%      its reported SEE of +/-0.03,
%   5) a calculation input contains NaN, or
%   6) a calculated pressure is NaN or Inf.
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if it makes a ratio undefined, the
% resulting NaN or Inf is retained and reported.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%
% The FIRST column of the Opx table is treated as an identifier ("data
% code") displayed in the selection dialog. Orthopyroxene oxide columns are
% read using common names such as SiO2, TiO2, Al2O3, FeO or FeOt, MnO, MgO,
% CaO, Na2O, K2O, and Cr2O3. Columns with a "_value" suffix are also
% recognized.
%
% Liquid compositions are obtained using liquid.readLiquidExcel. The
% Equation (3) calculation uses liquid cation fractions normalized to the
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
% Lessel and Putirka (2015), Equation (3), pp. 2165-2166:
%
%   P(GPa) =
%     -5.050
%     + 35.05*XTiO2_Liq
%     + 6.458*XAlO1.5_Liq
%     - 10.67*XCaO_Liq
%     + 1.438*Mg#_Liq
%     - 20.74*Ti_Opx
%     - 6.188*Al(total)_Opx
%     + 0.01915*D(MgO)_Opx/Liq
%     + 1.111*[Al(total)_Opx/XAlO1.5_Liq
%              + Na_Opx/XNaO0.5_Liq
%              + Si_Opx/XSiO2_Liq]
%
% where:
%
%   Mg#_Liq        = XMgO_Liq/(XMgO_Liq + XFeO_Liq)
%   D(MgO)_Opx/Liq = Mg_Opx/XMgO_Liq
%
% Orthopyroxene cations are calculated on a six-oxygen basis, and Liquid
% components are cation fractions, following Putirka (2008), as specified
% by Lessel and Putirka (2015; pp. 2164-2166).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = LesselPutirka2015OpxLiqbaro(rawdata_struct, T_degreeC)
%   results = LesselPutirka2015OpxLiqbaro(..., 'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing an Opx table
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Name-value option:
%   'LiquidRow'    : positive integer scalar selecting a row in the loaded
%                    Liquid table. Default [] uses row 1.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Orthopyroxene-Liquid pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('LesselPutirka2015OpxLiqbaro requires (rawdata_struct, T_degreeC).');
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

%% 1) Retrieve Orthopyroxene and Liquid datasets
disp('=== Step 1: Preparing Orthopyroxene and Liquid datasets ===');

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end

dataset_opx = rawdata_struct.Opx;
if height(dataset_opx) < 1
    error('rawdata_struct.Opx is empty.');
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
disp('=== Preparing Orthopyroxene and Liquid datasets has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a preallocated cell buffer and
% concatenate once after the interactive loop. The buffer size does not
% change on each loop iteration.
disp('=== Step 2: Preparing output container ===');

maxResultBlocks = max(1024, 16 * height(dataset_opx));
resultBlocks = cell(maxResultBlocks, 1);
nResultBlocks = 0;

% Broad ranges of the combined experimental database. These are used only
% as warning envelopes and are not strict Equation (3)-specific limits.
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
disp('=== Step 3: Selecting a data code from the list (Orthopyroxene) ===');

while true
    % ----- Orthopyroxene selection -----
    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Orthopyroxene data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Orthopyroxene selected: ' char(string(selectedCode_opx))]);

    % ----- Calculation -----
    disp('=== Step 4: Checking calculation inputs ===');

    selectedData_opx = dataset_opx(selectedIdx_opx, :);

    % Check NaN only in variables that enter cation normalization,
    % component calculation, temperature metadata, or Equation (3). Liquid
    % F and Cl are intentionally excluded from this check.
    nanInputNames = findNaNInputs( ...
        selectedData_opx, selectedData_liq, T_degreeC);

    % Reject Inf and finite negative values in calculation inputs. NaN and
    % zero are retained. Liquid F and Cl are not calculation inputs here.
    validateNonNegativeInputs(selectedData_opx, selectedData_liq);

    disp('=== Step 5: Calculating the pressure ===');
    row = calcPressure( ...
        selectedData_opx, selectedData_liq, T_degreeC, MWinfo);

    % Repeat identifiers for all temperatures in the current calculation.
    nRows = height(row);
    row.dataCode_opx = repmat(string(selectedCode_opx), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_opx', 'dataRow_liq'}, 'Before', 1);

    % Store one block per selected pair without resizing the result buffer.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > maxResultBlocks
        error(['Maximum number of interactive selections (%d) exceeded. ' ...
               'Restart the function to continue.'], maxResultBlocks);
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_opx)) ' & Liquid row ' ...
            num2str(selectedIdx_liq) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_opx)) ' & Liquid row ' ...
            num2str(selectedIdx_liq) ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Print the model-specific limitation once per function call.
    if ~modelCautionIssued
        fprintf(2, ...
            ['CAUTION: Lessel and Putirka (2015) Equation (3) was calibrated ' ...
             'specifically for martian igneous compositions (pp. 2163-2166). ' ...
             'Application to terrestrial or compositionally dissimilar systems is ' ...
             'extrapolation. Orthopyroxene-Liquid equilibrium must be checked ' ...
             'independently (Eq. 12, p. 2167). Equation (3) contains no ' ...
             'temperature term, so pressure is repeated unchanged across the ' ...
             'input temperature vector.\n']);
        modelCautionIssued = true;
    end

    % Input temperature is common to all selected pairs, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideExperimental) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the 950-1540 degreeC range ' ...
             'of the combined experimental database of Lessel and Putirka ' ...
             '(2015; pp. 2164-2165). This is not a strict Equation (3)-specific ' ...
             'range, and Equation (3) itself is temperature-independent. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
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
             'a strict Equation (3)-specific range. %d of %d finite pressure ' ...
             'point(s) are outside the range; calculated finite range = ' ...
             '%.4g-%.4g kbar for %s & Liquid row %d.\n'], ...
            sum(pressureOutsideExperimental), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_opx)), ...
            selectedIdx_liq);
    end

    % Warn when finite liquid compositions lie outside the broad combined
    % experimental database ranges listed in the paper.
    printLiquidCompositionWarnings( ...
        row, selectedCode_opx, selectedIdx_liq);

    % Compare observed Fe-Mg exchange with the composition-dependent
    % equilibrium relation of Equation (12). This is diagnostic only.
    printKDEquilibriumWarning(row, selectedCode_opx, selectedIdx_liq);

    % List the exact calculation-input names containing NaN. For vector
    % temperature input, the NaN element indices are included. F and Cl are
    % excluded from this warning by design.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s & ' ...
             'Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated where they enter ' ...
             'Equation (3); the calculated pressure may remain NaN. Because ' ...
             'Equation (3) is temperature-independent, NaN in T_degreeC does not ' ...
             'by itself force P to NaN. Liquid F and Cl were excluded from this ' ...
             'check and from cationTotal_liq.\n'], ...
            char(string(selectedCode_opx)), ...
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
            char(string(selectedCode_opx)), ...
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
            char(string(selectedCode_opx)), ...
            selectedIdx_liq, ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    % Report non-finite component ratios explicitly. These values are not
    % replaced or clipped.
    invalidRatio = ~isfinite(row.ratioTerm_Eq3) | ...
        ~isfinite(row.MgNumber_liq) | ~isfinite(row.DMgO_OpxLiq);
    if any(invalidRatio)
        fprintf(2, ...
            ['WARNING: One or more Equation (3) ratio terms are non-finite for ' ...
             '%s & Liquid row %d (%d of %d points). NaN or Inf values were ' ...
             'retained and propagated.\n'], ...
            char(string(selectedCode_opx)), ...
            selectedIdx_liq, ...
            sum(invalidRatio), ...
            numel(row.ratioTerm_Eq3));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Orthopyroxene selection (same Liquid row)?', ...
        'LesselPutirka2015OpxLiqbaro', ...
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
function nanInputNames = findNaNInputs(data_opx, data_liq, T_degreeC)
% findNaNInputs
% Return names of Equation (3) calculation inputs containing NaN. NaN values
% do not cause an error and are not replaced by zero. Liquid F and Cl are
% intentionally excluded.

maxNames = 64;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

% Opx oxides used for six-oxygen normalization and Equation (3).
opxOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO','Na2O','K2O','Cr2O3'};
for i = 1:numel(opxOxides)
    oxide = opxOxides{i};
    columnName = findOxideColumn(data_opx.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDouble(data_opx.(columnName));
        if isnan(value)
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = ...
                "Opx." + string(columnName);
        end
    end
end

% Use FeO when that column exists; otherwise use FeOt.
feColumnName = findOxideColumn(data_opx.Properties.VariableNames, 'FeO');
if isempty(feColumnName)
    feColumnName = findOxideColumn(data_opx.Properties.VariableNames, 'FeOt');
end
if ~isempty(feColumnName)
    value = toScalarDouble(data_opx.(feColumnName));
    if isnan(value)
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Opx." + string(feColumnName);
    end
end

% Liquid oxides included in cationTotal_liq. F and Cl are excluded.
liqOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO','Na2O','K2O', ...
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

% Use FeO when present, otherwise FeOt, for the Liquid Fe term.
liqFeColumnName = findOxideColumn(data_liq.Properties.VariableNames, 'FeO');
if isempty(liqFeColumnName)
    liqFeColumnName = findOxideColumn(data_liq.Properties.VariableNames, 'FeOt');
end
if ~isempty(liqFeColumnName)
    value = toScalarDouble(data_liq.(liqFeColumnName));
    if isnan(value)
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = ...
            "Liquid." + string(liqFeColumnName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_opx, data_liq)
% validateNonNegativeInputs
% Reject Inf and finite negative values in raw oxide variables used by the
% calculation. Zero and NaN are intentionally allowed and retained. Liquid
% F and Cl are not included because they are excluded from cationTotal_liq.

maxNames = 64;
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

opxOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO','Na2O','K2O','Cr2O3'};
for i = 1:numel(opxOxides)
    oxide = opxOxides{i};
    columnName = findOxideColumn(data_opx.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDouble(data_opx.(columnName));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = ...
                "Opx." + string(columnName);
        end
    end
end

feColumnName = findOxideColumn(data_opx.Properties.VariableNames, 'FeO');
if isempty(feColumnName)
    feColumnName = findOxideColumn(data_opx.Properties.VariableNames, 'FeOt');
end
if ~isempty(feColumnName)
    value = toScalarDouble(data_opx.(feColumnName));
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Opx." + string(feColumnName);
    end
end

liqOxides = {'SiO2','TiO2','Al2O3','MnO','MgO','CaO','Na2O','K2O', ...
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

liqFeColumnName = findOxideColumn(data_liq.Properties.VariableNames, 'FeO');
if isempty(liqFeColumnName)
    liqFeColumnName = findOxideColumn(data_liq.Properties.VariableNames, 'FeOt');
end
if ~isempty(liqFeColumnName)
    value = toScalarDouble(data_liq.(liqFeColumnName));
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Liquid." + string(liqFeColumnName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['LesselPutirka2015OpxLiqbaro: calculation inputs must be ' ...
           'non-negative. NaN is allowed, but Inf and finite negative ' ...
           'value(s) are prohibited. Invalid input(s): ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_opx, data_liq, T_degreeC, MWinfo)
% calcPressure
% Compute pressure for one Orthopyroxene row and one Liquid row. Equation
% (3) is temperature-independent, so the scalar pressure is repeated for
% each input temperature. NaN values propagate through the calculation.
%
% Inputs:
%   data_opx    : 1-row Orthopyroxene table
%   data_liq    : 1-row Liquid table
%   T_degreeC   : scalar or vector temperature in degreeC
%   MWinfo      : molar-weight and cation-number structure
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Extract and normalize one-row Opx and Liquid data. Existing NaN values are
% preserved. Missing optional columns use zero only for backward
% compatibility.
opx = prepareOpxRow(data_opx, MWinfo);
liq = extractLiquidRow(data_liq, MWinfo);

% Equation (3) composition terms. No finite-value guard is used: NaN and Inf
% remain available for diagnostic output.
MgNumberScalar = liq.XMgO ./ (liq.XMgO + liq.XFeO);
DMgOScalar = opx.XMg ./ liq.XMgO;
ratioTermScalar = ...
    (opx.XAl ./ liq.XAlO1_5) + ...
    (opx.XNa ./ liq.XNaO0_5) + ...
    (opx.XSi ./ liq.XSiO2);

% Lessel and Putirka (2015) Equation (3), pp. 2165-2166.
P_GPa_scalar = ...
    -5.050 ...
    + 35.05 .* liq.XTiO2 ...
    + 6.458 .* liq.XAlO1_5 ...
    - 10.67 .* liq.XCaO ...
    + 1.438 .* MgNumberScalar ...
    - 20.74 .* opx.XTi ...
    - 6.188 .* opx.XAl ...
    + 0.01915 .* DMgOScalar ...
    + 1.111 .* ratioTermScalar;

P_kbar_scalar = P_GPa_scalar .* 10;

% Fe-Mg exchange diagnostics from Equation (12), p. 2167. The observed KD
% uses the Fe and Mg terms available to this calculation and does not add an
% independent Kress-Carmichael redox correction.
KD_observed_scalar = ...
    (opx.XFe ./ opx.XMg) ./ (liq.XFeO ./ liq.XMgO);
KD_predicted_scalar = 0.32 - 0.05 .* liq.K2O;
KD_deviation_scalar = KD_observed_scalar - KD_predicted_scalar;

% Expand all selected-pair scalars to the temperature-vector length.
XSi_opx = repmat(opx.XSi, nT, 1);
XTi_opx = repmat(opx.XTi, nT, 1);
XAl_opx = repmat(opx.XAl, nT, 1);
XFe_opx = repmat(opx.XFe, nT, 1);
XMn_opx = repmat(opx.XMn, nT, 1);
XMg_opx = repmat(opx.XMg, nT, 1);
XCa_opx = repmat(opx.XCa, nT, 1);
XNa_opx = repmat(opx.XNa, nT, 1);
XK_opx = repmat(opx.XK, nT, 1);
XCr_opx = repmat(opx.XCr, nT, 1);
cationSum_opx = repmat(opx.cationSum, nT, 1);
XAlIV_opx = repmat(opx.XAlIV, nT, 1);
XAlVI_opx = repmat(opx.XAlVI, nT, 1);

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

MgNumber_liq = repmat(MgNumberScalar, nT, 1);
DMgO_OpxLiq = repmat(DMgOScalar, nT, 1);
ratioTerm_Eq3 = repmat(ratioTermScalar, nT, 1);
P_GPa = repmat(P_GPa_scalar, nT, 1);
P_kbar = repmat(P_kbar_scalar, nT, 1);

KD_FeMg_OpxLiq_observed = repmat(KD_observed_scalar, nT, 1);
KD_FeMg_OpxLiq_predicted = repmat(KD_predicted_scalar, nT, 1);
KD_FeMg_OpxLiq_deviation = repmat(KD_deviation_scalar, nT, 1);

% Applicability and diagnostic flags. Combined-database ranges are warning
% envelopes only, not strict Equation (3)-specific calibration limits.
isWithinEquationDomain = ...
    isfinite(cationTotal_liq) & cationTotal_liq > 0 & ...
    isfinite(XMgO_liq) & XMgO_liq > 0 & ...
    isfinite(XAlO1_5_liq) & XAlO1_5_liq > 0 & ...
    isfinite(XNaO0_5_liq) & XNaO0_5_liq > 0 & ...
    isfinite(XSiO2_liq) & XSiO2_liq > 0 & ...
    isfinite(XMgO_liq + XFeO_liq) & (XMgO_liq + XFeO_liq) > 0;

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

isWithinKDEquilibriumSEE = ...
    isfinite(KD_FeMg_OpxLiq_deviation) & ...
    abs(KD_FeMg_OpxLiq_deviation) <= 0.03;

% Pack outputs using pre-sized vectors of equal height.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.XSi_opx = XSi_opx;
row.XTi_opx = XTi_opx;
row.XAl_opx = XAl_opx;
row.XFe_opx = XFe_opx;
row.XMn_opx = XMn_opx;
row.XMg_opx = XMg_opx;
row.XCa_opx = XCa_opx;
row.XNa_opx = XNa_opx;
row.XK_opx = XK_opx;
row.XCr_opx = XCr_opx;
row.cationSum_opx = cationSum_opx;
row.XAlIV_opx = XAlIV_opx;
row.XAlVI_opx = XAlVI_opx;

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

row.MgNumber_liq = MgNumber_liq;
row.DMgO_OpxLiq = DMgO_OpxLiq;
row.ratioTerm_Eq3 = ratioTerm_Eq3;

row.KD_FeMg_OpxLiq_observed = KD_FeMg_OpxLiq_observed;
row.KD_FeMg_OpxLiq_predicted = KD_FeMg_OpxLiq_predicted;
row.KD_FeMg_OpxLiq_deviation = KD_FeMg_OpxLiq_deviation;
row.KD_FeMg_OpxLiq_SEE = repmat(0.03, nT, 1);

% Primary pressure variables used by the common barometer launchers.
row.P_GPa = P_GPa;
row.P_kbar = P_kbar;

% Backward-compatible aliases retained from the original implementation.
row.PEq3_GPa = P_GPa;
row.PEq3_kbar = P_kbar;

row.P_RMSE_calibration_GPa = repmat(0.20, nT, 1);
row.P_RMSE_test_GPa = repmat(0.21, nT, 1);
row.P_RMSE_calibration_kbar = repmat(2.0, nT, 1);
row.P_RMSE_test_kbar = repmat(2.1, nT, 1);

row.isWithinEquationDomain = isWithinEquationDomain;
row.isWithinExperimentalTRange = isWithinExperimentalTRange;
row.isWithinExperimentalPRange = isWithinExperimentalPRange;
row.isWithinBroadLiquidRange = isWithinBroadLiquidRange;
row.isWithinKDEquilibriumSEE = isWithinKDEquilibriumSEE;

% Backward-compatible/general aliases.
row.isRecommended_T_range = isWithinExperimentalTRange;
row.isApplicable_composition = isWithinBroadLiquidRange;

end

function opx = prepareOpxRow(data_opx, MWinfo)
% prepareOpxRow
% Extract one-row Orthopyroxene oxide data and calculate six-oxygen cations.
% Existing NaN values are retained.

if height(data_opx) ~= 1
    error('Orthopyroxene input must be a 1-row table.');
end

SiO2 = getMineralOxRequired(data_opx, 'SiO2', 'Orthopyroxene');
TiO2 = getMineralOxOptional(data_opx, 'TiO2', 0);
Al2O3 = getMineralOxRequired(data_opx, 'Al2O3', 'Orthopyroxene');
MnO = getMineralOxOptional(data_opx, 'MnO', 0);
MgO = getMineralOxRequired(data_opx, 'MgO', 'Orthopyroxene');
CaO = getMineralOxOptional(data_opx, 'CaO', 0);
Na2O = getMineralOxOptional(data_opx, 'Na2O', 0);
K2O = getMineralOxOptional(data_opx, 'K2O', 0);
Cr2O3 = getMineralOxOptional(data_opx, 'Cr2O3', 0);

% Prefer FeO when the column exists. If FeO exists but is NaN, retain that
% NaN rather than silently replacing it with FeOt. Use FeOt only when the
% FeO column is absent.
feOColumn = findOxideColumn(data_opx.Properties.VariableNames, 'FeO');
if ~isempty(feOColumn)
    FeO = toScalarDouble(data_opx.(feOColumn));
else
    feOtColumn = findOxideColumn(data_opx.Properties.VariableNames, 'FeOt');
    if isempty(feOtColumn)
        error('Orthopyroxene table must contain FeO or FeOt.');
    end
    FeO = toScalarDouble(data_opx.(feOtColumn));
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

% Prefer FeO when present. If FeO is absent, use FeOt as the FeO-equivalent
% input. A present NaN remains NaN.
liqFeOColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeO');
if ~isempty(liqFeOColumn)
    liq.FeO = toScalarDouble(data_liq.(liqFeOColumn));
else
    liqFeOtColumn = findOxideColumn(data_liq.Properties.VariableNames, 'FeOt');
    if isempty(liqFeOtColumn)
        liq.FeO = 0;
    else
        liq.FeO = toScalarDouble(data_liq.(liqFeOtColumn));
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

function printLiquidCompositionWarnings(row, selectedCode_opx, selectedIdx_liq)
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
         '(3)-specific limits.\n'], ...
        char(string(selectedCode_opx)), ...
        selectedIdx_liq, ...
        char(strjoin(messages, '; ')));
end

end

function printKDEquilibriumWarning(row, selectedCode_opx, selectedIdx_liq)
% printKDEquilibriumWarning
% Compare the observed Fe-Mg exchange coefficient with Equation (12) on
% p. 2167. The +/-0.03 SEE is used as a diagnostic warning envelope only.

KD_observed = row.KD_FeMg_OpxLiq_observed(1);
KD_predicted = row.KD_FeMg_OpxLiq_predicted(1);
KD_deviation = row.KD_FeMg_OpxLiq_deviation(1);

if isfinite(KD_deviation) && abs(KD_deviation) > 0.03
    fprintf(2, ...
        ['WARNING: Opx-Liquid Fe-Mg exchange is outside the +/-0.03 SEE of ' ...
         'Lessel and Putirka (2015) Equation (12), p. 2167, for %s & Liquid ' ...
         'row %d: KD_observed = %.4g, KD_predicted = %.4g, difference = %.4g. ' ...
         'This is a diagnostic equilibrium check, not an automatic rejection ' ...
         'criterion, and liquid Fe3+/Fe2+ treatment must also be considered.\n'], ...
        char(string(selectedCode_opx)), ...
        selectedIdx_liq, ...
        KD_observed, ...
        KD_predicted, ...
        KD_deviation);
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
