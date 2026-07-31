function results = Molina2015Liq(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Amphibole/Molina2015Liq.m
% Tested with MATLAB R2024b
%
% Liquid-only thermometer
% Molina, J.F., Moreno, J.A., Castro, A., Rodriguez, C., Fershtater, G.B. (2015)
% Lithos, 232, 286-305
% DOI: https://doi.org/10.1016/j.lithos.2015.06.027
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function selects one mineral analysis for traceability and one liquid
% analysis, then calculates temperature using the Mg-in-liquid thermometer
% of Molina et al. (2015). The selected mineral composition is not used in
% the thermometer equation.
%
% The function accepts either a scalar or vector P_kbar input. One output
% row is returned for each pressure value so that the function can be called
% from both startThermoCalc_fixedP and startThermoCalc_rangeP. Pressure is
% retained in the output for plotting and traceability, but it is not used
% in the published liquid-only thermometer equation.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another mineral-liquid selection and stores all
% result blocks in a single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Molina et al. (2015) calibrated the Mg-in-liquid thermometer using liquids
% from experimental amphibole-glass assemblages. The principal calibration
% and application constraints are:
%
%   Temperature : 800-1100 degreeC
%                 The selected experimental database listed in Table 3 spans
%                 830-1100 degreeC, whereas the stated application range for
%                 the thermometers is 800-1100 degreeC.
%
%   Pressure    : 0.8-20 kbar experimental database range
%                 Pressure was not statistically significant in the fitted
%                 thermometer and therefore does not occur in the equation.
%                 Use outside 0.8-20 kbar nevertheless represents pressure
%                 extrapolation beyond the experimental database.
%
%   Rock type   : alkaline and subalkaline igneous systems
%
%   Phase context:
%                 The calibration liquids coexisted with calcic amphibole
%                 having CaM4/(CaM4 + NaM4) > 0.75. Although this liquid-only
%                 equation does not require amphibole composition as an
%                 input, it was calibrated from amphibole-bearing,
%                 multiply-saturated experimental liquids.
%
%   Precision   : approximately +/-37 degreeC for the calibration dataset
%                 and +/-42 degreeC for the independent test dataset (1 s).
%
%   Selected liquid-composition ranges (Table 3, p. 289):
%     XMg_liq = 0.0009773-0.1411817
%     XCa_liq = 0.0142063-0.1465581
%     XAl_liq = 0.1465827-0.2305199
%     ln[XCa_liq/(XCa_liq + XAl_liq)] = -2.523202 to -0.7690551
%
% The experimental database and data-selection procedure are described on
% pp. 287-289; Table 3 on p. 289 reports the selected amphibole-glass P-T and
% compositional ranges. Glass compositions and the alkaline/subalkaline
% coverage are described on pp. 293-294. The liquid-only thermometer and its
% statistical calibration are presented on pp. 296-298. Application to
% natural volcanic and plutonic rocks is discussed on pp. 303-304, and the
% principal restrictions are summarized in the Conclusions on p. 304.
%
% Natural applications require a liquid composition representative of the
% melt of interest. Matrix glass or glass-inclusion compositions may be used
% where appropriate. Whole-rock compositions should not be treated as liquid
% compositions without considering crystal accumulation, fractionation,
% mixing, alteration, and whether the inferred melt was amphibole-saturated.
% Comparison with an independent thermometer is recommended when possible.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 0.8-20 kbar,
%   2) a finite calculated temperature is outside 800-1100 degreeC,
%   3) a calculation input contains NaN, or
%   4) the calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain at least one of the following mineral tables:
%   rawdata_struct.Amphibole
%   rawdata_struct.Cpx
%   rawdata_struct.Opx
%   rawdata_struct.Olivine
%
% The mineral table is used only for interactive row selection and output
% traceability. Its FIRST column is treated as the displayed data code.
%
% Liquid data are read using liquid.readLiquidExcel(). Oxide columns are
% matched without regard to spaces, underscores, or hyphens, and either the
% oxide name or the oxide name followed by "value" is accepted.
%
% Liquid normalization uses the following cations:
%   Si, Ti, Al, Cr, Ni, Fe, Mn, Mg, Ca, Na, K, P
%
% Corresponding oxide inputs used in the calculation are:
%   SiO2, TiO2, Al2O3, FeO (or FeOt if FeO is absent), Fe2O3, MnO,
%   MgO, CaO, Na2O, K2O, Cr2O3, NiO, and P2O5.
%
% A missing optional oxide column is treated as zero, preserving the behavior
% of the original implementation. In contrast, a column that is present but
% contains NaN remains NaN; it is never replaced by zero. NaN therefore
% propagates through SumCat_liq, cation fractions, logarithms, and temperature.
%
% All finite liquid-composition values read by this function must be >= 0.
% Negative finite values stop the calculation. NaN values are retained and
% reported by non-stopping fprintf warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Cation mole fractions are calculated on an anhydrous basis:
%
%   SumCat_liq = Si + Ti + Al + Cr + Ni + Fe + Mn + Mg
%                + Ca + Na + K + P
%
%   XN_liq = N_liq / SumCat_liq
%
% The liquid-only thermometer is:
%
%   T(degreeC) = 107 * ln(XMg_liq)
%              - 108 * ln[XCa_liq / (XCa_liq + XAl_liq)]
%              + 1184
%
% Natural logarithms are used. Pressure is not used in this equation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Molina2015Liq(rawdata_struct, P_kbar)
%   results = Molina2015Liq(rawdata_struct, P_kbar, 'LiquidRow', n)
%
% Inputs:
%   rawdata_struct : struct containing at least one supported mineral table
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Name-value option:
%   'LiquidRow'    : positive integer scalar or [] (default [])
%                    If empty, row 1 of the selected liquid dataset is used.
%
% Output:
%   results : table containing one row per pressure value for every selected
%             mineral-liquid calculation. Both Molina-specific temperature
%             names and standard T_degreeC/T_K aliases are included.
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or an
% invalid pressure vector.
if nargin < 2
    error('Molina2015Liq requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve mineral dataset and parse options
% The mineral dataset is retained only for user selection and traceability.
disp('=== Step 1: Preparing mineral and liquid datasets ===');

[mineralField, dataset_min] = localPickMineralDataset(rawdata_struct);

ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == floor(x)));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

% Molecular weights and cation numbers are supplied by the project helper.
MWinfo = liquid.getMolarWeights();

% Load the liquid dataset using the standard project interface.
[liqAll, metaLiq] = liquid.readLiquidExcel();
if ~istable(liqAll) || isempty(liqAll)
    error('Selected liquid dataset must be a non-empty table.');
end

disp('=== Preparing mineral and liquid datasets has been finished ===');

%% 2) Initialize output container and calibration limits
% Each calculation is stored as one table block. Repeated concatenation of
% the full output table inside the interactive loop is avoided.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Calibration/application limits reported by Molina et al. (2015).
calibrationT_min_degC = 800;
calibrationT_max_degC = 1100;
calibrationP_min_kbar = 0.8;
calibrationP_max_kbar = 20;

% Pressure is common to all calculations in this function call, so the
% pressure warning is printed only once.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop and calculation
% The same liquid dataset remains loaded while the user selects additional
% mineral rows for traceability.
disp(['=== Step 3: Selecting a data code from the list (' mineralField ') ===']);

while true
    % ----- Mineral selection for traceability -----
    dataCodes_min = dataset_min{:, 1};
    listStrings_min = cellstr(string(dataCodes_min));

    [selectedIdx_min, ok] = listdlg( ...
        'PromptString', ['Please select the ' mineralField ' data you would like to use:'], ...
        'SelectionMode', 'single', ...
        'ListString', listStrings_min, ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_min)
        disp('Selection canceled');
        break;
    end

    selectedCodeString = string(listStrings_min{selectedIdx_min});
    disp([mineralField ' selected: ' char(selectedCodeString)]);

    % ----- Liquid-row selection -----
    disp('=== Step 4: Selecting the liquid row ===');

    if isempty(liquidRowOpt)
        idxLiq = 1;
        if height(liqAll) > 1
            fprintf(2, ...
                ['WARNING: The selected liquid dataset contains %d rows. ' ...
                 'LiquidRow was not specified, so row 1 is being used.\n'], ...
                height(liqAll));
        end
    else
        idxLiq = liquidRowOpt;
        if idxLiq > height(liqAll)
            error(['Requested LiquidRow (%d) exceeds the number of rows in ' ...
                   'the selected liquid dataset (%d).'], idxLiq, height(liqAll));
        end
    end

    disp(['Liquid selected: Row ' num2str(idxLiq)]);

    % ----- Prepare selected rows and liquid inputs -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_min = dataset_min(selectedIdx_min, :);
    selectedData_liq = liqAll(idxLiq, :);

    % Extract oxide values once. Present NaN values remain NaN; missing
    % optional columns retain the original zero-default behavior.
    [liqInputs, sourceNames] = extractLiquidInputs(selectedData_liq);

    % Identify NaN inputs before calculation so their exact source columns
    % can be printed after the result without interrupting the calculation.
    nanInputNames = findNaNInputs(liqInputs, sourceNames);

    % Negative finite oxide values are physically invalid and are prohibited.
    % NaN is intentionally allowed and will propagate through the equation.
    validateNonnegativeInputs(liqInputs, sourceNames);

    % The pressure vector is passed intact. calcTemp returns one row for every
    % pressure value even though pressure is not used by the thermometer.
    row = calcTemp(selectedData_min, liqInputs, sourceNames, P_kbar, MWinfo);
    nRows = height(row);

    % Store identifiers and liquid metadata for traceability.
    row.dataCode_mineral = repmat(selectedCodeString, nRows, 1);
    row.mineralType = repmat(string(mineralField), nRows, 1);
    row.dataRow_liq = repmat(idxLiq, nRows, 1);
    row = localAttachLiquidIDs(row, selectedData_liq);
    row = movevars(row, ...
        {'mineralType', 'dataCode_mineral', 'dataRow_liq'}, 'Before', 1);

    % Store this calculation as one buffered table block. The buffer grows
    % geometrically only when full, rather than changing on every iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % ----- Display calculated temperature -----
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if nRows == 1
        disp([char(selectedCodeString) ' & Liquid row ' num2str(idxLiq) ...
            ': ' num2str(row.T_Molina2015Liq_C) ' degreeC']);
    else
        disp([char(selectedCodeString) ' & Liquid row ' num2str(idxLiq) ...
            ': ' num2str(row.T_Molina2015Liq_C(1)) ' to ' ...
            num2str(row.T_Molina2015Liq_C(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside the 0.8-20 kbar
    % experimental database range. The equation itself is pressure-independent.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental database range ' ...
             'of Molina et al. (2015): 0.8-20 kbar. Pressure is not used in ' ...
             'the liquid-only equation, but %d of %d pressure point(s) are ' ...
             'outside the experimental range; input range = %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite temperature lies outside the stated 800-1100
    % degreeC application range. NaN and Inf are handled separately below.
    finiteTemperature = isfinite(row.T_Molina2015Liq_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_Molina2015Liq_C < calibrationT_min_degC | ...
         row.T_Molina2015Liq_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_Molina2015Liq_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the stated application ' ...
             'range of Molina et al. (2015): 800-1100 degreeC. %d of %d finite ' ...
             'temperature point(s) are outside the range; calculated finite ' ...
             'range = %.4g-%.4g degreeC for %s and Liquid row %d.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(selectedCodeString), ...
            idxLiq);
    end

    % Print the names of all calculation inputs that contained NaN. These
    % values were retained as NaN and were not replaced by zero.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the liquid thermometer input(s) for %s ' ...
             'and Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated through the calculation.\n'], ...
            char(selectedCodeString), ...
            idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report NaN/Inf output values caused by missing values, zero
    % logarithm arguments, zero cation sums, or other numerical conditions.
    invalidTemperature = ~isfinite(row.T_Molina2015Liq_C);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s and ' ...
             'Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(selectedCodeString), ...
            idxLiq, ...
            sum(invalidTemperature), ...
            numel(row.T_Molina2015Liq_C), ...
            sum(isnan(row.T_Molina2015Liq_C)), ...
            sum(isinf(row.T_Molina2015Liq_C)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another mineral selection.
    userAction = questdlg( ...
        'Continue with another data selection (same Liquid dataset)?', ...
        'Molina2015Liq', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered result blocks once after the selection loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

% Preserve liquid-file metadata used by downstream project functions.
results.Properties.UserData = struct('liquid', metaLiq);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function [mineralField, dataset_min] = localPickMineralDataset(rawdata_struct)
% localPickMineralDataset
% Return the first supported mineral table found in rawdata_struct. The
% mineral composition is not used in the liquid-only thermometer.

candidateFields = {'Amphibole', 'Cpx', 'Opx', 'Olivine'};
found = false;
mineralField = '';
dataset_min = table();

for i = 1:numel(candidateFields)
    fieldName = candidateFields{i};
    if isfield(rawdata_struct, fieldName) && istable(rawdata_struct.(fieldName))
        mineralField = fieldName;
        dataset_min = rawdata_struct.(fieldName);
        found = true;
        break;
    end
end

if ~found
    error(['rawdata_struct must contain at least one mineral table for ' ...
           'interactive selection: rawdata_struct.Amphibole, ' ...
           'rawdata_struct.Cpx, rawdata_struct.Opx, or ' ...
           'rawdata_struct.Olivine.']);
end
if isempty(dataset_min)
    error('The selected mineral table is empty.');
end

end

function [liqInputs, sourceNames] = extractLiquidInputs(data_liq)
% extractLiquidInputs
% Extract scalar liquid oxide values from one liquid-table row. If a column
% exists and contains NaN, NaN is preserved. If an optional column is absent,
% the supplied default value is used.

liqInputs = struct();
sourceNames = struct();

[liqInputs.SiO2,  sourceNames.SiO2]  = localGetLiqOxOptional(data_liq, 'SiO2',  0);
[liqInputs.TiO2,  sourceNames.TiO2]  = localGetLiqOxOptional(data_liq, 'TiO2',  0);
[liqInputs.Al2O3, sourceNames.Al2O3] = localGetLiqOxOptional(data_liq, 'Al2O3', 0);

% Prefer FeO. Use FeOt only when an FeO column is absent. A present FeO
% value of NaN remains NaN and is not replaced by FeOt or zero.
[liqInputs.FeO, sourceNames.FeO] = localGetLiqOxOptional(data_liq, 'FeO', NaN);
if strlength(sourceNames.FeO) == 0
    [liqInputs.FeO, sourceNames.FeO] = ...
        localGetLiqOxOptional(data_liq, 'FeOt', 0);
end

[liqInputs.MnO,   sourceNames.MnO]   = localGetLiqOxOptional(data_liq, 'MnO',   0);
[liqInputs.MgO,   sourceNames.MgO]   = localGetLiqOxOptional(data_liq, 'MgO',   0);
[liqInputs.CaO,   sourceNames.CaO]   = localGetLiqOxOptional(data_liq, 'CaO',   0);
[liqInputs.Na2O,  sourceNames.Na2O]  = localGetLiqOxOptional(data_liq, 'Na2O',  0);
[liqInputs.K2O,   sourceNames.K2O]   = localGetLiqOxOptional(data_liq, 'K2O',   0);
[liqInputs.V2O3,  sourceNames.V2O3]  = localGetLiqOxOptional(data_liq, 'V2O3',  0);
[liqInputs.Cr2O3, sourceNames.Cr2O3] = localGetLiqOxOptional(data_liq, 'Cr2O3', 0);
[liqInputs.NiO,   sourceNames.NiO]   = localGetLiqOxOptional(data_liq, 'NiO',   0);
[liqInputs.P2O5,  sourceNames.P2O5]  = localGetLiqOxOptional(data_liq, 'P2O5',  0);
[liqInputs.SO3,   sourceNames.SO3]   = localGetLiqOxOptional(data_liq, 'SO3',   0);
[liqInputs.F,     sourceNames.F]     = localGetLiqOxOptional(data_liq, 'F',     0);
[liqInputs.Cl,    sourceNames.Cl]    = localGetLiqOxOptional(data_liq, 'Cl',    0);
[liqInputs.Fe2O3, sourceNames.Fe2O3] = localGetLiqOxOptional(data_liq, 'Fe2O3', 0);
[liqInputs.H2O,   sourceNames.H2O]   = localGetLiqOxOptional(data_liq, 'H2O',   0);

end

function nanInputNames = findNaNInputs(liqInputs, sourceNames)
% findNaNInputs
% Return exact source-column labels for NaN values used in the thermometer
% calculation. Missing columns that were assigned zero are not listed.

calculationFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'Cr2O3', 'NiO', 'P2O5', 'Fe2O3'};

nameBuffer = strings(numel(calculationFields), 1);
nNames = 0;

for i = 1:numel(calculationFields)
    fieldName = calculationFields{i};
    fieldValue = liqInputs.(fieldName);

    if isnan(fieldValue)
        nNames = nNames + 1;
        sourceName = sourceNames.(fieldName);
        if strlength(sourceName) == 0
            sourceName = string(fieldName) + " (column missing)";
        end
        nameBuffer(nNames) = "Liquid." + sourceName + " = NaN";
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonnegativeInputs(liqInputs, sourceNames)
% validateNonnegativeInputs
% Stop when a finite liquid input is negative or when an input is Inf. NaN
% is intentionally allowed so that it propagates and is reported later.

allFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', ...
    'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', 'F', ...
    'Cl', 'Fe2O3', 'H2O'};

invalidBuffer = strings(numel(allFields), 1);
nInvalid = 0;

for i = 1:numel(allFields)
    fieldName = allFields{i};
    fieldValue = liqInputs.(fieldName);

    if isinf(fieldValue) || (isfinite(fieldValue) && fieldValue < 0)
        nInvalid = nInvalid + 1;
        sourceName = sourceNames.(fieldName);
        if strlength(sourceName) == 0
            sourceName = string(fieldName);
        end
        invalidBuffer(nInvalid) = ...
            "Liquid." + sourceName + " = " + string(fieldValue);
    end
end

if nInvalid > 0
    invalidInputNames = invalidBuffer(1:nInvalid);
    error(['Molina2015Liq: liquid-composition values must be finite or ' ...
           'NaN and must not be negative. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_min, liqInputs, sourceNames, P_kbar, MWinfo)
% calcTemp
% Calculate the Molina et al. (2015) liquid-only temperature. One output row
% is returned per pressure value, although pressure is not used in the
% thermometer equation.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();

% Store pressure for compatibility with fixed-P and range-P launchers.
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.P_used_in_equation = false(nP, 1);

% Expand scalar liquid inputs once to the pressure-vector length. This keeps
% every output variable at a fixed and consistent height.
SiO2  = repmat(liqInputs.SiO2,  nP, 1);
TiO2  = repmat(liqInputs.TiO2,  nP, 1);
Al2O3 = repmat(liqInputs.Al2O3, nP, 1);
FeO   = repmat(liqInputs.FeO,   nP, 1);
MnO   = repmat(liqInputs.MnO,   nP, 1);
MgO   = repmat(liqInputs.MgO,   nP, 1);
CaO   = repmat(liqInputs.CaO,   nP, 1);
Na2O  = repmat(liqInputs.Na2O,  nP, 1);
K2O   = repmat(liqInputs.K2O,   nP, 1);
V2O3  = repmat(liqInputs.V2O3,  nP, 1);
Cr2O3 = repmat(liqInputs.Cr2O3, nP, 1);
NiO   = repmat(liqInputs.NiO,   nP, 1);
P2O5  = repmat(liqInputs.P2O5,  nP, 1);
SO3   = repmat(liqInputs.SO3,   nP, 1);
F     = repmat(liqInputs.F,     nP, 1);
Cl    = repmat(liqInputs.Cl,    nP, 1);
Fe2O3 = repmat(liqInputs.Fe2O3, nP, 1);
H2O   = repmat(liqInputs.H2O,   nP, 1);

% Store original liquid inputs for inspection and reproducibility.
row.SiO2_liq = SiO2;
row.TiO2_liq = TiO2;
row.Al2O3_liq = Al2O3;
row.FeO_liq = FeO;
row.Fe_input_column = repmat(sourceNames.FeO, nP, 1);
row.MnO_liq = MnO;
row.MgO_liq = MgO;
row.CaO_liq = CaO;
row.Na2O_liq = Na2O;
row.K2O_liq = K2O;
row.V2O3_liq = V2O3;
row.Cr2O3_liq = Cr2O3;
row.NiO_liq = NiO;
row.P2O5_liq = P2O5;
row.SO3_liq = SO3;
row.F_liq = F;
row.Cl_liq = Cl;
row.Fe2O3_liq = Fe2O3;
row.H2O_liq = H2O;

% Convert oxide wt% to anhydrous cation moles. Any NaN input remains NaN
% through the arithmetic and therefore makes SumCat_liq and temperature NaN.
n_Si = SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
n_Ti = TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
n_Al = Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n_Cr = Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
n_Ni = NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;
n_Fe = FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO ...
     + Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;
n_Mn = MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
n_Mg = MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
n_Ca = CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
n_Na = Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
n_K = K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
n_P = P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;

SumCat_liq = n_Si + n_Ti + n_Al + n_Cr + n_Ni + n_Fe + n_Mn ...
    + n_Mg + n_Ca + n_Na + n_K + n_P;

% Cation mole fractions. No NaN-to-zero conversion or early replacement is
% performed. Zero denominators and zero logarithm arguments naturally yield
% NaN or Inf, which are retained and reported by the caller.
XMg_liq = n_Mg ./ SumCat_liq;
XCa_liq = n_Ca ./ SumCat_liq;
XAl_liq = n_Al ./ SumCat_liq;

ln_XMg_liq = log(XMg_liq);
ln_XCa_over_XCaPlusXAl = log(XCa_liq ./ (XCa_liq + XAl_liq));

% Molina et al. (2015) Mg-in-liquid thermometer.
T_Molina2015Liq_C = 107 .* ln_XMg_liq ...
    - 108 .* ln_XCa_over_XCaPlusXAl + 1184;
T_Molina2015Liq_K = T_Molina2015Liq_C + 273.15;

% Store cation moles, mole fractions, logarithm terms, and temperatures.
row.nSi_liq = n_Si;
row.nTi_liq = n_Ti;
row.nAl_liq = n_Al;
row.nCr_liq = n_Cr;
row.nNi_liq = n_Ni;
row.nFe_liq = n_Fe;
row.nMn_liq = n_Mn;
row.nMg_liq = n_Mg;
row.nCa_liq = n_Ca;
row.nNa_liq = n_Na;
row.nK_liq = n_K;
row.nP_liq = n_P;
row.SumCat_liq = SumCat_liq;

row.XMg_liq = XMg_liq;
row.XCa_liq = XCa_liq;
row.XAl_liq = XAl_liq;
row.ln_XMg_liq = ln_XMg_liq;
row.ln_XCa_over_XCaPlusXAl = ln_XCa_over_XCaPlusXAl;

row.T_Molina2015Liq_C = T_Molina2015Liq_C;
row.T_Molina2015Liq_K = T_Molina2015Liq_K;

% Standard aliases used by downstream plotting/export routines.
row.T_degreeC = T_Molina2015Liq_C;
row.T_K = T_Molina2015Liq_K;
row.T_deg = T_Molina2015Liq_C;

% The mineral row is not used in the equation but is retained for interface
% traceability and consistency with other interactive thermometers.
row.trace_hasMineralRow = true(nP, 1);
row.trace_mineralRowHeight = repmat(height(data_min), nP, 1);

end

function row = localAttachLiquidIDs(row, data_liq)
% localAttachLiquidIDs
% Attach common liquid identifiers and expand scalar metadata to the output
% table height.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    rawValue = data_liq.('Index');
    if isnumeric(rawValue) || islogical(rawValue)
        row.liq_Index = repmat(rawValue(1), nRows, 1);
    else
        row.liq_Index = repmat(localToScalarString(rawValue), nRows, 1);
    end
end
if any(strcmp(variableNames, 'Experiment'))
    rawValue = data_liq.('Experiment');
    row.liq_Experiment = repmat(localToScalarString(rawValue), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    rawValue = data_liq.('Citation');
    row.liq_Citation = repmat(localToScalarString(rawValue), nRows, 1);
end

end

function [value, sourceName] = localGetLiqOxOptional(data_liq, oxide, defaultValue)
% localGetLiqOxOptional
% Return one scalar oxide value. A missing column receives defaultValue. A
% present NaN remains NaN and is never replaced with defaultValue.

columnName = localFindOxideColumn(data_liq.Properties.VariableNames, oxide);

if isempty(columnName)
    value = defaultValue;
    sourceName = "";
    return;
end

rawValue = data_liq.(columnName);
value = localToScalarDouble(rawValue);
sourceName = string(columnName);

end

function columnName = localFindOxideColumn(variableNames, oxide)
% localFindOxideColumn
% Match an oxide column after removing spaces, underscores, and hyphens.

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
    idx = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(idx)
        columnName = variableNames{idx};
        return;
    end
end

end

function value = localToScalarString(rawValue)
% localToScalarString
% Convert the first table entry to one string without truncating char rows.

value = string(missing);

if isempty(rawValue)
    return;
end

if ischar(rawValue)
    value = string(strtrim(rawValue(1, :)));
    return;
end

if iscell(rawValue)
    if isempty(rawValue{1})
        return;
    end
    value = string(rawValue{1});
    return;
end

value = string(rawValue(1));

end

function value = localToScalarDouble(rawValue)
% localToScalarDouble
% Convert the first table entry to a scalar double. Missing, empty, or
% unparseable entries become NaN so that missingness is preserved.

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
    value = str2double(rawValue(1));
    return;
end

if ischar(rawValue)
    if isempty(strtrim(rawValue))
        return;
    end
    value = str2double(string(rawValue));
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
    if isstring(firstValue)
        if ismissing(firstValue(1))
            return;
        end
        value = str2double(firstValue(1));
        return;
    end
    if ischar(firstValue)
        if isempty(strtrim(firstValue))
            return;
        end
        value = str2double(string(firstValue));
        return;
    end
end

end
