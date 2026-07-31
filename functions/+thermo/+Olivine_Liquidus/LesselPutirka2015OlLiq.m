function results = LesselPutirka2015OlLiq(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Olivine_Liquidus/LesselPutirka2015OlLiq.m
% Tested with MATLAB R2024b
%
% Olivine-Liquid thermometer
% Lessel, J. and Putirka, K.D. (2015)
% American Mineralogist, 100, 2163-2171
% DOI: https://doi.org/10.2138/am-2015-4732
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis with one Liquid
% analysis and calculates temperature using Equation (5) of Lessel and
% Putirka (2015).
%
% The function is compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP. P_kbar may therefore be supplied as either a
% scalar or a vector. Equation (5) itself contains no pressure term, so the
% calculated temperature is repeated for every supplied pressure value.
% Pressure is retained in the output table for downstream compatibility.
%
% The Liquid dataset is loaded using liquid.readLiquidExcel(). The current
% Liquid-selection workflow is retained so that downstream Liquid-related
% processing can be revised separately.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Olivine-Liquid pair, and stores all
% result blocks in a preallocated cell buffer before one final concatenation.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Lessel and Putirka (2015) calibrated Equation (5) specifically for
% olivine-liquid equilibria in high-Fe, relatively low-Al martian basaltic
% compositions. It should not be applied uncritically to ordinary terrestrial
% basaltic systems.
%
% The complete experimental database used in the study spans:
%
%   Temperature : 950-1540 degreeC
%   Pressure    : approximately 1 atm-2.3 GPa
%   Liquid SiO2 : 40.2-66.16 wt%
%   Liquid MgO  : 0.62-24.52 wt%
%   Liquid FeO  : 2.80-30.2 wt%
%   Liquid Al2O3: 2.97-20.5 wt%
%   Total alkali: 0.19-6.77 wt%
%
% These database limits are described in the Methodology section and Table 1
% on pp. 2164-2165. The olivine-bearing experimental studies listed in
% Table 1 span approximately:
%
%   Temperature : 960-1540 degreeC
%   Pressure    : approximately 1 atm-2.3 GPa
%
% Lessel and Putirka (2015) do not provide a separate table containing the
% exact minimum and maximum P-T-composition values of only the 115
% olivine-liquid pairs used for Equation (5). Therefore, this implementation
% uses the outer P-T limits of the olivine-bearing experimental studies as a
% practical warning range, rather than as a strict equation-specific limit.
%
% Equation (5) and its calibration statistics are presented on p. 2166.
% The calibration dataset contains 95 olivine-liquid pairs and reproduces
% temperature with RMSE = 48 K. The independent test dataset contains
% 20 pairs and predicts temperature with RMSE = 38 K.
%
% Application requires an equilibrium olivine-liquid pair. Lessel and
% Putirka (2015) recommend an Fe-Mg exchange coefficient near:
%
%   KD(Fe-Mg)Ol-Liq = 0.36
%
% for martian compositions (pp. 2166-2167). The inferred KD depends partly
% on the model used to calculate Fe3+/Fe2+ in the liquid; the authors use
% Kress and Carmichael (1991), Equation (7).
%
% Hydrous and Cl-rich experiments were included in the calibration, but the
% authors explicitly note that olivine-liquid thermometers are sensitive to
% volatile contents (p. 2165; Figure 6 on p. 2167). Results for strongly
% hydrous or Cl-rich liquids should therefore be interpreted cautiously.
%
% Whole-rock compositions should be used as liquid compositions only when
% they plausibly approximate a primitive or parental melt. Cumulates,
% fractionated rocks, mixed magmas, and rocks affected by crystal addition
% or removal may not represent an equilibrium liquid composition.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside approximately 1 atm-2.3 GPa,
%   2) a finite calculated temperature is outside 960-1540 degreeC,
%   3) a used composition value is NaN, or
%   4) a non-finite temperature is calculated.
%
% NaN values are retained as missing values. They are not converted to zero.
% If a used oxide column exists and contains NaN, that NaN is propagated
% through normalization and temperature calculation. A missing optional
% oxide column is excluded from the relevant normalization sum rather than
% being represented as a measured zero.
%
% All finite composition values actually used in normalization or Equation
% (5) must be strictly greater than zero. A finite zero or negative value
% stops the calculation with an error. NaN is deliberately allowed so that
% it remains missing and propagates to the output.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%
% The FIRST column of the Olivine table is treated as an identifier
% ("data code") displayed in the selection dialog.
%
% Required Olivine variables:
%   SiO2
%   MgO
%   FeO or FeOt
%
% Optional Olivine variables included in normalization when their columns
% exist:
%   TiO2, Al2O3, MnO, CaO, NiO
%
% Required Liquid variables:
%   SiO2
%   TiO2
%   MgO
%   Na2O
%   FeO or FeOt
%
% Optional Liquid variables included in cation-fraction normalization when
% their columns exist:
%   Al2O3, MnO, CaO, K2O, V2O3, Cr2O3, NiO, P2O5, SO3, Fe2O3
%
% F and Cl may be present in the input table, but they are excluded from the
% cation-fraction normalization denominator because they are anions.
%
% FeOt is used only when the FeO column is absent. If an FeO column exists
% but its selected value is NaN, that NaN is retained and FeOt is not used
% as a replacement.
%
% Molecular weights and cation numbers are loaded from:
%   functions/+liquid/Cationuli.xlsx
% through:
%   liquid.getMolarWeights()
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Lessel and Putirka (2015), Equation (5):
%
%   1 / T(K) =
%       6.529e-4
%     - 6.425e-4 * XTiO2_liq
%     - 1.049e-3 * XMgO_liq
%     - 4.206e-5 * ln(XFeO_liq)
%     - 4.121e-5 * ln(XNaO0.5_liq)
%     + 2.047e-3 * (XNaO0.5_liq)^2
%     - 8.807e-4 * (XSiO2_Ol)^2
%     + 2.299e-6 * DMgO_OlLiq
%
%   where:
%     DMgO_OlLiq = XMgO_Ol / XMgO_liq
%
% IMPORTANT IMPLEMENTATION NOTE
% In this implementation, XSiO2_Ol and XMgO_Ol are treated as oxide
% component fractions in olivine, calculated from oxide molar amounts,
% rather than as 4-oxygen-basis cation numbers. This preserves the behavior
% of the original supplied implementation.
%
% Notes:
% - Natural logarithm is used.
% - Liquid components are cation fractions on an anhydrous basis.
% - Olivine cations are also reported on a 4-oxygen basis for reference.
% - P_kbar is accepted for interface compatibility but is not used in Eq. 5.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = LesselPutirka2015OlLiq(rawdata_struct, P_kbar)
%   results = LesselPutirka2015OlLiq(..., 'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing an Olivine table
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Name-value option:
%   'LiquidRow'    : scalar positive integer or [] (default [])
%                    If empty, row 1 of the selected Liquid dataset is used.
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Olivine-Liquid pair. The output variable set is
%             intended to remain stable for downstream processing.
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('LesselPutirka2015OlLiq requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end
if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end

P_kbar = P_kbar(:);

%% Options
ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && fix(x) == x));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve required datasets
disp('=== Step 1: Preparing Olivine and Liquid datasets ===');

dataset_ol = rawdata_struct.Olivine;
MWinfo = liquid.getMolarWeights();

% Retain the existing upstream Liquid-dataset selection workflow.
[liqAll, metaLiq] = liquid.readLiquidExcel();

if isempty(liqAll) || ~istable(liqAll)
    error('Selected liquid dataset is empty or invalid.');
end

disp('=== Preparing Olivine and Liquid datasets has been finished ===');

%% 2) Initialize output container
% Each result is stored as one table block. Repeated concatenation of the
% complete output table inside the loop is avoided.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Practical warning limits based on the outer limits of the olivine-bearing
% experimental studies listed in Table 1 of Lessel and Putirka (2015).
calibrationT_min_degC = 960;
calibrationT_max_degC = 1540;
calibrationP_min_GPa = 0.0001;
calibrationP_max_GPa = 2.3;

P_GPa_input = P_kbar ./ 10;
pressureOutsideCalibration = ...
    P_GPa_input < calibrationP_min_GPa | ...
    P_GPa_input > calibrationP_max_GPa;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Olivine) ===');

while true
    % ----- Olivine selection -----
    dataCodes_ol = dataset_ol{:, 1};

    [selectedIdx_ol, ok] = listdlg( ...
        'PromptString', 'Please select the Olivine data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_ol, ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_ol)
        disp('Selection canceled');
        break;
    end

    selectedCode_ol = dataCodes_ol(selectedIdx_ol);
    disp(['Olivine selected: ' char(string(selectedCode_ol))]);

    % ----- Liquid selection -----
    disp('=== Step 4: Selecting the Liquid row ===');

    if isempty(liquidRowOpt)
        if height(liqAll) > 1
            fprintf(2, ...
                ['WARNING: The selected Liquid dataset contains %d rows. ' ...
                 'Row 1 is used by the current upstream Liquid workflow.\n'], ...
                height(liqAll));
        end
        selectedIdx_liq = 1;
    else
        selectedIdx_liq = liquidRowOpt;
        if selectedIdx_liq > height(liqAll)
            error(['Requested LiquidRow (%d) exceeds the number of rows ' ...
                   'in the selected Liquid dataset (%d).'], ...
                selectedIdx_liq, height(liqAll));
        end
    end

    disp(['Liquid selected: Row ' num2str(selectedIdx_liq)]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_liq = liqAll(selectedIdx_liq, :);

    % Check only values that are actually used. Existing NaN values are
    % retained and later propagated through the calculation.
    nanInputNames = findNaNInputs(selectedData_ol, selectedData_liq);

    % Finite values used in the calculation must be strictly positive.
    % NaN values are allowed and deliberately excluded from this check.
    validatePositiveInputs(selectedData_ol, selectedData_liq);

    row = calcTemp(selectedData_ol, selectedData_liq, P_kbar, MWinfo);

    nRows = height(row);

    % Store identifiers for traceability.
    row.dataCode_ol = repmat(string(selectedCode_ol), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = localAttachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_ol', 'dataRow_liq'}, 'Before', 1);

    % Store the result as one buffered table block.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_ol)) ' & Liquid row ' ...
            num2str(selectedIdx_liq) ': ' num2str(row.TEq5_C) ' degreeC']);
    else
        disp([char(string(selectedCode_ol)) ' & Liquid row ' ...
            num2str(selectedIdx_liq) ': ' num2str(row.TEq5_C(1)) ...
            ' to ' num2str(row.TEq5_C(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside the practical pressure
    % range represented by the olivine-bearing experiments.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the practical experimental range ' ...
             'represented by the olivine-bearing studies used by Lessel and ' ...
             'Putirka (2015): approximately 1 atm-2.3 GPa. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g GPa.\n' ...
             '         Equation (5) has no pressure term; this warning only ' ...
             'indicates extrapolation beyond the experimental P range.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_GPa_input), ...
            min(P_GPa_input), ...
            max(P_GPa_input));
        pressureWarningIssued = true;
    end

    % Warn when any finite temperature lies outside the practical
    % temperature range represented by the olivine-bearing experiments.
    finiteTemperature = isfinite(row.TEq5_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.TEq5_C < calibrationT_min_degC | ...
         row.TEq5_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.TEq5_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the practical experimental ' ...
             'range represented by the olivine-bearing studies used by Lessel and ' ...
             'Putirka (2015): 960-1540 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g-%.4g degreeC for %s & Liquid row %d.\n' ...
             '         This is an outer study range because the paper does not ' ...
             'report a separate exact P-T range for only the Equation (5) pairs.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_ol)), ...
            selectedIdx_liq);
    end

    % Report NaN inputs without stopping the calculation.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & ' ...
             'Liquid row %d: %s.\n' ...
             '         NaN was retained as a missing value and propagated through ' ...
             'the calculation; it was not replaced with zero.\n'], ...
            char(string(selectedCode_ol)), ...
            selectedIdx_liq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.TEq5_C);

    % Retain non-finite results in the output and report them.
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & ' ...
             'Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_ol)), ...
            selectedIdx_liq, ...
            sum(invalidTemperature), ...
            numel(row.TEq5_C), ...
            sum(isnan(row.TEq5_C)), ...
            sum(isinf(row.TEq5_C)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection (same Liquid dataset)?', ...
        'LesselPutirka2015OlLiq', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered table blocks once.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct('liquid', metaLiq);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_olivine, data_liquid)
% findNaNInputs
% Return names of used input variables that contain NaN. Missing optional
% columns are excluded from normalization and are not reported as NaN.

olivineOxides = {'SiO2', 'TiO2', 'Al2O3', 'MnO', 'MgO', 'CaO', 'NiO'};
liquidOxides = {'SiO2', 'TiO2', 'Al2O3', 'MnO', 'MgO', 'CaO', ...
    'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', 'Fe2O3'};

maxNames = numel(olivineOxides) + numel(liquidOxides) + 2;
nameBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(olivineOxides)
    oxide = olivineOxides{i};
    variableName = localFindOxideColumn( ...
        data_olivine.Properties.VariableNames, oxide);

    if ~isempty(variableName) && localValueIsNaN(data_olivine.(variableName))
        nNames = nNames + 1;
        nameBuffer(nNames) = "Olivine." + string(variableName);
    end
end

[~, olFeVariable] = localGetFeValue(data_olivine);
if localValueIsNaN(data_olivine.(olFeVariable))
    nNames = nNames + 1;
    nameBuffer(nNames) = "Olivine." + string(olFeVariable);
end

for i = 1:numel(liquidOxides)
    oxide = liquidOxides{i};
    variableName = localFindOxideColumn( ...
        data_liquid.Properties.VariableNames, oxide);

    if ~isempty(variableName) && localValueIsNaN(data_liquid.(variableName))
        nNames = nNames + 1;
        nameBuffer(nNames) = "Liquid." + string(variableName);
    end
end

[~, liqFeVariable] = localGetFeValue(data_liquid);
if localValueIsNaN(data_liquid.(liqFeVariable))
    nNames = nNames + 1;
    nameBuffer(nNames) = "Liquid." + string(liqFeVariable);
end

nanInputNames = nameBuffer(1:nNames);

end

function validatePositiveInputs(data_olivine, data_liquid)
% validatePositiveInputs
% Stop when a finite used composition value is zero or negative. NaN is
% intentionally allowed so that it propagates through the calculation.

requiredOlivine = {'SiO2', 'MgO'};
requiredLiquid = {'SiO2', 'TiO2', 'MgO', 'Na2O'};

optionalOlivine = {'TiO2', 'Al2O3', 'MnO', 'CaO', 'NiO'};
optionalLiquid = {'Al2O3', 'MnO', 'CaO', 'K2O', 'V2O3', ...
    'Cr2O3', 'NiO', 'P2O5', 'SO3', 'Fe2O3'};

maxInvalid = numel(requiredOlivine) + numel(requiredLiquid) + ...
    numel(optionalOlivine) + numel(optionalLiquid) + 2;
invalidBuffer = strings(maxInvalid, 1);
nInvalid = 0;

% Validate required Olivine columns.
for i = 1:numel(requiredOlivine)
    oxide = requiredOlivine{i};
    [value, variableName] = localGetRequiredOxide(data_olivine, oxide);

    if isinf(value) || (isfinite(value) && value <= 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Olivine." + string(variableName);
    end
end

% Validate the Olivine Fe column actually used.
[olFeValue, olFeVariable] = localGetFeValue(data_olivine);
if isinf(olFeValue) || (isfinite(olFeValue) && olFeValue <= 0)
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "Olivine." + string(olFeVariable);
end

% Validate present optional Olivine columns.
for i = 1:numel(optionalOlivine)
    oxide = optionalOlivine{i};
    [value, isPresent, variableName] = localGetOptionalOxide(data_olivine, oxide);

    if isPresent && (isinf(value) || (isfinite(value) && value <= 0))
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Olivine." + string(variableName);
    end
end

% Validate required Liquid columns.
for i = 1:numel(requiredLiquid)
    oxide = requiredLiquid{i};
    [value, variableName] = localGetRequiredOxide(data_liquid, oxide);

    if isinf(value) || (isfinite(value) && value <= 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Liquid." + string(variableName);
    end
end

% Validate the Liquid Fe column actually used.
[liqFeValue, liqFeVariable] = localGetFeValue(data_liquid);
if isinf(liqFeValue) || (isfinite(liqFeValue) && liqFeValue <= 0)
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "Liquid." + string(liqFeVariable);
end

% Validate present optional Liquid columns. F and Cl are intentionally not
% included because they are excluded from cation-fraction normalization.
for i = 1:numel(optionalLiquid)
    oxide = optionalLiquid{i};
    [value, isPresent, variableName] = localGetOptionalOxide(data_liquid, oxide);

    if isPresent && (isinf(value) || (isfinite(value) && value <= 0))
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Liquid." + string(variableName);
    end
end

if nInvalid > 0
    invalidInputNames = invalidBuffer(1:nInvalid);
    error(['LesselPutirka2015OlLiq: used composition values must be > 0. ' ...
           'Zero, negative, or infinite value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_ol, data_liq, P_kbar, MWinfo)
% calcTemp
% Compute Equation (5) for one Olivine row and one Liquid row. P_kbar may
% contain one or more values. Because Equation (5) is pressure independent,
% composition-derived quantities and temperature are repeated for each
% pressure value.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% ----- Olivine components -----
ol = localPrepareOlRow(data_ol, MWinfo);

row.XSi_ol = repmat(ol.XSi, nP, 1);
row.XTi_ol = repmat(ol.XTi, nP, 1);
row.XAl_ol = repmat(ol.XAl, nP, 1);
row.XFe_ol = repmat(ol.XFe, nP, 1);
row.XMn_ol = repmat(ol.XMn, nP, 1);
row.XMg_ol = repmat(ol.XMg, nP, 1);
row.XCa_ol = repmat(ol.XCa, nP, 1);
row.XNi_ol = repmat(ol.XNi, nP, 1);
row.cationSum_ol = repmat(ol.cationSum, nP, 1);

row.XSiO2_ol = repmat(ol.XSiO2_comp, nP, 1);
row.XMgO_ol = repmat(ol.XMgO_comp, nP, 1);
row.oxideSum_ol = repmat(ol.oxideSum, nP, 1);

% ----- Liquid components -----
liq = localPrepareLiquidRow(data_liq, MWinfo);

row.SiO2_liq = repmat(liq.SiO2, nP, 1);
row.TiO2_liq = repmat(liq.TiO2, nP, 1);
row.Al2O3_liq = repmat(liq.Al2O3, nP, 1);
row.FeO_liq = repmat(liq.FeO, nP, 1);
row.MnO_liq = repmat(liq.MnO, nP, 1);
row.MgO_liq = repmat(liq.MgO, nP, 1);
row.CaO_liq = repmat(liq.CaO, nP, 1);
row.Na2O_liq = repmat(liq.Na2O, nP, 1);
row.K2O_liq = repmat(liq.K2O, nP, 1);
row.Cr2O3_liq = repmat(liq.Cr2O3, nP, 1);
row.NiO_liq = repmat(liq.NiO, nP, 1);
row.Fe2O3_liq = repmat(liq.Fe2O3, nP, 1);

row.XSiO2_liq = repmat(liq.XSiO2, nP, 1);
row.XTiO2_liq = repmat(liq.XTiO2, nP, 1);
row.XAlO1_5_liq = repmat(liq.XAlO1_5, nP, 1);
row.XFeO_liq = repmat(liq.XFeO, nP, 1);
row.XMnO_liq = repmat(liq.XMnO, nP, 1);
row.XMgO_liq = repmat(liq.XMgO, nP, 1);
row.XCaO_liq = repmat(liq.XCaO, nP, 1);
row.XNaO0_5_liq = repmat(liq.XNaO0_5, nP, 1);
row.XKO0_5_liq = repmat(liq.XKO0_5, nP, 1);
row.XNiO_liq = repmat(liq.XNiO, nP, 1);

% ----- Derived Olivine/Liquid terms -----
XSiO2_Ol = ol.XSiO2_comp;
XMgO_Ol = ol.XMgO_comp;
DMgO_OlLiq = XMgO_Ol ./ liq.XMgO;

row.XSiO2_Ol = repmat(XSiO2_Ol, nP, 1);
row.XMgO_Ol = repmat(XMgO_Ol, nP, 1);
row.DMgO_OlLiq = repmat(DMgO_OlLiq, nP, 1);

% NaN values propagate naturally through this equation.
invT = ...
    6.529e-4 ...
    - 6.425e-4 .* liq.XTiO2 ...
    - 1.049e-3 .* liq.XMgO ...
    - 4.206e-5 .* log(liq.XFeO) ...
    - 4.121e-5 .* log(liq.XNaO0_5) ...
    + 2.047e-3 .* (liq.XNaO0_5 .^ 2) ...
    - 8.807e-4 .* (XSiO2_Ol .^ 2) ...
    + 2.299e-6 .* DMgO_OlLiq;

row.invT_Eq5 = repmat(invT, nP, 1);

% A non-positive finite reciprocal temperature is not physically usable.
% It is retained as an intermediate value, while the final temperature is
% returned as NaN and reported by the non-finite-result warning.
if isfinite(invT) && invT <= 0
    TEq5_K = NaN;
    TEq5_C = NaN;
else
    TEq5_K = 1 ./ invT;
    TEq5_C = TEq5_K - 273.15;
end

row.TEq5_K = repmat(TEq5_K, nP, 1);
row.TEq5_C = repmat(TEq5_C, nP, 1);

end

function ol = localPrepareOlRow(data_ol, MWinfo)
% localPrepareOlRow
% Convert one Olivine analysis into 4-oxygen-basis cations and the
% oxide-component fractions used by the supplied Equation (5)
% implementation. Missing optional columns are excluded from sums, while
% existing NaN values are retained and propagate.

[SiO2, ~] = localGetRequiredOxide(data_ol, 'SiO2');
[MgO, ~] = localGetRequiredOxide(data_ol, 'MgO');
[FeO, ~] = localGetFeValue(data_ol);

[TiO2, hasTiO2] = localGetOptionalOxide(data_ol, 'TiO2');
[Al2O3, hasAl2O3] = localGetOptionalOxide(data_ol, 'Al2O3');
[MnO, hasMnO] = localGetOptionalOxide(data_ol, 'MnO');
[CaO, hasCaO] = localGetOptionalOxide(data_ol, 'CaO');
[NiO, hasNiO] = localGetOptionalOxide(data_ol, 'NiO');

molSiO2 = SiO2 ./ MWinfo.MW.SiO2;
molMgO = MgO ./ MWinfo.MW.MgO;
molFeO = FeO ./ MWinfo.MW.FeO;

molTiO2 = localComponentForSum(TiO2, hasTiO2, MWinfo.MW.TiO2);
molAl2O3 = localComponentForSum(Al2O3, hasAl2O3, MWinfo.MW.Al2O3);
molMnO = localComponentForSum(MnO, hasMnO, MWinfo.MW.MnO);
molCaO = localComponentForSum(CaO, hasCaO, MWinfo.MW.CaO);
molNiO = localComponentForSum(NiO, hasNiO, MWinfo.MW.NiO);

oxySum = ...
    2 .* molSiO2 + ...
    2 .* molTiO2 + ...
    3 .* molAl2O3 + ...
    molFeO + ...
    molMnO + ...
    molMgO + ...
    molCaO + ...
    molNiO;

if isfinite(oxySum) && oxySum <= 0
    error('LesselPutirka2015OlLiq: invalid finite Olivine oxygen sum.');
end

ORF = 4 ./ oxySum;

XSi = molSiO2 .* ORF;
XTi = localReportedComponent(molTiO2, hasTiO2, ORF, 1);
XAl = localReportedComponent(molAl2O3, hasAl2O3, ORF, 2);
XFe = molFeO .* ORF;
XMn = localReportedComponent(molMnO, hasMnO, ORF, 1);
XMg = molMgO .* ORF;
XCa = localReportedComponent(molCaO, hasCaO, ORF, 1);
XNi = localReportedComponent(molNiO, hasNiO, ORF, 1);

cationSum = XSi + localContributionForReportedSum(XTi, hasTiO2) + ...
    localContributionForReportedSum(XAl, hasAl2O3) + XFe + ...
    localContributionForReportedSum(XMn, hasMnO) + XMg + ...
    localContributionForReportedSum(XCa, hasCaO) + ...
    localContributionForReportedSum(XNi, hasNiO);

oxideSum = molSiO2 + molTiO2 + molAl2O3 + molFeO + ...
    molMnO + molMgO + molCaO + molNiO;

if isfinite(oxideSum) && oxideSum <= 0
    error('LesselPutirka2015OlLiq: invalid finite Olivine oxide sum.');
end

ol = struct();
ol.XSi = XSi;
ol.XTi = XTi;
ol.XAl = XAl;
ol.XFe = XFe;
ol.XMn = XMn;
ol.XMg = XMg;
ol.XCa = XCa;
ol.XNi = XNi;
ol.cationSum = cationSum;
ol.oxideSum = oxideSum;
ol.XSiO2_comp = molSiO2 ./ oxideSum;
ol.XMgO_comp = molMgO ./ oxideSum;

end

function liq = localPrepareLiquidRow(data_liq, MWinfo)
% localPrepareLiquidRow
% Calculate Liquid cation fractions. Existing NaN values are retained.
% Missing optional oxide columns are excluded from the denominator. F and Cl
% are excluded because they are anions.

[SiO2, ~] = localGetRequiredOxide(data_liq, 'SiO2');
[TiO2, ~] = localGetRequiredOxide(data_liq, 'TiO2');
[MgO, ~] = localGetRequiredOxide(data_liq, 'MgO');
[Na2O, ~] = localGetRequiredOxide(data_liq, 'Na2O');
[FeO, ~] = localGetFeValue(data_liq);

[Al2O3, hasAl2O3] = localGetOptionalOxide(data_liq, 'Al2O3');
[MnO, hasMnO] = localGetOptionalOxide(data_liq, 'MnO');
[CaO, hasCaO] = localGetOptionalOxide(data_liq, 'CaO');
[K2O, hasK2O] = localGetOptionalOxide(data_liq, 'K2O');
[V2O3, hasV2O3] = localGetOptionalOxide(data_liq, 'V2O3');
[Cr2O3, hasCr2O3] = localGetOptionalOxide(data_liq, 'Cr2O3');
[NiO, hasNiO] = localGetOptionalOxide(data_liq, 'NiO');
[P2O5, hasP2O5] = localGetOptionalOxide(data_liq, 'P2O5');
[SO3, hasSO3] = localGetOptionalOxide(data_liq, 'SO3');
[Fe2O3, hasFe2O3] = localGetOptionalOxide(data_liq, 'Fe2O3');

nSiO2 = SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
nTiO2 = TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
nFeO = FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
nMgO = MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
nNa2O = Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;

nAl2O3 = localCationComponentForSum( ...
    Al2O3, hasAl2O3, MWinfo.Cat.Al2O3, MWinfo.MW.Al2O3);
nMnO = localCationComponentForSum( ...
    MnO, hasMnO, MWinfo.Cat.MnO, MWinfo.MW.MnO);
nCaO = localCationComponentForSum( ...
    CaO, hasCaO, MWinfo.Cat.CaO, MWinfo.MW.CaO);
nK2O = localCationComponentForSum( ...
    K2O, hasK2O, MWinfo.Cat.K2O, MWinfo.MW.K2O);
nV2O3 = localCationComponentForSum( ...
    V2O3, hasV2O3, MWinfo.Cat.V2O3, MWinfo.MW.V2O3);
nCr2O3 = localCationComponentForSum( ...
    Cr2O3, hasCr2O3, MWinfo.Cat.Cr2O3, MWinfo.MW.Cr2O3);
nNiO = localCationComponentForSum( ...
    NiO, hasNiO, MWinfo.Cat.NiO, MWinfo.MW.NiO);
nP2O5 = localCationComponentForSum( ...
    P2O5, hasP2O5, MWinfo.Cat.P2O5, MWinfo.MW.P2O5);
nSO3 = localCationComponentForSum( ...
    SO3, hasSO3, MWinfo.Cat.SO3, MWinfo.MW.SO3);
nFe2O3 = localCationComponentForSum( ...
    Fe2O3, hasFe2O3, MWinfo.Cat.Fe2O3, MWinfo.MW.Fe2O3);

% F and Cl are intentionally absent from this sum.
totalCations = nSiO2 + nTiO2 + nAl2O3 + nFeO + nMnO + nMgO + ...
    nCaO + nNa2O + nK2O + nV2O3 + nCr2O3 + nNiO + ...
    nP2O5 + nSO3 + nFe2O3;

if isfinite(totalCations) && totalCations <= 0
    error('LesselPutirka2015OlLiq: invalid finite Liquid cation sum.');
end

liq = struct();

% Raw wt% values. Missing optional columns are reported as NaN.
liq.SiO2 = SiO2;
liq.TiO2 = TiO2;
liq.Al2O3 = Al2O3;
liq.FeO = FeO;
liq.MnO = MnO;
liq.MgO = MgO;
liq.CaO = CaO;
liq.Na2O = Na2O;
liq.K2O = K2O;
liq.Cr2O3 = Cr2O3;
liq.NiO = NiO;
liq.Fe2O3 = Fe2O3;

% Cation fractions. Missing optional columns are reported as NaN.
liq.XSiO2 = nSiO2 ./ totalCations;
liq.XTiO2 = nTiO2 ./ totalCations;
liq.XAlO1_5 = localFractionOrNaN(nAl2O3, hasAl2O3, totalCations);
liq.XFeO = nFeO ./ totalCations;
liq.XMnO = localFractionOrNaN(nMnO, hasMnO, totalCations);
liq.XMgO = nMgO ./ totalCations;
liq.XCaO = localFractionOrNaN(nCaO, hasCaO, totalCations);
liq.XNaO0_5 = nNa2O ./ totalCations;
liq.XKO0_5 = localFractionOrNaN(nK2O, hasK2O, totalCations);
liq.XNiO = localFractionOrNaN(nNiO, hasNiO, totalCations);

end

function [value, variableName] = localGetRequiredOxide(data_tbl, oxide)
% localGetRequiredOxide
% Read a required oxide value. NaN is retained. A missing column is an error.

variableName = localFindOxideColumn(data_tbl.Properties.VariableNames, oxide);

if isempty(variableName)
    error('Selected row must contain variable: %s', oxide);
end

value = localToScalarDouble(data_tbl.(variableName));

end

function [value, isPresent, variableName] = localGetOptionalOxide(data_tbl, oxide)
% localGetOptionalOxide
% Read an optional oxide. Missing columns are marked absent and excluded
% from normalization; existing NaN values are retained.

variableName = localFindOxideColumn(data_tbl.Properties.VariableNames, oxide);
isPresent = ~isempty(variableName);

if isPresent
    value = localToScalarDouble(data_tbl.(variableName));
else
    value = NaN;
end

end

function [value, variableName] = localGetFeValue(data_tbl)
% localGetFeValue
% Use FeO when the FeO column exists. Use FeOt only when FeO is absent.
% An existing FeO value of NaN is retained and is not replaced by FeOt.

variableName = localFindOxideColumn(data_tbl.Properties.VariableNames, 'FeO');

if isempty(variableName)
    variableName = localFindOxideColumn(data_tbl.Properties.VariableNames, 'FeOt');
end

if isempty(variableName)
    error('Selected row must contain FeO or FeOt.');
end

value = localToScalarDouble(data_tbl.(variableName));

end

function row = localAttachLiquidIDs(row, data_liq)
% localAttachLiquidIDs
% Attach available Liquid identifiers and repeat them to match row height.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    value = localToScalarDouble(data_liq.('Index'));
    row.liq_Index = repmat(value, nRows, 1);
end
if any(strcmp(variableNames, 'Experiment'))
    value = string(data_liq.('Experiment'));
    row.liq_Experiment = repmat(value(1), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    value = string(data_liq.('Citation'));
    row.liq_Citation = repmat(value(1), nRows, 1);
end

end

function variableName = localFindOxideColumn(variableNames, oxide)
% localFindOxideColumn
% Match oxide columns while ignoring spaces, underscores, and hyphens.
% Both "SiO2" and "SiO2Value"-style names are accepted.

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
variableName = '';

for i = 1:numel(targets)
    index = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(index)
        variableName = variableNames{index};
        return
    end
end

end

function value = localToScalarDouble(rawValue)
% localToScalarDouble
% Convert the first table entry to a scalar double. Missing, empty, or
% non-convertible entries become NaN and are not replaced with zero.

value = NaN;

if isempty(rawValue)
    return
end

if isnumeric(rawValue)
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

if iscell(rawValue)
    if isempty(rawValue{1})
        return
    end

    firstValue = rawValue{1};

    if isnumeric(firstValue)
        value = double(firstValue(1));
    elseif isstring(firstValue)
        if ~ismissing(firstValue(1))
            value = str2double(firstValue(1));
        end
    elseif ischar(firstValue)
        value = str2double(string(firstValue));
    end
end

end

function tf = localValueIsNaN(rawValue)
% localValueIsNaN
% Return true when the first value is missing, empty, or converts to NaN.

tf = isnan(localToScalarDouble(rawValue));

end

function component = localComponentForSum(value, isPresent, molecularWeight)
% localComponentForSum
% Return zero only for an absent optional column so that it is excluded from
% the sum. If the column exists, its value, including NaN, is retained.

if isPresent
    component = value ./ molecularWeight;
else
    component = 0;
end

end

function component = localCationComponentForSum( ...
        value, isPresent, cationNumber, molecularWeight)
% localCationComponentForSum
% Return a cation amount for a present optional column. An absent column is
% excluded from the sum; an existing NaN propagates.

if isPresent
    component = value .* cationNumber ./ molecularWeight;
else
    component = 0;
end

end

function value = localReportedComponent(component, isPresent, ORF, multiplier)
% localReportedComponent
% Report NaN for an absent optional component; otherwise calculate its
% normalized cation value.

if isPresent
    value = multiplier .* component .* ORF;
else
    value = NaN;
end

end

function contribution = localContributionForReportedSum(value, isPresent)
% localContributionForReportedSum
% Exclude an absent optional component from the reported cation sum while
% preserving an existing NaN from a present component.

if isPresent
    contribution = value;
else
    contribution = 0;
end

end

function fraction = localFractionOrNaN(component, isPresent, total)
% localFractionOrNaN
% Report NaN for an absent optional component; otherwise calculate its
% cation fraction, retaining any existing NaN.

if isPresent
    fraction = component ./ total;
else
    fraction = NaN;
end

end
