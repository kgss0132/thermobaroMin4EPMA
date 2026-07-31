function results = Henry2005(rawdata_struct, P_kbar)
% functions/+thermo/+Mica/Henry2005.m
% Tested with MATLAB R2024b
%
% Ti-in-biotite thermometer
% Henry, D.J., Guidotti, C.V., Thomson, J.A. (2005)
% American Mineralogist, 90, 316–328
% DOI: https://doi.org/10.2138/am.2005.1498
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Biotite analysis from
% rawdata_struct.Mica and calculates temperature using the empirical
% Ti-in-biotite thermometer of Henry et al. (2005).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Biotite analysis and stores the result
% as a table block. All result blocks are concatenated only once after the
% interactive loop has finished.
%
% Both startThermoCalc_fixedP and startThermoCalc_rangeP are supported.
% P_kbar may be a finite non-negative scalar or vector. Because pressure does
% not appear explicitly in the published thermometer equation, the same
% calculated temperature is repeated for every supplied pressure value, while
% the corresponding pressure values are retained in the output table.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Henry et al. (2005) calibrated the Ti-saturation surface using 529 natural
% biotite analyses from graphitic, peraluminous metapelites containing a
% Ti-saturating mineral, specifically ilmenite and/or rutile. The principal
% calibration conditions are:
%
%   Temperature : 480–800 degreeC
%   Pressure    : 4–6 kbar (direct calibration data)
%   XMg         : 0.275–1.000
%   Ti          : 0.04–0.60 apfu, normalized to 22 O
%   Rock type   : graphitic, peraluminous metapelites
%   Assemblage  : ilmenite and/or rutile must be present to approach
%                 Ti-saturation in biotite
%
% The calibration data, 22 O normalization, and assumption that all Fe is
% treated as Fe2+ are described on p. 318. The fitted Ti-saturation surface,
% calibration ranges, estimated precision, and thermometer equation are
% given on p. 319 (Eq. 2). Natural applications at approximately 3–6 kbar
% are discussed on pp. 319–321, but 4–6 kbar is the direct calibration range.
%
% Reported precision varies across the calibration surface:
%   480–600 degreeC : approximately +/-24 degreeC
%   600–700 degreeC : approximately +/-23 degreeC
%   700–800 degreeC : approximately +/-12 degreeC
%
% These values describe the precision of the empirical surface fit and do
% not represent the complete absolute accuracy of every application.
%
% IMPORTANT APPLICATION LIMITATIONS
% - The thermometer should be applied to graphitic, peraluminous metapelites
%   containing ilmenite and/or rutile. Application to rocks lacking a
%   Ti-saturating mineral commonly underestimates temperature, in some cases
%   by approximately 60–80 degreeC (pp. 320–321).
% - Non-graphitic metapelites may yield lower or substantially more scattered
%   temperature estimates because fluid composition, oxidation state, and Ti
%   substitution mechanisms differ from the calibration assemblages (p. 321).
% - Metaluminous, low-Al, mafic, or igneous biotites are outside the intended
%   calibration system and can contain systematically different Ti contents
%   (pp. 322–326).
% - Chloritized, green retrograde, compositionally zoned, or locally
%   re-equilibrated biotite may record temperatures below the peak condition.
%   Intergranular and intragranular Ti heterogeneity should be checked before
%   interpreting the result (pp. 320–322).
% - Pressure affects Ti solubility in biotite even though Eq. 2 has no explicit
%   pressure term. The absence of a pressure term must not be interpreted as
%   general pressure independence outside the calibration interval
%   (discussion on pp. 316–317).
% - The input Fe_cation_apfu should represent total Fe expressed as Fe2+ in a
%   biotite formula normalized to 22 O, following the published calibration.
%   A separately supplied Fe3_cation_apfu value is retained for reporting but
%   is not used in XMg or in the thermometer equation.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 4–6 kbar,
%   2) finite XMg is outside 0.275–1.000,
%   3) finite Ti is outside 0.04–0.60 apfu,
%   4) finite calculated temperature is outside 480–800 degreeC,
%   5) a required thermometer input is NaN, or
%   6) the calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Mica : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must contain
% normalized cation data.
%
% Required variables used in the thermometer:
%   Mg_cation_apfu
%   Fe_cation_apfu       % total Fe expressed as Fe2+, 22 O basis
%   Ti_cation_apfu       % Ti apfu, 22 O basis
%
% Optional variables retained in the output:
%   Al_cation_apfu
%   Si_cation_apfu
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ca_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Missing optional variables are stored as NaN, not zero. NaN values in the
% required variables are retained, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Finite negative values are not
% allowed in any mineral-composition variable read by this function.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Henry et al. (2005) fitted the Ti-saturation surface as:
%
%   ln(Ti) = a + b*T^3 + c*(XMg)^3
%
% Rearranged for temperature (Henry et al., 2005, Eq. 2, p. 319):
%
%   T(degreeC) = nthroot([ln(Ti) - a - c*(XMg)^3] / b, 3)
%
% where:
%   Ti  = Ti in biotite (apfu normalized to 22 O)
%   XMg = Mg / (Mg + Fetotal), with all Fe treated as Fe2+
%   a   = -2.3594
%   b   = 4.6482e-9
%   c   = -1.7283
%
% nthroot is used instead of the fractional power operator so that a negative
% cubic argument remains a real value. NaN and Inf values are not replaced by
% zero and remain visible in the output for diagnosis.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Henry2005(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing a Mica table (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Biotite analysis
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values. Pressure NaN/Inf is not accepted because it is an
% independent launcher input rather than a mineral-analysis missing value.
if nargin < 2
    error('Henry2005 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation dataset
% Extract the required Mica table. The table itself is not modified.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_bt = rawdata_struct.Mica;

if height(dataset_bt) == 0
    error('rawdata_struct.Mica must contain at least one analysis row.');
end

requiredVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', 'Ti_cation_apfu'};
missingRequiredVariables = requiredVariables( ...
    ~ismember(requiredVariables, dataset_bt.Properties.VariableNames));

if ~isempty(missingRequiredVariables)
    error(['rawdata_struct.Mica is missing required variable(s): ' ...
        char(strjoin(string(missingRequiredVariables), ', ')) '.']);
end

% Prepare display strings once instead of rebuilding them in every loop.
dataCodes_bt = dataset_bt{:, 1};
dataCodeList_bt = cellstr(string(dataCodes_bt));

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated concatenation of the full results table inside the loop is avoided.
% The buffer follows the Ballhaus1991 implementation and doubles only when
% its current capacity is exhausted.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct empirical calibration limits reported by Henry et al. (2005).
calibrationT_min_degC = 480;
calibrationT_max_degC = 800;
calibrationP_min_kbar = 4;
calibrationP_max_kbar = 6;
calibrationXMg_min = 0.275;
calibrationXMg_max = 1.000;
calibrationTi_min_apfu = 0.04;
calibrationTi_max_apfu = 0.60;

% Pressure is common to all selected biotite analyses in this function call,
% so the pressure-range warning is printed only once.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels the dialog or selects Finish.
disp('=== Step 3: Selecting a data code from the list (Biotite) ===');

while true
    % ----- Biotite selection -----
    [selectedIdx_bt, ok] = listdlg( ...
        'PromptString', 'Please select the Biotite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodeList_bt, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_bt)
        disp('Selection canceled');
        break;
    end

    selectedCode_bt = dataCodes_bt(selectedIdx_bt);
    disp(['Biotite selected: ' char(string(selectedCode_bt))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    selectedData_bt = dataset_bt(selectedIdx_bt, :);

    % Identify NaN values in variables that are actually used by the
    % thermometer. NaN does not stop the calculation.
    nanInputNames = findNaNInputs(selectedData_bt);

    % Finite negative values are prohibited. Zero and NaN are retained so
    % their mathematical consequences remain visible and diagnosable.
    validateNonNegativeInputs(selectedData_bt);

    row = calcTemp(selectedData_bt, P_kbar);

    % Store the selected identifier once for each pressure row.
    row.dataCode_bt = repmat(string(selectedCode_bt), height(row), 1);
    row = movevars(row, {'dataCode_bt'}, 'Before', 1);

    % Store this calculation as one buffered table block.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_bt)) ...
            ': Henry2005 = ' num2str(row.THenry2005_C) ' degreeC']);
    else
        disp([char(string(selectedCode_bt)) ...
            ': Henry2005 = ' num2str(row.THenry2005_C(1)) ...
            ' to ' num2str(row.THenry2005_C(end)) ' degreeC']);
    end

    % Warn once when any supplied pressure lies outside the direct empirical
    % calibration range of 4–6 kbar. Calculation is not stopped.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the direct calibration range ' ...
             'of Henry et al. (2005): 4–6 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when finite XMg lies outside the fitted composition range.
    finiteXMg = isfinite(row.XMg);
    XMgOutsideCalibration = finiteXMg & ...
        (row.XMg < calibrationXMg_min | row.XMg > calibrationXMg_max);

    if any(XMgOutsideCalibration)
        finiteValues = row.XMg(finiteXMg);
        fprintf(2, ...
            ['WARNING: XMg is outside the calibration range of Henry et al. ' ...
             '(2005): 0.275–1.000. %d of %d finite point(s) are outside ' ...
             'the range; finite XMg range = %.4g–%.4g for %s.\n'], ...
            sum(XMgOutsideCalibration), ...
            sum(finiteXMg), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_bt)));
    end

    % Warn when finite Ti lies outside the fitted composition range.
    finiteTi = isfinite(row.Ti_apfu_22O);
    TiOutsideCalibration = finiteTi & ...
        (row.Ti_apfu_22O < calibrationTi_min_apfu | ...
         row.Ti_apfu_22O > calibrationTi_max_apfu);

    if any(TiOutsideCalibration)
        finiteValues = row.Ti_apfu_22O(finiteTi);
        fprintf(2, ...
            ['WARNING: Ti in biotite is outside the calibration range of ' ...
             'Henry et al. (2005): 0.04–0.60 apfu on a 22 O basis. ' ...
             '%d of %d finite point(s) are outside the range; ' ...
             'finite Ti range = %.4g–%.4g apfu for %s.\n'], ...
            sum(TiOutsideCalibration), ...
            sum(finiteTi), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_bt)));
    end

    % Warn when any finite calculated temperature lies outside 480–800 degreeC.
    finiteTemperature = isfinite(row.THenry2005_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.THenry2005_C < calibrationT_min_degC | ...
         row.THenry2005_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.THenry2005_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the calibration ' ...
             'range of Henry et al. (2005): 480–800 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_bt)));
    end

    % Display the exact required input variables that contained NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'temperature may therefore be NaN.\n'], ...
            char(string(selectedCode_bt)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report any NaN/Inf result instead of replacing it or stopping.
    invalidTemperature = ~isfinite(row.THenry2005_C);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_bt)), ...
            sum(invalidTemperature), ...
            numel(row.THenry2005_C), ...
            sum(isnan(row.THenry2005_C)), ...
            sum(isinf(row.THenry2005_C)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another biotite analysis.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Henry2005', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after all selections finish.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_biotite)
% findNaNInputs
% Return names of required thermometer variables containing NaN. The output
% array has a fixed maximum size and does not grow during the loop.

requiredVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', 'Ti_cation_apfu'};
qualifiedNames = "Mica." + string(requiredVariables(:));
containsNaN = false(numel(requiredVariables), 1);

for i = 1:numel(requiredVariables)
    variableValue = data_biotite.(requiredVariables{i});
    containsNaN(i) = any(isnan(variableValue(:)));
end

nanInputNames = qualifiedNames(containsNaN);

end

function validateNonNegativeInputs(data_biotite)
% validateNonNegativeInputs
% Stop calculation when a finite mineral-composition value is negative or
% when a supplied value is nonnumeric/non-scalar/Inf. NaN and zero are
% intentionally allowed so that they remain visible in the calculation.

variablesToCheck = { ...
    'Mg_cation_apfu', 'Fe_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Si_cation_apfu', 'Fe3_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu', 'K_cation_apfu', ...
    'Na_cation_apfu'};

variableExists = ismember(variablesToCheck, ...
    data_biotite.Properties.VariableNames);
invalidShapeOrType = false(numel(variablesToCheck), 1);
containsInf = false(numel(variablesToCheck), 1);
containsNegative = false(numel(variablesToCheck), 1);

for i = 1:numel(variablesToCheck)
    if ~variableExists(i)
        continue;
    end

    variableValue = data_biotite.(variablesToCheck{i});
    invalidShapeOrType(i) = ~isnumeric(variableValue) || ~isscalar(variableValue);

    if ~invalidShapeOrType(i)
        containsInf(i) = isinf(variableValue);
        containsNegative(i) = isfinite(variableValue) && variableValue < 0;
    end
end

if any(invalidShapeOrType)
    invalidNames = string(variablesToCheck(invalidShapeOrType));
    error(['Henry2005: mineral-composition variables must contain one ' ...
        'numeric scalar per selected row. Invalid variable(s): ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

if any(containsInf)
    invalidNames = string(variablesToCheck(containsInf));
    error(['Henry2005: Inf is not permitted in mineral-composition ' ...
        'variables. Invalid variable(s): ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

if any(containsNegative)
    invalidNames = string(variablesToCheck(containsNegative));
    error(['Henry2005: mineral-composition values must be >= 0 or NaN. ' ...
        'Negative finite value(s) were found in: ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_biotite, P_kbar)
% calcTemp
% Calculate the Henry et al. (2005) Ti-in-biotite temperature for one
% selected biotite row and every supplied pressure value.
%
% Pressure is retained for interface consistency and calibration screening,
% but it is not used explicitly in Henry et al. (2005, Eq. 2).

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();

% --- Henry et al. (2005) empirical constants ---
a_const = -2.3594;
b_const = 4.6482e-9;
c_const = -1.7283;

% --- Extract one-row mineral data ---
bt = prepareMineralRow(data_biotite, 'Mica');

% Replicate pressure-independent composition values to match the supplied
% pressure vector. NaN remains NaN and is never converted to zero.
Mg_bt = repmat(bt.Mg, nP, 1);
Fe_total_as_Fe2_bt = repmat(bt.Fe, nP, 1);
Ti_bt = repmat(bt.Ti, nP, 1);

Al_bt = repmat(bt.Al, nP, 1);
Si_bt = repmat(bt.Si, nP, 1);
Fe3_bt = repmat(bt.Fe3, nP, 1);
Mn_bt = repmat(bt.Mn, nP, 1);
Ca_bt = repmat(bt.Ca, nP, 1);
K_bt = repmat(bt.K, nP, 1);
Na_bt = repmat(bt.Na, nP, 1);

% --- Core thermometer parameters ---
XMg = Mg_bt ./ (Mg_bt + Fe_total_as_Fe2_bt);
Ti_apfu_22O = Ti_bt;

lnTi = log(Ti_apfu_22O);
cubic_argument = ...
    (lnTi - a_const - c_const .* (XMg .^ 3)) ./ b_const;

% Use the real cube root. NaN and Inf propagate without replacement.
THenry2005_C = nthroot(cubic_argument, 3);
T_K = THenry2005_C + 273.15;

% --- Calibration flags ---
is_XMg_in_calibration_range = isfinite(XMg) & ...
    XMg >= 0.275 & XMg <= 1.000;
is_Ti_in_calibration_range = isfinite(Ti_apfu_22O) & ...
    Ti_apfu_22O >= 0.04 & Ti_apfu_22O <= 0.60;
is_P_in_calibration_range = P_kbar >= 4 & P_kbar <= 6;
is_T_in_calibration_range = isfinite(THenry2005_C) & ...
    THenry2005_C >= 480 & THenry2005_C <= 800;

% --- Pack outputs ---
row.P_kbar = P_kbar;

row.Mg_bt = Mg_bt;
row.Fe_total_as_Fe2_bt = Fe_total_as_Fe2_bt;
row.Fe2_bt = Fe_total_as_Fe2_bt; % Backward-compatible alias.
row.Fe3_bt = Fe3_bt;
row.Al_bt = Al_bt;
row.Si_bt = Si_bt;
row.Mn_bt = Mn_bt;
row.Ca_bt = Ca_bt;
row.Ti_bt = Ti_bt;
row.K_bt = K_bt;
row.Na_bt = Na_bt;

row.XMg = XMg;
row.Ti_apfu_22O = Ti_apfu_22O;
row.lnTi = lnTi;
row.cubic_argument = cubic_argument;

row.T_K = T_K;
row.T_deg = THenry2005_C;
row.THenry2005_C = THenry2005_C;

row.is_P_in_calibration_range = is_P_in_calibration_range;
row.is_XMg_in_calibration_range = is_XMg_in_calibration_range;
row.is_Ti_in_calibration_range = is_Ti_in_calibration_range;
row.is_T_in_calibration_range = is_T_in_calibration_range;

end

function mineral = prepareMineralRow(data_tbl, mineralLabel)
% prepareMineralRow
% Extract cation values from one selected biotite row. Required and optional
% NaN values are retained. Missing optional variables are represented by NaN.

if height(data_tbl) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();

mineral.Mg = getRequiredVar(data_tbl, 'Mg_cation_apfu', mineralLabel);
mineral.Fe = getRequiredVar(data_tbl, 'Fe_cation_apfu', mineralLabel);
mineral.Ti = getRequiredVar(data_tbl, 'Ti_cation_apfu', mineralLabel);

mineral.Al = getOptionalVar(data_tbl, 'Al_cation_apfu', mineralLabel);
mineral.Si = getOptionalVar(data_tbl, 'Si_cation_apfu', mineralLabel);
mineral.Fe3 = getOptionalVar(data_tbl, 'Fe3_cation_apfu', mineralLabel);
mineral.Mn = getOptionalVar(data_tbl, 'Mn_cation_apfu', mineralLabel);
mineral.Ca = getOptionalVar(data_tbl, 'Ca_cation_apfu', mineralLabel);
mineral.K = getOptionalVar(data_tbl, 'K_cation_apfu', mineralLabel);
mineral.Na = getOptionalVar(data_tbl, 'Na_cation_apfu', mineralLabel);

end

function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read one required scalar value. NaN is allowed and retained; Inf and
% negative finite values are rejected.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must contain one numeric scalar.', mineralLabel, varName);
end
if isinf(value)
    error('%s.%s must not contain Inf.', mineralLabel, varName);
end
if isfinite(value) && value < 0
    error('%s.%s must be >= 0 or NaN.', mineralLabel, varName);
end

end

function value = getOptionalVar(tbl, varName, mineralLabel)
% getOptionalVar
% Read one optional scalar value. Missing optional variables are represented
% by NaN so that missing information is not silently converted to zero.

if ~ismember(varName, tbl.Properties.VariableNames)
    value = NaN;
    return;
end

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must contain one numeric scalar.', mineralLabel, varName);
end
if isinf(value)
    error('%s.%s must not contain Inf.', mineralLabel, varName);
end
if isfinite(value) && value < 0
    error('%s.%s must be >= 0 or NaN.', mineralLabel, varName);
end

end
