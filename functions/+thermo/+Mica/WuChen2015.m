function results = WuChen2015(rawdata_struct, P_kbar)
% functions/+thermo/+Mica/WuChen2015.m
% Tested with MATLAB R2024b
%
% Revised Ti-in-biotite thermometer for rutile- or ilmenite-bearing
% crustal metapelites
% Wu, C.-M. and Chen, H.-X. (2015)
% Science Bulletin, 60, 116–121
% DOI: https://doi.org/10.1007/s11434-014-0674-y
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Biotite analysis from
% rawdata_struct.Mica and calculates temperature using the empirical
% Ti-in-biotite thermometer of Wu and Chen (2015).
%
% The function is designed for repeated calculations. Each calculation is
% stored temporarily as one table block, and all blocks are concatenated only
% once after the interactive loop has finished.
%
% Both startThermoCalc_fixedP and startThermoCalc_rangeP are supported.
% P_kbar may be a finite non-negative numeric scalar or vector. One output row
% is returned for every supplied pressure value for each selected analysis.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Wu and Chen (2015) revised the Ti-in-biotite thermometer using 334 natural
% rutile- and/or ilmenite-bearing metapelitic samples collected worldwide.
% Twenty-four samples were discarded because of possible disequilibrium, and
% 310 samples were retained for calibration (p. 117).
%
% The reported calibration ranges are:
%
%   Temperature : 450–840 degreeC
%   Pressure    : 0.1–1.9 GPa (1–19 kbar)
%   XTi         : 0.02–0.14
%   XFe         : 0.19–0.55
%   XMg         : 0.23–0.67
%   Rock type   : TiO2-saturated crustal metapelites
%   Ti phases   : rutile and/or ilmenite
%   Facies      : greenschist to granulite facies
%
% The calibration dataset and compositional ranges are described on p. 117.
% The empirical thermometer equation is Eq. 3 on p. 117. Error analysis is
% presented on pp. 117–118. Applications to prograde and inverted terrains
% and thermal contact aureoles are discussed on pp. 118–120. The empirical
% nature, pressure term, and remaining theoretical uncertainty are discussed
% on p. 120. The final recommended applicability is summarized on p. 120.
%
% IMPORTANT APPLICATION LIMITATIONS
% - The thermometer is intended for TiO2-saturated metapelites containing
%   rutile and/or ilmenite. Biotite from Ti-undersaturated rocks, igneous
%   rocks, mafic rocks, hydrothermal assemblages, or mantle phlogopite is
%   outside the direct calibration system.
% - The equation is empirical and lacks a strict thermodynamic and
%   crystallographic basis. Extrapolation beyond the calibration P-T and
%   compositional ranges is therefore not recommended (pp. 117 and 120).
% - The calibration assumes Fe3+ = 11.6 mol% of total Fe in biotite. This
%   implementation follows that calibration convention for every analysis:
%   Fe_cation_apfu must represent total Fe on an 11 O basis, Fe3+ is set to
%   0.116*Fetotal, and Fe2+ is set to 0.884*Fetotal. A separately supplied
%   Fe3_cation_apfu value is retained for reporting only and is not used in
%   the thermometer calculation.
% - Mg, total Fe, Al, Si, and Ti must all be normalized consistently to
%   11 oxygen atoms per formula unit. Al and Si are required because AlVI is
%   included in the denominator of XTi, XFe, and XMg.
% - The 310 calibration P-T values were established using garnet-biotite
%   thermometry and GASP barometry. Possible disequilibrium samples were
%   excluded, so altered, zoned, retrogressed, or texturally unrelated
%   biotite should be treated cautiously (p. 117).
% - Seventy-two percent of calibration samples were reproduced within
%   +/-50 degreeC of the reference garnet-biotite temperatures. The authors
%   estimated an approximate random error of +/-65 degreeC for natural
%   applications, but the absolute error could not be evaluated because of
%   limited experimental data (pp. 117–118).
% - An input-pressure error of +/-0.2 GPa can produce approximately
%   +/-14–25 degreeC temperature error. Analytical errors of +/-2 percent in
%   Ti, Fe, and Mg produce smaller estimated temperature effects (p. 117).
% - The authors reported two applications in which this and an earlier
%   Ti-in-biotite thermometer indicated an unexpected temperature decrease
%   between adjacent metamorphic zones; the causes remained unknown
%   (p. 120). Geological consistency should therefore be checked.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 0.1–1.9 GPa,
%   2) finite calculated temperature is outside 450–840 degreeC,
%   3) finite XTi, XFe, or XMg is outside the calibration range,
%   4) a required thermometer input is NaN,
%   5) AlIV, AlVI, the site-fraction denominator, or a logarithm argument is
%      invalid for calculation, or
%   6) the calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Mica : table
%
% The FIRST column is treated as an identifier ("data code") displayed in
% the selection dialog. The remaining columns must contain cation data
% normalized consistently to 11 O.
%
% Required variables used directly in the thermometer:
%   Mg_cation_apfu
%   Fe_cation_apfu       % total Fe, not Fe2+ only
%   Ti_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%
% Optional variables retained in the output:
%   Fe3_cation_apfu      % reporting only; not used in the equation
%   Mn_cation_apfu
%   Ca_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Missing optional variables are stored as NaN, not zero. NaN values in the
% required variables are retained, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Finite negative values are not
% allowed in any mineral-composition variable read by this function. Zero is
% retained; if it makes a logarithm or division undefined, the corresponding
% result remains NaN and is reported without stopping the calculation.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Wu and Chen (2015), Eq. 3 (p. 117):
%
%   ln[T(degreeC)] = 6.313
%                  + 0.224*ln(XTi)
%                  - 0.288*ln(XFe)
%                  - 0.449*ln(XMg)
%                  + 0.15*P(GPa)
%
% where:
%
%   Xj = j / (Fe2 + Mg + AlVI + Ti)
%
% and this implementation uses:
%
%   Fe3 = 0.116*Fetotal
%   Fe2 = 0.884*Fetotal
%   AlIV = 4 - Si
%   AlVI = Altotal - AlIV
%
% Temperature is obtained as:
%
%   T(degreeC) = exp{ln[T(degreeC)]}
%
% -------------------------------------------------------------------------
% Syntax:
%   results = WuChen2015(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing the Mica table described above
%   P_kbar         : pressure in kbar (finite non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Biotite analysis
%

%% Input validation
% Basic checks prevent silent failures due to missing arguments or invalid
% launcher pressure values.
if nargin < 2
    error('WuChen2015 requires (rawdata_struct, P_kbar).');
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

requiredVariables = { ...
    'Mg_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Al_cation_apfu', ...
    'Si_cation_apfu'};

missingRequiredVariables = requiredVariables( ...
    ~ismember(requiredVariables, dataset_bt.Properties.VariableNames));

if ~isempty(missingRequiredVariables)
    error(['rawdata_struct.Mica is missing required variable(s): ' ...
        char(strjoin(string(missingRequiredVariables), ', ')) '.']);
end

% Prepare dialog strings once instead of rebuilding them in each iteration.
dataCodes_bt = dataset_bt{:, 1};
dataCodeList_bt = cellstr(string(dataCodes_bt));

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation is stored as one table block. Repeated concatenation of
% the full results table inside the interactive loop is avoided.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Calibration limits from Wu and Chen (2015, pp. 116–120).
calibrationT_min_degC = 450;
calibrationT_max_degC = 840;
calibrationP_min_GPa = 0.1;
calibrationP_max_GPa = 1.9;
calibrationXTi_min = 0.02;
calibrationXTi_max = 0.14;
calibrationXFe_min = 0.19;
calibrationXFe_max = 0.55;
calibrationXMg_min = 0.23;
calibrationXMg_max = 0.67;

P_GPa_input = P_kbar ./ 10;
pressureOutsideCalibration = ...
    P_GPa_input < calibrationP_min_GPa | ...
    P_GPa_input > calibrationP_max_GPa;
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

    % Identify NaN values in variables used directly by the thermometer.
    % NaN values do not stop the calculation and are not replaced with zero.
    nanInputNames = findNaNInputs(selectedData_bt);

    % Reject finite negative input values. Zero and NaN are retained so that
    % undefined derived quantities remain visible and diagnosable.
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
            ': WuChen2015 = ' num2str(row.TWuChen2015_C) ' degreeC']);
    else
        disp([char(string(selectedCode_bt)) ...
            ': WuChen2015 = ' num2str(row.TWuChen2015_C(1)) ...
            ' to ' num2str(row.TWuChen2015_C(end)) ' degreeC']);
    end

    % Warn once when any supplied pressure lies outside 0.1–1.9 GPa.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the calibration range of ' ...
             'Wu and Chen (2015): 0.1–1.9 GPa (1–19 kbar). %d of %d ' ...
             'pressure point(s) are outside the range; input range = ' ...
             '%.4g–%.4g GPa.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_GPa_input), ...
            min(P_GPa_input), ...
            max(P_GPa_input));
        pressureWarningIssued = true;
    end

    % Warn when finite calculated temperatures lie outside 450–840 degreeC.
    finiteTemperature = isfinite(row.TWuChen2015_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.TWuChen2015_C < calibrationT_min_degC | ...
         row.TWuChen2015_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.TWuChen2015_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the calibration ' ...
             'range of Wu and Chen (2015): 450–840 degreeC. %d of %d ' ...
             'finite temperature point(s) are outside the range; calculated ' ...
             'finite range = %.4g–%.4g degreeC for %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_bt)));
    end

    % Warn when finite XTi lies outside 0.02–0.14.
    printRangeWarning( ...
        row.XTi, calibrationXTi_min, calibrationXTi_max, ...
        'XTi', selectedCode_bt);

    % Warn when finite XFe lies outside 0.19–0.55.
    printRangeWarning( ...
        row.XFe, calibrationXFe_min, calibrationXFe_max, ...
        'XFe', selectedCode_bt);

    % Warn when finite XMg lies outside 0.23–0.67.
    printRangeWarning( ...
        row.XMg, calibrationXMg_min, calibrationXMg_max, ...
        'XMg', selectedCode_bt);

    % Display the exact required input variables that contained NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'temperature may therefore be NaN.\n'], ...
            char(string(selectedCode_bt)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Report invalid derived quantities used to form XTi, XFe, and XMg.
    invalidDerived = ...
        ~isfinite(row.AlIV_bt) | row.AlIV_bt < 0 | ...
        ~isfinite(row.AlVI_bt) | row.AlVI_bt < 0 | ...
        ~isfinite(row.denominator_X) | row.denominator_X <= 0;

    if any(invalidDerived)
        fprintf(2, ...
            ['WARNING: Wu and Chen (2015) site fractions could not be ' ...
             'evaluated from a valid finite AlIV, AlVI, and positive ' ...
             'denominator for %s (%d of %d points). The affected X values ' ...
             'and temperatures remain NaN.\n'], ...
            char(string(selectedCode_bt)), ...
            sum(invalidDerived), ...
            numel(invalidDerived));
    end

    % Report non-positive or non-finite logarithm arguments.
    invalidLogArguments = ...
        ~isfinite(row.XTi) | row.XTi <= 0 | ...
        ~isfinite(row.XFe) | row.XFe <= 0 | ...
        ~isfinite(row.XMg) | row.XMg <= 0;

    if any(invalidLogArguments)
        fprintf(2, ...
            ['WARNING: One or more logarithm arguments (XTi, XFe, XMg) are ' ...
             'non-positive or non-finite for %s (%d of %d points). The ' ...
             'affected logarithms and temperatures remain NaN.\n'], ...
            char(string(selectedCode_bt)), ...
            sum(invalidLogArguments), ...
            numel(invalidLogArguments));
    end

    % Retain and report NaN/Inf temperature results without stopping.
    invalidTemperature = ~isfinite(row.TWuChen2015_C);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for ' ...
             '%s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_bt)), ...
            sum(invalidTemperature), ...
            numel(row.TWuChen2015_C), ...
            sum(isnan(row.TWuChen2015_C)), ...
            sum(isinf(row.TWuChen2015_C)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another biotite analysis.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'WuChen2015', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered result blocks only once after all selections finish.
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
% is selected from a fixed-size logical mask and does not grow in the loop.

requiredVariables = { ...
    'Mg_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Al_cation_apfu', ...
    'Si_cation_apfu'};

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
% Stop calculation when a supplied mica-composition value is negative,
% nonnumeric, nonscalar, or Inf. NaN and zero are intentionally allowed.

variablesToCheck = { ...
    'Mg_cation_apfu', 'Fe_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Si_cation_apfu', 'Fe3_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu', 'K_cation_apfu', ...
    'Na_cation_apfu'};

variableExists = ismember(variablesToCheck, ...
    data_biotite.Properties.VariableNames);
invalidNames = strings(numel(variablesToCheck), 1);
invalidMask = false(numel(variablesToCheck), 1);

for i = 1:numel(variablesToCheck)
    if ~variableExists(i)
        continue;
    end

    variableName = variablesToCheck{i};
    variableValue = data_biotite.(variableName);
    isInvalid = ~isnumeric(variableValue) || ~isscalar(variableValue);

    if ~isInvalid
        isInvalid = isinf(variableValue) || ...
            (isfinite(variableValue) && variableValue < 0);
    end

    if isInvalid
        invalidNames(i) = "Mica." + string(variableName);
        invalidMask(i) = true;
    end
end

if any(invalidMask)
    invalidNames = invalidNames(invalidMask);
    error(['WuChen2015: mica-composition values must be numeric scalars ' ...
        'that are >= 0 or NaN; Inf is not permitted. Invalid variable(s): ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_biotite, P_kbar)
% calcTemp
% Calculate the Wu and Chen (2015) Ti-in-biotite temperature for one selected
% biotite row and every supplied pressure value.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();

% --- Wu and Chen (2015), Eq. 3 coefficients ---
a0 = 6.313;
aTi = 0.224;
aFe = -0.288;
aMg = -0.449;
aP = 0.15;

% --- Extract one-row mineral data ---
bt = prepareMineralRow(data_biotite, 'Mica');

% Replicate pressure-independent compositions to match the pressure vector.
Mg_bt = repmat(bt.Mg, nP, 1);
Fe_total_bt = repmat(bt.Fe_total, nP, 1);
Ti_bt = repmat(bt.Ti, nP, 1);
Al_bt = repmat(bt.Al, nP, 1);
Si_bt = repmat(bt.Si, nP, 1);
Fe3_input_bt = repmat(bt.Fe3_input, nP, 1);
Mn_bt = repmat(bt.Mn, nP, 1);
Ca_bt = repmat(bt.Ca, nP, 1);
K_bt = repmat(bt.K, nP, 1);
Na_bt = repmat(bt.Na, nP, 1);

% --- Iron partitioning used in the original calibration ---
Fe3_assumed_bt = 0.116 .* Fe_total_bt;
Fe2_assumed_bt = 0.884 .* Fe_total_bt;

% --- Al site allocation on the 11 O basis ---
AlIV_bt = 4 - Si_bt;
AlVI_bt = Al_bt - AlIV_bt;

% --- Site-fraction denominator and fractions ---
denominator_X = Fe2_assumed_bt + Mg_bt + AlVI_bt + Ti_bt;

XTi = NaN(nP, 1);
XFe = NaN(nP, 1);
XMg = NaN(nP, 1);

validSiteFractions = ...
    isfinite(AlIV_bt) & AlIV_bt >= 0 & ...
    isfinite(AlVI_bt) & AlVI_bt >= 0 & ...
    isfinite(denominator_X) & denominator_X > 0;

XTi(validSiteFractions) = ...
    Ti_bt(validSiteFractions) ./ denominator_X(validSiteFractions);
XFe(validSiteFractions) = ...
    Fe2_assumed_bt(validSiteFractions) ./ denominator_X(validSiteFractions);
XMg(validSiteFractions) = ...
    Mg_bt(validSiteFractions) ./ denominator_X(validSiteFractions);

% --- Logarithms and Eq. 3 ---
lnXTi = NaN(nP, 1);
lnXFe = NaN(nP, 1);
lnXMg = NaN(nP, 1);

validLogArguments = ...
    isfinite(XTi) & XTi > 0 & ...
    isfinite(XFe) & XFe > 0 & ...
    isfinite(XMg) & XMg > 0;

lnXTi(validLogArguments) = log(XTi(validLogArguments));
lnXFe(validLogArguments) = log(XFe(validLogArguments));
lnXMg(validLogArguments) = log(XMg(validLogArguments));

lnT_C = NaN(nP, 1);
lnT_C(validLogArguments) = ...
    a0 ...
    + aTi .* lnXTi(validLogArguments) ...
    + aFe .* lnXFe(validLogArguments) ...
    + aMg .* lnXMg(validLogArguments) ...
    + aP .* P_GPa(validLogArguments);

TWuChen2015_C = NaN(nP, 1);
validLnT = isfinite(lnT_C);
TWuChen2015_C(validLnT) = exp(lnT_C(validLnT));

invalidFiniteTemperature = ...
    isfinite(TWuChen2015_C) & TWuChen2015_C <= 0;
TWuChen2015_C(invalidFiniteTemperature) = NaN;
T_K = TWuChen2015_C + 273.15;

% --- Calibration-range flags ---
is_P_in_calibration_range = ...
    P_GPa >= 0.1 & P_GPa <= 1.9;
is_T_in_calibration_range = isfinite(TWuChen2015_C) & ...
    TWuChen2015_C >= 450 & TWuChen2015_C <= 840;
is_XTi_in_calibration_range = isfinite(XTi) & ...
    XTi >= 0.02 & XTi <= 0.14;
is_XFe_in_calibration_range = isfinite(XFe) & ...
    XFe >= 0.19 & XFe <= 0.55;
is_XMg_in_calibration_range = isfinite(XMg) & ...
    XMg >= 0.23 & XMg <= 0.67;

% --- Pack outputs ---
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

row.Mg_bt = Mg_bt;
row.Fe_total_bt = Fe_total_bt;
row.Fe_bt = Fe_total_bt;
row.Fe2_bt = Fe2_assumed_bt;
row.Fe3_bt = Fe3_assumed_bt;
row.Fe3_input_bt = Fe3_input_bt;
row.Fe3_assumption_used = true(nP, 1);
row.Al_bt = Al_bt;
row.AlIV_bt = AlIV_bt;
row.AlVI_bt = AlVI_bt;
row.Si_bt = Si_bt;
row.Mn_bt = Mn_bt;
row.Ca_bt = Ca_bt;
row.Ti_bt = Ti_bt;
row.K_bt = K_bt;
row.Na_bt = Na_bt;

row.denominator_X = denominator_X;
row.XTi = XTi;
row.XFe = XFe;
row.XMg = XMg;
row.lnXTi = lnXTi;
row.lnXFe = lnXFe;
row.lnXMg = lnXMg;
row.lnT_C = lnT_C;

row.T_K = T_K;
row.T_deg = TWuChen2015_C;
row.TWuChen2015_C = TWuChen2015_C;

row.is_P_in_calibration_range = is_P_in_calibration_range;
row.is_T_in_calibration_range = is_T_in_calibration_range;
row.is_XTi_in_calibration_range = is_XTi_in_calibration_range;
row.is_XFe_in_calibration_range = is_XFe_in_calibration_range;
row.is_XMg_in_calibration_range = is_XMg_in_calibration_range;

end

function printRangeWarning(values, minimumValue, maximumValue, ...
        variableLabel, selectedCode)
% printRangeWarning
% Print one non-stopping calibration-range warning for a finite composition
% variable when one or more points lie outside the specified range.

finiteMask = isfinite(values);
outsideRange = finiteMask & ...
    (values < minimumValue | values > maximumValue);

if any(outsideRange)
    finiteValues = values(finiteMask);
    fprintf(2, ...
        ['WARNING: %s is outside the Wu and Chen (2015) calibration ' ...
         'range %.4g–%.4g. %d of %d finite point(s) are outside; ' ...
         'finite range = %.4g–%.4g for %s.\n'], ...
        variableLabel, ...
        minimumValue, ...
        maximumValue, ...
        sum(outsideRange), ...
        sum(finiteMask), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(string(selectedCode)));
end

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
mineral.Fe_total = getRequiredVar( ...
    data_tbl, 'Fe_cation_apfu', mineralLabel);
mineral.Ti = getRequiredVar(data_tbl, 'Ti_cation_apfu', mineralLabel);
mineral.Al = getRequiredVar(data_tbl, 'Al_cation_apfu', mineralLabel);
mineral.Si = getRequiredVar(data_tbl, 'Si_cation_apfu', mineralLabel);

mineral.Fe3_input = getOptionalVar( ...
    data_tbl, 'Fe3_cation_apfu', mineralLabel);
mineral.Mn = getOptionalVar(data_tbl, 'Mn_cation_apfu', mineralLabel);
mineral.Ca = getOptionalVar(data_tbl, 'Ca_cation_apfu', mineralLabel);
mineral.K = getOptionalVar(data_tbl, 'K_cation_apfu', mineralLabel);
mineral.Na = getOptionalVar(data_tbl, 'Na_cation_apfu', mineralLabel);

end

function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read one required scalar value. NaN is allowed and retained; Inf and finite
% negative values are rejected.

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
% Read one optional scalar value. A missing variable is stored as NaN so that
% missing information is not silently converted to zero.

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
