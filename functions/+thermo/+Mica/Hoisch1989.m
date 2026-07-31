function results = Hoisch1989(rawdata_struct, P_kbar)
% functions/+thermo/+Mica/Hoisch1989.m
% Tested with MATLAB R2024b
%
% Muscovite-biotite thermometer
% Hoisch, T.D. (1989)
% American Mineralogist, 74, 565–572
% DOI: None
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Muscovite analysis and one Biotite
% analysis from rawdata_struct.Mica and calculates temperature using the
% Hoisch (1989) muscovite-biotite Mg-Tschermak exchange thermometer.
%
% Muscovite and biotite are selected independently from the same Mica table.
% The first column is used only as the displayed data code; the selected rows
% must therefore correspond to muscovite and biotite analyses identified by
% the user.
%
% The function is designed for repeated calculations. Each calculation is
% stored temporarily as one table block, and all blocks are concatenated only
% once after the interactive loop has finished.
%
% Both startThermoCalc_fixedP and startThermoCalc_rangeP are supported.
% P_kbar may be a finite non-negative numeric scalar or vector. One output row
% is returned for every supplied pressure value for each selected Ms-Bt pair.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Hoisch (1989) empirically calibrated exchange of the Mg-Tschermak
% component between coexisting muscovite and biotite. The original dataset
% contained 43 pelitic rocks with the assemblage:
%
%   Muscovite + Biotite + Quartz + Plagioclase + Garnet + Al2SiO5
%
% Four discrepant samples were removed from the final regression, leaving
% 39 calibration samples. The final calibration conditions are:
%
%   Temperature : 450–700 degreeC
%   Pressure    : 2000–9600 bar (2.0–9.6 kbar)
%   Rock type   : Al-rich pelitic metamorphic rocks
%   Assemblage  : Ms + Bt + Qz + Pl + Grt + Al2SiO5
%   Regression  : multiple correlation coefficient = 0.92
%   Residual SD : approximately 22 K
%
% Calibration temperatures and pressures were obtained by simultaneous
% application of garnet-biotite thermometry and GASP barometry (pp. 566–568).
% The regression, final thermometer equation, pressure sensitivity, and fit
% statistics are given on pp. 568–570; the implemented thermometer is Eq. 10
% on p. 569. Restrictions on use and the calibration-composition ranges are
% given on pp. 570–571 and in Table 4 on p. 571.
%
% IMPORTANT APPLICATION LIMITATIONS
% - Hoisch (1989) states that applications should be restricted to micas
%   compositionally similar to the calibration dataset. The calibration micas
%   are Al rich and occupy relatively narrow compositional ranges because
%   they coexist with Al2SiO5 in similar low-variance assemblages
%   (pp. 570–571).
% - Application outside 450–700 degreeC or 2.0–9.6 kbar, or to conditions
%   unrepresented in the dataset such as blueschist-facies metamorphism,
%   carries higher uncertainty (p. 570).
% - Pressure sensitivity is weak but not zero. For the end-member reaction,
%   a 1 kbar pressure uncertainty corresponds to approximately 14.8 K in
%   temperature (Eq. 11 and discussion on p. 569).
% - Muscovite and biotite must be interpreted as an equilibrium pair.
%   Polymetamorphic, zoned, altered, retrogressed, or texturally unrelated
%   mica generations may give misleading temperatures.
% - The original mica analyses were normalized to 11 anhydrous oxygens, and
%   all Fe was treated as Fe2+ (p. 566). Input cation values should follow the
%   same normalization and Fe convention for consistency with the calibration.
% - The approximately 22 K residual standard deviation describes scatter of
%   the regression and is not the complete absolute uncertainty of an
%   individual application. Systematic uncertainties in the thermometers and
%   barometer used to establish the calibration P-T values also contribute
%   to total uncertainty (pp. 569–570).
%
% Table 4 composition ranges for the 39 final-regression samples:
%   Kideal_R1             : 0.108–0.449
%   XMg_B - XAlVI_B       : 0.165–0.318
%   XMg_B                 : 0.319–0.417
%   XAlVI_B               : 0.089–0.170
%   XTi_B                 : 0.0240–0.0732
%   XFe_B                 : 0.339–0.443
%   XMg_M                 : 0.0085–0.0431
%
% XTi_B and XFe_B are used only for calibration-composition screening, not
% directly in Eq. 10. If Fe_cation_apfu or Ti_cation_apfu is missing, the
% corresponding output and range flag remain NaN rather than being replaced
% with zero.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 2.0–9.6 kbar,
%   2) finite calculated temperature is outside 450–700 degreeC,
%   3) a finite composition index is outside a Table 4 range,
%   4) a required thermometer input is NaN,
%   5) a derived quantity required for Kideal is non-positive or non-finite,
%   6) the equation denominator is zero or non-finite, or
%   7) the calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Mica : table containing both muscovite and biotite rows
%
% The FIRST column is treated as an identifier ("data code") displayed in
% the selection dialogs. The remaining columns must contain cation data
% normalized consistently with the original 11-anhydrous-oxygen basis.
%
% Required variables used directly in the thermometer:
%   Mg_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%
% Optional variables retained in the output and/or used for screening:
%   Fe_cation_apfu       % total Fe treated as Fe2+ in the original study
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Ti_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Missing optional variables are stored as NaN, not zero. NaN values in the
% required variables are retained, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Finite negative values are not
% allowed in any mineral-composition variable read by this function. Zero is
% retained; if it makes a logarithm or division undefined, the result remains
% NaN/Inf and is reported without replacing it.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Site fractions for muscovite (M) and biotite (B):
%
%   XMg_M   = Mg_M / 2
%   XMg_B   = Mg_B / 3
%   XAlVI_M = (Al_M + Si_M - 4) / 2
%   XAlVI_B = (Al_B + Si_B - 4) / 3
%   XTi_B   = Ti_B / 3
%   XFe_B   = Fe_B / 3
%
% The ideal equilibrium constant is (Hoisch, 1989, p. 568):
%
%   Kideal_R1 = 27 * (Mg_M / AlVI_M_apfu) / ...
%                    (Mg_B / AlVI_B_apfu)
%
% where:
%   AlVI_M_apfu = Al_M + Si_M - 4
%   AlVI_B_apfu = Al_B + Si_B - 4
%
% The factor 27 is required by the activity expressions in the original
% calibration and is included in this implementation.
%
% Hoisch (1989), Eq. 10:
%
%   T(K) =
%     [500.110 + 0.0147890*P(bar)
%      - 878.745*(XMg_B - XAlVI_B)
%      - 4532.67*XMg_M*(XMg_M - 2)]
%     / [1 + 0.0237527*R*ln(Kideal_R1)]
%
% where:
%   R      = 8.3144 J mol^-1 K^-1
%   P(bar) = 1000 * P_kbar
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Hoisch1989(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing the Mica table described above
%   P_kbar         : pressure in kbar (finite non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Muscovite-Biotite pair
%

%% Input validation
% Basic checks prevent silent failures due to missing arguments or invalid
% launcher pressure values.
if nargin < 2
    error('Hoisch1989 requires (rawdata_struct, P_kbar).');
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
% Muscovite and biotite are selected independently from the same Mica table.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_mica = rawdata_struct.Mica;

if height(dataset_mica) == 0
    error('rawdata_struct.Mica must contain at least one analysis row.');
end

requiredVariables = {'Mg_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu'};
missingRequiredVariables = requiredVariables( ...
    ~ismember(requiredVariables, dataset_mica.Properties.VariableNames));

if ~isempty(missingRequiredVariables)
    error(['rawdata_struct.Mica is missing required variable(s): ' ...
        char(strjoin(string(missingRequiredVariables), ', ')) '.']);
end

% Prepare dialog strings once instead of recreating them in every iteration.
dataCodes_mica = dataset_mica{:, 1};
dataCodeList_mica = cellstr(string(dataCodes_mica));

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Repeated concatenation of the full results table inside the loop is avoided.
% Each selected pair is stored as a table block and concatenated once after
% the interactive calculations finish.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Calibration limits from Hoisch (1989, pp. 570–571 and Table 4).
calibrationT_min_degC = 450;
calibrationT_max_degC = 700;
calibrationP_min_kbar = 2.0;
calibrationP_max_kbar = 9.6;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or selects
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Muscovite) ===');

while true
    % ----- Muscovite selection -----
    [selectedIdx_ms, ok] = listdlg( ...
        'PromptString', 'Please select the Muscovite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodeList_mica, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_ms)
        disp('Selection canceled');
        break;
    end

    selectedCode_ms = dataCodes_mica(selectedIdx_ms);
    disp(['Muscovite selected: ' char(string(selectedCode_ms))]);

    % ----- Biotite selection -----
    disp('=== Step 4: Selecting a data code from the list (Biotite) ===');

    [selectedIdx_bt, ok] = listdlg( ...
        'PromptString', 'Please select the Biotite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodeList_mica, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_bt)
        disp('Selection canceled');
        break;
    end

    selectedCode_bt = dataCodes_mica(selectedIdx_bt);
    disp(['Biotite selected: ' char(string(selectedCode_bt))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_ms = dataset_mica(selectedIdx_ms, :);
    selectedData_bt = dataset_mica(selectedIdx_bt, :);

    % Identify NaN values in variables used directly by Eq. 10. NaN values
    % do not stop the calculation and are not replaced with zero.
    nanInputNames = findNaNInputs(selectedData_ms, selectedData_bt);

    % Reject finite negative input values. Zero and NaN are retained so that
    % undefined derived quantities remain visible and diagnosable.
    validateNonNegativeInputs(selectedData_ms, selectedData_bt);

    row = calcTemp(selectedData_ms, selectedData_bt, P_kbar);

    % Store selected identifiers once for each pressure row.
    row.dataCode_ms = repmat(string(selectedCode_ms), height(row), 1);
    row.dataCode_bt = repmat(string(selectedCode_bt), height(row), 1);
    row = movevars(row, {'dataCode_ms', 'dataCode_bt'}, 'Before', 1);

    % Store the result as one buffered table block.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_ms)) ' & ' char(string(selectedCode_bt)) ...
            ': Hoisch1989 = ' num2str(row.THoisch1989_C) ' degreeC']);
    else
        disp([char(string(selectedCode_ms)) ' & ' char(string(selectedCode_bt)) ...
            ': Hoisch1989 = ' num2str(row.THoisch1989_C(1)) ...
            ' to ' num2str(row.THoisch1989_C(end)) ' degreeC']);
    end

    % Warn once when any supplied pressure lies outside 2.0–9.6 kbar.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the calibration range of ' ...
             'Hoisch (1989): 2.0–9.6 kbar. %d of %d pressure point(s) ' ...
             'are outside the range; input range = %.4g–%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when finite calculated temperatures lie outside 450–700 degreeC.
    finiteTemperature = isfinite(row.THoisch1989_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.THoisch1989_C < calibrationT_min_degC | ...
         row.THoisch1989_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.THoisch1989_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the calibration ' ...
             'range of Hoisch (1989): 450–700 degreeC. %d of %d finite ' ...
             'temperature point(s) are outside the range; calculated finite ' ...
             'range = %.4g–%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_ms)), ...
            char(string(selectedCode_bt)));
    end

    % Print Table 4 composition-range warnings for finite values.
    printCompositionRangeWarnings( ...
        row, selectedCode_ms, selectedCode_bt);

    % Display exact required thermometer inputs that contained NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'temperature may therefore be NaN.\n'], ...
            char(string(selectedCode_ms)), ...
            char(string(selectedCode_bt)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Report invalid derived quantities needed to calculate Kideal.
    invalidKInputs = ...
        ~isfinite(row.Mg_ms) | row.Mg_ms <= 0 | ...
        ~isfinite(row.Mg_bt) | row.Mg_bt <= 0 | ...
        ~isfinite(row.XAlVI_M_apfu) | row.XAlVI_M_apfu <= 0 | ...
        ~isfinite(row.XAlVI_B_apfu) | row.XAlVI_B_apfu <= 0;

    if any(invalidKInputs)
        fprintf(2, ...
            ['WARNING: Hoisch (1989) Kideal could not be evaluated from a ' ...
             'positive finite Mg and octahedral-Al combination for %s & %s ' ...
             '(%d of %d points). The affected Kideal, lnK, and temperature ' ...
             'values remain NaN.\n'], ...
            char(string(selectedCode_ms)), ...
            char(string(selectedCode_bt)), ...
            sum(invalidKInputs), ...
            numel(invalidKInputs));
    end

    % Report a zero or non-finite Eq. 10 denominator separately.
    invalidDenominator = ...
        ~isfinite(row.denominator) | row.denominator == 0;

    if any(invalidDenominator)
        fprintf(2, ...
            ['WARNING: The Hoisch (1989) Eq. 10 denominator is zero or ' ...
             'non-finite for %s & %s (%d of %d points). The corresponding ' ...
             'temperature remains NaN.\n'], ...
            char(string(selectedCode_ms)), ...
            char(string(selectedCode_bt)), ...
            sum(invalidDenominator), ...
            numel(invalidDenominator));
    end

    % Retain and report NaN/Inf temperature results without stopping.
    invalidTemperature = ~isfinite(row.THoisch1989_C);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for ' ...
             '%s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_ms)), ...
            char(string(selectedCode_bt)), ...
            sum(invalidTemperature), ...
            numel(row.THoisch1989_C), ...
            sum(isnan(row.THoisch1989_C)), ...
            sum(isinf(row.THoisch1989_C)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another mica pair.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Hoisch1989', ...
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
function nanInputNames = findNaNInputs(data_muscovite, data_biotite)
% findNaNInputs
% Return names of required Eq. 10 input variables containing NaN. The output
% is selected from a fixed-size logical mask and does not grow in the loop.

requiredVariables = {'Mg_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu'};
qualifiedNames = [ ...
    "Muscovite." + string(requiredVariables(:)); ...
    "Biotite." + string(requiredVariables(:))];
containsNaN = false(2 * numel(requiredVariables), 1);

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    muscoviteValue = data_muscovite.(variableName);
    biotiteValue = data_biotite.(variableName);
    containsNaN(i) = any(isnan(muscoviteValue(:)));
    containsNaN(numel(requiredVariables) + i) = ...
        any(isnan(biotiteValue(:)));
end

nanInputNames = qualifiedNames(containsNaN);

end

function validateNonNegativeInputs(data_muscovite, data_biotite)
% validateNonNegativeInputs
% Stop calculation when a supplied mica-composition value is negative,
% nonnumeric, nonscalar, or Inf. NaN and zero are intentionally allowed.

variablesToCheck = { ...
    'Mg_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu', ...
    'Fe_cation_apfu', 'Fe3_cation_apfu', 'Mn_cation_apfu', ...
    'Ca_cation_apfu', 'Ti_cation_apfu', 'K_cation_apfu', ...
    'Na_cation_apfu'};

invalidNames = strings(2 * numel(variablesToCheck), 1);
invalidMask = false(2 * numel(variablesToCheck), 1);
writeIndex = 0;

for mineralIndex = 1:2
    if mineralIndex == 1
        dataTable = data_muscovite;
        mineralName = "Muscovite";
    else
        dataTable = data_biotite;
        mineralName = "Biotite";
    end

    for i = 1:numel(variablesToCheck)
        variableName = variablesToCheck{i};
        if ~ismember(variableName, dataTable.Properties.VariableNames)
            continue;
        end

        variableValue = dataTable.(variableName);
        isInvalid = ~isnumeric(variableValue) || ~isscalar(variableValue);

        if ~isInvalid
            isInvalid = isinf(variableValue) || ...
                (isfinite(variableValue) && variableValue < 0);
        end

        if isInvalid
            writeIndex = writeIndex + 1;
            invalidNames(writeIndex) = mineralName + "." + string(variableName);
            invalidMask(writeIndex) = true;
        end
    end
end

if any(invalidMask)
    invalidNames = invalidNames(invalidMask);
    error(['Hoisch1989: mica-composition values must be numeric scalars ' ...
        'that are >= 0 or NaN; Inf is not permitted. Invalid variable(s): ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_muscovite, data_biotite, P_kbar)
% calcTemp
% Calculate the Hoisch (1989) muscovite-biotite temperature for one selected
% mica pair and every supplied pressure value.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

row = table();

% --- Constant used by Hoisch (1989) ---
R_J = 8.3144;

% --- Extract one-row mica data ---
ms = prepareMineralRow(data_muscovite, 'Muscovite');
bt = prepareMineralRow(data_biotite, 'Biotite');

% Replicate pressure-independent compositions to match the pressure vector.
Mg_ms = repmat(ms.Mg, nP, 1);
Al_ms = repmat(ms.Al, nP, 1);
Si_ms = repmat(ms.Si, nP, 1);
Fe_ms = repmat(ms.Fe, nP, 1);
Fe3_ms = repmat(ms.Fe3, nP, 1);
Mn_ms = repmat(ms.Mn, nP, 1);
Ca_ms = repmat(ms.Ca, nP, 1);
Ti_ms = repmat(ms.Ti, nP, 1);
K_ms = repmat(ms.K, nP, 1);
Na_ms = repmat(ms.Na, nP, 1);

Mg_bt = repmat(bt.Mg, nP, 1);
Al_bt = repmat(bt.Al, nP, 1);
Si_bt = repmat(bt.Si, nP, 1);
Fe_bt = repmat(bt.Fe, nP, 1);
Fe3_bt = repmat(bt.Fe3, nP, 1);
Mn_bt = repmat(bt.Mn, nP, 1);
Ca_bt = repmat(bt.Ca, nP, 1);
Ti_bt = repmat(bt.Ti, nP, 1);
K_bt = repmat(bt.K, nP, 1);
Na_bt = repmat(bt.Na, nP, 1);

% --- Site-based quantities ---
XMg_M = Mg_ms ./ 2;
XMg_B = Mg_bt ./ 3;

XAlVI_M_apfu = Al_ms + Si_ms - 4;
XAlVI_B_apfu = Al_bt + Si_bt - 4;

XAlVI_M = XAlVI_M_apfu ./ 2;
XAlVI_B = XAlVI_B_apfu ./ 3;

XTi_B = Ti_bt ./ 3;
XFe_B = Fe_bt ./ 3;
XMgB_minus_XAlVIB = XMg_B - XAlVI_B;

% --- Ideal equilibrium constant ---
% The factor 27 follows Hoisch (1989, p. 568) and must not be omitted.
Kideal_R1 = NaN(nP, 1);
validKInputs = ...
    isfinite(Mg_ms) & Mg_ms > 0 & ...
    isfinite(Mg_bt) & Mg_bt > 0 & ...
    isfinite(XAlVI_M_apfu) & XAlVI_M_apfu > 0 & ...
    isfinite(XAlVI_B_apfu) & XAlVI_B_apfu > 0;

Kideal_R1(validKInputs) = 27 .* ...
    (Mg_ms(validKInputs) ./ XAlVI_M_apfu(validKInputs)) ./ ...
    (Mg_bt(validKInputs) ./ XAlVI_B_apfu(validKInputs));

lnK = NaN(nP, 1);
validK = isfinite(Kideal_R1) & Kideal_R1 > 0;
lnK(validK) = log(Kideal_R1(validK));

% --- Hoisch (1989), Eq. 10 ---
numerator = 500.110 + 0.0147890 .* P_bar ...
    - 878.745 .* XMgB_minus_XAlVIB ...
    - 4532.67 .* (XMg_M .* (XMg_M - 2));

denominator = 1 + 0.0237527 .* R_J .* lnK;

THoisch1989_K = NaN(nP, 1);
validEquation = isfinite(numerator) & isfinite(denominator) & ...
    denominator ~= 0;
THoisch1989_K(validEquation) = ...
    numerator(validEquation) ./ denominator(validEquation);

% Non-positive Kelvin values are not treated as physical temperatures.
nonPositiveKelvin = isfinite(THoisch1989_K) & THoisch1989_K <= 0;
THoisch1989_K(nonPositiveKelvin) = NaN;
THoisch1989_C = THoisch1989_K - 273.15;

% --- Calibration-range flags ---
is_P_in_calibration_range = P_kbar >= 2.0 & P_kbar <= 9.6;
is_T_in_calibration_range = isfinite(THoisch1989_C) & ...
    THoisch1989_C >= 450 & THoisch1989_C <= 700;

is_Kideal_in_calibration_range = isfinite(Kideal_R1) & ...
    Kideal_R1 >= 0.108 & Kideal_R1 <= 0.449;
is_XMgB_minus_XAlVIB_in_calibration_range = ...
    isfinite(XMgB_minus_XAlVIB) & ...
    XMgB_minus_XAlVIB >= 0.165 & XMgB_minus_XAlVIB <= 0.318;
is_XMg_B_in_calibration_range = isfinite(XMg_B) & ...
    XMg_B >= 0.319 & XMg_B <= 0.417;
is_XAlVI_B_in_calibration_range = isfinite(XAlVI_B) & ...
    XAlVI_B >= 0.089 & XAlVI_B <= 0.170;
is_XTi_B_in_calibration_range = isfinite(XTi_B) & ...
    XTi_B >= 0.0240 & XTi_B <= 0.0732;
is_XFe_B_in_calibration_range = isfinite(XFe_B) & ...
    XFe_B >= 0.339 & XFe_B <= 0.443;
is_XMg_M_in_calibration_range = isfinite(XMg_M) & ...
    XMg_M >= 0.0085 & XMg_M <= 0.0431;

% --- Pack outputs ---
row.P_kbar = P_kbar;
row.P_bar = P_bar;
row.R_J_molK = repmat(R_J, nP, 1);

row.Mg_ms = Mg_ms;
row.Al_ms = Al_ms;
row.Si_ms = Si_ms;
row.Fe_ms = Fe_ms;
row.Fe2_ms = Fe_ms; % Backward-compatible alias.
row.Fe3_ms = Fe3_ms;
row.Mn_ms = Mn_ms;
row.Ca_ms = Ca_ms;
row.Ti_ms = Ti_ms;
row.K_ms = K_ms;
row.Na_ms = Na_ms;

row.Mg_bt = Mg_bt;
row.Al_bt = Al_bt;
row.Si_bt = Si_bt;
row.Fe_bt = Fe_bt;
row.Fe2_bt = Fe_bt; % Backward-compatible alias.
row.Fe3_bt = Fe3_bt;
row.Mn_bt = Mn_bt;
row.Ca_bt = Ca_bt;
row.Ti_bt = Ti_bt;
row.K_bt = K_bt;
row.Na_bt = Na_bt;

row.XMg_M = XMg_M;
row.XMg_B = XMg_B;
row.XAlVI_M_apfu = XAlVI_M_apfu;
row.XAlVI_B_apfu = XAlVI_B_apfu;
row.XAlVI_M = XAlVI_M;
row.XAlVI_B = XAlVI_B;
row.XTi_B = XTi_B;
row.XFe_B = XFe_B;
row.XMgB_minus_XAlVIB = XMgB_minus_XAlVIB;

row.Kideal_R1 = Kideal_R1;
row.K_exchange = Kideal_R1; % Backward-compatible alias.
row.lnK = lnK;
row.numerator = numerator;
row.denominator = denominator;

row.T_K = THoisch1989_K;
row.T_deg = THoisch1989_C;
row.THoisch1989_K = THoisch1989_K;
row.THoisch1989_C = THoisch1989_C;

row.is_P_in_calibration_range = is_P_in_calibration_range;
row.is_T_in_calibration_range = is_T_in_calibration_range;
row.is_Kideal_in_calibration_range = ...
    is_Kideal_in_calibration_range;
row.is_XMgB_minus_XAlVIB_in_calibration_range = ...
    is_XMgB_minus_XAlVIB_in_calibration_range;
row.is_XMg_B_in_calibration_range = ...
    is_XMg_B_in_calibration_range;
row.is_XAlVI_B_in_calibration_range = ...
    is_XAlVI_B_in_calibration_range;
row.is_XTi_B_in_calibration_range = ...
    is_XTi_B_in_calibration_range;
row.is_XFe_B_in_calibration_range = ...
    is_XFe_B_in_calibration_range;
row.is_XMg_M_in_calibration_range = ...
    is_XMg_M_in_calibration_range;

end

function printCompositionRangeWarnings(row, selectedCode_ms, selectedCode_bt)
% printCompositionRangeWarnings
% Print non-stopping warnings for finite composition indices outside the
% final-regression ranges in Hoisch (1989, Table 4, p. 571).

variableNames = { ...
    'Kideal_R1', ...
    'XMgB_minus_XAlVIB', ...
    'XMg_B', ...
    'XAlVI_B', ...
    'XTi_B', ...
    'XFe_B', ...
    'XMg_M'};

rangeMinimum = [0.108, 0.165, 0.319, 0.089, 0.0240, 0.339, 0.0085];
rangeMaximum = [0.449, 0.318, 0.417, 0.170, 0.0732, 0.443, 0.0431];
rangeLabels = { ...
    'Kideal_R1', ...
    'XMg_B - XAlVI_B', ...
    'XMg_B', ...
    'XAlVI_B', ...
    'XTi_B', ...
    'XFe_B', ...
    'XMg_M'};

for i = 1:numel(variableNames)
    values = row.(variableNames{i});
    finiteValuesMask = isfinite(values);
    outsideRange = finiteValuesMask & ...
        (values < rangeMinimum(i) | values > rangeMaximum(i));

    if any(outsideRange)
        finiteValues = values(finiteValuesMask);
        fprintf(2, ...
            ['WARNING: %s is outside the Hoisch (1989) Table 4 calibration ' ...
             'range %.4g–%.4g. %d of %d finite point(s) are outside; ' ...
             'finite range = %.4g–%.4g for %s & %s.\n'], ...
            rangeLabels{i}, ...
            rangeMinimum(i), ...
            rangeMaximum(i), ...
            sum(outsideRange), ...
            sum(finiteValuesMask), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_ms)), ...
            char(string(selectedCode_bt)));
    end
end

end

function mineral = prepareMineralRow(data_tbl, mineralLabel)
% prepareMineralRow
% Extract cation values from one selected mica row. Required and optional NaN
% values are retained. Missing optional variables are represented by NaN.

if height(data_tbl) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();

mineral.Mg = getRequiredVar(data_tbl, 'Mg_cation_apfu', mineralLabel);
mineral.Al = getRequiredVar(data_tbl, 'Al_cation_apfu', mineralLabel);
mineral.Si = getRequiredVar(data_tbl, 'Si_cation_apfu', mineralLabel);

mineral.Fe = getOptionalVar(data_tbl, 'Fe_cation_apfu', mineralLabel);
mineral.Fe3 = getOptionalVar(data_tbl, 'Fe3_cation_apfu', mineralLabel);
mineral.Mn = getOptionalVar(data_tbl, 'Mn_cation_apfu', mineralLabel);
mineral.Ca = getOptionalVar(data_tbl, 'Ca_cation_apfu', mineralLabel);
mineral.Ti = getOptionalVar(data_tbl, 'Ti_cation_apfu', mineralLabel);
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
