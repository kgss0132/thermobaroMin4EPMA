function results = Ashchepkov2006(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet/Ashchepkov2006.m
% Tested with MATLAB R2024b
%
% Empirical single-garnet thermometer for mantle peridotites
% Ashchepkov, I.V. (2006)
% Russian Geology and Geophysics, 47, 1071–1085
% DOI: None
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Garnet analysis at a time and
% calculates temperatures using the empirical single-garnet thermometers
% proposed by Ashchepkov (2006).
%
% Implemented temperatures:
%   Eq. (1): Opx-thermometer-calibrated temperature
%   Eq. (2): Cpx-thermobarometry-calibrated temperature
%   Eq. (3): Gar-Cpx-thermometer-calibrated temperature
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Garnet analysis and appends results
% into a single output table.
%
% In the original paper, the thermometer equations are intended to be used
% iteratively together with garnet barometers. In this implementation, the
% user-supplied pressure P_kbar is used directly.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION DOMAIN AND APPLICATION NOTES
%
% Ashchepkov (2006) did not define a simple rectangular experimental
% calibration range in temperature and pressure. The thermometer equations
% were obtained by statistical regression against temperatures calculated
% with Opx-, Cpx-, and Gar-Cpx-based thermometers, using heterogeneous
% natural and experimental datasets for mantle peridotites, kimberlite
% xenoliths and concentrates, diamond inclusions, and related mantle
% assemblages. The data sources and calibration strategy are described on
% pp. 1071–1072; Eqs. (1)–(3) and their rotation corrections are presented
% on p. 1072.
%
% Important application notes from the original paper:
%   1) The method is primarily intended for Cr-bearing pyrope garnets from
%      mantle peridotites, kimberlite xenoliths/concentrates, and related
%      mantle assemblages, rather than general crustal metamorphic garnets
%      (pp. 1071–1072 and 1078–1082).
%   2) Results depend strongly on the analytical quality of low-abundance
%      TiO2 and Na2O, and Eq. (2) also depends on MnO. Oxide concentrations
%      in wt% are used directly rather than cation apfu (p. 1072).
%   3) Eq. (1) shows increased dispersion in the high-temperature part of
%      the comparison, whereas Eq. (2) shows increased dispersion in the
%      low-temperature part. Eq. (3) may agree better in some Cpx-bearing
%      assemblages (pp. 1072–1074).
%   4) Regression polynomials tend to smooth natural geotherms. The final
%      correction used here was tied empirically to the Udachnaya Opx
%      geotherm and is not an independent thermodynamic correction (p. 1076).
%   5) Multistage melt percolation, incomplete equilibration, reheating,
%      garnet zoning, or relict garnet may produce unreliable temperatures
%      and/or overestimated pressures. The preferred equation may differ
%      between depleted and nondepleted peridotites (pp. 1081–1082).
%   6) The original computer implementation iterated temperature and
%      pressure until successive temperatures differed by 1 degreeC and
%      rejected analyses with unacceptable KD values (pp. 1078–1079).
%
% Because no strict numerical T-P calibration limits were stated, the
% following values are used only as non-stopping empirical screening limits,
% based on the plotted comparison domains in Figs. 1–4 (pp. 1073–1079):
%
%   Temperature screening domain : approximately 400–1500 degreeC
%   Pressure screening domain    : approximately 0–80 kbar
%
% These are NOT formal experimental calibration limits. A warning outside
% these plotted domains means extrapolation beyond the empirical comparison
% field shown in the paper. In addition, 40–50 kbar is flagged separately
% because this interval showed especially large deviations and was excluded
% from part of the calibration of coupled garnet barometer Eq. (6) (p. 1077).
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) input pressure is outside the plotted 0–80 kbar comparison domain,
%   2) input pressure is exactly 0 kbar, where ln(KD)/P is undefined,
%   3) input pressure is within the 40–50 kbar caution interval,
%   4) a finite corrected temperature is outside 400–1500 degreeC,
%   5) a required thermometer input contains NaN or zero, or
%   6) a calculated corrected temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% Expected Garnet oxide variables (wt%):
%   SiO2, TiO2, Al2O3, Cr2O3, FeO, MnO, MgO, CaO, Na2O
%
% Common aliases such as TiO2_wt, TiO2_wtpct, and TiO2_wtpercent are
% accepted. FeOt and FeOt_wt are also accepted as FeO aliases.
%
% Missing optional variables that were previously replaced by zero are now
% represented by NaN. All NaN values are retained as missing values,
% propagated through calculations that use them, and reported by
% non-stopping warnings. NaN is never converted to zero.
%
% All finite oxide values must be greater than or equal to zero. Negative
% finite values and Inf values stop the calculation with an error. Zero is
% retained, but may produce KD = 0, division by zero, or a non-finite
% temperature; these outcomes are reported by fprintf warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% 1) Equation (1), calibrated against Opx thermometry
%     KD1 = MgO*TiO2 / ((CaO + MgO)^2 * FeO * Al2O3)
%     x1  = ln(KD1) / P_kbar
%
%     T1_raw(°C) = 5272.5*x1^3 + 10265*x1^2 + 6472*x1 + 2113
%     T1_rot(°C) = (T1_raw - 725)*1.20 + 923
%
% 2) Equation (2), calibrated against Cpx thermobarometry
%     KD2 = Na2O*MnO*TiO2 / ((CaO + MgO) * FeO * Al2O3)
%     x2  = ln(KD2) / P_kbar
%
%     T2_raw(°C) = 362.05*x2^3 + 1880.4*x2^2
%                  + 2659.6*x2 + 1695.5
%     T2_rot(°C) = (T2_raw - 800)*1.35 + 1023
%
% 3) Equation (3), calibrated against Gar-Cpx thermometry
%     The printed paper does not explicitly redefine KD for Eq. (3).
%     This implementation uses KD2, following the placement of Eq. (3)
%     immediately after the Cpx-based KD discussion.
%
%     T3_raw(°C) = 1700 + 3607*x2 + 3138*x2^2
%
% 4) Universal empirical correction applied to Eqs. (1)–(3)
%     Tfinal = T - 100*(T - 1500)/750 - 125
%              + 200*(1223 - T)/T + 2*TiO2
%
% IMPORTANT IMPLEMENTATION NOTES
% - Oxide concentrations are used in wt%, not cation apfu.
% - Pressure P is supplied and used in kbar.
% - The printed form of Eq. (2) is typographically ambiguous. This function
%   implements the polynomial in x = ln(KD)/P shown above.
% - Raw, rotation-corrected, and final corrected temperatures are retained
%   in the output for transparency and reproducibility.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Ashchepkov2006(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing a Garnet table (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet analysis. The output variable set is
%             intended to remain stable for downstream processing.
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values. Both fixed-pressure and pressure-range launchers
% are supported by accepting a numeric scalar or vector.
if nargin < 2
    error('Ashchepkov2006 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve oxide dataset
% Extract the required Garnet table from the input struct. The table is not
% modified here; aliases and missing-value handling are resolved only after
% the user selects a row.
disp('=== Step 1: Preparing oxide dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end

dataset_gt = rawdata_struct.Garnet;

disp('=== Preparing oxide dataset has been finished ===');

%% 2) Initialize output container
% Each result is stored temporarily as one table block. Repeated table
% concatenation inside the interactive loop is avoided because it forces
% MATLAB to repeatedly reallocate and copy the entire output table.
%
% The cell buffer is preallocated and doubled only when necessary. After the
% loop finishes, all blocks are concatenated once with vertcat.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Numerical screening limits based on the empirical comparison domains
% plotted in Ashchepkov (2006), not formal experimental calibration limits.
screeningT_min_degC = 400;
screeningT_max_degC = 1500;
screeningP_min_kbar = 0;
screeningP_max_kbar = 80;
pressureCaution_min_kbar = 40;
pressureCaution_max_kbar = 50;

% Pressure is common to every Garnet selection in this function call, so
% pressure-related messages are printed only once after the first result.
pressureOutsideScreening = ...
    P_kbar < screeningP_min_kbar | P_kbar > screeningP_max_kbar;
pressureAtZero = P_kbar == 0;
pressureInCautionInterval = ...
    P_kbar >= pressureCaution_min_kbar & ...
    P_kbar <= pressureCaution_max_kbar;
pressureMessagesIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–4) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
    % Assumption: the first table column stores the identifier displayed to
    % the user. Conversion to string makes numeric, categorical, cellstr,
    % and string identifiers usable in listdlg.
    dataCodes_gt = dataset_gt{:, 1};

    [selectedIdx_gt, ok] = listdlg( ...
        'PromptString', 'Please select the Garnet data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_gt)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_gt)
        disp('Selection canceled');
        break;
    end

    selectedCode_gt = dataCodes_gt(selectedIdx_gt);
    selectedCodeText = char(string(selectedCode_gt));
    disp(['Garnet selected: ' selectedCodeText]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperatures ===');

    selectedData_gt = dataset_gt(selectedIdx_gt, :);

    % Resolve accepted aliases once and retain NaN values. Optional missing
    % columns are represented by NaN rather than being silently replaced by
    % zero.
    garnet = prepareGarnetRow(selectedData_gt);

    % Check only variables used by at least one implemented thermometer.
    % NaN and zero do not stop calculation; they are reported after results.
    nanInputNames = findNaNInputs(garnet);
    zeroInputNames = findZeroInputs(garnet);

    % Negative finite values and Inf values are prohibited. NaN values are
    % deliberately allowed so that they propagate through the calculation.
    validateNonNegativeInputs(garnet);

    row = calcTemp(garnet, P_kbar);

    % Store the selected identifier for every pressure point.
    row.dataCode_garnet = ...
        repmat(string(selectedCode_gt), height(row), 1);
    row = movevars(row, {'dataCode_garnet'}, 'Before', 1);

    % Store this result as one table block. If the preallocated cell buffer
    % is full, double its capacity rather than growing by one element during
    % every loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo final corrected temperatures for immediate feedback. Raw and
    % rotation-corrected values remain available in the output table.
    disp('--------------------------------------------------');
    disp('=== Temperatures were calculated: ===');
    disp([selectedCodeText ': Eq1(final) = ' ...
        formatTemperatureRange(row.T1_final_deg) ' degreeC']);
    disp([selectedCodeText ': Eq2(final) = ' ...
        formatTemperatureRange(row.T2_final_deg) ' degreeC']);
    disp([selectedCodeText ': Eq3(final) = ' ...
        formatTemperatureRange(row.T3_final_deg) ' degreeC']);

    % Print pressure messages only once because P_kbar is shared by all
    % Garnet selections in this function call.
    if ~pressureMessagesIssued
        if any(pressureOutsideScreening)
            fprintf(2, ...
                ['WARNING: Input pressure is outside the approximate plotted comparison ' ...
                 'domain of Ashchepkov (2006): 0–80 kbar. This is not a formal ' ...
                 'experimental calibration limit, but values outside it are extrapolations. ' ...
                 '%d of %d pressure point(s) are outside; input range = %.4g–%.4g kbar.\n'], ...
                sum(pressureOutsideScreening), ...
                numel(P_kbar), ...
                min(P_kbar), ...
                max(P_kbar));
        end

        if any(pressureAtZero)
            fprintf(2, ...
                ['WARNING: Input pressure includes 0 kbar (%d of %d point(s)). ' ...
                 'The term ln(KD)/P is undefined at P = 0, so corresponding ' ...
                 'temperature results are expected to be NaN or Inf.\n'], ...
                sum(pressureAtZero), ...
                numel(P_kbar));
        end

        if any(pressureInCautionInterval)
            fprintf(2, ...
                ['WARNING: Input pressure includes the 40–50 kbar caution interval ' ...
                 'identified by Ashchepkov (2006). This interval showed especially ' ...
                 'large deviations and was excluded from part of the calibration of ' ...
                 'the coupled garnet barometer Eq. (6). %d of %d pressure point(s) ' ...
                 'fall in this interval.\n'], ...
                sum(pressureInCautionInterval), ...
                numel(P_kbar));
        end

        pressureMessagesIssued = true;
    end

    % Warn when finite final corrected temperatures lie outside the plotted
    % 400–1500 degreeC comparison domain. Each equation is checked separately
    % because the three calibrations have different scatter characteristics.
    printTemperatureScreeningWarning( ...
        row.T1_final_deg, 'Eq. (1)', selectedCodeText, ...
        screeningT_min_degC, screeningT_max_degC);
    printTemperatureScreeningWarning( ...
        row.T2_final_deg, 'Eq. (2)', selectedCodeText, ...
        screeningT_min_degC, screeningT_max_degC);
    printTemperatureScreeningWarning( ...
        row.T3_final_deg, 'Eq. (3)', selectedCodeText, ...
        screeningT_min_degC, screeningT_max_degC);

    % Print non-stopping input warnings immediately after the result. NaN is
    % retained and zero is not converted or rejected.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in thermometer input(s) for Garnet %s: %s.\n' ...
             '         NaN values were retained, and affected temperature results may remain NaN.\n'], ...
            selectedCodeText, ...
            char(strjoin(nanInputNames, ', ')));
    end

    if ~isempty(zeroInputNames)
        fprintf(2, ...
            ['WARNING: Zero was found in thermometer input(s) for Garnet %s: %s.\n' ...
             '         Zero was retained, but logarithms or divisions may produce NaN or Inf.\n'], ...
            selectedCodeText, ...
            char(strjoin(zeroInputNames, ', ')));
    end

    % Retain and report non-finite final outputs instead of stopping.
    printNonFiniteTemperatureWarning( ...
        row.T1_final_deg, 'Eq. (1)', selectedCodeText);
    printNonFiniteTemperatureWarning( ...
        row.T2_final_deg, 'Eq. (2)', selectedCodeText);
    printNonFiniteTemperatureWarning( ...
        row.T3_final_deg, 'Eq. (3)', selectedCodeText);

    disp('--------------------------------------------------');

    % Ask whether to repeat with another Garnet analysis.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Ashchepkov2006', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after all selections have been
% completed. Return an empty table when no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(garnet)
% findNaNInputs
% Return canonical names of thermometer inputs containing NaN. Variables
% not used in Eqs. (1)–(3), such as SiO2 and Cr2O3, are not included.

variableNames = {'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O'};
nameBuffer = strings(numel(variableNames), 1);
nNames = 0;

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    variableValue = garnet.(variableName);

    if any(isnan(variableValue(:)))
        nNames = nNames + 1;
        nameBuffer(nNames) = "Garnet." + string(variableName);
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function zeroInputNames = findZeroInputs(garnet)
% findZeroInputs
% Return canonical names of thermometer inputs containing finite zero.
% Zero is allowed by the input validator but can make KD or ln(KD) invalid.

variableNames = {'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O'};
nameBuffer = strings(numel(variableNames), 1);
nNames = 0;

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    variableValue = garnet.(variableName);

    if any(isfinite(variableValue(:)) & variableValue(:) == 0)
        nNames = nNames + 1;
        nameBuffer(nNames) = "Garnet." + string(variableName);
    end
end

zeroInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(garnet)
% validateNonNegativeInputs
% Stop calculation when any oxide value is negative or infinite. NaN and
% finite zero are intentionally allowed and handled by non-stopping warnings.

variableNames = fieldnames(garnet);
nameBuffer = strings(numel(variableNames), 1);
nNames = 0;

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    variableValue = garnet.(variableName);

    if ~isnumeric(variableValue) || ~isscalar(variableValue) || ...
            any(isinf(variableValue(:))) || ...
            any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nNames = nNames + 1;
        nameBuffer(nNames) = "Garnet." + string(variableName);
    end
end

invalidInputNames = nameBuffer(1:nNames);

if ~isempty(invalidInputNames)
    error(['Ashchepkov2006: oxide values must be numeric scalars that are ' ...
           'NaN or finite and non-negative. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(garnet, P_kbar)
% calcTemp
% Compute Ashchepkov (2006) temperatures for one Garnet composition at one
% or more user-supplied pressures.
%
% Inputs:
%   garnet : scalar struct containing Garnet oxide concentrations in wt%
%   P_kbar : pressure in kbar; finite non-negative column vector
%
% Output:
%   row : table containing one row per pressure value, including oxide
%         inputs, KD values, intermediate variables, and temperatures.

P = P_kbar(:);
nP = numel(P);

row = table();

% Replicate the selected Garnet composition once across all pressure points.
% NaN and zero values are retained exactly as supplied.
SiO2  = repmat(garnet.SiO2,  nP, 1);
TiO2  = repmat(garnet.TiO2,  nP, 1);
Al2O3 = repmat(garnet.Al2O3, nP, 1);
Cr2O3 = repmat(garnet.Cr2O3, nP, 1);
FeO   = repmat(garnet.FeO,   nP, 1);
MnO   = repmat(garnet.MnO,   nP, 1);
MgO   = repmat(garnet.MgO,   nP, 1);
CaO   = repmat(garnet.CaO,   nP, 1);
Na2O  = repmat(garnet.Na2O,  nP, 1);

% ---------------------------------------------------------------------
% Equation (1): Opx-thermometer-calibrated temperature
% ---------------------------------------------------------------------
den1 = (CaO + MgO).^2 .* FeO .* Al2O3;
num1 = MgO .* TiO2;
KD1 = num1 ./ den1;
x1 = log(KD1) ./ P;

T1_raw_deg = ...
    5272.5 .* x1.^3 + ...
    10265  .* x1.^2 + ...
    6472   .* x1 + ...
    2113;

% Rotation correction described on p. 1072.
T1_rot_deg = (T1_raw_deg - 725) .* 1.20 + 923;

% Universal empirical geotherm correction described on p. 1076.
T1_final_deg = applyUniversalTempCorrection(T1_rot_deg, TiO2);
isValidEq1 = isfinite(T1_final_deg);

% ---------------------------------------------------------------------
% Equation (2): Cpx-thermobarometry-calibrated temperature
% ---------------------------------------------------------------------
den2 = (CaO + MgO) .* FeO .* Al2O3;
num2 = Na2O .* MnO .* TiO2;
KD2 = num2 ./ den2;
x2 = log(KD2) ./ P;

T2_raw_deg = ...
    362.05 .* x2.^3 + ...
    1880.4 .* x2.^2 + ...
    2659.6 .* x2 + ...
    1695.5;

% Rotation correction described on p. 1072.
T2_rot_deg = (T2_raw_deg - 800) .* 1.35 + 1023;

% Universal empirical geotherm correction described on p. 1076.
T2_final_deg = applyUniversalTempCorrection(T2_rot_deg, TiO2);
isValidEq2 = isfinite(T2_final_deg);

% ---------------------------------------------------------------------
% Equation (3): Gar-Cpx-thermometer-calibrated temperature
% ---------------------------------------------------------------------
x3 = x2;
T3_raw_deg = ...
    1700 + ...
    3607 .* x3 + ...
    3138 .* x3.^2;

T3_final_deg = applyUniversalTempCorrection(T3_raw_deg, TiO2);
isValidEq3 = isfinite(T3_final_deg);

% --- Pack outputs ---
% Every output variable has nP rows, allowing the same function to support
% both fixed-pressure and pressure-range launchers.
row.P_kbar = P;

row.SiO2_gt = SiO2;
row.TiO2_gt = TiO2;
row.Al2O3_gt = Al2O3;
row.Cr2O3_gt = Cr2O3;
row.FeO_gt = FeO;
row.MnO_gt = MnO;
row.MgO_gt = MgO;
row.CaO_gt = CaO;
row.Na2O_gt = Na2O;

row.KD1 = KD1;
row.lnKD1_over_P = x1;
row.T1_raw_deg = T1_raw_deg;
row.T1_rot_deg = T1_rot_deg;
row.T1_final_deg = T1_final_deg;
row.isValidEq1 = isValidEq1;

row.KD2 = KD2;
row.lnKD2_over_P = x2;
row.T2_raw_deg = T2_raw_deg;
row.T2_rot_deg = T2_rot_deg;
row.T2_final_deg = T2_final_deg;
row.isValidEq2 = isValidEq2;

row.T3_raw_deg = T3_raw_deg;
row.T3_final_deg = T3_final_deg;
row.isValidEq3 = isValidEq3;

end

function Tcorr = applyUniversalTempCorrection(T_in, TiO2_wt)
% applyUniversalTempCorrection
% Apply the empirical correction described by Ashchepkov (2006, p. 1076):
%
%   Tcorr = T - 100*(T - 1500)/750 - 125
%           + 200*(1223 - T)/T + 2*TiO2
%
% T and Tcorr are in degreeC; TiO2 is in wt%. The calculation is fully
% vectorized. NaN propagates naturally. Non-finite T inputs and T = 0 are
% returned as NaN rather than stopping the calculation.

Tcorr = ...
    T_in ...
    - 100 .* (T_in - 1500) ./ 750 ...
    - 125 ...
    + 200 .* (1223 - T_in) ./ T_in ...
    + 2 .* TiO2_wt;

invalidTemperature = ~isfinite(T_in) | T_in == 0;
Tcorr(invalidTemperature) = NaN;

end

function garnet = prepareGarnetRow(data_gt)
% prepareGarnetRow
% Extract one selected Garnet row and resolve accepted oxide aliases. Values
% are not converted to zero. Optional missing variables are returned as NaN.

if height(data_gt) ~= 1
    error('Garnet input must be a 1-row table.');
end

garnet = struct();

garnet.SiO2 = getOxideWtOrError( ...
    data_gt, {'SiO2','SiO2_wt','SiO2_wtpct','SiO2_wtpercent'}, 'Garnet');
garnet.TiO2 = getOxideWtOrNaN( ...
    data_gt, {'TiO2','TiO2_wt','TiO2_wtpct','TiO2_wtpercent'});
garnet.Al2O3 = getOxideWtOrError( ...
    data_gt, {'Al2O3','Al2O3_wt','Al2O3_wtpct','Al2O3_wtpercent'}, 'Garnet');
garnet.Cr2O3 = getOxideWtOrNaN( ...
    data_gt, {'Cr2O3','Cr2O3_wt','Cr2O3_wtpct','Cr2O3_wtpercent'});
garnet.FeO = getOxideWtOrError( ...
    data_gt, {'FeO','FeO_wt','FeO_wtpct','FeO_wtpercent','FeOt','FeOt_wt'}, 'Garnet');
garnet.MnO = getOxideWtOrNaN( ...
    data_gt, {'MnO','MnO_wt','MnO_wtpct','MnO_wtpercent'});
garnet.MgO = getOxideWtOrError( ...
    data_gt, {'MgO','MgO_wt','MgO_wtpct','MgO_wtpercent'}, 'Garnet');
garnet.CaO = getOxideWtOrError( ...
    data_gt, {'CaO','CaO_wt','CaO_wtpct','CaO_wtpercent'}, 'Garnet');
garnet.Na2O = getOxideWtOrNaN( ...
    data_gt, {'Na2O','Na2O_wt','Na2O_wtpct','Na2O_wtpercent'});

end

function value = getOxideWtOrError(tbl, candidateNames, mineralLabel)
% getOxideWtOrError
% Return the first matching oxide variable. Missing required variables stop
% calculation. NaN is allowed and retained; Inf and negative values are
% checked later by validateNonNegativeInputs.

varName = findFirstExistingVar(tbl, candidateNames);

if isempty(varName)
    error('%s table must contain one of these variables: %s', ...
        mineralLabel, strjoin(candidateNames, ', '));
end

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s must be a numeric scalar in a 1-row table.', varName);
end

end

function value = getOxideWtOrNaN(tbl, candidateNames)
% getOxideWtOrNaN
% Return the first matching optional oxide variable. If no accepted alias is
% present, return NaN. This replaces the previous zero-substitution behavior.

varName = findFirstExistingVar(tbl, candidateNames);

if isempty(varName)
    value = NaN;
    return;
end

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s must be a numeric scalar in a 1-row table.', varName);
end

end

function varName = findFirstExistingVar(tbl, candidateNames)
% findFirstExistingVar
% Return the first accepted variable name present in the input table.

varName = '';

for i = 1:numel(candidateNames)
    if ismember(candidateNames{i}, tbl.Properties.VariableNames)
        varName = candidateNames{i};
        return;
    end
end

end

function txt = formatTemperatureRange(values)
% formatTemperatureRange
% Format scalar or vector temperatures for the command-window summary.

if isempty(values)
    txt = 'NaN';
elseif isscalar(values)
    txt = numberOrNaN(values(1));
else
    txt = [numberOrNaN(values(1)) ' to ' numberOrNaN(values(end))];
end

end

function txt = numberOrNaN(value)
% numberOrNaN
% Convert one numeric value to text while preserving NaN/Inf labels.

if isnan(value)
    txt = 'NaN';
elseif isinf(value)
    if value > 0
        txt = 'Inf';
    else
        txt = '-Inf';
    end
else
    txt = num2str(value);
end

end

function printTemperatureScreeningWarning( ...
        temperatures, equationLabel, dataCode, Tmin, Tmax)
% printTemperatureScreeningWarning
% Print a non-stopping warning when finite corrected temperatures lie
% outside the empirical comparison domain plotted by Ashchepkov (2006).

finiteTemperature = isfinite(temperatures);
outsideScreening = finiteTemperature & ...
    (temperatures < Tmin | temperatures > Tmax);

if any(outsideScreening)
    finiteValues = temperatures(finiteTemperature);

    fprintf(2, ...
        ['WARNING: %s corrected temperature is outside the approximate plotted ' ...
         'comparison domain of Ashchepkov (2006): %.4g–%.4g degreeC. ' ...
         'This is not a formal experimental calibration limit. %d of %d finite ' ...
         'temperature point(s) are outside; calculated finite range = ' ...
         '%.4g–%.4g degreeC for Garnet %s.\n'], ...
        equationLabel, ...
        Tmin, ...
        Tmax, ...
        sum(outsideScreening), ...
        sum(finiteTemperature), ...
        min(finiteValues), ...
        max(finiteValues), ...
        dataCode);
end

end

function printNonFiniteTemperatureWarning(temperatures, equationLabel, dataCode)
% printNonFiniteTemperatureWarning
% Report retained NaN/Inf corrected temperatures without stopping.

invalidTemperature = ~isfinite(temperatures);

if any(invalidTemperature)
    fprintf(2, ...
        ['WARNING: Non-finite %s corrected temperature values were calculated ' ...
         'for Garnet %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and calculation has not been stopped.\n'], ...
        equationLabel, ...
        dataCode, ...
        sum(invalidTemperature), ...
        numel(temperatures), ...
        sum(isnan(temperatures)), ...
        sum(isinf(temperatures)));
end

end
