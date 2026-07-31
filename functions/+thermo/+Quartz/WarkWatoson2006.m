function results = WarkWatoson2006(rawdata_struct, P_kbar)
% functions/+thermo/+Quartz/WarkWatoson2006.m
% Tested with MATLAB R2024b
%
% Ti-in-Quartz geothermometer (TitaniQ)
% Wark, D.A. and Watson, E.B. (2006)
% Contributions to Mineralogy and Petrology, 152, 743-754
% DOI: https://doi.org/10.1007/s00410-006-0132-1
%
% NOTE ON FILE NAME
% The original function name supplied for this project is
% "WarkWatoson2006" (with "Watoson"). This spelling is retained in the
% function and file name for compatibility with existing launcher settings.
%
% -------------------------------------------------------------------------
% OVERVIEW
%
% This function interactively selects one Quartz analysis and calculates
% temperature using the original Wark and Watson (2006) Ti-in-quartz
% geothermometer (TitaniQ).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Quartz analysis, one output row is returned for every
% pressure supplied in P_kbar. The original Wark and Watson (2006) equation
% contains no pressure term, so all pressure rows for one selected Quartz
% analysis contain the same calculated temperature. This interface is
% compatible with both startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% TiO2 activity, aTiO2, is entered once through a dialog before Quartz
% selection. The same activity value is used for all selected analyses and
% all supplied pressure values in one function call.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Wark and Watson (2006) calibrated TitaniQ using 13 experiments in which
% quartz crystallized in equilibrium with rutile from aqueous fluid, HF
% fluid, or hydrous silicate melt. The direct experimental calibration is:
%
%   Temperature : 600-1000 degreeC
%   Pressure    : 1.0 GPa = 10 kbar only
%   Quartz Ti   : approximately 22-526 ppm elemental Ti by weight
%   TiO2 state  : rutile saturated, aTiO2 = 1
%
% These conditions are stated in the abstract on p. 743, listed in Table 1
% on p. 744, and described in the Experimental methods section on
% pp. 744-745. The calibration relation is Equation (6) on p. 747, and the
% temperature forms are Equations (7) and (8) on pp. 747 and 749.
%
% IMPORTANT APPLICATION LIMITATIONS
%
%   1) The thermometer directly records quartz crystallization or
%      re-equilibration temperature only when the measured Ti concentration
%      represents equilibrium Ti solubility at the stage of interest.
%      Post-crystallization diffusion, recrystallization, dissolution-
%      regrowth, or rutile exsolution can modify the preserved Ti
%      distribution (Ti diffusion discussion, pp. 749-750).
%
%   2) Use aTiO2 = 1 only when the analyzed Quartz domain crystallized or
%      re-equilibrated in equilibrium with rutile. Rutile elsewhere in the
%      rock does not by itself demonstrate that the selected Quartz domain
%      was rutile saturated at the stage being investigated
%      (pp. 748-750 and Summary, p. 753).
%
%   3) Rutile-absent systems can be treated only when aTiO2 is independently
%      constrained. In rutile-absent systems, 0 < aTiO2 < 1. If aTiO2 = 1
%      is assumed despite rutile absence, the calculated result is a minimum
%      equilibration temperature, not necessarily the true temperature
%      (pp. 748-750; Equation 8 on p. 749).
%
%   4) Errors in aTiO2 propagate directly into temperature. Wark and Watson
%      (2006) show that an activity uncertainty of approximately +/-0.1 can
%      produce errors of approximately +/-20 degreeC at lower temperatures
%      and approximately 30-35 degreeC near 1000 degreeC
%      (Fig. 6 and discussion on p. 749).
%
%   5) The equation has no pressure term and was directly calibrated only
%      at 10 kbar. The authors argued from natural and experimental
%      comparisons that the pressure effect appeared small
%      (Pressure effect discussion, p. 748), but the original calibration
%      did not experimentally establish pressure independence over a broad
%      pressure range. Results calculated at pressures other than 10 kbar
%      are therefore extrapolations of the original calibration.
%
%   6) Later polybaric Ti-in-quartz calibrations demonstrate measurable
%      pressure dependence. For applications over broad pressure ranges,
%      a pressure-dependent calibration such as Osborne et al. (2022)
%      should be considered in addition to the original Wark-Watson model.
%
%   7) Electron-microprobe analyses near Ti-bearing minerals can be
%      contaminated by secondary X-ray fluorescence. Wark and Watson (2006)
%      selected Quartz crystals with no visible rutile within approximately
%      100 micrometres. Rutile below the polished surface can also increase
%      apparent Ti concentrations (Analytical procedures and Fig. 3,
%      pp. 745-747; Summary, p. 753).
%
%   8) Quartz containing exsolved rutile needles requires special treatment.
%      A point analysis avoiding rutile records Ti remaining in the quartz
%      host, whereas a broad-beam or area analysis including exsolved rutile
%      may be used to approximate pre-exsolution bulk Ti. These approaches
%      record different quantities and must not be interpreted
%      interchangeably (UHT application, pp. 752-753).
%
%   9) Ti diffusion can modify zoning. Approximate characteristic transport
%      distances discussed by the authors are approximately 2 micrometres at
%      500 degreeC, 125 micrometres at 700 degreeC, and 2 mm at 900 degreeC
%      over 1 million years (pp. 749-750). Quartz cores and rims may
%      therefore record different stages of a thermal history.
%
%  10) The direct experimental temperature range is 600-1000 degreeC.
%      The authors suggested practical application down to approximately
%      400 degreeC when sub-ppm Ti can be measured by SIMS, but 400-600
%      degreeC is an extrapolation below the direct experimental
%      calibration. EPMA is generally practical from approximately
%      600 degreeC upward (abstract, p. 743; Applications and Summary,
%      pp. 750 and 753).
%
%  11) Reported model accuracy based only on regression parameters is
%      approximately +/-5 degreeC, but the uncertainty of a natural
%      application also depends on Ti analytical accuracy, aTiO2,
%      equilibrium interpretation, pressure, diffusion, and secondary
%      fluorescence (p. 747 and Fig. 5, p. 748).
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure differs from the direct calibration pressure of
%      10 kbar,
%   2) a finite calculated temperature lies outside 600-1000 degreeC,
%   3) equivalent finite elemental Ti lies outside approximately
%      22-526 ppm,
%   4) aTiO2 is greater than 1,
%   5) an explicitly stored calculation input is NaN,
%   6) the Ti conversion, activity, logarithm, denominator, or temperature
%      term is invalid, or
%   7) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
%
% rawdata_struct must contain:
%   rawdata_struct.Quartz : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% Required Quartz variable:
%   Ti_cation_apfu
%
% Ti_cation_apfu must represent Ti atoms per formula unit normalized to
% O = 2. It is not a wt% value and is not multiplied directly by 1e4.
%
% The published Wark-Watson equation requires elemental Ti concentration in
% ppm by weight. This implementation converts Ti_cation_apfu to equivalent
% elemental Ti ppm using a binary TiO2-SiO2 mass relation:
%
%   Ti_ppm =
%       [X_TiO2 * M_Ti] /
%       [X_TiO2*M_TiO2 + (1-X_TiO2)*M_SiO2] * 1e6
%
% where X_TiO2 is represented by Ti_cation_apfu on the O = 2 normalized
% quartz basis.
%
% Existing NaN values are retained and propagated; they are never replaced
% by zero. Finite Ti_cation_apfu and entered aTiO2 values must not be
% negative. Negative finite values and Inf stop the calculation. A value of
% zero is retained as an input, but the logarithm is undefined and the
% resulting temperature is returned as NaN with non-stopping diagnostics.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and from liquid NaN-input warnings is not
% applicable to this function.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Wark and Watson (2006), Equation (6), p. 747:
%
%   log10(Ti_ppm) = 5.69 - 3765 / T(K)
%
% For rutile-undersaturated systems, Equation (8), p. 749:
%
%                        3765
%   T(K) = ---------------------------------
%           5.69 - log10(Ti_ppm / aTiO2)
%
%   T(degreeC) = T(K) - 273.15
%
% where:
%   Ti_ppm = elemental Ti concentration in quartz, ppm by weight
%   aTiO2  = TiO2 activity relative to rutile saturation
%
% Pressure is accepted and stored for interface compatibility, but is not
% used explicitly in the original Wark and Watson (2006) equation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = WarkWatoson2006(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing a Quartz table
%   P_kbar         : finite non-negative numeric scalar or vector
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Quartz analysis. T_K, T_degreeC, and T_deg are supplied for
%             downstream launcher and plotting compatibility.
%

%% Input validation
if nargin < 2
    error('WarkWatoson2006 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve Quartz dataset
disp('=== Step 1: Preparing Quartz dataset ===');

if ~isfield(rawdata_struct, 'Quartz') || ~istable(rawdata_struct.Quartz)
    error('rawdata_struct must contain table: rawdata_struct.Quartz');
end
if isempty(rawdata_struct.Quartz)
    error('rawdata_struct.Quartz is empty.');
end

dataset_qtz = rawdata_struct.Quartz;

if ~ismember('Ti_cation_apfu', dataset_qtz.Properties.VariableNames)
    error('Quartz table must contain variable: Ti_cation_apfu');
end

disp('=== Preparing Quartz dataset has been finished ===');

%% 2) Initialize output container and screening settings
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct experimental calibration limits.
calibrationT_min_degreeC = 600;
calibrationT_max_degreeC = 1000;
calibrationP_kbar = 10;

% Approximate elemental Ti range in the 13 direct calibration experiments.
screeningTi_min_ppm = 22;
screeningTi_max_ppm = 526;

% The original calibration used one pressure only. Use a small numerical
% tolerance so values intended to be exactly 10 kbar are not misclassified
% because of floating-point representation.
pressureTolerance_kbar = 1e-9;
pressureOutsideCalibration = ...
    abs(P_kbar - calibrationP_kbar) > pressureTolerance_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3) Enter TiO2 activity
disp('=== Step 3: Entering TiO2 activity ===');

prompt = { ...
    ['Enter TiO2 activity (aTiO2):' newline ...
     'Use aTiO2 = 1 only when the selected Quartz domain' newline ...
     'crystallized or re-equilibrated in equilibrium with rutile.' newline ...
     'NaN is allowed and will be retained and reported.']};

dlgtitle = 'TiO2 activity input';
dims = [4 72];
definput = {'1'};

answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    error('TiO2 activity input canceled.');
end

aTiO2 = str2double(answer{1});

% NaN is deliberately allowed and propagated. Negative finite values and
% Inf are prohibited. Zero is retained and diagnosed after calculation.
if isinf(aTiO2) || (isfinite(aTiO2) && aTiO2 < 0)
    error('aTiO2 must not be negative or infinite. NaN and zero are allowed.');
end

disp(['TiO2 activity (aTiO2) = ' num2str(aTiO2)]);
disp(['NOTE: use aTiO2 = 1 only when rutile saturation and equilibrium ' ...
      'with the selected Quartz domain are independently justified.']);

if isfinite(aTiO2) && aTiO2 > 1
    fprintf(2, ...
        ['CAUTION: Entered aTiO2 = %.6g is greater than 1. ' ...
         'For the pure-rutile saturation reference, aTiO2 is normally ' ...
         'expected to be less than or equal to 1. The value has been ' ...
         'retained.\n'], ...
        aTiO2);
end

%% 4-5) Interactive Quartz selection and calculation
dataCodes_qtz = dataset_qtz{:, 1};
displayCodes_qtz = cellstr(string(dataCodes_qtz));

disp('=== Step 4: Selecting a data code from the list (Quartz) ===');

while true
    [selectedIdx_qtz, ok] = listdlg( ...
        'PromptString', ...
            'Please select the Quartz data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_qtz, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_qtz)
        disp('Selection canceled');
        break;
    end

    selectedCode_qtz = string(dataCodes_qtz(selectedIdx_qtz));
    selectedData_qtz = dataset_qtz(selectedIdx_qtz, :);
    disp(['Quartz selected: ' char(selectedCode_qtz)]);

    disp('=== Step 5: Calculating the temperature ===');

    % Report explicitly stored NaN inputs without replacing them by zero.
    nanInputNames = findNaNInputs(selectedData_qtz, aTiO2);

    % Stop only for negative finite values or Inf. NaN and zero are retained.
    validateNonNegativeInputs(selectedData_qtz, aTiO2);

    row = calcTemp(selectedData_qtz, P_kbar, aTiO2);

    nRows = height(row);
    row.dataCode_qtz = repmat(selectedCode_qtz, nRows, 1);
    row = movevars(row, 'dataCode_qtz', 'Before', 1);

    % Store one completed table block. Repeated enlargement of the complete
    % output table inside the interactive loop is avoided.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = ...
            [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % ----- Immediate result display -----
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    printTemperatureSummary(selectedCode_qtz, row.T_degreeC);

    % The direct experimental calibration was conducted only at 10 kbar.
    % Print the pressure warning only once for the complete function call.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure differs from the single direct ' ...
             'experimental calibration pressure of Wark and Watson ' ...
             '(2006): 10 kbar (1.0 GPa; abstract, p. 743; Experimental ' ...
             'methods, pp. 744-745). %d of %d pressure point(s) differ ' ...
             'from 10 kbar; input range = %.6g-%.6g kbar. The original ' ...
             'TitaniQ equation contains no pressure term, so these results ' ...
             'are pressure extrapolations and have been retained.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    printTemperatureCalibrationWarning( ...
        row.T_degreeC, ...
        calibrationT_min_degreeC, ...
        calibrationT_max_degreeC, ...
        selectedCode_qtz);

    printTiScreeningWarning( ...
        row.Ti_ppm, ...
        screeningTi_min_ppm, ...
        screeningTi_max_ppm, ...
        selectedCode_qtz);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Wark-Watson TitaniQ ' ...
             'calculation input(s) for %s: %s.\n' ...
             '         Existing NaN values were retained and propagated; ' ...
             'the calculated temperature may be NaN.\n'], ...
            char(selectedCode_qtz), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTerms = findInvalidEquationTerms(row);

    if ~isempty(invalidTerms)
        fprintf(2, ...
            ['WARNING: Invalid Wark-Watson TitaniQ composition, activity, ' ...
             'logarithm, denominator, or temperature term(s) were found ' ...
             'for %s: %s.\n' ...
             '         Affected values were retained as NaN, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(selectedCode_qtz), ...
            char(strjoin(invalidTerms, ', ')));
    end

    printNonFiniteTemperatureWarning(row.T_degreeC, selectedCode_qtz);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'WarkWatoson2006', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered result blocks only once after the selection loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'primaryTemperatureEquation', ...
        'Wark and Watson (2006) TitaniQ Equations (7) and (8)', ...
    'doi', ...
        'https://doi.org/10.1007/s00410-006-0132-1', ...
    'experimentalCalibrationTemperatureRange_degreeC', ...
        [calibrationT_min_degreeC, calibrationT_max_degreeC], ...
    'directExperimentalCalibrationPressure_kbar', ...
        calibrationP_kbar, ...
    'approximateCalibrationTiRange_ppm', ...
        [screeningTi_min_ppm, screeningTi_max_ppm], ...
    'pressureUsedInEquation', ...
        false, ...
    'aTiO2Required', ...
        true, ...
    'rutileSaturationActivity', ...
        1, ...
    'practicalLowerTemperatureWithSIMS_degreeC', ...
        400, ...
    'fileNameSpellingNote', ...
        'WarkWatoson spelling retained for project compatibility');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_qtz, aTiO2)
% findNaNInputs
% Return names of explicitly stored calculation inputs that are NaN. NaN is
% reported but does not stop the calculation.

nameBuffer = strings(2, 1);
nNames = 0;

value = data_qtz.Ti_cation_apfu;
validateScalarVariable(value, 'Quartz', 'Ti_cation_apfu');

if isnan(value)
    nNames = nNames + 1;
    nameBuffer(nNames) = "Quartz.Ti_cation_apfu";
end

if isnan(aTiO2)
    nNames = nNames + 1;
    nameBuffer(nNames) = "aTiO2";
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_qtz, aTiO2)
% validateNonNegativeInputs
% Stop when a stored Ti input or entered activity is negative or infinite.
% Zero and NaN are deliberately allowed and handled by non-stopping
% diagnostics after calculation.

value = data_qtz.Ti_cation_apfu;
validateScalarVariable(value, 'Quartz', 'Ti_cation_apfu');

invalidNames = strings(2, 1);
nInvalid = 0;

if isinf(value) || (isfinite(value) && value < 0)
    nInvalid = nInvalid + 1;
    invalidNames(nInvalid) = "Quartz.Ti_cation_apfu";
end

if isinf(aTiO2) || (isfinite(aTiO2) && aTiO2 < 0)
    nInvalid = nInvalid + 1;
    invalidNames(nInvalid) = "aTiO2";
end

if nInvalid > 0
    invalidNames = invalidNames(1:nInvalid);

    error([ ...
        'WarkWatoson2006: calculation inputs must not be negative or ' ...
        'infinite. Invalid value(s) were found in: ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_qtz, P_kbar, aTiO2)
% calcTemp
% Calculate Wark and Watson (2006) TitaniQ temperature for one selected
% Quartz row and repeat the pressure-independent result for every supplied
% pressure. Existing NaN values and invalid derived terms are retained as
% NaN.
%
% Ti_cation_apfu is interpreted as Ti atoms per formula unit normalized to
% O = 2 and converted to elemental Ti ppm by weight before applying the
% published TitaniQ equation.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

% --- Constants ---
M_Ti = 47.867;
M_TiO2 = 79.866;
M_SiO2 = 60.0843;

X_TiO2_qtz_scalar = data_qtz.Ti_cation_apfu;
validateScalarVariable( ...
    X_TiO2_qtz_scalar, 'Quartz', 'Ti_cation_apfu');

% Convert O = 2 normalized Ti apfu to equivalent elemental Ti ppm using a
% binary TiO2-SiO2 mass relation.
Ti_ppm_scalar = NaN;
Ti_conversion_valid_scalar = false;

if isnan(X_TiO2_qtz_scalar)
    Ti_conversion_valid_scalar = false;
elseif isfinite(X_TiO2_qtz_scalar) && ...
        X_TiO2_qtz_scalar >= 0 && X_TiO2_qtz_scalar <= 1

    formulaMass = ...
        X_TiO2_qtz_scalar .* M_TiO2 + ...
        (1 - X_TiO2_qtz_scalar) .* M_SiO2;

    if isfinite(formulaMass) && formulaMass > 0
        Ti_ppm_scalar = ...
            (X_TiO2_qtz_scalar .* M_Ti ./ formulaMass) .* 1e6;

        Ti_conversion_valid_scalar = ...
            isfinite(Ti_ppm_scalar) && Ti_ppm_scalar >= 0;
    end
end

% Repeat scalar Quartz/activity values for each supplied pressure.
X_TiO2_qtz = repmat(X_TiO2_qtz_scalar, nP, 1);
Ti_ppm = repmat(Ti_ppm_scalar, nP, 1);
aTiO2_vector = repmat(aTiO2, nP, 1);
Ti_conversion_valid = repmat(Ti_conversion_valid_scalar, nP, 1);

% Initialize all derived terms with NaN so that invalid inputs remain
% diagnosable and are never silently converted to zero.
Ti_over_aTiO2 = NaN(nP, 1);
log_term = NaN(nP, 1);
denominator = NaN(nP, 1);
T_raw_K = NaN(nP, 1);
T_raw_degreeC = NaN(nP, 1);
T_K = NaN(nP, 1);
T_degreeC = NaN(nP, 1);

validTi = ...
    isfinite(Ti_ppm) & ...
    Ti_ppm > 0;

validActivity = ...
    isfinite(aTiO2_vector) & ...
    aTiO2_vector > 0;

validRatioInputs = validTi & validActivity;

Ti_over_aTiO2(validRatioInputs) = ...
    Ti_ppm(validRatioInputs) ./ aTiO2_vector(validRatioInputs);

validRatio = ...
    validRatioInputs & ...
    isfinite(Ti_over_aTiO2) & ...
    Ti_over_aTiO2 > 0;

log_term(validRatio) = log10(Ti_over_aTiO2(validRatio));
denominator(validRatio) = 5.69 - log_term(validRatio);

validDenominator = ...
    validRatio & ...
    isfinite(denominator) & ...
    abs(denominator) > 1e-12;

T_raw_K(validDenominator) = ...
    3765 ./ denominator(validDenominator);
T_raw_degreeC(validDenominator) = ...
    T_raw_K(validDenominator) - 273.15;

% Non-positive Kelvin is physically invalid. Preserve the raw solution for
% diagnosis, but return NaN as the accepted temperature.
validTemperature = ...
    isfinite(T_raw_K) & T_raw_K > 0;

T_K(validTemperature) = T_raw_K(validTemperature);
T_degreeC(validTemperature) = T_raw_degreeC(validTemperature);

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PressureUsedInEquation = false(nP, 1);
row.aTiO2 = aTiO2_vector;
row.X_TiO2_qtz = X_TiO2_qtz;
row.Ti_cation_apfu = X_TiO2_qtz;
row.Ti_ppm = Ti_ppm;
row.Ti_conversion_valid = Ti_conversion_valid;
row.Ti_over_aTiO2 = Ti_over_aTiO2;
row.log_term = log_term;
row.denominator = denominator;
row.T_raw_K = T_raw_K;
row.T_raw_degreeC = T_raw_degreeC;
row.T_K = T_K;
row.T_degreeC = T_degreeC;
row.T_deg = T_degreeC;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid conversion, activity, logarithm, denominator, and
% temperature terms. Repeated pressure rows are summarized by term name.

termBuffer = strings(9, 1);
nTerms = 0;

if any(~row.Ti_conversion_valid)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Ti_cation_apfu-to-Ti_ppm conversion";
end

if any(~isfinite(row.X_TiO2_qtz) | ...
        row.X_TiO2_qtz <= 0 | row.X_TiO2_qtz > 1)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = ...
        "X_TiO2_qtz (must satisfy 0 < X <= 1)";
end

if any(~isfinite(row.Ti_ppm) | row.Ti_ppm <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Ti_ppm (> 0 required)";
end

if any(~isfinite(row.aTiO2) | row.aTiO2 <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "aTiO2 (> 0 required)";
end

if any(~isfinite(row.Ti_over_aTiO2) | row.Ti_over_aTiO2 <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Ti_ppm/aTiO2";
end

if any(~isfinite(row.log_term))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "log10(Ti_ppm/aTiO2)";
end

if any(~isfinite(row.denominator) | abs(row.denominator) <= 1e-12)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "TitaniQ denominator";
end

if any(~isfinite(row.T_raw_K) | row.T_raw_K <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "raw T_K";
end

if any(~isfinite(row.T_K) | row.T_K <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "accepted T_K";
end

invalidTerms = unique(termBuffer(1:nTerms), 'stable');

end

function printTemperatureSummary(selectedCode_qtz, temperatureValues)
% printTemperatureSummary
% Display one temperature or the first-to-last values for a pressure vector.

label = char(selectedCode_qtz);

if numel(temperatureValues) == 1
    disp([label ': ' num2str(temperatureValues) ' degreeC']);
else
    disp([label ': ' ...
        num2str(temperatureValues(1)) ' to ' ...
        num2str(temperatureValues(end)) ' degreeC']);
end

end

function printTemperatureCalibrationWarning( ...
        temperatureValues, minimumTemperature, maximumTemperature, ...
        selectedCode_qtz)
% printTemperatureCalibrationWarning
% Warn when finite temperatures lie outside the direct experimental
% calibration range. Results are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);

    fprintf(2, ...
        ['WARNING: Calculated Wark-Watson TitaniQ temperature is ' ...
         'outside the direct experimental calibration range of ' ...
         '%.4g-%.4g degreeC (Wark and Watson, 2006; abstract, p. 743; ' ...
         'Table 1, p. 744). %d of %d finite point(s) are outside; ' ...
         'calculated finite range = %.6g-%.6g degreeC for %s. ' ...
         'The result has been retained as an extrapolation.\n'], ...
        minimumTemperature, ...
        maximumTemperature, ...
        sum(outsideMask), ...
        sum(finiteMask), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(selectedCode_qtz));
end

end

function printTiScreeningWarning( ...
        TiValues, minimumTi, maximumTi, selectedCode_qtz)
% printTiScreeningWarning
% Warn when finite positive elemental Ti lies outside the approximate range
% of the direct calibration experiments. Results are retained.

finitePositiveMask = isfinite(TiValues) & TiValues > 0;
outsideMask = finitePositiveMask & ...
    (TiValues < minimumTi | TiValues > maximumTi);

if any(outsideMask)
    finiteValues = TiValues(finitePositiveMask);

    fprintf(2, ...
        ['CAUTION: Equivalent elemental Ti concentration is outside the ' ...
         'approximate %.6g-%.6g ppm range of the 13 direct calibration ' ...
         'experiments in Wark and Watson (2006, Table 1, p. 744). ' ...
         'Calculated finite input range = %.6g-%.6g ppm for %s. ' ...
         'This is a composition-screening caution, not a strict rejection. ' ...
         'The result has been retained.\n'], ...
        minimumTi, ...
        maximumTi, ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(selectedCode_qtz));
end

end

function printNonFiniteTemperatureWarning( ...
        temperatureValues, selectedCode_qtz)
% printNonFiniteTemperatureWarning
% Report NaN and Inf temperatures without modifying or deleting them.

invalidMask = ~isfinite(temperatureValues);

if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Wark-Watson TitaniQ temperature values ' ...
         'were calculated for %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the ' ...
         'calculation has not been stopped.\n'], ...
        char(selectedCode_qtz), ...
        sum(invalidMask), ...
        numel(temperatureValues), ...
        sum(isnan(temperatureValues)), ...
        sum(isinf(temperatureValues)));
end

end

function validateScalarVariable(value, mineralLabel, variableName)
% validateScalarVariable
% Require one numeric scalar from the selected table row. NaN is allowed;
% negative finite values and Inf are handled separately.

if ~isnumeric(value) || ~isscalar(value)
    error('%s variable %s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end
