function results = Zhang2020(rawdata_struct, P_kbar)
% functions/+thermo/+Quartz/Zhang2020.m
% Tested with MATLAB R2024b
%
% Ti-in-Quartz thermobarometer / geothermometer
% Zhang, C., Li, X., Almeev, R.R., Horn, I., Behrens, H., and Holtz, F. (2020)
% Earth and Planetary Science Letters, 538, 116213
% DOI: https://doi.org/10.1016/j.epsl.2020.116213
%
% -------------------------------------------------------------------------
% OVERVIEW
%
% This function interactively selects one Quartz analysis and calculates
% temperature using the pressure-dependent Ti-in-quartz model of
% Zhang et al. (2020), Equation (5).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Quartz analysis, one output row is returned for every
% pressure supplied in P_kbar. This interface is compatible with both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% TiO2 activity relative to rutile saturation, aTiO2, is entered once
% through a dialog before Quartz selection. The same activity value is used
% for all selected analyses and all supplied pressure values in one
% function call.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Zhang et al. (2020) calibrated the Ti-in-quartz relation by crystallizing
% quartz from hydrous high-silica aluminosilicate melts at crustal
% pressures. The direct experimental calibration conditions are:
%
%   Temperature : 700-900 degreeC
%   Pressure    : 0.5-4 kbar
%   Quartz Ti   : approximately 97-694 ppm elemental Ti for the
%                 17 rutile-saturated experiments used to define the
%                 rutile-saturated Ti-in-quartz model
%   TiO2 state  : rutile saturated for the primary calibration,
%                 corresponding to aTiO2 = 1
%   Melt type   : hydrous high-silica, rhyolitic, strongly peraluminous
%                 aluminosilicate melt
%
% The overall calibration range is stated in the abstract on p. 1.
% Starting materials and intended melt compositions are described on p. 2.
% Experimental conditions and measured Quartz Ti concentrations are listed
% in Table 1 on p. 3. Analytical procedures are described on pp. 3-4.
% Attainment of near-equilibrium conditions is discussed on pp. 5-7.
% The rutile-saturated regression is Equation (3) on p. 7, and the
% activity-corrected form implemented here is Equation (5) on p. 9.
%
% IMPORTANT APPLICATION LIMITATIONS
%
%   1) The model was calibrated for Quartz crystallized from hydrous
%      high-silica melt, principally to investigate storage conditions of
%      rhyolitic and other silicic magmas at shallow crustal depths
%      (abstract and Introduction, pp. 1-2). It was not calibrated as a
%      general thermometer for all metamorphic, hydrothermal, mafic, or
%      ultramafic Quartz-bearing systems.
%
%   2) The direct pressure range is only 0.5-4 kbar. Although Equation (5)
%      contains a non-linear P^0.2 term, the exponent and coefficients were
%      fitted using low-pressure experiments. Calculations below 0.5 kbar
%      or above 4 kbar are extrapolations of the direct calibration
%      (Sections 4.2-4.3, pp. 6-8).
%
%   3) Use aTiO2 = 1 only when the analyzed Quartz domain crystallized or
%      re-equilibrated at rutile saturation. The occurrence of rutile
%      elsewhere in a sample does not by itself demonstrate equilibrium
%      with the selected Quartz domain at the stage of interest
%      (Sections 4.4-4.5, pp. 8-10).
%
%   4) Rutile-undersaturated systems can be treated using Equation (5) only
%      when aTiO2 is independently constrained. The activity correction
%      assumes ideal TiO2 activity behavior in the silicate melt
%      (activity coefficient approximately 1). If this assumption is not
%      appropriate, calculated temperatures may be inaccurate
%      (Section 4.5, p. 9).
%
%   5) The analyzed Quartz Ti must represent dissolved Ti in the Quartz
%      lattice. Zhang et al. (2020) excluded experiments using the
%      highest-Ti starting glasses because the synthesized Quartz contained
%      rutile inclusions that prevented reliable measurement of dissolved
%      Quartz Ti (Starting materials, p. 2).
%
%   6) Natural Quartz analyses containing rutile inclusions, rutile
%      lamellae, Ti-oxide contamination, glass inclusions, or mixed
%      analytical volumes can yield artificially high Quartz Ti and
%      overestimated temperatures. BSE, CL, or time-resolved analytical
%      signals should be used to verify that the measured signal represents
%      the intended Quartz domain (pp. 3-5).
%
%   7) The experimental Quartz was compositionally homogeneous in CL
%      images and the authors used multiple textural, compositional, and
%      activity comparisons to argue for near-equilibrium Ti partitioning
%      between Quartz and melt (Sections 3.1 and 4.1, pp. 4-7). Natural
%      applications require equivalent evidence that the analyzed Quartz
%      domain records equilibrium at the temperature of interest.
%
%   8) Analytical precision is important. The reported Ti detection limits
%      were approximately 38 ppm by EPMA and approximately 10 ppm by
%      fs-LA-ICP-MS (Analytical methods, pp. 3-4). Low-Ti natural Quartz
%      near an analytical detection limit can produce large relative
%      uncertainty in calculated temperature.
%
%   9) The experimental starting melts were high-silica, Qz-Ab-Or-rich,
%      strongly peraluminous compositions with an aluminum saturation index
%      of approximately 1.6-1.8 (Starting materials, p. 2). Application to
%      substantially different melt compositions should be regarded as an
%      extrapolation even when P and T lie inside the numerical calibration
%      range.
%
%  10) Ti concentration alone cannot determine both pressure and
%      temperature. This function calculates temperature from independently
%      supplied pressure and aTiO2. If pressure is unknown, another
%      independent temperature or equilibrium constraint is required
%      (Section 4.5, pp. 9-10).
%
%  11) The approximate 97-694 ppm range is the envelope of the
%      rutile-saturated experimental Quartz data used for the direct model.
%      It is not an independent rectangular compositional calibration
%      because equilibrium Quartz Ti changes with both pressure and
%      temperature. Values outside this envelope are therefore flagged as
%      cautions rather than rejected automatically.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside 0.5-4 kbar,
%   2) a finite calculated temperature is outside 700-900 degreeC,
%   3) equivalent finite elemental Ti is outside approximately
%      97-694 ppm,
%   4) aTiO2 is greater than 1,
%   5) an explicitly stored calculation input is NaN,
%   6) the Ti conversion, activity, logarithm, pressure term, denominator,
%      or temperature term is invalid, or
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
% The published Zhang et al. (2020) equation requires elemental Ti
% concentration in Quartz in ppm by weight. This implementation converts
% Ti_cation_apfu to equivalent elemental Ti ppm using the binary
% TiO2-SiO2 mass relation:
%
%   Ti_ppm =
%       [X_TiO2 * M_Ti] /
%       [X_TiO2*M_TiO2 + (1-X_TiO2)*M_SiO2] * 1e6
%
% where X_TiO2 is represented by Ti_cation_apfu on the O = 2 normalized
% Quartz basis.
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
% Zhang et al. (2020), Equation (5), p. 9:
%
%   log10(C_Ti_Qtz / aTiO2) =
%       5.3226 - 1948.4/T - 981.4*P^0.2/T
%
% Rearranged to calculate temperature:
%
%          1948.4 + 981.4*P^0.2
%   T(K) = ----------------------
%          5.3226 - log10(C_Ti_Qtz/aTiO2)
%
%   T(degreeC) = T(K) - 273.15
%
% where:
%   C_Ti_Qtz : elemental Ti concentration in Quartz, ppm by weight
%   aTiO2    : TiO2 activity relative to rutile saturation
%   P        : pressure in kbar
%   T        : temperature in K
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Zhang2020(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing a Quartz table
%   P_kbar         : finite non-negative numeric scalar or vector
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Quartz analysis. T_K, T_degreeC, T_degC, and T_deg are
%             supplied for downstream launcher and plotting compatibility.
%

%% Input validation
if nargin < 2
    error('Zhang2020 requires (rawdata_struct, P_kbar).');
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

%% 2) Initialize output container and calibration screening settings
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct experimental calibration limits.
calibrationT_min_degreeC = 700;
calibrationT_max_degreeC = 900;
calibrationP_min_kbar = 0.5;
calibrationP_max_kbar = 4;

% Approximate elemental Ti envelope of the 17 rutile-saturated experiments.
screeningTi_min_ppm = 97;
screeningTi_max_ppm = 694;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3) Enter TiO2 activity
disp('=== Step 3: Entering TiO2 activity ===');

prompt = { ...
    ['Enter TiO2 activity relative to rutile saturation (aTiO2):' newline ...
     'Use aTiO2 = 1 only when the selected Quartz domain' newline ...
     'crystallized or re-equilibrated at rutile saturation.' newline ...
     'NaN is allowed and will be retained and reported.']};

dlgtitle = 'TiO2 activity input';
dims = [4 74];
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

if numel(P_kbar) == 1
    disp(['Pressure P = ' num2str(P_kbar) ' kbar']);
else
    disp(['Pressure P = ' num2str(P_kbar(1)) ' to ' ...
        num2str(P_kbar(end)) ' kbar']);
end

disp(['TiO2 activity (aTiO2) = ' num2str(aTiO2)]);
disp(['NOTE: use aTiO2 = 1 only when rutile saturation and equilibrium ' ...
      'with the selected Quartz domain are independently justified.']);

if isfinite(aTiO2) && aTiO2 > 1
    fprintf(2, ...
        ['CAUTION: Entered aTiO2 = %.6g is greater than 1. ' ...
         'For activity relative to rutile saturation, aTiO2 is normally ' ...
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

    % Pressure is common to every selected Quartz analysis in this function
    % call, so the pressure warning is printed only once.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the direct experimental ' ...
             'calibration range of Zhang et al. (2020): 0.5-4 kbar ' ...
             '(abstract, p. 1; Table 1, p. 3). %d of %d pressure ' ...
             'point(s) are outside; input range = %.6g-%.6g kbar. ' ...
             'The supplied pressures and calculated results have been ' ...
             'retained as extrapolations.\n'], ...
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
            ['WARNING: NaN was found in the Zhang et al. (2020) ' ...
             'calculation input(s) for %s: %s.\n' ...
             '         Existing NaN values were retained and propagated; ' ...
             'the calculated temperature may be NaN.\n'], ...
            char(selectedCode_qtz), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTerms = findInvalidEquationTerms(row);

    if ~isempty(invalidTerms)
        fprintf(2, ...
            ['WARNING: Invalid Zhang et al. (2020) composition, activity, ' ...
             'pressure, logarithm, denominator, or temperature term(s) ' ...
             'were found for %s: %s.\n' ...
             '         Affected values were retained as NaN, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(selectedCode_qtz), ...
            char(strjoin(invalidTerms, ', ')));
    end

    printNonFiniteTemperatureWarning(row.T_degreeC, selectedCode_qtz);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Zhang2020', ...
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
        'Zhang et al. (2020) Equation (5)', ...
    'doi', ...
        'https://doi.org/10.1016/j.epsl.2020.116213', ...
    'experimentalCalibrationTemperatureRange_degreeC', ...
        [calibrationT_min_degreeC, calibrationT_max_degreeC], ...
    'experimentalCalibrationPressureRange_kbar', ...
        [calibrationP_min_kbar, calibrationP_max_kbar], ...
    'approximateRutileSaturatedQuartzTiRange_ppm', ...
        [screeningTi_min_ppm, screeningTi_max_ppm], ...
    'pressureUsedInEquation', ...
        true, ...
    'aTiO2Required', ...
        true, ...
    'rutileSaturationActivity', ...
        1, ...
    'calibrationSystem', ...
        'hydrous high-silica rhyolitic peraluminous melt', ...
    'liquidCompositionUsed', ...
        false);

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
        'Zhang2020: calculation inputs must not be negative or ' ...
        'infinite. Invalid value(s) were found in: ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_qtz, P_kbar, aTiO2)
% calcTemp
% Calculate Zhang et al. (2020) Equation (5) temperature for one selected
% Quartz row and every supplied pressure. Existing NaN values and invalid
% derived terms are retained as NaN.
%
% Ti_cation_apfu is interpreted as Ti atoms per formula unit normalized to
% O = 2 and converted to elemental Ti ppm by weight before applying the
% published equation.

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

% Repeat scalar Quartz/activity values for each pressure.
X_TiO2_qtz = repmat(X_TiO2_qtz_scalar, nP, 1);
Ti_ppm = repmat(Ti_ppm_scalar, nP, 1);
aTiO2_vector = repmat(aTiO2, nP, 1);
Ti_conversion_valid = repmat(Ti_conversion_valid_scalar, nP, 1);

% Initialize all derived terms with NaN so invalid inputs remain
% diagnosable and are never silently converted to zero.
Ti_over_aTiO2 = NaN(nP, 1);
log_term = NaN(nP, 1);
P_power_0p2 = NaN(nP, 1);
numerator = NaN(nP, 1);
denominator = NaN(nP, 1);
T_raw_K = NaN(nP, 1);
T_raw_degreeC = NaN(nP, 1);
T_K = NaN(nP, 1);
T_degreeC = NaN(nP, 1);

% P_kbar was validated as finite and non-negative before entering calcTemp.
P_power_0p2(:) = P_kbar .^ 0.2;
numerator(:) = 1948.4 + 981.4 .* P_power_0p2;

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
denominator(validRatio) = 5.3226 - log_term(validRatio);

validDenominator = ...
    validRatio & ...
    isfinite(P_power_0p2) & ...
    isfinite(numerator) & ...
    isfinite(denominator) & ...
    abs(denominator) > 1e-12;

T_raw_K(validDenominator) = ...
    numerator(validDenominator) ./ denominator(validDenominator);
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
row.aTiO2 = aTiO2_vector;
row.X_TiO2_qtz = X_TiO2_qtz;
row.Ti_cation_apfu = X_TiO2_qtz;
row.Ti_ppm = Ti_ppm;
row.Ti_conversion_valid = Ti_conversion_valid;
row.Ti_over_aTiO2 = Ti_over_aTiO2;
row.log_term = log_term;
row.P_power_0p2 = P_power_0p2;
row.numerator = numerator;
row.denominator = denominator;
row.T_raw_K = T_raw_K;
row.T_raw_degreeC = T_raw_degreeC;
row.T_K = T_K;
row.T_degreeC = T_degreeC;
row.T_degC = T_degreeC;
row.T_deg = T_degreeC;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid conversion, activity, pressure, logarithm, denominator,
% and temperature terms. Repeated pressure rows are summarized by term name.

termBuffer = strings(10, 1);
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

if any(~isfinite(row.P_power_0p2))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "P_kbar^0.2";
end

if any(~isfinite(row.numerator))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Equation (5) numerator";
end

if any(~isfinite(row.denominator) | abs(row.denominator) <= 1e-12)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Equation (5) denominator";
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
        ['WARNING: Calculated Zhang et al. (2020) temperature is outside ' ...
         'the direct experimental calibration range of %.4g-%.4g degreeC ' ...
         '(abstract, p. 1; Table 1, p. 3). %d of %d finite point(s) ' ...
         'are outside; calculated finite range = %.6g-%.6g degreeC ' ...
         'for %s. The result has been retained as an extrapolation.\n'], ...
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
% of the rutile-saturated experiments used for the direct model.

finitePositiveMask = isfinite(TiValues) & TiValues > 0;
outsideMask = finitePositiveMask & ...
    (TiValues < minimumTi | TiValues > maximumTi);

if any(outsideMask)
    finiteValues = TiValues(finitePositiveMask);

    fprintf(2, ...
        ['CAUTION: Equivalent elemental Ti concentration is outside the ' ...
         'approximate %.6g-%.6g ppm envelope of the 17 rutile-saturated ' ...
         'experiments used for the Zhang et al. (2020) model ' ...
         '(Table 1, p. 3). Calculated finite input range = ' ...
         '%.6g-%.6g ppm for %s. This is a composition-screening caution, ' ...
         'not a strict rejection. The result has been retained.\n'], ...
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
        ['WARNING: Non-finite Zhang et al. (2020) temperature values ' ...
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
