function results = Osborne2022(rawdata_struct, P_kbar)
% functions/+thermo/+Quartz/Osborne2022.m
% Tested with MATLAB R2024b
%
% Expanded-volume (Vex) Ti-in-Quartz thermobarometer / geothermometer
% Osborne, Z.R., Thomas, J.B., Nachlas, W.O., Angel, R.J.,
% Hoff, C.M., and Watson, E.B. (2022)
% Contributions to Mineralogy and Petrology, 177, Article 31, 1-21
% DOI: https://doi.org/10.1007/s00410-022-01896-8
%
% -------------------------------------------------------------------------
% OVERVIEW
%
% This function interactively selects one Quartz analysis and calculates
% temperature using the preferred expanded-volume (Vex) Ti-in-quartz
% solubility model of Osborne et al. (2022).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Quartz analysis, one output row is returned for every
% pressure supplied in P_kbar. This interface is compatible with both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
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
% Osborne et al. (2022) developed the preferred Vex model from a collective
% dataset of 89 experiments in which quartz co-crystallized with rutile.
% The experiments span much of the alpha- and beta-quartz stability fields:
%
%   Temperature : 550-1050 degreeC
%   Pressure    : 2-30 kbar
%   Quartz field: 57 alpha-quartz and 32 beta-quartz experiments
%   Quartz Ti   : approximately 8-1971 ppm elemental Ti across the
%                 collective rutile-saturated calibration dataset
%
% The overall P-T range is stated in the abstract and summarized in the
% Discussion and Conclusions (pp. 1, 4, and 18). Experimental conditions
% and Ti concentrations are listed in Table 1 on pp. 5-6. The preferred
% Vex model and fit parameters are presented in Equation (7) and Table 3 on
% pp. 9-11. The temperature form used here is Equation (8) on p. 11.
%
% IMPORTANT APPLICATION LIMITATIONS
%
%   1) The model describes equilibrium Ti solubility in quartz. The analyzed
%      Quartz domain must have co-crystallized or re-equilibrated with the
%      minerals and medium used to constrain pressure and TiO2 activity.
%      Quartz and associated minerals must not have been substantially
%      modified by later diffusion, exsolution, recrystallization, or other
%      post-crystallization processes (application discussion, pp. 11-14).
%
%   2) Use aTiO2 = 1 only when the analyzed Quartz domain crystallized or
%      re-equilibrated at rutile saturation. The occurrence of rutile
%      elsewhere in the rock does not by itself demonstrate equilibrium
%      with the selected Quartz domain. The activity is referenced to pure
%      rutile at the standard state of 1 bar and 25 degreeC (pp. 7 and
%      11-14).
%
%   3) Rutile-free systems can be treated only when aTiO2 is independently
%      constrained. Osborne et al. (2022) tested sub-unity activities using
%      quartz-titanite-wollastonite experiments (Table 5, pp. 16-18).
%      Activity estimates from melt or saturation models may introduce
%      substantial uncertainty and can produce inaccurate P-T estimates
%      when aTiO2 is overestimated (discussion, p. 14).
%
%   4) Ti concentration alone cannot independently determine both pressure
%      and temperature. This implementation calculates temperature from an
%      independently supplied pressure and aTiO2. The model can instead be
%      combined with another independent equilibrium, such as Zr-in-rutile
%      or quartz-in-garnet, to obtain a unique P-T intersection
%      (application discussion, pp. 11-14).
%
%   5) The model applies to alpha- and beta-quartz. Osborne et al. (2022)
%      found no significant discontinuity in Ti solubility across the
%      alpha-beta quartz transition (pp. 7-9). The equation must not be
%      applied to coesite; a separate Ti-in-coesite model is required.
%
%   6) The high-pressure limit approaches the quartz-coesite boundary.
%      Model isopleths curve strongly near that boundary, and uncertainties
%      generally increase with increasing pressure and temperature
%      (pp. 11-13). Extrapolation beyond 30 kbar or outside the quartz
%      stability field is not recommended.
%
%   7) The model-fit uncertainty is not the complete uncertainty of a
%      natural-rock result. For example, quartz with 100 ppm Ti at
%      aTiO2 = 1 and 7.49 kbar gives approximately 700 +/- 21 degreeC from
%      fit-parameter uncertainty alone (p. 12). Natural applications must
%      additionally consider analytical error, pressure uncertainty,
%      aTiO2 uncertainty, and textural/equilibrium interpretation.
%
%   8) High-pressure quartz commonly contains low Ti concentrations, making
%      analytical accuracy especially important. Osborne et al. (2022)
%      obtained average EPMA detection limits of approximately 5-8 ppm Ti
%      using long-count, high-current trace-element analyses (p. 3).
%
%   9) To minimize secondary X-ray fluorescence from rutile inclusions, the
%      authors analyzed quartz more than 50 micrometres from rutile
%      inclusions (Electron probe microanalysis section, p. 3). Analyses too
%      close to rutile can overestimate quartz Ti and calculated temperature.
%
%  10) The collective calibration contains approximately 8-1971 ppm Ti,
%      but this is a dataset envelope rather than an independent rectangular
%      composition range. Equilibrium Ti varies strongly with both pressure
%      and temperature. Ti values outside this envelope are therefore
%      flagged as cautions rather than rejected automatically.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure lies outside 2-30 kbar,
%   2) a finite calculated temperature lies outside 550-1050 degreeC,
%   3) equivalent finite elemental Ti lies outside approximately
%      8-1971 ppm,
%   4) aTiO2 is greater than 1,
%   5) an explicitly stored calculation input is NaN,
%   6) the Ti mole fraction, activity, logarithm, denominator, or
%      temperature term is invalid, or
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
% O = 2. For the dilute binary TiO2-SiO2 substitution used by Osborne et al.
% (2022), this value is used directly as X_TiO2_quartz in the preferred Vex
% model. It is NOT treated as Ti wt% and is NOT multiplied by 1e4.
%
% For traceability and range screening, an equivalent elemental Ti
% concentration in ppm is calculated from X_TiO2_quartz using the binary
% TiO2-SiO2 mass relation corresponding to Appendix 1 on p. 18.
%
% Existing NaN values are retained and propagated; they are never replaced
% by zero. Finite Ti_cation_apfu values and entered aTiO2 values must not be
% negative. Negative finite values and Inf stop the calculation. A value of
% zero is retained as an input, but the logarithm/activity expression is
% undefined and the resulting temperature is returned as NaN with
% non-stopping diagnostic messages.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and from liquid NaN-input warnings is not
% applicable to this function.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Osborne et al. (2022), preferred Vex model, with:
%   R = 0.0083145 kJ mol^-1 K^-1
%   T = temperature in K
%   P = pressure in kbar
%   X_TiO2_quartz = mole fraction of TiO2 in quartz
%   aTiO2 = TiO2 activity referenced to pure rutile
%
%   R*T*ln(X_TiO2_quartz) =
%       -55.287
%       + P*(-2.625 + 0.0403*P)
%       + R*T*ln(aTiO2)
%
% Rearranged to calculate temperature:
%
%                  -55.287 + P*(-2.625 + 0.0403*P)
%   T(K) = ---------------------------------------------------------
%          R*[ln(X_TiO2_quartz) - ln(aTiO2)]
%
%   T(degreeC) = T(K) - 273.15
%
% The sign convention above reproduces the worked example on p. 12:
% approximately 100 ppm Ti, P = 7.49 kbar, and aTiO2 = 1 give
% approximately 700 degreeC.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Osborne2022(rawdata_struct, P_kbar)
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
    error('Osborne2022 requires (rawdata_struct, P_kbar).');
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

% Experimental calibration limits of the collective dataset.
calibrationT_min_degreeC = 550;
calibrationT_max_degreeC = 1050;
calibrationP_min_kbar = 2;
calibrationP_max_kbar = 30;

% Approximate elemental Ti envelope of the collective calibration dataset.
screeningTi_min_ppm = 8;
screeningTi_max_ppm = 1971;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3) Enter TiO2 activity
disp('=== Step 3: Entering TiO2 activity ===');

prompt = { ...
    ['Enter TiO2 activity (aTiO2):' newline ...
     'Use aTiO2 = 1 only when the selected Quartz domain' newline ...
     'co-crystallized or equilibrated with rutile.' newline ...
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
         'For the pure-rutile standard state, aTiO2 is normally expected ' ...
         'to be less than or equal to 1. The value has been retained.\n'], ...
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
            ['WARNING: Input pressure is outside the experimental ' ...
             'calibration range of Osborne et al. (2022): 2-30 kbar ' ...
             '(abstract, p. 1; Table 1, pp. 5-6; Conclusions, p. 18). ' ...
             '%d of %d pressure point(s) are outside; input range = ' ...
             '%.6g-%.6g kbar. The supplied pressures and calculated ' ...
             'results have been retained.\n'], ...
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
        row.Ti_ppm_equiv, ...
        screeningTi_min_ppm, ...
        screeningTi_max_ppm, ...
        selectedCode_qtz);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Osborne et al. (2022) ' ...
             'calculation input(s) for %s: %s.\n' ...
             '         Existing NaN values were retained and propagated; ' ...
             'the calculated temperature may be NaN.\n'], ...
            char(selectedCode_qtz), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTerms = findInvalidEquationTerms(row);

    if ~isempty(invalidTerms)
        fprintf(2, ...
            ['WARNING: Invalid Osborne Vex-model composition, activity, ' ...
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
        'Osborne2022', ...
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
        'Osborne et al. (2022) preferred expanded-volume Vex model', ...
    'doi', ...
        'https://doi.org/10.1007/s00410-022-01896-8', ...
    'experimentalCalibrationTemperatureRange_degreeC', ...
        [calibrationT_min_degreeC, calibrationT_max_degreeC], ...
    'experimentalCalibrationPressureRange_kbar', ...
        [calibrationP_min_kbar, calibrationP_max_kbar], ...
    'approximateCalibrationTiRange_ppm', ...
        [screeningTi_min_ppm, screeningTi_max_ppm], ...
    'pressureUsedInEquation', ...
        true, ...
    'aTiO2Required', ...
        true, ...
    'rutileSaturationActivity', ...
        1, ...
    'quartzPhasesIncluded', ...
        'alpha quartz and beta quartz', ...
    'coesiteApplicable', ...
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
        'Osborne2022: calculation inputs must not be negative or ' ...
        'infinite. Invalid value(s) were found in: ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_qtz, P_kbar, aTiO2)
% calcTemp
% Calculate the Osborne et al. (2022) preferred Vex-model temperature for
% one selected Quartz row and every supplied pressure. Existing NaN values
% and invalid derived terms are retained as NaN.
%
% Ti_cation_apfu is interpreted as Ti atoms per formula unit normalized to
% O = 2 and used directly as X_TiO2_quartz.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

% --- Constants and preferred Vex fit parameters ---
R = 0.0083145;  % kJ mol^-1 K^-1
A0 = -55.287;
cP = -2.625;
dP = 0.0403;

M_Ti = 47.867;
M_TiO2 = 79.866;
M_SiO2 = 60.0843;

X_TiO2_qtz_scalar = data_qtz.Ti_cation_apfu;
validateScalarVariable( ...
    X_TiO2_qtz_scalar, 'Quartz', 'Ti_cation_apfu');

% Convert X_TiO2 to an equivalent elemental Ti concentration in ppm using
% the binary TiO2-SiO2 mass relation corresponding to Appendix 1.
Ti_ppm_equiv_scalar = NaN;
compositionConversionValid = false;

if isnan(X_TiO2_qtz_scalar)
    compositionConversionValid = false;
elseif isfinite(X_TiO2_qtz_scalar) && ...
        X_TiO2_qtz_scalar >= 0 && X_TiO2_qtz_scalar <= 1

    formulaMass = ...
        X_TiO2_qtz_scalar .* M_TiO2 + ...
        (1 - X_TiO2_qtz_scalar) .* M_SiO2;

    if isfinite(formulaMass) && formulaMass > 0
        Ti_ppm_equiv_scalar = ...
            (X_TiO2_qtz_scalar .* M_Ti ./ formulaMass) .* 1e6;
        compositionConversionValid = ...
            isfinite(Ti_ppm_equiv_scalar) && Ti_ppm_equiv_scalar >= 0;
    end
end

% Repeat scalar Quartz/activity inputs for each pressure value.
X_TiO2_qtz = repmat(X_TiO2_qtz_scalar, nP, 1);
Ti_ppm_equiv = repmat(Ti_ppm_equiv_scalar, nP, 1);
aTiO2_vector = repmat(aTiO2, nP, 1);
compositionConversionValidVector = ...
    repmat(compositionConversionValid, nP, 1);

R_vector = repmat(R, nP, 1);
A0_vector = repmat(A0, nP, 1);
cP_vector = repmat(cP, nP, 1);
dP_vector = repmat(dP, nP, 1);

% Initialize derived arrays with NaN so invalid inputs remain diagnosable
% and no input NaN is silently replaced by zero.
lnX_TiO2_qtz = NaN(nP, 1);
ln_aTiO2 = NaN(nP, 1);
numerator = NaN(nP, 1);
denominator = NaN(nP, 1);
T_raw_K = NaN(nP, 1);
T_raw_degreeC = NaN(nP, 1);
T_K = NaN(nP, 1);
T_degreeC = NaN(nP, 1);

% The pressure polynomial is defined for every validated non-negative
% pressure, independent of whether composition/activity inputs are valid.
numerator(:) = ...
    A0 + P_kbar .* (cP + dP .* P_kbar);

validX = ...
    isfinite(X_TiO2_qtz) & ...
    X_TiO2_qtz > 0 & ...
    X_TiO2_qtz <= 1;

validActivity = ...
    isfinite(aTiO2_vector) & ...
    aTiO2_vector > 0;

lnX_TiO2_qtz(validX) = log(X_TiO2_qtz(validX));
ln_aTiO2(validActivity) = log(aTiO2_vector(validActivity));

validLogTerms = ...
    validX & validActivity & ...
    isfinite(lnX_TiO2_qtz) & ...
    isfinite(ln_aTiO2);

denominator(validLogTerms) = ...
    R .* (lnX_TiO2_qtz(validLogTerms) - ...
          ln_aTiO2(validLogTerms));

validDenominator = ...
    validLogTerms & ...
    isfinite(denominator) & ...
    abs(denominator) > 1e-12;

T_raw_K(validDenominator) = ...
    numerator(validDenominator) ./ denominator(validDenominator);
T_raw_degreeC(validDenominator) = ...
    T_raw_K(validDenominator) - 273.15;

% Non-positive Kelvin is physically invalid. Preserve the raw solution for
% diagnostics, but return NaN as the accepted temperature.
validTemperature = ...
    isfinite(T_raw_K) & T_raw_K > 0;

T_K(validTemperature) = T_raw_K(validTemperature);
T_degreeC(validTemperature) = T_raw_degreeC(validTemperature);

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.aTiO2 = aTiO2_vector;
row.R_kJ_molK = R_vector;
row.A0 = A0_vector;
row.cP = cP_vector;
row.dP = dP_vector;
row.X_TiO2_qtz = X_TiO2_qtz;
row.Ti_cation_apfu = X_TiO2_qtz;
row.Ti_ppm_equiv = Ti_ppm_equiv;
row.Ti_conversion_valid = compositionConversionValidVector;
row.lnX_TiO2_qtz = lnX_TiO2_qtz;
row.ln_aTiO2 = ln_aTiO2;
row.numerator = numerator;
row.denominator = denominator;
row.T_raw_K = T_raw_K;
row.T_raw_degreeC = T_raw_degreeC;
row.T_K = T_K;
row.T_degreeC = T_degreeC;
row.T_deg = T_degreeC;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid composition, activity, logarithm, denominator, and
% temperature terms. Repeated pressure rows are summarized by term name.

termBuffer = strings(9, 1);
nTerms = 0;

if any(~row.Ti_conversion_valid)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "X_TiO2-to-equivalent-Ti conversion";
end

if any(~isfinite(row.X_TiO2_qtz) | ...
        row.X_TiO2_qtz <= 0 | row.X_TiO2_qtz > 1)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = ...
        "X_TiO2_qtz (must satisfy 0 < X <= 1 for logarithm)";
end

if any(~isfinite(row.aTiO2) | row.aTiO2 <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "aTiO2 (> 0 required for logarithm)";
end

if any(~isfinite(row.lnX_TiO2_qtz))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "ln(X_TiO2_qtz)";
end

if any(~isfinite(row.ln_aTiO2))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "ln(aTiO2)";
end

if any(~isfinite(row.numerator))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Vex-model numerator";
end

if any(~isfinite(row.denominator) | abs(row.denominator) <= 1e-12)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Vex-model denominator";
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
% Warn when finite temperatures lie outside the experimental calibration
% range. Results are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);

    fprintf(2, ...
        ['WARNING: Calculated Osborne Vex-model temperature is outside ' ...
         'the experimental calibration range of %.4g-%.4g degreeC ' ...
         '(Osborne et al., 2022; abstract, p. 1; Table 1, pp. 5-6; ' ...
         'Conclusions, p. 18). %d of %d finite point(s) are outside; ' ...
         'calculated finite range = %.6g-%.6g degreeC for %s. ' ...
         'The result has been retained.\n'], ...
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
% Warn when finite positive equivalent elemental Ti lies outside the
% approximate collective calibration-data envelope. Results are retained.

finitePositiveMask = isfinite(TiValues) & TiValues > 0;
outsideMask = finitePositiveMask & ...
    (TiValues < minimumTi | TiValues > maximumTi);

if any(outsideMask)
    finiteValues = TiValues(finitePositiveMask);

    fprintf(2, ...
        ['CAUTION: Equivalent elemental Ti concentration is outside the ' ...
         'approximate %.6g-%.6g ppm envelope represented by the ' ...
         'collective rutile-saturated calibration dataset of Osborne ' ...
         'et al. (2022, Table 1, pp. 5-6). Calculated finite input ' ...
         'range = %.6g-%.6g ppm for %s. This is a composition-screening ' ...
         'caution, not a strict rejection. The result has been retained.\n'], ...
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
        ['WARNING: Non-finite Osborne Vex-model temperature values were ' ...
         'calculated for %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
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
