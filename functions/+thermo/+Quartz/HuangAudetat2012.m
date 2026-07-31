function results = HuangAudetat2012(rawdata_struct, P_kbar)
% functions/+thermo/+Quartz/HuangAudetat2012.m
% Tested with MATLAB R2024b
%
% Ti-in-Quartz thermometer / thermobarometer (TitaniQ)
% Huang, R. and Audetat, A. (2012)
% Geochimica et Cosmochimica Acta, 84, 75-89
% DOI: https://doi.org/10.1016/j.gca.2012.01.009
%
% -------------------------------------------------------------------------
% OVERVIEW
%
% This function interactively selects one Quartz analysis and calculates
% temperature using the Huang and Audetat (2012) TitaniQ calibration.
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
% Huang and Audetat (2012) grew synthetic quartz by dissolution and
% reprecipitation in rutile-bearing H2O (+/-NaCl) fluids. The experimental
% calibration conditions summarized in the abstract are:
%
%   Temperature : 600-800 degreeC
%   Pressure    : 1-10 kbar
%   TiO2 state  : rutile-bearing experiments; aTiO2 = 1 for the calibration
%   Fluid       : H2O (+/-NaCl)
%
% The experimental range, growth-rate dependence, natural-sample test, and
% principal application warnings are summarized in the abstract on p. 75.
% The detailed experimental methods, results, growth-rate discussion,
% natural-system evaluation, and summary occupy the remainder of the paper
% on pp. 76-88.
%
% Huang and Audetat (2012) also tested the calibration using igneous quartz
% from five intrusive and three volcanic magma systems with independently
% constrained conditions of:
%
%   Temperature : 675-780 degreeC
%   Pressure    : 0.8-2.7 kbar
%
% These natural-system test conditions are reported in the abstract on
% p. 75. They are a narrower validation field than the full experimental
% calibration range and are not used here as hard rejection limits.
%
% IMPORTANT APPLICATION LIMITATIONS
%
%   1) Quartz growth rate strongly affects Ti incorporation. Huang and
%      Audetat (2012) observed an increase of as much as a factor of 2.5 as
%      growth rate varied from approximately 4 to 110 micrometre/day. Their
%      calibration uses the most slowly grown quartz as the closest
%      approximation to equilibrium (abstract, p. 75; detailed discussion,
%      pp. 76-88). Rapid or strongly disequilibrium growth can therefore
%      produce excessive Ti contents and overestimated temperatures.
%
%   2) The authors conclude that the calibration should not be applied to
%      quartz grown from hydrothermal fluids because natural hydrothermal
%      growth rates may be highly variable. They consider the calibration
%      more likely to work for igneous quartz, while noting that TiO2
%      solubility/activity models still require refinement (abstract,
%      p. 75; detailed discussion and conclusions, pp. 76-88).
%
%   3) aTiO2 must be constrained independently. Use aTiO2 = 1 only when the
%      analyzed quartz grew at rutile saturation and rutile was in relevant
%      equilibrium with the quartz-forming system. The mere presence of
%      rutile elsewhere in the rock does not by itself demonstrate that the
%      selected Quartz domain formed at aTiO2 = 1.
%
%   4) In rutile-absent systems, uncertainty in the chosen aTiO2 model can
%      materially change the calculated temperature. Huang and Audetat
%      (2012) explicitly note that the natural-sample comparison depends on
%      the model used to calculate aTiO2 (abstract, p. 75; detailed
%      natural-system discussion, pp. 76-88).
%
%   5) Pressure and aTiO2 must be supplied independently to solve the
%      published relation for temperature. A single measured Ti
%      concentration cannot independently determine T, P, and aTiO2.
%
%   6) Select and analyze a texturally appropriate Quartz domain. Avoid
%      combining measured Ti from inherited cores, rapidly grown rims,
%      recrystallized domains, hydrothermal overgrowths, fractures, or
%      altered regions unless those domains are demonstrably relevant to
%      the event whose temperature is being estimated.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure lies outside 1-10 kbar,
%   2) a finite calculated temperature lies outside 600-800 degreeC,
%   3) aTiO2 is greater than 1,
%   4) an explicitly stored calculation input is NaN,
%   5) a Ti conversion, activity ratio, logarithm, denominator, or
%      temperature term is invalid, or
%   6) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
%
% rawdata_struct must contain:
%   rawdata_struct.Quartz : table
%
% The FIRST column is treated as an identifier ("data code") displayed in
% the selection dialog.
%
% The Quartz table must contain at least ONE of the following Ti variables.
% If more than one is present, the first available variable in this order is
% used for all selected rows:
%
%   1) Ti_ppm
%      Elemental Ti concentration in microgram/gram (= ppm). This is the
%      preferred and most direct input for the published equation.
%
%   2) TiO2_wtpercent
%      TiO2 concentration in wt%. It is converted to elemental Ti ppm by:
%
%        Ti_ppm = TiO2_wtpercent * (M_Ti / M_TiO2) * 1e4
%
%   3) Ti_cation_apfu
%      Ti cations per formula unit normalized to two oxygens. It is
%      converted to elemental Ti ppm using an ideal quartz approximation:
%
%        quartz formula = (Si_(1-Ti) Ti_Ti) O2
%
%      This conversion replaces the scientifically incorrect legacy
%      conversion Ti_ppm = Ti_cation_apfu * 1e4. Direct Ti_ppm or measured
%      TiO2_wtpercent is preferred whenever available.
%
% Existing NaN values are retained and propagated; they are never replaced
% by zero. All finite stored Ti inputs and the entered aTiO2 must be greater
% than or equal to zero. Negative finite values and Inf stop the
% calculation. A value of zero is retained, but produces an invalid
% logarithm/activity term and therefore a NaN temperature with a
% non-stopping diagnostic message.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and from liquid NaN-input warnings is not
% applicable to this function.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Huang and Audetat (2012), with T in K, P in kbar, and elemental Ti in ppm:
%
%   log10(Ti_ppm) =
%       -2794.3 / T
%       -660.53 * P^0.35 / T
%       +5.6459
%
% Including TiO2 activity:
%
%   log10(Ti_ppm / aTiO2) =
%       -(2794.3 + 660.53 * P^0.35) / T
%       +5.6459
%
% Rearranged to calculate temperature:
%
%                 2794.3 + 660.53 * P^0.35
%   T(K) = ------------------------------------------------
%          5.6459 - log10(Ti_ppm / aTiO2)
%
%   T(degreeC) = T(K) - 273.15
%
% IMPORTANT:
% The signs above correct the sign errors in the original local script.
% Increasing pressure at constant Ti_ppm and aTiO2 increases the calculated
% temperature in this formulation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = HuangAudetat2012(rawdata_struct, P_kbar)
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
    error('HuangAudetat2012 requires (rawdata_struct, P_kbar).');
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
tiInputSource = determineTiInputSource(dataset_qtz);

disp(['Ti input source: Quartz.' char(tiInputSource)]);
disp('=== Preparing Quartz dataset has been finished ===');

%% 2) Initialize output container and calibration screening settings
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Experimental calibration limits summarized by Huang and Audetat (2012,
% abstract, p. 75).
calibrationT_min_degreeC = 600;
calibrationT_max_degreeC = 800;
calibrationP_min_kbar = 1;
calibrationP_max_kbar = 10;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3) Enter TiO2 activity
disp('=== Step 3: Entering TiO2 activity ===');

prompt = { ...
    ['Enter TiO2 activity (aTiO2):' newline ...
     'Use aTiO2 = 1 only for quartz grown at rutile saturation.' newline ...
     'NaN is allowed and will be retained and reported.']};

dlgtitle = 'TiO2 activity input';
dims = [3 70];
definput = {'1'};

answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    error('TiO2 activity input canceled.');
end

aTiO2 = str2double(answer{1});

% NaN is deliberately allowed and propagated. Negative finite values and
% Inf are prohibited. Zero is retained and diagnosed after calculation.
if isinf(aTiO2) || (isfinite(aTiO2) && aTiO2 < 0)
    error('aTiO2 must not be negative or infinite. NaN is allowed.');
end

disp(['TiO2 activity (aTiO2) = ' num2str(aTiO2)]);
disp('NOTE: use aTiO2 = 1 only when rutile saturation is independently justified.');

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
    nanInputNames = findNaNInputs( ...
        selectedData_qtz, tiInputSource, aTiO2);

    % Stop only for negative finite values or Inf. NaN and zero are retained.
    validateNonNegativeInputs( ...
        selectedData_qtz, tiInputSource, aTiO2);

    row = calcTemp( ...
        selectedData_qtz, P_kbar, aTiO2, tiInputSource);

    nRows = height(row);
    row.dataCode_qtz = repmat(selectedCode_qtz, nRows, 1);
    row = movevars(row, 'dataCode_qtz', 'Before', 1);

    % Store one completed table block. The full output table is not enlarged
    % on every interactive iteration.
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
             'calibration range of Huang and Audetat (2012): 1-10 kbar ' ...
             '(abstract, p. 75). %d of %d pressure point(s) are outside; ' ...
             'input range = %.6g-%.6g kbar. The supplied pressure values ' ...
             'and calculated results have been retained.\n'], ...
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

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Huang-Audetat calculation ' ...
             'input(s) for %s: %s.\n' ...
             '         Existing NaN values were retained and propagated; ' ...
             'the calculated temperature may be NaN.\n'], ...
            char(selectedCode_qtz), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTerms = findInvalidEquationTerms(row);

    if ~isempty(invalidTerms)
        fprintf(2, ...
            ['WARNING: Invalid Huang-Audetat Ti conversion, activity, ' ...
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
        'HuangAudetat2012', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate the buffered result blocks only once after all selections have
% been completed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'primaryTemperatureEquation', ...
        'Huang and Audetat (2012) TitaniQ calibration', ...
    'doi', ...
        'https://doi.org/10.1016/j.gca.2012.01.009', ...
    'tiInputSource', ...
        char(tiInputSource), ...
    'experimentalCalibrationTemperatureRange_degreeC', ...
        [calibrationT_min_degreeC, calibrationT_max_degreeC], ...
    'experimentalCalibrationPressureRange_kbar', ...
        [calibrationP_min_kbar, calibrationP_max_kbar], ...
    'naturalIgneousTestTemperatureRange_degreeC', ...
        [675, 780], ...
    'naturalIgneousTestPressureRange_kbar', ...
        [0.8, 2.7], ...
    'hydrothermalQuartzRecommended', ...
        false, ...
    'pressureUsedInEquation', ...
        true, ...
    'aTiO2Required', ...
        true);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function tiInputSource = determineTiInputSource(dataset_qtz)
% determineTiInputSource
% Select one Ti input column for the complete function call. Direct elemental
% Ti ppm is preferred, followed by TiO2 wt%, then normalized Ti apfu.

candidateVariables = { ...
    'Ti_ppm', ...
    'TiO2_wtpercent', ...
    'Ti_cation_apfu'};

for i = 1:numel(candidateVariables)
    variableName = candidateVariables{i};

    if ismember(variableName, dataset_qtz.Properties.VariableNames)
        tiInputSource = string(variableName);
        return;
    end
end

error([ ...
    'Quartz table must contain at least one Ti input variable: ' ...
    'Ti_ppm, TiO2_wtpercent, or Ti_cation_apfu.']);

end

function nanInputNames = findNaNInputs( ...
        data_qtz, tiInputSource, aTiO2)
% findNaNInputs
% Return names of explicitly stored calculation inputs that are NaN. NaN is
% reported but does not stop the calculation.

nameBuffer = strings(2, 1);
nNames = 0;

value = data_qtz.(char(tiInputSource));
validateScalarVariable(value, 'Quartz', char(tiInputSource));

if isnan(value)
    nNames = nNames + 1;
    nameBuffer(nNames) = "Quartz." + tiInputSource;
end

if isnan(aTiO2)
    nNames = nNames + 1;
    nameBuffer(nNames) = "aTiO2";
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs( ...
        data_qtz, tiInputSource, aTiO2)
% validateNonNegativeInputs
% Stop when a stored Ti input or entered activity is negative or infinite.
% Zero and NaN are deliberately allowed and handled by non-stopping
% diagnostics after calculation.

value = data_qtz.(char(tiInputSource));
validateScalarVariable(value, 'Quartz', char(tiInputSource));

invalidNames = strings(2, 1);
nInvalid = 0;

if isinf(value) || (isfinite(value) && value < 0)
    nInvalid = nInvalid + 1;
    invalidNames(nInvalid) = "Quartz." + tiInputSource;
end

if isinf(aTiO2) || (isfinite(aTiO2) && aTiO2 < 0)
    nInvalid = nInvalid + 1;
    invalidNames(nInvalid) = "aTiO2";
end

if nInvalid > 0
    invalidNames = invalidNames(1:nInvalid);

    error([ ...
        'HuangAudetat2012: calculation inputs must not be negative ' ...
        'or infinite. Invalid value(s) were found in: ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_qtz, P_kbar, aTiO2, tiInputSource)
% calcTemp
% Calculate the Huang and Audetat (2012) TitaniQ temperature for one selected
% Quartz row and every supplied pressure. Existing NaN values are retained.
% Invalid derived terms become NaN rather than stopping the calculation.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

% Atomic/molecular masses used only when conversion from TiO2 wt% or Ti apfu
% is required.
M_Ti = 47.867;
M_Si = 28.0855;
M_O = 15.999;
M_TiO2 = M_Ti + 2 .* M_O;

[Ti_input_value, Ti_ppm_scalar, conversionValid, ...
    Ti_cation_apfu_scalar, TiO2_wtpercent_scalar] = ...
    convertTiToPpm( ...
        data_qtz, tiInputSource, M_Ti, M_Si, M_O, M_TiO2);

% Repeat scalar Quartz/activity inputs for each pressure value.
Ti_input_value_vector = repmat(Ti_input_value, nP, 1);
Ti_ppm = repmat(Ti_ppm_scalar, nP, 1);
aTiO2_vector = repmat(aTiO2, nP, 1);
conversionValidVector = repmat(conversionValid, nP, 1);

% Initialize all derived arrays with NaN so invalid inputs are retained as
% diagnosable NaN results without conditional changes in array size.
Ti_activity_ratio = NaN(nP, 1);
log_term = NaN(nP, 1);
pressure_term = NaN(nP, 1);
numerator = NaN(nP, 1);
denominator = NaN(nP, 1);
T_raw_K = NaN(nP, 1);
T_raw_degreeC = NaN(nP, 1);
T_K = NaN(nP, 1);
T_degreeC = NaN(nP, 1);

% Pressure is finite and non-negative because it was validated by the main
% function, so this term is always defined.
pressure_term(:) = P_kbar .^ 0.35;
numerator(:) = 2794.3 + 660.53 .* pressure_term;

% Calculate only where elemental Ti, activity, and the Ti conversion are
% finite and strictly positive. Zero is retained as an input but leads to an
% invalid logarithm and therefore remains NaN in the derived fields.
validLogInput = ...
    conversionValidVector & ...
    isfinite(Ti_ppm) & Ti_ppm > 0 & ...
    isfinite(aTiO2_vector) & aTiO2_vector > 0;

Ti_activity_ratio(validLogInput) = ...
    Ti_ppm(validLogInput) ./ aTiO2_vector(validLogInput);

validRatio = ...
    validLogInput & ...
    isfinite(Ti_activity_ratio) & Ti_activity_ratio > 0;

log_term(validRatio) = log10(Ti_activity_ratio(validRatio));
denominator(validRatio) = 5.6459 - log_term(validRatio);

validDenominator = ...
    validRatio & ...
    isfinite(denominator) & abs(denominator) > 1e-12;

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
row.Ti_input_source = repmat(tiInputSource, nP, 1);
row.Ti_input_value = Ti_input_value_vector;
row.Ti_cation_apfu = repmat(Ti_cation_apfu_scalar, nP, 1);
row.TiO2_wtpercent = repmat(TiO2_wtpercent_scalar, nP, 1);
row.Ti_ppm = Ti_ppm;
row.Ti_conversion_valid = conversionValidVector;
row.Ti_activity_ratio = Ti_activity_ratio;
row.log10_Ti_over_aTiO2 = log_term;
row.pressure_term_P_kbar_power_0p35 = pressure_term;
row.numerator = numerator;
row.denominator = denominator;
row.T_raw_K = T_raw_K;
row.T_raw_degreeC = T_raw_degreeC;
row.T_K = T_K;
row.T_degreeC = T_degreeC;
row.T_deg = T_degreeC;

end

function [Ti_input_value, Ti_ppm, conversionValid, ...
        Ti_cation_apfu, TiO2_wtpercent] = convertTiToPpm( ...
        data_qtz, tiInputSource, M_Ti, M_Si, M_O, M_TiO2)
% convertTiToPpm
% Convert the selected Ti input source to elemental Ti ppm. NaN is retained.
% The Ti-apfu conversion assumes an ideal two-oxygen quartz formula with
% Si + Ti = 1 cation per formula unit.

sourceName = char(tiInputSource);
Ti_input_value = data_qtz.(sourceName);
validateScalarVariable(Ti_input_value, 'Quartz', sourceName);

Ti_ppm = NaN;
Ti_cation_apfu = NaN;
TiO2_wtpercent = NaN;
conversionValid = false;

if strcmp(sourceName, 'Ti_ppm')
    Ti_ppm = Ti_input_value;

    if isnan(Ti_input_value)
        conversionValid = false;
    elseif isfinite(Ti_input_value) && Ti_input_value >= 0
        conversionValid = true;
    end

elseif strcmp(sourceName, 'TiO2_wtpercent')
    TiO2_wtpercent = Ti_input_value;

    if isnan(Ti_input_value)
        conversionValid = false;
    elseif isfinite(Ti_input_value) && Ti_input_value >= 0
        Ti_ppm = Ti_input_value .* (M_Ti ./ M_TiO2) .* 1e4;
        conversionValid = isfinite(Ti_ppm) && Ti_ppm >= 0;
    end

elseif strcmp(sourceName, 'Ti_cation_apfu')
    Ti_cation_apfu = Ti_input_value;

    if isnan(Ti_input_value)
        conversionValid = false;
    elseif isfinite(Ti_input_value) && ...
            Ti_input_value >= 0 && Ti_input_value <= 1

        idealSi_apfu = 1 - Ti_input_value;
        formulaMass = ...
            idealSi_apfu .* M_Si + ...
            Ti_input_value .* M_Ti + ...
            2 .* M_O;

        if isfinite(formulaMass) && formulaMass > 0
            Ti_ppm = ...
                (Ti_input_value .* M_Ti ./ formulaMass) .* 1e6;
            conversionValid = isfinite(Ti_ppm) && Ti_ppm >= 0;
        end
    end

else
    error('Unsupported Ti input source: %s', sourceName);
end

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid Ti conversion, activity, logarithm, denominator, and
% temperature terms. Repeated pressure rows are summarized by term name.

termBuffer = strings(10, 1);
nTerms = 0;

if any(~row.Ti_conversion_valid)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Ti-to-ppm conversion";
end

if any(~isfinite(row.Ti_ppm) | row.Ti_ppm <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "elemental Ti_ppm (> 0 required for logarithm)";
end

if any(~isfinite(row.aTiO2) | row.aTiO2 <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "aTiO2 (> 0 required for activity ratio)";
end

if any(~isfinite(row.Ti_activity_ratio) | row.Ti_activity_ratio <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Ti_ppm/aTiO2";
end

if any(~isfinite(row.log10_Ti_over_aTiO2))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "log10(Ti_ppm/aTiO2)";
end

if any(~isfinite(row.numerator))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "temperature numerator";
end

if any(~isfinite(row.denominator) | abs(row.denominator) <= 1e-12)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "temperature denominator";
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
        ['WARNING: Calculated Huang-Audetat temperature is outside the ' ...
         'experimental calibration range of 600-800 degreeC ' ...
         '(abstract, p. 75). %d of %d finite point(s) are outside; ' ...
         'calculated finite range = %.6g-%.6g degreeC for %s. ' ...
         'The result has been retained.\n'], ...
        sum(outsideMask), ...
        sum(finiteMask), ...
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
        ['WARNING: Non-finite Huang-Audetat temperature values were ' ...
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
