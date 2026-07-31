function results = KawasakiOsanai2008(rawdata_struct, P_kbar)
% functions/+thermo/+Quartz/KawasakiOsanai2008.m
% Tested with MATLAB R2024b
%
% Empirical Ti-in-Quartz geothermometer
% Kawasaki, T. and Osanai, Y. (2008)
% Geological Society, London, Special Publications, 308, 419-430
% DOI: https://doi.org/10.1144/SP308.21
%
% -------------------------------------------------------------------------
% OVERVIEW
%
% This function interactively selects one Quartz analysis and calculates
% temperature using Equation (14) of Kawasaki and Osanai (2008).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Quartz analysis, one output row is returned for every
% pressure supplied in P_kbar. Equation (14) contains no explicit pressure
% term, so all pressure rows for one selected analysis contain the same
% calculated temperature. This interface is compatible with both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Kawasaki and Osanai (2008) derived Equation (14) mainly from natural
% quartz-rutile pairs in geologically and petrologically constrained
% metamorphic rocks. It is therefore a natural-sample empirical calibration,
% not a broad experimental calibration over a formally rectangular P-T
% field.
%
% APPROXIMATE NATURAL CALIBRATION-DATA ENVELOPE
%
%   Temperature : approximately 607-1120 degreeC
%                 This interval covers the independent peak and retrograde
%                 temperature estimates compiled in Table 2 on p. 427.
%
%   Pressure    : approximately 10-25 kbar for the principal calibration
%                 datasets discussed by the authors. On p. 429, Kawasaki
%                 and Osanai (2008) summarize approximately 16-25 kbar for
%                 the Sanbagawa quartz eclogite and approximately 10-15
%                 kbar for the East Antarctic ultrahigh-temperature
%                 granulites.
%
%   Quartz Ti   : approximately 0.00029-0.00255 Ti atoms per formula unit,
%                 normalized to O = 2, based on Table 2 on p. 427.
%
% IMPORTANT:
% These numerical intervals are practical envelopes represented by the
% natural calibration and comparison dataset. Kawasaki and Osanai (2008)
% did not define a single formal rectangular temperature-pressure-
% composition calibration field. The intervals are used here only for
% non-stopping screening warnings and must not be interpreted as strict
% validity boundaries.
%
% IMPORTANT APPLICATION LIMITATIONS
%
%   1) Equation (14) is intended for quartz coexisting and equilibrated with
%      rutile. Rutile saturation is a fundamental assumption of the
%      thermodynamic model and empirical calibration (pp. 419-421; Eq. 14
%      and discussion on p. 429). The presence of rutile elsewhere in the
%      rock does not by itself prove that the selected Quartz domain was in
%      equilibrium with rutile at the stage of interest.
%
%   2) The thermometer was developed primarily for high-temperature to
%      ultrahigh-temperature metamorphic rocks. Its use for low-pressure
%      igneous quartz, hydrothermal quartz, or other geological settings was
%      not established by the calibration dataset.
%
%   3) Pressure is not included explicitly in Equation (14). The authors
%      neglected the pressure effect because the metamorphic pressures of
%      the natural calibration datasets were restricted to a comparatively
%      narrow interval, with a difference of about 9 kbar at most (p. 429).
%      This must not be interpreted as proof that Ti solubility in quartz is
%      pressure independent over all geological pressures.
%
%   4) Quartz microstructure must be considered. High-Ti quartz in direct
%      contact with rutile may preserve peak conditions, whereas low-Ti
%      quartz containing exsolved rutile needles or recrystallized during
%      cooling may record a retrograde stage (pp. 426-428; Table 2, p. 427).
%      Do not combine unrelated cores, rims, recrystallized domains, or
%      exsolution-bearing regions without a textural interpretation.
%
%   5) Fine rutile lamellae or grains included in an electron-microprobe
%      analysis can increase the measured bulk Ti concentration. Point
%      analyses of dissolved Ti in quartz and area analyses intended to
%      reconstruct pre-exsolution bulk Ti represent different quantities
%      and must not be interpreted interchangeably (Table 1, pp. 422-423;
%      discussion on pp. 424-426).
%
%   6) To minimize secondary X-ray fluorescence from adjacent rutile, the
%      authors analyzed quartz at positions approximately 10 micrometres or
%      more from the rutile-quartz boundary (Chemical analyses, p. 422).
%
%   7) The reported electron-microprobe detection limit for Ti was about
%      200 ppm (p. 422). Low-Ti quartz near the analytical detection limit
%      can carry large relative uncertainty, which propagates strongly
%      through the logarithmic temperature equation.
%
%   8) The 1-atm, 1300 degreeC experiment produced cristobalite and was not
%      included in the least-squares calibration of Equation (14) (p. 429).
%      It must not be used to extend the quartz calibration to 1300 degreeC
%      or to atmospheric pressure.
%
%   9) Kawasaki and Osanai (2008) found that the Wark and Watson (2006)
%      TitaniQ calibration yielded temperatures approximately 200 degreeC
%      higher for their samples. The two calibrations are not
%      interchangeable because they use different variables and empirical
%      datasets (abstract, p. 419; Table 2, p. 427; discussion on p. 429).
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure lies outside approximately 10-25 kbar,
%   2) a finite calculated temperature lies outside approximately
%      607-1120 degreeC,
%   3) finite Ti_cation_apfu lies outside approximately
%      0.00029-0.00255,
%   4) an explicitly stored calculation input is NaN,
%   5) the logarithm, denominator, or temperature term is invalid, or
%   6) a calculated temperature is NaN or Inf.
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
% O = 2. It is used directly as X_TiO2_Qtz in Equation (14). Do not supply
% elemental Ti ppm, TiO2 ppm, or TiO2 wt% without first converting and
% normalizing the analysis to Ti cations per formula unit on an O = 2 basis.
%
% Existing NaN values are retained and propagated; they are never replaced
% by zero. Finite Ti_cation_apfu values must be greater than or equal to
% zero. Negative finite values and Inf stop the calculation. A value of zero
% is retained as an input, but the logarithm is undefined and the resulting
% temperature is returned as NaN with a non-stopping diagnostic message.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and from liquid NaN-input warnings is not
% applicable to this function.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Kawasaki and Osanai (2008), Equation (14), p. 429:
%
%                            -5895
%   T(K) = --------------------------------------------
%           ln(X_TiO2_Qtz) + 1.729
%
%   T(degreeC) = T(K) - 273
%
% where:
%
%   X_TiO2_Qtz = Ti atoms per formula unit normalized to O = 2
%
% IMPORTANT SIGN NOTE:
% The negative sign before 5895 is required. It follows from the Arrhenius
% relation and reproduces the temperatures in Table 2 on p. 427. Omitting
% this negative sign produces physically incorrect negative temperatures.
%
% Pressure is accepted and stored for interface compatibility, but is not
% used explicitly in Equation (14).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = KawasakiOsanai2008(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing a Quartz table
%   P_kbar         : finite non-negative numeric scalar or vector; stored in
%                    the output but not used explicitly by Equation (14)
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Quartz analysis. T_K, T_degreeC, and T_deg are supplied for
%             downstream launcher and plotting compatibility.
%

%% Input validation
if nargin < 2
    error('KawasakiOsanai2008 requires (rawdata_struct, P_kbar).');
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

% Approximate natural calibration-data envelopes, not a formal rectangular
% experimental calibration field.
screeningT_min_degreeC = 607;
screeningT_max_degreeC = 1120;
screeningP_min_kbar = 10;
screeningP_max_kbar = 25;
screeningTi_min_apfu = 0.00029;
screeningTi_max_apfu = 0.00255;

pressureOutsideScreening = ...
    P_kbar < screeningP_min_kbar | ...
    P_kbar > screeningP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive Quartz selection and calculation
dataCodes_qtz = dataset_qtz{:, 1};
displayCodes_qtz = cellstr(string(dataCodes_qtz));

disp('=== Step 3: Selecting a data code from the list (Quartz) ===');

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

    disp('=== Step 4: Calculating the temperature ===');

    % Report explicitly stored NaN inputs without replacing them by zero.
    nanInputNames = findNaNInputs(selectedData_qtz);

    % Stop only for negative finite values or Inf. NaN and zero are retained.
    validateNonNegativeInputs(selectedData_qtz);

    row = calcTemp(selectedData_qtz, P_kbar);

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

    % Pressure is not used in Equation (14), but the natural calibration
    % datasets occupy a restricted pressure range. Print this warning only
    % once for the complete function call.
    if any(pressureOutsideScreening) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate ' ...
             '10-25 kbar pressure envelope of the principal natural ' ...
             'calibration datasets discussed by Kawasaki and Osanai ' ...
             '(2008, p. 429). %d of %d pressure point(s) are outside; ' ...
             'input range = %.6g-%.6g kbar. This is a practical ' ...
             'source-data screening interval, not a formally stated ' ...
             'rectangular calibration range. Equation (14) contains no ' ...
             'explicit pressure term, and all supplied pressure values ' ...
             'and calculated results have been retained.\n'], ...
            sum(pressureOutsideScreening), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    printTemperatureScreeningWarning( ...
        row.T_degreeC, ...
        screeningT_min_degreeC, ...
        screeningT_max_degreeC, ...
        selectedCode_qtz);

    printTiScreeningWarning( ...
        row.X_TiO2_Qtz, ...
        screeningTi_min_apfu, ...
        screeningTi_max_apfu, ...
        selectedCode_qtz);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Kawasaki-Osanai Equation ' ...
             '(14) input(s) for %s: %s.\n' ...
             '         Existing NaN values were retained and propagated; ' ...
             'the calculated temperature may be NaN.\n'], ...
            char(selectedCode_qtz), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTerms = findInvalidEquationTerms(row);

    if ~isempty(invalidTerms)
        fprintf(2, ...
            ['WARNING: Invalid Kawasaki-Osanai logarithm, denominator, ' ...
             'or temperature term(s) were found for %s: %s.\n' ...
             '         Affected values were retained as NaN, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(selectedCode_qtz), ...
            char(strjoin(invalidTerms, ', ')));
    end

    printNonFiniteTemperatureWarning(row.T_degreeC, selectedCode_qtz);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'KawasakiOsanai2008', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after the selection loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'primaryTemperatureEquation', ...
        'Kawasaki and Osanai (2008) Equation (14)', ...
    'doi', ...
        'https://doi.org/10.1144/SP308.21', ...
    'calibrationType', ...
        'Natural-sample empirical quartz-rutile thermometer', ...
    'pressureUsedInEquation', ...
        false, ...
    'approximateNaturalDatasetTemperatureRange_degreeC', ...
        [screeningT_min_degreeC, screeningT_max_degreeC], ...
    'approximatePrincipalDatasetPressureRange_kbar', ...
        [screeningP_min_kbar, screeningP_max_kbar], ...
    'approximateNaturalDatasetTiRange_apfu_O2', ...
        [screeningTi_min_apfu, screeningTi_max_apfu], ...
    'rutileSaturationRequired', ...
        true, ...
    'rangeStatus', ...
        ['Practical natural source-dataset screening envelopes; not a ' ...
         'formal rectangular experimental calibration field']);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_qtz)
% findNaNInputs
% Return the names of explicitly stored calculation inputs that contain NaN.
% NaN values are reported but do not stop the calculation.

value = data_qtz.Ti_cation_apfu;
validateScalarVariable(value, 'Quartz', 'Ti_cation_apfu');

nameBuffer = strings(1, 1);
nNames = 0;

if isnan(value)
    nNames = 1;
    nameBuffer(nNames) = "Quartz.Ti_cation_apfu";
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_qtz)
% validateNonNegativeInputs
% Stop when the stored Ti input is negative or infinite. Zero and NaN are
% deliberately allowed and handled by non-stopping diagnostics after the
% calculation.

value = data_qtz.Ti_cation_apfu;
validateScalarVariable(value, 'Quartz', 'Ti_cation_apfu');

if isinf(value) || (isfinite(value) && value < 0)
    error([ ...
        'KawasakiOsanai2008: Quartz.Ti_cation_apfu must not be ' ...
        'negative or infinite. NaN and zero are allowed and retained.']);
end

end

function row = calcTemp(data_qtz, P_kbar)
% calcTemp
% Calculate Kawasaki and Osanai (2008) Equation (14) for one selected Quartz
% row and repeat the pressure-independent result for every supplied pressure.
% Existing NaN values and invalid derived terms are retained as NaN.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

X_TiO2_Qtz_scalar = data_qtz.Ti_cation_apfu;
validateScalarVariable( ...
    X_TiO2_Qtz_scalar, 'Quartz', 'Ti_cation_apfu');

% Initialize scalar derived terms with NaN. Zero and NaN Ti inputs therefore
% remain diagnosable without being replaced by zero or causing an early exit.
denominator_scalar = NaN;
T_raw_K_scalar = NaN;
T_raw_degreeC_scalar = NaN;
T_K_scalar = NaN;
T_degreeC_scalar = NaN;

if isfinite(X_TiO2_Qtz_scalar) && X_TiO2_Qtz_scalar > 0
    denominator_scalar = log(X_TiO2_Qtz_scalar) + 1.729;

    if isfinite(denominator_scalar) && abs(denominator_scalar) > 1e-12
        % Kawasaki and Osanai (2008), Equation (14). The negative sign is
        % required to reproduce Table 2 and the Arrhenius relation.
        T_raw_K_scalar = -5895 ./ denominator_scalar;
        T_raw_degreeC_scalar = T_raw_K_scalar - 273;

        % Non-positive Kelvin is physically invalid. Preserve the raw
        % solution for diagnosis, but return NaN as the accepted temperature.
        if isfinite(T_raw_K_scalar) && T_raw_K_scalar > 0
            T_K_scalar = T_raw_K_scalar;
            T_degreeC_scalar = T_raw_degreeC_scalar;
        end
    end
end

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PressureUsedInEquation = false(nP, 1);
row.PrimaryEquation = repmat("KawasakiOsanai2008_Eq14", nP, 1);
row.X_TiO2_Qtz = repmat(X_TiO2_Qtz_scalar, nP, 1);
row.denominator = repmat(denominator_scalar, nP, 1);
row.T_raw_K = repmat(T_raw_K_scalar, nP, 1);
row.T_raw_degreeC = repmat(T_raw_degreeC_scalar, nP, 1);
row.T_K_equiv = repmat(T_K_scalar, nP, 1);
row.T_K = repmat(T_K_scalar, nP, 1);
row.T_degreeC = repmat(T_degreeC_scalar, nP, 1);
row.T_deg = row.T_degreeC;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid logarithm, denominator, and temperature terms. Repeated
% pressure rows are summarized by term name.

termBuffer = strings(5, 1);
nTerms = 0;

if any(~isfinite(row.X_TiO2_Qtz) | row.X_TiO2_Qtz <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = ...
        "X_TiO2_Qtz (> 0 required for natural logarithm)";
end

if any(~isfinite(row.denominator) | abs(row.denominator) <= 1e-12)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Equation (14) denominator";
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

if isscalar(temperatureValues)
    disp([label ': ' num2str(temperatureValues) ' degreeC']);
else
    disp([label ': ' ...
        num2str(temperatureValues(1)) ' to ' ...
        num2str(temperatureValues(end)) ' degreeC']);
end

end

function printTemperatureScreeningWarning( ...
        temperatureValues, minimumTemperature, maximumTemperature, ...
        selectedCode_qtz)
% printTemperatureScreeningWarning
% Warn when finite temperatures lie outside the approximate natural-data
% envelope represented in Table 2. Results are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);

    fprintf(2, ...
        ['WARNING: Calculated Kawasaki-Osanai temperature is outside ' ...
         'the approximate %.4g-%.4g degreeC natural calibration-data ' ...
         'envelope represented in Table 2 of Kawasaki and Osanai ' ...
         '(2008, p. 427). %d of %d finite point(s) are outside; ' ...
         'calculated finite range = %.6g-%.6g degreeC for %s. This is ' ...
         'a practical source-data screening envelope, not a formally ' ...
         'stated rectangular calibration range. The result has been ' ...
         'retained.\n'], ...
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
% Warn when finite positive Ti apfu lies outside the approximate range of
% the natural calibration/comparison dataset in Table 2. Results are retained.

finitePositiveMask = isfinite(TiValues) & TiValues > 0;
outsideMask = finitePositiveMask & ...
    (TiValues < minimumTi | TiValues > maximumTi);

if any(outsideMask)
    finiteValues = TiValues(finitePositiveMask);

    fprintf(2, ...
        ['CAUTION: Quartz Ti_cation_apfu is outside the approximate ' ...
         '%.6g-%.6g range represented in Kawasaki and Osanai (2008, ' ...
         'Table 2, p. 427). Calculated finite input range = ' ...
         '%.6g-%.6g for %s. This is a composition-screening caution, ' ...
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
        ['WARNING: Non-finite Kawasaki-Osanai Equation (14) temperature ' ...
         'values were calculated for %s (%d of %d points; NaN: %d, ' ...
         'Inf: %d).\n' ...
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
