function results = Powell1985(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Clinopyroxene/Powell1985.m
% Tested with MATLAB R2024b
%
% Garnet-clinopyroxene Fe-Mg exchange geothermometer
% Powell, R. (1985)
% Journal of Metamorphic Geology, 3, 231-243
% DOI: https://doi.org/10.1111/j.1525-1314.1985.tb00319.x
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis (selected by the user from tables) and calculates
% temperature using the Powell (1985) garnet-clinopyroxene Fe-Mg exchange
% geothermometer.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Grt-Cpx pair. One output row is produced
% for every supplied pressure value, allowing the same function to operate
% with both startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Powell (1985) re-evaluated the experimental garnet-clinopyroxene data of
% Raheim and Green (1974) and Ellis and Green (1979) using regression
% diagnostics, robust regression, and jackknife uncertainty estimates. The
% preferred temperature-composition regression is based on 49 experiments at
% 30 kbar after the data screening and influential-outlier exclusions
% described on pp. 235-237 and 241 (Tables 1-3).
%
% The numerical envelope of the 30 kbar calibration data in Table 3 is:
%
%   Temperature       : 700-1500 degreeC
%   Pressure          : 30 kbar (direct preferred calibration)
%   Garnet XCa        : 0.056-0.459
%
% These limits are the envelope of the experiments used in the preferred
% regression (Table 3, pp. 238-239), not hard validity limits explicitly
% guaranteed by the author. Calculations outside these limits are
% extrapolations and are retained with non-stopping fprintf warnings.
%
% Most experiments are at 30 kbar, and the pressure spread is explicitly
% described as insufficient to determine the pressure coefficient by
% regression (pp. 235-236). The coefficient used in the preferred equation
% adopts V' = -10; its uncertainty is difficult to estimate, and the pressure
% contribution to uncertainty increases with distance from 30 kbar
% (pp. 241-242, equations 17-18). Powell (1985) does not provide a numerical
% lower and upper pressure application range. This implementation therefore
% warns whenever an input pressure differs from the directly calibrated
% pressure of 30 kbar, without inventing an unsupported pressure interval.
%
% Calibration data were screened to remove experiments affected by iron loss,
% disequilibrium, or large quoted KD uncertainty/range (p. 235). Influential
% experiments RG 4732 and EG 43A were also excluded (p. 237). Natural mineral
% pairs should therefore represent texturally corresponding, equilibrated
% garnet and clinopyroxene domains.
%
% Calibration uncertainty is composition- and temperature-dependent rather
% than one constant +/- value (Table 3 and equations 14-16, pp. 238-241).
% Additional uncertainty arises from ln(KD), microprobe analyses, and the
% adopted pressure term. Continued Fe-Mg re-equilibration during cooling may
% cause the calculated temperature not to represent the intended metamorphic
% event (p. 242).
%
% This implementation issues non-stopping warnings when:
%   1) input pressure differs from 30 kbar,
%   2) a finite calculated temperature is outside 700-1500 degreeC,
%   3) finite garnet XCa is outside 0.056-0.459, or
%   4) a required calculation input or calculated temperature is non-finite.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet          : table
% or
%   rawdata_struct.Grt             : table
%
% and
%
%   rawdata_struct.Clinopyroxene   : table
% or
%   rawdata_struct.Cpx             : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns contain normalized
% cation data as apfu.
%
% Required Garnet variables:
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%
% Required Clinopyroxene variables:
%   Fe_cation_apfu
%   Mg_cation_apfu
%
% Optional variables in both tables:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Si_cation_apfu
%   Al_cation_apfu
%
% Optional Clinopyroxene variables:
%   Ca_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%   Ti_cation_apfu
%
% IMPORTANT Fe convention:
% To preserve the numerical convention of the original implementation, this
% function uses
%
%   Fe_used = Fe_cation_apfu + Fe3_cation_apfu
%
% in the exchange expression. This is appropriate only when
% Fe_cation_apfu and Fe3_cation_apfu are separate quantities. If
% Fe_cation_apfu already contains total Fe, Fe3_cation_apfu must not be added
% again. Powell (1985) calibrates Fe-Mg exchange but does not provide a
% general ferric-iron recalculation procedure for natural analyses.
%
% Garnet Mn is included in the denominator of XCa_g. An absent Mn column
% retains the legacy assumption Mn = 0; an explicitly supplied NaN remains
% NaN and propagates into XCa_g and temperature. An absent Fe3 column likewise
% defaults to zero, whereas explicit Fe3 = NaN is preserved. Other missing
% optional output variables are stored as NaN.
%
% All supplied finite cation values must be greater than or equal to zero.
% Negative values and Inf values stop the calculation. NaN values do not stop
% the calculation and are never converted to zero. NaN in a thermometer input
% propagates to the calculated temperature and generates a non-stopping
% fprintf warning.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
%   T(K) = (2790 + 10*P_kbar + 3140*XCa_g) / (1.735 + ln(KD))
%
% where
%
%   KD    = (Fe_used/Mg)_garnet / (Fe_used/Mg)_clinopyroxene
%   XCa_g = Ca_g / (Fe_used_g + Mg_g + Mn_g + Ca_g)
%
% The preferred equation is equation 17 on p. 242. Pressure is supplied in
% kbar and is used directly in the equation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Powell1985(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet/Grt and
%                    Clinopyroxene/Cpx tables
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Grt-Cpx pair
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('Powell1985 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || ~isreal(P_kbar) || isempty(P_kbar) || ...
        ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the accepted Garnet/Grt and Clinopyroxene/Cpx table aliases.
disp('=== Step 1: Preparing cation dataset ===');

if isfield(rawdata_struct, 'Garnet') && istable(rawdata_struct.Garnet)
    dataset_gt = rawdata_struct.Garnet;
elseif isfield(rawdata_struct, 'Grt') && istable(rawdata_struct.Grt)
    dataset_gt = rawdata_struct.Grt;
else
    error(['rawdata_struct must contain garnet table as either ' ...
           'rawdata_struct.Garnet or rawdata_struct.Grt']);
end

if isfield(rawdata_struct, 'Clinopyroxene') && ...
        istable(rawdata_struct.Clinopyroxene)
    dataset_cpx = rawdata_struct.Clinopyroxene;
elseif isfield(rawdata_struct, 'Cpx') && istable(rawdata_struct.Cpx)
    dataset_cpx = rawdata_struct.Cpx;
else
    error(['rawdata_struct must contain clinopyroxene table as either ' ...
           'rawdata_struct.Clinopyroxene or rawdata_struct.Cpx']);
end

requireVariables(dataset_gt, ...
    {'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu'}, 'Garnet');
requireVariables(dataset_cpx, ...
    {'Fe_cation_apfu', 'Mg_cation_apfu'}, 'Clinopyroxene');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Result table blocks are buffered and concatenated once after all selections,
% avoiding repeated full-table reallocation on every loop iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degC = 700;
calibrationT_max_degC = 1500;
directCalibrationP_kbar = 30;
calibrationXCa_min = 0.056;
calibrationXCa_max = 0.459;

pressureOutsideDirectCalibration = P_kbar ~= directCalibrationP_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
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
    disp(['Garnet selected: ' char(string(selectedCode_gt))]);

    % ----- Clinopyroxene selection -----
    disp('=== Step 4: Selecting a data code from the list (Clinopyroxene) ===');

    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', ...
        'Please select the Clinopyroxene data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIdx_cpx);
    disp(['Clinopyroxene selected: ' char(string(selectedCode_cpx))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_gt = dataset_gt(selectedIdx_gt, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    % Negative and infinite cation values are invalid. Zero and NaN are
    % accepted; mathematically undefined zero ratios return NaN.
    validateNonnegativeInputs(selectedData_gt, selectedData_cpx);

    % Missing or NaN thermometer inputs are recorded for a warning. They are
    % not replaced and do not stop the calculation.
    missingOrNaNInputNames = ...
        findMissingOrNaNInputs(selectedData_gt, selectedData_cpx);

    row = calcTemp(selectedData_gt, selectedData_cpx, P_kbar);

    row.dataCode_gt = repmat(string(selectedCode_gt), height(row), 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_gt', 'dataCode_cpx'}, 'Before', 1);

    % Store one table block per selected pair. Capacity is expanded
    % geometrically only when necessary instead of resizing results each loop.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature immediately, including NaN if present.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_gt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_gt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    % The preferred regression is directly calibrated at 30 kbar. The paper
    % supplies no defensible lower/upper pressure limits, so every pressure
    % differing from 30 kbar is identified as a pressure extrapolation.
    if any(pressureOutsideDirectCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure differs from the directly calibrated ' ...
             'pressure of Powell (1985): 30 kbar (pp. 235-236, 241-242). ' ...
             '%d of %d pressure point(s) differ from 30 kbar; input range = ' ...
             '%.4g-%.4g kbar. The pressure coefficient in equation 17 was ' ...
             'adopted rather than well constrained by regression, and the ' ...
             'paper gives no numerical lower/upper pressure application range.\n'], ...
            sum(pressureOutsideDirectCalibration), numel(P_kbar), ...
            min(P_kbar), max(P_kbar));
        pressureWarningIssued = true;
    end

    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the experimental ' ...
             'calibration-data range of Powell (1985): 700-1500 degreeC ' ...
             '(Table 3, pp. 238-239). %d of %d finite temperature point(s) ' ...
             'are outside the range; calculated finite range = %.4g-%.4g ' ...
             'degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_gt)), char(string(selectedCode_cpx)));
    end

    xCaValue = row.XCa_g(1);
    if isfinite(xCaValue) && ...
            (xCaValue < calibrationXCa_min || ...
             xCaValue > calibrationXCa_max)
        fprintf(2, ...
            ['WARNING: Garnet XCa is outside the calibration-data range of ' ...
             'Powell (1985): 0.056-0.459 (Table 3, pp. 238-239). ' ...
             'Calculated XCa_g = %.4g for %s & %s; the result is a ' ...
             'compositional extrapolation.\n'], ...
            xCaValue, char(string(selectedCode_gt)), ...
            char(string(selectedCode_cpx)));
    end

    if ~isempty(missingOrNaNInputNames)
        fprintf(2, ...
            ['WARNING: Missing or NaN thermometer input(s) were found for ' ...
             '%s & %s: %s.\n' ...
             '         NaN was retained and propagated; it was not treated as zero.\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_cpx)), ...
            char(strjoin(missingOrNaNInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), numel(row.T_deg), ...
            sum(isnan(row.T_deg)), sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Powell1985', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered table blocks once after the selection loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function requireVariables(tbl, requiredVariables, mineralLabel)
% requireVariables
% Confirm that the table contains every required temperature input variable.

missingVariables = strings(numel(requiredVariables), 1);
nMissing = 0;

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    if ~ismember(variableName, tbl.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingVariables(nMissing) = string(variableName);
    end
end

if nMissing > 0
    missingVariables = missingVariables(1:nMissing);
    error('%s table is missing required variable(s): %s.', ...
        mineralLabel, char(strjoin(missingVariables, ', ')));
end

end


function inputNames = findMissingOrNaNInputs(data_garnet, data_clinopyroxene)
% findMissingOrNaNInputs
% Return names of missing or NaN variables used directly by the thermometer.
% The output buffer is preallocated so its size does not grow in the loop.

garnetVariables = {'Fe_cation_apfu', 'Fe3_cation_apfu', ...
    'Mg_cation_apfu', 'Ca_cation_apfu', 'Mn_cation_apfu'};
clinopyroxeneVariables = {'Fe_cation_apfu', 'Fe3_cation_apfu', ...
    'Mg_cation_apfu'};

maximumNames = numel(garnetVariables) + numel(clinopyroxeneVariables);
inputNames = strings(maximumNames, 1);
nNames = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ~ismember(variableName, data_garnet.Properties.VariableNames)
        % Missing Fe3 and Mn retain the legacy default of zero. Other
        % thermometer variables are required and were checked earlier.
        if ~ismember(variableName, ...
                {'Fe3_cation_apfu', 'Mn_cation_apfu'})
            nNames = nNames + 1;
            inputNames(nNames) = ...
                "Garnet." + string(variableName) + " (missing)";
        end
    else
        variableValue = data_garnet.(variableName);
        if any(isnan(variableValue(:)))
            nNames = nNames + 1;
            inputNames(nNames) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(clinopyroxeneVariables)
    variableName = clinopyroxeneVariables{i};
    if ~ismember(variableName, data_clinopyroxene.Properties.VariableNames)
        % A missing Fe3 column retains the legacy default of zero.
        if ~strcmp(variableName, 'Fe3_cation_apfu')
            nNames = nNames + 1;
            inputNames(nNames) = ...
                "Clinopyroxene." + string(variableName) + " (missing)";
        end
    else
        variableValue = data_clinopyroxene.(variableName);
        if any(isnan(variableValue(:)))
            nNames = nNames + 1;
            inputNames(nNames) = ...
                "Clinopyroxene." + string(variableName);
        end
    end
end

inputNames = inputNames(1:nNames);

end


function validateNonnegativeInputs(data_garnet, data_clinopyroxene)
% validateNonnegativeInputs
% Stop when a supplied cation value is negative or Inf. Zero and NaN are
% accepted. NaN in thermometer variables propagates into the temperature.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu', 'Mn_cation_apfu', 'Ca_cation_apfu', ...
    'Si_cation_apfu', 'Al_cation_apfu'};
clinopyroxeneVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu', 'Mn_cation_apfu', 'Ca_cation_apfu', ...
    'Na_cation_apfu', 'K_cation_apfu', 'Ti_cation_apfu', ...
    'Si_cation_apfu', 'Al_cation_apfu'};

maximumNames = numel(garnetVariables) + numel(clinopyroxeneVariables);
invalidInputNames = strings(maximumNames, 1);
nInvalid = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if ~isnumeric(variableValue) || ~isreal(variableValue) || ...
                any(isinf(variableValue(:))) || ...
                any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nInvalid = nInvalid + 1;
            invalidInputNames(nInvalid) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(clinopyroxeneVariables)
    variableName = clinopyroxeneVariables{i};
    if ismember(variableName, data_clinopyroxene.Properties.VariableNames)
        variableValue = data_clinopyroxene.(variableName);
        if ~isnumeric(variableValue) || ~isreal(variableValue) || ...
                any(isinf(variableValue(:))) || ...
                any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nInvalid = nInvalid + 1;
            invalidInputNames(nInvalid) = ...
                "Clinopyroxene." + string(variableName);
        end
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNames(1:nInvalid);
    error(['Powell1985: cation inputs must be finite or NaN and must ' ...
           'be >= 0. Negative or infinite value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function row = calcTemp(data_garnet, data_clinopyroxene, P_kbar)
% calcTemp
% Compute temperature for one garnet-clinopyroxene pair and return one table
% row per supplied pressure value.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

gt = prepareGarnetRow(data_garnet);
cpx = prepareClinopyroxeneRow(data_clinopyroxene);
param = calcPowell1985Params(gt, cpx);

if isfinite(param.KD) && param.KD > 0
    lnKD = log(param.KD);
else
    lnKD = NaN;
end

denominator = 1.735 + lnKD;
canCalculate = isfinite(denominator) && denominator > 0 && ...
    isfinite(param.XCa_g);

if canCalculate
    T_K = (2790 + 10 .* P_kbar + 3140 .* param.XCa_g) ./ denominator;
    T_deg = T_K - 273.15;
else
    T_K = NaN(nP, 1);
    T_deg = NaN(nP, 1);
end

if isfinite(gt.Fe_used) && isfinite(gt.Mg) && gt.Fe_used + gt.Mg > 0
    Mg_number_gt = gt.Mg ./ (gt.Fe_used + gt.Mg);
else
    Mg_number_gt = NaN;
end

if isfinite(cpx.Fe_used) && isfinite(cpx.Mg) && cpx.Fe_used + cpx.Mg > 0
    Mg_number_cpx = cpx.Mg ./ (cpx.Fe_used + cpx.Mg);
else
    Mg_number_cpx = NaN;
end

if isfinite(cpx.Ca) && isfinite(cpx.Na) && cpx.Ca + cpx.Na > 0
    Ca_fraction_cpx = cpx.Ca ./ (cpx.Ca + cpx.Na);
else
    Ca_fraction_cpx = NaN;
end

is_positive_FeMg_garnet = isfinite(gt.Fe_used) && isfinite(gt.Mg) && ...
    gt.Fe_used > 0 && gt.Mg > 0;
is_positive_FeMg_clinopyroxene = ...
    isfinite(cpx.Fe_used) && isfinite(cpx.Mg) && ...
    cpx.Fe_used > 0 && cpx.Mg > 0;
is_positive_KD = isfinite(param.KD) && param.KD > 0;
is_positive_denominator = isfinite(denominator) && denominator > 0;
is_XCa_g_reasonable = isfinite(param.XCa_g) && ...
    param.XCa_g >= 0 && param.XCa_g <= 1;
is_XMn_g_reasonable = isfinite(param.XMn_g) && ...
    param.XMn_g >= 0 && param.XMn_g <= 1;
is_XCa_within_calibration = isfinite(param.XCa_g) && ...
    param.XCa_g >= 0.056 && param.XCa_g <= 0.459;

is_temperature_within_calibration = isfinite(T_deg) & ...
    T_deg >= 700 & T_deg <= 1500;
is_direct_calibration_pressure = P_kbar == 30;

% This numerical flag cannot determine whether the selected mineral domains
% were in equilibrium or whether cooling re-equilibration has occurred.
recommendedScalar = is_positive_FeMg_garnet && ...
    is_positive_FeMg_clinopyroxene && is_positive_KD && ...
    is_positive_denominator && is_XCa_within_calibration;
recommended = repmat(recommendedScalar, nP, 1) & ...
    is_temperature_within_calibration & is_direct_calibration_pressure;

row = table();
row.P_kbar = P_kbar;

% Repeat scalar composition and thermometer terms for every pressure row.
row.Fe2_g = repmat(gt.Fe_used, nP, 1);
row.Fe_raw_g = repmat(gt.Fe_raw, nP, 1);
row.Fe3_g = repmat(gt.Fe3, nP, 1);
row.Mg_g = repmat(gt.Mg, nP, 1);
row.Mn_g = repmat(gt.Mn, nP, 1);
row.Ca_g = repmat(gt.Ca, nP, 1);
row.Si_g = repmat(gt.Si, nP, 1);
row.Al_g = repmat(gt.Al, nP, 1);

row.Fe2_cpx = repmat(cpx.Fe_used, nP, 1);
row.Fe_raw_cpx = repmat(cpx.Fe_raw, nP, 1);
row.Fe3_cpx = repmat(cpx.Fe3, nP, 1);
row.Mg_cpx = repmat(cpx.Mg, nP, 1);
row.Mn_cpx = repmat(cpx.Mn, nP, 1);
row.Ca_cpx = repmat(cpx.Ca, nP, 1);
row.Na_cpx = repmat(cpx.Na, nP, 1);
row.K_cpx = repmat(cpx.K, nP, 1);
row.Ti_cpx = repmat(cpx.Ti, nP, 1);
row.Si_cpx = repmat(cpx.Si, nP, 1);
row.Al_cpx = repmat(cpx.Al, nP, 1);

row.FeMg_garnet = repmat(param.FeMg_garnet, nP, 1);
row.FeMg_clinopyroxene = repmat(param.FeMg_clinopyroxene, nP, 1);
row.KD = repmat(param.KD, nP, 1);
row.lnKD = repmat(lnKD, nP, 1);
row.XCa_g = repmat(param.XCa_g, nP, 1);
row.XMn_g = repmat(param.XMn_g, nP, 1);

row.Mg_number_garnet = repmat(Mg_number_gt, nP, 1);
row.Mg_number_clinopyroxene = repmat(Mg_number_cpx, nP, 1);
row.Ca_fraction_clinopyroxene_over_CaNa = ...
    repmat(Ca_fraction_cpx, nP, 1);

row.T_K = T_K;
row.T_deg = T_deg;

row.is_positive_FeMg_garnet = ...
    repmat(is_positive_FeMg_garnet, nP, 1);
row.is_positive_FeMg_clinopyroxene = ...
    repmat(is_positive_FeMg_clinopyroxene, nP, 1);
row.is_positive_KD = repmat(is_positive_KD, nP, 1);
row.is_positive_denominator = repmat(is_positive_denominator, nP, 1);
row.is_XCa_g_reasonable = repmat(is_XCa_g_reasonable, nP, 1);
row.is_XMn_g_reasonable = repmat(is_XMn_g_reasonable, nP, 1);
row.is_XCa_within_calibration = ...
    repmat(is_XCa_within_calibration, nP, 1);
row.is_temperature_within_calibration = ...
    is_temperature_within_calibration;
row.is_direct_calibration_pressure = is_direct_calibration_pressure;
row.recommended_by_Powell1985 = recommended;

end


function gt = prepareGarnetRow(data_garnet)
% prepareGarnetRow
% Read one garnet row while preserving all explicitly supplied NaN values.

if height(data_garnet) ~= 1
    error('Garnet input must be a 1-row table.');
end

gt = struct();
[gt.Fe_raw, gt.Fe3, gt.Fe_used] = getFeUsed(data_garnet, 'Garnet');
gt.Mg = getRequiredVar(data_garnet, 'Mg_cation_apfu', 'Garnet');
gt.Ca = getRequiredVar(data_garnet, 'Ca_cation_apfu', 'Garnet');

% Mn enters XCa_g and XMn_g. An absent Mn column retains the legacy
% assumption Mn = 0; an explicitly supplied NaN remains NaN.
gt.Mn = getOptionalVar(data_garnet, 'Mn_cation_apfu', 0, 'Garnet');

% Output-only optional variables are NaN when absent.
gt.Si = getOptionalVar(data_garnet, 'Si_cation_apfu', NaN, 'Garnet');
gt.Al = getOptionalVar(data_garnet, 'Al_cation_apfu', NaN, 'Garnet');

end


function cpx = prepareClinopyroxeneRow(data_clinopyroxene)
% prepareClinopyroxeneRow
% Read one clinopyroxene row while preserving explicitly supplied NaN values.

if height(data_clinopyroxene) ~= 1
    error('Clinopyroxene input must be a 1-row table.');
end

cpx = struct();
[cpx.Fe_raw, cpx.Fe3, cpx.Fe_used] = ...
    getFeUsed(data_clinopyroxene, 'Clinopyroxene');
cpx.Mg = getRequiredVar( ...
    data_clinopyroxene, 'Mg_cation_apfu', 'Clinopyroxene');

% Output-only optional variables are NaN when absent.
cpx.Mn = getOptionalVar( ...
    data_clinopyroxene, 'Mn_cation_apfu', NaN, 'Clinopyroxene');
cpx.Si = getOptionalVar( ...
    data_clinopyroxene, 'Si_cation_apfu', NaN, 'Clinopyroxene');
cpx.Al = getOptionalVar( ...
    data_clinopyroxene, 'Al_cation_apfu', NaN, 'Clinopyroxene');
cpx.Ti = getOptionalVar( ...
    data_clinopyroxene, 'Ti_cation_apfu', NaN, 'Clinopyroxene');
cpx.Ca = getOptionalVar( ...
    data_clinopyroxene, 'Ca_cation_apfu', NaN, 'Clinopyroxene');
cpx.Na = getOptionalVar( ...
    data_clinopyroxene, 'Na_cation_apfu', NaN, 'Clinopyroxene');
cpx.K = getOptionalVar( ...
    data_clinopyroxene, 'K_cation_apfu', NaN, 'Clinopyroxene');

end


function param = calcPowell1985Params(gt, cpx)
% calcPowell1985Params
% Calculate KD and garnet XCa used by equation 17 on p. 242.

param = struct();

garnetDivalentSum = gt.Fe_used + gt.Mg + gt.Mn + gt.Ca;
if isfinite(garnetDivalentSum) && garnetDivalentSum > 0
    param.XCa_g = gt.Ca ./ garnetDivalentSum;
    param.XMn_g = gt.Mn ./ garnetDivalentSum;
else
    param.XCa_g = NaN;
    param.XMn_g = NaN;
end

if isfinite(gt.Fe_used) && isfinite(gt.Mg) && gt.Mg > 0
    param.FeMg_garnet = gt.Fe_used ./ gt.Mg;
else
    param.FeMg_garnet = NaN;
end

if isfinite(cpx.Fe_used) && isfinite(cpx.Mg) && cpx.Mg > 0
    param.FeMg_clinopyroxene = cpx.Fe_used ./ cpx.Mg;
else
    param.FeMg_clinopyroxene = NaN;
end

if isfinite(param.FeMg_garnet) && ...
        isfinite(param.FeMg_clinopyroxene) && ...
        param.FeMg_clinopyroxene > 0
    param.KD = param.FeMg_garnet ./ param.FeMg_clinopyroxene;
else
    param.KD = NaN;
end

end


function [Fe_raw, Fe3, Fe_used] = getFeUsed(tbl, mineralLabel)
% getFeUsed
% Read Fe without replacing explicitly supplied NaN values. A missing Fe3
% column retains the legacy default Fe3 = 0.

Fe_raw = getRequiredVar(tbl, 'Fe_cation_apfu', mineralLabel);
Fe3 = getOptionalVar(tbl, 'Fe3_cation_apfu', 0, mineralLabel);
Fe_used = Fe_raw + Fe3;

end


function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read one required numeric scalar. NaN is allowed and preserved.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('Variable %s in %s table must be a numeric scalar.', ...
        varName, mineralLabel);
end
if isinf(value) || (isfinite(value) && value < 0)
    error('Variable %s in %s table must be NaN or a finite value >= 0.', ...
        varName, mineralLabel);
end

end


function value = getOptionalVar(tbl, varName, defaultValue, mineralLabel)
% getOptionalVar
% Read one optional numeric scalar. If the column is absent, use the supplied
% default. If the column exists and contains NaN, preserve NaN.

if ismember(varName, tbl.Properties.VariableNames)
    value = tbl.(varName);
else
    value = defaultValue;
end

if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('Variable %s in %s table must be a numeric scalar.', ...
        varName, mineralLabel);
end
if isinf(value) || (isfinite(value) && value < 0)
    error('Variable %s in %s table must be NaN or a finite value >= 0.', ...
        varName, mineralLabel);
end

end
