function results = Perchuk1985GrtAmph(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Amphibole/Perchuk1985GrtAmph.m
% Tested with MATLAB R2024b
%
% Garnet-amphibole Fe-Mg exchange geothermometer
% Perchuk, L.L. et al. (1985)
% Journal of Metamorphic Geology, 3, 265-310
% DOI: https://doi.org/10.1111/j.1525-1314.1985.tb00321.x
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Amphibole
% analysis (selected by the user from tables) and calculates temperature
% using the Perchuk et al. (1985) garnet-amphibole Fe-Mg exchange
% geothermometer.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Grt-Amp pair. One output row is produced
% for every supplied pressure value, allowing the same function to operate
% with both startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% The amphibole-garnet thermometer presented by Perchuk et al. (1985) is a
% slightly modified form of calibrations developed in earlier work. The first
% calibration was based on an amphibole-plagioclase thermometer; the 1000
% degreeC isotherm was subsequently constrained using amphibole-garnet pairs
% from peridotite inclusions in alkaline basalt, and some lower-temperature
% isotherms were modified (Appendix, p. 308; Fig. A1).
%
% Perchuk et al. (1985) explicitly state that Fe-Mg distribution between
% garnet and amphibole can be treated as ideal above 500 degreeC (p. 308).
% Fig. A1 extends to a 1000 degreeC isotherm. The paper does not define 1000
% degreeC as a strict upper applicability limit, but calculations above the
% highest constrained isotherm are extrapolations. This implementation
% therefore issues non-stopping warnings for finite temperatures:
%
%   Temperature <= 500 degreeC : outside the stated ideal-distribution
%                                condition
%   Temperature > 1000 degreeC : above the highest explicitly constrained
%                                isotherm discussed on p. 308
%
% A 470 degreeC garnet-amphibole rim temperature is reported in the main text
% as an example of retrograde re-equilibration (p. 295). It should not be
% interpreted as extending the calibration below the >500 degreeC condition.
%
% No numerical pressure calibration or application range is reported for this
% thermometer, and the equation contains no pressure term (p. 308). P_kbar is
% accepted for interface compatibility and traceability but is not used in the
% temperature equation. Because no source-based pressure limits exist, the
% function prints an fprintf notice instead of inventing a pressure cutoff.
%
% The calibration also assumes (p. 308):
%   1) garnet XCa differs very little from approximately 0.212,
%   2) Mn content in garnet is very low, and
%   3) amphibole composition varies over a wide range.
%
% The paper gives no numerical tolerance around XCa = 0.212 and no numerical
% upper limit for garnet Mn. Consequently, XCa_g and XMn_g are reported when
% the necessary data are available, but no arbitrary pass/fail threshold is
% imposed. Users must assess whether their garnet composition is sufficiently
% similar to the calibration compositions.
%
% Domain selection is critical. In the Sut-31 metabasite, core and rim
% compositions give approximately 750-470 degreeC, and the authors conclude
% that mineral thermometry mainly records retrograde Fe-Mg exchange (p. 295).
% Compare texturally corresponding, equilibrated domains (for example,
% core-core or contacting rim-contacting rim); calculated T may represent a
% closure or retrograde re-equilibration temperature rather than peak T.
%
% Relative errors for the four metabasite thermometers presented in the
% Appendix are stated not to exceed 10% (p. 310). This is a collective relative
% error estimate, not a separately defined 1-sigma or 95% uncertainty for the
% garnet-amphibole equation.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet    : table
% or
%   rawdata_struct.Grt       : table
%
% and
%
%   rawdata_struct.Amphibole : table
% or
%   rawdata_struct.Amp       : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns contain normalized
% cation data as apfu.
%
% Required variables in both Garnet and Amphibole tables:
%   Fe_cation_apfu          % Fe2+ used in the exchange equation
%   Mg_cation_apfu
%
% Optional variables in both tables:
%   Fe3_cation_apfu         % retained for output; not added to Fe2+
%   Mn_cation_apfu
%   Si_cation_apfu
%   Al_cation_apfu
%
% Optional variables in Garnet table:
%   Ca_cation_apfu
%
% Optional variables in Amphibole table:
%   Ca_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%   Ti_cation_apfu
%
% IMPORTANT Fe convention:
% The published exchange is between Fe2+ and Mg, and the average amphibole
% formula on p. 308 distinguishes Fe2+ from Fe3+. This implementation therefore
% uses Fe_cation_apfu as Fe2+ and does not add Fe3_cation_apfu. If the input
% Fe_cation_apfu contains total Fe rather than Fe2+, it must be converted to
% Fe2+ before this function is called.
%
% All finite cation values must be greater than or equal to zero. Negative
% values and Inf values stop the calculation. NaN values do not stop the
% calculation and are never converted to zero. NaN in Fe2+ or Mg propagates to
% the calculated temperature and generates a non-stopping fprintf warning.
% Missing optional variables and explicit NaN values are stored as NaN.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
%   T(K) = 3330 / (ln(KD) + 2.333)
%
% where
%
%   KD = (Fe2+/Mg)_garnet / (Fe2+/Mg)_amphibole
%
% The equation and its assumptions are given in the Appendix on p. 308.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Perchuk1985GrtAmph(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet/Grt and Amphibole/Amp tables
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Grt-Amp pair
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('Perchuk1985GrtAmph requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the accepted Garnet/Grt and Amphibole/Amp table aliases.
disp('=== Step 1: Preparing cation dataset ===');

if isfield(rawdata_struct, 'Garnet') && istable(rawdata_struct.Garnet)
    dataset_gt = rawdata_struct.Garnet;
elseif isfield(rawdata_struct, 'Grt') && istable(rawdata_struct.Grt)
    dataset_gt = rawdata_struct.Grt;
else
    error(['rawdata_struct must contain garnet table as either ' ...
           'rawdata_struct.Garnet or rawdata_struct.Grt']);
end

if isfield(rawdata_struct, 'Amphibole') && istable(rawdata_struct.Amphibole)
    dataset_amp = rawdata_struct.Amphibole;
elseif isfield(rawdata_struct, 'Amp') && istable(rawdata_struct.Amp)
    dataset_amp = rawdata_struct.Amp;
else
    error(['rawdata_struct must contain amphibole table as either ' ...
           'rawdata_struct.Amphibole or rawdata_struct.Amp']);
end

requireVariables(dataset_gt, ...
    {'Fe_cation_apfu', 'Mg_cation_apfu'}, 'Garnet');
requireVariables(dataset_amp, ...
    {'Fe_cation_apfu', 'Mg_cation_apfu'}, 'Amphibole');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Result table blocks are buffered and concatenated once after all selections,
% avoiding repeated full-table reallocation on every loop iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

idealDistributionT_min_degC = 500;
highestConstrainedT_degC = 1000;
pressureNoticeIssued = false;

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

    % ----- Amphibole selection -----
    disp('=== Step 4: Selecting a data code from the list (Amphibole) ===');

    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, ok] = listdlg( ...
        'PromptString', 'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_amp)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_gt = dataset_gt(selectedIdx_gt, :);
    selectedData_amp = dataset_amp(selectedIdx_amp, :);

    % NaN in required Fe2+ or Mg is recorded for a warning but is not replaced
    % and does not stop the calculation.
    nanInputNames = findNaNInputs(selectedData_gt, selectedData_amp);

    % Negative and infinite cation values are invalid. Zero and NaN are
    % accepted; mathematically undefined zero ratios return NaN.
    validateNonnegativeInputs(selectedData_gt, selectedData_amp);

    row = calcTemp(selectedData_gt, selectedData_amp, P_kbar);

    row.dataCode_gt = repmat(string(selectedCode_gt), height(row), 1);
    row.dataCode_amp = repmat(string(selectedCode_amp), height(row), 1);
    row = movevars(row, {'dataCode_gt', 'dataCode_amp'}, 'Before', 1);

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
        disp([char(string(selectedCode_gt)) ' & ' char(string(selectedCode_amp)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_gt)) ' & ' char(string(selectedCode_amp)) ...
            ': ' num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ' degreeC']);
    end

    % The paper gives no numerical pressure range. Print this once so that a
    % supplied pressure vector is not mistaken for a pressure correction.
    if ~pressureNoticeIssued
        fprintf(2, ...
            ['NOTICE: Perchuk et al. (1985) provide no numerical pressure ' ...
             'calibration or application range for the garnet-amphibole ' ...
             'thermometer, and the equation has no pressure term (p. 308). ' ...
             'Input pressure %.4g-%.4g kbar is stored for traceability but ' ...
             'is not used or range-validated.\n'], min(P_kbar), max(P_kbar));
        pressureNoticeIssued = true;
    end

    finiteTemperature = isfinite(row.T_deg);

    % At and below 500 degreeC the ideal-distribution assumption stated by the
    % authors is not satisfied.
    temperatureAtOrBelowIdealLimit = ...
        finiteTemperature & row.T_deg <= idealDistributionT_min_degC;
    if any(temperatureAtOrBelowIdealLimit)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is at or below 500 degreeC, ' ...
             'outside the condition under which Perchuk et al. (1985) state ' ...
             'that garnet-amphibole Fe-Mg distribution can be treated as ' ...
             'ideal (p. 308). %d of %d finite temperature point(s) are ' ...
             'outside this condition; calculated finite range = %.4g-%.4g ' ...
             'degreeC for %s & %s.\n'], ...
            sum(temperatureAtOrBelowIdealLimit), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)));
    end

    % Above 1000 degreeC the calculation lies beyond the highest isotherm
    % explicitly constrained and discussed by the authors.
    temperatureAboveHighestIsotherm = ...
        finiteTemperature & row.T_deg > highestConstrainedT_degC;
    if any(temperatureAboveHighestIsotherm)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is above 1000 degreeC, the ' ...
             'highest explicitly constrained isotherm discussed for the ' ...
             'Perchuk et al. (1985) garnet-amphibole thermometer (p. 308). ' ...
             '%d of %d finite temperature point(s) are extrapolated; ' ...
             'calculated finite range = %.4g-%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureAboveHighestIsotherm), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)));
    end

    % The paper gives composition targets but no numerical tolerances. Report
    % the calculated values without inventing an automatic pass/fail boundary.
    xCaValue = row.XCa_g(1);
    xMnValue = row.XMn_g(1);
    if isfinite(xCaValue) && isfinite(xMnValue)
        fprintf(2, ...
            ['CAUTION: Perchuk et al. (1985) calibrated the thermometer for ' ...
             'garnet XCa close to 0.212 and very low garnet Mn (p. 308). ' ...
             'For %s & %s, calculated XCa_g = %.4g and XMn_g = %.4g. ' ...
             'No numerical tolerances are provided in the paper; assess ' ...
             'compositional applicability manually.\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)), ...
            xCaValue, xMnValue);
    else
        fprintf(2, ...
            ['CAUTION: Perchuk et al. (1985) calibrated the thermometer for ' ...
             'garnet XCa close to 0.212 and very low garnet Mn (p. 308), but ' ...
             'XCa_g and/or XMn_g could not be evaluated for %s & %s because ' ...
             'Ca and/or Mn data are missing or NaN.\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; it was not treated as zero.\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)), ...
            sum(invalidTemperature), numel(row.T_deg), ...
            sum(isnan(row.T_deg)), sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Perchuk1985GrtAmph', ...
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


function nanInputNames = findNaNInputs(data_garnet, data_amphibole)
% findNaNInputs
% Return required Fe2+/Mg input names containing NaN. The output buffer is
% preallocated so its size does not grow on each loop iteration.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
amphiboleVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

maximumNames = numel(garnetVariables) + numel(amphiboleVariables);
nanInputNames = strings(maximumNames, 1);
nNames = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    variableValue = data_garnet.(variableName);
    if any(isnan(variableValue(:)))
        nNames = nNames + 1;
        nanInputNames(nNames) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(amphiboleVariables)
    variableName = amphiboleVariables{i};
    variableValue = data_amphibole.(variableName);
    if any(isnan(variableValue(:)))
        nNames = nNames + 1;
        nanInputNames(nNames) = "Amphibole." + string(variableName);
    end
end

nanInputNames = nanInputNames(1:nNames);

end


function validateNonnegativeInputs(data_garnet, data_amphibole)
% validateNonnegativeInputs
% Stop when a supplied cation value is negative or Inf. Zero and NaN are
% accepted. NaN in required variables propagates into the temperature.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu', 'Mn_cation_apfu', 'Ca_cation_apfu', ...
    'Si_cation_apfu', 'Al_cation_apfu'};
amphiboleVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu', 'Mn_cation_apfu', 'Ca_cation_apfu', ...
    'Na_cation_apfu', 'K_cation_apfu', 'Ti_cation_apfu', ...
    'Si_cation_apfu', 'Al_cation_apfu'};

maximumNames = numel(garnetVariables) + numel(amphiboleVariables);
invalidInputNames = strings(maximumNames, 1);
nInvalid = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if ~isnumeric(variableValue) || any(isinf(variableValue(:))) || ...
                any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nInvalid = nInvalid + 1;
            invalidInputNames(nInvalid) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(amphiboleVariables)
    variableName = amphiboleVariables{i};
    if ismember(variableName, data_amphibole.Properties.VariableNames)
        variableValue = data_amphibole.(variableName);
        if ~isnumeric(variableValue) || any(isinf(variableValue(:))) || ...
                any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nInvalid = nInvalid + 1;
            invalidInputNames(nInvalid) = "Amphibole." + string(variableName);
        end
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNames(1:nInvalid);
    error(['Perchuk1985GrtAmph: cation inputs must be finite or NaN ' ...
           'and must be >= 0. Negative or infinite value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function row = calcTemp(data_garnet, data_amphibole, P_kbar)
% calcTemp
% Compute temperature for one garnet-amphibole pair and return one table row
% per supplied pressure. Pressure is stored but not used by the equation.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

gt = prepareGarnetRow(data_garnet);
amp = prepareAmphiboleRow(data_amphibole);
param = calcPerchuk1985Params(gt, amp);

if isfinite(param.KD) && param.KD > 0
    lnKD = log(param.KD);
else
    lnKD = NaN;
end

denominator = lnKD + 2.333;
canCalculate = isfinite(denominator) && denominator > 0;

if canCalculate
    T_K_scalar = 3330 ./ denominator;
    T_deg_scalar = T_K_scalar - 273.15;
else
    T_K_scalar = NaN;
    T_deg_scalar = NaN;
end

if isfinite(gt.Fe2) && isfinite(gt.Mg) && gt.Fe2 + gt.Mg > 0
    Mg_number_gt = gt.Mg ./ (gt.Fe2 + gt.Mg);
else
    Mg_number_gt = NaN;
end

if isfinite(amp.Fe2) && isfinite(amp.Mg) && amp.Fe2 + amp.Mg > 0
    Mg_number_amp = amp.Mg ./ (amp.Fe2 + amp.Mg);
else
    Mg_number_amp = NaN;
end

garnetDivalentSum = gt.Fe2 + gt.Mg + gt.Mn + gt.Ca;
if isfinite(garnetDivalentSum) && garnetDivalentSum > 0
    XCa_g = gt.Ca ./ garnetDivalentSum;
    XMn_g = gt.Mn ./ garnetDivalentSum;
else
    XCa_g = NaN;
    XMn_g = NaN;
end

if isfinite(amp.Ca) && isfinite(amp.Na) && amp.Ca + amp.Na > 0
    Ca_fraction_amp = amp.Ca ./ (amp.Ca + amp.Na);
else
    Ca_fraction_amp = NaN;
end

is_positive_FeMg_garnet = isfinite(gt.Fe2) && isfinite(gt.Mg) && ...
    gt.Fe2 > 0 && gt.Mg > 0;
is_positive_FeMg_amphibole = isfinite(amp.Fe2) && isfinite(amp.Mg) && ...
    amp.Fe2 > 0 && amp.Mg > 0;
is_positive_KD = isfinite(param.KD) && param.KD > 0;
is_positive_denominator = isfinite(denominator) && denominator > 0;
is_above_500C_scalar = isfinite(T_deg_scalar) && T_deg_scalar > 500;
is_at_or_below_1000C_scalar = isfinite(T_deg_scalar) && T_deg_scalar <= 1000;

% This numerical flag cannot evaluate whether XCa is sufficiently close to
% 0.212, Mn is sufficiently low, or the selected domains were equilibrated.
recommendedScalar = is_positive_FeMg_garnet && ...
    is_positive_FeMg_amphibole && is_positive_KD && ...
    is_positive_denominator && is_above_500C_scalar && ...
    is_at_or_below_1000C_scalar;

row = table();
row.P_kbar = P_kbar;

% Repeat scalar composition and thermometer terms for every pressure row.
row.Fe2_g = repmat(gt.Fe2, nP, 1);
row.Fe_raw_g = repmat(gt.Fe2, nP, 1);
row.Fe3_g = repmat(gt.Fe3, nP, 1);
row.Mg_g = repmat(gt.Mg, nP, 1);
row.Mn_g = repmat(gt.Mn, nP, 1);
row.Ca_g = repmat(gt.Ca, nP, 1);
row.Si_g = repmat(gt.Si, nP, 1);
row.Al_g = repmat(gt.Al, nP, 1);

row.Fe2_amp = repmat(amp.Fe2, nP, 1);
row.Fe_raw_amp = repmat(amp.Fe2, nP, 1);
row.Fe3_amp = repmat(amp.Fe3, nP, 1);
row.Mg_amp = repmat(amp.Mg, nP, 1);
row.Mn_amp = repmat(amp.Mn, nP, 1);
row.Ca_amp = repmat(amp.Ca, nP, 1);
row.Na_amp = repmat(amp.Na, nP, 1);
row.K_amp = repmat(amp.K, nP, 1);
row.Ti_amp = repmat(amp.Ti, nP, 1);
row.Si_amp = repmat(amp.Si, nP, 1);
row.Al_amp = repmat(amp.Al, nP, 1);

row.FeMg_garnet = repmat(param.FeMg_garnet, nP, 1);
row.FeMg_amphibole = repmat(param.FeMg_amphibole, nP, 1);
row.KD = repmat(param.KD, nP, 1);
row.lnKD = repmat(lnKD, nP, 1);

row.Mg_number_garnet = repmat(Mg_number_gt, nP, 1);
row.Mg_number_amphibole = repmat(Mg_number_amp, nP, 1);
row.XCa_g = repmat(XCa_g, nP, 1);
row.XMn_g = repmat(XMn_g, nP, 1);
row.Ca_fraction_amphibole_over_CaNa = repmat(Ca_fraction_amp, nP, 1);

row.T_K = repmat(T_K_scalar, nP, 1);
row.T_deg = repmat(T_deg_scalar, nP, 1);

row.is_positive_FeMg_garnet = repmat(is_positive_FeMg_garnet, nP, 1);
row.is_positive_FeMg_amphibole = repmat(is_positive_FeMg_amphibole, nP, 1);
row.is_positive_KD = repmat(is_positive_KD, nP, 1);
row.is_positive_denominator = repmat(is_positive_denominator, nP, 1);
row.is_above_500C = repmat(is_above_500C_scalar, nP, 1);
row.is_at_or_below_1000C = repmat(is_at_or_below_1000C_scalar, nP, 1);
row.is_within_figure_temperature_interval = ...
    row.is_above_500C & row.is_at_or_below_1000C;
row.recommended_by_Perchuk1985 = repmat(recommendedScalar, nP, 1);

end


function gt = prepareGarnetRow(data_garnet)
% prepareGarnetRow
% Read one garnet row while preserving all explicitly supplied NaN values.

if height(data_garnet) ~= 1
    error('Garnet input must be a 1-row table.');
end

gt = struct();
gt.Fe2 = getRequiredVar(data_garnet, 'Fe_cation_apfu', 'Garnet');
gt.Mg = getRequiredVar(data_garnet, 'Mg_cation_apfu', 'Garnet');
gt.Fe3 = getOptionalVar(data_garnet, 'Fe3_cation_apfu', 'Garnet');
gt.Mn = getOptionalVar(data_garnet, 'Mn_cation_apfu', 'Garnet');
gt.Ca = getOptionalVar(data_garnet, 'Ca_cation_apfu', 'Garnet');
gt.Si = getOptionalVar(data_garnet, 'Si_cation_apfu', 'Garnet');
gt.Al = getOptionalVar(data_garnet, 'Al_cation_apfu', 'Garnet');

end


function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Read one amphibole row while preserving all explicitly supplied NaN values.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();
amp.Fe2 = getRequiredVar(data_amphibole, 'Fe_cation_apfu', 'Amphibole');
amp.Mg = getRequiredVar(data_amphibole, 'Mg_cation_apfu', 'Amphibole');
amp.Fe3 = getOptionalVar(data_amphibole, 'Fe3_cation_apfu', 'Amphibole');
amp.Mn = getOptionalVar(data_amphibole, 'Mn_cation_apfu', 'Amphibole');
amp.Si = getOptionalVar(data_amphibole, 'Si_cation_apfu', 'Amphibole');
amp.Al = getOptionalVar(data_amphibole, 'Al_cation_apfu', 'Amphibole');
amp.Ti = getOptionalVar(data_amphibole, 'Ti_cation_apfu', 'Amphibole');
amp.Ca = getOptionalVar(data_amphibole, 'Ca_cation_apfu', 'Amphibole');
amp.Na = getOptionalVar(data_amphibole, 'Na_cation_apfu', 'Amphibole');
amp.K = getOptionalVar(data_amphibole, 'K_cation_apfu', 'Amphibole');

end


function param = calcPerchuk1985Params(gt, amp)
% calcPerchuk1985Params
% Calculate the Fe2+-Mg distribution coefficient used on p. 308.

param = struct();

if isfinite(gt.Fe2) && isfinite(gt.Mg) && gt.Mg > 0
    param.FeMg_garnet = gt.Fe2 ./ gt.Mg;
else
    param.FeMg_garnet = NaN;
end

if isfinite(amp.Fe2) && isfinite(amp.Mg) && amp.Mg > 0
    param.FeMg_amphibole = amp.Fe2 ./ amp.Mg;
else
    param.FeMg_amphibole = NaN;
end

if isfinite(param.FeMg_garnet) && ...
        isfinite(param.FeMg_amphibole) && param.FeMg_amphibole > 0
    param.KD = param.FeMg_garnet ./ param.FeMg_amphibole;
else
    param.KD = NaN;
end

end


function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read one required numeric scalar. NaN is allowed and preserved.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);
if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s in %s table must be a numeric scalar.', ...
        varName, mineralLabel);
end
if isinf(value) || (isfinite(value) && value < 0)
    error('Variable %s in %s table must be NaN or a finite value >= 0.', ...
        varName, mineralLabel);
end

end


function value = getOptionalVar(tbl, varName, mineralLabel)
% getOptionalVar
% Read one optional numeric scalar. Missing variables and explicitly supplied
% NaN values are returned as NaN and are never converted to zero.

if ismember(varName, tbl.Properties.VariableNames)
    value = tbl.(varName);
else
    value = NaN;
end

if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s in %s table must be a numeric scalar.', ...
        varName, mineralLabel);
end
if isinf(value) || (isfinite(value) && value < 0)
    error('Variable %s in %s table must be NaN or a finite value >= 0.', ...
        varName, mineralLabel);
end

end
