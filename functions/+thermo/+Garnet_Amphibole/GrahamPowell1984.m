function results = GrahamPowell1984(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Amphibole/GrahamPowell1984.m
% Tested with MATLAB R2024b
%
% Garnet-hornblende Fe-Mg exchange geothermometer
% Graham, C.M. and Powell, R. (1984)
% Journal of Metamorphic Geology, 2, 13-31
% DOI: https://doi.org/10.1111/j.1525-1314.1984.tb00282.x
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Amphibole
% analysis (selected by the user from tables) and calculates temperature
% using the Graham and Powell (1984) garnet-hornblende Fe-Mg exchange
% geothermometer.
%
% The original calibration is for garnet + common hornblende pairs. In this
% implementation, compositions are read from an Amphibole (or Amp) table and
% used as the hornblende-equivalent phase. The user must independently verify
% that the selected amphibole is a common hornblende and that the selected
% garnet and hornblende attained Fe-Mg exchange equilibrium.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Grt-Amp pair. One output row is produced
% for every supplied pressure value, allowing the same function to operate
% with both startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Graham and Powell (1984) calibrated the thermometer empirically against
% garnet-clinopyroxene temperatures for garnet + hornblende + clinopyroxene
% assemblages. Table 1 (p. 16), together with the outlier discussion on
% pp. 18-19, indicates the following approximate span of the retained primary
% calibration data:
%
%   Temperature : approximately 599-920 degreeC
%                 (natural calibrants are mainly approximately 599-815
%                  degreeC; the high-temperature end includes an experimental
%                  datum)
%   Pressure    : approximately 5-18 kbar
%   XCa_g       : approximately 0.165-0.418
%
% These numerical spans describe the retained calibration data; they are not
% all stated as strict hard limits by the authors. In particular, the fitted
% pressure coefficient was indistinguishable from zero, so the pressure term
% was omitted from the final equation (p. 17). P_kbar is therefore retained
% for interface compatibility and traceability but is not used to calculate T.
%
% The authors explicitly recommend application under the following conditions:
%   1) Temperature below approximately 850 degreeC because higher-temperature
%      natural assemblages may record Fe-Mg closure rather than the thermal
%      maximum (p. 20; conclusion on p. 28).
%   2) Mn-poor garnet, in practice XMn_g < 0.1 (p. 20; p. 28).
%   3) Rocks metamorphosed at low oxygen activity because unusually ferric
%      hornblende may yield overestimated temperatures (pp. 16-17, 20, 28).
%   4) Common hornblende compositions. Amphiboles with Na(M4) < 0.09 or high
%      Ca(M4) were excluded from calibration, and the upper Na(M4) limit was
%      uncertain (p. 17).
%   5) Texturally and chemically equilibrated garnet-hornblende pairs.
%      Eclogitic garnet amphibolites commonly failed this requirement and
%      should be treated with particular caution (p. 22; p. 28).
%
% Independent tests reported approximately 520-610 degreeC for Dalradian
% garnet amphibolites (pp. 20-21), and the Pelona Schist application extended
% to approximately 480-650 degreeC (p. 28). Therefore, approximately 599
% degreeC is the lower end of the retained primary calibration data, not an
% explicit lower applicability limit proposed by the authors.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside approximately 5-18 kbar,
%   2) a finite calculated temperature is below approximately 599 degreeC,
%   3) a finite calculated temperature is at or above 850 degreeC,
%   4) XCa_g is outside approximately 0.165-0.418,
%   5) XMn_g is at or above 0.1,
%   6) a calculation input or calculated temperature is NaN/Inf.
%
% The numerical flags and warnings cannot determine oxygen activity, mineral
% identity, textural equilibrium, or whether amphibole formed later than
% garnet. These conditions must be evaluated petrographically and chemically.
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
% Required Garnet variables:
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%
% Required Amphibole variables:
%   Fe_cation_apfu
%   Mg_cation_apfu
%
% Optional calculation variables:
%   Garnet and Amphibole: Fe3_cation_apfu
%   Garnet              : Mn_cation_apfu
%
% If an optional calculation variable is absent, it is assumed to be zero to
% preserve compatibility with the original input format. If the variable is
% present but its value is NaN, NaN is retained and propagated; it is never
% converted to zero.
%
% Optional variables retained for output:
%   Garnet    : Si_cation_apfu, Al_cation_apfu
%   Amphibole : Mn_cation_apfu, Si_cation_apfu, Al_cation_apfu,
%               Ti_cation_apfu, Ca_cation_apfu, Na_cation_apfu,
%               K_cation_apfu
%
% Missing optional output-only variables are stored as NaN. All finite cation
% values must be greater than or equal to zero. Negative values and Inf values
% stop the calculation. NaN values do not stop the calculation; they propagate
% into the result and generate non-stopping fprintf warnings.
%
% IMPORTANT Fe convention:
% Graham and Powell (1984) treated all hornblende Fe as Fe2+ because reliable
% Fe3+/Fe2+ recalculation was unavailable (pp. 16-17). This implementation uses
%
%   Fe_used = Fe_cation_apfu + Fe3_cation_apfu
%
% and treats Fe_used as Fe2+ in the exchange expression. This is appropriate
% only when Fe_cation_apfu and Fe3_cation_apfu are separate quantities. If
% Fe_cation_apfu already contains total Fe, Fe3_cation_apfu must not be added
% again.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
%   T(K) = (2880 + 3280 * XCa_g) / (ln(KD) + 2.426)
%
% where
%
%   KD    = (Fe/Mg)_garnet / (Fe/Mg)_hornblende
%   XCa_g = Ca_g / (Fe_g + Mg_g + Mn_g + Ca_g)
%   XMn_g = Mn_g / (Fe_g + Mg_g + Mn_g + Ca_g)
%
% The equation is presented on p. 20 of Graham and Powell (1984).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = GrahamPowell1984(rawdata_struct, P_kbar)
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
    error('GrahamPowell1984 requires (rawdata_struct, P_kbar).');
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
    {'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu'}, 'Garnet');
requireVariables(dataset_amp, ...
    {'Fe_cation_apfu', 'Mg_cation_apfu'}, 'Amphibole');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Repeated concatenation of the full results table inside the interactive
% loop is avoided. Result table blocks are buffered and concatenated once
% after all selections have been completed.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Approximate retained primary-calibration spans and the recommended upper
% application limit reported by Graham and Powell (1984).
calibrationT_min_degC = 599;
recommendedT_max_degC = 850;
calibrationP_min_kbar = 5;
calibrationP_max_kbar = 18;
calibrationXCa_min = 0.165;
calibrationXCa_max = 0.418;
recommendedXMn_max = 0.1;

% Pressure is common to all selected mineral pairs in this function call,
% so the pressure warning is printed only once after the first calculation.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
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

    % NaN is recorded for a warning but does not stop the calculation.
    nanInputNames = findNaNInputs(selectedData_gt, selectedData_amp);

    % Negative and infinite cation values are invalid. Zero is accepted;
    % invalid zero denominators are converted to NaN by calcTemp.
    validateNonnegativeInputs(selectedData_gt, selectedData_amp);

    row = calcTemp(selectedData_gt, selectedData_amp, P_kbar);

    % Store the selected identifiers for every pressure row.
    row.dataCode_gt = repmat(string(selectedCode_gt), height(row), 1);
    row.dataCode_amp = repmat(string(selectedCode_amp), height(row), 1);
    row = movevars(row, {'dataCode_gt', 'dataCode_amp'}, 'Before', 1);

    % Buffer the result block. Capacity is expanded geometrically only when
    % needed, avoiding full table reallocation on every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature immediately.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_gt)) ' & ' char(string(selectedCode_amp)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_gt)) ' & ' char(string(selectedCode_amp)) ...
            ': ' num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ' degreeC']);
    end

    % Warn once when pressure is outside the approximate span represented by
    % the retained primary calibration data. The published equation itself
    % has no pressure term, and calculation is not stopped.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate pressure span ' ...
             'of the retained Graham and Powell (1984) calibration data: ' ...
             '5-18 kbar. %d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar. The equation is pressure-independent, ' ...
             'so calculation was continued.\n'], ...
            sum(pressureOutsideCalibration), numel(P_kbar), ...
            min(P_kbar), max(P_kbar));
        pressureWarningIssued = true;
    end

    finiteTemperature = isfinite(row.T_deg);

    % The lower warning identifies extrapolation below the retained primary
    % calibration dataset. It is not a strict lower applicability cutoff.
    temperatureBelowPrimaryCalibration = ...
        finiteTemperature & row.T_deg < calibrationT_min_degC;
    if any(temperatureBelowPrimaryCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is below the approximate lower end ' ...
             'of the retained primary calibration data of Graham and Powell ' ...
             '(1984): approximately 599 degreeC. %d of %d finite temperature ' ...
             'point(s) are below this value; calculated finite range = ' ...
             '%.4g-%.4g degreeC for %s & %s. The authors reported tests and ' ...
             'applications down to approximately 480-520 degreeC, so this is ' ...
             'an extrapolation warning rather than a hard rejection.\n'], ...
            sum(temperatureBelowPrimaryCalibration), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)));
    end

    % The upper warning follows the authors' explicit recommendation that the
    % thermometer be used below approximately 850 degreeC.
    temperatureAtOrAboveRecommendedMaximum = ...
        finiteTemperature & row.T_deg >= recommendedT_max_degC;
    if any(temperatureAtOrAboveRecommendedMaximum)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is at or above the recommended ' ...
             'application limit of Graham and Powell (1984): approximately ' ...
             '850 degreeC. %d of %d finite temperature point(s) are outside ' ...
             'the recommended range; calculated finite range = %.4g-%.4g ' ...
             'degreeC for %s & %s. Fe-Mg closure effects may prevent recovery ' ...
             'of the thermal maximum.\n'], ...
            sum(temperatureAtOrAboveRecommendedMaximum), ...
            sum(finiteTemperature), min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)));
    end

    % Report compositional extrapolation without stopping calculation.
    finiteXCa = isfinite(row.XCa_g);
    xCaOutsideCalibration = finiteXCa & ...
        (row.XCa_g < calibrationXCa_min | row.XCa_g > calibrationXCa_max);
    if any(xCaOutsideCalibration)
        fprintf(2, ...
            ['WARNING: Garnet XCa is outside the approximate range represented ' ...
             'by the retained Graham and Powell (1984) calibration data: ' ...
             '0.165-0.418. Calculated XCa_g = %.4g for %s & %s.\n'], ...
            row.XCa_g(find(finiteXCa, 1, 'first')), ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)));
    end

    finiteXMn = isfinite(row.XMn_g);
    xMnOutsideRecommendation = finiteXMn & row.XMn_g >= recommendedXMn_max;
    if any(xMnOutsideRecommendation)
        fprintf(2, ...
            ['WARNING: Garnet is outside the Mn-poor application criterion of ' ...
             'Graham and Powell (1984): XMn_g < 0.1. Calculated XMn_g = %.4g ' ...
             'for %s & %s.\n'], ...
            row.XMn_g(find(finiteXMn, 1, 'first')), ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)));
    end

    % Print a non-stopping warning when any calculation input contains NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; it was not treated as zero.\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report any NaN/Inf temperature caused by an input NaN, a
    % zero denominator, or another mathematically undefined operation.
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
        'GrahamPowell1984', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered blocks once. Return an empty table if no calculation
% was performed.
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
% Confirm that the input table contains every required calculation variable.

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
% Return names of calculation variables containing NaN without throwing an
% error. The output buffer is preallocated to avoid growth inside the loop.

garnetVariables = {'Fe_cation_apfu', 'Fe3_cation_apfu', ...
    'Mg_cation_apfu', 'Mn_cation_apfu', 'Ca_cation_apfu'};
amphiboleVariables = {'Fe_cation_apfu', 'Fe3_cation_apfu', ...
    'Mg_cation_apfu'};

maximumNames = numel(garnetVariables) + numel(amphiboleVariables);
nanInputNames = strings(maximumNames, 1);
nNames = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if any(isnan(variableValue(:)))
            nNames = nNames + 1;
            nanInputNames(nNames) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(amphiboleVariables)
    variableName = amphiboleVariables{i};
    if ismember(variableName, data_amphibole.Properties.VariableNames)
        variableValue = data_amphibole.(variableName);
        if any(isnan(variableValue(:)))
            nNames = nNames + 1;
            nanInputNames(nNames) = "Amphibole." + string(variableName);
        end
    end
end

nanInputNames = nanInputNames(1:nNames);

end

function validateNonnegativeInputs(data_garnet, data_amphibole)
% validateNonnegativeInputs
% Stop when a finite calculation input is negative or when an input is Inf.
% Zero and NaN are allowed. NaN is propagated and reported by fprintf.

garnetVariables = {'Fe_cation_apfu', 'Fe3_cation_apfu', ...
    'Mg_cation_apfu', 'Mn_cation_apfu', 'Ca_cation_apfu'};
amphiboleVariables = {'Fe_cation_apfu', 'Fe3_cation_apfu', ...
    'Mg_cation_apfu'};

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
    error(['GrahamPowell1984: calculation inputs must be finite or NaN ' ...
           'and must be >= 0. Negative or infinite value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_garnet, data_amphibole, P_kbar)
% calcTemp
% Compute temperature for one garnet-amphibole pair and return one table row
% per supplied pressure. Pressure is stored but is not used by the equation.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

gt = prepareGarnetRow(data_garnet);
amp = prepareAmphiboleRow(data_amphibole);
param = calcGrahamPowell1984Params(gt, amp);

% Guard only mathematically undefined exchange expressions. Input NaN is
% retained in the structures and therefore produces a NaN temperature.
if isfinite(param.KD) && param.KD > 0
    lnKD = log(param.KD);
else
    lnKD = NaN;
end

denominator = lnKD + 2.426;
canCalculate = isfinite(param.XCa_g) && isfinite(denominator) && ...
    denominator > 0;

if canCalculate
    T_K_scalar = (2880 + 3280 .* param.XCa_g) ./ denominator;
    T_deg_scalar = T_K_scalar - 273.15;
else
    T_K_scalar = NaN;
    T_deg_scalar = NaN;
end

if isfinite(gt.Mg) && isfinite(gt.Fe_used) && gt.Mg + gt.Fe_used > 0
    Mg_number_gt = gt.Mg ./ (gt.Mg + gt.Fe_used);
else
    Mg_number_gt = NaN;
end

if isfinite(amp.Mg) && isfinite(amp.Fe_used) && amp.Mg + amp.Fe_used > 0
    Mg_number_amp = amp.Mg ./ (amp.Mg + amp.Fe_used);
else
    Mg_number_amp = NaN;
end

if isfinite(amp.Ca) && isfinite(amp.Na) && amp.Ca + amp.Na > 0
    Ca_fraction_amp = amp.Ca ./ (amp.Ca + amp.Na);
else
    Ca_fraction_amp = NaN;
end

is_positive_FeMg_garnet = isfinite(gt.Fe_used) && isfinite(gt.Mg) && ...
    gt.Fe_used > 0 && gt.Mg > 0;
is_positive_FeMg_amphibole = isfinite(amp.Fe_used) && isfinite(amp.Mg) && ...
    amp.Fe_used > 0 && amp.Mg > 0;
is_positive_KD = isfinite(param.KD) && param.KD > 0;
is_positive_denominator = isfinite(denominator) && denominator > 0;
is_XCa_g_reasonable = isfinite(param.XCa_g) && ...
    param.XCa_g >= 0 && param.XCa_g <= 1;
is_XMn_g_low = isfinite(param.XMn_g) && param.XMn_g < 0.1;
is_below_850C_scalar = isfinite(T_deg_scalar) && T_deg_scalar < 850;

% This legacy flag contains only numerical checks available to the function.
% It cannot assess low oxygen activity, hornblende identity, or equilibrium.
recommendedScalar = is_positive_FeMg_garnet && ...
    is_positive_FeMg_amphibole && is_positive_KD && ...
    is_positive_denominator && is_XCa_g_reasonable && ...
    is_XMn_g_low && is_below_850C_scalar;

row = table();
row.P_kbar = P_kbar;

% Repeat scalar composition terms so that every output variable has nP rows.
row.Fe2_g = repmat(gt.Fe_used, nP, 1);
row.Fe_raw_g = repmat(gt.Fe_raw, nP, 1);
row.Fe3_g = repmat(gt.Fe3, nP, 1);
row.Mg_g = repmat(gt.Mg, nP, 1);
row.Mn_g = repmat(gt.Mn, nP, 1);
row.Ca_g = repmat(gt.Ca, nP, 1);
row.Si_g = repmat(gt.Si, nP, 1);
row.Al_g = repmat(gt.Al, nP, 1);

row.Fe2_amp = repmat(amp.Fe_used, nP, 1);
row.Fe_raw_amp = repmat(amp.Fe_raw, nP, 1);
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
row.XCa_g = repmat(param.XCa_g, nP, 1);
row.XMn_g = repmat(param.XMn_g, nP, 1);

row.Mg_number_garnet = repmat(Mg_number_gt, nP, 1);
row.Mg_number_amphibole = repmat(Mg_number_amp, nP, 1);
row.Ca_fraction_amphibole_over_CaNa = repmat(Ca_fraction_amp, nP, 1);

row.T_K = repmat(T_K_scalar, nP, 1);
row.T_deg = repmat(T_deg_scalar, nP, 1);

row.is_positive_FeMg_garnet = repmat(is_positive_FeMg_garnet, nP, 1);
row.is_positive_FeMg_amphibole = repmat(is_positive_FeMg_amphibole, nP, 1);
row.is_positive_KD = repmat(is_positive_KD, nP, 1);
row.is_positive_denominator = repmat(is_positive_denominator, nP, 1);
row.is_XCa_g_reasonable = repmat(is_XCa_g_reasonable, nP, 1);
row.is_XMn_g_low = repmat(is_XMn_g_low, nP, 1);
row.is_below_850C = repmat(is_below_850C_scalar, nP, 1);
row.is_within_primary_calibration_T = ...
    row.T_deg >= 599 & row.T_deg <= 920;
row.is_within_primary_calibration_P = ...
    P_kbar >= 5 & P_kbar <= 18;
row.is_within_primary_calibration_XCa = ...
    isfinite(row.XCa_g) & row.XCa_g >= 0.165 & row.XCa_g <= 0.418;
row.recommended_by_GrahamPowell1984 = repmat(recommendedScalar, nP, 1);

end

function gt = prepareGarnetRow(data_garnet)
% prepareGarnetRow
% Read one garnet row while preserving any explicitly supplied NaN values.

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

function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Read one amphibole row while preserving any explicitly supplied NaN values.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();
[amp.Fe_raw, amp.Fe3, amp.Fe_used] = ...
    getFeUsed(data_amphibole, 'Amphibole');
amp.Mg = getRequiredVar(data_amphibole, 'Mg_cation_apfu', 'Amphibole');

% Output-only optional variables are NaN when absent.
amp.Mn = getOptionalVar(data_amphibole, 'Mn_cation_apfu', NaN, 'Amphibole');
amp.Si = getOptionalVar(data_amphibole, 'Si_cation_apfu', NaN, 'Amphibole');
amp.Al = getOptionalVar(data_amphibole, 'Al_cation_apfu', NaN, 'Amphibole');
amp.Ti = getOptionalVar(data_amphibole, 'Ti_cation_apfu', NaN, 'Amphibole');
amp.Ca = getOptionalVar(data_amphibole, 'Ca_cation_apfu', NaN, 'Amphibole');
amp.Na = getOptionalVar(data_amphibole, 'Na_cation_apfu', NaN, 'Amphibole');
amp.K = getOptionalVar(data_amphibole, 'K_cation_apfu', NaN, 'Amphibole');

end

function param = calcGrahamPowell1984Params(gt, amp)
% calcGrahamPowell1984Params
% Calculate the exchange coefficient and garnet Ca/Mn fractions.

param = struct();

sumDivalent_g = gt.Fe_used + gt.Mg + gt.Mn + gt.Ca;
if isfinite(sumDivalent_g) && sumDivalent_g > 0
    param.XCa_g = gt.Ca ./ sumDivalent_g;
    param.XMn_g = gt.Mn ./ sumDivalent_g;
else
    param.XCa_g = NaN;
    param.XMn_g = NaN;
end

if isfinite(gt.Fe_used) && isfinite(gt.Mg) && gt.Mg > 0
    param.FeMg_garnet = gt.Fe_used ./ gt.Mg;
else
    param.FeMg_garnet = NaN;
end

if isfinite(amp.Fe_used) && isfinite(amp.Mg) && amp.Mg > 0
    param.FeMg_amphibole = amp.Fe_used ./ amp.Mg;
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

function [Fe_raw, Fe3, Fe_used] = getFeUsed(tbl, mineralLabel)
% getFeUsed
% Read Fe without replacing explicitly supplied NaN values.

Fe_raw = getRequiredVar(tbl, 'Fe_cation_apfu', mineralLabel);
Fe3 = getOptionalVar(tbl, 'Fe3_cation_apfu', 0, mineralLabel);
Fe_used = Fe_raw + Fe3;

end

function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read one required scalar numeric value. NaN is allowed and preserved.

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

function value = getOptionalVar(tbl, varName, defaultValue, mineralLabel)
% getOptionalVar
% Read one optional scalar numeric value. If the column is absent, use the
% supplied default. If the column exists and contains NaN, preserve NaN.

if ismember(varName, tbl.Properties.VariableNames)
    value = tbl.(varName);
else
    value = defaultValue;
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
