function results = Ravna2000(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Amphibole/Ravna2000.m
% Tested with MATLAB R2024b
%
% Garnet-hornblende Fe2+-Mg exchange geothermometer
% Ravna, E.K. (2000)
% Lithos, 53, 265-277
% DOI: https://doi.org/10.1016/S0024-4937(00)00029-3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Amphibole
% analysis (selected by the user from tables) and calculates temperature
% using the Ravna (2000) garnet-hornblende Fe2+-Mg exchange geothermometer.
%
% The term hornblende in Ravna (2000) comprises the calcic amphiboles
% pargasite, hastingsite, tschermakite, hornblende, and edenite (p. 271).
% An arbitrary amphibole analysis must not automatically be assumed to lie
% within this calibration domain.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Grt-Amp pair. One output row is produced
% for every supplied pressure value, allowing the same function to operate
% with both startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Ravna (2000) calibrated the thermometer by multiple regression of 65
% accepted garnet-hornblende pairs: 22 experimental pairs from basaltic to
% intermediate compositions and 43 natural Grt-Cpx-Hbl assemblages from
% intermediate to basaltic rocks (abstract, p. 265; p. 274).
%
% The combined calibration-data envelope reported in the abstract is:
%
%   Temperature       : 515-1025 degreeC
%   Pressure          : 5-16 kbar
%
% The experimental subset covers 700-1025 degreeC and 10-15 kbar, whereas
% the natural subset covers 515-810 degreeC and 5-14 kbar (p. 270). The
% combined ranges are data envelopes, not separately demonstrated hard
% validity limits. Applications outside them are retained with non-stopping
% fprintf warnings.
%
% The garnet-composition envelope read from the accepted calibration data in
% Tables 1 and 2 is approximately:
%
%   Garnet XCa              : 0.154-0.473
%   Garnet XMn              : 0.004-0.158
%   Garnet XCa + XMn        : 0.167-0.516
%
% Ravna (2000) explicitly states that the geothermometer should not be
% applied to mineral compositions beyond the ranges used in the calibration
% (Tables 1-2; conclusion, p. 275). These numerical table envelopes are
% therefore used for non-stopping compositional warnings.
%
% Adding a pressure correction did not improve the regression statistics, so
% pressure was not included in the final equation (p. 274). P_kbar is stored
% for traceability but does not affect calculated temperature. This absence
% of a pressure term should not be interpreted as validation outside the
% 5-16 kbar calibration-data envelope.
%
% Ravna (2000) used Fe2+/Mg in the hornblende M1-M3 sites, not bulk
% amphibole Fe/Mg. Hornblende structural formulae and ferric/ferrous iron
% were recalculated using the empirical procedure of Schumacher (1991), with
% cations allocated among T, M1-M3, M4, and A sites (p. 271). The present
% function retains the simplified site-allocation algorithm of the supplied
% implementation; it follows the same allocation sequence but is not a full
% independent implementation of the Schumacher (1991) recalculation.
%
% The calibration excluded amphiboles with Ca_M4 below 1.50, analyses giving
% Ca_M4 above 2.00, analyses giving negative calculated Fe3+, ultrabasic
% experimental systems, replacement textures, and disequilibrium or suspect
% analyses (pp. 271-273). Pairs deviating irregularly by more than about
% +/-90 degreeC from the main regression trend were also excluded (p. 273).
% Compare texturally corresponding, equilibrated garnet and hornblende
% domains and avoid applying the thermometer to ultrabasic systems.
%
% This implementation issues non-stopping warnings when:
%   1) input pressure is outside 5-16 kbar,
%   2) a finite calculated temperature is outside 515-1025 degreeC,
%   3) garnet XCa, XMn, or XCa+XMn is outside the Table 1-2 envelope,
%   4) estimated Ca_M4 is outside 1.50-2.00 apfu,
%   5) the simplified amphibole site allocation is invalid, or
%   6) a thermometer input or calculated temperature is non-finite.
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
%   Fe_cation_apfu          % Fe2+ used in the exchange equation
%   Mg_cation_apfu
%   Ca_cation_apfu
%
% Required Amphibole variables:
%   Fe_cation_apfu          % Fe2+ before simplified site allocation
%   Mg_cation_apfu
%   Si_cation_apfu
%   Al_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%
% Optional Garnet variables:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Si_cation_apfu
%   Al_cation_apfu
%
% Optional Amphibole variables:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   K_cation_apfu
%   Ti_cation_apfu
%   Cr_cation_apfu
%
% Missing optional variables retain the legacy default of zero. If an
% optional variable is present and explicitly contains NaN, that NaN is never
% converted to zero. NaN in a variable used by the thermometer or site
% allocation propagates to the calculated temperature and generates a
% non-stopping fprintf warning.
%
% All supplied finite cation values must be greater than or equal to zero.
% Negative values and Inf values stop the calculation. Zero and NaN are
% accepted at input; mathematically undefined ratios or site allocations
% produce NaN and a non-stopping warning.
%
% IMPORTANT Fe and site convention:
% Fe_cation_apfu must represent Fe2+ and Fe3_cation_apfu must represent Fe3+.
% If Fe_cation_apfu contains total Fe, it must be separated into Fe2+ and Fe3+
% before this function is called. The original calibration used recalculated
% structural formulae rather than assuming all amphibole Fe was Fe2+ (p. 271).
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
%   T_K = (1504 + 1784*(XCa_g + XMn_g)) / (ln(KD) + 0.720)
%   T_degreeC = T_K - 273.15
%
% where
%
%   KD    = (Fe2+/Mg)_garnet / (Fe2+/Mg)_hornblende_M1-M3
%   XCa_g = Ca_g / (Ca_g + Mn_g + Fe2_g + Mg_g)
%   XMn_g = Mn_g / (Ca_g + Mn_g + Fe2_g + Mg_g)
%
% The published expression prints subtraction of 273 degreeC (p. 274). This
% implementation uses 273.15 for internally consistent Kelvin-to-degreeC
% conversion, preserving the convention of the supplied implementation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Ravna2000(rawdata_struct, P_kbar)
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
if nargin < 2
    error('Ravna2000 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || ~isreal(P_kbar) || isempty(P_kbar) || ...
        ~isvector(P_kbar) || any(~isfinite(P_kbar(:))) || ...
        any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
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
    {'Fe_cation_apfu', 'Mg_cation_apfu', 'Si_cation_apfu', ...
     'Al_cation_apfu', 'Ca_cation_apfu', 'Na_cation_apfu'}, ...
    'Amphibole');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Result table blocks are buffered and concatenated once after all selections,
% avoiding repeated full-table reallocation on every loop iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degC = 515;
calibrationT_max_degC = 1025;
calibrationP_min_kbar = 5;
calibrationP_max_kbar = 16;
calibrationXCa_min = 0.154;
calibrationXCa_max = 0.473;
calibrationXMn_min = 0.004;
calibrationXMn_max = 0.158;
calibrationXCaMn_min = 0.167;
calibrationXCaMn_max = 0.516;
calibrationCaM4_min = 1.50;
calibrationCaM4_max = 2.00;

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
        'PromptString', ...
        'Please select the Amphibole data you would like to use:', ...
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

    % Negative and infinite cation values are invalid. Zero and NaN are
    % accepted; undefined ratios or site allocations return NaN.
    validateNonnegativeInputs(selectedData_gt, selectedData_amp);

    % Record explicitly supplied NaN values used by the thermometer or the
    % simplified site allocation. Missing optional columns default to zero.
    nanInputNames = findNaNInputs(selectedData_gt, selectedData_amp);

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
            ': ' num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the calibration-data range ' ...
             'of Ravna (2000): 5-16 kbar (abstract, p. 265). %d of %d ' ...
             'pressure point(s) are outside the range; input range = ' ...
             '%.4g-%.4g kbar. Pressure is stored but is not used in the ' ...
             'published temperature equation (p. 274).\n'], ...
            sum(pressureOutsideCalibration), numel(P_kbar), ...
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
            ['WARNING: Calculated temperature is outside the combined ' ...
             'calibration-data range of Ravna (2000): 515-1025 degreeC ' ...
             '(abstract, p. 265; data subsets on p. 270). %d of %d finite ' ...
             'temperature point(s) are outside the range; calculated finite ' ...
             'range = %.4g-%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)));
    end

    xCaValue = row.XCa_g(1);
    xMnValue = row.XMn_g(1);
    xCaMnValue = row.XCa_plus_XMn_g(1);
    compositionOutsideCalibration = ...
        (isfinite(xCaValue) && ...
         (xCaValue < calibrationXCa_min || xCaValue > calibrationXCa_max)) || ...
        (isfinite(xMnValue) && ...
         (xMnValue < calibrationXMn_min || xMnValue > calibrationXMn_max)) || ...
        (isfinite(xCaMnValue) && ...
         (xCaMnValue < calibrationXCaMn_min || ...
          xCaMnValue > calibrationXCaMn_max));

    if compositionOutsideCalibration
        fprintf(2, ...
            ['WARNING: Garnet composition is outside the Ravna (2000) ' ...
             'calibration-data envelope (Tables 1-2, pp. 266, 268-269): ' ...
             'XCa = 0.154-0.473, XMn = 0.004-0.158, and XCa+XMn = ' ...
             '0.167-0.516. Calculated values for %s & %s are XCa_g = ' ...
             '%.4g, XMn_g = %.4g, and XCa_g+XMn_g = %.4g. Ravna (2000) ' ...
             'states that compositions beyond the calibration ranges ' ...
             'should not be used (p. 275).\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)), ...
            xCaValue, xMnValue, xCaMnValue);
    end

    caM4Value = row.Ca_B_amp(1);
    rawCaValue = row.Ca_amp(1);
    caM4OutsideCalibration = ...
        (isfinite(caM4Value) && caM4Value < calibrationCaM4_min) || ...
        (isfinite(rawCaValue) && rawCaValue > calibrationCaM4_max);

    if caM4OutsideCalibration
        fprintf(2, ...
            ['WARNING: Estimated amphibole Ca_M4 is outside the Ravna (2000) ' ...
             'calibration criterion of 1.50-2.00 apfu (pp. 271-273). ' ...
             'Simplified allocated Ca_M4 = %.4g and input amphibole Ca = ' ...
             '%.4g apfu for %s & %s.\n'], ...
            caM4Value, rawCaValue, char(string(selectedCode_gt)), ...
            char(string(selectedCode_amp)));
    end

    if ~all(row.is_C_site_capacity_ok & row.is_B_site_capacity_ok)
        fprintf(2, ...
            ['WARNING: The simplified amphibole site allocation is invalid ' ...
             'for %s & %s. Ravna (2000) used recalculated hornblende ' ...
             'structural formulae and M1-M3 Fe2+/Mg (p. 271). The affected ' ...
             'temperature values remain NaN.\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_amp)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer/site-allocation ' ...
             'input(s) for %s & %s: %s.\n' ...
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
        'Ravna2000', ...
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
% Confirm that the table contains every required thermometer input variable.

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
% Return explicitly supplied NaN names used by the thermometer or simplified
% site allocation. The output buffer is preallocated.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
amphiboleVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Si_cation_apfu', 'Al_cation_apfu', 'Ca_cation_apfu', ...
    'Na_cation_apfu', 'Fe3_cation_apfu', 'Mn_cation_apfu', ...
    'Ti_cation_apfu', 'Cr_cation_apfu'};

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
% Stop when a supplied cation value is negative, infinite, complex, or
% nonnumeric. Zero and NaN are accepted.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Fe3_cation_apfu', 'Mn_cation_apfu', ...
    'Si_cation_apfu', 'Al_cation_apfu'};
amphiboleVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Si_cation_apfu', 'Al_cation_apfu', 'Ca_cation_apfu', ...
    'Na_cation_apfu', 'Fe3_cation_apfu', 'Mn_cation_apfu', ...
    'K_cation_apfu', 'Ti_cation_apfu', 'Cr_cation_apfu'};

maximumNames = numel(garnetVariables) + numel(amphiboleVariables);
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

for i = 1:numel(amphiboleVariables)
    variableName = amphiboleVariables{i};
    if ismember(variableName, data_amphibole.Properties.VariableNames)
        variableValue = data_amphibole.(variableName);
        if ~isnumeric(variableValue) || ~isreal(variableValue) || ...
                any(isinf(variableValue(:))) || ...
                any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nInvalid = nInvalid + 1;
            invalidInputNames(nInvalid) = ...
                "Amphibole." + string(variableName);
        end
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNames(1:nInvalid);
    error(['Ravna2000: cation inputs must be finite or NaN and must ' ...
           'be >= 0. Negative, infinite, complex, or nonnumeric value(s) ' ...
           'were found in: ' char(strjoin(invalidInputNames, ', ')) '.']);
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
param = calcRavna2000Params(gt, amp);

if isfinite(param.KD) && param.KD > 0
    lnKD = log(param.KD);
else
    lnKD = NaN;
end

denominator = lnKD + 0.720;
canCalculate = isfinite(denominator) && denominator > 0 && ...
    isfinite(param.XCa_g) && isfinite(param.XMn_g) && ...
    amp.is_C_site_valid && amp.is_B_site_valid;

if canCalculate
    T_K_scalar = ...
        (1504 + 1784 .* (param.XCa_g + param.XMn_g)) ./ denominator;
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

if isfinite(amp.Fe2_C) && isfinite(amp.Mg_C) && ...
        amp.Fe2_C + amp.Mg_C > 0
    Mg_number_amp_C = amp.Mg_C ./ (amp.Fe2_C + amp.Mg_C);
else
    Mg_number_amp_C = NaN;
end

if isfinite(amp.Ca_B) && isfinite(amp.Na_B) && ...
        amp.Ca_B + amp.Na_B > 0
    Ca_fraction_amp_B = amp.Ca_B ./ (amp.Ca_B + amp.Na_B);
else
    Ca_fraction_amp_B = NaN;
end

is_positive_FeMg_garnet = isfinite(gt.Fe2) && isfinite(gt.Mg) && ...
    gt.Fe2 > 0 && gt.Mg > 0;
is_positive_FeMg_amphibole_C = ...
    isfinite(amp.Fe2_C) && isfinite(amp.Mg_C) && ...
    amp.Fe2_C > 0 && amp.Mg_C > 0;
is_positive_KD = isfinite(param.KD) && param.KD > 0;
is_positive_denominator = isfinite(denominator) && denominator > 0;
is_XCa_g_reasonable = isfinite(param.XCa_g) && ...
    param.XCa_g >= 0 && param.XCa_g <= 1;
is_XMn_g_reasonable = isfinite(param.XMn_g) && ...
    param.XMn_g >= 0 && param.XMn_g <= 1;

xCaMn = param.XCa_g + param.XMn_g;
is_XCa_within_calibration = isfinite(param.XCa_g) && ...
    param.XCa_g >= 0.154 && param.XCa_g <= 0.473;
is_XMn_within_calibration = isfinite(param.XMn_g) && ...
    param.XMn_g >= 0.004 && param.XMn_g <= 0.158;
is_XCaMn_within_calibration = isfinite(xCaMn) && ...
    xCaMn >= 0.167 && xCaMn <= 0.516;
is_CaM4_within_calibration = isfinite(amp.Ca_B) && ...
    isfinite(amp.Ca) && amp.Ca_B >= 1.50 && amp.Ca <= 2.00;
is_common_calcic_amphibole_like = is_CaM4_within_calibration;
is_C_site_capacity_ok = amp.is_C_site_valid;
is_B_site_capacity_ok = amp.is_B_site_valid;

is_pressure_within_calibration = P_kbar >= 5 & P_kbar <= 16;
is_temperature_within_calibration = repmat( ...
    isfinite(T_deg_scalar) && T_deg_scalar >= 515 && T_deg_scalar <= 1025, ...
    nP, 1);

recommendedScalar = is_positive_FeMg_garnet && ...
    is_positive_FeMg_amphibole_C && is_positive_KD && ...
    is_positive_denominator && is_XCa_within_calibration && ...
    is_XMn_within_calibration && is_XCaMn_within_calibration && ...
    is_CaM4_within_calibration && is_C_site_capacity_ok && ...
    is_B_site_capacity_ok;
recommended = repmat(recommendedScalar, nP, 1) & ...
    is_pressure_within_calibration & is_temperature_within_calibration;

row = table();
row.P_kbar = P_kbar;

% Repeat scalar compositions, allocations, and thermometer terms for every
% pressure row. Temperature is repeated because the equation has no P term.
row.Fe2_g = repmat(gt.Fe2, nP, 1);
row.Fe3_g = repmat(gt.Fe3, nP, 1);
row.Fe_total_g = repmat(gt.Fe_total, nP, 1);
row.Mg_g = repmat(gt.Mg, nP, 1);
row.Mn_g = repmat(gt.Mn, nP, 1);
row.Ca_g = repmat(gt.Ca, nP, 1);
row.Si_g = repmat(gt.Si, nP, 1);
row.Al_g = repmat(gt.Al, nP, 1);

row.Fe2_amp = repmat(amp.Fe2, nP, 1);
row.Fe3_amp = repmat(amp.Fe3, nP, 1);
row.Fe_total_amp = repmat(amp.Fe_total, nP, 1);
row.Mg_amp = repmat(amp.Mg, nP, 1);
row.Mn_amp = repmat(amp.Mn, nP, 1);
row.Ca_amp = repmat(amp.Ca, nP, 1);
row.Na_amp = repmat(amp.Na, nP, 1);
row.K_amp = repmat(amp.K, nP, 1);
row.Ti_amp = repmat(amp.Ti, nP, 1);
row.Cr_amp = repmat(amp.Cr, nP, 1);
row.Si_amp = repmat(amp.Si, nP, 1);
row.Al_amp = repmat(amp.Al, nP, 1);

row.Al_T_amp = repmat(amp.Al_T, nP, 1);
row.Al_C_amp = repmat(amp.Al_C, nP, 1);
row.Fe2_C_amp = repmat(amp.Fe2_C, nP, 1);
row.Mg_C_amp = repmat(amp.Mg_C, nP, 1);
row.Mn_C_amp = repmat(amp.Mn_C, nP, 1);
row.Ca_B_amp = repmat(amp.Ca_B, nP, 1);
row.Na_B_amp = repmat(amp.Na_B, nP, 1);
row.Fe2_B_amp = repmat(amp.Fe2_B, nP, 1);
row.Mn_B_amp = repmat(amp.Mn_B, nP, 1);
row.Na_A_amp = repmat(amp.Na_A, nP, 1);
row.K_A_amp = repmat(amp.K_A, nP, 1);

row.C_site_sum_amp = repmat(amp.C_site_sum, nP, 1);
row.B_site_sum_amp = repmat(amp.B_site_sum, nP, 1);
row.A_site_sum_amp = repmat(amp.A_site_sum, nP, 1);

row.FeMg_garnet = repmat(param.FeMg_garnet, nP, 1);
row.FeMg_amphibole_M1M3 = ...
    repmat(param.FeMg_amphibole_M1M3, nP, 1);
row.KD = repmat(param.KD, nP, 1);
row.lnKD = repmat(lnKD, nP, 1);
row.XCa_g = repmat(param.XCa_g, nP, 1);
row.XMn_g = repmat(param.XMn_g, nP, 1);
row.XCa_plus_XMn_g = repmat(xCaMn, nP, 1);

row.Mg_number_garnet = repmat(Mg_number_gt, nP, 1);
row.Mg_number_amphibole_C = repmat(Mg_number_amp_C, nP, 1);
row.Ca_fraction_amphibole_B_over_CaNa = ...
    repmat(Ca_fraction_amp_B, nP, 1);

row.T_K = repmat(T_K_scalar, nP, 1);
row.T_deg = repmat(T_deg_scalar, nP, 1);

row.is_positive_FeMg_garnet = ...
    repmat(is_positive_FeMg_garnet, nP, 1);
row.is_positive_FeMg_amphibole_M1M3 = ...
    repmat(is_positive_FeMg_amphibole_C, nP, 1);
row.is_positive_KD = repmat(is_positive_KD, nP, 1);
row.is_positive_denominator = repmat(is_positive_denominator, nP, 1);
row.is_XCa_g_reasonable = repmat(is_XCa_g_reasonable, nP, 1);
row.is_XMn_g_reasonable = repmat(is_XMn_g_reasonable, nP, 1);
row.is_common_calcic_amphibole_like = ...
    repmat(is_common_calcic_amphibole_like, nP, 1);
row.is_C_site_capacity_ok = repmat(is_C_site_capacity_ok, nP, 1);
row.is_B_site_capacity_ok = repmat(is_B_site_capacity_ok, nP, 1);
row.is_XCa_within_calibration = ...
    repmat(is_XCa_within_calibration, nP, 1);
row.is_XMn_within_calibration = ...
    repmat(is_XMn_within_calibration, nP, 1);
row.is_XCaMn_within_calibration = ...
    repmat(is_XCaMn_within_calibration, nP, 1);
row.is_CaM4_within_calibration = ...
    repmat(is_CaM4_within_calibration, nP, 1);
row.is_pressure_within_calibration = is_pressure_within_calibration;
row.is_temperature_within_calibration = ...
    is_temperature_within_calibration;
row.recommended_by_Ravna2000 = recommended;

end


function gt = prepareGarnetRow(data_garnet)
% prepareGarnetRow
% Read one garnet row while preserving explicitly supplied NaN values.

if height(data_garnet) ~= 1
    error('Garnet input must be a 1-row table.');
end

gt = struct();
gt.Fe2 = getRequiredVar(data_garnet, 'Fe_cation_apfu', 'Garnet');
gt.Fe3 = getOptionalVar(data_garnet, 'Fe3_cation_apfu', 0, 'Garnet');
gt.Fe_total = gt.Fe2 + gt.Fe3;
gt.Mg = getRequiredVar(data_garnet, 'Mg_cation_apfu', 'Garnet');
gt.Ca = getRequiredVar(data_garnet, 'Ca_cation_apfu', 'Garnet');
gt.Mn = getOptionalVar(data_garnet, 'Mn_cation_apfu', 0, 'Garnet');
gt.Si = getOptionalVar(data_garnet, 'Si_cation_apfu', 0, 'Garnet');
gt.Al = getOptionalVar(data_garnet, 'Al_cation_apfu', 0, 'Garnet');

end


function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Read one amphibole row and perform the simplified site allocation while
% preserving explicitly supplied NaN values.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();
amp.Fe2 = getRequiredVar( ...
    data_amphibole, 'Fe_cation_apfu', 'Amphibole');
amp.Fe3 = getOptionalVar( ...
    data_amphibole, 'Fe3_cation_apfu', 0, 'Amphibole');
amp.Fe_total = amp.Fe2 + amp.Fe3;
amp.Mg = getRequiredVar( ...
    data_amphibole, 'Mg_cation_apfu', 'Amphibole');
amp.Si = getRequiredVar( ...
    data_amphibole, 'Si_cation_apfu', 'Amphibole');
amp.Al = getRequiredVar( ...
    data_amphibole, 'Al_cation_apfu', 'Amphibole');
amp.Ca = getRequiredVar( ...
    data_amphibole, 'Ca_cation_apfu', 'Amphibole');
amp.Na = getRequiredVar( ...
    data_amphibole, 'Na_cation_apfu', 'Amphibole');
amp.Mn = getOptionalVar( ...
    data_amphibole, 'Mn_cation_apfu', 0, 'Amphibole');
amp.K = getOptionalVar( ...
    data_amphibole, 'K_cation_apfu', 0, 'Amphibole');
amp.Ti = getOptionalVar( ...
    data_amphibole, 'Ti_cation_apfu', 0, 'Amphibole');
amp.Cr = getOptionalVar( ...
    data_amphibole, 'Cr_cation_apfu', 0, 'Amphibole');

amp = allocateAmphiboleSites(amp);

end


function param = calcRavna2000Params(gt, amp)
% calcRavna2000Params
% Calculate garnet fractions and the M1-M3 Fe2+-Mg distribution coefficient.

param = struct();

sumDivalent_g = gt.Ca + gt.Mn + gt.Fe2 + gt.Mg;
if isfinite(sumDivalent_g) && sumDivalent_g > 0
    param.XCa_g = gt.Ca ./ sumDivalent_g;
    param.XMn_g = gt.Mn ./ sumDivalent_g;
else
    param.XCa_g = NaN;
    param.XMn_g = NaN;
end

if isfinite(gt.Fe2) && isfinite(gt.Mg) && gt.Mg > 0
    param.FeMg_garnet = gt.Fe2 ./ gt.Mg;
else
    param.FeMg_garnet = NaN;
end

if isfinite(amp.Fe2_C) && isfinite(amp.Mg_C) && amp.Mg_C > 0
    param.FeMg_amphibole_M1M3 = amp.Fe2_C ./ amp.Mg_C;
else
    param.FeMg_amphibole_M1M3 = NaN;
end

if isfinite(param.FeMg_garnet) && ...
        isfinite(param.FeMg_amphibole_M1M3) && ...
        param.FeMg_amphibole_M1M3 > 0
    param.KD = param.FeMg_garnet ./ param.FeMg_amphibole_M1M3;
else
    param.KD = NaN;
end

end


function amp = allocateAmphiboleSites(amp)
% allocateAmphiboleSites
% Simplified allocation retained from the supplied Ravna2000 implementation.
% NaN in a cation used here invalidates the allocation rather than being
% silently omitted by min/max operations.

T_capacity = 8.0;
C_capacity = 5.0;
B_capacity = 2.0;

allocationInputs = [amp.Fe2, amp.Fe3, amp.Mg, amp.Si, amp.Al, ...
    amp.Ca, amp.Na, amp.Mn, amp.Ti, amp.Cr];
if any(isnan(allocationInputs))
    amp.Al_T = NaN;
    amp.Al_C = NaN;
    amp.Fe2_C = NaN;
    amp.Mg_C = NaN;
    amp.Mn_C = NaN;
    amp.Ca_B = NaN;
    amp.Na_B = NaN;
    amp.Fe2_B = NaN;
    amp.Mn_B = NaN;
    amp.Na_A = NaN;
    amp.K_A = amp.K;
    amp.C_site_sum = NaN;
    amp.B_site_sum = NaN;
    amp.A_site_sum = NaN;
    amp.is_C_site_valid = false;
    amp.is_B_site_valid = false;
    return
end

T_deficit = max(0, T_capacity - amp.Si);
amp.Al_T = min(amp.Al, T_deficit);
amp.Al_C = max(0, amp.Al - amp.Al_T);

C_fixed = amp.Al_C + amp.Ti + amp.Cr + amp.Fe3 + amp.Mg;
C_remaining = C_capacity - C_fixed;

if C_remaining < -1e-8
    amp.Fe2_C = NaN;
    amp.Mg_C = amp.Mg;
    amp.Mn_C = NaN;
    amp.Fe2_B = NaN;
    amp.Mn_B = NaN;
    amp.Ca_B = min(amp.Ca, B_capacity);
    amp.Na_B = NaN;
    amp.Na_A = NaN;
    amp.K_A = amp.K;
    amp.C_site_sum = C_fixed;
    amp.B_site_sum = NaN;
    amp.A_site_sum = NaN;
    amp.is_C_site_valid = false;
    amp.is_B_site_valid = false;
    return
end

C_remaining = max(0, C_remaining);
amp.Fe2_C = min(amp.Fe2, C_remaining);
remainingAfterFe2 = C_remaining - amp.Fe2_C;
amp.Mn_C = min(amp.Mn, remainingAfterFe2);

remainingFe2 = amp.Fe2 - amp.Fe2_C;
remainingMn = amp.Mn - amp.Mn_C;

Ca_for_B = min(amp.Ca, B_capacity);
B_remaining = B_capacity - Ca_for_B;
Na_for_B = min(amp.Na, B_remaining);
B_remaining = B_remaining - Na_for_B;
Fe2_for_B = min(remainingFe2, B_remaining);
B_remaining = B_remaining - Fe2_for_B;
Mn_for_B = min(remainingMn, B_remaining);

amp.Ca_B = Ca_for_B;
amp.Na_B = Na_for_B;
amp.Fe2_B = Fe2_for_B;
amp.Mn_B = Mn_for_B;
amp.Na_A = max(0, amp.Na - amp.Na_B);
amp.K_A = amp.K;

amp.Mg_C = amp.Mg;
amp.C_site_sum = amp.Al_C + amp.Ti + amp.Cr + amp.Fe3 + ...
    amp.Mg_C + amp.Fe2_C + amp.Mn_C;
amp.B_site_sum = amp.Ca_B + amp.Na_B + amp.Fe2_B + amp.Mn_B;
amp.A_site_sum = amp.Na_A + amp.K_A;

amp.is_C_site_valid = ...
    amp.C_site_sum <= C_capacity + 1e-6;
amp.is_B_site_valid = ...
    amp.B_site_sum <= B_capacity + 1e-6;

end


function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read one required numeric scalar. NaN is allowed and preserved.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('Variable %s in %s table must be a real numeric scalar.', ...
        varName, mineralLabel);
end
if isinf(value) || (isfinite(value) && value < 0)
    error('Variable %s in %s table must be NaN or a finite value >= 0.', ...
        varName, mineralLabel);
end

end


function value = getOptionalVar(tbl, varName, defaultValue, mineralLabel)
% getOptionalVar
% Use the supplied default only if an optional column is absent. Explicit NaN
% is preserved and is never converted to zero.

if ismember(varName, tbl.Properties.VariableNames)
    value = tbl.(varName);
else
    value = defaultValue;
end

if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('Variable %s in %s table must be a real numeric scalar.', ...
        varName, mineralLabel);
end
if isinf(value) || (isfinite(value) && value < 0)
    error('Variable %s in %s table must be NaN or a finite value >= 0.', ...
        varName, mineralLabel);
end

end
