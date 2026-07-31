function results = Wu2018(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/Wu2018.m
% Tested with MATLAB R2024b
%
% Empirical Fe-Mg exchange thermometer between Garnet and Muscovite
% Wu, C.-M. (2018)
% Journal of Earth Science, 29, 977-988
% DOI: https://doi.org/10.1007/s12583-018-0851-z
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% analysis (selected independently by the user from tables) and calculates
% temperature using the Wu (2018) garnet-muscovite Fe-Mg exchange
% thermometer.
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Garnet-Mica pair, the output contains one row per pressure
% value. This allows the same function to be called by both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Garnet-Mica pair and appends the new
% table block to the final output.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Wu (2018) calibrated the thermometer using 161 natural Al2SiO5-bearing
% metapelitic samples over the following ranges:
%
%   Temperature               : 460-760 degreeC
%   Pressure                  : 1-12 kbar
%   Garnet composition        : XFe = 0.53-0.91
%                               XMg = 0.04-0.27
%                               XCa = 0.02-0.22
%                               XMn = 0.00-0.29
%   Muscovite composition     : Fe   = 0.03-0.21 apfu
%                               Mg   = 0.02-0.32 apfu
%                               AlVI = 1.62-1.96 apfu
%                               (11-oxygen basis)
%   Intended application      : equilibrated crustal metapelites, including
%                               CaO-saturated and CaO-undersaturated rocks
%
% The P-T and mineral-composition ranges are summarized in the abstract on
% p. 977 and described with the calibration data on p. 978. The garnet and
% muscovite activity-model assumptions are described on p. 979, and the
% thermometer formulation (Eq. 14) is presented on p. 980.
%
% Important application cautions stated or demonstrated in the article:
% - This is an empirical calibration based on natural metapelites, not a
%   direct experimental calibration. Calibration P-T values were obtained
%   using the Holdaway (2000) Garnet-Biotite thermometer and Holdaway (2001)
%   GASP barometer, so their systematic uncertainties can be inherited by
%   the Wu (2018) thermometer (p. 978).
% - Suitable calibration samples were required to show equilibrium textures,
%   lack retrogressive metamorphism, and high-quality mineral analyses. For
%   growth-zoned garnet, the garnet rim and corresponding matrix-mineral rim
%   compositions were used. Disequilibrium assemblages were excluded
%   (p. 978).
% - The garnet activity model treats Fe3+ and Cr3+ as negligible and uses the
%   divalent Fe-Mg-Ca-Mn components. Fe_cation_apfu is therefore treated as
%   Fe2+ in this implementation; Fe3_cation_apfu is retained for reference
%   but is not added to the thermometer Fe term (p. 979).
% - The paper reports that the garnet activity model is invalid for an
%   MnO-rich garnet beyond the calibrated XMn limit and that an unexpectedly
%   high temperature was obtained. Application above XMn = 0.29 is therefore
%   especially unreliable (p. 984).
% - Post-peak Fe-Mg re-equilibration between garnet and muscovite, or
%   formation of retrograde muscovite, can yield temperatures below the peak
%   metamorphic temperature (p. 985). Pair analyses representing the same
%   equilibrium stage.
% - If Al2SiO5 is absent, Wu (2018) applies the thermometer using an assumed
%   pressure. An input-pressure error of +/-1 kbar produces approximately
%   +/-10-17 degreeC temperature error (pp. 981 and 984).
% - The thermometer reproduced the reference Garnet-Biotite temperatures
%   within approximately +/-50 degreeC for 88% of calibrants, but the total
%   random error is estimated at approximately +/-60 degreeC. Its precision
%   cannot be established directly because phase-equilibrium experimental
%   data are unavailable (pp. 980-981; conclusions on p. 986).
%
% Wu (2018) does not provide a separate P-T application range wider than the
% empirical calibration. This implementation therefore uses 460-760 degreeC
% and 1-12 kbar as the non-stopping warning limits. It also reports finite
% mineral compositions outside the ranges represented by the calibrants.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Mica   : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following variable names (normalized cations):
%
%   Garnet table variables used in the calculation:
%     Fe_cation_apfu       % Fe2+ in garnet (assumed)
%     Mg_cation_apfu
%
%   Mica table variables used in the calculation:
%     Fe_cation_apfu       % Fe2+ in mica (assumed)
%     Mg_cation_apfu
%     Al_cation_apfu
%     Si_cation_apfu
%
% Optional variables used in the calculation when present:
%   Garnet : Mn_cation_apfu, Ca_cation_apfu
%
% Optional variables retained in the output table when present:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Ti_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Missing optional variables are represented by zero. An optional variable
% that exists and contains NaN remains NaN; it is never converted to zero.
% AlVI in mica is calculated on an 11-oxygen basis as Al + Si - 4.
%
% All finite mineral-composition values used by the thermometer must be
% non-negative. Negative and infinite values stop the calculation. Zero is
% allowed, although it may produce a non-finite exchange coefficient or
% temperature. NaN values are retained, propagated through the calculation,
% and reported by non-stopping fprintf messages.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Wu, 2018, Eq. 14, p. 980)
%
% This implementation uses the Wu (2018) garnet-muscovite thermometer:
%
%   T(K) = [33648.0 + [2.0 - 0.7*(Feb - Mgb)]*P(bar)
%           - 0.7*Fec + 0.7*Mgc
%           - 142260.0*(X_Fe_Ms - X_Mg_Ms)
%           + 68610.0*X_Al_Ms]
%          /
%          [100.0 + 0.7*(R*ln(K_D1) + Fea - Mga)]
%
% where
%   K_D1 = (X_alm^Grt * X_Mg-cel^Ms) / (X_prp^Grt * X_Fe-cel^Ms)
%
% and, using the definitions in Wu (2018) Table 1,
%   X_alm^Grt    = (X_Fe^Grt)^3
%   X_prp^Grt    = (X_Mg^Grt)^3
%   X_Mg-cel^Ms  = 4*X_Mg^Ms*X_Al^Ms
%   X_Fe-cel^Ms  = 4*X_Fe^Ms*X_Al^Ms
%
% so that
%   ln(K_D1) = 3*ln(X_Fe^Grt / X_Mg^Grt) + ln(X_Mg^Ms / X_Fe^Ms)
%
% Garnet activity terms Fea, Feb, Fec, Mga, Mgb, Mgc are from:
% Wu, Zhang and Ren (2004)
% Journal of Petrology, 45, 1907-1921, Appendix A4-A9.
%
% NOTES
% - In the Wu (2018) formulation, Fe3+ is not used explicitly in the
%   thermometer equation or in the garnet/muscovite mole-fraction terms.
% - Therefore, Fe_cation_apfu is treated as Fe2+ for the thermometer.
% - Optional Fe3_cation_apfu is retained in output for reference only.
% - P_kbar is converted internally to P_bar because Eq. (14) uses P(bar).
% - Zero is allowed as requested, but logarithms or divisions involving zero
%   can produce NaN or Inf; these results are retained and reported.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Wu2018(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Mica tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Mica pair. The output variable set is
%             intended to remain stable for downstream processing.

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('Wu2018 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        ~isreal(P_kbar) || any(~isfinite(P_kbar(:))) || ...
        any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative real numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The tables are not
% modified here; only the relevant columns are read during calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_grt = rawdata_struct.Garnet;
dataset_mica = rawdata_struct.Mica;

requiredVariables_grt = {'Fe_cation_apfu', 'Mg_cation_apfu'};
requiredVariables_mica = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Al_cation_apfu', 'Si_cation_apfu'};

validateRequiredVariables(dataset_grt, requiredVariables_grt, 'Garnet');
validateRequiredVariables(dataset_mica, requiredVariables_mica, 'Mica');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated table concatenation inside the loop is avoided because it forces
% MATLAB to repeatedly reallocate and copy the entire results table.
%
% The cell buffer is allocated in blocks and expanded only when necessary.
% After the interactive loop finishes, all blocks are concatenated once.
disp('=== Step 2: Preparing output container ===');

bufferBlockSize = 16;
resultBlocks = cell(bufferBlockSize, 1);
nResultBlocks = 0;

% Empirical calibration limits reported by Wu (2018).
calibrationT_min_degC = 460;
calibrationT_max_degC = 760;
calibrationP_min_kbar = 1;
calibrationP_max_kbar = 12;

% Pressure is common to all selected mineral pairs in this function call.
% The pressure warning is therefore printed only once, after the first
% completed calculation.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
    dataCodes_grt = dataset_grt{:, 1};

    [selectedIdx_grt, ok] = listdlg( ...
        'PromptString', 'Please select the Garnet data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_grt)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_grt)
        disp('Selection canceled');
        break;
    end

    selectedCode_grt = dataCodes_grt(selectedIdx_grt);
    disp(['Garnet selected: ' char(string(selectedCode_grt))]);

    % ----- Mica selection -----
    disp('=== Step 4: Selecting a data code from the list (Mica) ===');

    dataCodes_mica = dataset_mica{:, 1};

    [selectedIdx_mica, ok] = listdlg( ...
        'PromptString', 'Please select the Mica data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_mica)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_mica)
        disp('Selection canceled');
        break;
    end

    selectedCode_mica = dataCodes_mica(selectedIdx_mica);
    disp(['Mica selected: ' char(string(selectedCode_mica))]);

    % ----- Calculation -----
    % Garnet and mica are selected independently; do not assume row indices
    % correspond between the two datasets.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_mica = dataset_mica(selectedIdx_mica, :);

    % Check all raw variables actually used by the thermometer. The
    % calculation is intentionally allowed to continue when NaN is present.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_mica);

    % Negative, infinite, nonnumeric, nonreal, and nonscalar calculation
    % inputs stop the calculation. NaN and zero remain allowed.
    validateNonNegativeInputs(selectedData_grt, selectedData_mica);

    row = calcTemp(selectedData_grt, selectedData_mica, P_kbar);

    % Store the selected identifiers for every pressure point.
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_mica = repmat(string(selectedCode_mica), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_mica'}, 'Before', 1);

    % Store this result as one table block. Expand the buffer only in fixed
    % blocks, rather than growing the results table on every loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(bufferBlockSize, 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': Wu2018 = ' num2str(row.TWu2018_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': Wu2018 = ' num2str(row.TWu2018_C(1)) ' to ' ...
            num2str(row.TWu2018_C(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside the empirical
    % calibration range of 1-12 kbar. The calculation is not stopped.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the empirical calibration ' ...
             'range of Wu (2018): 1-12 kbar. %d of %d pressure point(s) ' ...
             'are outside the range; input range = %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the empirical
    % calibration range of 460-760 degreeC. NaN and Inf are handled below.
    finiteTemperature = isfinite(row.TWu2018_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.TWu2018_C < calibrationT_min_degC | ...
         row.TWu2018_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.TWu2018_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the empirical ' ...
             'calibration range of Wu (2018): 460-760 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g-%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)));
    end

    % Report mineral compositions outside the finite ranges represented by
    % the empirical calibration dataset. NaN inputs are handled separately.
    printCompositionRangeMessages(row, selectedCode_grt, selectedCode_mica);

    % Print a non-stopping warning when an input used by the thermometer
    % contains NaN. fprintf is used so the message remains visible even when
    % MATLAB warnings have been disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, and NaN values were not converted to zero.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report NaN/Inf results caused by missing data, zero
    % denominators, logarithms at zero, or other non-finite operations.
    invalidTemperature = ~isfinite(row.TWu2018_C);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite Wu2018 temperature values were calculated ' ...
             'for %s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            sum(invalidTemperature), ...
            numel(row.TWu2018_C), ...
            sum(isnan(row.TWu2018_C)), ...
            sum(isinf(row.TWu2018_C)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Wu2018', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate the buffered table blocks only once after all selections have
% been completed. Return an empty table when no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(dataTable, requiredVariables, mineralLabel)
% validateRequiredVariables
% Confirm that every variable required by the thermometer exists before the
% user enters the interactive selection loop.

missingVariables = requiredVariables(~ismember( ...
    requiredVariables, dataTable.Properties.VariableNames));

if ~isempty(missingVariables)
    error('%s table must contain variable(s): %s', ...
        mineralLabel, char(strjoin(string(missingVariables), ', ')));
end

end

function nanInputNames = findNaNInputs(data_garnet, data_mica)
% findNaNInputs
% Return the names of raw thermometer inputs that contain NaN. The output
% buffer is preallocated so its size does not grow inside the checking loops.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Al_cation_apfu', 'Si_cation_apfu'};

nanInputNames = strings(numel(garnetVariables) + numel(micaVariables), 1);
nNaNInputs = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if isnumeric(variableValue) && any(isnan(variableValue(:)))
            nNaNInputs = nNaNInputs + 1;
            nanInputNames(nNaNInputs) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
    if ismember(variableName, data_mica.Properties.VariableNames)
        variableValue = data_mica.(variableName);
        if isnumeric(variableValue) && any(isnan(variableValue(:)))
            nNaNInputs = nNaNInputs + 1;
            nanInputNames(nNaNInputs) = "Mica." + string(variableName);
        end
    end
end

nanInputNames = nanInputNames(1:nNaNInputs);

end

function validateNonNegativeInputs(data_garnet, data_mica)
% validateNonNegativeInputs
% Stop when a raw input used by the thermometer is negative, infinite,
% nonnumeric, nonreal, or nonscalar. NaN and zero are intentionally allowed.
% The error-name buffer is preallocated and does not grow inside the loops.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Al_cation_apfu', 'Si_cation_apfu'};

invalidInputNames = strings(numel(garnetVariables) + numel(micaVariables), 1);
nInvalidInputs = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if isInvalidCationValue(variableValue)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputNames(nInvalidInputs) = ...
                "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
    if ismember(variableName, data_mica.Properties.VariableNames)
        variableValue = data_mica.(variableName);
        if isInvalidCationValue(variableValue)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputNames(nInvalidInputs) = ...
                "Mica." + string(variableName);
        end
    end
end

invalidInputNames = invalidInputNames(1:nInvalidInputs);

if ~isempty(invalidInputNames)
    error(['Wu2018: thermometer inputs must be real, scalar, ' ...
           'non-negative, and not Inf. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function invalid = isInvalidCationValue(value)
% isInvalidCationValue
% NaN and zero are valid at this stage. Negative, infinite, nonnumeric,
% nonreal, and nonscalar values are invalid.

invalid = ~isnumeric(value) || ~isscalar(value) || ~isreal(value);

if ~invalid
    invalid = isinf(value) || (isfinite(value) && value < 0);
end

end

function row = calcTemp(data_grt, data_mica, P_kbar)
% calcTemp
% Compute Wu (2018) Eq. (14) for one garnet row, one mica row, and every
% pressure supplied in P_kbar.
%
% Inputs:
%   data_grt  : 1-row table containing garnet apfu cations
%   data_mica : 1-row table containing mica apfu cations
%   P_kbar    : pressure in kbar (column vector after input validation)
%
% Output:
%   row : table with one row per pressure value, including mineral inputs,
%         compositional terms, activity terms, and calculated temperature.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

row = table();

% --- Physical constant ---
R_J = 8.314462618;

% --- Prepare one-row mineral data ---
grt = prepareMineralRow(data_grt, 'Garnet');
mica = prepareMineralRow(data_mica, 'Mica');

% AlVI is a thermometer input derived from the required mica Al and Si
% analyses. NaN is retained; a finite negative value is rejected.
AlVI_mica_scalar = mica.Al + mica.Si - 4.0;

if isfinite(AlVI_mica_scalar) && AlVI_mica_scalar < 0
    error(['Wu2018: calculated Mica.AlVI is negative. Check the ' ...
           'Al and Si cation normalization (the calculation assumes an ' ...
           '11-oxygen mica basis).']);
end

% --- Replicate scalar mineral inputs for every pressure point ---
Fe2_grt = repmat(grt.Fe2, nP, 1);
Fe3_grt = repmat(grt.Fe3, nP, 1);
FeUsed_grt = repmat(grt.FeUsed, nP, 1);
Mg_grt = repmat(grt.Mg, nP, 1);
Mn_grt = repmat(grt.Mn, nP, 1);
Ca_grt = repmat(grt.Ca, nP, 1);
Al_grt = repmat(grt.Al, nP, 1);
Si_grt = repmat(grt.Si, nP, 1);

Fe2_mica = repmat(mica.Fe2, nP, 1);
Fe3_mica = repmat(mica.Fe3, nP, 1);
FeUsed_mica = repmat(mica.FeUsed, nP, 1);
Mg_mica = repmat(mica.Mg, nP, 1);
Mn_mica = repmat(mica.Mn, nP, 1);
Ca_mica = repmat(mica.Ca, nP, 1);
Ti_mica = repmat(mica.Ti, nP, 1);
Al_mica = repmat(mica.Al, nP, 1);
Si_mica = repmat(mica.Si, nP, 1);
K_mica = repmat(mica.K, nP, 1);
Na_mica = repmat(mica.Na, nP, 1);
AlVI_mica = repmat(AlVI_mica_scalar, nP, 1);

% --- Garnet site fractions (Wu, 2018, Table 1; Fe2+ only) ---
sum_grt_div = Fe2_grt + Mg_grt + Ca_grt + Mn_grt;

XFe_grt = Fe2_grt ./ sum_grt_div;
XMg_grt = Mg_grt ./ sum_grt_div;
XCa_grt = Ca_grt ./ sum_grt_div;
XMn_grt = Mn_grt ./ sum_grt_div;

% --- Mica octahedral fractions (Wu, 2018, Table 1; Fe2+ only) ---
sum_ms_div = Fe2_mica + Mg_mica + AlVI_mica;

XFe_ms = Fe2_mica ./ sum_ms_div;
XMg_ms = Mg_mica ./ sum_ms_div;
XAl_ms = AlVI_mica ./ sum_ms_div;

% --- ln(K_D1) from Wu (2018), Table 1 ---
lnK_D1 = 3 .* log(XFe_grt ./ XMg_grt) + log(XMg_ms ./ XFe_ms);

% --- Garnet activity polynomial terms from Wu et al. (2004) Appendix ---
act = calcWu2004bGarnetActivityTerms( ...
    XFe_grt, XMg_grt, XCa_grt, XMn_grt);

Fea = act.Fea;
Feb = act.Feb;
Fec = act.Fec;
Mga = act.Mga;
Mgb = act.Mgb;
Mgc = act.Mgc;

% --- Wu (2018), Eq. (14) ---
% The published equation uses pressure in bar. P_kbar is therefore
% converted to P_bar before evaluating the pressure term.
num = 33648.0 ...
    + (2.0 - 0.7 .* (Feb - Mgb)) .* P_bar ...
    - 0.7 .* Fec ...
    + 0.7 .* Mgc ...
    - 142260.0 .* (XFe_ms - XMg_ms) ...
    + 68610.0 .* XAl_ms;

den = 100.0 + 0.7 .* (R_J .* lnK_D1 + Fea - Mga);

% Direct calculation intentionally retains NaN and Inf results.
TWu2018_K = num ./ den;
TWu2018_C = TWu2018_K - 273.15;

% --- Pack pressure ---
row.P_kbar = P_kbar;
row.P_bar = P_bar;

% --- Pack garnet inputs ---
row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.FeUsed_grt = FeUsed_grt;
row.Mg_grt = Mg_grt;
row.Mn_grt = Mn_grt;
row.Ca_grt = Ca_grt;
row.Al_grt = Al_grt;
row.Si_grt = Si_grt;

% --- Pack mica inputs ---
row.Fe2_mica = Fe2_mica;
row.Fe3_mica = Fe3_mica;
row.FeUsed_mica = FeUsed_mica;
row.Mg_mica = Mg_mica;
row.Mn_mica = Mn_mica;
row.Ca_mica = Ca_mica;
row.Ti_mica = Ti_mica;
row.Al_mica = Al_mica;
row.Si_mica = Si_mica;
row.K_mica = K_mica;
row.Na_mica = Na_mica;

% --- Pack compositional and activity terms ---
row.sumDiv_grt = sum_grt_div;
row.XFe_grt = XFe_grt;
row.XMg_grt = XMg_grt;
row.XCa_grt = XCa_grt;
row.XMn_grt = XMn_grt;

row.AlVI_mica = AlVI_mica;
row.sumDiv_mica = sum_ms_div;
row.XFe_mica = XFe_ms;
row.XMg_mica = XMg_ms;
row.XAl_mica = XAl_ms;

row.lnK_D1 = lnK_D1;

row.Fea = Fea;
row.Feb = Feb;
row.Fec = Fec;
row.Mga = Mga;
row.Mgb = Mgb;
row.Mgc = Mgc;

row.num_Wu2018 = num;
row.den_Wu2018 = den;

row.TWu2018_K = TWu2018_K;
row.TWu2018_C = TWu2018_C;

end

function mineral = prepareMineralRow(dataTable, mineralLabel)
% prepareMineralRow
% Extract one-row mineral cation data. Missing optional variables are set to
% zero; existing NaN values are retained and never converted to zero.

if height(dataTable) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();

mineral.Fe2 = getVarOrError(dataTable, 'Fe_cation_apfu', mineralLabel);
mineral.Mg = getVarOrError(dataTable, 'Mg_cation_apfu', mineralLabel);

mineral.Fe3 = getVarOrZero(dataTable, 'Fe3_cation_apfu', mineralLabel);
mineral.Mn = getVarOrZero(dataTable, 'Mn_cation_apfu', mineralLabel);
mineral.Ca = getVarOrZero(dataTable, 'Ca_cation_apfu', mineralLabel);
mineral.Ti = getVarOrZero(dataTable, 'Ti_cation_apfu', mineralLabel);
mineral.Al = getVarOrZero(dataTable, 'Al_cation_apfu', mineralLabel);
mineral.Si = getVarOrZero(dataTable, 'Si_cation_apfu', mineralLabel);
mineral.K = getVarOrZero(dataTable, 'K_cation_apfu', mineralLabel);
mineral.Na = getVarOrZero(dataTable, 'Na_cation_apfu', mineralLabel);

% Eq. (14) and Table 1 use Fe2+. Fe3+ is retained separately but is not
% added to the value recorded as FeUsed.
mineral.FeUsed = mineral.Fe2;

end

function value = getVarOrError(dataTable, variableName, mineralLabel)
% getVarOrError
% Read a required scalar cation value. NaN and zero are retained; negative,
% infinite, nonnumeric, nonreal, and nonscalar values are rejected.

if ~ismember(variableName, dataTable.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = dataTable.(variableName);

if isInvalidCationValue(value)
    error('%s contains an invalid value for %s.', mineralLabel, variableName);
end

end

function value = getVarOrZero(dataTable, variableName, mineralLabel)
% getVarOrZero
% Read an optional scalar cation value. A missing variable is set to zero,
% whereas an existing NaN value remains NaN.

if ismember(variableName, dataTable.Properties.VariableNames)
    value = dataTable.(variableName);

    if isInvalidCationValue(value)
        error('%s contains an invalid value for %s.', ...
            mineralLabel, variableName);
    end
else
    value = 0;
end

end

function printCompositionRangeMessages(row, selectedCode_grt, selectedCode_mica)
% printCompositionRangeMessages
% Print non-stopping messages when finite mineral compositions fall outside
% the ranges represented by the Wu (2018) calibration dataset.

XFe_grt = row.XFe_grt(1);
XMg_grt = row.XMg_grt(1);
XCa_grt = row.XCa_grt(1);
XMn_grt = row.XMn_grt(1);
Fe_mica = row.FeUsed_mica(1);
Mg_mica = row.Mg_mica(1);
AlVI_mica = row.AlVI_mica(1);

if isfinite(XFe_grt) && (XFe_grt < 0.53 || XFe_grt > 0.91)
    fprintf(2, ...
        ['WARNING: Garnet XFe = %.4g is outside the Wu (2018) ' ...
         'calibration-composition range 0.53-0.91 for %s & %s.\n'], ...
        XFe_grt, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(XMg_grt) && (XMg_grt < 0.04 || XMg_grt > 0.27)
    fprintf(2, ...
        ['WARNING: Garnet XMg = %.4g is outside the Wu (2018) ' ...
         'calibration-composition range 0.04-0.27 for %s & %s.\n'], ...
        XMg_grt, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(XCa_grt) && (XCa_grt < 0.02 || XCa_grt > 0.22)
    fprintf(2, ...
        ['WARNING: Garnet XCa = %.4g is outside the Wu (2018) ' ...
         'calibration-composition range 0.02-0.22 for %s & %s.\n'], ...
        XCa_grt, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(XMn_grt) && (XMn_grt < 0.00 || XMn_grt > 0.29)
    fprintf(2, ...
        ['WARNING: Garnet XMn = %.4g is outside the Wu (2018) ' ...
         'calibration-composition range 0.00-0.29 for %s & %s. ' ...
         'Wu (2018, p. 984) reports invalid garnet activity-model behavior ' ...
         'and unexpectedly high temperatures for MnO-rich garnet beyond ' ...
         'the calibrated XMn limit.\n'], ...
        XMn_grt, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(Fe_mica) && (Fe_mica < 0.03 || Fe_mica > 0.21)
    fprintf(2, ...
        ['WARNING: Muscovite Fe = %.4g apfu is outside the Wu (2018) ' ...
         'calibration-composition range 0.03-0.21 apfu ' ...
         '(11-oxygen basis) for %s & %s.\n'], ...
        Fe_mica, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(Mg_mica) && (Mg_mica < 0.02 || Mg_mica > 0.32)
    fprintf(2, ...
        ['WARNING: Muscovite Mg = %.4g apfu is outside the Wu (2018) ' ...
         'calibration-composition range 0.02-0.32 apfu ' ...
         '(11-oxygen basis) for %s & %s.\n'], ...
        Mg_mica, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(AlVI_mica) && (AlVI_mica < 1.62 || AlVI_mica > 1.96)
    fprintf(2, ...
        ['WARNING: Muscovite AlVI = %.4g apfu is outside the Wu (2018) ' ...
         'calibration-composition range 1.62-1.96 apfu ' ...
         '(11-oxygen basis) for %s & %s.\n'], ...
        AlVI_mica, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

end

function act = calcWu2004bGarnetActivityTerms(XFe, XMg, XCa, XMn)
% calcWu2004bGarnetActivityTerms
% Garnet activity polynomial expressions from Wu et al. (2004),
% Journal of Petrology, 45, 1907-1921, Appendix A4-A9.

act = struct();

% ---------------------------------------------------------------------
% Fea : Appendix A4
% ---------------------------------------------------------------------
Fea = ...
    59.93 .* (XMg.^2) ...
    - 9.67 .* (XCa.^2) ...
    - 11.006 .* XFe .* XMg ...
    + 0.674 .* XFe .* XCa ...
    - 95.725 .* XMg .* XCa ...
    - 22.765 .* XMg .* XMn ...
    - 46.665 .* XCa .* XMn ...
    + 11.006 .* (XFe.^2) .* XMg ...
    - 0.674 .* (XFe.^2) .* XCa ...
    - 11.986 .* (XMg.^2) .* XFe ...
    + 37.966 .* (XMg.^2) .* XCa ...
    + 4.602 .* (XMg.^2) .* XMn ...
    + 1.934 .* (XCa.^2) .* XFe ...
    + 25.746 .* (XCa.^2) .* XMg ...
    + 4.602 .* (XMn.^2) .* XMg ...
    + 19.145 .* XFe .* XMg .* XCa ...
    + 4.553 .* XFe .* XMg .* XMn ...
    + 9.333 .* XFe .* XCa .* XMn ...
    + 77.876 .* XMg .* XCa .* XMn;

% ---------------------------------------------------------------------
% Feb : Appendix A5
% ---------------------------------------------------------------------
Feb = ...
    -0.034 .* (XMg.^2) ...
    + 0.135 .* (XCa.^2) ...
    + 0.024 .* (XMn.^2) ...
    + 0.1 .* XFe .* XMg ...
    + 0.08 .* XFe .* XCa ...
    + 0.048 .* XFe .* XMn ...
    - 0.0515 .* XMg .* XCa ...
    + 0.007 .* XMg .* XMn ...
    + 0.1765 .* XCa .* XMn ...
    - 0.1 .* (XFe.^2) .* XMg ...
    - 0.08 .* (XFe.^2) .* XCa ...
    - 0.048 .* (XFe.^2) .* XMn ...
    + 0.068 .* (XMg.^2) .* XFe ...
    - 0.136 .* (XMg.^2) .* XCa ...
    - 0.076 .* (XMg.^2) .* XMn ...
    - 0.27 .* (XCa.^2) .* XFe ...
    - 0.28 .* (XCa.^2) .* XMg ...
    - 0.13 .* (XCa.^2) .* XMn ...
    - 0.048 .* (XMn.^2) .* XFe ...
    - 0.076 .* (XMn.^2) .* XMg ...
    - 0.13 .* (XMn.^2) .* XCa ...
    + 0.103 .* XFe .* XMg .* XCa ...
    - 0.14 .* XFe .* XMg .* XMn ...
    - 0.353 .* XFe .* XCa .* XMn ...
    - 0.414 .* XMg .* XCa .* XMn;

% ---------------------------------------------------------------------
% Fec : Appendix A6
% ---------------------------------------------------------------------
Fec = ...
    -5672.0 .* (XMg.^2) ...
    + 19932.0 .* (XCa.^2) ...
    + 1617.0 .* (XMn.^2) ...
    + 23244.0 .* XFe .* XMg ...
    - 2608.0 .* XFe .* XCa ...
    + 3234.0 .* XFe .* XMn ...
    + 31326.5 .* XMg .* XCa ...
    + 45841.0 .* XMg .* XMn ...
    + 12356.0 .* XCa .* XMn ...
    - 23244.0 .* (XFe.^2) .* XMg ...
    + 2608.0 .* (XFe.^2) .* XCa ...
    - 3234.0 .* (XFe.^2) .* XMn ...
    + 11344.0 .* (XMg.^2) .* XFe ...
    - 132228.0 .* (XMg.^2) .* XCa ...
    - 82498.0 .* (XMg.^2) .* XMn ...
    - 39864.0 .* (XCa.^2) .* XFe ...
    - 51518.0 .* (XCa.^2) .* XMg ...
    - 2850.0 .* (XCa.^2) .* XMn ...
    - 3234.0 .* (XMn.^2) .* XFe ...
    - 82498.0 .* (XMn.^2) .* XMg ...
    - 2850.0 .* (XMn.^2) .* XCa ...
    - 62653.0 .* XFe .* XMg .* XCa ...
    - 91682.0 .* XFe .* XMg .* XMn ...
    - 24712.0 .* XFe .* XCa .* XMn ...
    - 177221.0 .* XMg .* XCa .* XMn;

% ---------------------------------------------------------------------
% Mga : Appendix A7
% ---------------------------------------------------------------------
Mga = ...
    55.03 .* (XFe.^2) ...
    - 128.73 .* (XCa.^2) ...
    - 23.01 .* (XMn.^2) ...
    + 11.986 .* XFe .* XMg ...
    - 95.725 .* XFe .* XCa ...
    - 22.765 .* XFe .* XMn ...
    - 37.966 .* XMg .* XCa ...
    - 4.602 .* XMg .* XMn ...
    - 38.938 .* XCa .* XMn ...
    + 11.006 .* (XFe.^2) .* XMg ...
    - 0.674 .* (XFe.^2) .* XCa ...
    - 11.986 .* (XMg.^2) .* XFe ...
    + 37.966 .* (XMg.^2) .* XCa ...
    + 4.602 .* (XMg.^2) .* XMn ...
    + 1.934 .* (XCa.^2) .* XFe ...
    + 25.746 .* (XCa.^2) .* XMg ...
    + 4.602 .* (XMn.^2) .* XMg ...
    + 19.145 .* XFe .* XMg .* XCa ...
    + 4.553 .* XFe .* XMg .* XMn ...
    + 9.333 .* XFe .* XCa .* XMn ...
    + 77.876 .* XMg .* XCa .* XMn;

% ---------------------------------------------------------------------
% Mgb : Appendix A8
% ---------------------------------------------------------------------
Mgb = ...
    0.05 .* (XFe.^2) ...
    + 0.14 .* (XCa.^2) ...
    + 0.038 .* (XMn.^2) ...
    - 0.068 .* XFe .* XMg ...
    - 0.0515 .* XFe .* XCa ...
    + 0.007 .* XFe .* XMn ...
    - 0.136 .* XMg .* XCa ...
    + 0.076 .* XMg .* XMn ...
    + 0.207 .* XCa .* XMn ...
    - 0.1 .* (XFe.^2) .* XMg ...
    - 0.08 .* (XFe.^2) .* XCa ...
    - 0.048 .* (XFe.^2) .* XMn ...
    + 0.068 .* (XMg.^2) .* XFe ...
    - 0.136 .* (XMg.^2) .* XCa ...
    - 0.076 .* (XMg.^2) .* XMn ...
    - 0.27 .* (XCa.^2) .* XFe ...
    - 0.28 .* (XCa.^2) .* XMg ...
    - 0.13 .* (XCa.^2) .* XMn ...
    - 0.048 .* (XMn.^2) .* XFe ...
    - 0.076 .* (XMn.^2) .* XMg ...
    - 0.13 .* (XMn.^2) .* XCa ...
    + 0.103 .* XFe .* XMg .* XCa ...
    - 0.14 .* XFe .* XMg .* XMn ...
    - 0.353 .* XFe .* XCa .* XMn ...
    - 0.414 .* XMg .* XCa .* XMn;

% ---------------------------------------------------------------------
% Mgc : Appendix A9
% ---------------------------------------------------------------------
Mgc = ...
    11622.0 .* (XFe.^2) ...
    + 25759.0 .* (XCa.^2) ...
    + 41249.0 .* (XMn.^2) ...
    - 11344.0 .* XFe .* XMg ...
    + 31326.5 .* XFe .* XCa ...
    + 45841.0 .* XFe .* XMn ...
    + 132228.0 .* XMg .* XCa ...
    + 82498.0 .* XMg .* XMn ...
    + 88610.5 .* XCa .* XMn ...
    - 23244.0 .* (XFe.^2) .* XMg ...
    + 2608.0 .* (XFe.^2) .* XCa ...
    - 3234.0 .* (XFe.^2) .* XMn ...
    + 11344.0 .* (XMg.^2) .* XFe ...
    - 132228.0 .* (XMg.^2) .* XCa ...
    - 82498.0 .* (XMg.^2) .* XMn ...
    - 39864.0 .* (XCa.^2) .* XFe ...
    - 51518.0 .* (XCa.^2) .* XMg ...
    - 2850.0 .* (XCa.^2) .* XMn ...
    - 3234.0 .* (XMn.^2) .* XFe ...
    - 82498.0 .* (XMn.^2) .* XMg ...
    - 2850.0 .* (XMn.^2) .* XCa ...
    - 62653.0 .* XFe .* XMg .* XCa ...
    - 91682.0 .* XFe .* XMg .* XMn ...
    - 24712.0 .* XFe .* XCa .* XMn ...
    - 177221.0 .* XMg .* XCa .* XMn;

act.Fea = Fea;
act.Feb = Feb;
act.Fec = Fec;
act.Mga = Mga;
act.Mgb = Mgb;
act.Mgc = Mgc;

end
