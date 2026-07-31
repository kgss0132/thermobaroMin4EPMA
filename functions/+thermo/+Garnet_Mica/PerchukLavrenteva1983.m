function results = PerchukLavrenteva1983(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/PerchukLavrenteva1983.m
% Tested with MATLAB R2024b
%
% Garnet-biotite Fe-Mg exchange thermometer
% Perchuk, L.L. and Lavrent'eva, I.V. (1983)
% Experimental Investigation of Exchange Equilibria in the System
% Cordierite-Garnet-Biotite. In: Kinetics and Equilibrium in Mineral
% Reactions, Advances in Physical Geochemistry, vol. 3, pp. 199-239.
% DOI: https://doi.org/10.1007/978-1-4612-5587-1_7
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% analysis (selected by the user from tables) and calculates temperature
% using the Perchuk & Lavrent'eva (1983) garnet-biotite Fe-Mg exchange
% thermometer (their Eq. 31).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Garnet-Mica pair, and appends results
% into a single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Perchuk & Lavrent'eva (1983) report the following conditions:
%
%   Experimental-method range for the chapter : 550-1000 degreeC
%                                                and 5-7 kbar
%   Grt-Bt calibration data in Tables 13-14    : 575-950 degreeC
%                                                at approximately 6 kbar
%   Overall composition range in Table 13      : N_Mg(Bt) = 22.4-83 mol%
%                                                N_Mg(Grt) = 4.9-76.5 mol%
%   Intended application                       : equilibrated Grt-Bt pairs
%                                                in metapelites
%
% The general experimental conditions are described on p. 202. The Grt-Bt
% results and the 575-950 degreeC, 6 kbar calibration dataset are presented
% on pp. 226-234 (especially Tables 13-14 on pp. 233-234). Equation (31) is
% presented on p. 235.
%
% Important application cautions stated or demonstrated in the chapter:
% - The authors treat Fe-Mg distribution as ideal above 600 degreeC; the
%   lowest-temperature part of the calibration therefore requires particular
%   caution (pp. 202 and 227-228).
% - Equation (31) is based on a limited and uneven N_Mg range. The authors
%   caution that K_D may depend on composition in very Fe-rich or very
%   Mg-rich systems and that experiments at different Bi Al contents and
%   pressures are still needed (p. 235).
% - Newly synthesized garnets were close to binary pyrope-almandine and
%   contained very little Ca and Mn. Application to Ca- or Mn-rich garnet is
%   therefore an extrapolation beyond the best-constrained compositions
%   (p. 226).
% - The pressure correction is uncertain because the chapter discusses two
%   estimates of DeltaV: -0.0246 and -0.0577 cal/bar per exchanging atom
%   (pp. 232-235). This implementation uses -0.0246 cal/bar, the value from
%   Hewitt and Wones (1975) adopted in common implementations of Eq. (31).
%   Pressure extrapolation away from 6 kbar should therefore be treated with
%   caution.
% - Fe-Mg diffusion is slow and both experimental and natural minerals may
%   preserve zoning or metastable compositions. Pair spatially adjacent
%   domains that represent the same equilibrium stage (for example,
%   core-core or rim-rim), rather than unrelated grain averages
%   (pp. 207-208 and 237-238).
% - The authors report unrealistically low temperatures for some
%   phlogopite-garnet granulite assemblages and basalt-hosted inclusions
%   (p. 235).
%
% This implementation therefore issues non-stopping warnings when:
%   1) input pressure is outside 5-7 kbar, or
%   2) a finite calculated temperature is outside 575-950 degreeC.
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
%
% Optional variables retained in the output table when present:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Fe3_cation_apfu is retained for reference but is not added to Fe2+ in the
% Fe-Mg exchange coefficient. Missing optional variables are represented by
% zero. An optional variable that exists and contains NaN remains NaN.
%
% All finite mineral-composition values must be non-negative. Negative and
% infinite values stop the calculation. Zero is allowed as requested, even
% though it may produce a non-finite exchange coefficient or temperature.
% NaN values in required thermometer inputs are retained as missing values,
% propagated through the calculation, and reported by non-stopping warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Perchuk & Lavrent'eva, 1983, Eq. 31, p. 235)
%
%   K_D = (Fe/Mg)_Grt / (Fe/Mg)_Bt
%       = (Fe_Grt * Mg_Bt) / (Mg_Grt * Fe_Bt)
%
%   T(K) = [7843.7 + DeltaV * (P_bar - 6000)] ...
%          / [1.987 * ln(K_D) + 5.699]
%
%   T(degreeC) = T(K) - 273.15
%
% where:
%   P_bar  = P_kbar * 1000
%   DeltaV = -0.0246 cal/bar per exchanging atom
%
% The 6000-bar reference pressure is the pressure at which the experimental
% Grt-Bt regression was obtained. No separate Mn correction is added to
% Eq. (31).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = PerchukLavrenteva1983(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Mica tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Mica pair. The output variable set is
%             intended to remain stable for downstream processing.
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('PerchukLavrenteva1983 requires (rawdata_struct, P_kbar).');
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

requiredVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
validateRequiredVariables(dataset_grt, requiredVariables, 'Garnet');
validateRequiredVariables(dataset_mica, requiredVariables, 'Mica');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated table concatenation inside the loop is avoided because it forces
% MATLAB to repeatedly reallocate and copy the entire results table.
%
% The cell buffer is preallocated and doubled only when necessary. After the
% interactive loop finishes, all blocks are concatenated once with vertcat.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Calibration limits used for the non-stopping range messages.
calibrationT_min_degC = 575;
calibrationT_max_degC = 950;
calibrationP_min_kbar = 5;
calibrationP_max_kbar = 7;

% Pressure is common to all selected mineral pairs in this function call.
% The warning is therefore printed only once, after the first calculation.
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
    % Assumption: the first column stores an identifier displayed to the user.
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

    % Check only the variables actually used in Eq. (31). The calculation is
    % intentionally allowed to continue when NaN is present.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_mica);

    % Negative and infinite required inputs stop the calculation. NaN and
    % zero are intentionally allowed and are handled after calculation.
    validateNonNegativeInputs(selectedData_grt, selectedData_mica);

    row = calcTemp(selectedData_grt, selectedData_mica, P_kbar);

    % Store the selected identifiers for each pressure point.
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_mica = repmat(string(selectedCode_mica), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_mica'}, 'Before', 1);

    % Store this result as one table block. If the preallocated buffer is
    % full, double its capacity rather than growing the results table during
    % every loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': ' num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside 5-7 kbar. Equation (31)
    % is referenced to 6 kbar, and the pressure correction remains uncertain.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental pressure range ' ...
             'reported by Perchuk & Lavrent''eva (1983): 5-7 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar. Eq. (31) is referenced to 6 kbar, ' ...
             'and its pressure correction is uncertain.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the Grt-Bt
    % calibration dataset range of 575-950 degreeC. NaN and Inf are handled
    % separately below.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the experimental Grt-Bt ' ...
             'calibration range of Perchuk & Lavrent''eva (1983): ' ...
             '575-950 degreeC. %d of %d finite temperature point(s) are ' ...
             'outside the range; calculated finite range = %.4g-%.4g degreeC ' ...
             'for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)));
    end

    % Print a non-stopping warning immediately after the temperature result
    % when a required Eq. (31) input contains NaN. fprintf is used so that
    % the message remains visible even when MATLAB warnings are disabled.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, and NaN values were not converted to zero.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.T_deg);

    % Retain and report NaN/Inf results caused by NaN, zero denominators, or
    % another mathematically non-finite operation.
    if any(invalidTemperature)
        fprintf(2, ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
                    '(%d of %d points; NaN: %d, Inf: %d).\n' ...
                    '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'PerchukLavrenteva1983', ...
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
% Confirm that every variable required by Eq. (31) exists before the user
% enters the interactive selection loop.

missingVariables = requiredVariables(~ismember( ...
    requiredVariables, dataTable.Properties.VariableNames));

if ~isempty(missingVariables)
    error('%s table must contain variable(s): %s', ...
        mineralLabel, char(strjoin(string(missingVariables), ', ')));
end

end

function nanInputNames = findNaNInputs(data_garnet, data_mica)
% findNaNInputs
% Return the names of required Eq. (31) input variables that contain NaN.
% The output buffer is preallocated so that its size does not grow during
% the variable-checking loops.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

nanInputNames = strings(numel(garnetVariables) + numel(micaVariables), 1);
nNaNInputs = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    variableValue = data_garnet.(variableName);
    if isnumeric(variableValue) && any(isnan(variableValue(:)))
        nNaNInputs = nNaNInputs + 1;
        nanInputNames(nNaNInputs) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
    variableValue = data_mica.(variableName);
    if isnumeric(variableValue) && any(isnan(variableValue(:)))
        nNaNInputs = nNaNInputs + 1;
        nanInputNames(nNaNInputs) = "Mica." + string(variableName);
    end
end

nanInputNames = nanInputNames(1:nNaNInputs);

end

function validateNonNegativeInputs(data_garnet, data_mica)
% validateNonNegativeInputs
% Stop the calculation when a required Eq. (31) input is negative, infinite,
% nonnumeric, nonreal, or nonscalar. NaN and zero are intentionally allowed.
% The error-name buffer is preallocated and does not grow inside the loops.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

invalidInputNames = strings(numel(garnetVariables) + numel(micaVariables), 1);
nInvalidInputs = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    variableValue = data_garnet.(variableName);
    if isInvalidCationValue(variableValue)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputNames(nInvalidInputs) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
    variableValue = data_mica.(variableName);
    if isInvalidCationValue(variableValue)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputNames(nInvalidInputs) = "Mica." + string(variableName);
    end
end

invalidInputNames = invalidInputNames(1:nInvalidInputs);

if ~isempty(invalidInputNames)
    error(['PerchukLavrenteva1983: required mineral-composition values ' ...
           'must be real, scalar, non-negative, and not Inf. Invalid value(s) ' ...
           'were found in: ' char(strjoin(invalidInputNames, ', ')) '.']);
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

function row = calcTemp(data_garnet, data_mica, P_kbar)
% calcTemp
% Compute Perchuk & Lavrent'eva (1983) Eq. (31) for one garnet row, one
% mica row, and every pressure supplied in P_kbar.
%
% Inputs:
%   data_garnet : 1-row table containing garnet apfu cations
%   data_mica   : 1-row table containing mica apfu cations
%   P_kbar      : pressure in kbar (column vector after input validation)
%
% Output:
%   row : table with one row per pressure value, including mineral inputs,
%         exchange variables, constants, and calculated temperatures.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

row = table();

% --- Published constants used in Eq. (31) ---
R_cal = 1.987;
deltaV_cal_bar = -0.0246;
referenceP_bar = 6000;

% --- Prepare one-row mineral data ---
grt = prepareMineralRow(data_garnet, 'Garnet');
mica = prepareMineralRow(data_mica, 'Mica');

% --- Replicate scalar mineral inputs for every pressure point ---
Fe2_grt = repmat(grt.Fe2, nP, 1);
Fe3_grt = repmat(grt.Fe3, nP, 1);
Mg_grt = repmat(grt.Mg, nP, 1);
Mn_grt = repmat(grt.Mn, nP, 1);
Ca_grt = repmat(grt.Ca, nP, 1);
Al_grt = repmat(grt.Al, nP, 1);
Si_grt = repmat(grt.Si, nP, 1);

Fe2_mica = repmat(mica.Fe2, nP, 1);
Fe3_mica = repmat(mica.Fe3, nP, 1);
Mg_mica = repmat(mica.Mg, nP, 1);
Mn_mica = repmat(mica.Mn, nP, 1);
Ca_mica = repmat(mica.Ca, nP, 1);
Ti_mica = repmat(mica.Ti, nP, 1);
Al_mica = repmat(mica.Al, nP, 1);
Si_mica = repmat(mica.Si, nP, 1);
K_mica = repmat(mica.K, nP, 1);
Na_mica = repmat(mica.Na, nP, 1);

% FeUsed records the Fe value actually used by Eq. (31). The exchange
% reaction is formulated with Fe2+; Fe3+ is retained separately but is not
% included in FeUsed.
FeUsed_grt = Fe2_grt;
FeUsed_mica = Fe2_mica;

% --- Composition indices and Fe-Mg exchange coefficient ---
FeMg_grt = FeUsed_grt ./ Mg_grt;
FeMg_mica = FeUsed_mica ./ Mg_mica;

KD2 = FeMg_grt ./ FeMg_mica;
lnKD2 = log(KD2);

NMg_grt = 100 .* Mg_grt ./ (Mg_grt + FeUsed_grt + Mn_grt);
NMg_mica = 100 .* Mg_mica ./ (Mg_mica + FeUsed_mica + Mn_mica);
XGr_Mn = Mn_grt ./ (Mn_grt + FeUsed_grt + Mg_grt);

% --- Solve temperature using Eq. (31) ---
equationNumerator = 7843.7 + deltaV_cal_bar .* ...
    (P_bar - referenceP_bar);
equationDenominator = R_cal .* lnKD2 + 5.699;

T_K = equationNumerator ./ equationDenominator;
T_deg = T_K - 273.15;

% --- Pack pressure and equation constants ---
row.P_kbar = P_kbar;
row.P_bar = P_bar;
row.R_cal = repmat(R_cal, nP, 1);
row.deltaV_cal_bar = repmat(deltaV_cal_bar, nP, 1);
row.referenceP_bar = repmat(referenceP_bar, nP, 1);

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

% --- Pack calculated indices and temperatures ---
row.FeMg_grt = FeMg_grt;
row.FeMg_mica = FeMg_mica;
row.KD2 = KD2;
row.lnKD2 = lnKD2;
row.NMg_grt = NMg_grt;
row.NMg_mica = NMg_mica;
row.XGr_Mn = XGr_Mn;
row.equationNumerator = equationNumerator;
row.equationDenominator = equationDenominator;
row.T_K = T_K;
row.T_deg = T_deg;

% Retain the former output name as an alias for backward compatibility.
row.T_C = T_deg;

end

function mineral = prepareMineralRow(dataTable, mineralLabel)
% prepareMineralRow
% Extract one-row mineral cation data. Missing optional variables are set to
% zero; existing NaN values are retained and are never converted to zero.

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

% FeUsed is defined explicitly as the ferrous Fe used by Eq. (31).
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
