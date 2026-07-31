function results = RidolfiRenzulli2012(rawdata_struct, P_kbar)
% functions/+thermo/+Amphibole/RidolfiRenzulli2012.m
% Tested with MATLAB R2024b
%
% Single-amphibole thermometer
% Ridolfi, F. & Renzulli, A. (2012)
% Contributions to Mineralogy and Petrology, 163, 877–895
% DOI: https://doi.org/10.1007/s00410-011-0704-6
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Amphibole analysis at a time and
% calculates temperature using Eq. (2) of Ridolfi and Renzulli (2012).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Amphibole analysis and appends the
% result into a single output table.
%
% Pressure is NOT estimated internally from amphibole composition in this
% implementation. Instead, the user-supplied pressure P_kbar is used. A
% scalar pressure produces one result row, whereas a pressure vector produces
% one result row for each pressure value. This allows the function to operate
% with both startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Ridolfi and Renzulli (2012) calibrated Eq. (2) using 61 Mg-rich calcic
% amphiboles synthesized in equilibrium crystallization experiments over the
% following range:
%
%   Temperature : 800–1130 degreeC
%   Pressure    : 130–2200 MPa (1.3–22 kbar; 0.13–2.2 GPa)
%   Amphibole   : Mg-rich calcic amphiboles from calc-alkaline and alkaline
%                 magmas, with compositions broadly matching Fig. 1 and
%                 Table 2
%   Formula     : total Si, Ti, Al, Fe, Mg, Ca, Na, and K calculated on a
%                 13-cation basis
%
% The experimental range and data-selection criteria are described on
% pp. 879–884, Table 1, Table 2, and Fig. 3. Equation (2) is given in Table 3
% on p. 887. Its calibration statistics are shown in Fig. 5g and discussed
% on pp. 888–890: standard error of estimate = 23.5 degreeC and maximum
% absolute calibration error = approximately 50 degreeC.
%
% IMPORTANT APPLICATION NOTES:
% - The valid experimental P–T region is not a simple rectangle. The selected
%   experiments define a horn-shaped field in Fig. 3a. Ridolfi and Renzulli
%   (2012) recommend checking that calculated P–T conditions fall within, or
%   near the boundary of, this field after accounting for expected uncertainty
%   (pp. 891–892).
% - Natural amphibole compositions should broadly match the experimental
%   amphiboles in Fig. 1 and Table 2 before the equation is applied
%   (pp. 891–892).
% - The method is intended for magmatic amphibole phenocrysts, including
%   crystals with dehydration or breakdown rims, and for euhedral,
%   compositionally homogeneous amphiboles from sub-volcanic or plutonic
%   rocks (pp. 891–892).
% - The method should not be applied to hydrothermal vein amphiboles,
%   microlites, rapidly grown or quenched amphibole zones, or clearly
%   disequilibrium compositions (pp. 891–892).
% - The original internally consistent procedure first estimates pressure
%   from amphibole composition using Eqs. (1a–e), then applies Eq. (2).
%   This implementation instead accepts an independently supplied pressure.
% - Results outside the experimental P–T or compositional fields may have
%   larger and unpredictable uncertainties (p. 891 and final remarks,
%   p. 892).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 1.3–22 kbar,
%   2) a finite calculated temperature is outside 800–1130 degreeC,
%   3) a required thermometer input contains NaN, or
%   4) the calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following normalized cation variables on a 13-cation basis:
%
%   Si_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu       % total Fe
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%
% All finite values in the required amphibole-composition variables above
% must be greater than or equal to zero. Negative values are prohibited.
% NaN values are retained as missing values, propagated through the
% calculation, and reported by non-stopping fprintf warnings. NaN is never
% replaced by zero.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Ridolfi and Renzulli (2012), Eq. (2):
%
%   T (degreeC) =
%       17098
%     - 1322.3 * Si
%     - 1035.1 * Ti
%     - 1208.2 * Al
%     - 1230.4 * Fe
%     - 1152.9 * Mg
%     -  130.40 * Ca
%     +  200.54 * Na
%     +   29.408 * K
%     +   24.410 * ln(P_MPa)
%
% where:
%   - Si, Ti, Al, Fe, Mg, Ca, Na, and K are amphibole total cations
%     calculated on a 13-cation basis
%   - Fe is total Fe
%   - P_MPa is pressure in MPa
%
% Pressure conversion:
%   P_MPa = 100 * P_kbar
%
% Because the equation contains ln(P_MPa), every pressure value must be
% strictly greater than zero.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = RidolfiRenzulli2012(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing an Amphibole table (see above)
%   P_kbar         : pressure in kbar (finite positive numeric scalar or
%                    vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Amphibole analysis. The output variable set is
%             intended to remain stable for downstream processing.
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% mathematically invalid pressure values.
if nargin < 2
    error('RidolfiRenzulli2012 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) <= 0)
    error('P_kbar must be a finite positive numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation dataset
% Extract the required table from the input struct. The table is not modified;
% only the required columns are read during calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Amphibole') || ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end

dataset_amp = rawdata_struct.Amphibole;

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

% Experimental calibration limits reported by Ridolfi and Renzulli (2012).
calibrationT_min_degC = 800;
calibrationT_max_degC = 1130;
calibrationP_min_kbar = 1.3;
calibrationP_max_kbar = 22;

% Pressure is common to all selected Amphibole analyses in this function
% call. The warning is therefore printed only once, after the first
% calculation.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % ----- Amphibole selection -----
    % Assumption: the first column stores an identifier to display to the user.
    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, ok] = listdlg( ...
        'PromptString', 'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_amp)), ...
        'ListSize', [320 320]);

    % If the user cancels, exit the loop gracefully.
    if ~ok || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    % Store the selected input row so that required thermometer inputs can be
    % checked explicitly without modifying the source table.
    selectedData_amp = dataset_amp(selectedIdx_amp, :);

    % Check only the variables actually used in Eq. (2). Calculation is
    % intentionally allowed to continue when NaN is present; the NaN values
    % remain unchanged and a warning is printed after the result.
    nanInputNames = findNaNInputs(selectedData_amp);

    % All finite cation inputs must be non-negative. NaN values are excluded
    % from this check so that they propagate through the calculation.
    validateNonnegativeInputs(selectedData_amp);

    row = calcTemp(selectedData_amp, P_kbar);

    % Store the user-selected identifier for traceability.
    row.dataCode_amphibole = ...
        repmat(string(selectedCode_amp), height(row), 1);

    % Move the identifier to the front for readability.
    row = movevars(row, {'dataCode_amphibole'}, 'Before', 1);

    % Store this result as one table block. If the preallocated buffer is
    % full, double its capacity. This occurs only occasionally and avoids
    % reallocating the full results table after every calculation.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_amp)) ': T = ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_amp)) ': T = ' ...
            num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside the experimental
    % calibration range of 1.3–22 kbar. The calculation is not stopped.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental calibration range ' ...
             'of Ridolfi and Renzulli (2012): 1.3–22 kbar ' ...
             '(130–2200 MPa; 0.13–2.2 GPa). ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the
    % experimental calibration range of 800–1130 degreeC. NaN and Inf are
    % handled separately by the non-finite-result warning below.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the experimental calibration ' ...
             'range of Ridolfi and Renzulli (2012): 800–1130 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_amp)));
    end

    % Print a non-stopping warning immediately after the temperature result
    % when any required thermometer input contains NaN. fprintf is used
    % instead of warning so that the message remains visible even when MATLAB
    % warnings have been disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         The calculation was continued, NaN was not replaced by zero, ' ...
             'and the calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_amp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.T_deg);

    % Retain a result-based check for NaN/Inf values caused by an explicitly
    % NaN input or by any other numerical issue. The result remains in the
    % output table and calculation is not stopped.
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_amp)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'RidolfiRenzulli2012', ...
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
function nanInputNames = findNaNInputs(data_amphibole)
% findNaNInputs
% Return the names of required thermometer input variables containing NaN.
% This function does not throw an error for NaN values; it only prepares a
% warning message for the calling function.

amphiboleVariables = { ...
    'Si_cation_apfu'; ...
    'Ti_cation_apfu'; ...
    'Al_cation_apfu'; ...
    'Fe_cation_apfu'; ...
    'Mg_cation_apfu'; ...
    'Ca_cation_apfu'; ...
    'Na_cation_apfu'; ...
    'K_cation_apfu'};

nVariables = numel(amphiboleVariables);
hasNaN = false(nVariables, 1);

for i = 1:nVariables
    variableName = amphiboleVariables{i};

    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        error('Amphibole table must contain variable: %s', variableName);
    end

    variableValue = data_amphibole.(variableName);
    hasNaN(i) = any(isnan(variableValue(:)));
end

nanInputNames = "Amphibole." + string(amphiboleVariables(hasNaN));

end

function validateNonnegativeInputs(data_amphibole)
% validateNonnegativeInputs
% Stop the calculation when a finite required amphibole-composition value is
% negative or when an infinite value is present. Zero is allowed. NaN is
% intentionally allowed so that it remains NaN, propagates through Eq. (2),
% and is reported by non-stopping warning messages.

amphiboleVariables = { ...
    'Si_cation_apfu'; ...
    'Ti_cation_apfu'; ...
    'Al_cation_apfu'; ...
    'Fe_cation_apfu'; ...
    'Mg_cation_apfu'; ...
    'Ca_cation_apfu'; ...
    'Na_cation_apfu'; ...
    'K_cation_apfu'};

nVariables = numel(amphiboleVariables);
hasInvalidValue = false(nVariables, 1);

for i = 1:nVariables
    variableName = amphiboleVariables{i};

    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        error('Amphibole table must contain variable: %s', variableName);
    end

    variableValue = data_amphibole.(variableName);
    hasInvalidValue(i) = any(isinf(variableValue(:)) | ...
        (isfinite(variableValue(:)) & variableValue(:) < 0));
end

if any(hasInvalidValue)
    invalidInputNames = ...
        "Amphibole." + string(amphiboleVariables(hasInvalidValue));

    error(['RidolfiRenzulli2012: required amphibole-composition values ' ...
           'must be non-negative or NaN. Negative or infinite value(s) were ' ...
           'found in: ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_amphibole, P_kbar)
% calcTemp
% Compute Ridolfi and Renzulli (2012) Eq. (2) for one Amphibole row and a
% scalar or vector of user-supplied pressure values.
%
% Inputs:
%   data_amphibole : 1-row table containing 13-cation-basis amphibole data
%   P_kbar         : pressure in kbar; converted internally to MPa
%
% Output:
%   row : table containing one row per pressure value, input cations,
%         pressure variables, calculated temperature, and calibration flags.

P_kbar = P_kbar(:);
P_MPa = 100 .* P_kbar;
lnP_MPa = log(P_MPa);
nP = numel(P_kbar);

% Extract one Amphibole composition without replacing NaN values.
amp = prepareAmphiboleRow(data_amphibole);

% Replicate the selected Amphibole composition to match the number of
% pressure points. This produces one output row for each pressure value.
Si = repmat(amp.Si, nP, 1);
Ti = repmat(amp.Ti, nP, 1);
Al = repmat(amp.Al, nP, 1);
Fe = repmat(amp.FeT, nP, 1);
Mg = repmat(amp.Mg, nP, 1);
Ca = repmat(amp.Ca, nP, 1);
Na = repmat(amp.Na, nP, 1);
K  = repmat(amp.K, nP, 1);

% --- Ridolfi and Renzulli (2012), Eq. (2) ---
T_deg = ...
    17098 ...
    - 1322.3 .* Si ...
    - 1035.1 .* Ti ...
    - 1208.2 .* Al ...
    - 1230.4 .* Fe ...
    - 1152.9 .* Mg ...
    - 130.40 .* Ca ...
    + 200.54 .* Na ...
    + 29.408 .* K ...
    + 24.410 .* lnP_MPa;

% Calibration flags use the reported experimental limits. The original
% P–T field is horn-shaped rather than rectangular; these flags therefore
% indicate only whether T and P separately fall within their reported ranges.
isWithinCalibration_T = isfinite(T_deg) & T_deg >= 800 & T_deg <= 1130;
isWithinCalibration_P = P_kbar >= 1.3 & P_kbar <= 22;
isWithinCalibration_TP = isWithinCalibration_T & isWithinCalibration_P;

% --- Pack outputs ---
row = table();

row.P_kbar = P_kbar;
row.P_MPa = P_MPa;
row.lnP_MPa = lnP_MPa;

row.Si_amp = Si;
row.Ti_amp = Ti;
row.Al_amp = Al;
row.Fe_amp = Fe;
row.Mg_amp = Mg;
row.Ca_amp = Ca;
row.Na_amp = Na;
row.K_amp = K;

row.T_deg = T_deg;
row.isWithinCalibration_T = isWithinCalibration_T;
row.isWithinCalibration_P = isWithinCalibration_P;
row.isWithinCalibration_TP = isWithinCalibration_TP;

end

function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Extract one Amphibole analysis from a 1-row table. NaN values are retained
% exactly as supplied and are never converted to zero.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();

amp.Si  = getRequiredVariable(data_amphibole, 'Si_cation_apfu', 'Amphibole');
amp.Ti  = getRequiredVariable(data_amphibole, 'Ti_cation_apfu', 'Amphibole');
amp.Al  = getRequiredVariable(data_amphibole, 'Al_cation_apfu', 'Amphibole');
amp.FeT = getRequiredVariable(data_amphibole, 'Fe_cation_apfu', 'Amphibole');
amp.Mg  = getRequiredVariable(data_amphibole, 'Mg_cation_apfu', 'Amphibole');
amp.Ca  = getRequiredVariable(data_amphibole, 'Ca_cation_apfu', 'Amphibole');
amp.Na  = getRequiredVariable(data_amphibole, 'Na_cation_apfu', 'Amphibole');
amp.K   = getRequiredVariable(data_amphibole, 'K_cation_apfu', 'Amphibole');

end

function value = getRequiredVariable(tbl, variableName, mineralLabel)
% getRequiredVariable
% Read one required scalar variable from a 1-row mineral table. NaN is
% accepted and returned unchanged. Negative and infinite values are checked
% separately by validateNonnegativeInputs.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);

if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s must be a numeric scalar in a 1-row table.', ...
        variableName);
end

end
