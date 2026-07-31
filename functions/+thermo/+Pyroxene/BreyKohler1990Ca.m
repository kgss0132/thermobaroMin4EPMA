function results = BreyKohler1990Ca(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/BreyKohler1990Ca.m
% Tested with MATLAB R2024b
%
% Ca-in-Opx thermometer for peridotitic systems
% Brey, G.P., Kohler, T. (1990)
% Journal of Petrology, 31(6), 1353–1378
% DOI: https://doi.org/10.1093/petrology/31.6.1353
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Orthopyroxene (Opx) analysis and
% calculates temperature using the Ca-in-Opx thermometer of Brey & Kohler
% (1990).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Opx temperature and stores the results
% in a single output table.
%
% Both scalar and vector pressure inputs are accepted. For each selected Opx
% analysis, one output row is returned for every pressure value. This allows
% the function to be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Brey & Kohler (1990) fitted the Ca-in-Opx thermometer to reversed
% experiments in the CaO-MgO-SiO2 (CMS) system and tested the resulting
% equation against experiments in a natural peridotitic composition. The
% thermometer is presented as equation (10) on p. 1366.
%
% Experimental and validation limits relevant to this implementation are:
%
%   Temperature : approximately 900–1400 degreeC for the natural-system
%                 validation experiments shown in Fig. 7a
%   Pressure    : 2–60 kbar for the underlying CMS experimental data;
%                 10–60 kbar for the natural-system experiments
%   Composition : Mg-rich peridotitic orthopyroxene buffered by coexisting
%                 clinopyroxene
%
% The CMS and natural-system pressure intervals are reported in the
% discussion of the pyroxene calibrations on pp. 1362–1366. The natural-
% system temperature span is shown in Fig. 7a on p. 1365. Equation (10)
% reproduced the CMS experiments to approximately +/-26 degreeC (1 sigma)
% and the natural-system experiments to approximately +/-19 degreeC
% (1 sigma) (p. 1366).
%
% Although only an Opx analysis is supplied to this function, the method is
% not independent of mineral assemblage. Ca in orthopyroxene is interpreted
% as being buffered by coexisting clinopyroxene. Application to an Opx that
% did not coexist or equilibrate with Cpx may therefore be inappropriate.
%
% Application to natural garnet peridotites produced systematic differences
% from the BKN two-pyroxene thermometer. These differences correlate with Na
% in Opx, and Brey & Kohler (1990) were unable to derive a satisfactory Na
% correction from the natural mineral compositions. They therefore stated
% that application of this thermometer to natural rocks must be viewed with
% caution (pp. 1365–1367; Fig. 7b–c).
%
% The selected Opx should be texturally and chemically equilibrated with the
% relevant Cpx generation. Altered grains, exsolved host compositions,
% analyses from different generations, and grains that re-equilibrated at
% different times may yield geologically misleading temperatures. Results
% should preferably be compared with an independent two-pyroxene
% thermometer, especially for Na-rich Opx or natural mantle xenoliths.
%
% For conservative application to natural peridotitic samples, this
% implementation issues non-stopping warnings when:
%   1) input pressure is outside 10–60 kbar,
%   2) a finite calculated temperature is outside 900–1400 degreeC,
%   3) the required Ca input contains NaN, or
%   4) a calculated temperature is NaN or Inf.
%
% The narrower 10–60 kbar pressure warning reflects the natural-system
% validation range. Pressures from 2 to <10 kbar fall within the underlying
% CMS experimental range but outside the natural-system validation range.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following normalized cation variable, preferably on a 6-oxygen basis:
%
%   Opx table variable:
%     Ca_cation_apfu
%
% Important NaN and sign handling:
% - A measured NaN in Ca_cation_apfu is retained as NaN and propagated
%   through the calculation. It is never replaced by zero.
% - Finite Ca_cation_apfu values must be >= 0. Negative finite values stop
%   the calculation with an error.
% - Ca_cation_apfu = 0 is outside the logarithmic domain of the thermometer.
%   It is retained in the output, while TA, T_K, and T_deg are returned as
%   NaN and reported by non-stopping fprintf warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
%   TA = -ln(Ca_in_opx) + 1.843
%   TB = 6425 + 26.4 * P_kbar
%   T(K) = TB / TA
%   T(degreeC) = T(K) - 273.15
%
% where
%   Ca_in_opx = Ca_cation_apfu
%
% Notes:
% - Pressure is supplied and used directly in kbar.
% - Temperature is calculated in Kelvin and converted to degree Celsius.
% - R is retained in the output for compatibility and traceability, although
%   it is not used explicitly in equation (10).
% - Invalid logarithmic-domain values are converted to NaN rather than
%   stopping the calculation or returning a misleading finite temperature.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = BreyKohler1990Ca(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing an Opx table (see above)
%   P_kbar         : pressure in kbar; finite non-negative numeric scalar or
%                    vector
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Opx analysis
%

%% Input validation
% Basic argument checks prevent silent failures caused by missing inputs or
% invalid pressure values.
if nargin < 2
    error('BreyKohler1990Ca requires (rawdata_struct, P_kbar).');
end

if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end

if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

% Convert pressure to a column vector so all downstream table variables have
% a consistent orientation in both fixed-P and range-P modes.
P_kbar = P_kbar(:);

%% 1) Retrieve cation dataset
% Extract the required table from the input struct. The source table is not
% modified.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end

dataset_opx = rawdata_struct.Opx;

% Fail early when the required table variable is missing.
requiredVariables = {'Ca_cation_apfu'};
validateRequiredVariables(dataset_opx, requiredVariables, 'Opx');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each selected Opx analysis produces one table block. Blocks are buffered in
% a preallocated cell array and concatenated only once after the interactive
% loop, avoiding repeated reallocation of the full results table on every
% iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Conservative natural-system validation limits used for warning messages.
calibrationT_min_degC = 900;
calibrationT_max_degC = 1400;
calibrationP_min_kbar = 10;
calibrationP_max_kbar = 60;

% Pressure is common to all Opx analyses selected during this function call,
% so its range warning is printed only once.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–4) Interactive selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Opx) ===');

while true
    % ----- Opx selection -----
    % The first table column is treated as the identifier shown to the user.
    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Orthopyroxene (Opx) data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Opx selected: ' char(string(selectedCode_opx))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    selectedData_opx = dataset_opx(selectedIdx_opx, :);

    % Report NaN without interrupting the calculation. Negative finite values
    % are rejected, while zero is retained and handled as an invalid log
    % domain in calcTemp.
    validateNonNegativeInputs(selectedData_opx);
    nanInputNames = findNaNInputs(selectedData_opx);

    row = calcTemp(selectedData_opx, P_kbar);

    % Store the user-selected identifier for every pressure row.
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row = movevars(row, 'dataCode_opx', 'Before', 1);

    % Store this result as one table block. If the preallocated buffer is
    % full, double its capacity. This occurs only occasionally and avoids
    % reallocating the full results table on every calculation.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_opx)) ': ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_opx)) ': ' ...
            num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside the conservative
    % natural-system validation range of 10–60 kbar. The calculation is not
    % stopped. The message also notes the broader CMS experimental range.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the natural-system validation ' ...
             'range of Brey & Kohler (1990): 10–60 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g kbar. The underlying CMS experimental ' ...
             'calibration extends approximately from 2 to 60 kbar, but ' ...
             '2 to <10 kbar is outside the natural-system validation range.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the
    % approximate natural-system validation range. NaN and Inf are handled
    % separately below.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximate ' ...
             'natural-system validation range of Brey & Kohler (1990): ' ...
             '900–1400 degreeC. %d of %d finite temperature point(s) are ' ...
             'outside the range; calculated finite range = %.4g–%.4g ' ...
             'degreeC for %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_opx)));
    end

    % Print the exact thermometer input names containing NaN. fprintf is
    % used instead of warning so the message remains visible even when
    % MATLAB warning display has been disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         The NaN value(s) were retained and propagated; they were not replaced by zero.\n'], ...
            char(string(selectedCode_opx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % A zero Ca value is not negative, but it lies outside the logarithmic
    % domain. Report it explicitly without stopping the calculation.
    if any(row.Ca_opx == 0)
        fprintf(2, ...
            ['WARNING: Opx.Ca_cation_apfu is zero for %s. ' ...
             'The logarithm is undefined at zero, so TA, T_K, and T_deg ' ...
             'were retained as NaN.\n'], ...
            char(string(selectedCode_opx)));
    end

    % Retain and report all NaN/Inf temperature values, whether caused by a
    % NaN input, zero Ca, or another invalid mathematical domain.
    invalidTemperature = ~isfinite(row.T_deg);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_opx)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    % Ask whether another Opx analysis should be calculated.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'BreyKohler1990Ca', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered blocks once. Return an empty table if the user
% canceled before completing a calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(dataset, requiredVariables, mineralLabel)
% validateRequiredVariables
% Confirm that every required table variable exists before interactive
% selection begins.

missingVariables = strings(numel(requiredVariables), 1);
nMissing = 0;

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    if ~ismember(variableName, dataset.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingVariables(nMissing) = string(variableName);
    end
end

if nMissing > 0
    missingVariables = missingVariables(1:nMissing);
    error('%s table is missing required variable(s): %s', ...
        mineralLabel, char(strjoin(missingVariables, ', ')));
end

end

function nanInputNames = findNaNInputs(data_opx)
% findNaNInputs
% Return the names of required thermometer input variables containing NaN.
% This function never stops the calculation.

thermometerVariables = {'Ca_cation_apfu'};
nanInputNamesBuffer = strings(numel(thermometerVariables), 1);
nFound = 0;

for i = 1:numel(thermometerVariables)
    variableName = thermometerVariables{i};
    variableValue = data_opx.(variableName);

    if any(isnan(variableValue(:)))
        nFound = nFound + 1;
        nanInputNamesBuffer(nFound) = "Opx." + string(variableName);
    end
end

nanInputNames = nanInputNamesBuffer(1:nFound);

end

function validateNonNegativeInputs(data_opx)
% validateNonNegativeInputs
% Stop when a finite thermometer input is negative. Zero and NaN are
% intentionally allowed: zero is handled as an invalid logarithmic domain,
% and NaN is retained and reported by non-stopping warnings.

thermometerVariables = {'Ca_cation_apfu'};
invalidInputNamesBuffer = strings(numel(thermometerVariables), 1);
nInvalid = 0;

for i = 1:numel(thermometerVariables)
    variableName = thermometerVariables{i};
    variableValue = data_opx.(variableName);

    if ~isnumeric(variableValue) || ~isscalar(variableValue) || ...
            ~isreal(variableValue)
        error('Opx variable %s must be a real numeric scalar in a 1-row table.', ...
            variableName);
    end

    if isinf(variableValue)
        error('Opx variable %s must not be Inf.', variableName);
    end

    if isfinite(variableValue) && variableValue < 0
        nInvalid = nInvalid + 1;
        invalidInputNamesBuffer(nInvalid) = "Opx." + string(variableName);
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNamesBuffer(1:nInvalid);
    error(['BreyKohler1990Ca: thermometer input values must be >= 0. ' ...
           'Negative finite value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, P_kbar)
% calcTemp
% Calculate the Ca-in-Opx temperature for one Opx row and one or more
% pressures. One output row is returned for every pressure value.
%
% Inputs:
%   data_opx : 1-row Opx table containing Ca_cation_apfu
%   P_kbar   : pressure in kbar; scalar or vector converted to a column
%
% Output:
%   row : table containing pressure, input Ca, intermediate terms, and
%         calculated temperatures

if height(data_opx) ~= 1
    error('Opx input must be a 1-row table.');
end

P_kbar = P_kbar(:);
nP = numel(P_kbar);

% Physical constant retained for compatibility and traceability. It is not
% used explicitly in equation (10).
R = 8.314462618;

% Preserve the measured scalar Ca value, including NaN, and replicate it for
% every pressure point.
Ca_opx = repmat(data_opx.Ca_cation_apfu, nP, 1);

% Pressure-dependent numerator term.
TB = 6425 + 26.4 .* P_kbar;

% The logarithm is defined only for finite Ca > 0. Initialize TA and the
% temperature outputs as NaN so zero or NaN inputs remain NaN without being
% replaced by zero or causing the calculation to stop.
TA = NaN(nP, 1);
T_K = NaN(nP, 1);
T_deg = NaN(nP, 1);

validLogInput = isfinite(Ca_opx) & Ca_opx > 0;
TA(validLogInput) = -log(Ca_opx(validLogInput)) + 1.843;

% Calculate temperature only where the logarithmic input and denominator are
% valid. A zero TA would otherwise generate Inf and is therefore retained as
% NaN for explicit reporting by the caller.
validTemperatureCalculation = validLogInput & ...
    isfinite(TA) & TA ~= 0 & isfinite(TB);

T_K(validTemperatureCalculation) = ...
    TB(validTemperatureCalculation) ./ TA(validTemperatureCalculation);
T_deg(validTemperatureCalculation) = ...
    T_K(validTemperatureCalculation) - 273.15;

% Build the complete table in one operation with stable variable sizes.
row = table( ...
    P_kbar, ...
    repmat(R, nP, 1), ...
    Ca_opx, ...
    TA, ...
    TB, ...
    T_K, ...
    T_deg, ...
    'VariableNames', { ...
        'P_kbar', ...
        'R', ...
        'Ca_opx', ...
        'TA', ...
        'TB', ...
        'T_K', ...
        'T_deg'});

end
