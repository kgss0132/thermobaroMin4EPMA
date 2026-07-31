function results = Cathelineau1988(rawdata_struct, P_kbar)
% functions/+thermo/+Chlorite/Cathelineau1988.m
% Tested with MATLAB R2024b
%
% Empirical Al(IV)-in-chlorite thermometer
% Cathelineau, M. (1988)
% Clay Minerals, 23, 471–485
% DOI: https://doi.org/10.1180/claymin.1988.023.4.13
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Chlorite analysis at a time and
% calculates temperature using the Cathelineau (1988) Al(IV)-in-chlorite
% geothermometer.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Chlorite analysis and stores all
% results in a single output table.
%
% Both a scalar pressure and a pressure vector are accepted. Therefore, this
% function can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. Pressure is retained in the output for
% traceability, although it is not explicitly used in the thermometer
% equation.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Cathelineau (1988) revised the empirical relationship between tetrahedral
% Al in Chlorite and crystallization temperature using authigenic and
% hydrothermal Chlorite data, principally from the Los Azufres and Salton
% Sea geothermal systems.
%
% Temperature ranges reported in the paper should be distinguished as:
%
%   Abstract calibration range       : approximately 150–300 degreeC
%                                      (abstract, p. 471)
%   Direct regression-data range     : approximately 130–310 degreeC
%                                      (Results, pp. 474 and 476)
%   Broader application discussion   : approximately 100–350 degreeC
%                                      (Table 4 and discussion, pp. 481–482)
%
% The direct regression is therefore treated here as 130–310 degreeC. The
% broader 100–350 degreeC interval is retained as an extended discussion
% range rather than as an equally well constrained calibration range.
%
% Cathelineau (1988) does not provide a quantitative pressure calibration
% range, and pressure does not appear in the thermometer equation. An input
% pressure can therefore be stored for interface compatibility and
% traceability, but pressure-range validity cannot be assessed numerically.
% This implementation prints a non-stopping fprintf warning explaining this
% limitation after the first calculation.
%
% The calibration materials were authigenic minerals formed by water–rock
% interaction in quartz-bearing environments with relatively high
% water/rock ratios. The paper cautions against uncritical application to:
%
%   - closed systems whose Chlorite behaviour is not well constrained,
%   - detrital or inherited Chlorite,
%   - incompletely altered Chlorite formed from earlier sheet silicates,
%   - mechanical or structural mixtures with other clay minerals at scales
%     smaller than approximately 1 micrometre,
%   - unusual Fe-, Zn-, or Cr-rich Chlorite compositions, and
%   - analyses affected by later thermal re-equilibration.
%
% These limitations are discussed in the section "Use and misuse of the
% geothermometric relationships" on pp. 482–483. The calculated temperature
% should be interpreted as the crystallization or later re-equilibration
% temperature recorded by Chlorite, not automatically as the peak
% temperature of the entire rock history.
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) a finite calculated temperature is outside 130–310 degreeC,
%   2) a finite calculated temperature is outside the broader 100–350
%      degreeC discussion range,
%   3) a required or auxiliary calculation input contains NaN,
%   4) a non-finite temperature is calculated, or
%   5) pressure is supplied, because no quantitative pressure calibration
%      range is defined in the paper.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain either:
%   rawdata_struct.Chlorite : table
% or
%   rawdata_struct.Chl      : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% The remaining columns should contain normalized Chlorite cation data on a
% 14-oxygen basis (half-formula basis), or values directly usable as apfu on
% that basis.
%
% Required variables:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu      % treated as total Fe and assumed Fe2+
%   Mg_cation_apfu
%
% Optional variable used in auxiliary compositional calculations:
%   Mn_cation_apfu      % assigned zero only when the column is absent
%
% Optional trace variables stored for traceability:
%   Ti_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%   Cr_cation_apfu
%
% All finite cation values must be greater than or equal to zero. Negative
% values are prohibited. NaN values are retained as missing values, are
% never replaced by zero when the corresponding column exists, propagate
% through outputs that depend on them, and are reported by non-stopping
% fprintf warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Cathelineau (1988) revised the Chlorite Al(IV)-temperature calibration to:
%
%   T(degreeC) = -61.92 + 321.98 * Al_IV
%
% where, on the 14-oxygen basis:
%
%   Al_IV = 4 - Si
%
% This implementation retains the practical compositional constraint used
% in the original script:
%
%   Al_IV = min(Al_total, max(0, 4 - Si))
%
% Auxiliary site quantities are calculated as:
%
%   Al_VI = Al_total - Al_IV
%   Sum_VI = Al_VI + Fe2 + Mg + Mn
%   VAC = 6 - Sum_VI
%
% The paper discusses octahedral vacancy and substitution trends, but the
% thermometer output implemented here is only the revised Al(IV)-based
% relationship.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Cathelineau1988(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing a Chlorite or Chl table
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Chlorite analysis
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('Cathelineau1988 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation dataset
% Extract the Chlorite table from either accepted field name.
disp('=== Step 1: Preparing cation dataset ===');

if isfield(rawdata_struct, 'Chlorite') && istable(rawdata_struct.Chlorite)
    dataset_chl = rawdata_struct.Chlorite;
elseif isfield(rawdata_struct, 'Chl') && istable(rawdata_struct.Chl)
    dataset_chl = rawdata_struct.Chl;
else
    error(['rawdata_struct must contain a Chlorite table as either ' ...
           'rawdata_struct.Chlorite or rawdata_struct.Chl']);
end

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated concatenation of the complete results table inside the loop is
% avoided because it repeatedly reallocates and copies the table.
%
% The cell buffer is preallocated and doubled only when necessary. After the
% interactive loop finishes, all result blocks are concatenated once.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Temperature limits distinguished in Cathelineau (1988).
coreCalibrationT_min_degC = 130;
coreCalibrationT_max_degC = 310;
extendedDiscussionT_min_degC = 100;
extendedDiscussionT_max_degC = 350;

% The paper defines no numerical pressure calibration range. This caution is
% therefore printed once after the first completed calculation.
pressureCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–4) Interactive selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Chlorite) ===');

while true
    % ----- Chlorite selection -----
    % Assumption: the first column stores an identifier displayed to the user.
    dataCodes_chl = dataset_chl{:, 1};

    [selectedIdx_chl, ok] = listdlg( ...
        'PromptString', 'Please select the Chlorite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_chl)), ...
        'ListSize', [320 320]);

    % If the user cancels, exit the loop gracefully.
    if ~ok || isempty(selectedIdx_chl)
        disp('Selection canceled');
        break;
    end

    selectedCode_chl = dataCodes_chl(selectedIdx_chl);
    disp(['Chlorite selected: ' char(string(selectedCode_chl))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    selectedData_chl = dataset_chl(selectedIdx_chl, :);

    % Identify NaN values without replacing them or interrupting calculation.
    nanInputNames = findNaNInputs(selectedData_chl);

    % Negative finite values and infinite values are prohibited. Zero and
    % NaN are allowed; NaN values propagate to outputs that depend on them.
    validateNonNegativeInputs(selectedData_chl);

    row = calcTemp(selectedData_chl, P_kbar);

    % Store the selected identifier for every pressure row.
    row.dataCode_chl = repmat(string(selectedCode_chl), height(row), 1);
    row = movevars(row, {'dataCode_chl'}, 'Before', 1);

    % Store this result as one table block. The buffer grows only when its
    % current capacity is exhausted, rather than on every loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_chl)) ': T_Cathelineau1988 = ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_chl)) ': T_Cathelineau1988 = ' ...
            num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ...
            ' degreeC']);
    end

    % The paper gives no quantitative pressure calibration range, and the
    % equation contains no pressure term. Print this limitation once without
    % stopping calculation.
    if ~pressureCautionIssued
        fprintf(2, ...
            ['WARNING: Cathelineau (1988) does not define a quantitative ' ...
             'pressure calibration range, and pressure is not used in the ' ...
             'thermometer equation. The supplied pressure is stored only for ' ...
             'traceability; input range = %.4g–%.4g kbar. Pressure-range ' ...
             'validity cannot be evaluated from this calibration.\n'], ...
            min(P_kbar), ...
            max(P_kbar));
        pressureCautionIssued = true;
    end

    % Warn first when any finite temperature lies outside the broader
    % 100–350 degreeC discussion range. This is a stronger extrapolation than
    % merely falling outside the direct 130–310 degreeC regression range.
    finiteTemperature = isfinite(row.T_deg);
    outsideExtendedDiscussion = finiteTemperature & ...
        (row.T_deg < extendedDiscussionT_min_degC | ...
         row.T_deg > extendedDiscussionT_max_degC);

    outsideCoreCalibration = finiteTemperature & ...
        (row.T_deg < coreCalibrationT_min_degC | ...
         row.T_deg > coreCalibrationT_max_degC);

    if any(outsideExtendedDiscussion)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside both the direct ' ...
             'regression range (approximately 130–310 degreeC) and the broader ' ...
             '100–350 degreeC application discussion of Cathelineau (1988). ' ...
             '%d of %d finite temperature point(s) are outside 100–350 degreeC; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s.\n'], ...
            sum(outsideExtendedDiscussion), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    elseif any(outsideCoreCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the direct regression ' ...
             'range of Cathelineau (1988): approximately 130–310 degreeC, ' ...
             'although it remains within the broader 100–350 degreeC range ' ...
             'discussed in the paper. %d of %d finite temperature point(s) ' ...
             'are outside 130–310 degreeC; calculated finite range = ' ...
             '%.4g–%.4g degreeC for %s.\n'], ...
            sum(outsideCoreCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    end

    % Print a non-stopping warning when an input contains NaN. Existing NaN
    % values are retained and are never converted to zero.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Cathelineau (1988) input(s) for ' ...
             '%s: %s.\n' ...
             '         NaN values were retained and were not replaced by zero. ' ...
             'Outputs that depend on these values may therefore be NaN.\n'], ...
            char(string(selectedCode_chl)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain a result-based check for NaN/Inf caused by missing inputs or
    % mathematically undefined intermediate values.
    invalidTemperature = ~isfinite(row.T_deg);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_chl)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Cathelineau1988', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after all selections have
% finished. Return an empty table when no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_chl)
% findNaNInputs
% Return the names of required or auxiliary calculation variables that
% contain NaN. This function does not convert NaN to zero and does not throw
% an error for missing values.

variableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu'};

displayNames = "Chlorite." + string(variableNames(:));
nanMask = false(numel(variableNames), 1);

for i = 1:numel(variableNames)
    variableName = variableNames{i};

    if ismember(variableName, data_chl.Properties.VariableNames)
        variableValue = data_chl.(variableName);
        nanMask(i) = any(isnan(variableValue(:)));
    end
end

nanInputNames = displayNames(nanMask);

end

function validateNonNegativeInputs(data_chl)
% validateNonNegativeInputs
% Stop the calculation when a finite cation value is negative or when an
% infinite value is present. Zero is allowed. NaN is intentionally allowed
% so that it remains missing and is reported by fprintf.

variableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu', ...
    'Cr_cation_apfu'};

displayNames = "Chlorite." + string(variableNames(:));
negativeMask = false(numel(variableNames), 1);
infiniteMask = false(numel(variableNames), 1);
nonNumericMask = false(numel(variableNames), 1);

for i = 1:numel(variableNames)
    variableName = variableNames{i};

    if ismember(variableName, data_chl.Properties.VariableNames)
        variableValue = data_chl.(variableName);

        if ~isnumeric(variableValue)
            nonNumericMask(i) = true;
            continue;
        end

        negativeMask(i) = any(isfinite(variableValue(:)) & variableValue(:) < 0);
        infiniteMask(i) = any(isinf(variableValue(:)));
    end
end

if any(nonNumericMask)
    error(['Cathelineau1988: cation variables must be numeric. ' ...
           'Non-numeric variable(s): ' ...
           char(strjoin(displayNames(nonNumericMask), ', ')) '.']);
end

if any(negativeMask)
    error(['Cathelineau1988: cation values must be greater than or ' ...
           'equal to zero. Negative value(s) were found in: ' ...
           char(strjoin(displayNames(negativeMask), ', ')) '.']);
end

if any(infiniteMask)
    error(['Cathelineau1988: infinite cation value(s) are not allowed in: ' ...
           char(strjoin(displayNames(infiniteMask), ', ')) '.']);
end

end

function row = calcTemp(data_chl, P_kbar)
% calcTemp
% Compute a temperature estimate for one Chlorite analysis and return one
% output row for each supplied pressure value.
%
% Pressure is not explicitly used in the Cathelineau (1988) equation.
% Therefore, composition-derived and temperature outputs are repeated for
% every pressure value, while P_kbar records the pressure supplied by the
% fixed-pressure or pressure-range workflow.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

% --- Prepare Chlorite cation row ---
chl = prepareChloriteRow(data_chl);

% --- Approximate site allocation on a 14-oxygen basis ---
site = calcChloriteSites(chl);

% --- Cathelineau (1988) Al(IV)-in-Chlorite thermometer ---
T_deg_scalar = -61.92 + 321.98 .* site.Al_IV;
T_K_scalar = T_deg_scalar + 273.15;

% --- Auxiliary compositional indices ---
% NaN and zero denominators are deliberately allowed to propagate naturally.
Mg_number_scalar = chl.Mg ./ (chl.Mg + chl.Fe2);
Fe_ratio_scalar = chl.Fe2 ./ (chl.Fe2 + chl.Mg + chl.Mn);

% --- Calibration / sanity flags ---
is_14O_tetra_reasonable_scalar = ...
    isfinite(chl.Si) && chl.Si >= 2.0 && chl.Si <= 4.0;

is_AlIV_reasonable_scalar = ...
    isfinite(site.Al_IV) && site.Al_IV >= 0.0 && site.Al_IV <= 2.0;

is_VAC_reasonable_scalar = ...
    isfinite(site.VAC) && site.VAC >= -0.25 && site.VAC <= 1.50;

is_in_core_calibration_T_range_scalar = ...
    isfinite(T_deg_scalar) && T_deg_scalar >= 130 && T_deg_scalar <= 310;

is_in_extended_discussion_T_range_scalar = ...
    isfinite(T_deg_scalar) && T_deg_scalar >= 100 && T_deg_scalar <= 350;

AlIV_at_130C = (130 + 61.92) ./ 321.98;
AlIV_at_310C = (310 + 61.92) ./ 321.98;

is_in_core_calibration_AlIV_range_scalar = ...
    isfinite(site.Al_IV) && ...
    site.Al_IV >= AlIV_at_130C && site.Al_IV <= AlIV_at_310C;

% A quantitative pressure calibration range is not defined in the paper.
pressure_calibration_defined = false(nP, 1);

% The conservative recommendation flag includes the direct temperature
% regression range as well as basic compositional sanity checks.
recommended_by_Cathelineau1988 = ...
    repmat(is_14O_tetra_reasonable_scalar, nP, 1) & ...
    repmat(is_AlIV_reasonable_scalar, nP, 1) & ...
    repmat(is_VAC_reasonable_scalar, nP, 1) & ...
    repmat(is_in_core_calibration_T_range_scalar, nP, 1);

% --- Pack outputs ---
% Every scalar composition-derived result is repeated to match the number of
% pressure values and to keep all table variables at a stable height.
row = table();

row.P_kbar = P_kbar;

row.Si = repmat(chl.Si, nP, 1);
row.Al = repmat(chl.Al, nP, 1);
row.Fe2 = repmat(chl.Fe2, nP, 1);
row.Mg = repmat(chl.Mg, nP, 1);
row.Mn = repmat(chl.Mn, nP, 1);
row.Ti = repmat(chl.Ti, nP, 1);
row.Ca = repmat(chl.Ca, nP, 1);
row.Na = repmat(chl.Na, nP, 1);
row.K = repmat(chl.K, nP, 1);
row.Cr = repmat(chl.Cr, nP, 1);

row.Al_IV = repmat(site.Al_IV, nP, 1);
row.Al_VI = repmat(site.Al_VI, nP, 1);
row.Sum_VI = repmat(site.Sum_VI, nP, 1);
row.VAC = repmat(site.VAC, nP, 1);

row.Mg_number = repmat(Mg_number_scalar, nP, 1);
row.Fe_ratio_Fe_over_FeMgMn = repmat(Fe_ratio_scalar, nP, 1);

row.T_deg = repmat(T_deg_scalar, nP, 1);
row.T_K = repmat(T_K_scalar, nP, 1);

row.is_14O_tetra_reasonable = ...
    repmat(is_14O_tetra_reasonable_scalar, nP, 1);
row.is_AlIV_reasonable = ...
    repmat(is_AlIV_reasonable_scalar, nP, 1);
row.is_VAC_reasonable = ...
    repmat(is_VAC_reasonable_scalar, nP, 1);
row.is_in_core_calibration_T_range = ...
    repmat(is_in_core_calibration_T_range_scalar, nP, 1);
row.is_in_extended_discussion_T_range = ...
    repmat(is_in_extended_discussion_T_range_scalar, nP, 1);
row.is_in_core_calibration_AlIV_range = ...
    repmat(is_in_core_calibration_AlIV_range_scalar, nP, 1);
row.pressure_calibration_defined = pressure_calibration_defined;
row.recommended_by_Cathelineau1988 = recommended_by_Cathelineau1988;

end

function chl = prepareChloriteRow(data_chl)
% prepareChloriteRow
% Extract one Chlorite analysis from a 1-row table.
%
% Required values may be finite or NaN. Optional variables are assigned zero
% only when the corresponding column is absent. If an optional column exists
% and contains NaN, that NaN is retained unchanged.

if height(data_chl) ~= 1
    error('Chlorite input must be a 1-row table.');
end

chl = struct();

chl.Si  = getVarOrError(data_chl, 'Si_cation_apfu', 'Chlorite');
chl.Al  = getVarOrError(data_chl, 'Al_cation_apfu', 'Chlorite');
chl.FeT = getVarOrError(data_chl, 'Fe_cation_apfu', 'Chlorite');
chl.Mg  = getVarOrError(data_chl, 'Mg_cation_apfu', 'Chlorite');
chl.Mn  = getVarOrZeroIfMissing(data_chl, 'Mn_cation_apfu');

% Optional trace constituents are retained for traceability.
chl.Ti  = getVarOrZeroIfMissing(data_chl, 'Ti_cation_apfu');
chl.Ca  = getVarOrZeroIfMissing(data_chl, 'Ca_cation_apfu');
chl.Na  = getVarOrZeroIfMissing(data_chl, 'Na_cation_apfu');
chl.K   = getVarOrZeroIfMissing(data_chl, 'K_cation_apfu');
chl.Cr  = getVarOrZeroIfMissing(data_chl, 'Cr_cation_apfu');

% Following the empirical Chlorite treatment, all Fe is treated as Fe2+.
chl.Fe2 = chl.FeT;

end

function site = calcChloriteSites(chl)
% calcChloriteSites
% Calculate the approximate Chlorite site quantities used by the original
% script on a 14-oxygen basis.
%
%   Al_IV = min(Al_total, max(0, 4 - Si))
%   Al_VI = Al_total - Al_IV
%   Sum_VI = Al_VI + Fe2 + Mg + Mn
%   VAC = 6 - Sum_VI
%
% NaN inputs remain NaN. Consistency errors are raised only for finite,
% mathematically invalid derived values.

site = struct();

% Explicitly preserve NaN in the thermometer variables. This avoids relying
% on version-dependent min/max handling of missing values.
if isnan(chl.Si) || isnan(chl.Al)
    site.Al_IV = NaN;
    site.Al_VI = NaN;
else
    site.Al_IV = min(chl.Al, max(0, 4 - chl.Si));
    site.Al_VI = chl.Al - site.Al_IV;
end

if isfinite(site.Al_VI) && site.Al_VI < -1e-10
    error('Negative octahedral Al calculated. Check Chlorite normalization.');
end

if isfinite(site.Al_VI)
    site.Al_VI = max(0, site.Al_VI);
end

site.Sum_VI = site.Al_VI + chl.Fe2 + chl.Mg + chl.Mn;
site.VAC = 6 - site.Sum_VI;

end

function value = getVarOrError(tbl, varName, mineralLabel)
% getVarOrError
% Retrieve a required numeric scalar. NaN is allowed and retained. Infinite
% and negative values are prohibited.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s must be a numeric scalar in a 1-row table.', varName);
end

if isinf(value)
    error('Variable %s contains an infinite value.', varName);
end

if isfinite(value) && value < 0
    error('Variable %s contains a negative value.', varName);
end

end

function value = getVarOrZeroIfMissing(tbl, varName)
% getVarOrZeroIfMissing
% Retrieve an optional numeric scalar. Assign zero only when the column is
% absent. If the column exists and contains NaN, retain NaN unchanged.

if ismember(varName, tbl.Properties.VariableNames)
    value = tbl.(varName);

    if ~isnumeric(value) || ~isscalar(value)
        error('Variable %s must be a numeric scalar in a 1-row table.', varName);
    end

    if isinf(value)
        error('Variable %s contains an infinite value.', varName);
    end

    if isfinite(value) && value < 0
        error('Variable %s contains a negative value.', varName);
    end
else
    value = 0;
end

end
