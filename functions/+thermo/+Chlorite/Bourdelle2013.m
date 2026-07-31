function results = Bourdelle2013(rawdata_struct, P_kbar)
% functions/+thermo/+Chlorite/Bourdelle2013.m
% Tested with MATLAB R2024b
%
% Semi-empirical single-chlorite thermometer
% Bourdelle, F., Parra, T., Chopin, C., and Beyssac, O. (2013)
% Contributions to Mineralogy and Petrology, 165, 723–735
% DOI: https://doi.org/10.1007/s00410-012-0832-7
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Chlorite analysis at a time and
% calculates temperature using the Bourdelle et al. (2013) chlorite
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
% Bourdelle et al. (2013) calibrated the thermometer using 161 published
% Chlorite analyses for which independently constrained formation
% temperatures and pressures were available. The selected calibration data
% cover:
%
%   Temperature : 50–350 degreeC
%   Pressure    : less than 4 kbar
%   Material    : authigenic, hydrothermal, diagenetic, or low-grade
%                 metamorphic Chlorite; analyses identified as detrital
%                 were excluded
%   Assemblage  : quartz-bearing samples only
%
% The calibration-data selection criteria are described on pp. 724–725.
% The linear calibration equation is presented on p. 729. The limitations
% and recommended application range are reiterated on p. 734.
%
% The implementation assumes:
%
%   a_SiO2 = 1       (quartz-bearing system)
%   a_H2O  = 1
%   Fe_total = Fe2+
%   pressure effects are negligible within P < 4 kbar
%
% The assumption a_H2O = 1 is regarded as reasonable for many diagenetic to
% low-grade metamorphic systems, but the authors note that it may be
% questionable locally (p. 729). Application to quartz-free samples,
% water-undersaturated systems, strongly saline or non-aqueous fluids,
% detrital/inherited Chlorite, mixed-layer material, or contaminated
% analyses should therefore be treated cautiously.
%
% To reduce contamination and mixed-layer effects, the calibration dataset
% was screened using K2O + Na2O + CaO < 1 wt% (pp. 724–725 and 730–731).
% The present function receives normalized cation data and cannot apply this
% oxide-based screening automatically; the user should inspect the original
% analyses before calculation.
%
% The authors report that scatter is commonly lower at low temperature
% (approximately +/-30 degreeC) and larger above 300 degreeC
% (approximately +/-50–60 degreeC), although many test analyses agreed with
% independent temperatures within about 20 degreeC (pp. 729–730 and 734).
%
% Use outside the calibration range, especially at T > 350 degreeC or
% P >= 4 kbar, is not recommended because assumptions concerning pressure,
% reaction volume, heat capacity, and water activity may no longer be valid
% (p. 734).
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) any input pressure is greater than or equal to 4 kbar,
%   2) any finite calculated temperature is outside 50–350 degreeC,
%   3) a required calculation input contains NaN, or
%   4) a non-finite temperature is calculated.
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
% Required variables used in the thermometer:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu      % treated as total Fe and assumed Fe2+
%   Mg_cation_apfu
%
% Optional variable used in the site-allocation calculation:
%   Mn_cation_apfu      % assumed zero only when the column is absent
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
% through the calculation, and are reported by non-stopping warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Bourdelle et al. (2013) use a semi-ordered Chlorite site-allocation model
% and the following reaction:
%
%   Mg-ChlS + 3 Mg-Sud = 3 Mg-Am + 7 Qtz + 4 H2O
%
% The equilibrium constant is:
%
%   logK = 3*log10(a_MgAm) - log10(a_MgChlS)
%          - 3*log10(a_MgSud)
%
% with ideal end-member activities:
%
%   a_MgChlS = (X_Si_T2)^2 * (X_Mg_M23)^4 * (X_Mg_M14)^2
%
%   a_MgAm   = (X_Al_T2)^2 * (X_Mg_M23)^4 * (X_Al_M14)^2
%
%   a_MgSud  = 256 * X_Si_T2 * X_Al_T2 * X_Al_M14 * X_h_M14 ...
%                    * (X_Mg_M23)^2 * (X_Al_M23)^2
%
% The principal thermometer is Eq. (7) on p. 729:
%
%   logK = 9400 / T(K) + 23.40
%
% rearranged as:
%
%   T(K) = 9400 / (logK - 23.40)
%
%   T(degreeC) = 9400 / (logK - 23.40) - 273.15
%
% Bourdelle et al. (2013) also provide a quadratic regression, Eq. (8), for
% improved fitting near 300–350 degreeC. They explicitly state that it has
% no physical basis and should not be used outside 150–350 degreeC (p. 729).
% It is retained here only as an auxiliary output.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Bourdelle2013(rawdata_struct, P_kbar)
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
    error('Bourdelle2013 requires (rawdata_struct, P_kbar).');
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
% Repeated concatenation of the full results table inside the loop is
% avoided because it repeatedly reallocates and copies the table.
%
% The cell buffer is preallocated and doubled only when necessary. After the
% interactive loop finishes, all result blocks are concatenated once.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Calibration limits reported by Bourdelle et al. (2013).
calibrationT_min_degC = 50;
calibrationT_max_degC = 350;
calibrationP_max_kbar = 4;

% Pressure is common to all selected analyses in this function call.
% Therefore, the pressure warning is printed only once.
pressureOutsideCalibration = P_kbar >= calibrationP_max_kbar;
pressureWarningIssued = false;

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

    % Check only variables that are used in the thermometer calculation.
    % NaN is intentionally allowed and retained.
    nanInputNames = findNaNInputs(selectedData_chl);

    % Negative finite values and infinite values are prohibited. Zero and
    % NaN are allowed so that mathematically undefined cases can propagate
    % to NaN and be reported without stopping the overall workflow.
    validateNonNegativeInputs(selectedData_chl);

    row = calcTemp(selectedData_chl, P_kbar);

    % Store the selected identifier for every pressure row.
    row.dataCode_chl = repmat(string(selectedCode_chl), height(row), 1);
    row = movevars(row, {'dataCode_chl'}, 'Before', 1);

    % Store this result as one table block.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_chl)) ': T_Bourdelle2013 = ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_chl)) ': T_Bourdelle2013 = ' ...
            num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ...
            ' degreeC']);
    end

    % Warn once when any input pressure is outside the calibration range.
    % The calculation is not stopped.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the calibration range of ' ...
             'Bourdelle et al. (2013): P < 4 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside 50–350 degreeC.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the calibration ' ...
             'range of Bourdelle et al. (2013): 50–350 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    end

    % Print a non-stopping warning immediately after the result when a
    % calculation input contains NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         NaN values were retained and were not replaced by zero. ' ...
             'The calculated temperature may therefore be NaN.\n'], ...
            char(string(selectedCode_chl)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain a result-based check for NaN/Inf caused by either missing inputs
    % or mathematically undefined intermediate values.
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
        'Bourdelle2013', ...
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
% Return the names of thermometer input variables that contain NaN.
% This function never converts NaN to zero and does not throw an error.

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
% Stop the calculation when a finite cation value used by the thermometer
% is negative, or when an infinite value is present. Zero is allowed. NaN is
% intentionally allowed so that it propagates and is reported by fprintf.

variableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu'};

displayNames = "Chlorite." + string(variableNames(:));
negativeMask = false(numel(variableNames), 1);
infiniteMask = false(numel(variableNames), 1);

for i = 1:numel(variableNames)
    variableName = variableNames{i};

    if ismember(variableName, data_chl.Properties.VariableNames)
        variableValue = data_chl.(variableName);

        if ~isnumeric(variableValue)
            error('Variable %s must be numeric.', variableName);
        end

        negativeMask(i) = any(isfinite(variableValue(:)) & variableValue(:) < 0);
        infiniteMask(i) = any(isinf(variableValue(:)));
    end
end

if any(negativeMask)
    error(['Bourdelle2013: cation values used in the thermometer must ' ...
           'be greater than or equal to zero. Negative value(s) were found in: ' ...
           char(strjoin(displayNames(negativeMask), ', ')) '.']);
end

if any(infiniteMask)
    error(['Bourdelle2013: infinite cation value(s) are not allowed in: ' ...
           char(strjoin(displayNames(infiniteMask), ', ')) '.']);
end

end

function row = calcTemp(data_chl, P_kbar)
% calcTemp
% Compute a temperature estimate for one Chlorite analysis and return one
% output row for each supplied pressure value.
%
% Pressure is not explicitly used in Bourdelle et al. (2013) Eq. (7).
% Therefore, compositional and temperature outputs are repeated for each
% pressure value, while P_kbar records the pressure used by the calling
% fixed-pressure or pressure-range workflow.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

% --- Prepare Chlorite cation row ---
chl = prepareChloriteRow(data_chl);

% --- Site allocation on a 14-oxygen basis ---
site = calcChloriteSites(chl);

% --- Ideal activities for Bourdelle et al. (2013) ---
act = calcIdealActivities(site);

% --- Bourdelle et al. (2013) logK and linear thermometer Eq. (7) ---
logK = 3 .* log10(act.a_MgAm) ...
     - 1 .* log10(act.a_MgChlS) ...
     - 3 .* log10(act.a_MgSud);

% Preserve NaN/non-finite values. A non-finite or non-positive denominator
% is represented as NaN rather than stopping the calculation.
if isfinite(logK) && logK > 23.40
    T_K_scalar = 9400 ./ (logK - 23.40);
    T_deg_scalar = T_K_scalar - 273.15;
else
    T_K_scalar = NaN;
    T_deg_scalar = NaN;
end

% --- Auxiliary quadratic fit, Eq. (8) ---
T_quad_K_scalar = calcQuadraticTK(logK);
T_quad_deg_scalar = T_quad_K_scalar - 273.15;

% --- Simple compositional indices ---
Mg_number_scalar = chl.Mg ./ (chl.Mg + chl.Fe2);
Fe_ratio_scalar = chl.Fe2 ./ (chl.Fe2 + chl.Mg + chl.Mn);

% --- Calibration / sanity flags ---
is_14O_tetra_reasonable_scalar = ...
    isfinite(chl.Si) && chl.Si >= 2.0 && chl.Si <= 4.0;

is_AlIV_reasonable_scalar = ...
    isfinite(site.Al_IV) && site.Al_IV >= 0.0 && site.Al_IV <= 2.0;

is_VAC_reasonable_scalar = ...
    isfinite(site.VAC) && site.VAC >= -0.10 && site.VAC <= 1.50;

is_positive_sitefractions_for_logK_scalar = ...
    isfinite(site.X_Si_T2)  && site.X_Si_T2  > 0 && ...
    isfinite(site.X_Al_T2)  && site.X_Al_T2  > 0 && ...
    isfinite(site.X_Mg_M23) && site.X_Mg_M23 > 0 && ...
    isfinite(site.X_Mg_M14) && site.X_Mg_M14 > 0 && ...
    isfinite(site.X_Al_M23) && site.X_Al_M23 > 0 && ...
    isfinite(site.X_Al_M14) && site.X_Al_M14 > 0 && ...
    isfinite(site.X_h_M14)  && site.X_h_M14  > 0;

is_logK_valid_scalar = isfinite(logK) && logK > 23.40;

is_in_recommended_T_range_scalar = ...
    isfinite(T_deg_scalar) && ...
    T_deg_scalar >= 50 && T_deg_scalar <= 350;

is_in_recommended_P_range = P_kbar < 4;

% Include both the temperature and pressure calibration checks.
recommended_by_Bourdelle2013 = ...
    repmat(is_14O_tetra_reasonable_scalar, nP, 1) & ...
    repmat(is_AlIV_reasonable_scalar, nP, 1) & ...
    repmat(is_VAC_reasonable_scalar, nP, 1) & ...
    repmat(is_positive_sitefractions_for_logK_scalar, nP, 1) & ...
    repmat(is_logK_valid_scalar, nP, 1) & ...
    repmat(is_in_recommended_T_range_scalar, nP, 1) & ...
    is_in_recommended_P_range;

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
row.VAC = repmat(site.VAC, nP, 1);
row.Sum_VI = repmat(site.Sum_VI, nP, 1);

row.Si_T2 = repmat(site.Si_T2, nP, 1);
row.Al_T2 = repmat(site.Al_T2, nP, 1);
row.Al_M14 = repmat(site.Al_M14, nP, 1);
row.Al_M23 = repmat(site.Al_M23, nP, 1);
row.Mg_M14 = repmat(site.Mg_M14, nP, 1);
row.Mg_M23 = repmat(site.Mg_M23, nP, 1);
row.FeLike_M14 = repmat(site.FeLike_M14, nP, 1);
row.FeLike_M23 = repmat(site.FeLike_M23, nP, 1);
row.h_M14 = repmat(site.h_M14, nP, 1);

row.X_Si_T2 = repmat(site.X_Si_T2, nP, 1);
row.X_Al_T2 = repmat(site.X_Al_T2, nP, 1);
row.X_Al_M14 = repmat(site.X_Al_M14, nP, 1);
row.X_Al_M23 = repmat(site.X_Al_M23, nP, 1);
row.X_Mg_M14 = repmat(site.X_Mg_M14, nP, 1);
row.X_Mg_M23 = repmat(site.X_Mg_M23, nP, 1);
row.X_FeLike_M14 = repmat(site.X_FeLike_M14, nP, 1);
row.X_FeLike_M23 = repmat(site.X_FeLike_M23, nP, 1);
row.X_h_M14 = repmat(site.X_h_M14, nP, 1);

row.a_MgChlS = repmat(act.a_MgChlS, nP, 1);
row.a_FeChlS = repmat(act.a_FeChlS, nP, 1);
row.a_MgAm = repmat(act.a_MgAm, nP, 1);
row.a_FeAm = repmat(act.a_FeAm, nP, 1);
row.a_MgSud = repmat(act.a_MgSud, nP, 1);
row.a_FeSud = repmat(act.a_FeSud, nP, 1);

row.logK = repmat(logK, nP, 1);

row.Mg_number = repmat(Mg_number_scalar, nP, 1);
row.Fe_ratio_Fe_over_FeMgMn = repmat(Fe_ratio_scalar, nP, 1);

row.T_deg = repmat(T_deg_scalar, nP, 1);
row.T_K = repmat(T_K_scalar, nP, 1);
row.T_quad_deg = repmat(T_quad_deg_scalar, nP, 1);
row.T_quad_K = repmat(T_quad_K_scalar, nP, 1);
row.deltaT_linear_minus_quadratic = ...
    repmat(T_deg_scalar - T_quad_deg_scalar, nP, 1);

row.is_14O_tetra_reasonable = ...
    repmat(is_14O_tetra_reasonable_scalar, nP, 1);
row.is_AlIV_reasonable = ...
    repmat(is_AlIV_reasonable_scalar, nP, 1);
row.is_VAC_reasonable = ...
    repmat(is_VAC_reasonable_scalar, nP, 1);
row.is_positive_sitefractions_for_logK = ...
    repmat(is_positive_sitefractions_for_logK_scalar, nP, 1);
row.is_logK_valid = repmat(is_logK_valid_scalar, nP, 1);
row.is_in_recommended_P_range = is_in_recommended_P_range;
row.is_in_recommended_T_range = ...
    repmat(is_in_recommended_T_range_scalar, nP, 1);
row.recommended_by_Bourdelle2013 = recommended_by_Bourdelle2013;

end

function chl = prepareChloriteRow(data_chl)
% prepareChloriteRow
% Extract one Chlorite analysis from a 1-row table.
%
% Required values may be finite or NaN. Optional variables are assigned zero
% only when the corresponding column is absent. If an optional column exists
% and contains NaN, that NaN is retained.

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

% Following Bourdelle et al. (2013), all Fe is treated as Fe2+.
chl.Fe2 = chl.FeT;

end

function site = calcChloriteSites(chl)
% calcChloriteSites
% Perform the semi-ordered Chlorite site allocation on a 14-oxygen basis
% following the practical framework described by Bourdelle et al. (2013).
%
% Site groups:
%   T2     : 2 tetrahedral positions
%   M2+M3  : 4 octahedral positions
%   M1+M4  : 2 octahedral positions
%
% Assumptions:
% - Al_IV is restricted to T2
% - vacancies are restricted to M1+M4
% - Al_VI fills M1+M4 first, then excess Al_VI goes to M2+M3
% - Mg-Fe-Mn fills M2+M3 first, then the remainder goes to M1+M4
% - all Fe is treated as Fe2+
%
% If a calculation input is NaN, derived site occupancies remain NaN.
% Consistency errors are raised only when the relevant derived values are
% finite and mathematically invalid.

calculationInputs = [chl.Si, chl.Al, chl.Fe2, chl.Mg, chl.Mn];

% Guarantee explicit NaN propagation without relying on min/max missing-data
% behavior. If any value used in the site calculation is NaN, return NaN for
% every derived site quantity and continue to the non-stopping warnings.
if any(isnan(calculationInputs))
    site = makeNaNSiteStruct();
    return;
end

site = struct();

% --- Tetrahedral allocation ---
site.Al_IV = min(chl.Al, max(0, 4 - chl.Si));
site.Al_VI = chl.Al - site.Al_IV;

if isfinite(site.Al_VI) && site.Al_VI < -1e-10
    error('Negative octahedral Al calculated. Check Chlorite normalization.');
end

if isfinite(site.Al_VI)
    site.Al_VI = max(0, site.Al_VI);
end

% --- Raw octahedral sum and vacancy ---
site.Sum_VI = site.Al_VI + chl.Fe2 + chl.Mg + chl.Mn;
site.VAC = 6 - site.Sum_VI;

% Small finite negative values from rounding are clipped for site-fraction use.
VAC_eff = site.VAC;

if isfinite(VAC_eff) && VAC_eff < 0 && VAC_eff > -1e-6
    VAC_eff = 0;
end

if isfinite(VAC_eff) && VAC_eff < 0
    error(['Negative octahedral vacancy calculated. ' ...
           'Check Chlorite cation normalization on the 14-oxygen basis.']);
end

if isfinite(VAC_eff) && VAC_eff > 2
    error('Octahedral vacancy exceeds the M1+M4 site capacity.');
end

% --- T2 site group ---
site.Si_T2 = chl.Si - 2;
site.Al_T2 = site.Al_IV;

if (isfinite(site.Si_T2) && site.Si_T2 < -1e-8) || ...
        (isfinite(site.Al_T2) && site.Al_T2 < -1e-8)
    error('Negative T2 site occupancy calculated.');
end

sumT2_occupancy = site.Si_T2 + site.Al_T2;

if isfinite(sumT2_occupancy) && abs(sumT2_occupancy - 2) > 1e-6
    error('T2 site occupancy does not sum to 2. Check Chlorite normalization.');
end

% --- Octahedral divalent pool ---
Div_total = chl.Mg + chl.Fe2 + chl.Mn;
Mg_ratio = chl.Mg ./ Div_total;
FeLike_ratio = (chl.Fe2 + chl.Mn) ./ Div_total;

% --- M1+M4 and M2+M3 capacities ---
M14_capacity = 2 - VAC_eff;
M23_capacity = 4;

if (isfinite(M14_capacity) && M14_capacity < -1e-10) || ...
        (isfinite(M23_capacity) && M23_capacity < -1e-10)
    error('Invalid octahedral site capacities.');
end

% --- Al_VI fills M1+M4 first, then M2+M3 ---
site.Al_M14 = min(site.Al_VI, M14_capacity);
site.Al_M23 = site.Al_VI - site.Al_M14;

if isfinite(site.Al_M23) && site.Al_M23 < -1e-10
    error('Negative Al_M23 calculated.');
end

if isfinite(site.Al_M23)
    site.Al_M23 = max(0, site.Al_M23);
end

if isfinite(site.Al_M23) && site.Al_M23 > M23_capacity + 1e-10
    error('Al_M23 exceeds M2+M3 site capacity.');
end

% --- Fe-Mg-Mn fills M2+M3 first, then M1+M4 ---
Div_M23_capacity = M23_capacity - site.Al_M23;
Div_M14_capacity = M14_capacity - site.Al_M14;

if (isfinite(Div_M23_capacity) && Div_M23_capacity < -1e-10) || ...
        (isfinite(Div_M14_capacity) && Div_M14_capacity < -1e-10)
    error('Negative divalent-cation site capacity calculated.');
end

site.Div_M23 = min(Div_total, Div_M23_capacity);
site.Div_M14 = Div_total - site.Div_M23;

if isfinite(site.Div_M14) && isfinite(Div_M14_capacity) && ...
        site.Div_M14 > Div_M14_capacity + 1e-8
    error(['Calculated divalent cations exceed available M1+M4 capacity. ' ...
           'Check Chlorite cation normalization.']);
end

if isfinite(site.Div_M14)
    site.Div_M14 = max(0, site.Div_M14);
end

% --- Split the divalent pool using a common Mg/(Fe+Mg+Mn) ratio ---
site.Mg_M23 = site.Div_M23 .* Mg_ratio;
site.Mg_M14 = site.Div_M14 .* Mg_ratio;

site.FeLike_M23 = site.Div_M23 .* FeLike_ratio;
site.FeLike_M14 = site.Div_M14 .* FeLike_ratio;

site.h_M14 = VAC_eff;

% --- Convert occupancies to site fractions ---
site.X_Si_T2 = site.Si_T2 ./ 2;
site.X_Al_T2 = site.Al_T2 ./ 2;

site.X_Al_M23 = site.Al_M23 ./ 4;
site.X_Mg_M23 = site.Mg_M23 ./ 4;
site.X_FeLike_M23 = site.FeLike_M23 ./ 4;

site.X_Al_M14 = site.Al_M14 ./ 2;
site.X_Mg_M14 = site.Mg_M14 ./ 2;
site.X_FeLike_M14 = site.FeLike_M14 ./ 2;
site.X_h_M14 = site.h_M14 ./ 2;

% --- Final consistency checks ---
sumT2 = site.X_Si_T2 + site.X_Al_T2;
sumM23 = site.X_Al_M23 + site.X_Mg_M23 + site.X_FeLike_M23;
sumM14 = site.X_Al_M14 + site.X_Mg_M14 + ...
    site.X_FeLike_M14 + site.X_h_M14;

if isfinite(sumT2) && abs(sumT2 - 1) > 1e-6
    error('T2 site fractions do not sum to 1.');
end
if isfinite(sumM23) && abs(sumM23 - 1) > 1e-6
    error('M2+M3 site fractions do not sum to 1.');
end
if isfinite(sumM14) && abs(sumM14 - 1) > 1e-6
    error('M1+M4 site fractions do not sum to 1.');
end

end

function site = makeNaNSiteStruct()
% makeNaNSiteStruct
% Create a Chlorite site-allocation structure populated entirely with NaN.
% This guarantees that missing calculation inputs remain missing throughout
% all intermediate and final outputs.

fieldNames = { ...
    'Al_IV'; ...
    'Al_VI'; ...
    'Sum_VI'; ...
    'VAC'; ...
    'Si_T2'; ...
    'Al_T2'; ...
    'Al_M14'; ...
    'Al_M23'; ...
    'Div_M23'; ...
    'Div_M14'; ...
    'Mg_M23'; ...
    'Mg_M14'; ...
    'FeLike_M23'; ...
    'FeLike_M14'; ...
    'h_M14'; ...
    'X_Si_T2'; ...
    'X_Al_T2'; ...
    'X_Al_M23'; ...
    'X_Mg_M23'; ...
    'X_FeLike_M23'; ...
    'X_Al_M14'; ...
    'X_Mg_M14'; ...
    'X_FeLike_M14'; ...
    'X_h_M14'};

site = cell2struct(repmat({NaN}, numel(fieldNames), 1), fieldNames, 1);

end

function act = calcIdealActivities(site)
% calcIdealActivities
% Calculate ideal end-member activities after Bourdelle et al. (2013).
%
% Zero activity is allowed to produce a non-finite logarithm and ultimately
% a NaN temperature. NaN activity is retained. Negative finite activity is
% treated as a mathematically invalid internal result.

act = struct();

act.a_MgChlS = (site.X_Si_T2) .^ 2 ...
             .* (site.X_Mg_M23) .^ 4 ...
             .* (site.X_Mg_M14) .^ 2;

act.a_FeChlS = (site.X_Si_T2) .^ 2 ...
             .* (site.X_FeLike_M23) .^ 4 ...
             .* (site.X_FeLike_M14) .^ 2;

act.a_MgAm = (site.X_Al_T2) .^ 2 ...
           .* (site.X_Mg_M23) .^ 4 ...
           .* (site.X_Al_M14) .^ 2;

act.a_FeAm = (site.X_Al_T2) .^ 2 ...
           .* (site.X_FeLike_M23) .^ 4 ...
           .* (site.X_Al_M14) .^ 2;

act.a_MgSud = 256 ...
            .* (site.X_Si_T2) ...
            .* (site.X_Al_T2) ...
            .* (site.X_Al_M14) ...
            .* (site.X_h_M14) ...
            .* (site.X_Mg_M23) .^ 2 ...
            .* (site.X_Al_M23) .^ 2;

act.a_FeSud = 256 ...
            .* (site.X_Si_T2) ...
            .* (site.X_Al_T2) ...
            .* (site.X_Al_M14) ...
            .* (site.X_h_M14) ...
            .* (site.X_FeLike_M23) .^ 2 ...
            .* (site.X_Al_M23) .^ 2;

activityNames = fieldnames(act);
negativeMask = false(numel(activityNames), 1);

for i = 1:numel(activityNames)
    value = act.(activityNames{i});
    negativeMask(i) = isfinite(value) && value < 0;
end

if any(negativeMask)
    invalidNames = string(activityNames(negativeMask));
    error(['Negative finite activity encountered in Bourdelle2013: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function T_K = calcQuadraticTK(logK)
% calcQuadraticTK
% Auxiliary inversion of Bourdelle et al. (2013) Eq. (8):
%
%   logK = 11185729 / T(K)^2 - 56598 / T(K) + 72.3
%
% Rearranged:
%
%   (logK - 72.3) T^2 + 56598 T - 11185729 = 0
%
% The equation has no physical basis and should not be used outside
% 150–350 degreeC. It is returned only as an auxiliary output.

if ~isscalar(logK) || ~isfinite(logK)
    T_K = NaN;
    return;
end

a = logK - 72.3;
b = 56598;
c = -11185729;

if abs(a) < 1e-12
    T_K = NaN;
    return;
end

disc = b.^2 - 4 .* a .* c;

if disc < 0
    T_K = NaN;
    return;
end

root1 = (-b + sqrt(disc)) ./ (2 .* a);
root2 = (-b - sqrt(disc)) ./ (2 .* a);

candidates = [root1, root2];
validCandidate = isfinite(candidates) & candidates > 0;
candidates = candidates(validCandidate);

if isempty(candidates)
    T_K = NaN;
    return;
end

% Prefer a root within the stated 150–350 degreeC range.
inIntendedRange = candidates >= 423.15 & candidates <= 623.15;

if any(inIntendedRange)
    T_K = candidates(find(inIntendedRange, 1, 'first'));
else
    T_K = candidates(1);
end

end

function value = getVarOrError(tbl, varName, mineralLabel)
% getVarOrError
% Retrieve a required numeric scalar. NaN is allowed and retained, whereas
% infinite and negative values are rejected elsewhere.

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
