function results = Jowett1991(rawdata_struct, P_kbar)
% functions/+thermo/+Chlorite/Jowett1991.m
% Tested with MATLAB R2024b
%
% Empirical Fe-Mg-corrected single-Chlorite thermometer
% Jowett, E.C. (1991)
% Fitting Iron and Magnesium into the Hydrothermal Chlorite Geothermometer
% GAC/MAC/SEG Joint Annual Meeting, Program with Abstracts, 16, A62
% DOI: None
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Chlorite analysis at a time and
% calculates temperature using the Fe-Mg-corrected hydrothermal Chlorite
% geothermometer proposed by Jowett (1991).
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
% PUBLICATION STATUS:
% The thermometer was presented by Jowett (1991). The longer manuscript used
% here was accepted by Economic Geology after revision but was not published.
% The March 2016 author note states that the requested statistical
% verification of the isotherm slope coefficient a = 0.10 was not completed
% (2016 note, p. 1). Results should therefore be interpreted more cautiously
% than those from a fully published and independently reproduced calibration.
%
% The thermometer was field-calibrated using Chlorites from the Los Azufres
% and Salton Sea geothermal fields. The isothermal Fe-Mg correction and
% temperature equations are presented on pp. 7–8:
%
%   Al_IVC = Al_IVM + 0.10 * Fe/(Fe+Mg)
%   T(degreeC) = 318.5 * Al_IVC - 68.7
%
% The directly illustrated calibration/application ranges are (Fig. 5, p. 8):
%
%   Temperature    : 150–325 degreeC
%   Fe/(Fe+Mg)     : 0.2–0.6
%
% The conclusions describe application over 150–325 degreeC and
% Fe/(Fe+Mg) < 0.6 (p. 12). This implementation conservatively uses the
% directly illustrated 0.2–0.6 interval for the composition-range flag and
% warning.
%
% No quantitative pressure term or formal pressure calibration range is
% defined by Jowett (1991). In the application-limits discussion, Jowett
% cites the recommendation of Kranidiotis and MacLean (1987) that the method
% be used in low-pressure (<100 MPa, approximately <1 kbar),
% non-metamorphosed environments (p. 11). This is a cited prior-study
% recommendation rather than a Jowett calibration limit. The implementation
% therefore retains all pressures but prints a caution when any input
% pressure exceeds 1 kbar.
%
% Detailed petrography, composition, and paragenesis must be established
% before applying the thermometer with confidence (abstract, p. 2;
% applications, pp. 9–10; limits and conclusions, pp. 11–12). In
% particular:
%
%   - The analysed Chlorite must represent the geological event whose
%     temperature is sought.
%   - Relict material from minerals replaced by Chlorite may modify the
%     measured composition (p. 11).
%   - Chlorite/corrensite intergrowths and corrensite can produce strongly
%     underestimated and inconsistent temperatures. Pure corrensite should
%     not be treated as calibrated; mixed Chlorite/corrensite requires
%     caution and further calibration (pp. 11–12).
%   - The Si < 3.3 apfu and Ca < 0.07 apfu criteria cited from Shau et al.
%     (1990) are useful screening indicators for relatively pure Chlorite
%     (pp. 11–12), but they do not replace mineralogical confirmation.
%   - Al saturation is not imposed as a mandatory condition because the Los
%     Azufres calibration Chlorites were themselves Al-unsaturated on a Hey
%     diagram (p. 11).
%   - Structural formulae and the coefficient 0.10 are defined for a
%     14-oxygen half-formula basis (abstract, p. 2; equation 4, p. 8).
%   - The method uses the EPMA total-Fe term in Fe/(Fe+Mg) and does not
%     provide an independent Fe3+ correction. Oxidized Chlorites with
%     appreciable Fe3+ therefore require additional caution.
%   - Multiple analyses should be evaluated because application examples
%     show spreads of several tens of degrees within individual systems
%     (pp. 9–11).
%
% Two anomalous Salton Sea points near approximately 190 and 280–289 degreeC
% were excluded when defining the isothermal correction (pp. 7–8). This
% exclusion emphasizes that Chlorite composition is not controlled by
% temperature alone.
%
% This implementation issues non-stopping fprintf messages when:
%   1) a finite calculated temperature is outside 150–325 degreeC,
%   2) finite Fe/(Fe+Mg) is outside 0.2–0.6,
%   3) input pressure exceeds the cited low-pressure recommendation of
%      approximately 1 kbar, or when the absence of a formal Jowett pressure
%      calibration range must be stated,
%   4) Si >= 3.3 apfu or Ca >= 0.07 apfu suggests possible
%      corrensite/mixed-layer influence,
%   5) a calculation input contains NaN, or
%   6) a non-finite temperature is calculated.
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
% Required variables:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu      % total Fe used in Fe/(Fe+Mg)
%   Mg_cation_apfu
%
% Optional calculation variable:
%   Mn_cation_apfu      % used only in an auxiliary Fe/(Fe+Mg+Mn) index
%
% Optional trace/screening variables:
%   Ti_cation_apfu
%   Ca_cation_apfu      % used for the Ca < 0.07 Chlorite-purity screen
%   Na_cation_apfu
%   K_cation_apfu
%   Cr_cation_apfu
%
% All finite cation values must be greater than or equal to zero. Negative
% values are prohibited. NaN values are retained as missing values and are
% never replaced by zero when the corresponding column exists. Optional
% variables are assigned zero only when their columns are absent.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% On a 14-oxygen half-formula basis:
%
%   Al_IVM = min(Al_total, max(0, 4 - Si))
%
%   FeMg_ratio = Fe_total / (Fe_total + Mg)
%
%   Al_IVC = Al_IVM + 0.10 * FeMg_ratio
%
%   T(degreeC) = 318.5 * Al_IVC - 68.7
%
%   T(K) = T(degreeC) + 273.15
%
% Pressure is not used in the equation. A scalar or vector P_kbar is accepted
% only for compatibility with the fixed-pressure and pressure-range
% workflows and is stored in the output table.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Jowett1991(rawdata_struct, P_kbar)
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
% Basic checks prevent silent failures caused by missing arguments or invalid
% pressure values.
if nargin < 2
    error('Jowett1991 requires (rawdata_struct, P_kbar).');
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
% Store each calculation as one table block. This avoids repeatedly growing
% and copying the complete results table during the interactive loop.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Directly illustrated Jowett application ranges.
calibrationT_min_degC = 150;
calibrationT_max_degC = 325;
calibrationFeMg_min = 0.2;
calibrationFeMg_max = 0.6;

% Jowett gives no formal pressure calibration. The 1-kbar value is a
% low-pressure recommendation cited from Kranidiotis and MacLean (1987).
citedLowPressureLimit_kbar = 1.0;
pressureCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–4) Interactive selection loop + calculation
% Continue until the user cancels the selection or chooses Finish.
disp('=== Step 3: Selecting a data code from the list (Chlorite) ===');

while true
    % ----- Chlorite selection -----
    dataCodes_chl = dataset_chl{:, 1};

    [selectedIdx_chl, ok] = listdlg( ...
        'PromptString', 'Please select the Chlorite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_chl)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_chl)
        disp('Selection canceled');
        break;
    end

    selectedCode_chl = dataCodes_chl(selectedIdx_chl);
    disp(['Chlorite selected: ' char(string(selectedCode_chl))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    selectedData_chl = dataset_chl(selectedIdx_chl, :);

    % Identify existing NaN inputs without replacing them or stopping.
    nanInputNames = findNaNInputs(selectedData_chl);

    % Reject finite negative and infinite cation values. Zero and NaN are
    % allowed; NaN propagates through dependent outputs.
    validateNonNegativeInputs(selectedData_chl);

    row = calcTemp(selectedData_chl, P_kbar);

    % Add the selected identifier to every pressure row.
    row.dataCode_chl = repmat(string(selectedCode_chl), height(row), 1);
    row = movevars(row, {'dataCode_chl'}, 'Before', 1);

    % Add this table block to the preallocated buffer.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_chl)) ': T_Jowett1991 = ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_chl)) ': T_Jowett1991 = ' ...
            num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ...
            ' degreeC']);
    end

    % Explain the pressure limitation once per function call. The calculation
    % is never stopped because the equation contains no pressure term.
    if ~pressureCautionIssued
        pressureAboveCitedLimit = P_kbar > citedLowPressureLimit_kbar;

        if any(pressureAboveCitedLimit)
            fprintf(2, ...
                ['WARNING: %d of %d input pressure point(s) exceed 1 kbar ' ...
                 '(100 MPa); input range = %.4g–%.4g kbar. Jowett (1991, ' ...
                 'p. 11) cites a prior recommendation for low-pressure ' ...
                 '(<100 MPa), non-metamorphosed environments. This is not ' ...
                 'a formal Jowett pressure calibration range, and pressure ' ...
                 'is not used in the equation.\n'], ...
                sum(pressureAboveCitedLimit), ...
                numel(P_kbar), ...
                min(P_kbar), ...
                max(P_kbar));
        else
            fprintf(2, ...
                ['WARNING: Jowett (1991) does not define a formal numerical ' ...
                 'pressure calibration range, and pressure is not used in ' ...
                 'the equation. All supplied pressures are <=1 kbar, which ' ...
                 'is consistent with the low-pressure recommendation cited ' ...
                 'on p. 11; input range = %.4g–%.4g kbar.\n'], ...
                min(P_kbar), ...
                max(P_kbar));
        end

        pressureCautionIssued = true;
    end

    % Warn when a finite temperature lies outside the directly illustrated
    % 150–325 degreeC range.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the Jowett (1991) ' ...
             'field-calibration/application range of 150–325 degreeC ' ...
             '(Fig. 5, p. 8). %d of %d finite point(s) are outside; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    end

    % Warn when Fe/(Fe+Mg) lies outside the directly illustrated 0.2–0.6
    % composition interval.
    finiteFeMg = isfinite(row.Fe_ratio_Fe_over_FeMg);
    FeMgOutsideCalibration = finiteFeMg & ...
        (row.Fe_ratio_Fe_over_FeMg < calibrationFeMg_min | ...
         row.Fe_ratio_Fe_over_FeMg > calibrationFeMg_max);

    if any(FeMgOutsideCalibration)
        finiteValues = row.Fe_ratio_Fe_over_FeMg(finiteFeMg);
        fprintf(2, ...
            ['WARNING: Fe/(Fe+Mg) is outside the directly illustrated ' ...
             'Jowett (1991) range of 0.2–0.6 (Fig. 5, p. 8). ' ...
             'Calculated finite range = %.4g–%.4g for %s.\n'], ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    end

    % Screen for possible corrensite or mixed-layer influence using the
    % criteria cited by Jowett from Shau et al. (1990).
    finiteSi = isfinite(row.Si);
    finiteCa = isfinite(row.Ca);
    possibleCorrensiteInfluence = ...
        (finiteSi & row.Si >= 3.3) | ...
        (finiteCa & row.Ca >= 0.07);

    if any(possibleCorrensiteInfluence)
        fprintf(2, ...
            ['WARNING: Si >= 3.3 apfu and/or Ca >= 0.07 apfu was found for ' ...
             '%s. Jowett (1991, pp. 11–12) cites these as indicators of ' ...
             'possible corrensite or mixed Chlorite/corrensite influence. ' ...
             'Pure corrensite is not calibrated, and mineralogical and ' ...
             'petrographic confirmation is required.\n'], ...
            char(string(selectedCode_chl)));
    end

    % Report NaN inputs after the displayed result. NaN remains unchanged.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Jowett (1991) calculation ' ...
             'input(s) for %s: %s.\n' ...
             '         NaN values were retained and were not replaced by ' ...
             'zero. Dependent outputs may therefore be NaN.\n'], ...
            char(string(selectedCode_chl)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Report non-finite calculated temperatures without stopping.
    invalidTemperature = ~isfinite(row.T_deg);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite Jowett (1991) temperature values were ' ...
             'calculated for %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
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
        'Jowett1991', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all result blocks once after the loop.
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
% Return names of existing variables used in the thermometer or screening
% calculations that contain NaN.

variableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ca_cation_apfu'};

displayNames = "Chlorite." + string(variableNames(:));
nanMask = false(numel(variableNames), 1);

for i = 1:numel(variableNames)
    variableName = variableNames{i};

    if ismember(variableName, data_chl.Properties.VariableNames)
        variableValue = data_chl.(variableName);

        if isnumeric(variableValue)
            nanMask(i) = any(isnan(variableValue(:)));
        end
    end
end

nanInputNames = displayNames(nanMask);

end

function validateNonNegativeInputs(data_chl)
% validateNonNegativeInputs
% Reject finite negative values, infinite values, and non-numeric cation
% columns. Zero and NaN are allowed.

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
    error(['Jowett1991: cation variables must be numeric. ' ...
           'Non-numeric variable(s): ' ...
           char(strjoin(displayNames(nonNumericMask), ', ')) '.']);
end

if any(negativeMask)
    error(['Jowett1991: cation values must be greater than or equal to ' ...
           'zero. Negative value(s) were found in: ' ...
           char(strjoin(displayNames(negativeMask), ', ')) '.']);
end

if any(infiniteMask)
    error(['Jowett1991: infinite cation value(s) are not allowed in: ' ...
           char(strjoin(displayNames(infiniteMask), ', ')) '.']);
end

end

function row = calcTemp(data_chl, P_kbar)
% calcTemp
% Compute one composition-derived temperature and repeat all scalar outputs
% for each supplied pressure value.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

% --- Prepare Chlorite data and site allocation ---
chl = prepareChloriteRow(data_chl);
site = calcChloriteSites(chl);

% --- Jowett Fe-Mg correction and temperature ---
% Division by zero and NaN values propagate naturally to NaN.
FeMg_ratio_scalar = chl.FeT ./ (chl.FeT + chl.Mg);
Mg_number_scalar = chl.Mg ./ (chl.Mg + chl.FeT);
Fe_ratio_Fe_over_FeMgMn_scalar = ...
    chl.FeT ./ (chl.FeT + chl.Mg + chl.Mn);

Al_IVC_scalar = site.Al_IV + 0.10 .* FeMg_ratio_scalar;
T_deg_scalar = 318.5 .* Al_IVC_scalar - 68.7;
T_K_scalar = T_deg_scalar + 273.15;

% --- Applicability and sanity flags ---
is_14O_tetra_reasonable_scalar = ...
    isfinite(chl.Si) && chl.Si >= 2.0 && chl.Si <= 4.0;

is_AlIV_reasonable_scalar = ...
    isfinite(site.Al_IV) && site.Al_IV >= 0.0 && site.Al_IV <= 2.0;

is_FeMg_ratio_reasonable_scalar = ...
    isfinite(FeMg_ratio_scalar) && ...
    FeMg_ratio_scalar >= 0.0 && FeMg_ratio_scalar <= 1.0;

is_in_Jowett_T_range_scalar = ...
    isfinite(T_deg_scalar) && ...
    T_deg_scalar >= 150 && T_deg_scalar <= 325;

is_in_Jowett_FeMg_range_scalar = ...
    isfinite(FeMg_ratio_scalar) && ...
    FeMg_ratio_scalar >= 0.2 && FeMg_ratio_scalar <= 0.6;

is_Si_pureChlorite_screen_scalar = ...
    isfinite(chl.Si) && chl.Si < 3.3;

is_Ca_pureChlorite_screen_scalar = ...
    isfinite(chl.Ca) && chl.Ca < 0.07;

recommended_by_Jowett1991 = ...
    repmat(is_14O_tetra_reasonable_scalar, nP, 1) & ...
    repmat(is_AlIV_reasonable_scalar, nP, 1) & ...
    repmat(is_FeMg_ratio_reasonable_scalar, nP, 1) & ...
    repmat(is_in_Jowett_T_range_scalar, nP, 1) & ...
    repmat(is_in_Jowett_FeMg_range_scalar, nP, 1) & ...
    repmat(is_Si_pureChlorite_screen_scalar, nP, 1) & ...
    repmat(is_Ca_pureChlorite_screen_scalar, nP, 1);

% --- Pack outputs ---
row = table();

row.P_kbar = P_kbar;

row.Si = repmat(chl.Si, nP, 1);
row.Al = repmat(chl.Al, nP, 1);
row.Fe_total = repmat(chl.FeT, nP, 1);
row.Fe2 = repmat(chl.FeT, nP, 1);
row.Fe2_assumed = repmat(chl.FeT, nP, 1);
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

row.Fe_ratio_Fe_over_FeMg = repmat(FeMg_ratio_scalar, nP, 1);
row.Mg_number = repmat(Mg_number_scalar, nP, 1);
row.Fe_ratio_Fe_over_FeMgMn = ...
    repmat(Fe_ratio_Fe_over_FeMgMn_scalar, nP, 1);

row.Al_IVC = repmat(Al_IVC_scalar, nP, 1);

row.T_deg = repmat(T_deg_scalar, nP, 1);
row.T_K = repmat(T_K_scalar, nP, 1);

row.is_14O_tetra_reasonable = ...
    repmat(is_14O_tetra_reasonable_scalar, nP, 1);
row.is_AlIV_reasonable = ...
    repmat(is_AlIV_reasonable_scalar, nP, 1);
row.is_FeMg_ratio_reasonable = ...
    repmat(is_FeMg_ratio_reasonable_scalar, nP, 1);
row.is_in_Jowett_T_range = ...
    repmat(is_in_Jowett_T_range_scalar, nP, 1);
row.is_in_Jowett_FeMg_range = ...
    repmat(is_in_Jowett_FeMg_range_scalar, nP, 1);
row.is_Si_pureChlorite_screen = ...
    repmat(is_Si_pureChlorite_screen_scalar, nP, 1);
row.is_Ca_pureChlorite_screen = ...
    repmat(is_Ca_pureChlorite_screen_scalar, nP, 1);
row.pressure_calibration_defined = false(nP, 1);
row.is_within_cited_lowP_recommendation = ...
    P_kbar <= 1.0;
row.recommended_by_Jowett1991 = recommended_by_Jowett1991;

end

function chl = prepareChloriteRow(data_chl)
% prepareChloriteRow
% Extract one Chlorite analysis from a 1-row table. Existing NaN values are
% retained. Optional variables are assigned zero only when their columns are
% absent.

if height(data_chl) ~= 1
    error('Chlorite input must be a 1-row table.');
end

chl = struct();

chl.Si  = getVarOrError(data_chl, 'Si_cation_apfu', 'Chlorite');
chl.Al  = getVarOrError(data_chl, 'Al_cation_apfu', 'Chlorite');
chl.FeT = getVarOrError(data_chl, 'Fe_cation_apfu', 'Chlorite');
chl.Mg  = getVarOrError(data_chl, 'Mg_cation_apfu', 'Chlorite');
chl.Mn  = getVarOrZeroIfMissing(data_chl, 'Mn_cation_apfu');

chl.Ti = getVarOrZeroIfMissing(data_chl, 'Ti_cation_apfu');
chl.Ca = getVarOrZeroIfMissing(data_chl, 'Ca_cation_apfu');
chl.Na = getVarOrZeroIfMissing(data_chl, 'Na_cation_apfu');
chl.K  = getVarOrZeroIfMissing(data_chl, 'K_cation_apfu');
chl.Cr = getVarOrZeroIfMissing(data_chl, 'Cr_cation_apfu');

end

function site = calcChloriteSites(chl)
% calcChloriteSites
% Calculate approximate site quantities on a 14-oxygen basis. NaN values
% remain NaN and do not trigger replacement with zero.

site = struct( ...
    'Al_IV', NaN, ...
    'Al_VI', NaN, ...
    'Sum_VI', NaN, ...
    'VAC', NaN);

% Preserve NaN explicitly because min/max handling can vary with options and
% MATLAB versions.
if isnan(chl.Si) || isnan(chl.Al)
    return;
end

site.Al_IV = min(chl.Al, max(0, 4 - chl.Si));
site.Al_VI = chl.Al - site.Al_IV;

% A finite negative Al_VI indicates an internally inconsistent structural
% formula. Return NaN outputs without stopping the complete workflow.
if ~isfinite(site.Al_VI) || site.Al_VI < -1e-10
    site.Al_IV = NaN;
    site.Al_VI = NaN;
    return;
end

site.Al_VI = max(0, site.Al_VI);
site.Sum_VI = site.Al_VI + chl.FeT + chl.Mg + chl.Mn;
site.VAC = 6 - site.Sum_VI;

end

function value = getVarOrError(tbl, varName, mineralLabel)
% getVarOrError
% Retrieve a required numeric scalar. NaN is retained. Infinite and finite
% negative values are prohibited.

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
% absent. Existing NaN values remain NaN.

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
