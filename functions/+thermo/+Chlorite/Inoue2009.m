function results = Inoue2009(rawdata_struct, P_kbar)
% functions/+thermo/+Chlorite/Inoue2009.m
% Tested with MATLAB R2024b
%
% Semi-empirical single-Chlorite thermometer
% Inoue, A., Meunier, A., Patrier-Mas, P., Rigault, C.,
% Beaufort, D. and Vieillard, P. (2009)
% Clays and Clay Minerals, 57, 371–382
% DOI: https://doi.org/10.1346/CCMN.2009.0570309
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Chlorite analysis at a time and
% calculates temperature using the Inoue et al. (2009) low-temperature
% trioctahedral-Chlorite geothermometer.
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
% Inoue et al. (2009) developed the thermometer from natural
% low-temperature trioctahedral Chlorites from diagenetic, hydrothermal, and
% low-grade metamorphic settings, including Niger, Rouez, and St Martin
% (paper overview and data discussion, pp. 371–382).
%
% The paper does not define a single strict numerical calibration interval
% equivalent to an experimental calibration range. This implementation uses
% 100–350 degreeC only as a practical low-temperature screening interval.
% It must not be interpreted as a formally stated experimental calibration
% range of Inoue et al. (2009). Temperatures outside this interval are
% retained but reported by non-stopping fprintf warnings.
%
% A quantitative pressure calibration range is not defined, and pressure
% does not appear explicitly in the thermometer equation. Input pressure is
% therefore stored for interface compatibility and traceability, but its
% validity cannot be assessed numerically from this calibration.
%
% The most important application constraint is the ferric-iron correction.
% Inoue et al. (2009) used Chlorites for which Fe3+/SigmaFe information was
% available and showed that assuming all Fe as Fe2+ can cause substantial
% temperature overestimation. The relevant Fe-valence treatment and
% thermometer comparisons are discussed throughout the paper, particularly
% in the methods, results, and discussion sections (approximately
% pp. 374–381).
%
% Consequently:
%
%   - Fe3_cation_apfu should be supplied independently whenever possible.
%   - If the Fe3 column is absent, this implementation assigns NaN rather
%     than zero, so that an apparently precise but unsupported temperature is
%     not calculated.
%   - If the Fe3 column exists and contains NaN, the NaN is retained.
%   - Fe3 must not exceed total Fe.
%
% The thermometer is based on a Chlorite + quartz + water reaction and
% assumes:
%
%   a_quartz = 1
%   a_H2O    = 1
%
% It is therefore most appropriate for quartz-bearing or silica-saturated,
% water-rich systems. Application to quartz-absent systems or fluids with
% substantially reduced water activity requires caution.
%
% Structural formulae must be supplied on a 14-oxygen half-formula basis.
% The model requires positive site fractions and positive ideal activities
% for the logarithmic equilibrium expression. Analyses affected by mixed
% layering, intergrowths, inherited material, incomplete equilibration, or
% unusual compositions may yield non-evaluable activities and NaN
% temperatures.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) a finite temperature is outside the practical 100–350 degreeC
%      low-temperature screening interval,
%   2) no quantitative pressure calibration range can be evaluated,
%   3) Fe3_cation_apfu is absent,
%   4) a calculation input contains NaN,
%   5) the site-activity expression is not evaluable, or
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
%   Fe_cation_apfu       % total Fe
%   Mg_cation_apfu
%
% Strongly recommended calculation variable:
%   Fe3_cation_apfu      % Fe3+; absent column is treated as unknown (NaN)
%
% Optional calculation variable:
%   Mn_cation_apfu       % treated as Mg-equivalent; zero only if absent
%
% Optional trace variables:
%   Ti_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%   Cr_cation_apfu
%
% All finite cation values must be greater than or equal to zero. Negative
% values are prohibited. NaN values are retained as missing values and are
% never replaced by zero when the corresponding column exists.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Reaction:
%
%   3 Crdp + Afch = 3 Sud + 7 Qtz + 4 H2O
%
% where:
%   Crdp = corundophilite component
%   Afch = Al-free Chlorite component
%   Sud  = sudoite component
%
% Equilibrium variable:
%
%   x = log10(K)
%     = 4*log10(a_Crdp_ideal)
%       - 3*log10(a_Sud_ideal)
%       - log10(a_Afch_ideal)
%
% Temperature:
%
%   T(K) = 1 / (0.00293 - 5.13e-4*x + 3.904e-5*x^2)
%   T(degreeC) = T(K) - 273.15
%
% Ideal activities:
%
%   a_Afch = (X_Mg_oct)^6 * (X_Si_tet)^2
%
%   a_Crdp = 45.563 * (X_Mg_oct)^4 * (X_Al_oct)^2
%                    * (X_Al_tet)^2
%
%   a_Sud  = 1728 * (X_Mg_oct)^2 * (X_Al_oct)^3
%                   * X_Vac_oct * X_Si_tet * X_Al_tet
%
% Site allocation on a 14-oxygen basis:
%
%   Al_IV = min(Al_total, max(0, 4 - Si))
%   Al_VI = Al_total - Al_IV
%   Fe2   = Fe_total - Fe3
%   Mg_equiv = Mg + Mn
%   Sum_VI = Al_VI + Fe2 + Fe3 + Mg_equiv
%   VAC = 6 - Sum_VI
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Inoue2009(rawdata_struct, P_kbar)
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
if nargin < 2
    error('Inoue2009 requires (rawdata_struct, P_kbar).');
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
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Practical screening interval; not a formally stated experimental range.
practicalT_min_degC = 100;
practicalT_max_degC = 350;

% No numerical pressure calibration range is defined in the paper.
pressureCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–4) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Chlorite) ===');

while true
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

    disp('=== Step 4: Calculating the temperature ===');

    selectedData_chl = dataset_chl(selectedIdx_chl, :);

    nanInputNames = findNaNInputs(selectedData_chl);
    validateNonNegativeInputs(selectedData_chl);

    row = calcTemp(selectedData_chl, P_kbar);

    row.dataCode_chl = repmat(string(selectedCode_chl), height(row), 1);
    row = movevars(row, {'dataCode_chl'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_chl)) ': T_Inoue2009 = ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_chl)) ': T_Inoue2009 = ' ...
            num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ...
            ' degreeC']);
    end

    if ~pressureCautionIssued
        fprintf(2, ...
            ['WARNING: Inoue et al. (2009) do not define a quantitative ' ...
             'pressure calibration range, and pressure is not used explicitly ' ...
             'in the thermometer equation. The supplied pressure is stored ' ...
             'only for traceability; input range = %.4g–%.4g kbar. ' ...
             'Pressure-range validity cannot be evaluated numerically.\n'], ...
            min(P_kbar), max(P_kbar));
        pressureCautionIssued = true;
    end

    finiteTemperature = isfinite(row.T_deg);
    outsidePracticalRange = finiteTemperature & ...
        (row.T_deg < practicalT_min_degC | ...
         row.T_deg > practicalT_max_degC);

    if any(outsidePracticalRange)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the practical ' ...
             '100–350 degreeC low-temperature screening interval used in ' ...
             'this implementation of Inoue et al. (2009). This interval is ' ...
             'not a formally stated experimental calibration range. ' ...
             '%d of %d finite point(s) are outside; calculated finite ' ...
             'range = %.4g–%.4g degreeC for %s.\n'], ...
            sum(outsidePracticalRange), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    end

    if ~row.has_Fe3_input(1)
        fprintf(2, ...
            ['WARNING: Fe3_cation_apfu is absent for %s. Inoue et al. ' ...
             '(2009) require ferric-iron information for reliable use. ' ...
             'Fe3 was treated as unknown (NaN), not as zero, and dependent ' ...
             'outputs may therefore be NaN.\n'], ...
            char(string(selectedCode_chl)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Inoue et al. (2009) input(s) ' ...
             'for %s: %s.\n' ...
             '         NaN values were retained and were not replaced by ' ...
             'zero. Dependent outputs may therefore be NaN.\n'], ...
            char(string(selectedCode_chl)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    if any(~row.is_logK_evaluable)
        fprintf(2, ...
            ['WARNING: The Inoue et al. (2009) logarithmic activity ' ...
             'expression was not evaluable for %s. Required site fractions ' ...
             'or ideal activities were non-positive or non-finite. The ' ...
             'temperature was retained as NaN and calculation continued.\n'], ...
            char(string(selectedCode_chl)));
    end

    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite Inoue et al. (2009) temperature values ' ...
             'were calculated for %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_chl)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Inoue2009', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

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
% Return names of existing calculation variables that contain NaN.

variableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Fe3_cation_apfu'};

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
% Reject negative finite values and infinite values. Zero and NaN are
% allowed; NaN is retained for propagation.

variableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Fe3_cation_apfu', ...
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
    error(['Inoue2009: cation variables must be numeric. ' ...
           'Non-numeric variable(s): ' ...
           char(strjoin(displayNames(nonNumericMask), ', ')) '.']);
end

if any(negativeMask)
    error(['Inoue2009: cation values must be greater than or equal to ' ...
           'zero. Negative value(s) were found in: ' ...
           char(strjoin(displayNames(negativeMask), ', ')) '.']);
end

if any(infiniteMask)
    error(['Inoue2009: infinite cation value(s) are not allowed in: ' ...
           char(strjoin(displayNames(infiniteMask), ', ')) '.']);
end

end

function row = calcTemp(data_chl, P_kbar)
% calcTemp
% Compute one composition-derived temperature and repeat the output for each
% supplied pressure value.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

chl = prepareChloriteRow(data_chl);
site = calcChloriteSites(chl);
act = calcIdealActivities(site);

is_logK_evaluable_scalar = ...
    isfinite(act.a_Afch_ideal) && act.a_Afch_ideal > 0 && ...
    isfinite(act.a_Crdp_ideal) && act.a_Crdp_ideal > 0 && ...
    isfinite(act.a_Sud_ideal)  && act.a_Sud_ideal  > 0;

if is_logK_evaluable_scalar
    x_logK7_scalar = 4 .* log10(act.a_Crdp_ideal) ...
                   - 3 .* log10(act.a_Sud_ideal) ...
                   - log10(act.a_Afch_ideal);

    denominator = 0.00293 ...
                - 5.13e-4 .* x_logK7_scalar ...
                + 3.904e-5 .* x_logK7_scalar.^2;

    if isfinite(denominator) && denominator > 0
        T_K_scalar = 1 ./ denominator;
        T_deg_scalar = T_K_scalar - 273.15;
    else
        x_logK7_scalar = NaN;
        T_K_scalar = NaN;
        T_deg_scalar = NaN;
        is_logK_evaluable_scalar = false;
    end
else
    x_logK7_scalar = NaN;
    T_K_scalar = NaN;
    T_deg_scalar = NaN;
end

Mg_number_scalar = chl.Mg_equiv ./ (chl.Mg_equiv + chl.Fe2);
Fe_ratio_scalar = chl.Fe2 ./ (chl.Fe2 + chl.Mg + chl.Mn);
Fe3_ratio_scalar = chl.Fe3 ./ chl.FeT;

is_14O_tetra_reasonable_scalar = ...
    isfinite(chl.Si) && chl.Si >= 2.0 && chl.Si <= 4.0;
is_AlIV_reasonable_scalar = ...
    isfinite(site.Al_IV) && site.Al_IV >= 0.0 && site.Al_IV <= 2.0;
is_VAC_nonnegative_scalar = ...
    isfinite(site.VAC) && site.VAC >= 0.0;
is_VAC_positive_scalar = ...
    isfinite(site.VAC) && site.VAC > 0.0;

is_XSi_tet_reasonable_scalar = ...
    isfinite(site.X_Si_tet) && site.X_Si_tet >= 0.0 && site.X_Si_tet <= 1.0;
is_XAl_tet_reasonable_scalar = ...
    isfinite(site.X_Al_tet) && site.X_Al_tet >= 0.0 && site.X_Al_tet <= 1.0;
is_XMg_oct_reasonable_scalar = ...
    isfinite(site.X_Mg_oct) && site.X_Mg_oct >= 0.0 && site.X_Mg_oct <= 1.0;
is_XAl_oct_reasonable_scalar = ...
    isfinite(site.X_Al_oct) && site.X_Al_oct >= 0.0 && site.X_Al_oct <= 1.0;

isfinite_temperature_scalar = isfinite(T_deg_scalar) && isreal(T_deg_scalar);
is_in_practical_lowT_range_scalar = ...
    isfinite(T_deg_scalar) && T_deg_scalar >= 100 && T_deg_scalar <= 350;

recommended_by_Inoue2009 = ...
    repmat(chl.has_Fe3_input, nP, 1) & ...
    repmat(is_14O_tetra_reasonable_scalar, nP, 1) & ...
    repmat(is_AlIV_reasonable_scalar, nP, 1) & ...
    repmat(is_VAC_positive_scalar, nP, 1) & ...
    repmat(is_XSi_tet_reasonable_scalar, nP, 1) & ...
    repmat(is_XAl_tet_reasonable_scalar, nP, 1) & ...
    repmat(is_XMg_oct_reasonable_scalar, nP, 1) & ...
    repmat(is_XAl_oct_reasonable_scalar, nP, 1) & ...
    repmat(is_logK_evaluable_scalar, nP, 1) & ...
    repmat(isfinite_temperature_scalar, nP, 1) & ...
    repmat(is_in_practical_lowT_range_scalar, nP, 1);

row = table();

row.P_kbar = P_kbar;

row.Si = repmat(chl.Si, nP, 1);
row.Al = repmat(chl.Al, nP, 1);
row.FeT = repmat(chl.FeT, nP, 1);
row.Fe2 = repmat(chl.Fe2, nP, 1);
row.Fe3 = repmat(chl.Fe3, nP, 1);
row.Mg = repmat(chl.Mg, nP, 1);
row.Mn = repmat(chl.Mn, nP, 1);
row.Mg_equiv_for_thermometer = repmat(chl.Mg_equiv, nP, 1);
row.Ti = repmat(chl.Ti, nP, 1);
row.Ca = repmat(chl.Ca, nP, 1);
row.Na = repmat(chl.Na, nP, 1);
row.K = repmat(chl.K, nP, 1);
row.Cr = repmat(chl.Cr, nP, 1);

row.Al_IV = repmat(site.Al_IV, nP, 1);
row.Al_VI = repmat(site.Al_VI, nP, 1);
row.Sum_VI = repmat(site.Sum_VI, nP, 1);
row.VAC = repmat(site.VAC, nP, 1);

row.X_Si_tet = repmat(site.X_Si_tet, nP, 1);
row.X_Al_tet = repmat(site.X_Al_tet, nP, 1);
row.X_Mg_oct = repmat(site.X_Mg_oct, nP, 1);
row.X_Al_oct = repmat(site.X_Al_oct, nP, 1);
row.X_Vac_oct = repmat(site.X_Vac_oct, nP, 1);

row.a_Afch_ideal = repmat(act.a_Afch_ideal, nP, 1);
row.a_Crdp_ideal = repmat(act.a_Crdp_ideal, nP, 1);
row.a_Sud_ideal = repmat(act.a_Sud_ideal, nP, 1);
row.x_logK7 = repmat(x_logK7_scalar, nP, 1);

row.Mg_number = repmat(Mg_number_scalar, nP, 1);
row.Fe_ratio_Fe2_over_FeMgMn = repmat(Fe_ratio_scalar, nP, 1);
row.Fe3_ratio_Fe3_over_FeT = repmat(Fe3_ratio_scalar, nP, 1);

row.T_deg = repmat(T_deg_scalar, nP, 1);
row.T_K = repmat(T_K_scalar, nP, 1);

row.has_Fe3_input = repmat(chl.has_Fe3_input, nP, 1);
row.is_14O_tetra_reasonable = ...
    repmat(is_14O_tetra_reasonable_scalar, nP, 1);
row.is_AlIV_reasonable = repmat(is_AlIV_reasonable_scalar, nP, 1);
row.is_VAC_nonnegative = repmat(is_VAC_nonnegative_scalar, nP, 1);
row.is_VAC_positive = repmat(is_VAC_positive_scalar, nP, 1);
row.is_XSi_tet_reasonable = ...
    repmat(is_XSi_tet_reasonable_scalar, nP, 1);
row.is_XAl_tet_reasonable = ...
    repmat(is_XAl_tet_reasonable_scalar, nP, 1);
row.is_XMg_oct_reasonable = ...
    repmat(is_XMg_oct_reasonable_scalar, nP, 1);
row.is_XAl_oct_reasonable = ...
    repmat(is_XAl_oct_reasonable_scalar, nP, 1);
row.is_logK_evaluable = repmat(is_logK_evaluable_scalar, nP, 1);
row.isfinite_temperature = repmat(isfinite_temperature_scalar, nP, 1);
row.is_in_practical_lowT_range = ...
    repmat(is_in_practical_lowT_range_scalar, nP, 1);
row.pressure_calibration_defined = false(nP, 1);
row.recommended_by_Inoue2009 = recommended_by_Inoue2009;

end

function chl = prepareChloriteRow(data_chl)
% prepareChloriteRow
% Extract one Chlorite analysis. Existing NaN values are retained. Optional
% variables are assigned zero only when their columns are absent.

if height(data_chl) ~= 1
    error('Chlorite input must be a 1-row table.');
end

chl = struct();

chl.Si  = getVarOrError(data_chl, 'Si_cation_apfu', 'Chlorite');
chl.Al  = getVarOrError(data_chl, 'Al_cation_apfu', 'Chlorite');
chl.FeT = getVarOrError(data_chl, 'Fe_cation_apfu', 'Chlorite');
chl.Mg  = getVarOrError(data_chl, 'Mg_cation_apfu', 'Chlorite');
chl.Mn  = getVarOrZeroIfMissing(data_chl, 'Mn_cation_apfu');

if ismember('Fe3_cation_apfu', data_chl.Properties.VariableNames)
    chl.Fe3 = getVarOrError(data_chl, 'Fe3_cation_apfu', 'Chlorite');
    chl.has_Fe3_input = true;
else
    chl.Fe3 = NaN;
    chl.has_Fe3_input = false;
end

chl.Ti = getVarOrZeroIfMissing(data_chl, 'Ti_cation_apfu');
chl.Ca = getVarOrZeroIfMissing(data_chl, 'Ca_cation_apfu');
chl.Na = getVarOrZeroIfMissing(data_chl, 'Na_cation_apfu');
chl.K  = getVarOrZeroIfMissing(data_chl, 'K_cation_apfu');
chl.Cr = getVarOrZeroIfMissing(data_chl, 'Cr_cation_apfu');

if isfinite(chl.Fe3) && isfinite(chl.FeT) && chl.Fe3 > chl.FeT + 1e-10
    error('Fe3_cation_apfu cannot exceed Fe_cation_apfu.');
end

chl.Fe2 = chl.FeT - chl.Fe3;
chl.Mg_equiv = chl.Mg + chl.Mn;

end

function site = calcChloriteSites(chl)
% calcChloriteSites
% Calculate Chlorite site quantities on a 14-oxygen basis. Invalid or
% missing derived quantities are retained as NaN rather than stopping.

site = makeNaNSiteStruct();

if isnan(chl.Si) || isnan(chl.Al)
    return;
end

site.Al_IV = min(chl.Al, max(0, 4 - chl.Si));
site.Al_VI = chl.Al - site.Al_IV;

if ~isfinite(site.Al_VI) || site.Al_VI < -1e-10
    site = makeNaNSiteStruct();
    return;
end

site.Al_VI = max(0, site.Al_VI);
site.Sum_VI = site.Al_VI + chl.Fe2 + chl.Fe3 + chl.Mg_equiv;
site.VAC = 6 - site.Sum_VI;

if isfinite(site.VAC) && site.VAC > -1e-8 && site.VAC < 0
    site.VAC = 0;
end

site.X_Si_tet = (chl.Si - 2) ./ 2;
site.X_Al_tet = site.Al_IV ./ 2;
site.X_Mg_oct = chl.Mg_equiv ./ 6;
site.X_Al_oct = site.Al_VI ./ 6;
site.X_Vac_oct = site.VAC ./ 6;

end

function site = makeNaNSiteStruct()
% makeNaNSiteStruct
% Return a site-quantity structure whose fields are all NaN.

site = struct( ...
    'Al_IV', NaN, ...
    'Al_VI', NaN, ...
    'Sum_VI', NaN, ...
    'VAC', NaN, ...
    'X_Si_tet', NaN, ...
    'X_Al_tet', NaN, ...
    'X_Mg_oct', NaN, ...
    'X_Al_oct', NaN, ...
    'X_Vac_oct', NaN);

end

function act = calcIdealActivities(site)
% calcIdealActivities
% Calculate ideal activities. Return NaN values when logarithmic activity
% requirements are not satisfied.

act = struct( ...
    'a_Afch_ideal', NaN, ...
    'a_Crdp_ideal', NaN, ...
    'a_Sud_ideal', NaN);

isValid = ...
    isfinite(site.X_Si_tet) && site.X_Si_tet > 0 && ...
    isfinite(site.X_Al_tet) && site.X_Al_tet > 0 && ...
    isfinite(site.X_Mg_oct) && site.X_Mg_oct > 0 && ...
    isfinite(site.X_Al_oct) && site.X_Al_oct > 0 && ...
    isfinite(site.X_Vac_oct) && site.X_Vac_oct > 0;

if ~isValid
    return;
end

act.a_Afch_ideal = ...
    (site.X_Mg_oct.^6) .* (site.X_Si_tet.^2);

act.a_Crdp_ideal = ...
    45.563 .* (site.X_Mg_oct.^4) .* ...
    (site.X_Al_oct.^2) .* (site.X_Al_tet.^2);

act.a_Sud_ideal = ...
    1728 .* (site.X_Mg_oct.^2) .* ...
    (site.X_Al_oct.^3) .* site.X_Vac_oct .* ...
    site.X_Si_tet .* site.X_Al_tet;

end

function value = getVarOrError(tbl, varName, mineralLabel)
% getVarOrError
% Retrieve a required numeric scalar. NaN is retained; Inf and negative
% finite values are rejected.

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
% Retrieve an optional numeric scalar. Zero is assigned only if the column is
% absent; an existing NaN is retained.

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
