function results = Inoue2018(rawdata_struct, P_kbar)
% functions/+thermo/+Chlorite/Inoue2018.m
% Tested with MATLAB R2024b
%
% Single-Chlorite thermometer and redox estimator
% Inoue, A., Inoue, S. and Utada, M. (2018)
% Clay Minerals, 53, 143–158
% DOI: https://doi.org/10.1180/clm.2018.10
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Chlorite analysis at a time and
% applies the practical temperature and redox workflow described by
% Inoue et al. (2018).
%
% The implemented workflow is:
%
%   1) Calculate logK(Fe2+) assuming total Fe = Fe2+.
%   2) Convert logK(Fe2+) to temperature using Eq. (6).
%   3) Estimate Fe3+/SigmaFe by matching the temperature from Appendix
%      Eq. (A3) to the Eq. (6) temperature.
%   4) Estimate log fO2 using the Walshe-type redox relation.
%
% Both a scalar pressure and a pressure vector are accepted. Therefore, this
% function can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. Pressure is retained in the output for
% traceability, although it is not explicitly used in the temperature
% equation.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Inoue et al. (2018) is not an independent experimental calibration with a
% single formally stated temperature and pressure range. Equation (6) was
% fitted using the Niger, Rouez, and St Martin data previously treated with
% the Inoue et al. (2009) thermometer (Eq. 6 and Fig. 5, pp. 150–151).
%
% The Noboribetsu application produced average temperatures of approximately
% 205–299 degreeC (Table 1 and discussion, pp. 148–151), but this observed
% result is not a universal formal calibration interval.
%
% This implementation uses 100–350 degreeC only as a practical
% low-temperature screening interval. Temperatures outside this interval are
% retained and reported by non-stopping fprintf warnings. The interval must
% not be described as a formally stated experimental calibration range.
%
% A quantitative pressure calibration range is not defined for the
% temperature equation. Pressure effects on the redox calculations are
% discussed on pp. 145–146 and were found to be comparatively insensitive
% for the reactions considered. The Noboribetsu application assumed
% pressures near the liquid-vapour saturation curve, possibly below
% approximately 200 bar (p. 151), but this is a site-specific geological
% condition rather than a general calibration range.
%
% Important application constraints include:
%
%   - Structural formulae are calculated on a 14-oxygen half-formula basis
%     (Appendix, pp. 156–158).
%   - The temperature and activity model assumes quartz and water activities
%     of unity: a_quartz = 1 and a_H2O = 1 (Appendix, pp. 156–158).
%   - If a_H2O is less than one, the calculated log fO2 changes. The water
%     activity was unknown for the Noboribetsu system (p. 152).
%   - The estimated Fe3+/SigmaFe value is obtained through model matching
%     rather than direct Fe-valence measurement (pp. 150–151 and Appendix).
%   - Quantitative log fO2 estimation requires further work; the method is
%     most defensible for relative or qualitative redox comparisons
%     (concluding remarks, p. 154).
%   - The amount of Fe3+ accommodated by trioctahedral Chlorite is
%     structurally constrained and remains controversial, especially for
%     very Fe-poor analyses (pp. 153–154).
%   - Ti, Mn, alkali, and alkaline-earth elements are ignored in the
%     Appendix activity calculation (pp. 156–158). Unusual compositions rich
%     in these elements require caution.
%
% IMPORTANT IMPLEMENTATION LIMITATION:
% The original Appendix repartitions FeO* into FeO and Fe2O3 and then
% recalculates the 14-oxygen structural formula. The thermo.Chlorite workflow
% used here receives cation apfu directly. Therefore, this implementation
% keeps the supplied Si, Al, Mg, and total-Fe apfu fixed while partitioning
% total Fe. It is a practical cation-apfu approximation and is not identical
% to strict oxide-wt.% repartitioning and renormalization.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) a finite temperature is outside the practical 100–350 degreeC
%      screening interval,
%   2) no quantitative pressure calibration range can be evaluated,
%   3) a calculation input contains NaN,
%   4) the temperature or activity calculation is not evaluable,
%   5) Fe3+/SigmaFe matching does not converge adequately,
%   6) estimated Fe3+/SigmaFe is unusually high, or
%   7) a non-finite temperature, Fe3 ratio, or log fO2 is calculated.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain either:
%   rawdata_struct.Chlorite : table
% or
%   rawdata_struct.Chl      : table
%
% The FIRST column is treated as an identifier ("data code").
%
% Required variables:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu       % total Fe
%   Mg_cation_apfu
%
% Optional variables:
%   Mn_cation_apfu
%   Ti_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%   Cr_cation_apfu
%
% All finite cation values must be greater than or equal to zero. Negative
% values are prohibited. Existing NaN values are retained and are never
% replaced by zero. Optional variables are assigned zero only if their
% columns are absent.
%
% -------------------------------------------------------------------------
% IMPLEMENTED EQUATIONS
%
% [1] Inoue et al. (2018), Eq. (6), p. 150:
%
%   T(degreeC) = 63.83 + 50.41*x + 2.617*x^2 + 2.846*x^3
%                - 1.097*x^4 + 0.09285*x^5
%
%   x = logK(Fe2+)
%
% [2] Appendix Eq. (A2), pp. 156–158:
%
%   logK(Fe2+Fe3+) =
%       3*log10(a_Crdp) - 3*log10(a_Sud) - log10(a_Afch)
%
% [3] Appendix Eq. (A3), pp. 156–158:
%
%   T(degreeC) =
%       1/(0.00293 - 5.13e-4*logK + 3.904e-5*logK^2) - 273.15
%
% [4] Walshe-type redox expression:
%
%   logK1 = -1.926 + 6075.8/T(K)
%   log fO2 = 4*(log10(a6) - log10(a3) - logK1)
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Inoue2018(rawdata_struct, P_kbar)
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
    error('Inoue2018 requires (rawdata_struct, P_kbar).');
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

pressureCautionIssued = false;
redoxCautionIssued = false;

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
        disp([char(string(selectedCode_chl)) ': T_Inoue2018 = ' ...
            num2str(row.T_deg, '%.2f') ' degreeC']);
        disp(['Fe3+/SigmaFe = ' ...
            num2str(row.Fe3_ratio_est, '%.4f') ...
            ', logfO2 = ' num2str(row.logfO2, '%.3f')]);
    else
        disp([char(string(selectedCode_chl)) ': T_Inoue2018 = ' ...
            num2str(row.T_deg(1), '%.2f') ' to ' ...
            num2str(row.T_deg(end), '%.2f') ' degreeC']);
        disp(['Fe3+/SigmaFe = ' ...
            num2str(row.Fe3_ratio_est(1), '%.4f') ...
            ', logfO2 = ' num2str(row.logfO2(1), '%.3f')]);
    end

    if ~pressureCautionIssued
        fprintf(2, ...
            ['WARNING: Inoue et al. (2018) do not define a quantitative ' ...
             'pressure calibration range for the temperature equation. ' ...
             'Pressure is stored only for traceability; input range = ' ...
             '%.4g–%.4g kbar. The <200 bar condition discussed for ' ...
             'Noboribetsu is site-specific and is not used as a universal ' ...
             'calibration limit.\n'], ...
            min(P_kbar), max(P_kbar));
        pressureCautionIssued = true;
    end

    if ~redoxCautionIssued
        fprintf(2, ...
            ['WARNING: Inoue et al. (2018) state that quantitative log fO2 ' ...
             'estimation requires additional work. The calculated log fO2 ' ...
             'assumes ideal mixing, a_quartz = 1, and a_H2O = 1, and should ' ...
             'be interpreted primarily for relative or qualitative redox ' ...
             'comparisons.\n']);
        redoxCautionIssued = true;
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
             'this implementation of Inoue et al. (2018). This interval is ' ...
             'not a formally stated experimental calibration range. ' ...
             '%d of %d finite point(s) are outside; calculated finite ' ...
             'range = %.4g–%.4g degreeC for %s.\n'], ...
            sum(outsidePracticalRange), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Inoue et al. (2018) input(s) ' ...
             'for %s: %s.\n' ...
             '         NaN values were retained and were not replaced by ' ...
             'zero. Dependent outputs may therefore be NaN.\n'], ...
            char(string(selectedCode_chl)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    if any(~row.iter_solution_reliable)
        fprintf(2, ...
            ['WARNING: Fe3+/SigmaFe matching did not produce a reliable ' ...
             'solution for %s. The residual between Eq. (6) and Appendix ' ...
             'Eq. (A3) exceeds 0.5 degreeC or the calculation was not ' ...
             'evaluable. Derived Fe3 and log fO2 values may be NaN or ' ...
             'unreliable.\n'], ...
            char(string(selectedCode_chl)));
    end

    finiteFe3Ratio = isfinite(row.Fe3_ratio_est);
    highFe3Ratio = finiteFe3Ratio & row.Fe3_ratio_est > 0.60;
    if any(highFe3Ratio)
        fprintf(2, ...
            ['WARNING: Estimated Fe3+/SigmaFe exceeds 0.60 for %s. ' ...
             'Inoue et al. (2018, pp. 153–154) note that ferric-iron ' ...
             'capacity in trioctahedral Chlorite is structurally constrained ' ...
             'and remains controversial, particularly for Fe-poor ' ...
             'compositions.\n'], ...
            char(string(selectedCode_chl)));
    end

    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite Inoue et al. (2018) temperature values ' ...
             'were calculated for %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_chl)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    invalidFe3Ratio = ~isfinite(row.Fe3_ratio_est);
    if any(invalidFe3Ratio)
        fprintf(2, ...
            ['WARNING: Non-finite Fe3+/SigmaFe values were estimated for ' ...
             '%s. These values remain in the output table, and dependent ' ...
             'redox outputs may be NaN.\n'], ...
            char(string(selectedCode_chl)));
    end

    invalidLogfO2 = ~isfinite(row.logfO2);
    if any(invalidLogfO2)
        fprintf(2, ...
            ['WARNING: Non-finite log fO2 values were calculated for %s. ' ...
             'The values remain in the output table, and calculation has ' ...
             'not been stopped.\n'], ...
            char(string(selectedCode_chl)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Inoue2018', ...
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
% Return names of existing calculation or auxiliary variables containing NaN.

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
    error(['Inoue2018: cation variables must be numeric. ' ...
           'Non-numeric variable(s): ' ...
           char(strjoin(displayNames(nonNumericMask), ', ')) '.']);
end

if any(negativeMask)
    error(['Inoue2018: cation values must be greater than or equal to ' ...
           'zero. Negative value(s) were found in: ' ...
           char(strjoin(displayNames(negativeMask), ', ')) '.']);
end

if any(infiniteMask)
    error(['Inoue2018: infinite cation value(s) are not allowed in: ' ...
           char(strjoin(displayNames(infiniteMask), ', ')) '.']);
end

end

function row = calcTemp(data_chl, P_kbar)
% calcTemp
% Compute one composition-derived result and repeat it for each supplied
% pressure value.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

chl = prepareChloriteRow(data_chl);

site_fe2 = calcChloriteSites(chl, 0.0);
act_fe2 = calcIdealActivities(site_fe2);

logK_Fe2_scalar = calcLogK(act_fe2);
T_deg_eq6_scalar = calcTFromEq6(logK_Fe2_scalar);
T_K_eq6_scalar = T_deg_eq6_scalar + 273.15;

[Fe3_ratio_est_scalar, T_deg_A3_scalar, ...
    logK_Fe2Fe3_scalar, residual_deg_scalar, iterFlag_scalar] = ...
    estimateFe3Ratio(chl, T_deg_eq6_scalar);

T_K_A3_scalar = T_deg_A3_scalar + 273.15;

site_final = calcChloriteSites(chl, Fe3_ratio_est_scalar);
act_final = calcIdealActivities(site_final);

if isfinite(T_K_A3_scalar) && T_K_A3_scalar > 0
    logK1_scalar = -1.926 + 6075.8 ./ T_K_A3_scalar;
else
    logK1_scalar = NaN;
end

if isfinite(act_final.a6) && act_final.a6 > 0 && ...
        isfinite(act_final.a3) && act_final.a3 > 0 && ...
        isfinite(logK1_scalar)
    logfO2_scalar = 4 .* ...
        (log10(act_final.a6) - log10(act_final.a3) - logK1_scalar);
else
    logfO2_scalar = NaN;
end

Mg_number_scalar = chl.Mg ./ (chl.Mg + site_final.Fe2);
Fe_ratio_scalar = chl.FeT ./ (chl.FeT + chl.Mg + chl.Mn);

is_tetra_reasonable_scalar = ...
    isfinite(chl.Si) && chl.Si >= 2.0 && chl.Si <= 4.0;
is_AlIV_reasonable_scalar = ...
    isfinite(site_final.Al_IV) && ...
    site_final.Al_IV >= 0.0 && site_final.Al_IV <= 2.0;
is_oct_sum_reasonable_scalar = ...
    isfinite(site_final.Sum_oct) && ...
    site_final.Sum_oct >= 5.0 && site_final.Sum_oct <= 6.2;
is_vacancy_reasonable_scalar = ...
    isfinite(site_final.Vacancy) && ...
    site_final.Vacancy >= 0.0 && site_final.Vacancy <= 1.0;

has_positive_aAfch_scalar = ...
    isfinite(act_final.aAfch) && act_final.aAfch > 0;
has_positive_acrdp_scalar = ...
    isfinite(act_final.acrdp) && act_final.acrdp > 0;
has_positive_asud_scalar = ...
    isfinite(act_final.asud) && act_final.asud > 0;
has_positive_a3_scalar = ...
    isfinite(act_final.a3) && act_final.a3 > 0;
has_positive_a6_scalar = ...
    isfinite(act_final.a6) && act_final.a6 > 0;

is_in_practical_lowT_range_scalar = ...
    isfinite(T_deg_eq6_scalar) && ...
    T_deg_eq6_scalar >= 100 && T_deg_eq6_scalar <= 350;

recommended_by_Inoue2018 = ...
    repmat(is_tetra_reasonable_scalar, nP, 1) & ...
    repmat(is_AlIV_reasonable_scalar, nP, 1) & ...
    repmat(is_oct_sum_reasonable_scalar, nP, 1) & ...
    repmat(is_vacancy_reasonable_scalar, nP, 1) & ...
    repmat(has_positive_aAfch_scalar, nP, 1) & ...
    repmat(has_positive_acrdp_scalar, nP, 1) & ...
    repmat(has_positive_asud_scalar, nP, 1) & ...
    repmat(has_positive_a3_scalar, nP, 1) & ...
    repmat(has_positive_a6_scalar, nP, 1) & ...
    repmat(iterFlag_scalar, nP, 1) & ...
    repmat(is_in_practical_lowT_range_scalar, nP, 1);

row = table();

row.P_kbar = P_kbar;

row.Si = repmat(chl.Si, nP, 1);
row.Al = repmat(chl.Al, nP, 1);
row.Fe_total = repmat(chl.FeT, nP, 1);
row.Fe2 = repmat(site_final.Fe2, nP, 1);
row.Fe3 = repmat(site_final.Fe3, nP, 1);
row.Mg = repmat(chl.Mg, nP, 1);
row.Mn = repmat(chl.Mn, nP, 1);
row.Ti = repmat(chl.Ti, nP, 1);
row.Ca = repmat(chl.Ca, nP, 1);
row.Na = repmat(chl.Na, nP, 1);
row.K = repmat(chl.K, nP, 1);
row.Cr = repmat(chl.Cr, nP, 1);

row.Al_IV = repmat(site_final.Al_IV, nP, 1);
row.Al_VI = repmat(site_final.Al_VI, nP, 1);
row.Sum_oct = repmat(site_final.Sum_oct, nP, 1);
row.Vacancy = repmat(site_final.Vacancy, nP, 1);

row.XMg_oct = repmat(site_final.XMg_oct, nP, 1);
row.XAl_oct = repmat(site_final.XAl_oct, nP, 1);
row.XFe2_oct = repmat(site_final.XFe2_oct, nP, 1);
row.XFe3_oct = repmat(site_final.XFe3_oct, nP, 1);
row.Xvac_oct = repmat(site_final.Xvac_oct, nP, 1);
row.XSi_tet = repmat(site_final.XSi_tet, nP, 1);
row.XAl_tet = repmat(site_final.XAl_tet, nP, 1);

row.aAfch = repmat(act_final.aAfch, nP, 1);
row.acrdp = repmat(act_final.acrdp, nP, 1);
row.asud = repmat(act_final.asud, nP, 1);
row.a3 = repmat(act_final.a3, nP, 1);
row.a6 = repmat(act_final.a6, nP, 1);

row.logK_Fe2 = repmat(logK_Fe2_scalar, nP, 1);
row.logK_Fe2Fe3 = repmat(logK_Fe2Fe3_scalar, nP, 1);
row.logK1 = repmat(logK1_scalar, nP, 1);

row.T_deg = repmat(T_deg_eq6_scalar, nP, 1);
row.T_K = repmat(T_K_eq6_scalar, nP, 1);
row.T_deg_A3 = repmat(T_deg_A3_scalar, nP, 1);
row.T_K_A3 = repmat(T_K_A3_scalar, nP, 1);
row.deltaT_Eq6_minus_A3 = ...
    repmat(T_deg_eq6_scalar - T_deg_A3_scalar, nP, 1);

row.Fe3_ratio_est = repmat(Fe3_ratio_est_scalar, nP, 1);
row.Fe_ratio_Fe_over_FeMgMn = repmat(Fe_ratio_scalar, nP, 1);
row.Mg_number = repmat(Mg_number_scalar, nP, 1);
row.logfO2 = repmat(logfO2_scalar, nP, 1);

row.iter_residual_deg = repmat(residual_deg_scalar, nP, 1);
row.iter_solution_reliable = repmat(iterFlag_scalar, nP, 1);

row.is_tetra_reasonable = repmat(is_tetra_reasonable_scalar, nP, 1);
row.is_AlIV_reasonable = repmat(is_AlIV_reasonable_scalar, nP, 1);
row.is_oct_sum_reasonable = repmat(is_oct_sum_reasonable_scalar, nP, 1);
row.is_vacancy_reasonable = ...
    repmat(is_vacancy_reasonable_scalar, nP, 1);
row.has_positive_aAfch = repmat(has_positive_aAfch_scalar, nP, 1);
row.has_positive_acrdp = repmat(has_positive_acrdp_scalar, nP, 1);
row.has_positive_asud = repmat(has_positive_asud_scalar, nP, 1);
row.has_positive_a3 = repmat(has_positive_a3_scalar, nP, 1);
row.has_positive_a6 = repmat(has_positive_a6_scalar, nP, 1);
row.is_in_practical_lowT_range = ...
    repmat(is_in_practical_lowT_range_scalar, nP, 1);
row.pressure_calibration_defined = false(nP, 1);
row.recommended_by_Inoue2018 = recommended_by_Inoue2018;

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

chl.Ti = getVarOrZeroIfMissing(data_chl, 'Ti_cation_apfu');
chl.Ca = getVarOrZeroIfMissing(data_chl, 'Ca_cation_apfu');
chl.Na = getVarOrZeroIfMissing(data_chl, 'Na_cation_apfu');
chl.K  = getVarOrZeroIfMissing(data_chl, 'K_cation_apfu');
chl.Cr = getVarOrZeroIfMissing(data_chl, 'Cr_cation_apfu');

end

function site = calcChloriteSites(chl, Fe3_ratio)
% calcChloriteSites
% Practical cation-apfu site allocation. Missing or invalid derived
% quantities are returned as NaN rather than stopping the complete workflow.

site = makeNaNSiteStruct();

if ~isnumeric(Fe3_ratio) || ~isscalar(Fe3_ratio) || ...
        (~isnan(Fe3_ratio) && ...
         (~isfinite(Fe3_ratio) || Fe3_ratio < 0 || Fe3_ratio > 1))
    return;
end

if isnan(Fe3_ratio) || isnan(chl.Si) || isnan(chl.Al) || ...
        isnan(chl.FeT) || isnan(chl.Mg)
    return;
end

site.Al_IV = min(chl.Al, max(0, 4 - chl.Si));
site.Al_VI = chl.Al - site.Al_IV;

if ~isfinite(site.Al_VI) || site.Al_VI < -1e-10
    site = makeNaNSiteStruct();
    return;
end

site.Al_VI = max(0, site.Al_VI);
site.Fe3 = chl.FeT .* Fe3_ratio;
site.Fe2 = chl.FeT - site.Fe3;

if site.Fe2 < -1e-10 || site.Fe3 < -1e-10
    site = makeNaNSiteStruct();
    return;
end

site.Fe2 = max(0, site.Fe2);
site.Fe3 = max(0, site.Fe3);

site.Sum_oct = site.Al_VI + site.Fe2 + site.Fe3 + chl.Mg;
site.Vacancy = 6 - site.Sum_oct;

if isfinite(site.Vacancy) && ...
        site.Vacancy > -1e-8 && site.Vacancy < 0
    site.Vacancy = 0;
end

site.XMg_oct = chl.Mg ./ 6;
site.XAl_oct = site.Al_VI ./ 6;
site.XFe2_oct = site.Fe2 ./ 6;
site.XFe3_oct = site.Fe3 ./ 6;
site.Xvac_oct = site.Vacancy ./ 6;
site.XSi_tet = (chl.Si - 2) ./ 2;
site.XAl_tet = site.Al_IV ./ 2;

end

function site = makeNaNSiteStruct()
% makeNaNSiteStruct
% Return a final-site structure whose fields are all NaN.

site = struct( ...
    'Al_IV', NaN, ...
    'Al_VI', NaN, ...
    'Fe2', NaN, ...
    'Fe3', NaN, ...
    'Sum_oct', NaN, ...
    'Vacancy', NaN, ...
    'XMg_oct', NaN, ...
    'XAl_oct', NaN, ...
    'XFe2_oct', NaN, ...
    'XFe3_oct', NaN, ...
    'Xvac_oct', NaN, ...
    'XSi_tet', NaN, ...
    'XAl_tet', NaN);

end

function act = calcIdealActivities(site)
% calcIdealActivities
% Calculate ideal activities. Non-evaluable terms are retained as NaN.

act = struct( ...
    'aAfch', NaN, ...
    'acrdp', NaN, ...
    'asud', NaN, ...
    'a3', NaN, ...
    'a6', NaN);

baseValid = ...
    isfinite(site.XMg_oct) && site.XMg_oct > 0 && ...
    isfinite(site.XAl_oct) && site.XAl_oct > 0 && ...
    isfinite(site.XSi_tet) && site.XSi_tet > 0 && ...
    isfinite(site.XAl_tet) && site.XAl_tet > 0;

if baseValid
    act.aAfch = ...
        (site.XMg_oct.^6) .* (site.XSi_tet.^2);

    act.acrdp = ...
        45.563 .* (site.XMg_oct.^4) .* ...
        (site.XAl_oct.^2) .* (site.XAl_tet.^2);

    if isfinite(site.Xvac_oct) && site.Xvac_oct > 0
        act.asud = ...
            1728 .* (site.XMg_oct.^2) .* ...
            (site.XAl_oct.^3) .* site.Xvac_oct .* ...
            site.XSi_tet .* site.XAl_tet;
    end

    if isfinite(site.XFe2_oct) && site.XFe2_oct > 0
        act.a3 = ...
            59.720 .* (site.XFe2_oct.^5) .* ...
            site.XAl_oct .* site.XSi_tet .* site.XAl_tet;
    end
end

if isfinite(site.Fe3) && site.Fe3 > 0
    F = 28 ./ (28 + site.Fe3);
    act.a6 = site.Fe3 .* (1 - F);
end

end

function logK = calcLogK(act)
% calcLogK
% Calculate the Appendix equilibrium variable. Return NaN if any required
% activity is non-positive or non-finite.

if isfinite(act.acrdp) && act.acrdp > 0 && ...
        isfinite(act.asud) && act.asud > 0 && ...
        isfinite(act.aAfch) && act.aAfch > 0
    logK = 3 .* log10(act.acrdp) ...
         - 3 .* log10(act.asud) ...
         - log10(act.aAfch);
else
    logK = NaN;
end

end

function T_deg = calcTFromEq6(logK_Fe2)
% calcTFromEq6
% Inoue et al. (2018), Eq. (6). Return NaN for non-finite input.

if ~isfinite(logK_Fe2)
    T_deg = NaN;
    return;
end

x = logK_Fe2;

T_deg = 63.83 ...
      + 50.41 .* x ...
      + 2.617 .* x.^2 ...
      + 2.846 .* x.^3 ...
      - 1.097 .* x.^4 ...
      + 0.09285 .* x.^5;

if ~isfinite(T_deg) || ~isreal(T_deg)
    T_deg = NaN;
end

end

function T_deg = calcTFromA3(logK_Fe2Fe3)
% calcTFromA3
% Appendix Eq. (A3). Return NaN if the denominator is invalid.

if ~isfinite(logK_Fe2Fe3)
    T_deg = NaN;
    return;
end

denominator = 0.00293 ...
            - 5.13e-4 .* logK_Fe2Fe3 ...
            + 3.904e-5 .* logK_Fe2Fe3.^2;

if ~isfinite(denominator) || denominator <= 0
    T_deg = NaN;
    return;
end

T_deg = 1 ./ denominator - 273.15;

if ~isfinite(T_deg) || ~isreal(T_deg)
    T_deg = NaN;
end

end

function [Fe3_ratio_best, T_deg_best, logK_best, residual_deg, iterFlag] = ...
    estimateFe3Ratio(chl, T_target_deg)
% estimateFe3Ratio
% Estimate Fe3+/SigmaFe by minimizing the mismatch between Eq. (6) and
% Appendix Eq. (A3). Failures return NaN without stopping the main workflow.

Fe3_ratio_best = NaN;
T_deg_best = NaN;
logK_best = NaN;
residual_deg = NaN;
iterFlag = false;

if ~isfinite(T_target_deg) || ~isfinite(chl.FeT) || chl.FeT <= 0
    return;
end

upperBound = 0.95;
objective = @(ratio) localObjective(ratio, chl, T_target_deg);
options = optimset('TolX', 1e-6, 'Display', 'off');

try
    ratioCandidate = fminbnd(objective, 0.0, upperBound, options);
catch
    return;
end

[T_candidate, logK_candidate] = localTA3(ratioCandidate, chl);

if ~isfinite(T_candidate) || ~isfinite(logK_candidate)
    return;
end

Fe3_ratio_best = ratioCandidate;
T_deg_best = T_candidate;
logK_best = logK_candidate;
residual_deg = abs(T_deg_best - T_target_deg);
iterFlag = isfinite(residual_deg) && residual_deg <= 0.5;

end

function misfit = localObjective(Fe3_ratio, chl, T_target_deg)
% localObjective
% Objective function used by fminbnd. Invalid trial values receive a large
% finite penalty rather than throwing an error.

[T_deg_A3, ~] = localTA3(Fe3_ratio, chl);

if isfinite(T_deg_A3)
    misfit = abs(T_deg_A3 - T_target_deg);
else
    misfit = realmax('double');
end

end

function [T_deg_A3, logK_A3] = localTA3(Fe3_ratio, chl)
% localTA3
% Calculate Appendix Eq. (A3) temperature for one trial Fe3 ratio.

site = calcChloriteSites(chl, Fe3_ratio);
act = calcIdealActivities(site);
logK_A3 = calcLogK(act);
T_deg_A3 = calcTFromA3(logK_A3);

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
% Retrieve an optional numeric scalar. Zero is assigned only when the column
% is absent; an existing NaN is retained.

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
