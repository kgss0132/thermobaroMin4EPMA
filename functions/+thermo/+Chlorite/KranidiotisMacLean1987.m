function results = KranidiotisMacLean1987(rawdata_struct, P_kbar)
% functions/+thermo/+Chlorite/KranidiotisMacLean1987.m
% Tested with MATLAB R2024b
%
% Empirical Fe-Mg-corrected single-Chlorite thermometer
% Kranidiotis, P. and MacLean, W.H. (1987)
% Economic Geology, 82, 1898–1911
% DOI: https://doi.org/10.2113/gsecongeo.82.7.1898
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Chlorite analysis at a time and
% calculates temperature using the Fe-Mg-corrected Chlorite geothermometer
% described by Kranidiotis and MacLean (1987).
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
% FIELD-CALIBRATION RANGE AND APPLICATION NOTES
%
% Kranidiotis and MacLean (1987) did not establish an independent
% experimental calibration. They adopted the field-calibrated relationship
% of Cathelineau and Nieva (1985) and introduced an Fe/(Fe+Mg) correction
% for Chlorites on the Al-saturated boundary of the Chlorite solid-solution
% field. The temperature relation and Fe-Mg correction are discussed on
% p. 1909 and illustrated in Fig. 11 (p. 1909).
%
% The directly illustrated low-temperature field-calibration interval is
% approximately:
%
%   Temperature : 150–300 degreeC
%
% This is a field-based reference interval rather than a controlled
% experimental calibration range. Temperatures outside this interval are
% retained but reported by non-stopping fprintf warnings.
%
% PRESSURE LIMITATIONS:
% The effect of pressure on the position and shape of the Chlorite
% solid-solution field is stated to be unknown (p. 1909). The natural systems
% directly discussed in connection with the thermometer formed at low
% pressure:
%
%   Phelps Dodge : approximately 100–300 bar = 0.1–0.3 kbar
%   Los Azufres  : up to approximately 700 bar = 0.7 kbar
%
% High-pressure experimental data obtained at 2–6 kbar did not follow the
% field-calibration regression (p. 1909). The authors therefore recommend
% limiting use to Chlorites that clearly formed in low-pressure environments
% and were not recrystallized during metamorphism (p. 1909).
%
% This implementation uses 0.1–0.7 kbar only as the directly referenced
% natural low-pressure interval. It is not a formal experimental pressure
% calibration range. Input pressures outside 0.1–0.7 kbar are retained and
% reported by fprintf. Pressures of 2 kbar or higher receive an additional
% warning because the cited high-pressure experiments did not agree with the
% field regression.
%
% COMPOSITIONAL AND PETROGRAPHIC LIMITATIONS:
%
%   - The thermometer is intended primarily for Al-saturated Chlorites on
%     the Al-rich boundary of the solid-solution field. Such Chlorites occur
%     with an Al-rich phase such as sericite, albite, epidote, or clay
%     minerals (pp. 1904–1906 and 1909).
%   - Al-undersaturated Chlorites associated with talc or stilpnomelane were
%     recognized as a possible future geothermometer, but that boundary was
%     not calibrated by the equation implemented here (pp. 1909–1910).
%   - The analysed Chlorite must represent the hydrothermal event of
%     interest and must not have recrystallized or re-equilibrated during
%     later metamorphism (p. 1909).
%   - Coexisting Fe-Mg-Al minerals such as biotite, hornblende, talc, and
%     stilpnomelane can alter the relevant phase relations and require
%     careful petrographic interpretation (pp. 1905–1906 and 1910).
%   - Fe/(Fe+Mg) is influenced not only by temperature but also by fluid
%     composition, mixing of Mg-rich seawater and Fe-rich hydrothermal
%     fluid, water/rock ratio, oxygen fugacity, sulfur fugacity, and pH
%     (pp. 1908–1909).
%   - The Phelps Dodge Chlorites span Fe/(Fe+Mg) = 0.183–0.642
%     (Table 1 and discussion, pp. 1901–1906). Values outside this interval
%     are compositional extrapolations. Values above approximately 0.60
%     require additional caution because Fe3+ may contribute to charge
%     balance (p. 1906).
%   - Ferric iron was not measured directly by EPMA. The authors estimated
%     relatively small Fe3+ contents and inferred reducing conditions below
%     QFM (p. 1901). This implementation uses EPMA total Fe in the published
%     Fe/(Fe+Mg) correction and does not apply an independent Fe3+ correction.
%   - Application to Phelps Dodge gave approximately 290 degreeC, which the
%     authors considered reasonable but approximately 50–75 degreeC lower
%     than their geological model. They explicitly state that additional
%     calibration is needed (p. 1909).
%
% STRUCTURAL-FORMULA NORMALIZATION (IMPORTANT):
% The published equations use a full Chlorite formula normalized to eight
% tetrahedral cations (Si + Al_IV = 8). The thermo.Chlorite input used here
% is normalized to a 14-oxygen half-formula basis, for which
% Si + Al_IV = 4. The published equation must therefore be converted before
% it is applied to these inputs.
%
% Published full-formula equation (p. 1909):
%
%   Al_IV_corrected_full = Al_IV_full + 0.70*Fe/(Fe+Mg)
%   T(degreeC) = 106*Al_IV_corrected_full + 18
%
% Equivalent 14-oxygen half-formula equation implemented here:
%
%   Al_IV_corrected_14O = Al_IV_14O + 0.35*Fe/(Fe+Mg)
%   T(degreeC) = 212*Al_IV_corrected_14O + 18
%
% These two forms are algebraically equivalent because
% Al_IV_full = 2*Al_IV_14O. Applying the coefficients 106 and 0.70 directly
% to a 14-oxygen half-formula would not reproduce the published equation.
%
% This implementation issues non-stopping fprintf messages when:
%   1) a finite calculated temperature is outside 150–300 degreeC,
%   2) input pressure is outside the directly referenced 0.1–0.7 kbar
%      natural low-pressure interval,
%   3) input pressure is 2 kbar or higher,
%   4) Fe/(Fe+Mg) is outside the observed 0.183–0.642 interval,
%   5) Fe/(Fe+Mg) exceeds 0.60,
%   6) a calculation input contains NaN,
%   7) site allocation is not evaluable, or
%   8) a non-finite temperature is calculated.
%
% The script cannot determine mineral assemblage or recrystallization state
% from a single EPMA row. Al saturation, phase assemblage, paragenesis, and
% preservation of the hydrothermal Chlorite generation must be verified
% independently before interpreting the result.
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
% Required variables, normalized on a 14-oxygen half-formula basis:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu      % EPMA total Fe used in Fe/(Fe+Mg)
%   Mg_cation_apfu
%
% Optional calculation variable:
%   Mn_cation_apfu      % used only in an auxiliary Fe/(Fe+Mg+Mn) index
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
% never replaced by zero when the corresponding column exists. Optional
% variables are assigned zero only when their columns are absent.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% On a 14-oxygen half-formula basis:
%
%   Al_IV_14O = min(Al_total, max(0, 4 - Si))
%
%   FeMg_ratio = Fe_total/(Fe_total + Mg)
%
%   Al_IV_corrected_14O = Al_IV_14O + 0.35*FeMg_ratio
%
%   T(degreeC) = 212*Al_IV_corrected_14O + 18
%
% Equivalent full-formula quantities retained for traceability:
%
%   Al_IV_full = 2*Al_IV_14O
%   Al_IV_corrected_full = Al_IV_full + 0.70*FeMg_ratio
%   T(degreeC) = 106*Al_IV_corrected_full + 18
%
% Pressure is not used in the equation. A scalar or vector P_kbar is accepted
% for compatibility with fixed-pressure and pressure-range workflows and is
% stored in the output table.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = KranidiotisMacLean1987(rawdata_struct, P_kbar)
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
    error('KranidiotisMacLean1987 requires (rawdata_struct, P_kbar).');
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

initialBufferCapacity = max(16, height(dataset_chl));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degC = 150;
calibrationT_max_degC = 300;
referenceP_min_kbar = 0.1;
referenceP_max_kbar = 0.7;
highPressureCaution_kbar = 2.0;
observedFeMg_min = 0.183;
observedFeMg_max = 0.642;
Fe3CautionRatio = 0.60;

pressureCautionIssued = false;
applicationCautionIssued = false;

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
        disp([char(string(selectedCode_chl)) ...
            ': T_KranidiotisMacLean1987 = ' num2str(row.T_deg) ...
            ' degreeC, AlIV_corr_14O = ' ...
            num2str(row.Al_IV_corrected_14O)]);
    else
        disp([char(string(selectedCode_chl)) ...
            ': T_KranidiotisMacLean1987 = ' ...
            num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ...
            ' degreeC, AlIV_corr_14O = ' ...
            num2str(row.Al_IV_corrected_14O(1))]);
    end

    if ~pressureCautionIssued
        outsideReferencePressure = ...
            P_kbar < referenceP_min_kbar | P_kbar > referenceP_max_kbar;

        if any(outsideReferencePressure)
            fprintf(2, ...
                ['WARNING: %d of %d input pressure point(s) are outside ' ...
                 'the directly referenced natural low-pressure interval ' ...
                 'of approximately 0.1–0.7 kbar; input range = ' ...
                 '%.4g–%.4g kbar. This interval combines the Phelps Dodge ' ...
                 'and Los Azufres conditions discussed by Kranidiotis and ' ...
                 'MacLean (1987, p. 1909) and is not a formal experimental ' ...
                 'pressure calibration range. Pressure is not used in the ' ...
                 'thermometer equation.\n'], ...
                sum(outsideReferencePressure), ...
                numel(P_kbar), ...
                min(P_kbar), ...
                max(P_kbar));
        else
            fprintf(2, ...
                ['WARNING: Kranidiotis and MacLean (1987) do not define a ' ...
                 'formal pressure calibration range, and pressure is not ' ...
                 'used in the equation. All supplied pressures fall within ' ...
                 'the directly referenced natural low-pressure interval of ' ...
                 'approximately 0.1–0.7 kbar; input range = %.4g–%.4g kbar.\n'], ...
                min(P_kbar), ...
                max(P_kbar));
        end

        highPressurePoints = P_kbar >= highPressureCaution_kbar;
        if any(highPressurePoints)
            fprintf(2, ...
                ['WARNING: %d of %d pressure point(s) are >=2 kbar. ' ...
                 'Kranidiotis and MacLean (1987, p. 1909) report that ' ...
                 'high-pressure experiments at 2–6 kbar did not follow the ' ...
                 'field-calibration regression. Results at these pressures ' ...
                 'are strong extrapolations.\n'], ...
                sum(highPressurePoints), ...
                numel(P_kbar));
        end

        pressureCautionIssued = true;
    end

    if ~applicationCautionIssued
        fprintf(2, ...
            ['WARNING: This calculation cannot verify Al saturation, ' ...
             'mineral assemblage, paragenesis, or preservation from later ' ...
             'metamorphic recrystallization. Kranidiotis and MacLean ' ...
             '(1987, pp. 1909–1910) recommend application only to clearly ' ...
             'low-pressure, non-recrystallized, Al-saturated Chlorites.\n']);
        applicationCautionIssued = true;
    end

    finiteTemperature = isfinite(row.T_deg);
    outsideTemperatureRange = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(outsideTemperatureRange)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximately ' ...
             '150–300 degreeC low-temperature field-calibration interval ' ...
             'illustrated by Kranidiotis and MacLean (1987, Fig. 11, ' ...
             'p. 1909). %d of %d finite point(s) are outside; calculated ' ...
             'finite range = %.4g–%.4g degreeC for %s. The interval is ' ...
             'field based and is not a controlled experimental range.\n'], ...
            sum(outsideTemperatureRange), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    end

    finiteFeMg = isfinite(row.Fe_ratio_Fe_over_FeMg);
    outsideObservedFeMg = finiteFeMg & ...
        (row.Fe_ratio_Fe_over_FeMg < observedFeMg_min | ...
         row.Fe_ratio_Fe_over_FeMg > observedFeMg_max);

    if any(outsideObservedFeMg)
        finiteValues = row.Fe_ratio_Fe_over_FeMg(finiteFeMg);
        fprintf(2, ...
            ['WARNING: Fe/(Fe+Mg) is outside the observed Phelps Dodge ' ...
             'Chlorite range of 0.183–0.642 reported by Kranidiotis and ' ...
             'MacLean (1987, Table 1 and pp. 1905–1906). Calculated finite ' ...
             'range = %.4g–%.4g for %s.\n'], ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    end

    highFeMg = finiteFeMg & row.Fe_ratio_Fe_over_FeMg > Fe3CautionRatio;
    if any(highFeMg)
        fprintf(2, ...
            ['WARNING: Fe/(Fe+Mg) exceeds 0.60 for %s. Kranidiotis and ' ...
             'MacLean (1987, p. 1906) indicate that Fe3+ may contribute to ' ...
             'charge balance in the more Fe-rich Chlorites. The present ' ...
             'implementation uses EPMA total Fe and does not apply an ' ...
             'independent Fe3+ correction.\n'], ...
            char(string(selectedCode_chl)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Kranidiotis and MacLean (1987) ' ...
             'calculation input(s) for %s: %s.\n' ...
             '         NaN values were retained and were not replaced by ' ...
             'zero. Dependent outputs may therefore be NaN.\n'], ...
            char(string(selectedCode_chl)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    if any(~row.is_site_allocation_evaluable)
        fprintf(2, ...
            ['WARNING: Chlorite site allocation was not evaluable for %s. ' ...
             'The supplied Si and Al values are missing or internally ' ...
             'inconsistent with the 14-oxygen half-formula allocation. ' ...
             'Dependent outputs were retained as NaN and calculation ' ...
             'continued.\n'], ...
            char(string(selectedCode_chl)));
    end

    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite Kranidiotis and MacLean (1987) ' ...
             'temperature values were calculated for %s (%d of %d points; ' ...
             'NaN: %d, Inf: %d).\n' ...
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
        'KranidiotisMacLean1987', ...
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
% Return names of existing variables used in the thermometer calculation
% that contain NaN.

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
    error(['KranidiotisMacLean1987: cation variables must be numeric. ' ...
           'Non-numeric variable(s): ' ...
           char(strjoin(displayNames(nonNumericMask), ', ')) '.']);
end

if any(negativeMask)
    error(['KranidiotisMacLean1987: cation values must be greater than ' ...
           'or equal to zero. Negative value(s) were found in: ' ...
           char(strjoin(displayNames(negativeMask), ', ')) '.']);
end

if any(infiniteMask)
    error(['KranidiotisMacLean1987: infinite cation value(s) are not ' ...
           'allowed in: ' ...
           char(strjoin(displayNames(infiniteMask), ', ')) '.']);
end

end

function row = calcTemp(data_chl, P_kbar)
% calcTemp
% Compute one composition-derived temperature and repeat all scalar outputs
% for each supplied pressure value.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

chl = prepareChloriteRow(data_chl);
site = calcChloriteSites(chl);

FeMg_ratio_scalar = chl.FeT ./ (chl.FeT + chl.Mg);
Mg_number_scalar = chl.Mg ./ (chl.Mg + chl.FeT);
Fe_ratio_Fe_over_FeMgMn_scalar = ...
    chl.FeT ./ (chl.FeT + chl.Mg + chl.Mn);

% Convert the full-formula published coefficients to the 14-oxygen
% half-formula input basis.
Al_IV_corrected_14O_scalar = ...
    site.Al_IV + 0.35 .* FeMg_ratio_scalar;
Al_IV_full_scalar = 2 .* site.Al_IV;
Al_IV_corrected_full_scalar = ...
    Al_IV_full_scalar + 0.70 .* FeMg_ratio_scalar;

% Both equations below are algebraically equivalent. The 14-oxygen form is
% used as the primary implementation; the full-formula result is retained as
% an internal consistency check.
T_deg_scalar = 212 .* Al_IV_corrected_14O_scalar + 18;
T_deg_from_full_formula_scalar = ...
    106 .* Al_IV_corrected_full_scalar + 18;
T_K_scalar = T_deg_scalar + 273.15;

T_uncorrected_deg_scalar = 212 .* site.Al_IV + 18;
T_uncorrected_K_scalar = T_uncorrected_deg_scalar + 273.15;

normalization_consistent_scalar = ...
    isfinite(T_deg_scalar) && ...
    isfinite(T_deg_from_full_formula_scalar) && ...
    abs(T_deg_scalar - T_deg_from_full_formula_scalar) <= 1e-10;

is_14O_tetra_reasonable_scalar = ...
    isfinite(chl.Si) && chl.Si >= 2.0 && chl.Si <= 4.0;

is_AlIV_reasonable_scalar = ...
    isfinite(site.Al_IV) && site.Al_IV >= 0.0 && site.Al_IV <= 2.0;

is_FeMg_ratio_reasonable_scalar = ...
    isfinite(FeMg_ratio_scalar) && ...
    FeMg_ratio_scalar >= 0.0 && FeMg_ratio_scalar <= 1.0;

is_in_reference_T_range_scalar = ...
    isfinite(T_deg_scalar) && ...
    T_deg_scalar >= 150 && T_deg_scalar <= 300;

is_in_observed_FeMg_range_scalar = ...
    isfinite(FeMg_ratio_scalar) && ...
    FeMg_ratio_scalar >= 0.183 && FeMg_ratio_scalar <= 0.642;

is_in_reference_P_range = ...
    P_kbar >= 0.1 & P_kbar <= 0.7;

% This is a numerical screening flag only. The required Al-saturated phase
% assemblage and absence of later metamorphic recrystallization cannot be
% determined from the EPMA row and must be confirmed independently.
recommended_for_KranidiotisMacLean1987 = ...
    repmat(site.is_evaluable, nP, 1) & ...
    repmat(is_14O_tetra_reasonable_scalar, nP, 1) & ...
    repmat(is_AlIV_reasonable_scalar, nP, 1) & ...
    repmat(is_FeMg_ratio_reasonable_scalar, nP, 1) & ...
    repmat(is_in_reference_T_range_scalar, nP, 1) & ...
    repmat(is_in_observed_FeMg_range_scalar, nP, 1) & ...
    is_in_reference_P_range & ...
    repmat(normalization_consistent_scalar, nP, 1);

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
row.Al_IV_14O = repmat(site.Al_IV, nP, 1);
row.Al_IV_full_formula = repmat(Al_IV_full_scalar, nP, 1);
row.Al_VI = repmat(site.Al_VI, nP, 1);
row.Sum_VI = repmat(site.Sum_VI, nP, 1);
row.VAC = repmat(site.VAC, nP, 1);

row.Mg_number = repmat(Mg_number_scalar, nP, 1);
row.Fe_ratio_Fe_over_FeMg = repmat(FeMg_ratio_scalar, nP, 1);
row.Fe_ratio_Fe_over_FeMgMn = ...
    repmat(Fe_ratio_Fe_over_FeMgMn_scalar, nP, 1);

% Preserve the legacy output name, but define it explicitly on the 14-oxygen
% basis. Additional basis-specific columns remove ambiguity.
row.Al_IV_corrected = repmat(Al_IV_corrected_14O_scalar, nP, 1);
row.Al_IV_corrected_14O = ...
    repmat(Al_IV_corrected_14O_scalar, nP, 1);
row.Al_IV_corrected_full_formula = ...
    repmat(Al_IV_corrected_full_scalar, nP, 1);

row.T_deg = repmat(T_deg_scalar, nP, 1);
row.T_K = repmat(T_K_scalar, nP, 1);
row.T_deg_from_full_formula = ...
    repmat(T_deg_from_full_formula_scalar, nP, 1);

row.T_uncorrected_deg = repmat(T_uncorrected_deg_scalar, nP, 1);
row.T_uncorrected_K = repmat(T_uncorrected_K_scalar, nP, 1);
row.deltaT_corrected_minus_uncorrected = ...
    repmat(T_deg_scalar - T_uncorrected_deg_scalar, nP, 1);

row.is_site_allocation_evaluable = repmat(site.is_evaluable, nP, 1);
row.is_14O_tetra_reasonable = ...
    repmat(is_14O_tetra_reasonable_scalar, nP, 1);
row.is_AlIV_reasonable = repmat(is_AlIV_reasonable_scalar, nP, 1);
row.is_FeMg_ratio_reasonable = ...
    repmat(is_FeMg_ratio_reasonable_scalar, nP, 1);
row.is_in_KranidiotisMacLean_T_range = ...
    repmat(is_in_reference_T_range_scalar, nP, 1);
row.is_in_KranidiotisMacLean_FeMg_range = ...
    repmat(is_in_observed_FeMg_range_scalar, nP, 1);
row.pressure_calibration_defined = false(nP, 1);
row.is_in_directly_referenced_lowP_range = is_in_reference_P_range;
row.is_high_pressure_extrapolation = P_kbar >= 2.0;
row.normalization_conversion_consistent = ...
    repmat(normalization_consistent_scalar, nP, 1);
row.requires_Al_saturation_confirmation = true(nP, 1);
row.requires_no_metamorphic_recrystallization_confirmation = true(nP, 1);
row.recommended_for_KranidiotisMacLean1987 = ...
    recommended_for_KranidiotisMacLean1987;

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
% Calculate approximate site quantities on a 14-oxygen half-formula basis.
% Missing or internally inconsistent derived values are retained as NaN.

site = struct( ...
    'Al_IV', NaN, ...
    'Al_VI', NaN, ...
    'Sum_VI', NaN, ...
    'VAC', NaN, ...
    'is_evaluable', false);

if isnan(chl.Si) || isnan(chl.Al)
    return;
end

site.Al_IV = min(chl.Al, max(0, 4 - chl.Si));
site.Al_VI = chl.Al - site.Al_IV;

if ~isfinite(site.Al_VI) || site.Al_VI < -1e-10
    site.Al_IV = NaN;
    site.Al_VI = NaN;
    return;
end

site.Al_VI = max(0, site.Al_VI);
site.Sum_VI = site.Al_VI + chl.FeT + chl.Mg + chl.Mn;
site.VAC = 6 - site.Sum_VI;
site.is_evaluable = ...
    isfinite(site.Al_IV) && isfinite(site.Al_VI) && ...
    isfinite(site.Sum_VI) && isfinite(site.VAC);

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
