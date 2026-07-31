function results = NimisTaylor2000(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/NimisTaylor2000.m
% Tested with MATLAB R2024b
%
% Single-clinopyroxene thermometer for Cr-diopside from garnet peridotite
% Nimis, P. and Taylor, W.R. (2000), equation (17)
% Contributions to Mineralogy and Petrology, 139, 541-554
% DOI: https://doi.org/10.1007/s004100000156
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Clinopyroxene analysis and
% calculates temperature using the enstatite-in-clinopyroxene thermometer
% of Nimis and Taylor (2000), equation (17).
%
% The function accepts either a scalar pressure or a pressure vector and is
% therefore compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP. For each selected Cpx analysis, the output table
% contains one row per input pressure value.
%
% The function is designed for repeated calculations. Each selected-Cpx
% result is stored in a fixed-size preallocated cell buffer, and all result
% blocks are concatenated only once after the interactive loop finishes.
% No result array is enlarged inside the loop.
%
% This is a mineral-only thermometer and does not use a Liquid dataset.
% Therefore cationTotal_liq, F, and Cl handling is not applicable.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Nimis and Taylor (2000) developed this thermometer primarily for
% Cr-diopside derived from garnet peridotite. The principal assumption is
% that the analyzed Cpx originally equilibrated with orthopyroxene and
% garnet. Application to igneous phenocrysts, eclogitic Cpx, metasomatic Cpx,
% or Cpx of uncertain paragenesis may yield geologically meaningless
% apparent temperatures.
%
% The experimental compilation summarized in the abstract spans:
%
%   Pressure    : approximately 0-60 kbar
%   Temperature : approximately 850-1500 degreeC
%
% For natural-like garnet-bearing lherzolite systems, the principal range
% described on p. 544 is approximately:
%
%   Pressure    : 20-60 kbar
%   Temperature : 900-1400 degreeC
%
% Some simple-system tests shown in Figure 3 extend to approximately
% 75 kbar and 1600 degreeC. These extensions should not be interpreted as a
% general calibration range for compositionally complex natural Cpx.
%
% Equation (17) is presented on p. 546. It reproduces the full calibration
% dataset (n = 131) with a standard deviation of approximately 30 degreeC
% and correlation coefficient R = 0.981 (Table 2 and Figure 3, p. 547).
% Individual experimental subsets show standard deviations of approximately
% 25-32 degreeC; simple-system tests in Figure 3 are approximately
% 36-40 degreeC.
%
% Important application cautions stated or demonstrated in the paper:
%
%   1) The selected Cpx should have equilibrated with Opx and garnet.
%      Apparent P-T values from Cpx that did not equilibrate with these
%      phases may be meaningless (pp. 545-546 and 550-551).
%
%   2) Low-Al Cpx is not necessarily derived from garnet peridotite. Strongly
%      depleted or metasomatized spinel peridotites may also contain low-Al
%      Cr-diopside. Mineral paragenesis, texture, and where possible trace
%      elements should be used to assess origin (pp. 550-551).
%
%   3) The Fe, Ti, Al, Cr, Na, and K correction terms in equation (17) are
%      empirical. Application to compositions substantially different from
%      the calibration dataset may produce biased estimates (p. 546).
%
%   4) The formulation treats Fe_cation_apfu as total Fe and does not apply
%      a separate Fe3+ correction. Neglecting significant Fe3+ may lead to
%      underestimation of temperature; the authors considered this effect
%      generally small for natural garnet-peridotite Cpx (p. 549).
%
%   5) Natural-sample screening in the paper required Cpx cations normalized
%      to 6 oxygens and cation sums greater than 1.990 for both the T site
%      and the combined M1+M2 sites. Analyses failing these criteria were
%      rejected (p. 548). This implementation reports non-stopping warnings
%      when calculated site sums are below 1.990.
%
%   6) The authors restricted their natural Cr-diopside application dataset
%      to Cpx with Cr2O3 < 5 wt.% because the simplified Cr-activity model is
%      probably inaccurate for very kosmochlor-rich compositions (p. 548).
%      This limit cannot be checked directly from cation-apfu inputs alone.
%
% This implementation uses two levels of non-stopping fprintf warnings:
%
%   Principal natural-like range:
%     Pressure    : 20-60 kbar
%     Temperature : 900-1400 degreeC
%
%   Broader experimental compilation:
%     Pressure    : 0-60 kbar
%     Temperature : 850-1500 degreeC
%
% Values outside the principal range are extrapolative for natural-like
% garnet-lherzolite Cpx. Values outside the broader compilation represent
% stronger extrapolation.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. Cation variables should be normalized
% to a 6-oxygen basis.
%
% Required variables:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu         % total Fe used directly in equation (17)
%   Mg_cation_apfu
%   Ca_cation_apfu
%
% Optional variables:
%   Na_cation_apfu         % absent column is treated as 0
%   K_cation_apfu          % absent column is treated as 0
%   Ti_cation_apfu         % absent column is treated as 0
%   Cr_cation_apfu         % absent column is treated as 0
%   Mn_cation_apfu         % absent column is treated as 0; stored only
%   Fe3_cation_apfu        % absent column is treated as 0; not used in T
%
% If an optional column is absent, the documented zero assumption above is
% used. If a column is present and its selected value is NaN, the NaN is
% retained; it is never replaced by zero. NaN values in equation (17) inputs
% propagate to the temperature result and are reported using fprintf.
%
% All finite mineral-composition inputs must be greater than or equal to
% zero. Finite negative values and Inf are rejected. Zero is allowed as a
% raw input, but an invalid logarithm, denominator, or activity term produces
% a retained NaN temperature and a non-stopping warning.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Nimis and Taylor (2000), equation (17):
%
%   T(K) = (23166 + 39.28*P_kbar) / denominator
%
%   denominator = 13.25
%                 + 15.35*Ti
%                 + 4.50*Fe
%                 - 1.55*(Al + Cr - Na - K)
%                 + [ln(aEn_cpx)]^2
%
%   aEn_cpx = (1 - Ca - Na - K)
%              * [1 - (Al + Cr + Na + K)/2]
%
% where all compositional variables are Cpx cations per formula unit on a
% 6-oxygen basis, Fe is total Fe from Fe_cation_apfu, P is in kbar, and
% natural logarithms are used.
%
% The two factors forming aEn_cpx and aEn_cpx itself must be finite and
% strictly greater than zero. The thermometer denominator must also be
% finite and strictly greater than zero. Domain failures return NaN rather
% than stopping the complete calculation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = NimisTaylor2000(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing the Cpx table described above
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Cpx analysis. NaN and Inf results are retained.
%

%% Input validation
if nargin < 2
    error('NimisTaylor2000 requires (rawdata_struct, P_kbar).');
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

if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

dataset_cpx = rawdata_struct.Cpx;

requiredVariables = {'Si_cation_apfu', 'Al_cation_apfu', ...
    'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu'};

validateRequiredVariables(dataset_cpx, requiredVariables);

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% A fixed-size buffer avoids all loop-dependent array resizing. The limit is
% intentionally much larger than a normal interactive session. If the limit
% is reached, the completed results are returned without enlarging the array.
disp('=== Step 2: Preparing output container ===');

maxResultBlocks = max(1024, height(dataset_cpx));
resultBlocks = cell(maxResultBlocks, 1);
nResultBlocks = 0;

% Principal natural-like range from the garnet-lherzolite calibration.
principalP_min_kbar = 20;
principalP_max_kbar = 60;
principalT_min_degC = 900;
principalT_max_degC = 1400;

% Broader experimental compilation summarized in the paper.
broadP_min_kbar = 0;
broadP_max_kbar = 60;
broadT_min_degC = 850;
broadT_max_degC = 1500;

pressureOutsidePrincipal = P_kbar < principalP_min_kbar | ...
    P_kbar > principalP_max_kbar;
pressureOutsideBroad = P_kbar < broadP_min_kbar | ...
    P_kbar > broadP_max_kbar;
pressureWarningsIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive selection loop and calculation
disp('=== Step 3: Selecting a data code from the list (Cpx) ===');

while true
    % Prevent any result-buffer size change inside the loop.
    if nResultBlocks >= maxResultBlocks
        fprintf(2, ...
            ['WARNING: The fixed result-buffer limit of %d selections was ' ...
             'reached. The completed calculations will be returned without ' ...
             'enlarging the result array.\n'], maxResultBlocks);
        break;
    end

    % ----- Cpx selection -----
    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Cpx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIdx_cpx);
    disp(['Cpx selected: ' char(string(selectedCode_cpx))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    % Report NaN separately for equation inputs and auxiliary variables.
    [nanEquationInputs, nanAuxiliaryInputs] = findNaNInputs(selectedData_cpx);

    % Reject finite negative values and Inf. Zero and NaN are retained so
    % equation-domain failures can be returned as non-stopping NaN results.
    validateInputValues(selectedData_cpx);

    row = calcTemp(selectedData_cpx, P_kbar);

    % Store the selected identifier once per pressure row.
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_cpx'}, 'Before', 1);

    % Store one result block in the fixed preallocated buffer.
    nResultBlocks = nResultBlocks + 1;
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature range.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_cpx)) ': ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_cpx)) ': ' ...
            num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    % Pressure warnings are common to all selected Cpx rows, so print them
    % only once during this function call.
    if ~pressureWarningsIssued
        if any(pressureOutsidePrincipal)
            fprintf(2, ...
                ['WARNING: Input pressure is outside the principal natural-like ' ...
                 'garnet-lherzolite range of Nimis and Taylor (2000): ' ...
                 '20-60 kbar. %d of %d pressure point(s) are outside; ' ...
                 'input range = %.6g-%.6g kbar.\n'], ...
                sum(pressureOutsidePrincipal), numel(P_kbar), ...
                min(P_kbar), max(P_kbar));
        end

        if any(pressureOutsideBroad)
            fprintf(2, ...
                ['WARNING: Input pressure is outside the broader experimental ' ...
                 'compilation of Nimis and Taylor (2000): 0-60 kbar. ' ...
                 '%d of %d pressure point(s) are outside; input range = ' ...
                 '%.6g-%.6g kbar. This is strong pressure extrapolation.\n'], ...
                sum(pressureOutsideBroad), numel(P_kbar), ...
                min(P_kbar), max(P_kbar));
        end
        pressureWarningsIssued = true;
    end

    % Temperature warnings for the principal natural-like range.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsidePrincipal = finiteTemperature & ...
        (row.T_deg < principalT_min_degC | ...
         row.T_deg > principalT_max_degC);

    if any(temperatureOutsidePrincipal)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the principal ' ...
             'natural-like garnet-lherzolite range of Nimis and Taylor ' ...
             '(2000): 900-1400 degreeC. %d of %d finite point(s) are ' ...
             'outside; calculated finite range = %.6g-%.6g degreeC for %s.\n'], ...
            sum(temperatureOutsidePrincipal), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_cpx)));
    end

    % Strong-extrapolation warning outside the broader compilation.
    temperatureOutsideBroad = finiteTemperature & ...
        (row.T_deg < broadT_min_degC | row.T_deg > broadT_max_degC);

    if any(temperatureOutsideBroad)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the broader ' ...
             'experimental compilation of Nimis and Taylor (2000): ' ...
             '850-1500 degreeC. %d of %d finite point(s) are outside; ' ...
             'calculated finite range = %.6g-%.6g degreeC for %s. ' ...
             'This is strong temperature extrapolation.\n'], ...
            sum(temperatureOutsideBroad), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_cpx)));
    end

    % Report every equation input containing NaN. These values remain NaN.
    if ~isempty(nanEquationInputs)
        fprintf(2, ...
            ['WARNING: NaN was found in Nimis and Taylor (2000) equation ' ...
             '(17) input(s) for %s: %s.\n' ...
             '         NaN values were retained and propagated; they were ' ...
             'not replaced by zero.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanEquationInputs, ', ')));
    end

    % Auxiliary values do not directly enter Eq. (17), but their NaN values
    % are still reported for data-quality and traceability purposes.
    if ~isempty(nanAuxiliaryInputs)
        fprintf(2, ...
            ['WARNING: NaN was found in auxiliary Cpx input(s) for %s: %s.\n' ...
             '         These values were retained. They do not directly enter ' ...
             'equation (17), but may affect quality diagnostics or metadata.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanAuxiliaryInputs, ', ')));
    end

    % Report non-finite temperature results without removing them.
    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), numel(row.T_deg), ...
            sum(isnan(row.T_deg)), sum(isinf(row.T_deg)));
    end

    % Report mathematical-domain failures explicitly.
    invalidEquationDomain = ~row.equation_domain_valid;
    if any(invalidEquationDomain)
        fprintf(2, ...
            ['WARNING: Nimis and Taylor (2000) equation (17) was outside ' ...
             'its mathematical domain for %s at %d of %d pressure point(s). ' ...
             'Both enstatite-activity factors, aEn_cpx, and the thermometer ' ...
             'denominator must be finite and > 0. Corresponding T values ' ...
             'remain NaN.\n'], ...
            char(string(selectedCode_cpx)), ...
            sum(invalidEquationDomain), numel(invalidEquationDomain));
    end

    % Report the natural-sample site-sum screening used in the paper.
    siteSumInvalid = ...
        (isfinite(row.T_site_sum_cpx) & row.T_site_sum_cpx < 1.990) | ...
        (isfinite(row.M1M2_site_sum_cpx) & row.M1M2_site_sum_cpx < 1.990);

    if any(siteSumInvalid)
        fprintf(2, ...
            ['WARNING: Cpx site-sum screening is outside the analytical ' ...
             'quality criterion used by Nimis and Taylor (2000), p. 548: ' ...
             'T-site and M1+M2-site sums should each be > 1.990 on a ' ...
             '6-oxygen basis. T-site sum = %.6g; M1+M2 sum = %.6g for %s.\n'], ...
            row.T_site_sum_cpx(1), row.M1M2_site_sum_cpx(1), ...
            char(string(selectedCode_cpx)));
    end

    % Fe3+ is stored but is not used separately in Eq. (17).
    if isfinite(row.Fe3_cpx(1)) && row.Fe3_cpx(1) > 0
        fprintf(2, ...
            ['WARNING: Fe3_cation_apfu is present and > 0 for %s, but ' ...
             'Nimis and Taylor (2000) equation (17) uses Fe_cation_apfu ' ...
             'as total Fe and applies no separate Fe3+ correction. The paper ' ...
             'notes that significant Fe3+ may cause temperature ' ...
             'underestimation (p. 549).\n'], ...
            char(string(selectedCode_cpx)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'NimisTaylor2000', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all result blocks once after all selections are complete.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(tbl, requiredVariables)
% validateRequiredVariables
% Confirm that every required Cpx input column is present.

missingVariables = requiredVariables(~ismember(requiredVariables, ...
    tbl.Properties.VariableNames));

if ~isempty(missingVariables)
    error('Cpx table is missing required variable(s): %s', ...
        char(strjoin(string(missingVariables), ', ')));
end

end

function [nanEquationInputs, nanAuxiliaryInputs] = findNaNInputs(data_cpx)
% findNaNInputs
% Return names of present Cpx inputs containing NaN. Equation (17) inputs
% and auxiliary inputs are listed separately. Missing optional columns are
% not listed because they follow the documented zero assumption.

% Equation (17) directly uses Al, Fe, Ca, Na, K, Ti, and Cr.
equationVariables = {'Al_cation_apfu', 'Fe_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu', 'K_cation_apfu', ...
    'Ti_cation_apfu', 'Cr_cation_apfu'};

% Si, Mg, Mn, and Fe3 are stored and/or used for quality diagnostics but do
% not directly enter equation (17).
auxiliaryVariables = {'Si_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Fe3_cation_apfu'};

nanEquationBuffer = strings(numel(equationVariables), 1);
nEquationNaN = 0;

for i = 1:numel(equationVariables)
    variableName = equationVariables{i};
    if ismember(variableName, data_cpx.Properties.VariableNames)
        value = data_cpx.(variableName);
        if any(isnan(value(:)))
            nEquationNaN = nEquationNaN + 1;
            nanEquationBuffer(nEquationNaN) = "Cpx." + string(variableName);
        end
    end
end

nanAuxiliaryBuffer = strings(numel(auxiliaryVariables), 1);
nAuxiliaryNaN = 0;

for i = 1:numel(auxiliaryVariables)
    variableName = auxiliaryVariables{i};
    if ismember(variableName, data_cpx.Properties.VariableNames)
        value = data_cpx.(variableName);
        if any(isnan(value(:)))
            nAuxiliaryNaN = nAuxiliaryNaN + 1;
            nanAuxiliaryBuffer(nAuxiliaryNaN) = ...
                "Cpx." + string(variableName);
        end
    end
end

nanEquationInputs = nanEquationBuffer(1:nEquationNaN);
nanAuxiliaryInputs = nanAuxiliaryBuffer(1:nAuxiliaryNaN);

end

function validateInputValues(data_cpx)
% validateInputValues
% Reject finite negative values and Inf in every present Cpx input used or
% stored by this function. Zero and NaN are allowed.

variablesToCheck = {'Si_cation_apfu', 'Al_cation_apfu', ...
    'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu', ...
    'Na_cation_apfu', 'K_cation_apfu', 'Mn_cation_apfu', ...
    'Ti_cation_apfu', 'Cr_cation_apfu', 'Fe3_cation_apfu'};

invalidBuffer = strings(numel(variablesToCheck), 1);
nInvalid = 0;

for i = 1:numel(variablesToCheck)
    variableName = variablesToCheck{i};

    if ~ismember(variableName, data_cpx.Properties.VariableNames)
        continue;
    end

    value = data_cpx.(variableName);
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value)
        error('Cpx variable %s must be a real numeric scalar.', variableName);
    end

    if isinf(value) || (isfinite(value) && value < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Cpx." + string(variableName);
    end
end

if nInvalid > 0
    invalidNames = invalidBuffer(1:nInvalid);
    error(['NimisTaylor2000: finite mineral inputs must be greater ' ...
           'than or equal to zero, and Inf is not permitted. Invalid ' ...
           'value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_cpx, P_kbar)
% calcTemp
% Calculate Nimis and Taylor (2000) equation (17) for one selected Cpx row
% over a scalar or vector of pressures. One output row is returned for each
% pressure value.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

% --- Prepare Cpx cation row ---
cpx = prepareCpxRow(data_cpx, 'Cpx');

% --- Calculate enstatite activity and quality diagnostics ---
activity = calcAenCpx(cpx);

% Thermometer denominator is composition dependent and common to all input
% pressure points.
if activity.domainValid
    denominator_scalar = 13.25 ...
        + 15.35 .* cpx.Ti ...
        + 4.50 .* cpx.Fe ...
        - 1.55 .* (cpx.Al + cpx.Cr - cpx.Na - cpx.K) ...
        + activity.lnAEn.^2;
else
    denominator_scalar = NaN;
end

baseDomainValid = activity.domainValid && ...
    isfinite(denominator_scalar) && denominator_scalar > 0;

T_K = NaN(nP, 1);
T_deg = NaN(nP, 1);
equationDomainValid = false(nP, 1);

if baseDomainValid
    numerator = 23166 + 39.28 .* P_kbar;
    validNumerator = isfinite(numerator) & numerator > 0;
    equationDomainValid(validNumerator) = true;
    T_K(validNumerator) = numerator(validNumerator) ./ denominator_scalar;
    T_deg(validNumerator) = T_K(validNumerator) - 273.15;
end

% --- Pack outputs ---
row = table();
row.P_kbar = P_kbar;

row.Si_cpx = repmat(cpx.Si, nP, 1);
row.Al_cpx = repmat(cpx.Al, nP, 1);
row.Fe_cpx = repmat(cpx.Fe, nP, 1);
row.Fe3_cpx = repmat(cpx.Fe3, nP, 1);
row.Mg_cpx = repmat(cpx.Mg, nP, 1);
row.Ca_cpx = repmat(cpx.Ca, nP, 1);
row.Na_cpx = repmat(cpx.Na, nP, 1);
row.K_cpx = repmat(cpx.K, nP, 1);
row.Mn_cpx = repmat(cpx.Mn, nP, 1);
row.Ti_cpx = repmat(cpx.Ti, nP, 1);
row.Cr_cpx = repmat(cpx.Cr, nP, 1);

row.total_cation_sum_cpx = repmat(cpx.totalCationSum, nP, 1);
row.AlIV_cpx = repmat(cpx.AlIV, nP, 1);
row.T_site_sum_cpx = repmat(cpx.TSiteSum, nP, 1);
row.M1M2_site_sum_cpx = repmat(cpx.M1M2SiteSum, nP, 1);

row.term_M2_En = repmat(activity.termM2, nP, 1);
row.term_M1_En = repmat(activity.termM1, nP, 1);
row.aEn_cpx = repmat(activity.aEn, nP, 1);
row.ln_aEn_cpx = repmat(activity.lnAEn, nP, 1);
row.denom_T = repmat(denominator_scalar, nP, 1);
row.equation_domain_valid = equationDomainValid;

row.T_K = T_K;
row.T_deg = T_deg;

end

function cpx = prepareCpxRow(data_px, mineralLabel)
% prepareCpxRow
% Extract one selected Cpx row. Missing optional columns follow documented
% zero assumptions. Present NaN values are retained and propagated.

if height(data_px) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

cpx = struct();
cpx.Si = getRequiredVar(data_px, 'Si_cation_apfu', mineralLabel);
cpx.Al = getRequiredVar(data_px, 'Al_cation_apfu', mineralLabel);
cpx.Fe = getRequiredVar(data_px, 'Fe_cation_apfu', mineralLabel);
cpx.Mg = getRequiredVar(data_px, 'Mg_cation_apfu', mineralLabel);
cpx.Ca = getRequiredVar(data_px, 'Ca_cation_apfu', mineralLabel);

cpx.Na = getOptionalVar(data_px, 'Na_cation_apfu', 0);
cpx.K = getOptionalVar(data_px, 'K_cation_apfu', 0);
cpx.Mn = getOptionalVar(data_px, 'Mn_cation_apfu', 0);
cpx.Ti = getOptionalVar(data_px, 'Ti_cation_apfu', 0);
cpx.Cr = getOptionalVar(data_px, 'Cr_cation_apfu', 0);
cpx.Fe3 = getOptionalVar(data_px, 'Fe3_cation_apfu', 0);

% Total Fe is counted once. Fe3 is stored separately but is not added to Fe.
cpx.totalCationSum = cpx.Si + cpx.Al + cpx.Fe + cpx.Mg + ...
    cpx.Ca + cpx.Na + cpx.K + cpx.Mn + cpx.Ti + cpx.Cr;

% Approximate site sums for the analytical-quality screening described on
% p. 548. Tetrahedral Al is limited by both the T-site vacancy and total Al.
if isnan(cpx.Si) || isnan(cpx.Al)
    cpx.AlIV = NaN;
else
    cpx.AlIV = min(max(2 - cpx.Si, 0), cpx.Al);
end

cpx.TSiteSum = cpx.Si + cpx.AlIV;
cpx.M1M2SiteSum = cpx.totalCationSum - cpx.TSiteSum;

end

function activity = calcAenCpx(cpx)
% calcAenCpx
% Calculate the enstatite activity approximation used in equation (17),
% preserving NaN and returning an explicit mathematical-domain flag.

activity = struct();
activity.termM2 = 1 - cpx.Ca - cpx.Na - cpx.K;
activity.termM1 = 1 - (cpx.Al + cpx.Cr + cpx.Na + cpx.K) ./ 2;
activity.aEn = activity.termM2 .* activity.termM1;

activity.domainValid = ...
    isfinite(activity.termM2) && activity.termM2 > 0 && ...
    isfinite(activity.termM1) && activity.termM1 > 0 && ...
    isfinite(activity.aEn) && activity.aEn > 0;

if activity.domainValid
    activity.lnAEn = log(activity.aEn);
else
    activity.lnAEn = NaN;
end

end

function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read one required real numeric scalar without modifying NaN.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value)
    error('%s variable %s must be a real numeric scalar.', ...
        mineralLabel, varName);
end

end

function value = getOptionalVar(tbl, varName, defaultValue)
% getOptionalVar
% Return defaultValue only when the optional column is absent. A present NaN
% is retained and is never replaced by zero.

if ismember(varName, tbl.Properties.VariableNames)
    value = tbl.(varName);
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value)
        error('Variable %s must be a real numeric scalar.', varName);
    end
else
    value = defaultValue;
end

end
