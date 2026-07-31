function results = Perchuk1985GrtOpx(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/Perchuk1985GrtOpx.m
% Tested structurally against MATLAB R2024b syntax
%
% Garnet-orthopyroxene Fe-Mg exchange thermometer
% Perchuk, L.L. et al. (1985)
% Journal of Metamorphic Geology, 3, 265-310
% DOI: https://doi.org/10.1111/j.1525-1314.1985.tb00321.x
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Orthopyroxene analysis (selected independently from tables) and calculates
% temperature using the Perchuk et al. (1985) Grt-Opx Fe2+-Mg thermometer.
%
% A scalar pressure from startThermoCalc_fixedP or a pressure vector from
% startThermoCalc_rangeP can be supplied. For each selected Grt-Opx pair,
% the output contains one row for every pressure value.
%
% The function supports repeated mineral-pair selections. Each result is
% buffered as a table block and the blocks are concatenated only once after
% the interactive selection loop has finished.
%
% -------------------------------------------------------------------------
% CALIBRATION BASIS, APPLICATION DOMAIN, AND APPLICATION NOTES
%
% Equation (A23) is not presented with one explicit rectangular experimental
% calibration range in temperature and pressure. It was derived by combining
% internally consistent thermodynamic parameters for reactions (A18) and
% (A22), the latter constrained experimentally, with garnet and Opx activity
% models (pp. 303-306). Consequently, the numerical P-T limits below are the
% natural Aldan-granulite application domain illustrated in this paper, not
% formal experimental calibration bounds for Eq. (A23):
%
%   Illustrated temperature domain : approximately 590-840 degreeC
%                                    (pp. 289, 294-295; Table 8)
%
%   Illustrated pressure domain    : approximately 3-7.5 kbar
%                                    (pp. 289, 294-295; Table 8)
%
% A specifically described Grt-Opx core-rim example gives approximately
% 750 degreeC/7.3 kbar for the cores and 590 degreeC/3 kbar for the rims
% (sample Sut-2/2, p. 289). This implementation uses the broader illustrated
% 590-840 degreeC and 3-7.5 kbar intervals for non-stopping warnings and
% labels them as natural-application domains rather than calibration ranges.
%
% Important compositional limits and cautions in the paper:
%
%   1) The garnet activity model was simplified for Ca-bearing metapelitic
%      garnets. Equations (A4)-(A9) are stated to be valid only for
%      0 < XCa_Grt < 0.30 and XFe_Grt >= 0.40 (p. 303). This function prints
%      non-stopping warnings outside those intervals.
%
%   2) Experimental data were unavailable for Ca-Mn garnets. The Ca-Mn
%      interaction was assumed equal to Ca-Fe, and Mg-Mn was assumed ideal
%      (p. 303). Mn-rich garnets therefore require particular caution.
%
%   3) The Opx activity model treats Opx as a regular mixture of En, Fs, and
%      a fictive ortho-Al2O3 component. It was fitted to MAS and FMAS
%      experimental data; pressure dependence of Opx activities was
%      neglected (p. 304). Opx Al is used explicitly and must not silently
%      default to zero when unmeasured.
%
%   4) For garnets with high Cr2O3, characteristic of mantle xenoliths, an
%      additional term must be added to Eq. (A23) (p. 306). The paper gives
%      no numerical Cr threshold, and that correction is not implemented
%      here. If a positive Cr_cation_apfu value is present, a caution is
%      printed so the user can assess whether the high-Cr correction is
%      required.
%
%   5) At low to intermediate pressure the pressure term may sometimes be
%      neglected, but the authors state that Eq. (A23) should in most cases
%      be solved simultaneously with a geobarometric equation (p. 306).
%      P in Eq. (A23) is in bar; P_kbar is converted internally to bar.
%
%   6) KD represents ferrous Fe2+-Mg exchange among almandine-pyrope and
%      ferrosilite-enstatite components (reaction A21, pp. 305-306). This
%      implementation uses Fe_cation_apfu as Fe2+ and does not add
%      Fe3_cation_apfu to KD or the mineral mole fractions.
%
%   7) Garnet and other Fe-Mg minerals in the Aldan granulites commonly show
%      zoning produced or modified during retrogression. Core and rim pairs
%      can record markedly different P-T conditions, and much of the
%      calculated record reflects retrograde Fe-Mg exchange rather than the
%      peak conditions (pp. 289, 294-295). Analyses should be paired by
%      textural generation and position; unrelated cores and rims or
%      secondary reaction products should not be combined.
%
%   8) The paper does not report a universal statistical precision for
%      Eq. (A23). Precision values quoted for other thermometers in the
%      Appendix must not be assigned to this Grt-Opx thermometer (pp. 305-306).
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside the illustrated 3-7.5 kbar domain,
%   2) a finite calculated temperature is outside the illustrated
%      590-840 degreeC domain,
%   3) XCa_Grt is outside 0-0.30 or XFe_Grt is below 0.40,
%   4) positive garnet Cr is reported and the unimplemented Cr correction
%      may need consideration,
%   5) a thermometer input contains NaN, or
%   6) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain either:
%   rawdata_struct.Garnet : table
% or
%   rawdata_struct.Grt    : table
%
% and either:
%   rawdata_struct.Orthopyroxene : table
% or
%   rawdata_struct.Opx           : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns contain
% normalised cation data.
%
% Required thermometer variables:
%   Garnet:
%     Fe_cation_apfu       % Fe2+
%     Mg_cation_apfu
%     Mn_cation_apfu
%     Ca_cation_apfu
%
%   Orthopyroxene:
%     Fe_cation_apfu       % Fe2+
%     Mg_cation_apfu
%     Al_cation_apfu       % used as Al/2 in XFe_Opx and XEn_Opx
%
% Optional variables retained in the output when present:
%     Fe3_cation_apfu, Si_cation_apfu, Cr_cation_apfu
%   Garnet Al and Opx Mn/Ca are also retained when present.
%
% An optional output-only column that is absent defaults to zero. An
% explicitly stored NaN is retained and is never converted to zero.
% Required equation inputs never default to zero when a column is absent.
%
% All finite values actually used by the calculation must be >= 0. Negative
% or infinite values stop the calculation with an error. NaN values remain
% NaN, propagate through the calculation, and trigger non-stopping fprintf
% diagnostics. Zero is permitted as an input value, but zero Fe or Mg, an
% invalid mole-fraction denominator, or an invalid logarithm produces NaN;
% the NaN result is retained and reported without stopping the calculation.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Perchuk et al. (1985), Appendix Eq. (A23), pp. 305-306:
%
%       4766 + 2533*(XFe_Opx - XEn_Opx) + 0.023*P_bar
% T(K)= ---------------------------------------------------------
%       R*ln(KD) + 2.65 - 5.214*XMg_Grt + 5.704*XFe_Grt
%
%   T(degreeC) = T(K) - 273.15
%
% where:
%
%   KD       = (Fe2+/Mg)_Grt / (Fe2+/Mg)_Opx
%   R        = 1.987 cal K^-1 mol^-1
%   P_bar    = 1000 * P_kbar
%
%   XFe_Opx  = Fe2_Opx/(Fe2_Opx + Mg_Opx + Al_Opx/2)
%   XEn_Opx  = Mg_Opx /(Fe2_Opx + Mg_Opx + Al_Opx/2)
%
%   XFe_Grt  = Fe2_Grt/(Fe2_Grt + Mg_Grt + Mn_Grt + Ca_Grt)
%   XMg_Grt  = Mg_Grt /(Fe2_Grt + Mg_Grt + Mn_Grt + Ca_Grt)
%   XMn_Grt  = Mn_Grt /(Fe2_Grt + Mg_Grt + Mn_Grt + Ca_Grt)
%   XCa_Grt  = Ca_Grt /(Fe2_Grt + Mg_Grt + Mn_Grt + Ca_Grt)
%
% Fe3+ is deliberately excluded from the exchange and mole-fraction terms.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Perchuk1985GrtOpx(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet/Grt and Orthopyroxene/Opx
%                    tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Opx pair

%% Input validation
if nargin < 2
    error('Perchuk1985GrtOpx requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || ~isreal(P_kbar) || isempty(P_kbar) || ...
        ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
disp('=== Step 1: Preparing cation dataset ===');

if isfield(rawdata_struct, 'Garnet') && istable(rawdata_struct.Garnet)
    dataset_gt = rawdata_struct.Garnet;
elseif isfield(rawdata_struct, 'Grt') && istable(rawdata_struct.Grt)
    dataset_gt = rawdata_struct.Grt;
else
    error(['rawdata_struct must contain garnet table as either ' ...
        'rawdata_struct.Garnet or rawdata_struct.Grt']);
end

if isfield(rawdata_struct, 'Orthopyroxene') && ...
        istable(rawdata_struct.Orthopyroxene)
    dataset_opx = rawdata_struct.Orthopyroxene;
elseif isfield(rawdata_struct, 'Opx') && istable(rawdata_struct.Opx)
    dataset_opx = rawdata_struct.Opx;
else
    error(['rawdata_struct must contain orthopyroxene table as either ' ...
        'rawdata_struct.Orthopyroxene or rawdata_struct.Opx']);
end

validateRequiredVariables(dataset_gt, dataset_opx);

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each calculation as one table block. The cell buffer is preallocated
% and doubled only when its capacity is exhausted; the full output table is
% concatenated once after the interactive loop.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

illustratedT_min_degC = 590;
illustratedT_max_degC = 840;
illustratedP_min_kbar = 3;
illustratedP_max_kbar = 7.5;
modelXCa_min = 0;
modelXCa_max = 0.30;
modelXFe_min = 0.40;

pressureOutsideIllustratedDomain = ...
    P_kbar < illustratedP_min_kbar | P_kbar > illustratedP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
    dataCodes_gt = dataset_gt{:, 1};

    [selectedIdx_gt, ok] = listdlg( ...
        'PromptString', 'Please select the Garnet data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_gt)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_gt)
        disp('Selection canceled');
        break;
    end

    selectedCode_gt = dataCodes_gt(selectedIdx_gt);
    disp(['Garnet selected: ' char(string(selectedCode_gt))]);

    % ----- Orthopyroxene selection -----
    disp('=== Step 4: Selecting a data code from the list (Orthopyroxene) ===');

    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', ...
        'Please select the Orthopyroxene data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Orthopyroxene selected: ' char(string(selectedCode_opx))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_gt = dataset_gt(selectedIdx_gt, :);
    selectedData_opx = dataset_opx(selectedIdx_opx, :);

    validateNonNegativeInputs(selectedData_gt, selectedData_opx);
    nanInputNames = findNaNInputs(selectedData_gt, selectedData_opx);

    row = calcTemp(selectedData_gt, selectedData_opx, P_kbar);

    row.dataCode_gt = repmat(string(selectedCode_gt), height(row), 1);
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row = movevars(row, {'dataCode_gt', 'dataCode_opx'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature values.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_gt)) ' & ' ...
            char(string(selectedCode_opx)) ': ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_gt)) ' & ' ...
            char(string(selectedCode_opx)) ': ' ...
            num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    % The paper gives no formal rectangular P-T calibration range for A23.
    % This warning therefore refers explicitly to the illustrated natural
    % Aldan-granulite application domain and is printed only once per call.
    if any(pressureOutsideIllustratedDomain) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the natural Aldan-' ...
             'granulite application domain illustrated by Perchuk et al. ' ...
             '(1985): approximately 3-7.5 kbar. %d of %d pressure ' ...
             'point(s) are outside the domain; input range = ' ...
             '%.4g-%.4g kbar (pp. 289, 294-295; Table 8).\n' ...
             '         This is not a formal experimental calibration ' ...
             'range for Eq. (A23). In most cases the authors recommend ' ...
             'solving the thermometer with a geobarometer (p. 306).\n'], ...
            sum(pressureOutsideIllustratedDomain), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideIllustratedDomain = finiteTemperature & ...
        (row.T_deg < illustratedT_min_degC | ...
         row.T_deg > illustratedT_max_degC);

    if any(temperatureOutsideIllustratedDomain)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the natural ' ...
             'Aldan-granulite application domain illustrated by ' ...
             'Perchuk et al. (1985): approximately 590-840 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the ' ...
             'domain; calculated finite range = %.4g-%.4g degreeC for ' ...
             '%s & %s (pp. 289, 294-295; Table 8).\n' ...
             '         This is not a formal experimental calibration ' ...
             'range for Eq. (A23).\n'], ...
            sum(temperatureOutsideIllustratedDomain), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_gt)), ...
            char(string(selectedCode_opx)));
    end

    % Garnet activity-model limits stated on p. 303.
    finiteXCa = isfinite(row.XCa_g);
    xCaOutsideModel = finiteXCa & ...
        (row.XCa_g <= modelXCa_min | row.XCa_g >= modelXCa_max);
    if any(xCaOutsideModel)
        fprintf(2, ...
            ['WARNING: Garnet XCa is outside the stated validity interval ' ...
             'of the activity model used by Perchuk et al. (1985): ' ...
             '0 < XCa_Grt < 0.30. Input XCa_Grt = %.8g for %s & %s ' ...
             '(Eqs. A4-A9, p. 303).\n'], ...
            row.XCa_g(1), ...
            char(string(selectedCode_gt)), ...
            char(string(selectedCode_opx)));
    end

    finiteXFe = isfinite(row.XFe_g);
    xFeOutsideModel = finiteXFe & row.XFe_g < modelXFe_min;
    if any(xFeOutsideModel)
        fprintf(2, ...
            ['WARNING: Garnet XFe is below the stated validity limit of ' ...
             'the activity model used by Perchuk et al. (1985): ' ...
             'XFe_Grt >= 0.40. Input XFe_Grt = %.8g for %s & %s ' ...
             '(Eqs. A4-A9, p. 303).\n'], ...
            row.XFe_g(1), ...
            char(string(selectedCode_gt)), ...
            char(string(selectedCode_opx)));
    end

    % No numerical high-Cr threshold is supplied in the paper. Report any
    % positive stored garnet Cr value as a prompt for manual assessment.
    finitePositiveCr = isfinite(row.Cr_g) & row.Cr_g > 0;
    if any(finitePositiveCr)
        fprintf(2, ...
            ['CAUTION: Garnet Cr_cation_apfu = %.8g for %s & %s. ' ...
             'Perchuk et al. (1985) state that high-Cr2O3 garnets, as in ' ...
             'mantle xenoliths, require an additional term in Eq. (A23); ' ...
             'the paper gives no numerical Cr threshold and the correction ' ...
             'is not implemented here (p. 306).\n'], ...
            row.Cr_g(1), ...
            char(string(selectedCode_gt)), ...
            char(string(selectedCode_opx)));
    end

    % Report explicitly stored NaN values in variables used by the equation.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for ' ...
             '%s & %s: %s.\n' ...
             '         The calculation was continued, and NaN was not ' ...
             'replaced by zero.\n'], ...
            char(string(selectedCode_gt)), ...
            char(string(selectedCode_opx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Preserve non-finite results and display all seven equation inputs.
    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        nonFiniteCauses = findNonFiniteCauses(row);

        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for ' ...
             '%s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_gt)), ...
            char(string(selectedCode_opx)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));

        fprintf(2, ...
            ['         Thermometer inputs used: ' ...
             'Garnet.Fe_cation_apfu=%s, Garnet.Mg_cation_apfu=%s, ' ...
             'Garnet.Mn_cation_apfu=%s, Garnet.Ca_cation_apfu=%s, ' ...
             'Opx.Fe_cation_apfu=%s, Opx.Mg_cation_apfu=%s, ' ...
             'Opx.Al_cation_apfu=%s.\n'], ...
            formatNumericValue(row.Fe2_g(1)), ...
            formatNumericValue(row.Mg_g(1)), ...
            formatNumericValue(row.Mn_g(1)), ...
            formatNumericValue(row.Ca_g(1)), ...
            formatNumericValue(row.Fe2_opx(1)), ...
            formatNumericValue(row.Mg_opx(1)), ...
            formatNumericValue(row.Al_opx(1)));

        if isempty(nonFiniteCauses)
            fprintf(2, ...
                ['         No explicit NaN, zero, Inf, or invalid ' ...
                 'intermediate value was identified; inspect the stored ' ...
                 'intermediate variables.\n']);
        else
            fprintf(2, '         Identified cause(s): %s.\n', ...
                char(strjoin(nonFiniteCauses, '; ')));
        end
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Perchuk1985GrtOpx', ...
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
function validateRequiredVariables(dataset_gt, dataset_opx)
% validateRequiredVariables
% Verify all table columns used by Eq. (A23).

requiredGt = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
requiredOpx = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Al_cation_apfu'};

missingNames = strings(numel(requiredGt) + numel(requiredOpx), 1);
nMissing = 0;

for i = 1:numel(requiredGt)
    if ~ismember(requiredGt{i}, dataset_gt.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Garnet." + string(requiredGt{i});
    end
end

for i = 1:numel(requiredOpx)
    if ~ismember(requiredOpx{i}, dataset_opx.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Opx." + string(requiredOpx{i});
    end
end

if nMissing > 0
    missingNames = missingNames(1:nMissing);
    error(['Perchuk1985GrtOpx: required variable(s) are missing: ' ...
        char(strjoin(missingNames, ', ')) '.']);
end

end


function nanInputNames = findNaNInputs(data_gt, data_opx)
% findNaNInputs
% Return names and values of Eq. (A23) inputs containing stored NaN values.

gtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
opxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Al_cation_apfu'};

nanInputNames = strings(numel(gtVariables) + numel(opxVariables), 1);
nNaN = 0;

for i = 1:numel(gtVariables)
    variableName = gtVariables{i};
    variableValue = data_gt.(variableName);
    if any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanInputNames(nNaN) = "Garnet." + string(variableName) + "=NaN";
    end
end

for i = 1:numel(opxVariables)
    variableName = opxVariables{i};
    variableValue = data_opx.(variableName);
    if any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanInputNames(nNaN) = "Opx." + string(variableName) + "=NaN";
    end
end

nanInputNames = nanInputNames(1:nNaN);

end


function validateNonNegativeInputs(data_gt, data_opx)
% validateNonNegativeInputs
% Reject negative or infinite Eq. (A23) inputs. Explicit NaN is retained.

gtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
opxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Al_cation_apfu'};

invalidInputNames = strings(numel(gtVariables) + numel(opxVariables), 1);
nInvalid = 0;

for i = 1:numel(gtVariables)
    variableName = gtVariables{i};
    variableValue = data_gt.(variableName);
    if ~isnumeric(variableValue) || ~isreal(variableValue) || ...
            ~isscalar(variableValue)
        error('Garnet.%s must be a numeric scalar in the selected row.', ...
            variableName);
    end
    if isinf(variableValue) || variableValue < 0
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(opxVariables)
    variableName = opxVariables{i};
    variableValue = data_opx.(variableName);
    if ~isnumeric(variableValue) || ~isreal(variableValue) || ...
            ~isscalar(variableValue)
        error('Opx.%s must be a numeric scalar in the selected row.', ...
            variableName);
    end
    if isinf(variableValue) || variableValue < 0
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Opx." + string(variableName);
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNames(1:nInvalid);
    error(['Perchuk1985GrtOpx: thermometer inputs must be finite or ' ...
        'NaN and must be >= 0. Negative or infinite value(s) were found ' ...
        'in: ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function row = calcTemp(data_gt, data_opx, P_kbar)
% calcTemp
% Compute Perchuk et al. (1985) temperatures for one Grt-Opx pair and a
% scalar or vector of pressures.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_bar = P_bar;

R_cal = 1.987;
row.R_cal = repmat(R_cal, nP, 1);

% --- Extract and expand garnet cations ---
Fe2_g = repmat(getRequiredValue(data_gt, ...
    'Fe_cation_apfu', 'Garnet'), nP, 1);
Fe3_g = repmat(getOptionalValue(data_gt, ...
    'Fe3_cation_apfu', 0, 'Garnet'), nP, 1);
Mg_g = repmat(getRequiredValue(data_gt, ...
    'Mg_cation_apfu', 'Garnet'), nP, 1);
Mn_g = repmat(getRequiredValue(data_gt, ...
    'Mn_cation_apfu', 'Garnet'), nP, 1);
Ca_g = repmat(getRequiredValue(data_gt, ...
    'Ca_cation_apfu', 'Garnet'), nP, 1);
Si_g = repmat(getOptionalValue(data_gt, ...
    'Si_cation_apfu', 0, 'Garnet'), nP, 1);
Al_g = repmat(getOptionalValue(data_gt, ...
    'Al_cation_apfu', 0, 'Garnet'), nP, 1);
Cr_g = repmat(getOptionalValue(data_gt, ...
    'Cr_cation_apfu', 0, 'Garnet'), nP, 1);

% --- Extract and expand orthopyroxene cations ---
Fe2_opx = repmat(getRequiredValue(data_opx, ...
    'Fe_cation_apfu', 'Opx'), nP, 1);
Fe3_opx = repmat(getOptionalValue(data_opx, ...
    'Fe3_cation_apfu', 0, 'Opx'), nP, 1);
Mg_opx = repmat(getRequiredValue(data_opx, ...
    'Mg_cation_apfu', 'Opx'), nP, 1);
Mn_opx = repmat(getOptionalValue(data_opx, ...
    'Mn_cation_apfu', 0, 'Opx'), nP, 1);
Si_opx = repmat(getOptionalValue(data_opx, ...
    'Si_cation_apfu', 0, 'Opx'), nP, 1);
Al_opx = repmat(getRequiredValue(data_opx, ...
    'Al_cation_apfu', 'Opx'), nP, 1);
Ca_opx = repmat(getOptionalValue(data_opx, ...
    'Ca_cation_apfu', 0, 'Opx'), nP, 1);
Cr_opx = repmat(getOptionalValue(data_opx, ...
    'Cr_cation_apfu', 0, 'Opx'), nP, 1);

% --- Perchuk et al. (1985), Eq. (A23) ---
% Fe3+ is deliberately excluded from every equation term.
sumDivalent_g = Fe2_g + Mg_g + Mn_g + Ca_g;
XCa_g = Ca_g ./ sumDivalent_g;
XMn_g = Mn_g ./ sumDivalent_g;
XFe_g = Fe2_g ./ sumDivalent_g;
XMg_g = Mg_g ./ sumDivalent_g;

FeMg_garnet = Fe2_g ./ Mg_g;
FeMg_opx = Fe2_opx ./ Mg_opx;
KD = FeMg_garnet ./ FeMg_opx;
lnKD = log(KD);

sumFeMgAlHalf_opx = Fe2_opx + Mg_opx + Al_opx ./ 2;
XFe_opx = Fe2_opx ./ sumFeMgAlHalf_opx;
XEn_opx = Mg_opx ./ sumFeMgAlHalf_opx;
XAl_opx = (Al_opx ./ 2) ./ sumFeMgAlHalf_opx;

temperatureNumerator = ...
    4766 + 2533 .* (XFe_opx - XEn_opx) + 0.023 .* P_bar;
temperatureDenominator = ...
    R_cal .* lnKD + 2.65 - 5.214 .* XMg_g + 5.704 .* XFe_g;

calculationDomainValid = ...
    isfinite(Fe2_g) & Fe2_g > 0 ...
    & isfinite(Mg_g) & Mg_g > 0 ...
    & isfinite(Mn_g) & Mn_g >= 0 ...
    & isfinite(Ca_g) & Ca_g >= 0 ...
    & isfinite(Fe2_opx) & Fe2_opx > 0 ...
    & isfinite(Mg_opx) & Mg_opx > 0 ...
    & isfinite(Al_opx) & Al_opx >= 0 ...
    & isfinite(sumDivalent_g) & sumDivalent_g > 0 ...
    & isfinite(sumFeMgAlHalf_opx) & sumFeMgAlHalf_opx > 0 ...
    & isfinite(FeMg_garnet) & FeMg_garnet > 0 ...
    & isfinite(FeMg_opx) & FeMg_opx > 0 ...
    & isfinite(KD) & KD > 0 ...
    & isfinite(lnKD) ...
    & isfinite(XCa_g) & XCa_g >= 0 & XCa_g <= 1 ...
    & isfinite(XMn_g) & XMn_g >= 0 & XMn_g <= 1 ...
    & isfinite(XFe_g) & XFe_g >= 0 & XFe_g <= 1 ...
    & isfinite(XMg_g) & XMg_g >= 0 & XMg_g <= 1 ...
    & isfinite(XFe_opx) & XFe_opx >= 0 & XFe_opx <= 1 ...
    & isfinite(XEn_opx) & XEn_opx >= 0 & XEn_opx <= 1 ...
    & isfinite(XAl_opx) & XAl_opx >= 0 & XAl_opx <= 1 ...
    & isfinite(temperatureNumerator) ...
    & isfinite(temperatureDenominator) & temperatureDenominator > 0;

T_K = NaN(nP, 1);
T_K(calculationDomainValid) = ...
    temperatureNumerator(calculationDomainValid) ./ ...
    temperatureDenominator(calculationDomainValid);
T_deg = T_K - 273.15;

Mg_number_gt = Mg_g ./ (Mg_g + Fe2_g);
Mg_number_opx = Mg_opx ./ (Mg_opx + Fe2_opx);

% --- Pack outputs ---
row.Fe2_g = Fe2_g;
row.Fe_raw_g = Fe2_g;
row.Fe3_g = Fe3_g;
row.Mg_g = Mg_g;
row.Mn_g = Mn_g;
row.Ca_g = Ca_g;
row.Si_g = Si_g;
row.Al_g = Al_g;
row.Cr_g = Cr_g;

row.Fe2_opx = Fe2_opx;
row.Fe_raw_opx = Fe2_opx;
row.Fe3_opx = Fe3_opx;
row.Mg_opx = Mg_opx;
row.Mn_opx = Mn_opx;
row.Si_opx = Si_opx;
row.Al_opx = Al_opx;
row.Ca_opx = Ca_opx;
row.Cr_opx = Cr_opx;

row.FeMg_garnet = FeMg_garnet;
row.FeMg_opx = FeMg_opx;
row.KD = KD;
row.lnKD = lnKD;

row.sumDivalent_g = sumDivalent_g;
row.XCa_g = XCa_g;
row.XMn_g = XMn_g;
row.XFe_g = XFe_g;
row.XMg_g = XMg_g;

row.XFe_opx = XFe_opx;
row.XEn_opx = XEn_opx;
row.XAl_opx = XAl_opx;
row.opx_denominator_FeMgAlHalf = sumFeMgAlHalf_opx;

row.Mg_number_garnet = Mg_number_gt;
row.Mg_number_opx = Mg_number_opx;
row.temperatureNumerator = temperatureNumerator;
row.denominator = temperatureDenominator;
row.calculationDomainValid = calculationDomainValid;
row.T_K = T_K;
row.T_deg = T_deg;

row.is_positive_FeMg_garnet = ...
    isfinite(Fe2_g) & Fe2_g > 0 & isfinite(Mg_g) & Mg_g > 0;
row.is_positive_FeMg_opx = ...
    isfinite(Fe2_opx) & Fe2_opx > 0 & isfinite(Mg_opx) & Mg_opx > 0;
row.is_positive_KD = isfinite(KD) & KD > 0;
row.is_positive_denominator = ...
    isfinite(temperatureDenominator) & temperatureDenominator > 0;
row.is_XCa_g_reasonable = isfinite(XCa_g) & XCa_g >= 0 & XCa_g <= 1;
row.is_XMn_g_reasonable = isfinite(XMn_g) & XMn_g >= 0 & XMn_g <= 1;
row.is_XEn_opx_reasonable = ...
    isfinite(XEn_opx) & XEn_opx >= 0 & XEn_opx <= 1;
row.is_XFe_opx_reasonable = ...
    isfinite(XFe_opx) & XFe_opx >= 0 & XFe_opx <= 1;

end


function nonFiniteCauses = findNonFiniteCauses(row)
% findNonFiniteCauses
% Identify input and derived-variable conditions that can produce NaN or Inf.

maximumCauses = 28;
nonFiniteCauses = strings(maximumCauses, 1);
nCauses = 0;

inputValues = {row.Fe2_g, row.Mg_g, row.Mn_g, row.Ca_g, ...
    row.Fe2_opx, row.Mg_opx, row.Al_opx};
inputNames = {"Garnet.Fe_cation_apfu", "Garnet.Mg_cation_apfu", ...
    "Garnet.Mn_cation_apfu", "Garnet.Ca_cation_apfu", ...
    "Opx.Fe_cation_apfu", "Opx.Mg_cation_apfu", ...
    "Opx.Al_cation_apfu"};
zeroIsInvalid = [true, true, false, false, true, true, false];

for i = 1:numel(inputValues)
    value = inputValues{i};
    if any(isnan(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is NaN";
    elseif any(isinf(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is Inf";
    elseif zeroIsInvalid(i) && any(value(:) == 0)
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is zero";
    end
end

derivedValues = {row.sumDivalent_g, ...
    row.opx_denominator_FeMgAlHalf, row.FeMg_garnet, row.FeMg_opx, ...
    row.KD, row.lnKD, row.XCa_g, row.XMn_g, row.XFe_g, row.XMg_g, ...
    row.XFe_opx, row.XEn_opx, row.XAl_opx, ...
    row.temperatureNumerator, row.denominator};
derivedNames = {"garnet divalent-cation sum", ...
    "Opx Fe-Mg-Al/2 sum", "garnet Fe/Mg", "Opx Fe/Mg", ...
    "KD", "lnKD", "XCa_Grt", "XMn_Grt", "XFe_Grt", "XMg_Grt", ...
    "XFe_Opx", "XEn_Opx", "XAl_Opx", ...
    "temperature numerator", "temperature denominator"};

for i = 1:numel(derivedValues)
    value = derivedValues{i};
    if any(~isfinite(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is non-finite";
    elseif any(value(:) <= 0) && ...
            (derivedNames{i} == "garnet divalent-cation sum" || ...
             derivedNames{i} == "Opx Fe-Mg-Al/2 sum" || ...
             derivedNames{i} == "garnet Fe/Mg" || ...
             derivedNames{i} == "Opx Fe/Mg" || ...
             derivedNames{i} == "KD" || ...
             derivedNames{i} == "temperature denominator")
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is zero or negative";
    end
end

if any(~row.calculationDomainValid)
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = ...
        "Perchuk et al. (1985) calculation domain is invalid";
end

nonFiniteCauses = nonFiniteCauses(1:nCauses);

end


function textValue = formatNumericValue(value)
% formatNumericValue
% Format one scalar value for compact fprintf diagnostics.

if isnan(value)
    textValue = 'NaN';
elseif isinf(value) && value > 0
    textValue = 'Inf';
elseif isinf(value) && value < 0
    textValue = '-Inf';
else
    textValue = sprintf('%.8g', value);
end

end


function value = getRequiredValue(data_tbl, variableName, mineralLabel)
% getRequiredValue
% Extract a required scalar numeric value; an explicit NaN is unchanged.

if ~ismember(variableName, data_tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = data_tbl.(variableName);
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end


function value = getOptionalValue(data_tbl, variableName, defaultValue, mineralLabel)
% getOptionalValue
% Use defaultValue only when a column is absent. Explicit NaN is unchanged.

if ismember(variableName, data_tbl.Properties.VariableNames)
    value = data_tbl.(variableName);
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
        error('%s.%s must be a numeric scalar in the selected row.', ...
            mineralLabel, variableName);
    end
else
    value = defaultValue;
end

end
