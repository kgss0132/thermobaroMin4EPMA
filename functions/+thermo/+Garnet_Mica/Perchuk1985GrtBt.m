function results = Perchuk1985GrtBt(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/Perchuk1985GrtBt.m
% Tested with MATLAB R2024b
%
% Garnet-biotite Fe-Mg exchange geothermometer
% Perchuk, L.L. et al. (1985)
% Journal of Metamorphic Geology, 3, 265-310
% DOI: https://doi.org/10.1111/j.1525-1314.1985.tb00321.x
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Biotite
% analysis and calculates the standard and F-corrected temperatures of
% Perchuk et al. (1985), equations (A19) and (A20), respectively.
%
% The function accepts pressure as either a scalar or a vector. Therefore,
% it can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for every pressure
% value for each user-selected Garnet-Biotite pair.
%
% -------------------------------------------------------------------------
% CALIBRATION BASIS, APPLICATION RANGE, AND IMPORTANT CAUTIONS
%
% Perchuk et al. (1985) is primarily a petrological and thermodynamic study
% of Aldan-shield granulites, not a stand-alone experimental calibration
% paper. The experimental calibration is referred to Perchuk & Lavrent'eva
% (1983), and the 1985 paper does not state one formal numerical temperature
% or pressure calibration range (Appendix, p. 305).
%
% The following interval describes the main core-rim P-T applications in
% the Aldan granulites; it is NOT a formal calibration interval:
%
%   Direct application T : approximately 590-840 degreeC
%   Direct application P : approximately 3-7.5 kbar
%   Main application     : high-grade metapelites and granulites
%
% These values are illustrated by the core-rim applications and Table 8
% discussion on p. 289. An explicit F-corrected garnet-biotite result of
% 778 degreeC is reported on p. 295. Values outside the interval are treated
% here as extrapolation beyond the paper's directly demonstrated Aldan
% application, not as automatically invalid results.
%
% Important application cautions stated or implied by the original paper:
% - The garnet activity model used by the thermometer is stated to be valid
%   only for 0 < XCa_Gr < 0.3 and XFe_Gr >= 0.4 (Appendix, p. 303).
% - Equation (A19) has a reported statistical precision of approximately
%   +/-15 degreeC. This is a regression precision, not a guarantee of
%   absolute accuracy for natural samples (Appendix, p. 305).
% - Fluorine increases Mg in biotite at constant temperature. Ignoring F in
%   F-rich biotite produces an artificial temperature underestimate in the
%   standard garnet-biotite thermometer (pp. 279 and 305).
% - Equation (A20) defines XF as F/(F+OH), and n as Mg+Fe in the biotite
%   formula calculated per 22 negative charges (p. 305). In this function,
%   the input variable F_anion must therefore contain XF, not an unconverted
%   raw F apfu value, unless the user's normalization makes them identical.
% - The garnet Ca-Mn interaction was assumed equal to the Ca-Fe interaction
%   because experimental Ca-Mn garnet data were unavailable. Mn-rich
%   garnets require additional caution (p. 303).
% - Biotite was treated as an ideal multi-site solid solution. Strongly
%   non-ideal or unusual Ti-, Al-, or Fe3+-rich compositions may not be
%   represented adequately (p. 304).
% - Garnet and biotite commonly preserve diffusion zoning and retrograde
%   Fe-Mg exchange. Pair compositions from the same equilibrium stage, such
%   as contacting rim-rim or corresponding core-core compositions
%   (pp. 278-279 and 289-295).
% - Biotite inclusions in garnet may have re-equilibrated their Fe/Mg ratio
%   with the host even when Ti retains evidence of an earlier high-
%   temperature stage (p. 279).
%
% This implementation therefore issues non-stopping fprintf messages when:
%   1) pressure is outside the 3-7.5 kbar direct-application interval,
%   2) a finite standard or F-corrected temperature is outside the
%      590-840 degreeC direct-application interval,
%   3) finite garnet composition is outside the stated XCa-XFe model range,
%   4) supplied XF is outside 0-1,
%   5) a calculation input is NaN, or
%   6) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain one garnet table named:
%   rawdata_struct.Garnet or rawdata_struct.Grt
%
% and one biotite table named:
%   rawdata_struct.Mica, rawdata_struct.Bt, or rawdata_struct.Bio
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% Required Garnet variables:
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%
% Required Biotite variables:
%   Fe_cation_apfu
%   Mg_cation_apfu
%
% Optional calculation variables:
%   Garnet : Fe3_cation_apfu, Mn_cation_apfu
%   Biotite: Fe3_cation_apfu, F_anion
%
% Optional output-only variables retained when present:
%   Garnet : Si_cation_apfu, Al_cation_apfu
%   Biotite: Mn_cation_apfu, Si_cation_apfu, Al_cation_apfu,
%            Ti_cation_apfu, Na_cation_apfu, K_cation_apfu
%
% This implementation preserves the original project convention:
%   Fe_used = Fe_cation_apfu + Fe3_cation_apfu
%
% when a separate Fe3_cation_apfu column is present. This is appropriate
% only when Fe_cation_apfu and Fe3_cation_apfu are separate, non-overlapping
% components and the intended thermometer convention treats their sum as
% Fe. If Fe_cation_apfu already stores total Fe, do not add Fe3 again. The
% paper distinguishes Fe2+ and Fe3+ in its biotite formula (p. 304), so the
% user's Fe normalization and exchange convention must be checked.
%
% All finite mineral-composition values used in the calculation must be
% greater than or equal to zero. Negative finite values and Inf stop the
% calculation with an error. Zero is allowed, although it may generate a
% non-finite mathematical result that is retained and reported. A present
% NaN is never replaced by zero; it propagates through the calculation.
%
% Missing optional calculation columns retain the original convention:
%   missing Fe3_cation_apfu -> 0
%   missing Mn_cation_apfu in Garnet -> 0
%   missing F_anion in Biotite -> 0
%
% This missing-column behavior is distinct from a present NaN, which is
% retained. With XF = 0, equation (A20) reduces exactly to equation (A19).
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
%   KD = (Fe/Mg)_Gr / (Fe/Mg)_Bt
%
%   XCa_Gr = Ca_Gr / (Fe_Gr + Mg_Gr + Mn_Gr + Ca_Gr)
%   XFe_Gr = Fe_Gr / (Fe_Gr + Mg_Gr + Mn_Gr + Ca_Gr)
%
%   n = Mg_Bt + Fe_Bt
%   XF_over_n = XF_Bt / n
%   XF_Bt = F/(F+OH), supplied in F_anion
%
% Standard equation (A19):
%   T(K) = (3720 + 2871*XCa_Gr + 0.038*P_bar) ...
%          / (lnKD + 2.868 + 0.625*XCa_Gr)
%
% F-corrected equation (A20), using the user-confirmed coefficients:
%   T(K) = (3720 + 2871*XCa_Gr ...
%          + (0.038 - 0.02*XF_over_n)*P_bar ...
%          + 2469*XF_over_n) ...
%          / (lnKD + 2.868 + 0.625*XCa_Gr ...
%             - 3.27*XF_over_n)
%
% Pressure is supplied in kbar and converted internally:
%   P_bar = 1000 * P_kbar
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Perchuk1985GrtBt(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Biotite tables
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Biotite pair
%

%% Input validation
if nargin < 2
    error('Perchuk1985GrtBt requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
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

if isfield(rawdata_struct, 'Mica') && istable(rawdata_struct.Mica)
    dataset_bt = rawdata_struct.Mica;
elseif isfield(rawdata_struct, 'Bt') && istable(rawdata_struct.Bt)
    dataset_bt = rawdata_struct.Bt;
elseif isfield(rawdata_struct, 'Bio') && istable(rawdata_struct.Bio)
    dataset_bt = rawdata_struct.Bio;
else
    error(['rawdata_struct must contain biotite table as either ' ...
           'rawdata_struct.Mica, rawdata_struct.Bt, or rawdata_struct.Bio']);
end

requiredGarnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu'};
requiredBiotiteVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

validateRequiredVariables(dataset_gt, requiredGarnetVariables, 'Garnet');
validateRequiredVariables(dataset_bt, requiredBiotiteVariables, 'Biotite');
validateOptionalCalculationVariables(dataset_gt, ...
    {'Fe3_cation_apfu', 'Mn_cation_apfu'}, 'Garnet');
validateOptionalCalculationVariables(dataset_bt, ...
    {'Fe3_cation_apfu', 'F_anion'}, 'Biotite');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected pair as one table block and concatenate only once
% after the loop. This avoids reallocating the complete results table after
% every calculation.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

applicationT_min_degC = 590;
applicationT_max_degC = 840;
applicationP_min_kbar = 3;
applicationP_max_kbar = 7.5;

pressureOutsideApplication = ...
    P_kbar < applicationP_min_kbar | P_kbar > applicationP_max_kbar;
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

    % ----- Biotite selection -----
    disp('=== Step 4: Selecting a data code from the list (Biotite) ===');
    dataCodes_bt = dataset_bt{:, 1};

    [selectedIdx_bt, ok] = listdlg( ...
        'PromptString', 'Please select the Biotite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_bt)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_bt)
        disp('Selection canceled');
        break;
    end

    selectedCode_bt = dataCodes_bt(selectedIdx_bt);
    disp(['Biotite selected: ' char(string(selectedCode_bt))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_gt = dataset_gt(selectedIdx_gt, :);
    selectedData_bt = dataset_bt(selectedIdx_bt, :);

    nanInputNames = findNaNInputs(selectedData_gt, selectedData_bt);
    validateNonNegativeInputs(selectedData_gt, selectedData_bt);

    row = calcTemp(selectedData_gt, selectedData_bt, P_kbar);

    row.dataCode_gt = repmat(string(selectedCode_gt), height(row), 1);
    row.dataCode_bt = repmat(string(selectedCode_bt), height(row), 1);
    row = movevars(row, {'dataCode_gt', 'dataCode_bt'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo both temperature versions.
    disp('--------------------------------------------------');
    disp('=== Standard temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_gt)) ' & ' char(string(selectedCode_bt)) ...
            ': ' num2str(row.T_deg_standard) ' degreeC']);
    else
        disp([char(string(selectedCode_gt)) ' & ' char(string(selectedCode_bt)) ...
            ': ' num2str(row.T_deg_standard(1)) ' to ' ...
            num2str(row.T_deg_standard(end)) ' degreeC']);
    end

    disp('=== F-corrected temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_gt)) ' & ' char(string(selectedCode_bt)) ...
            ': ' num2str(row.T_deg_Fcorrected) ' degreeC']);
    else
        disp([char(string(selectedCode_gt)) ' & ' char(string(selectedCode_bt)) ...
            ': ' num2str(row.T_deg_Fcorrected(1)) ' to ' ...
            num2str(row.T_deg_Fcorrected(end)) ' degreeC']);
    end

    % Pressure warning: direct application interval, not formal calibration.
    if any(pressureOutsideApplication) && ~pressureWarningIssued
        fprintf(2, ...
            ['CAUTION: Input pressure is outside the directly demonstrated ' ...
             'Aldan application interval of Perchuk et al. (1985): ' ...
             'approximately 3-7.5 kbar (Table 8 discussion, p. 289). ' ...
             'The paper does not state a formal numerical calibration range. ' ...
             '%d of %d pressure point(s) are outside; input range = ' ...
             '%.4g-%.4g kbar. Calculation was continued.\n'], ...
            sum(pressureOutsideApplication), numel(P_kbar), ...
            min(P_kbar), max(P_kbar));
        pressureWarningIssued = true;
    end

    % Temperature warnings for the standard and F-corrected versions.
    printTemperatureRangeWarning(row.T_deg_standard, ...
        'standard equation (A19)', applicationT_min_degC, ...
        applicationT_max_degC, selectedCode_gt, selectedCode_bt);
    printTemperatureRangeWarning(row.T_deg_Fcorrected, ...
        'F-corrected equation (A20)', applicationT_min_degC, ...
        applicationT_max_degC, selectedCode_gt, selectedCode_bt);

    % Garnet composition-model range from Appendix p. 303.
    finiteComposition = isfinite(row.XCa_g) & isfinite(row.XFe_g);
    compositionOutsideModel = finiteComposition & ...
        ~(row.XCa_g > 0 & row.XCa_g < 0.3 & row.XFe_g >= 0.4);
    if any(compositionOutsideModel)
        fprintf(2, ...
            ['CAUTION: Garnet composition is outside the stated range of the ' ...
             'activity model used by Perchuk et al. (1985): ' ...
             '0 < XCa_Gr < 0.3 and XFe_Gr >= 0.4 (Appendix, p. 303). ' ...
             'XCa_Gr = %.4g and XFe_Gr = %.4g for %s & %s. ' ...
             'Calculation was continued.\n'], ...
            row.XCa_g(1), row.XFe_g(1), ...
            char(string(selectedCode_gt)), char(string(selectedCode_bt)));
    end

    % XF is a site fraction and should lie between zero and one.
    finiteXF = isfinite(row.XF_bt);
    invalidXFRange = finiteXF & (row.XF_bt < 0 | row.XF_bt > 1);
    if any(invalidXFRange)
        fprintf(2, ...
            ['CAUTION: F_anion is used as XF = F/(F+OH) in equation (A20) ' ...
             'and should lie within 0-1 (p. 305). Supplied XF = %.4g for ' ...
             '%s & %s. Confirm the F normalization. Calculation was continued.\n'], ...
            row.XF_bt(1), char(string(selectedCode_gt)), ...
            char(string(selectedCode_bt)));
    end

    % NaN inputs remain NaN and are reported without stopping.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; it was not replaced by zero.\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_bt)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidStandard = ~isfinite(row.T_deg_standard);
    invalidFcorrected = ~isfinite(row.T_deg_Fcorrected);
    if any(invalidStandard) || any(invalidFcorrected)
        fprintf(2, ...
            ['WARNING: Non-finite temperatures were calculated for %s & %s. ' ...
             'Standard: %d of %d (NaN: %d, Inf: %d); ' ...
             'F-corrected: %d of %d (NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_gt)), char(string(selectedCode_bt)), ...
            sum(invalidStandard), numel(row.T_deg_standard), ...
            sum(isnan(row.T_deg_standard)), sum(isinf(row.T_deg_standard)), ...
            sum(invalidFcorrected), numel(row.T_deg_Fcorrected), ...
            sum(isnan(row.T_deg_Fcorrected)), ...
            sum(isinf(row.T_deg_Fcorrected)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Perchuk1985GrtBt', ...
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
function validateRequiredVariables(dataset, requiredVariables, mineralLabel)
% validateRequiredVariables
% Confirm that required calculation columns exist and contain numeric data.

missingVariables = requiredVariables( ...
    ~ismember(requiredVariables, dataset.Properties.VariableNames));
if ~isempty(missingVariables)
    error('%s table is missing required variable(s): %s.', ...
        mineralLabel, char(strjoin(string(missingVariables), ', ')));
end

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    if ~isnumeric(dataset.(variableName))
        error('%s.%s must be numeric.', mineralLabel, variableName);
    end
end

end

function validateOptionalCalculationVariables(dataset, variableNames, mineralLabel)
% validateOptionalCalculationVariables
% Confirm that optional calculation columns are numeric when present.

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    if ismember(variableName, dataset.Properties.VariableNames) && ...
            ~isnumeric(dataset.(variableName))
        error('%s.%s must be numeric when present.', mineralLabel, variableName);
    end
end

end

function nanInputNames = findNaNInputs(data_garnet, data_biotite)
% findNaNInputs
% Return calculation input names containing NaN. The buffer is preallocated.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Fe3_cation_apfu', 'Mn_cation_apfu'};
biotiteVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu', 'F_anion'};

nameBuffer = strings(numel(garnetVariables) + numel(biotiteVariables), 1);
nNames = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if any(isnan(variableValue(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(biotiteVariables)
    variableName = biotiteVariables{i};
    if ismember(variableName, data_biotite.Properties.VariableNames)
        variableValue = data_biotite.(variableName);
        if any(isnan(variableValue(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Biotite." + string(variableName);
        end
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_garnet, data_biotite)
% validateNonNegativeInputs
% Reject negative finite values and Inf in variables used by either equation.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Fe3_cation_apfu', 'Mn_cation_apfu'};
biotiteVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu', 'F_anion'};

maxNames = numel(garnetVariables) + numel(biotiteVariables);
negativeNameBuffer = strings(maxNames, 1);
infiniteNameBuffer = strings(maxNames, 1);
nNegative = 0;
nInfinite = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nNegative = nNegative + 1;
            negativeNameBuffer(nNegative) = "Garnet." + string(variableName);
        end
        if any(isinf(variableValue(:)))
            nInfinite = nInfinite + 1;
            infiniteNameBuffer(nInfinite) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(biotiteVariables)
    variableName = biotiteVariables{i};
    if ismember(variableName, data_biotite.Properties.VariableNames)
        variableValue = data_biotite.(variableName);
        if any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nNegative = nNegative + 1;
            negativeNameBuffer(nNegative) = "Biotite." + string(variableName);
        end
        if any(isinf(variableValue(:)))
            nInfinite = nInfinite + 1;
            infiniteNameBuffer(nInfinite) = "Biotite." + string(variableName);
        end
    end
end

if nNegative > 0
    invalidNames = negativeNameBuffer(1:nNegative);
    error(['Perchuk1985GrtBt: calculation inputs must be >= 0. ' ...
           'Negative finite value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

if nInfinite > 0
    invalidNames = infiniteNameBuffer(1:nInfinite);
    error(['Perchuk1985GrtBt: Inf is not permitted in calculation ' ...
           'inputs. Inf was found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function printTemperatureRangeWarning(temperatureValues, equationLabel, ...
        minimumTemperature, maximumTemperature, garnetCode, biotiteCode)
% printTemperatureRangeWarning
% Report extrapolation beyond the paper's directly demonstrated interval.

finiteTemperature = isfinite(temperatureValues);
outsideApplication = finiteTemperature & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideApplication)
    finiteValues = temperatureValues(finiteTemperature);
    fprintf(2, ...
        ['CAUTION: Temperature from the %s is outside the directly ' ...
         'demonstrated Aldan application interval of Perchuk et al. (1985): ' ...
         'approximately 590-840 degreeC (Table 8 discussion, p. 289; ' ...
         'see also p. 295). The paper does not state a formal numerical ' ...
         'calibration range. %d of %d finite point(s) are outside; ' ...
         'calculated finite range = %.4g-%.4g degreeC for %s & %s. ' ...
         'Calculation was continued.\n'], ...
        equationLabel, sum(outsideApplication), sum(finiteTemperature), ...
        min(finiteValues), max(finiteValues), ...
        char(string(garnetCode)), char(string(biotiteCode)));
end

end

function row = calcTemp(data_garnet, data_biotite, P_kbar)
% calcTemp
% Calculate equations (A19) and (A20) for one Garnet-Biotite pair and one
% or more pressure values.

P_kbar = P_kbar(:);
P_bar = 1000 .* P_kbar;
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_bar = P_bar;

gt = prepareGarnetRow(data_garnet);
bt = prepareBiotiteRow(data_biotite);

sumDivalent_g_scalar = gt.Fe_used + gt.Mg + gt.Mn + gt.Ca;
XCa_g_scalar = gt.Ca ./ sumDivalent_g_scalar;
XFe_g_scalar = gt.Fe_used ./ sumDivalent_g_scalar;
XMn_g_scalar = gt.Mn ./ sumDivalent_g_scalar;

FeMg_garnet_scalar = gt.Fe_used ./ gt.Mg;
FeMg_biotite_scalar = bt.Fe_used ./ bt.Mg;
KD_scalar = FeMg_garnet_scalar ./ FeMg_biotite_scalar;
lnKD_scalar = log(KD_scalar);

n_bt_scalar = bt.Mg + bt.Fe_used;
XF_bt_scalar = bt.F;
XF_over_n_scalar = XF_bt_scalar ./ n_bt_scalar;

den_standard_scalar = lnKD_scalar + 2.868 + 0.625 .* XCa_g_scalar;
den_Fcorrected_scalar = den_standard_scalar - 3.27 .* XF_over_n_scalar;

num_standard = 3720 + 2871 .* XCa_g_scalar + 0.038 .* P_bar;
num_Fcorrected = 3720 + 2871 .* XCa_g_scalar ...
    + (0.038 - 0.02 .* XF_over_n_scalar) .* P_bar ...
    + 2469 .* XF_over_n_scalar;

if isfinite(den_standard_scalar) && den_standard_scalar > 0
    T_K_standard = num_standard ./ den_standard_scalar;
    T_deg_standard = T_K_standard - 273.15;
else
    T_K_standard = NaN(nP, 1);
    T_deg_standard = NaN(nP, 1);
end

if isfinite(den_Fcorrected_scalar) && den_Fcorrected_scalar > 0
    T_K_Fcorrected = num_Fcorrected ./ den_Fcorrected_scalar;
    T_deg_Fcorrected = T_K_Fcorrected - 273.15;
else
    T_K_Fcorrected = NaN(nP, 1);
    T_deg_Fcorrected = NaN(nP, 1);
end

Mg_number_gt_scalar = gt.Mg ./ (gt.Mg + gt.Fe_used);
Mg_number_bt_scalar = bt.Mg ./ (bt.Mg + bt.Fe_used);

% Preserve the original output names while repeating scalar composition
% values to match the pressure-vector height.
row.Fe2_g = repmat(gt.Fe_used, nP, 1);
row.Fe_raw_g = repmat(gt.Fe_raw, nP, 1);
row.Fe3_g = repmat(gt.Fe3, nP, 1);
row.Mg_g = repmat(gt.Mg, nP, 1);
row.Mn_g = repmat(gt.Mn, nP, 1);
row.Ca_g = repmat(gt.Ca, nP, 1);
row.Si_g = repmat(gt.Si, nP, 1);
row.Al_g = repmat(gt.Al, nP, 1);

row.Fe2_bt = repmat(bt.Fe_used, nP, 1);
row.Fe_raw_bt = repmat(bt.Fe_raw, nP, 1);
row.Fe3_bt = repmat(bt.Fe3, nP, 1);
row.Mg_bt = repmat(bt.Mg, nP, 1);
row.Mn_bt = repmat(bt.Mn, nP, 1);
row.Ti_bt = repmat(bt.Ti, nP, 1);
row.Na_bt = repmat(bt.Na, nP, 1);
row.K_bt = repmat(bt.K, nP, 1);
row.Si_bt = repmat(bt.Si, nP, 1);
row.Al_bt = repmat(bt.Al, nP, 1);
row.F_bt = repmat(bt.F, nP, 1);

row.FeMg_garnet = repmat(FeMg_garnet_scalar, nP, 1);
row.FeMg_biotite = repmat(FeMg_biotite_scalar, nP, 1);
row.KD = repmat(KD_scalar, nP, 1);
row.lnKD = repmat(lnKD_scalar, nP, 1);
row.XCa_g = repmat(XCa_g_scalar, nP, 1);
row.XFe_g = repmat(XFe_g_scalar, nP, 1);
row.XMn_g = repmat(XMn_g_scalar, nP, 1);
row.n_bt_MgFe = repmat(n_bt_scalar, nP, 1);
row.XF_bt = repmat(XF_bt_scalar, nP, 1);
row.XF_over_n = repmat(XF_over_n_scalar, nP, 1);

row.Mg_number_garnet = repmat(Mg_number_gt_scalar, nP, 1);
row.Mg_number_biotite = repmat(Mg_number_bt_scalar, nP, 1);

row.T_K_standard = T_K_standard;
row.T_deg_standard = T_deg_standard;
row.T_K_Fcorrected = T_K_Fcorrected;
row.T_deg_Fcorrected = T_deg_Fcorrected;

row.denominator_standard = repmat(den_standard_scalar, nP, 1);
row.denominator_Fcorrected = repmat(den_Fcorrected_scalar, nP, 1);

row.is_positive_FeMg_garnet = repmat( ...
    isfinite(gt.Fe_used) && isfinite(gt.Mg) && ...
    gt.Fe_used > 0 && gt.Mg > 0, nP, 1);
row.is_positive_FeMg_biotite = repmat( ...
    isfinite(bt.Fe_used) && isfinite(bt.Mg) && ...
    bt.Fe_used > 0 && bt.Mg > 0, nP, 1);
row.is_positive_KD = repmat(isfinite(KD_scalar) && KD_scalar > 0, nP, 1);
row.is_positive_denominator_standard = repmat( ...
    isfinite(den_standard_scalar) && den_standard_scalar > 0, nP, 1);
row.is_positive_denominator_Fcorrected = repmat( ...
    isfinite(den_Fcorrected_scalar) && den_Fcorrected_scalar > 0, nP, 1);
row.is_XCa_g_reasonable = repmat(isfinite(XCa_g_scalar) && ...
    XCa_g_scalar > 0 && XCa_g_scalar < 0.3, nP, 1);
row.is_XFe_g_inel_range = repmat(isfinite(XFe_g_scalar) && ...
    XFe_g_scalar >= 0.4, nP, 1);
row.is_XMn_g_low = repmat(isfinite(XMn_g_scalar) && ...
    XMn_g_scalar < 0.1, nP, 1);
row.has_F_data = repmat(isfinite(bt.F) && bt.F > 0, nP, 1);

end

function gt = prepareGarnetRow(dataTable)
% prepareGarnetRow
% Extract one garnet row while retaining present NaN values.

if height(dataTable) ~= 1
    error('Garnet input must be a 1-row table.');
end

gt = struct();
gt.Fe_raw = getRequiredScalar(dataTable, 'Fe_cation_apfu', 'Garnet');
gt.Fe3 = getOptionalScalar(dataTable, 'Fe3_cation_apfu', 'Garnet', 0);
gt.Fe_used = gt.Fe_raw + gt.Fe3;
gt.Mg = getRequiredScalar(dataTable, 'Mg_cation_apfu', 'Garnet');
gt.Ca = getRequiredScalar(dataTable, 'Ca_cation_apfu', 'Garnet');
gt.Mn = getOptionalScalar(dataTable, 'Mn_cation_apfu', 'Garnet', 0);
gt.Si = getOptionalScalar(dataTable, 'Si_cation_apfu', 'Garnet', NaN);
gt.Al = getOptionalScalar(dataTable, 'Al_cation_apfu', 'Garnet', NaN);

end

function bt = prepareBiotiteRow(dataTable)
% prepareBiotiteRow
% Extract one biotite row while retaining present NaN values.

if height(dataTable) ~= 1
    error('Biotite input must be a 1-row table.');
end

bt = struct();
bt.Fe_raw = getRequiredScalar(dataTable, 'Fe_cation_apfu', 'Biotite');
bt.Fe3 = getOptionalScalar(dataTable, 'Fe3_cation_apfu', 'Biotite', 0);
bt.Fe_used = bt.Fe_raw + bt.Fe3;
bt.Mg = getRequiredScalar(dataTable, 'Mg_cation_apfu', 'Biotite');
bt.F = getOptionalScalar(dataTable, 'F_anion', 'Biotite', 0);

bt.Mn = getOptionalScalar(dataTable, 'Mn_cation_apfu', 'Biotite', NaN);
bt.Si = getOptionalScalar(dataTable, 'Si_cation_apfu', 'Biotite', NaN);
bt.Al = getOptionalScalar(dataTable, 'Al_cation_apfu', 'Biotite', NaN);
bt.Ti = getOptionalScalar(dataTable, 'Ti_cation_apfu', 'Biotite', NaN);
bt.Na = getOptionalScalar(dataTable, 'Na_cation_apfu', 'Biotite', NaN);
bt.K = getOptionalScalar(dataTable, 'K_cation_apfu', 'Biotite', NaN);

end

function value = getRequiredScalar(dataTable, variableName, mineralLabel)
% getRequiredScalar
% Return a required numeric scalar without changing NaN to zero.

value = dataTable.(variableName);
if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end

function value = getOptionalScalar(dataTable, variableName, ...
        mineralLabel, missingValue)
% getOptionalScalar
% Return a present optional scalar unchanged. Use missingValue only when the
% column itself is absent; a present NaN remains NaN.

if ismember(variableName, dataTable.Properties.VariableNames)
    value = dataTable.(variableName);
    if ~isnumeric(value) || ~isscalar(value)
        error('%s.%s must be a numeric scalar in the selected row.', ...
            mineralLabel, variableName);
    end
else
    value = missingValue;
end

end
