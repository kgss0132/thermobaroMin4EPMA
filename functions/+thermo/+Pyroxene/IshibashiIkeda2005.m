function results = IshibashiIkeda2005(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/IshibashiIkeda2005.m
% Tested with MATLAB R2024b
%
% Revised two-pyroxene thermometer
% Ishibashi, H. and Ikeda, T. (2005)
% Japanese Magazine of Mineralogical and Petrological Sciences, 34, 186-194
% DOI: https://doi.org/10.2465/gkk.34.186
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis and one
% Clinopyroxene analysis (selected independently from tables) and calculates
% temperature using the revised two-pyroxene thermometer of Ishibashi and
% Ikeda (2005).
%
% The function accepts either a scalar pressure or a pressure vector. It is
% therefore compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP. For each selected Opx-Cpx pair, the output table
% contains one row per input pressure value.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Opx-Cpx pair, stores each result block
% temporarily, and concatenates all blocks only once at the end.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Ishibashi and Ikeda (2005) evaluated and calibrated the revised
% two-pyroxene thermometer using published equilibrium experiments over the
% following ranges:
%
%   CMS system:
%     Temperature : 900-1500 degreeC
%     Pressure    : 15-60 kbar
%     Reproduction: mean DeltaT = +1 degreeC; 1 sigma = 21 degreeC
%
%   Multicomponent system:
%     Temperature : 900-1400 degreeC
%     Pressure    : 10-60 kbar
%     XFe in Opx and Cpx: 0.04-0.14
%     Opx Al2O3   : < 8 wt.%
%     Cpx Al2O3   : < 9 wt.%
%     Cpx Na2O    : < 3 wt.%
%     Reproduction: mean DeltaT = -5 degreeC; 1 sigma = 35 degreeC
%
% The experimental datasets and their P-T ranges are described on
% pp. 188-191 and in Table 1. The multicomponent data-screening criteria are
% described on pp. 188-191: experiments containing pigeonite were excluded,
% and only pyroxene analyses with 4.00 +/- 0.02 total cations per 6 oxygens
% were retained. The revised formulation is developed on pp. 192-193 and is
% given as equation (13); its performance is summarized in Table 2 and
% Figure 8.
%
% The thermometer should be applied to texturally and chemically equilibrated
% coexisting Opx-Cpx pairs. Pairing different generations, cores with
% unrelated rims, exsolution lamellae with unsuitable host compositions, or
% altered/metasomatized analyses may yield geologically meaningless results.
%
% The multicomponent calibration is most directly relevant to natural rocks.
% This implementation therefore issues non-stopping warnings when:
%   1) input pressure is outside 10-60 kbar, or
%   2) a finite calculated temperature is outside 900-1400 degreeC.
%
% Additional non-stopping composition warnings are issued when finite XFe
% values lie outside 0.04-0.14 or when the calculated total cations lie
% outside 3.98-4.02 per 6 oxygens. The Al2O3 and Na2O limits above cannot be
% checked directly from cation-apfu inputs alone.
%
% The CMS-system range extends to 1500 degreeC, but use of that upper limit
% should be restricted to genuinely simple CMS compositions. Application
% outside the stated P-T-composition ranges is extrapolative.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%   rawdata_struct.Cpx : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must contain
% normalized pyroxene cations on a 6-oxygen basis.
%
% Required variables in both Opx and Cpx tables:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu         % total Fe
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   Mn_cation_apfu
%   Ti_cation_apfu
%   Cr_cation_apfu
%
% Optional variable:
%   Fe3_cation_apfu
%
% If Fe3_cation_apfu is absent, Fe3+ is assumed to be zero. If the column is
% present and its selected value is NaN, the NaN is retained and propagated;
% it is not replaced by zero.
%
% All finite mineral-composition values used by the thermometer must be
% greater than or equal to zero. Finite negative values are rejected. NaN
% values are retained as missing values, propagated through the calculation,
% and reported using non-stopping fprintf warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Ishibashi and Ikeda (2005), equation (13):
%
%   lnK = [3867 - 9.686*(P_kbar - 20)] / T
%         + 5.518e-3*T + 2.44*XFe_opx - 12.48
%
% where
%   K       = aEn_cpx / aEn_opx
%   aEn     = XMg_M1 * XMg_M2
%   XFe_opx = Fe2_opx / (Fe2_opx + Mg_opx)
%   P_kbar  = pressure in kbar
%   T       = temperature in Kelvin
%
% Site allocation follows the Wood and Banno (1973)-style approximation:
%   T site:
%     AlIV = max(0, 2 - Si)
%     AlVI = max(0, Al_total - AlIV)
%
%   M1 preferential occupants:
%     AlVI, Cr, Ti, Fe3+
%
%   M2 preferential occupants:
%     Ca, Na, Mn
%
%   Remaining M1 and M2 capacities are filled by Mg and Fe2+ in proportion
%   to their available amounts.
%
% The published equation is rearranged to a quadratic in T:
%
%   B*T^2 + (C - 12.48 - lnK)*T + A = 0
%
% where
%   A = 3867 - 9.686*(P_kbar - 20)
%   B = 5.518e-3
%   C = 2.44*XFe_opx
%
% When two finite positive roots fall between 500 and 2500 K, this
% implementation retains the behavior of the original script and selects
% the root closest to 1273.15 K. When no physically screened real root is
% available, T remains NaN and a non-stopping warning is printed.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = IshibashiIkeda2005(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Opx and Cpx tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Opx-Cpx pair. NaN and Inf results are retained.
%

%% Input validation
% Basic argument checks prevent silent failures from missing inputs or
% invalid pressure vectors.
if nargin < 2
    error('IshibashiIkeda2005 requires (rawdata_struct, P_kbar).');
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
% Extract required input tables without modifying their contents.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end
if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

dataset_opx = rawdata_struct.Opx;
dataset_cpx = rawdata_struct.Cpx;

requiredVariables = {'Si_cation_apfu', 'Al_cation_apfu', ...
    'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu', ...
    'Na_cation_apfu', 'Mn_cation_apfu', 'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

validateRequiredVariables(dataset_opx, requiredVariables, 'Opx');
validateRequiredVariables(dataset_cpx, requiredVariables, 'Cpx');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a cell buffer and concatenate once at
% the end. This avoids repeated reallocation of the complete output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Multicomponent experimental evaluation limits from Ishibashi and Ikeda
% (2005), used for non-stopping range warnings.
calibrationT_min_degC = 900;
calibrationT_max_degC = 1400;
calibrationP_min_kbar = 10;
calibrationP_max_kbar = 60;
calibrationXFe_min = 0.04;
calibrationXFe_max = 0.14;
calibrationSumCat_min = 3.98;
calibrationSumCat_max = 4.02;

% Pressure is common to every selected pair in this function call, so the
% pressure warning is printed only once after the first calculation.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a dialog or chooses Finish.
disp('=== Step 3: Selecting a data code from the list (Opx) ===');

while true
    % ----- Opx selection -----
    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Opx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Opx selected: ' char(string(selectedCode_opx))]);

    % ----- Cpx selection -----
    disp('=== Step 4: Selecting a data code from the list (Cpx) ===');

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
    % Opx and Cpx are selected independently; row indices are not assumed to
    % correspond between the two input tables.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_opx = dataset_opx(selectedIdx_opx, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    % Identify NaN inputs without stopping the calculation. Values remain NaN
    % and propagate naturally through site allocation and thermometry.
    nanInputNames = findNaNInputs(selectedData_opx, selectedData_cpx, ...
        requiredVariables);

    % Reject finite negative values and Inf. Zero and NaN are allowed here;
    % zero may subsequently produce a non-finite result, which is retained.
    validateMineralInputs(selectedData_opx, selectedData_cpx, ...
        requiredVariables);

    row = calcTemp(selectedData_opx, selectedData_cpx, P_kbar);

    % Store selected identifiers once per pressure row for traceability.
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_opx', 'dataCode_cpx'}, 'Before', 1);

    % Store this result block. The cell buffer is doubled only when full,
    % avoiding repeated growth of the complete results table.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_cpx)) ': ' ...
            num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_cpx)) ': ' ...
            num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    % Warn once when any pressure point lies outside 10-60 kbar.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the multicomponent experimental ' ...
             'evaluation range of Ishibashi and Ikeda (2005): 10-60 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), numel(P_kbar), ...
            min(P_kbar), max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when finite calculated temperatures lie outside 900-1400 degreeC.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the multicomponent ' ...
             'experimental evaluation range of Ishibashi and Ikeda (2005): ' ...
             '900-1400 degreeC. %d of %d finite temperature point(s) are ' ...
             'outside the range; calculated finite range = %.4g-%.4g degreeC ' ...
             'for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_opx)), char(string(selectedCode_cpx)));
    end

    % Warn when XFe values are outside the multicomponent evaluation range.
    finiteXFe_opx = isfinite(row.XFe_opx);
    finiteXFe_cpx = isfinite(row.XFe_cpx);
    xFeOutside_opx = finiteXFe_opx & ...
        (row.XFe_opx < calibrationXFe_min | ...
         row.XFe_opx > calibrationXFe_max);
    xFeOutside_cpx = finiteXFe_cpx & ...
        (row.XFe_cpx < calibrationXFe_min | ...
         row.XFe_cpx > calibrationXFe_max);

    if any(xFeOutside_opx) || any(xFeOutside_cpx)
        fprintf(2, ...
            ['WARNING: Pyroxene XFe is outside the multicomponent experimental ' ...
             'evaluation range of Ishibashi and Ikeda (2005): 0.04-0.14. ' ...
             'Opx XFe = %.6g; Cpx XFe = %.6g for %s & %s.\n'], ...
            row.XFe_opx(1), row.XFe_cpx(1), ...
            char(string(selectedCode_opx)), char(string(selectedCode_cpx)));
    end

    % Warn when total cations are outside 4.00 +/- 0.02 per 6 oxygens.
    sumCatOutside_opx = isfinite(row.SumCat_opx) & ...
        (row.SumCat_opx < calibrationSumCat_min | ...
         row.SumCat_opx > calibrationSumCat_max);
    sumCatOutside_cpx = isfinite(row.SumCat_cpx) & ...
        (row.SumCat_cpx < calibrationSumCat_min | ...
         row.SumCat_cpx > calibrationSumCat_max);

    if any(sumCatOutside_opx) || any(sumCatOutside_cpx)
        fprintf(2, ...
            ['WARNING: Pyroxene total cations are outside the data-screening ' ...
             'criterion used by Ishibashi and Ikeda (2005): 3.98-4.02 ' ...
             'cations per 6 oxygens. Opx total = %.6g; Cpx total = %.6g ' ...
             'for %s & %s.\n'], ...
            row.SumCat_opx(1), row.SumCat_cpx(1), ...
            char(string(selectedCode_opx)), char(string(selectedCode_cpx)));
    end

    % Report every selected thermometer input whose value was NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN values were retained and propagated; they were not ' ...
             'replaced by zero.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report non-finite output temperatures.
    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), numel(row.T_deg), ...
            sum(isnan(row.T_deg)), sum(isinf(row.T_deg)));
    end

    % Report pressure points for which the quadratic has no retained real
    % root, without stopping calculations for other pressure points.
    noSelectedRoot = row.root_selection_code == 0;
    if any(noSelectedRoot)
        fprintf(2, ...
            ['WARNING: No finite physically screened real temperature root was ' ...
             'retained for %s & %s at %d of %d pressure point(s). ' ...
             'The corresponding T values remain NaN.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            sum(noSelectedRoot), numel(noSelectedRoot));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another independently selected pair.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'IshibashiIkeda2005', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all result blocks once. Return an empty table if the user did
% not complete any calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(tbl, requiredVariables, mineralLabel)
% validateRequiredVariables
% Confirm that every required thermometer input column is present.

missingVariables = requiredVariables(~ismember(requiredVariables, ...
    tbl.Properties.VariableNames));

if ~isempty(missingVariables)
    error('%s table is missing required variable(s): %s', ...
        mineralLabel, char(strjoin(string(missingVariables), ', ')));
end

end

function nanInputNames = findNaNInputs(data_opx, data_cpx, requiredVariables)
% findNaNInputs
% Return the names of selected thermometer input variables containing NaN.
% The returned list is used only for fprintf warnings; NaN does not stop the
% calculation and is never replaced by zero.

maxEntries = 2 * (numel(requiredVariables) + 1);
nanBuffer = strings(maxEntries, 1);
nNan = 0;

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};

    value_opx = data_opx.(variableName);
    if any(isnan(value_opx(:)))
        nNan = nNan + 1;
        nanBuffer(nNan) = "Opx." + string(variableName);
    end

    value_cpx = data_cpx.(variableName);
    if any(isnan(value_cpx(:)))
        nNan = nNan + 1;
        nanBuffer(nNan) = "Cpx." + string(variableName);
    end
end

% Fe3+ is optional. Check it only when the column is actually present.
if ismember('Fe3_cation_apfu', data_opx.Properties.VariableNames)
    value_opx = data_opx.Fe3_cation_apfu;
    if any(isnan(value_opx(:)))
        nNan = nNan + 1;
        nanBuffer(nNan) = "Opx.Fe3_cation_apfu";
    end
end

if ismember('Fe3_cation_apfu', data_cpx.Properties.VariableNames)
    value_cpx = data_cpx.Fe3_cation_apfu;
    if any(isnan(value_cpx(:)))
        nNan = nNan + 1;
        nanBuffer(nNan) = "Cpx.Fe3_cation_apfu";
    end
end

nanInputNames = nanBuffer(1:nNan);

end

function validateMineralInputs(data_opx, data_cpx, requiredVariables)
% validateMineralInputs
% Reject finite negative values and Inf in all thermometer inputs. Zero is
% allowed. NaN is intentionally allowed so it can propagate into the result.

maxEntries = 2 * (numel(requiredVariables) + 1);
invalidBuffer = strings(maxEntries, 1);
nInvalid = 0;

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};

    value_opx = data_opx.(variableName);
    if ~isscalar(value_opx)
        error('Opx variable %s must be scalar in a selected 1-row table.', ...
            variableName);
    end
    if isinf(value_opx) || (isfinite(value_opx) && value_opx < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Opx." + string(variableName);
    end

    value_cpx = data_cpx.(variableName);
    if ~isscalar(value_cpx)
        error('Cpx variable %s must be scalar in a selected 1-row table.', ...
            variableName);
    end
    if isinf(value_cpx) || (isfinite(value_cpx) && value_cpx < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Cpx." + string(variableName);
    end
end

% Validate optional Fe3+ values only if their columns are present.
if ismember('Fe3_cation_apfu', data_opx.Properties.VariableNames)
    value_opx = data_opx.Fe3_cation_apfu;
    if ~isscalar(value_opx)
        error('Opx variable Fe3_cation_apfu must be scalar in a selected 1-row table.');
    end
    if isinf(value_opx) || (isfinite(value_opx) && value_opx < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Opx.Fe3_cation_apfu";
    end
end

if ismember('Fe3_cation_apfu', data_cpx.Properties.VariableNames)
    value_cpx = data_cpx.Fe3_cation_apfu;
    if ~isscalar(value_cpx)
        error('Cpx variable Fe3_cation_apfu must be scalar in a selected 1-row table.');
    end
    if isinf(value_cpx) || (isfinite(value_cpx) && value_cpx < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Cpx.Fe3_cation_apfu";
    end
end

if nInvalid > 0
    invalidNames = invalidBuffer(1:nInvalid);
    error(['IshibashiIkeda2005: finite thermometer inputs must be ' ...
           'greater than or equal to zero, and Inf is not permitted. ' ...
           'Invalid value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_cpx, P_kbar)
% calcTemp
% Compute temperatures for one selected Opx-Cpx pair over a scalar or vector
% of pressures. The returned table contains one row per pressure value.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

% --- Prepare pyroxene cation rows ---
opx = preparePyroxeneRow(data_opx, 'Opx');
cpx = preparePyroxeneRow(data_cpx, 'Cpx');

% --- Site allocation ---
site_opx = calcSiteFractions(opx, 'Opx');
site_cpx = calcSiteFractions(cpx, 'Cpx');

% --- Idealized activities of the enstatite component ---
aEn_opx_scalar = site_opx.XMg_M1 .* site_opx.XMg_M2;
aEn_cpx_scalar = site_cpx.XMg_M1 .* site_cpx.XMg_M2;

% Direct arithmetic is intentional. Zero and NaN values may produce NaN or
% Inf, which are retained rather than converted or used to stop the run.
K_scalar = aEn_cpx_scalar ./ aEn_opx_scalar;
lnK_scalar = log(K_scalar);

% --- Fe numbers ---
XFe_opx_scalar = opx.Fe2 ./ (opx.Fe2 + opx.Mg);
XFe_cpx_scalar = cpx.Fe2 ./ (cpx.Fe2 + cpx.Mg);

% --- Total cations per 6 oxygens ---
% Fe_total is counted once; Fe2 and Fe3 are not added separately.
SumCat_opx_scalar = opx.Si + opx.Al + opx.Fe_total + opx.Mg + ...
    opx.Ca + opx.Na + opx.Mn + opx.Ti + opx.Cr;
SumCat_cpx_scalar = cpx.Si + cpx.Al + cpx.Fe_total + cpx.Mg + ...
    cpx.Ca + cpx.Na + cpx.Mn + cpx.Ti + cpx.Cr;

% Replicate composition-dependent scalars to match the pressure-vector size.
lnK = repmat(lnK_scalar, nP, 1);
XFe_opx = repmat(XFe_opx_scalar, nP, 1);
XFe_cpx = repmat(XFe_cpx_scalar, nP, 1);

% --- Ishibashi and Ikeda (2005) equation coefficients ---
A = 3867 - 9.686 .* (P_kbar - 20);
B = 5.518e-3;
C = 2.44 .* XFe_opx;

quad_a = repmat(B, nP, 1);
quad_b = C - 12.48 - lnK;
quad_c = A;

discriminant = quad_b.^2 - 4 .* quad_a .* quad_c;

% Avoid complex roots. Negative or non-finite discriminants are represented
% by NaN roots and are reported later by non-stopping warnings.
sqrtDiscriminant = nan(nP, 1);
realDiscriminant = isfinite(discriminant) & discriminant >= 0;
sqrtDiscriminant(realDiscriminant) = sqrt(discriminant(realDiscriminant));

T1 = (-quad_b + sqrtDiscriminant) ./ (2 .* quad_a);
T2 = (-quad_b - sqrtDiscriminant) ./ (2 .* quad_a);

% Select one root independently for every pressure point.
% code 0: no retained root; code 1: root 1; code 2: root 2.
T_K = nan(nP, 1);
rootSelectionCode = zeros(nP, 1);

validRoot1 = isfinite(T1) & T1 > 500 & T1 < 2500;
validRoot2 = isfinite(T2) & T2 > 500 & T2 < 2500;

onlyRoot1 = validRoot1 & ~validRoot2;
onlyRoot2 = validRoot2 & ~validRoot1;
bothRoots = validRoot1 & validRoot2;

T_K(onlyRoot1) = T1(onlyRoot1);
rootSelectionCode(onlyRoot1) = 1;

T_K(onlyRoot2) = T2(onlyRoot2);
rootSelectionCode(onlyRoot2) = 2;

if any(bothRoots)
    distanceRoot1 = abs(T1(bothRoots) - 1273.15);
    distanceRoot2 = abs(T2(bothRoots) - 1273.15);
    chooseRoot1 = distanceRoot1 <= distanceRoot2;

    bothIndices = find(bothRoots);
    root1Indices = bothIndices(chooseRoot1);
    root2Indices = bothIndices(~chooseRoot1);

    T_K(root1Indices) = T1(root1Indices);
    rootSelectionCode(root1Indices) = 1;

    T_K(root2Indices) = T2(root2Indices);
    rootSelectionCode(root2Indices) = 2;
end

T_deg = T_K - 273.15;

% --- Residual check against equation (13) ---
eq_residual = lnK - (A ./ T_K + B .* T_K + C - 12.48);

% --- Pack outputs ---
row = table();
row.P_kbar = P_kbar;

row.Si_opx = repmat(opx.Si, nP, 1);
row.Al_opx = repmat(opx.Al, nP, 1);
row.Fe_total_opx = repmat(opx.Fe_total, nP, 1);
row.Fe2_opx = repmat(opx.Fe2, nP, 1);
row.Fe3_opx = repmat(opx.Fe3, nP, 1);
row.Mg_opx = repmat(opx.Mg, nP, 1);
row.Ca_opx = repmat(opx.Ca, nP, 1);
row.Na_opx = repmat(opx.Na, nP, 1);
row.Mn_opx = repmat(opx.Mn, nP, 1);
row.Ti_opx = repmat(opx.Ti, nP, 1);
row.Cr_opx = repmat(opx.Cr, nP, 1);
row.SumCat_opx = repmat(SumCat_opx_scalar, nP, 1);

row.Si_cpx = repmat(cpx.Si, nP, 1);
row.Al_cpx = repmat(cpx.Al, nP, 1);
row.Fe_total_cpx = repmat(cpx.Fe_total, nP, 1);
row.Fe2_cpx = repmat(cpx.Fe2, nP, 1);
row.Fe3_cpx = repmat(cpx.Fe3, nP, 1);
row.Mg_cpx = repmat(cpx.Mg, nP, 1);
row.Ca_cpx = repmat(cpx.Ca, nP, 1);
row.Na_cpx = repmat(cpx.Na, nP, 1);
row.Mn_cpx = repmat(cpx.Mn, nP, 1);
row.Ti_cpx = repmat(cpx.Ti, nP, 1);
row.Cr_cpx = repmat(cpx.Cr, nP, 1);
row.SumCat_cpx = repmat(SumCat_cpx_scalar, nP, 1);

row.AlIV_opx = repmat(site_opx.AlIV, nP, 1);
row.AlVI_opx = repmat(site_opx.AlVI, nP, 1);
row.XMg_M1_opx = repmat(site_opx.XMg_M1, nP, 1);
row.XMg_M2_opx = repmat(site_opx.XMg_M2, nP, 1);
row.M1_total_opx = repmat(site_opx.M1_total, nP, 1);
row.M2_total_opx = repmat(site_opx.M2_total, nP, 1);

row.AlIV_cpx = repmat(site_cpx.AlIV, nP, 1);
row.AlVI_cpx = repmat(site_cpx.AlVI, nP, 1);
row.XMg_M1_cpx = repmat(site_cpx.XMg_M1, nP, 1);
row.XMg_M2_cpx = repmat(site_cpx.XMg_M2, nP, 1);
row.M1_total_cpx = repmat(site_cpx.M1_total, nP, 1);
row.M2_total_cpx = repmat(site_cpx.M2_total, nP, 1);

row.aEn_opx = repmat(aEn_opx_scalar, nP, 1);
row.aEn_cpx = repmat(aEn_cpx_scalar, nP, 1);
row.K_En = repmat(K_scalar, nP, 1);
row.lnK = lnK;
row.XFe_opx = XFe_opx;
row.XFe_cpx = XFe_cpx;

row.coeff_A = A;
row.coeff_B = quad_a;
row.coeff_C = C;
row.quad_b = quad_b;
row.discriminant = discriminant;
row.T_K_root1 = T1;
row.T_K_root2 = T2;
row.root_selection_code = rootSelectionCode;
row.eq_residual = eq_residual;

row.T_K = T_K;
row.T_deg = T_deg;

end

function px = preparePyroxeneRow(data_px, mineralLabel)
% preparePyroxeneRow
% Extract one selected pyroxene row. Required values are read directly.
% Missing Fe3_cation_apfu is treated as zero; a present NaN is preserved.

if height(data_px) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

px = struct();
px.Si = getRequiredVar(data_px, 'Si_cation_apfu', mineralLabel);
px.Al = getRequiredVar(data_px, 'Al_cation_apfu', mineralLabel);
px.Fe_total = getRequiredVar(data_px, 'Fe_cation_apfu', mineralLabel);
px.Mg = getRequiredVar(data_px, 'Mg_cation_apfu', mineralLabel);
px.Ca = getRequiredVar(data_px, 'Ca_cation_apfu', mineralLabel);
px.Na = getRequiredVar(data_px, 'Na_cation_apfu', mineralLabel);
px.Mn = getRequiredVar(data_px, 'Mn_cation_apfu', mineralLabel);
px.Ti = getRequiredVar(data_px, 'Ti_cation_apfu', mineralLabel);
px.Cr = getRequiredVar(data_px, 'Cr_cation_apfu', mineralLabel);
px.Fe3 = getOptionalFe3(data_px);

% Fe3 > total Fe is physically inconsistent when both values are finite.
if isfinite(px.Fe3) && isfinite(px.Fe_total) && ...
        px.Fe3 > px.Fe_total + 1e-12
    error('%s has Fe3_cation_apfu > Fe_cation_apfu.', mineralLabel);
end

px.Fe2 = px.Fe_total - px.Fe3;
if isfinite(px.Fe2) && px.Fe2 < -1e-12
    error('%s has negative Fe2 after Fe_total - Fe3.', mineralLabel);
end
if isfinite(px.Fe2) && px.Fe2 < 0
    px.Fe2 = 0;
end

end

function site = calcSiteFractions(px, mineralLabel)
% calcSiteFractions
% Calculate Wood and Banno-style pyroxene site fractions while preserving
% NaN values. Finite site occupancies exceeding capacity are rejected.

site = struct();

% --- T site ---
if isnan(px.Si)
    site.AlIV = NaN;
else
    site.AlIV = max(0, 2 - px.Si);
end

if isnan(px.Al) || isnan(site.AlIV)
    site.AlVI = NaN;
else
    site.AlVI = max(0, px.Al - site.AlIV);
end

% --- Fixed M1 and M2 occupants ---
M1_fixed = site.AlVI + px.Cr + px.Ti + px.Fe3;
M2_fixed = px.Ca + px.Na + px.Mn;

if isfinite(M1_fixed) && M1_fixed > 1 + 1e-8
    error('%s calculated M1 fixed occupancy exceeds 1. Check cation normalization.', ...
        mineralLabel);
end
if isfinite(M2_fixed) && M2_fixed > 1 + 1e-8
    error('%s calculated M2 fixed occupancy exceeds 1. Check cation normalization.', ...
        mineralLabel);
end

if isnan(M1_fixed)
    M1_remaining = NaN;
else
    M1_remaining = max(0, 1 - M1_fixed);
end

if isnan(M2_fixed)
    M2_remaining = NaN;
else
    M2_remaining = max(0, 1 - M2_fixed);
end

MgFe_total = px.Mg + px.Fe2;
Mg_fraction = px.Mg ./ MgFe_total;
Fe_fraction = px.Fe2 ./ MgFe_total;

Mg_M1 = M1_remaining .* Mg_fraction;
Fe2_M1 = M1_remaining .* Fe_fraction;
Mg_M2 = M2_remaining .* Mg_fraction;
Fe2_M2 = M2_remaining .* Fe_fraction;

site.M1_total = M1_fixed + Mg_M1 + Fe2_M1;
site.M2_total = M2_fixed + Mg_M2 + Fe2_M2;

site.XMg_M1 = Mg_M1 ./ site.M1_total;
site.XMg_M2 = Mg_M2 ./ site.M2_total;

end

function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read one required scalar value without modifying NaN.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);
if ~isscalar(value)
    error('%s variable %s must be scalar in a selected 1-row table.', ...
        mineralLabel, varName);
end

end

function value = getOptionalFe3(tbl)
% getOptionalFe3
% Return Fe3_cation_apfu when present, preserving NaN. If the column is
% absent, return zero as the documented Fe3+ assumption.

if ismember('Fe3_cation_apfu', tbl.Properties.VariableNames)
    value = tbl.Fe3_cation_apfu;
    if ~isscalar(value)
        error('Variable Fe3_cation_apfu must be scalar in a selected 1-row table.');
    end
else
    value = 0;
end

end
