function results = LesselPutirka2015twoPx(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/LesselPutirka2015twoPx.m
% Tested with MATLAB R2024b
%
% Two-pyroxene thermometer for martian igneous compositions
% Lessel, J. and Putirka, K. (2015), equation (7)
% American Mineralogist, 100, 2163-2171
% DOI: https://doi.org/10.2138/am-2015-4732
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis and one
% Clinopyroxene analysis and calculates temperature using equation (7) of
% Lessel and Putirka (2015).
%
% The function accepts either a scalar pressure or a pressure vector and is
% therefore compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP. For each selected Opx-Cpx pair, the output table
% contains one row per input pressure value.
%
% The function is designed for repeated calculations. Each selected-pair
% result is stored temporarily in a preallocated cell buffer, and all result
% blocks are concatenated only once after the interactive loop finishes.
%
% This is a mineral-mineral thermometer and does not use a Liquid dataset.
% Therefore Liquid cationTotal_liq, F, and Cl handling is not applicable to
% this function.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Lessel and Putirka (2015) developed this thermometer from experiments on
% martian meteorites and martian-analog bulk compositions. Martian basaltic
% systems are generally richer in FeO and poorer in Al2O3 than common
% terrestrial basaltic systems. Application to ordinary terrestrial
% pyroxenes is therefore extrapolative unless their compositions closely
% resemble the calibration data.
%
% The complete experimental compilation used in the paper spans:
%
%   Pressure     : approximately 0.0001-2.3 GPa
%                  (approximately 0.001-23 kbar)
%   Temperature  : approximately 950-1540 degreeC
%   Liquid SiO2  : 40.2-66.16 wt.%
%   Liquid MgO   : 0.62-24.52 wt.%
%   Liquid FeO   : 2.80-30.2 wt.%
%   Liquid Al2O3 : 2.97-20.5 wt.%
%   Total alkalis: 0.19-6.77 wt.%
%
% These overall limits are reported in the Methodology section on
% pp. 2164-2165. They combine several mineral-liquid and mineral-mineral
% models and are not strict rectangular limits for equation (7) alone.
%
% From the source experiments in Table 1 that contain both Opx and Cpx, the
% approximate phase-specific source-study envelope is:
%
%   Pressure    : approximately 0.0001-2.3 GPa
%   Temperature : approximately 950-1440 degreeC
%
% This source-study envelope is used below for non-stopping warnings. It is
% not presented by the authors as a strict rectangular calibration boundary
% for the final equation (7) regression subset.
%
% Equation (7) was calibrated with 24 data and reproduced temperature with
% R2 = 0.99 and RMSE = 12 K. An independent test set of 5 data was predicted
% with R2 = 0.95 and RMSE = 37 K (Table 2 on p. 2164; equation (7) on
% p. 2166; Figure 8 on p. 2167). The test-set RMSE of approximately 37 K is
% the more conservative precision estimate for application.
%
% Pyroxene components are calculated from cations normalized to 6 oxygens,
% following the Putirka (2008)-style normative scheme used by Lessel and
% Putirka (2015) (pp. 2164-2166).
%
% The selected Opx and Cpx must represent an equilibrated pair. Lessel and
% Putirka (2015) discuss the Cpx-Opx Fe-Mg exchange coefficient on p. 2168.
% They do not recommend their temperature-dependent KD relation as an
% independent weak thermometer; instead, it may be used to assess whether a
% calculated geothermometric result is plausible. They also give a
% composition-dependent, temperature-independent check based on Cpx CaO
% (equations (14)-(16), p. 2168).
%
% Pairing different crystal generations, unrelated cores and rims,
% exsolution-modified analyses, altered pyroxenes, or pyroxenes that did not
% equilibrate with one another may yield geologically meaningless results.
%
% Equation (7) contains reciprocal terms. The following calculated
% quantities must therefore be finite and strictly greater than zero:
%
%   Na_Opx
%   XCaTs_Cpx
%   XFmAl2SiO6_Opx
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside the approximate two-pyroxene source-study
%      envelope of 0.0001-2.3 GPa, or
%   2) a finite calculated temperature is outside the approximate
%      two-pyroxene source-study envelope of 950-1440 degreeC.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%   rawdata_struct.Cpx : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must contain
% pyroxene cations normalized to a 6-oxygen basis.
%
% Required variables in both Opx and Cpx tables:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu
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
% The legacy input convention of the original script is retained:
%   - If Fe3_cation_apfu is absent, Fe_cation_apfu is used as FeTot.
%   - If Fe3_cation_apfu is present, Fe_cation_apfu is treated as Fe2+ and
%     FeTot is calculated as Fe_cation_apfu + Fe3_cation_apfu.
%
% For equation (7), Cpx Fe3+ is then recalculated using the Papike et al.
% (1974)-style stoichiometric estimate specified by Lessel and Putirka
% (2015). Users must ensure that their Fe_cation_apfu and Fe3_cation_apfu
% columns follow the convention above; otherwise Fe may be double-counted.
%
% All finite mineral-composition inputs must be greater than or equal to
% zero. Finite negative values and Inf are rejected. NaN values are retained
% as missing values, propagated through component and temperature
% calculations, and reported using non-stopping fprintf warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Lessel and Putirka (2015), equation (7):
%
%   1 / T(K) =
%       6.644e-4
%     - 2.757e-5 * P(GPa)
%     + 1.499e-3 * Na_Opx
%     - 1.640e-4 * Fe2_Cpx
%     + 6.664e-5 * (XDiHd_Cpx)^2
%     + 2.611e-4 * (FeTot_Opx)^2
%     - 6.602e-8 / XCaTs_Cpx
%     + 6.869e-8 / Na_Opx
%     - 1.166e-8 / XFmAl2SiO6_Opx
%
% where:
%   P(GPa)          = P_kbar / 10
%   Fe2_Cpx         = FeTot_Cpx - Fe3_calc_Cpx
%   Fe3_calc_Cpx    = AlIV_Cpx + Na_Cpx - AlVI_Cpx
%                     - Cr_Cpx - 2*Ti_Cpx
%   XJd_Cpx         = min(AlVI_Cpx, Na_Cpx)
%   XCaTs_Cpx       = AlVI_Cpx - XJd_Cpx
%   XCaTi_Cpx       = max[(AlIV_Cpx - XCaTs_Cpx)/2, 0]
%   XCrCaTs_Cpx     = Cr_Cpx / 2
%   XDiHd_Cpx       = Ca_Cpx - XCaTi_Cpx - XCaTs_Cpx - XCrCaTs_Cpx
%   AlIV_Opx        = max(2 - Si_Opx, 0)
%   AlVI_Opx        = max(Al_Opx - AlIV_Opx, 0)
%   XFmAl2SiO6_Opx  = min(AlIV_Opx, AlVI_Opx)
%
% Derived normative components that are slightly negative because of
% analytical or rounding effects are clamped to zero while preserving NaN.
% Raw negative inputs are never clamped and instead cause an error.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = LesselPutirka2015twoPx(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Opx and Cpx tables described above
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Opx-Cpx pair. NaN and Inf results are retained.
%

%% Input validation
if nargin < 2
    error('LesselPutirka2015twoPx requires (rawdata_struct, P_kbar).');
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
disp('=== Step 1: Preparing cation datasets ===');

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

disp('=== Preparing cation datasets has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a cell buffer and concatenate once at
% the end. This avoids resizing and copying the complete output table on
% every loop iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max(16, height(dataset_opx) * height(dataset_cpx));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Approximate two-pyroxene source-study envelope inferred from Table 1 of
% Lessel and Putirka (2015). These are warning limits rather than strict
% rectangular boundaries explicitly defined for equation (7).
applicationP_min_GPa = 0.0001;
applicationP_max_GPa = 2.3;
applicationT_min_degC = 950;
applicationT_max_degC = 1440;

P_GPa_input = P_kbar ./ 10;
pressureOutsideRange = P_GPa_input < applicationP_min_GPa | ...
    P_GPa_input > applicationP_max_GPa;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop and calculation
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
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_opx = dataset_opx(selectedIdx_opx, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    % Identify NaN values without stopping the calculation. NaN values are
    % retained and propagated through the calculation.
    nanInputNames = findNaNInputs(selectedData_opx, selectedData_cpx, ...
        requiredVariables);

    % Reject finite negative values and Inf. Zero and NaN are allowed here;
    % reciprocal-domain failures are handled as non-stopping NaN results.
    validateInputValues(selectedData_opx, selectedData_cpx, ...
        requiredVariables);

    row = calcTemp(selectedData_opx, selectedData_cpx, P_kbar);

    % Store selected identifiers once per pressure row for traceability.
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_opx', 'dataCode_cpx'}, 'Before', 1);

    % Store this result block. The buffer is enlarged only if the user makes
    % more selections than the preallocated number of possible table pairs.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_cpx)) ': ' ...
            num2str(row.TEq7_C) ' degreeC']);
    else
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_cpx)) ': ' ...
            num2str(row.TEq7_C(1)) ' to ' ...
            num2str(row.TEq7_C(end)) ' degreeC']);
    end

    % Warn once when input pressure lies outside the approximate source-study
    % envelope. The calculation continues for all pressure points.
    if any(pressureOutsideRange) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate two-pyroxene ' ...
             'source-study envelope of Lessel and Putirka (2015): ' ...
             '0.0001-2.3 GPa (0.001-23 kbar). %d of %d pressure point(s) ' ...
             'are outside the range; input range = %.6g-%.6g GPa.\n'], ...
            sum(pressureOutsideRange), numel(P_GPa_input), ...
            min(P_GPa_input), max(P_GPa_input));
        pressureWarningIssued = true;
    end

    % Warn when finite calculated temperatures lie outside the approximate
    % two-pyroxene source-study envelope.
    finiteTemperature = isfinite(row.TEq7_C);
    temperatureOutsideRange = finiteTemperature & ...
        (row.TEq7_C < applicationT_min_degC | ...
         row.TEq7_C > applicationT_max_degC);

    if any(temperatureOutsideRange)
        finiteValues = row.TEq7_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximate ' ...
             'two-pyroxene source-study envelope of Lessel and Putirka ' ...
             '(2015): 950-1440 degreeC. %d of %d finite temperature ' ...
             'point(s) are outside the range; calculated finite range = ' ...
             '%.6g-%.6g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideRange), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)));
    end

    % Report every selected thermometer input whose value was NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN values were retained and propagated; they were ' ...
             'not replaced by zero.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report non-finite output temperatures.
    invalidTemperature = ~isfinite(row.TEq7_C);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for ' ...
             '%s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), numel(row.TEq7_C), ...
            sum(isnan(row.TEq7_C)), sum(isinf(row.TEq7_C)));
    end

    % Report equation-domain failures separately so the cause of a NaN is
    % easier to identify.
    invalidEquationDomain = ~row.equation_domain_valid;
    if any(invalidEquationDomain)
        fprintf(2, ...
            ['WARNING: Lessel and Putirka (2015) equation (7) was outside ' ...
             'its mathematical domain for %s & %s at %d of %d pressure ' ...
             'point(s). Na_Opx, XCaTs_Cpx, and XFmAl2SiO6_Opx must be ' ...
             'finite and > 0; other required component terms and 1/T must ' ...
             'be finite, and 1/T must be > 0. Corresponding T values ' ...
             'remain NaN.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            sum(invalidEquationDomain), numel(invalidEquationDomain));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another independently selected pair.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'LesselPutirka2015twoPx', ...
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
% Confirm that every required thermometer-input column is present.

missingVariables = requiredVariables(~ismember(requiredVariables, ...
    tbl.Properties.VariableNames));

if ~isempty(missingVariables)
    error('%s table is missing required variable(s): %s', ...
        mineralLabel, char(strjoin(string(missingVariables), ', ')));
end

end

function nanInputNames = findNaNInputs(data_opx, data_cpx, requiredVariables)
% findNaNInputs
% Return names of selected thermometer inputs containing NaN. This function
% only prepares an fprintf warning; it does not replace NaN or stop the run.

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

% Fe3+ is optional and is checked only when its column is present.
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

function validateInputValues(data_opx, data_cpx, requiredVariables)
% validateInputValues
% Reject finite negative values and Inf in all raw thermometer inputs. Zero
% and NaN are allowed so domain failures can be returned as non-stopping NaN.

maxEntries = 2 * (numel(requiredVariables) + 1);
invalidBuffer = strings(maxEntries, 1);
nInvalid = 0;

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};

    value_opx = data_opx.(variableName);
    if ~isnumeric(value_opx) || ~isscalar(value_opx) || ~isreal(value_opx)
        error('Opx variable %s must be a real numeric scalar in a selected 1-row table.', ...
            variableName);
    end
    if isinf(value_opx) || (isfinite(value_opx) && value_opx < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Opx." + string(variableName);
    end

    value_cpx = data_cpx.(variableName);
    if ~isnumeric(value_cpx) || ~isscalar(value_cpx) || ~isreal(value_cpx)
        error('Cpx variable %s must be a real numeric scalar in a selected 1-row table.', ...
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
    if ~isnumeric(value_opx) || ~isscalar(value_opx) || ~isreal(value_opx)
        error('Opx variable Fe3_cation_apfu must be a real numeric scalar.');
    end
    if isinf(value_opx) || (isfinite(value_opx) && value_opx < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Opx.Fe3_cation_apfu";
    end
end

if ismember('Fe3_cation_apfu', data_cpx.Properties.VariableNames)
    value_cpx = data_cpx.Fe3_cation_apfu;
    if ~isnumeric(value_cpx) || ~isscalar(value_cpx) || ~isreal(value_cpx)
        error('Cpx variable Fe3_cation_apfu must be a real numeric scalar.');
    end
    if isinf(value_cpx) || (isfinite(value_cpx) && value_cpx < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Cpx.Fe3_cation_apfu";
    end
end

if nInvalid > 0
    invalidNames = invalidBuffer(1:nInvalid);
    error(['LesselPutirka2015twoPx: finite mineral inputs must be ' ...
           'greater than or equal to zero, and Inf is not permitted. ' ...
           'Invalid value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_cpx, P_kbar)
% calcTemp
% Calculate Lessel and Putirka (2015) equation (7) for one selected Opx-Cpx
% pair over a scalar or vector of pressures. One output row is returned for
% each pressure value.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

% --- Prepare pyroxene cation rows ---
opx = preparePyroxeneRow(data_opx, 'Opx');
cpx = preparePyroxeneRow(data_cpx, 'Cpx');

% --- Calculate Opx and Cpx component terms ---
opxTerms = calcOpxTerms(opx);
cpxTerms = calcCpxTerms(cpx);

% --- Equation (7) domain screening ---
baseDomainValid = ...
    isfinite(opx.Na) && opx.Na > 0 && ...
    isfinite(cpxTerms.XCaTs) && cpxTerms.XCaTs > 0 && ...
    isfinite(opxTerms.XFmAl2SiO6) && opxTerms.XFmAl2SiO6 > 0 && ...
    isfinite(cpxTerms.Fe2) && cpxTerms.Fe2 >= 0 && ...
    isfinite(cpxTerms.XDiHd) && cpxTerms.XDiHd >= 0 && ...
    isfinite(opx.FeTot) && opx.FeTot >= 0;

invT = NaN(nP, 1);
if baseDomainValid
    invT = ...
        6.644e-4 ...
        - 2.757e-5 .* P_GPa ...
        + 1.499e-3 .* opx.Na ...
        - 1.640e-4 .* cpxTerms.Fe2 ...
        + 6.664e-5 .* (cpxTerms.XDiHd .^ 2) ...
        + 2.611e-4 .* (opx.FeTot .^ 2) ...
        - 6.602e-8 ./ cpxTerms.XCaTs ...
        + 6.869e-8 ./ opx.Na ...
        - 1.166e-8 ./ opxTerms.XFmAl2SiO6;
end

validInvT = isfinite(invT) & invT > 0;
equationDomainValid = baseDomainValid & validInvT;

TEq7_K = NaN(nP, 1);
TEq7_C = NaN(nP, 1);
TEq7_K(validInvT) = 1 ./ invT(validInvT);
TEq7_C(validInvT) = TEq7_K(validInvT) - 273.15;

% --- Pack outputs ---
row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

row.Si_opx = repmat(opx.Si, nP, 1);
row.Al_opx = repmat(opx.Al, nP, 1);
row.Fe2_opx = repmat(opx.Fe2, nP, 1);
row.Fe3_opx = repmat(opx.Fe3, nP, 1);
row.FeTot_opx = repmat(opx.FeTot, nP, 1);
row.Mg_opx = repmat(opx.Mg, nP, 1);
row.Ca_opx = repmat(opx.Ca, nP, 1);
row.Na_opx = repmat(opx.Na, nP, 1);
row.Mn_opx = repmat(opx.Mn, nP, 1);
row.Ti_opx = repmat(opx.Ti, nP, 1);
row.Cr_opx = repmat(opx.Cr, nP, 1);

row.Si_cpx = repmat(cpx.Si, nP, 1);
row.Al_cpx = repmat(cpx.Al, nP, 1);
row.Fe2_input_cpx = repmat(cpx.Fe2, nP, 1);
row.Fe3_input_cpx = repmat(cpx.Fe3, nP, 1);
row.FeTot_cpx = repmat(cpx.FeTot, nP, 1);
row.Mg_cpx = repmat(cpx.Mg, nP, 1);
row.Ca_cpx = repmat(cpx.Ca, nP, 1);
row.Na_cpx = repmat(cpx.Na, nP, 1);
row.Mn_cpx = repmat(cpx.Mn, nP, 1);
row.Ti_cpx = repmat(cpx.Ti, nP, 1);
row.Cr_cpx = repmat(cpx.Cr, nP, 1);

row.AlIV_opx = repmat(opxTerms.AlIV, nP, 1);
row.AlVI_opx = repmat(opxTerms.AlVI, nP, 1);
row.XFmAl2SiO6_opx = repmat(opxTerms.XFmAl2SiO6, nP, 1);

row.AlIV_cpx = repmat(cpxTerms.AlIV, nP, 1);
row.AlVI_cpx = repmat(cpxTerms.AlVI, nP, 1);
row.Fe3_calc_cpx = repmat(cpxTerms.Fe3, nP, 1);
row.Fe2_calc_cpx = repmat(cpxTerms.Fe2, nP, 1);
row.XJd_cpx = repmat(cpxTerms.XJd, nP, 1);
row.XCaTs_cpx = repmat(cpxTerms.XCaTs, nP, 1);
row.XCaTi_cpx = repmat(cpxTerms.XCaTi, nP, 1);
row.XCrCaTs_cpx = repmat(cpxTerms.XCrCaTs, nP, 1);
row.XDiHd_cpx = repmat(cpxTerms.XDiHd, nP, 1);

row.invT_Eq7 = invT;
row.equation_domain_valid = equationDomainValid;
row.TEq7_K = TEq7_K;
row.TEq7_C = TEq7_C;

% Ballhaus-style standardized temperature aliases for downstream launchers.
row.T_K = row.TEq7_K;
row.T_deg = row.TEq7_C;

end

function px = preparePyroxeneRow(data_px, mineralLabel)
% preparePyroxeneRow
% Extract one selected pyroxene row while preserving NaN. The legacy Fe
% convention of the original script is retained and documented in the file
% header.

if height(data_px) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

px = struct();
px.Si = getRequiredVar(data_px, 'Si_cation_apfu', mineralLabel);
px.Al = getRequiredVar(data_px, 'Al_cation_apfu', mineralLabel);
px.Fe2 = getRequiredVar(data_px, 'Fe_cation_apfu', mineralLabel);
px.Mg = getRequiredVar(data_px, 'Mg_cation_apfu', mineralLabel);
px.Ca = getRequiredVar(data_px, 'Ca_cation_apfu', mineralLabel);
px.Na = getRequiredVar(data_px, 'Na_cation_apfu', mineralLabel);
px.Mn = getRequiredVar(data_px, 'Mn_cation_apfu', mineralLabel);
px.Ti = getRequiredVar(data_px, 'Ti_cation_apfu', mineralLabel);
px.Cr = getRequiredVar(data_px, 'Cr_cation_apfu', mineralLabel);

[px.Fe3, hasFe3Column] = getOptionalFe3(data_px);

if hasFe3Column
    px.FeTot = px.Fe2 + px.Fe3;
else
    px.FeTot = px.Fe2;
    px.Fe3 = 0;
end

end

function terms = calcOpxTerms(opx)
% calcOpxTerms
% Calculate Opx terms required by Lessel and Putirka (2015) equation (7)
% while preserving NaN values.

terms = struct();
terms.AlIV = clampDerivedNonnegative(2 - opx.Si);
terms.AlVI = clampDerivedNonnegative(opx.Al - terms.AlIV);

if isnan(terms.AlIV) || isnan(terms.AlVI)
    terms.XFmAl2SiO6 = NaN;
else
    terms.XFmAl2SiO6 = min(terms.AlIV, terms.AlVI);
end

end

function terms = calcCpxTerms(cpx)
% calcCpxTerms
% Calculate Cpx terms required by Lessel and Putirka (2015) equation (7).
% Finite negative normative components are clamped to zero, while NaN is
% retained. A negative calculated Fe2+ is treated as invalid and retained as
% NaN rather than being converted to zero.

terms = struct();
terms.AlIV = clampDerivedNonnegative(2 - cpx.Si);
terms.AlVI = clampDerivedNonnegative(cpx.Al - terms.AlIV);

terms.Fe3 = clampDerivedNonnegative( ...
    terms.AlIV + cpx.Na - terms.AlVI - cpx.Cr - 2 .* cpx.Ti);

Fe2_raw = cpx.FeTot - terms.Fe3;
if isnan(Fe2_raw)
    terms.Fe2 = NaN;
elseif isfinite(Fe2_raw) && Fe2_raw < 0
    terms.Fe2 = NaN;
else
    terms.Fe2 = Fe2_raw;
end

if isnan(terms.AlVI) || isnan(cpx.Na)
    terms.XJd = NaN;
else
    terms.XJd = clampDerivedNonnegative(min(terms.AlVI, cpx.Na));
end

terms.XCaTs = clampDerivedNonnegative(terms.AlVI - terms.XJd);
terms.XCaTi = clampDerivedNonnegative((terms.AlIV - terms.XCaTs) ./ 2);
terms.XCrCaTs = clampDerivedNonnegative(cpx.Cr ./ 2);
terms.XDiHd = clampDerivedNonnegative( ...
    cpx.Ca - terms.XCaTi - terms.XCaTs - terms.XCrCaTs);

end

function value = clampDerivedNonnegative(value)
% clampDerivedNonnegative
% Clamp a finite derived component to zero when negative while preserving
% NaN. This function applies only to normative derived components.

if isnan(value)
    return;
end
if isfinite(value) && value < 0
    value = 0;
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
    error('%s variable %s must be a real numeric scalar in a selected 1-row table.', ...
        mineralLabel, varName);
end

end

function [value, columnPresent] = getOptionalFe3(tbl)
% getOptionalFe3
% Return Fe3_cation_apfu when present, preserving NaN. If absent, return
% zero together with columnPresent = false.

columnPresent = ismember('Fe3_cation_apfu', tbl.Properties.VariableNames);
if columnPresent
    value = tbl.Fe3_cation_apfu;
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value)
        error('Variable Fe3_cation_apfu must be a real numeric scalar.');
    end
else
    value = 0;
end

end
