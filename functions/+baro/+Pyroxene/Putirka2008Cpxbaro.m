function results = Putirka2008Cpxbaro(rawdata_struct, T_degreeC)
% functions/+baro/+Pyroxene/Putirka2008Cpxbaro.m
% Compatibility target: MATLAB R2024b
%
% Clinopyroxene-only geobarometer, Equation (32a)
% Putirka, K.D. (2008)
% Reviews in Mineralogy and Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Clinopyroxene analysis and
% calculates pressure using Equation (32a) of Putirka (2008; p. 91).
%
% Equation (32a) uses temperature and Clinopyroxene composition only. No
% liquid composition or second mineral is required.
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Clinopyroxene analysis, one output
% row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka (2008) recalibrated the Nimis (1995) Cpx-only approach using
% experiments performed from approximately 0.001 to 80 kbar (p. 91).
% The reported calibration standard error of estimate is approximately
% +/-3.1 kbar for the anhydrous calibration data.
%
% Equation (32a) retains systematic error for hydrous experiments. Putirka
% (2008) therefore introduced Equation (32b), which includes liquid H2O, to
% reduce this hydrous bias. Equation (32a) should consequently be treated
% cautiously for Cpx crystallized from hydrous melts.
%
% Later testing by Neave and Putirka (2017) also found that Equation (32a)
% may overestimate low-to-moderate pressures by approximately 1-2 kbar,
% even for H2O-poor basaltic systems. Passing the numerical range checks in
% this implementation does not demonstrate that the calculated pressure is
% accurate for a particular natural sample.
%
% Temperature uncertainty propagates directly into pressure because
% Equation (32a) contains both T and ln(T). The supplied temperature should
% be obtained independently using a suitable thermometer and should refer
% to the same crystallization or equilibration event recorded by the Cpx.
%
% Clinopyroxene components must be calculated on a six-oxygen basis using
% the normative scheme summarized in Putirka (2008, Tables 2-3) and Putirka
% et al. (2003). The input cation data used here are therefore expected to
% already be normalized to six oxygens. A Cpx cation sum close to four is an
% important analytical and normalization quality check.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite calculated pressure is outside 0.001-80 kbar,
%   2) a required calculation input contains NaN,
%   3) XAl(VI) is non-positive or a derived component is invalid,
%   4) a calculated pressure is NaN, Inf, or negative, or
%   5) the Cpx cation sum differs noticeably from four.
%
% Calculations outside the calibration range are retained for diagnostic
% purposes but should be treated as extrapolations.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must contain
% Cpx cations normalized on a six-oxygen basis.
%
% Required Cpx variables:
%   Si_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu
%   Mn_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   Cr_cation_apfu
%
% Optional Cpx variables retained in the output when present:
%   K_cation_apfu
%   Fe3_cation_apfu
%
% Fe_cation_apfu is treated as the Fe term used in the original normative
% component calculation. Fe3_cation_apfu, when present, is retained only for
% traceability and is not subtracted from Fe_cation_apfu.
%
% Finite cation inputs used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if XAl(VI) becomes zero, ln[XAl(VI)]
% is undefined and the corresponding pressure remains NaN.
%
% No liquid composition is used by this barometer.
%
% -------------------------------------------------------------------------
% CLINOPYROXENE COMPONENT CALCULATION
%
% Cpx cations are assumed to be normalized to six oxygens. Define:
%
%   XAlIV = max(2 - XSi, 0)
%   XAlVI = max(XAl - XAlIV, 0)
%   XJd   = min(XAlVI, XNa)
%   XCaTs = max(XAlVI - XJd, 0)
%   XCaTi = max((XAlIV - XCaTs)/2, 0), when XAlIV > XCaTs
%   XCrCaTs = max(XCr/2, 0)
%   XFm   = XFe + XMg
%
% For the usual Ca-sufficient case:
%
%   XDiHd = max(XCa - XCaTi - XCaTs - XCrCaTs, 0)
%   XEnFs = max((XFm - XDiHd)/2, 0)
%
% For Ca-poor pyroxenes, the alternative component allocation described by
% Putirka et al. (2003) is applied before XEnFs is calculated.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION
%
% Putirka (2008), Equation (32a), p. 91:
%
%   P(kbar) =
%       3205
%     + 0.384*T(K)
%     - 518*ln[T(K)]
%     - 5.62*XMg_cpx
%     + 83.2*XNa_cpx
%     + 68.2*XDiHd_cpx
%     + 2.52*ln[XAlVI_cpx]
%     - 51.1*(XDiHd_cpx)^2
%     + 34.8*(XEnFs_cpx)^2
%
% Natural logarithms are used. Temperature is in Kelvin and pressure is
% returned in kbar.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008Cpxbaro(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing the Cpx table described above
%   T_degreeC      : non-negative numeric scalar or vector in degreeC.
%                    NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Cpx analysis.
%

%% Input validation
% Accept scalar or vector temperature input so that the fixed-temperature
% and temperature-range launchers use the same implementation.
if nargin < 2
    error('Putirka2008Cpxbaro requires (rawdata_struct, T_degreeC).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(T_degreeC) || isempty(T_degreeC) || ~isvector(T_degreeC) || ...
        any(isinf(T_degreeC(:))) || ...
        any(isfinite(T_degreeC(:)) & T_degreeC(:) < 0)
    error(['T_degreeC must be a non-negative numeric scalar or vector. ' ...
           'NaN is allowed and retained, but Inf and finite negative ' ...
           'values are prohibited.']);
end

T_degreeC = T_degreeC(:);

%% 1) Retrieve Clinopyroxene cation dataset
% Extract the required table from the input struct. The source table is not
% modified; selected rows are read during each calculation.
disp('=== Step 1: Preparing Clinopyroxene cation dataset ===');

if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

if height(rawdata_struct.Cpx) < 1
    error('rawdata_struct.Cpx is empty.');
end

dataset_cpx = rawdata_struct.Cpx;

requiredCpxVariables = {'Si_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Fe_cation_apfu', 'Mn_cation_apfu', ...
    'Mg_cation_apfu', 'Ca_cation_apfu', 'Na_cation_apfu', ...
    'Cr_cation_apfu'};

checkRequiredVariables(dataset_cpx, requiredCpxVariables, 'Cpx');

disp('=== Preparing Clinopyroxene cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-Cpx result in a fixed-size cell buffer and concatenate
% once after the interactive loop. This avoids repeated growth of the result
% table and avoids changing the cell-array size inside the loop.
disp('=== Step 2: Preparing output container ===');

maxResultBlocks = max(1024, height(dataset_cpx));
resultBlocks = cell(maxResultBlocks, 1);
nResultBlocks = 0;

% Published experimental pressure span used as a warning envelope.
calibrationP_min_kbar = 0.001;
calibrationP_max_kbar = 80;

modelCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive Cpx selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a completed calculation.
disp('=== Step 3: Selecting a data code from the list (Cpx) ===');

while true
    if nResultBlocks >= maxResultBlocks
        fprintf(2, ...
            ['WARNING: The fixed result-buffer limit of %d selections was ' ...
             'reached. Completed calculations will be returned without ' ...
             'enlarging the result array.\n'], ...
            maxResultBlocks);
        break;
    end

    % ----- Clinopyroxene selection -----
    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Clinopyroxene (Cpx) data:', ...
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
    disp('=== Step 4: Calculating the pressure ===');

    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    % List NaN only in variables used directly by Equation (32a) or its
    % component calculations. NaN values are retained and do not stop the
    % calculation.
    nanInputNames = findNaNInputs(selectedData_cpx, T_degreeC);

    % Reject only Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(selectedData_cpx);

    row = calcPressure(selectedData_cpx, T_degreeC);

    % Repeat the identifier for all temperatures in the current calculation.
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, 'dataCode_cpx', 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_cpx)) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_cpx)) ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Print the main model limitation once per function call.
    if ~modelCautionIssued
        fprintf(2, ...
            ['CAUTION: Putirka (2008) Equation (32a) preserves systematic ' ...
             'error for hydrous experiments (p. 91). It also has a ' ...
             'temperature-dependent pressure result and should be paired ' ...
             'with an independently appropriate thermometer.\n']);
        modelCautionIssued = true;
    end

    % Warn when finite calculated pressures fall outside the experimental
    % pressure interval used to recalibrate Equation (32a).
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the approximately ' ...
             '0.001-80 kbar experimental range used for Putirka (2008) ' ...
             'Equation (32a; p. 91). %d of %d finite pressure point(s) ' ...
             'are outside the range; calculated finite range = ' ...
             '%.4g-%.4g kbar for %s.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_cpx)));
    end

    % A cation sum near four is an important quality check for a six-oxygen
    % pyroxene formula. The +/-0.05 threshold is used only as an
    % implementation diagnostic and is not a published calibration limit.
    finiteCationSum = isfinite(row.cationSum_cpx);
    cationSumOutsideDiagnostic = finiteCationSum & ...
        abs(row.cationSum_cpx - 4) > 0.05;

    if any(cationSumOutsideDiagnostic)
        fprintf(2, ...
            ['WARNING: Cpx cation sum differs from 4 by more than 0.05 for ' ...
             '%s (cation sum = %.6g). Equation (32a) requires Cpx cations ' ...
             'on a six-oxygen basis. The threshold is a diagnostic, not a ' ...
             'published calibration limit.\n'], ...
            char(string(selectedCode_cpx)), ...
            row.cationSum_cpx(find(finiteCationSum, 1, 'first')));
    end

    % List the exact required input names containing NaN. For vector
    % temperature input, the NaN element indices are included.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Equation (32a) input(s) for ' ...
             '%s: %s.\n' ...
             '         NaN values were retained and propagated; the ' ...
             'calculated pressure may remain NaN.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Report derived component quantities outside their mathematical domain.
    invalidComponentNames = findInvalidComponentQuantities(row);
    if ~isempty(invalidComponentNames)
        fprintf(2, ...
            ['WARNING: Invalid Equation (32a) component or logarithm term ' ...
             'was found for %s: %s.\n' ...
             '         XAl(VI) must be finite and > 0 for ln[XAl(VI)]. ' ...
             'Affected pressure values remain NaN where applicable.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(invalidComponentNames, ', ')));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_cpx)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressure is retained for diagnostic purposes.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s ' ...
             '(%d of %d points). The values were retained for diagnostic ' ...
             'purposes.\n'], ...
            char(string(selectedCode_cpx)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Putirka2008Cpxbaro', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered result blocks once. If the user canceled before
% any calculation, return an empty table.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'primaryPressureEquation', 'Putirka2008 Eq. (32a)', ...
    'pressureUnit', 'kbar', ...
    'temperatureUnitInEquation', 'K', ...
    'calibrationSEE_kbar', 3.1, ...
    'calibrationPressureRange_kbar', [0.001 80]);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_cpx, T_degreeC)
% findNaNInputs
% Return names of Equation (32a) calculation inputs containing NaN. NaN
% values do not cause an error and are never replaced by zero.

cpxVariables = {'Si_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Fe_cation_apfu', 'Mn_cation_apfu', ...
    'Mg_cation_apfu', 'Ca_cation_apfu', 'Na_cation_apfu', ...
    'Cr_cation_apfu'};

maxNames = 1 + numel(cpxVariables);
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

for i = 1:numel(cpxVariables)
    variableName = cpxVariables{i};
    variableValue = data_cpx.(variableName);
    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Cpx." + string(variableName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_cpx)
% validateNonNegativeInputs
% Reject Inf and finite negative values in variables directly used by
% Equation (32a). Zero and NaN are intentionally allowed and retained.

cpxVariables = {'Si_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Fe_cation_apfu', 'Mn_cation_apfu', ...
    'Mg_cation_apfu', 'Ca_cation_apfu', 'Na_cation_apfu', ...
    'Cr_cation_apfu'};

invalidInputBuffer = strings(numel(cpxVariables), 1);
nInvalidInputs = 0;

for i = 1:numel(cpxVariables)
    variableName = cpxVariables{i};
    variableValue = data_cpx.(variableName);
    if any(isinf(variableValue(:)) | ...
            (isfinite(variableValue(:)) & variableValue(:) < 0))
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Cpx." + string(variableName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['Putirka2008Cpxbaro: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are ' ...
           'prohibited. Invalid input(s): ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_cpx, T_degreeC)
% calcPressure
% Compute Putirka (2008) Equation (32a) pressure for one Cpx row at one or
% more input temperatures. NaN values propagate through the calculation.
%
% Inputs:
%   data_cpx   : 1-row Cpx table containing six-oxygen cations
%   T_degreeC  : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Extract one-row Cpx cation data and calculate normative components.
cpx = prepareCpxRow(data_cpx);

% Equation (32a) requires ln[XAl(VI)]. Avoid complex values for a
% non-positive argument. A positive Inf is retained as Inf.
lnXAlVI_scalar = NaN;
if isfinite(cpx.XAlVI) && cpx.XAlVI > 0
    lnXAlVI_scalar = log(cpx.XAlVI);
elseif isinf(cpx.XAlVI) && cpx.XAlVI > 0
    lnXAlVI_scalar = Inf;
end

% Expand composition-dependent scalars to the temperature-vector length.
XSi_cpx = repmat(cpx.XSi, nT, 1);
XTi_cpx = repmat(cpx.XTi, nT, 1);
XAl_cpx = repmat(cpx.XAl, nT, 1);
XFe_cpx = repmat(cpx.XFe, nT, 1);
XMn_cpx = repmat(cpx.XMn, nT, 1);
XMg_cpx = repmat(cpx.XMg, nT, 1);
XCa_cpx = repmat(cpx.XCa, nT, 1);
XNa_cpx = repmat(cpx.XNa, nT, 1);
XK_cpx = repmat(cpx.XK, nT, 1);
XCr_cpx = repmat(cpx.XCr, nT, 1);
XFe3_cpx = repmat(cpx.XFe3, nT, 1);
cationSum_cpx = repmat(cpx.cationSum, nT, 1);

XAlIV_cpx = repmat(cpx.XAlIV, nT, 1);
XAlVI_cpx = repmat(cpx.XAlVI, nT, 1);
XJd_cpx = repmat(cpx.XJd, nT, 1);
XCaTs_cpx = repmat(cpx.XCaTs, nT, 1);
XCaTi_cpx = repmat(cpx.XCaTi, nT, 1);
XCrCaTs_cpx = repmat(cpx.XCrCaTs, nT, 1);
XDiHd_cpx = repmat(cpx.XDiHd, nT, 1);
XEnFs_cpx = repmat(cpx.XEnFs, nT, 1);
XFmCaTs_cpx = repmat(cpx.XFmCaTs, nT, 1);
XFmTi_cpx = repmat(cpx.XFmTi, nT, 1);
lnXAlVI_cpx = repmat(lnXAlVI_scalar, nT, 1);

% Putirka (2008), Equation (32a), p. 91.
% No finite-value guard is applied to the complete expression: NaN and Inf
% propagate and remain available for diagnosis in the output table.
P_kbar = ...
    3205 ...
    + 0.384 .* T_K ...
    - 518 .* log(T_K) ...
    - 5.62 .* XMg_cpx ...
    + 83.2 .* XNa_cpx ...
    + 68.2 .* XDiHd_cpx ...
    + 2.52 .* lnXAlVI_cpx ...
    - 51.1 .* (XDiHd_cpx .^ 2) ...
    + 34.8 .* (XEnFs_cpx .^ 2);

P_GPa = P_kbar ./ 10;

% Applicability and diagnostic flags. These do not prove equilibrium or
% applicability to a particular natural sample.
isWithinCalibrationPRange = ...
    isfinite(P_kbar) & P_kbar >= 0.001 & P_kbar <= 80;

isWithinEquationDomain = ...
    isfinite(T_K) & T_K > 0 & ...
    isfinite(XAlVI_cpx) & XAlVI_cpx > 0 & ...
    isfinite(XDiHd_cpx) & XDiHd_cpx >= 0 & ...
    isfinite(XEnFs_cpx) & XEnFs_cpx >= 0;

isCationSumNearFour = ...
    isfinite(cationSum_cpx) & abs(cationSum_cpx - 4) <= 0.05;

% Pack outputs using equal-length, pre-sized vectors.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.P_SEE_calibration_kbar = repmat(3.1, nT, 1);
row.P_SEE_calibration_GPa = repmat(0.31, nT, 1);

row.XSi_cpx = XSi_cpx;
row.XTi_cpx = XTi_cpx;
row.XAl_cpx = XAl_cpx;
row.XFe_cpx = XFe_cpx;
row.XMn_cpx = XMn_cpx;
row.XMg_cpx = XMg_cpx;
row.XCa_cpx = XCa_cpx;
row.XNa_cpx = XNa_cpx;
row.XK_cpx = XK_cpx;
row.XCr_cpx = XCr_cpx;
row.XFe3_cpx = XFe3_cpx;
row.cationSum_cpx = cationSum_cpx;

row.XAlIV_cpx = XAlIV_cpx;
row.XAlVI_cpx = XAlVI_cpx;
row.XJd_cpx = XJd_cpx;
row.XCaTs_cpx = XCaTs_cpx;
row.XCaTi_cpx = XCaTi_cpx;
row.XCrCaTs_cpx = XCrCaTs_cpx;
row.XDiHd_cpx = XDiHd_cpx;
row.XEnFs_cpx = XEnFs_cpx;
row.XFmCaTs_cpx = XFmCaTs_cpx;
row.XFmTi_cpx = XFmTi_cpx;
row.lnXAlVI_cpx = lnXAlVI_cpx;

row.isWithinCalibrationPRange = isWithinCalibrationPRange;
row.isWithinEquationDomain = isWithinEquationDomain;
row.isCationSumNearFour = isCationSumNearFour;

end

function cpx = prepareCpxRow(data_cpx)
% prepareCpxRow
% Extract one-row Cpx cation data and calculate the normative components
% used by Equation (32a). Existing NaN values are retained.

if height(data_cpx) ~= 1
    error('Cpx input must be a 1-row table.');
end

cpx = struct();
cpx.XSi = getVarOrError(data_cpx, 'Si_cation_apfu', 'Cpx');
cpx.XTi = getVarOrError(data_cpx, 'Ti_cation_apfu', 'Cpx');
cpx.XAl = getVarOrError(data_cpx, 'Al_cation_apfu', 'Cpx');
cpx.XFe = getVarOrError(data_cpx, 'Fe_cation_apfu', 'Cpx');
cpx.XMn = getVarOrError(data_cpx, 'Mn_cation_apfu', 'Cpx');
cpx.XMg = getVarOrError(data_cpx, 'Mg_cation_apfu', 'Cpx');
cpx.XCa = getVarOrError(data_cpx, 'Ca_cation_apfu', 'Cpx');
cpx.XNa = getVarOrError(data_cpx, 'Na_cation_apfu', 'Cpx');
cpx.XCr = getVarOrError(data_cpx, 'Cr_cation_apfu', 'Cpx');
cpx.XK = getVarOrDefault(data_cpx, 'K_cation_apfu', 0);
cpx.XFe3 = getVarOrDefault(data_cpx, 'Fe3_cation_apfu', NaN);

cpx.cationSum = cpx.XSi + cpx.XTi + cpx.XAl + cpx.XFe + ...
    cpx.XMn + cpx.XMg + cpx.XCa + cpx.XNa + cpx.XK + cpx.XCr;

% Tetrahedral and octahedral Al.
cpx.XAlIV = maxPreserveNaN(2 - cpx.XSi, 0);
cpx.XAlVI = maxPreserveNaN(cpx.XAl - cpx.XAlIV, 0);

% Step 1: Jadeite.
cpx.XJd = minPreserveNaN(cpx.XAlVI, cpx.XNa);
cpx.XJd = maxPreserveNaN(cpx.XJd, 0);

% Step 2: Ca-Tschermak component.
cpx.XCaTs = maxPreserveNaN(cpx.XAlVI - cpx.XJd, 0);

% Step 3: CaTi component.
if isnan(cpx.XAlIV) || isnan(cpx.XCaTs)
    cpx.XCaTi = NaN;
elseif cpx.XAlIV > cpx.XCaTs
    cpx.XCaTi = max((cpx.XAlIV - cpx.XCaTs) ./ 2, 0);
else
    cpx.XCaTi = 0;
end

% Step 4: Cr-bearing Ca-Tschermak component.
cpx.XCrCaTs = maxPreserveNaN(cpx.XCr ./ 2, 0);

% Steps 5-6: DiHd and EnFs, including the Ca-poor alternative.
XFm = cpx.XFe + cpx.XMg;

if any(isnan([cpx.XCa, cpx.XCaTs, cpx.XCaTi, ...
        cpx.XCrCaTs, XFm]))
    cpx.XDiHd = NaN;
    cpx.XEnFs = NaN;
    cpx.XFmCaTs = NaN;
    cpx.XFmTi = NaN;
elseif cpx.XCa >= (cpx.XCaTs + cpx.XCaTi)
    cpx.XDiHd = max( ...
        cpx.XCa - cpx.XCaTi - cpx.XCaTs - cpx.XCrCaTs, 0);
    cpx.XFmCaTs = 0;
    cpx.XFmTi = 0;
    cpx.XEnFs = max((XFm - cpx.XDiHd) ./ 2, 0);
else
    cpx.XDiHd = 0;
    cpx.XCaTs = cpx.XCa;

    VIAlex = max(cpx.XAlVI - cpx.XCaTs, 0);
    cpx.XFmCaTs = max(VIAlex - cpx.XCaTs, 0);
    cpx.XFmTi = max( ...
        (cpx.XAlIV - cpx.XCaTs - cpx.XFmCaTs) ./ 2, 0);
    cpx.XEnFs = max( ...
        (XFm - cpx.XFmCaTs - cpx.XFmTi) ./ 2, 0);
end

end

function invalidComponentNames = findInvalidComponentQuantities(row)
% findInvalidComponentQuantities
% Return names of composition-dependent quantities outside the equation or
% component domain. These values are identical in every temperature row, so
% the first row is sufficient for the component checks.

maxNames = 9;
invalidBuffer = strings(maxNames, 1);
nInvalid = 0;

if isempty(row)
    invalidComponentNames = strings(0, 1);
    return;
end

if ~isfinite(row.XAlIV_cpx(1)) || row.XAlIV_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XAlIV_cpx";
end
if ~isfinite(row.XAlVI_cpx(1)) || row.XAlVI_cpx(1) <= 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XAlVI_cpx";
end
if ~isfinite(row.lnXAlVI_cpx(1))
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "lnXAlVI_cpx";
end
if ~isfinite(row.XJd_cpx(1)) || row.XJd_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XJd_cpx";
end
if ~isfinite(row.XCaTs_cpx(1)) || row.XCaTs_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XCaTs_cpx";
end
if ~isfinite(row.XCaTi_cpx(1)) || row.XCaTi_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XCaTi_cpx";
end
if ~isfinite(row.XCrCaTs_cpx(1)) || row.XCrCaTs_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XCrCaTs_cpx";
end
if ~isfinite(row.XDiHd_cpx(1)) || row.XDiHd_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XDiHd_cpx";
end
if ~isfinite(row.XEnFs_cpx(1)) || row.XEnFs_cpx(1) < 0
    nInvalid = nInvalid + 1;
    invalidBuffer(nInvalid) = "XEnFs_cpx";
end

invalidComponentNames = invalidBuffer(1:nInvalid);

end

function checkRequiredVariables(tbl, requiredVariables, tableName)
% checkRequiredVariables
% Stop before the interactive loop when required calculation columns are
% absent. Existing NaN values within present columns are allowed.

missingVariables = ...
    requiredVariables(~ismember(requiredVariables, tbl.Properties.VariableNames));

if ~isempty(missingVariables)
    error('%s table must contain variable(s): %s', ...
        tableName, strjoin(missingVariables, ', '));
end

end

function value = getVarOrError(tbl, variableName, mineralLabel)
% getVarOrError
% Retrieve a required scalar table variable without altering NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);
if ~isscalar(value) || ~isnumeric(value) || ~isreal(value)
    error('Variable %s must be one real numeric scalar.', variableName);
end
value = double(value);

end

function value = getVarOrDefault(tbl, variableName, defaultValue)
% getVarOrDefault
% Retrieve an optional scalar variable. Only an absent column receives the
% stated default. An existing NaN value is retained.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);
    if ~isscalar(value) || ~isnumeric(value) || ~isreal(value)
        error('Variable %s must be one real numeric scalar.', variableName);
    end
    value = double(value);
else
    value = defaultValue;
end

end

function value = maxPreserveNaN(a, b)
% maxPreserveNaN
% Scalar maximum that explicitly preserves NaN.

if isnan(a) || isnan(b)
    value = NaN;
else
    value = max(a, b);
end

end

function value = minPreserveNaN(a, b)
% minPreserveNaN
% Scalar minimum that explicitly preserves NaN.

if isnan(a) || isnan(b)
    value = NaN;
else
    value = min(a, b);
end

end
