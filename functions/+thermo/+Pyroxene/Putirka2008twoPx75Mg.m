function results = Putirka2008twoPx75Mg(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/Putirka2008twoPx75Mg.m
% Tested with MATLAB R2024b
%
% Two-pyroxene thermometer for Mg-rich systems
% Putirka, K.D. (2008), Equation (37)
% Reviews in Mineralogy and Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis and one
% Clinopyroxene analysis and calculates temperature using the Mg-rich
% two-pyroxene thermometer of Putirka (2008), Equation (37).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Opx-Cpx pair, one output row is returned for every pressure
% value supplied in P_kbar. It is therefore compatible with both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka (2008) discusses two-pyroxene equilibrium and calibration on
% pp. 95-96. Equation (37) and its component definitions are given on p. 96.
% Equation (37) was calibrated using only experimental Opx-Cpx pairs having:
%
%   Mg#_cpx > 0.75
%
% Figure 10c on p. 95 reports the following statistics:
%   Calibration subset, Mg#_cpx > 0.75:
%     SEE = +/-38 degreeC, R^2 = 0.93, n = 390
%
% When Equation (37) is applied to all 487 experimental pairs, including
% lower-Mg# compositions, Figure 10c reports:
%     SEE = +/-60 degreeC, R^2 = 0.86, n = 487
%
% Equation (37) should therefore be used only with Mg-rich Clinopyroxene
% compositions satisfying Mg#_cpx > 0.75. A result is still retained when
% this condition is not met, but a non-stopping warning is printed.
%
% Putirka (2008) does not provide one explicit rectangular P-T calibration
% box for Equation (37). Figure 10 displays the two-pyroxene experimental
% comparison over approximately:
%
%   Temperature : 600-1800 degreeC
%   Pressure    : 0-50 kbar for the associated two-pyroxene dataset
%
% These graphical intervals are used here only as non-stopping screening
% ranges. They are not represented as strict universal calibration limits.
%
% Two-pyroxene equilibrium must be evaluated before interpreting the
% temperature. On p. 95, Putirka (2008) reports:
%
%   KD(Fe-Mg)cpx-opx = 1.09 +/- 0.14
%
% from 311 experiments and states that data at the outskirts of 3 standard
% deviations were excluded for calibration and testing. This implementation
% therefore reports, but does not stop for, pairs outside the approximate
% 3-sigma interval 0.67-1.51.
%
% Pyroxene components must be calculated on a 6-oxygen basis using the same
% allocation procedure as the calibration. The component calculation scheme
% is given in Tables 2 and 3 on pp. 88-89. Putirka (2008) notes that a good
% pyroxene analysis should have a total cation sum close to 4.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside the approximate 0-50 kbar Figure-10
%      experimental envelope,
%   2) a finite calculated temperature is outside the approximate
%      600-1800 degreeC Figure-10 comparison envelope,
%   3) Mg#_cpx is not greater than 0.75,
%   4) KD(Fe-Mg)cpx-opx is outside the approximate 3-sigma interval
%      0.67-1.51,
%   5) an explicitly stored calculation input is NaN,
%   6) a logarithm, component denominator, or Equation (37) denominator is
%      invalid, or
%   7) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%   rawdata_struct.Cpx : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% The pyroxene tables must contain normalized cations, preferably on a
% 6-oxygen basis.
%
% Required variables for both Opx and Cpx:
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
%   K_cation_apfu
%
% Na, Mn, Ti, and Cr are required because they enter Equation (37) directly
% or affect the Table-2/Table-3 component allocations used by the equation.
% K does not enter Equation (37) or those component terms; an absent K column
% is therefore assigned zero only for cation-sum and output traceability.
%
% Fe_cation_apfu is treated as total Fe on the 6-oxygen cation basis.
%
% An explicitly stored NaN in a calculation input is retained and propagated;
% it is never replaced by zero. All finite composition values must be greater
% than or equal to zero. Negative finite values and Inf stop the calculation.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and liquid NaN warnings is not applicable.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Putirka (2008), Equation (37), p. 96:
%
%   10^4 / T(degreeC) =
%       13.4
%     - 3.4  * ln(XEnFs_cpx / XEnFs_opx)
%     + 5.59 * ln(XMg_cpx)
%     - 8.8  * MgNum_cpx
%     + 23.85 * XMn_opx
%     + 6.48 * XFmAl2SiO6_opx
%     - 2.38 * XDi_cpx
%     - 0.044 * P_kbar
%
% where:
%   MgNum_cpx = XMg_cpx / (XMg_cpx + XFe_cpx)
%   XDi_cpx   = [XMg_cpx / (XMg_cpx + XFe_cpx + XMn_cpx)] * XDiHd_cpx
%
% Therefore:
%   T(degreeC) = 10000 / denominator
%   T(K)       = T(degreeC) + 273.15
%
% The EnFs, FmAl2SiO6, DiHd, and Di components are calculated following
% Putirka (2008), Tables 2 and 3, on a 6-oxygen basis.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008twoPx75Mg(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Opx and Cpx tables
%   P_kbar         : pressure in kbar; finite, non-negative numeric scalar
%                    or vector
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Opx-Cpx pair. T_degreeC and T_deg are standardized aliases for
%             downstream launcher and plotting compatibility.
%

%% Input validation
if nargin < 2
    error('Putirka2008twoPx75Mg requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve pyroxene datasets
disp('=== Step 1: Preparing pyroxene datasets ===');

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end
if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end
if isempty(rawdata_struct.Opx)
    error('rawdata_struct.Opx is empty.');
end
if isempty(rawdata_struct.Cpx)
    error('rawdata_struct.Cpx is empty.');
end

dataset_opx = rawdata_struct.Opx;
dataset_cpx = rawdata_struct.Cpx;

validateRequiredColumns(dataset_opx, 'Opx');
validateRequiredColumns(dataset_cpx, 'Cpx');

disp('=== Preparing pyroxene datasets has been finished ===');

%% 2) Initialize output container and applicability limits
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Approximate graphical experimental envelopes represented in Figure 10.
screeningP_min_kbar = 0;
screeningP_max_kbar = 50;
screeningT_min_degC = 600;
screeningT_max_degC = 1800;

pressureOutsideScreening = ...
    P_kbar < screeningP_min_kbar | P_kbar > screeningP_max_kbar;
pressureWarningIssued = false;

% Equation (37) calibration condition and experimental equilibrium criterion.
calibrationMgNum_cpx_min = 0.75;
KD_FeMg_mean = 1.09;
KD_FeMg_sd = 0.14;
KD_FeMg_3sigma_min = KD_FeMg_mean - 3 .* KD_FeMg_sd;
KD_FeMg_3sigma_max = KD_FeMg_mean + 3 .* KD_FeMg_sd;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Opx) ===');

dataCodes_opx = dataset_opx{:, 1};
displayCodes_opx = cellstr(string(dataCodes_opx));

dataCodes_cpx = dataset_cpx{:, 1};
displayCodes_cpx = cellstr(string(dataCodes_cpx));

while true
    % ----- Orthopyroxene selection -----
    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Opx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_opx, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = string(dataCodes_opx(selectedIdx_opx));
    selectedData_opx = dataset_opx(selectedIdx_opx, :);
    disp(['Opx selected: ' char(selectedCode_opx)]);

    % ----- Clinopyroxene selection -----
    disp('=== Step 4: Selecting a data code from the list (Cpx) ===');

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Cpx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_cpx, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = string(dataCodes_cpx(selectedIdx_cpx));
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);
    disp(['Cpx selected: ' char(selectedCode_cpx)]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    % Report explicitly stored NaN values only for variables that enter the
    % component and temperature calculation. K is not used by Equation (37).
    nanInputNames = findNaNInputs(selectedData_opx, selectedData_cpx);

    % Reject negative finite values and Inf. NaN and zero are retained.
    validateNonNegativeInputs(selectedData_opx, selectedData_cpx);

    row = calcTemp(selectedData_opx, selectedData_cpx, P_kbar);

    nRows = height(row);
    row.dataCode_opx = repmat(selectedCode_opx, nRows, 1);
    row.dataCode_cpx = repmat(selectedCode_cpx, nRows, 1);
    row = movevars(row, {'dataCode_opx', 'dataCode_cpx'}, 'Before', 1);

    % Store one completed result block. The complete output table is not
    % enlarged on every interactive iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    printTemperatureSummary( ...
        selectedCode_opx, selectedCode_cpx, row.T_deg);

    % Pressure is common to all selected pairs in this function call.
    if any(pressureOutsideScreening) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate %.4g-%.4g ' ...
             'kbar experimental envelope represented by the associated ' ...
             'two-pyroxene dataset in Putirka (2008), Figure 10. Putirka ' ...
             '(2008) does not state a strict rectangular pressure ' ...
             'calibration range for Equation (37). %d of %d pressure ' ...
             'point(s) are outside the graphical envelope; input range = ' ...
             '%.4g-%.4g kbar.\n'], ...
            screeningP_min_kbar, ...
            screeningP_max_kbar, ...
            sum(pressureOutsideScreening), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn for finite temperatures outside the Figure-10 display envelope.
    printTemperatureRangeWarning( ...
        row.T_deg, screeningT_min_degC, screeningT_max_degC, ...
        selectedCode_opx, selectedCode_cpx);

    % Equation (37) was calibrated only for Cpx Mg# > 0.75.
    if isfinite(row.MgNum_cpx(1)) && ...
            row.MgNum_cpx(1) <= calibrationMgNum_cpx_min
        fprintf(2, ...
            ['WARNING: Mg#_cpx = %.4g for %s & %s does not satisfy the ' ...
             'Equation (37) calibration condition Mg#_cpx > 0.75 ' ...
             '(Putirka, 2008, pp. 95-96). The calculated result is retained ' ...
             'but lies outside the intended compositional calibration ' ...
             'subset.\n'], ...
            row.MgNum_cpx(1), ...
            char(selectedCode_opx), ...
            char(selectedCode_cpx));
    end

    % Report pairs outside the approximate 3-sigma Fe-Mg equilibrium range.
    if isfinite(row.KD_FeMg_cpx_opx(1)) && ...
            (row.KD_FeMg_cpx_opx(1) < KD_FeMg_3sigma_min || ...
             row.KD_FeMg_cpx_opx(1) > KD_FeMg_3sigma_max)
        fprintf(2, ...
            ['WARNING: KD(Fe-Mg)cpx-opx = %.4g for %s & %s is outside the ' ...
             'approximate 3-sigma interval %.4g-%.4g derived from ' ...
             '1.09 +/- 0.14 (Putirka, 2008, p. 95). The selected pair may ' ...
             'not represent equilibrium. The result is retained.\n'], ...
            row.KD_FeMg_cpx_opx(1), ...
            char(selectedCode_opx), ...
            char(selectedCode_cpx), ...
            KD_FeMg_3sigma_min, ...
            KD_FeMg_3sigma_max);
    end

    % Report exact names of explicitly stored NaN calculation inputs.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & ' ...
             '%s: %s.\n' ...
             '         NaN values were retained and propagated; the ' ...
             'calculated temperature may be NaN.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_cpx), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Give a specific warning for invalid component and equation terms.
    invalidEquationTerms = findInvalidEquationTerms(row);
    if ~isempty(invalidEquationTerms)
        fprintf(2, ...
            ['WARNING: Invalid Equation (37) term(s) were found for %s & ' ...
             '%s: %s.\n' ...
             '         Required logarithm arguments and denominators must ' ...
             'be finite and > 0. Affected temperatures were retained as ' ...
             'NaN.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_cpx), ...
            char(strjoin(invalidEquationTerms, ', ')));
    end

    printNonFiniteTemperatureWarning( ...
        row.T_deg, selectedCode_opx, selectedCode_cpx);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Putirka2008twoPx75Mg', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after the interactive loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'primaryTemperatureEquation', 'Putirka2008 Eq. (37)', ...
    'temperatureUnitInEquation', 'degreeC', ...
    'calibrationMgNumCpx', '> 0.75', ...
    'KD_FeMg_cpx_opx_mean', KD_FeMg_mean, ...
    'KD_FeMg_cpx_opx_sd', KD_FeMg_sd);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_opx, data_cpx)
% findNaNInputs
% Return names of explicitly stored variables used by Equation (37) or its
% component calculations that contain NaN. K is excluded because it does not
% enter the equation or the component terms used by it.

activeVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

maxNames = 2 .* numel(activeVariables);
nanNamesBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(activeVariables)
    variableName = activeVariables{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        value = data_opx.(variableName);
        if isnumeric(value) && any(isnan(value(:)))
            nNames = nNames + 1;
            nanNamesBuffer(nNames) = "Opx." + string(variableName);
        end
    end

    if ismember(variableName, data_cpx.Properties.VariableNames)
        value = data_cpx.(variableName);
        if isnumeric(value) && any(isnan(value(:)))
            nNames = nNames + 1;
            nanNamesBuffer(nNames) = "Cpx." + string(variableName);
        end
    end
end

nanInputNames = nanNamesBuffer(1:nNames);

end

function validateNonNegativeInputs(data_opx, data_cpx)
% validateNonNegativeInputs
% Stop when a stored pyroxene-composition input is negative or infinite.
% Zero is allowed. NaN is deliberately allowed and propagated.

variablesToCheck = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

maxNames = 2 .* numel(variablesToCheck);
invalidNamesBuffer = strings(maxNames, 1);
nInvalid = 0;

for i = 1:numel(variablesToCheck)
    variableName = variablesToCheck{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        value = data_opx.(variableName);
        validateScalarVariable(value, 'Opx', variableName);
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            invalidNamesBuffer(nInvalid) = "Opx." + string(variableName);
        end
    end

    if ismember(variableName, data_cpx.Properties.VariableNames)
        value = data_cpx.(variableName);
        validateScalarVariable(value, 'Cpx', variableName);
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            invalidNamesBuffer(nInvalid) = "Cpx." + string(variableName);
        end
    end
end

if nInvalid > 0
    invalidInputNames = invalidNamesBuffer(1:nInvalid);
    error(['Putirka2008twoPx75Mg: pyroxene-composition inputs must not ' ...
           'be negative or infinite. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_cpx, P_kbar)
% calcTemp
% Calculate Putirka (2008) Equation (37) for one selected Opx-Cpx pair and
% a scalar or vector of pressures. Output rows correspond one-to-one with
% input pressure values. NaN inputs are retained and propagated.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

opx = preparePyroxeneRow(data_opx, 'Opx');
cpx = preparePyroxeneRow(data_cpx, 'Cpx');

comp_opx = calcOpxComponents(opx);
comp_cpx = calcCpxComponents(cpx);

% Required EnFs ratio and logarithm.
if isfinite(comp_cpx.XEnFs) && comp_cpx.XEnFs > 0 && ...
        isfinite(comp_opx.XEnFs) && comp_opx.XEnFs > 0
    XEnFs_ratio = comp_cpx.XEnFs ./ comp_opx.XEnFs;
else
    XEnFs_ratio = NaN;
end
ln_XEnFs_ratio = safeLogPositive(XEnFs_ratio);

% Equation (37) contains ln(XMg_cpx).
ln_XMg_cpx = safeLogPositive(cpx.Mg);

% Clinopyroxene Mg# used both by Equation (37) and its calibration criterion.
MgNum_cpx_denominator = cpx.Mg + cpx.Fe;
if isfinite(cpx.Mg) && isfinite(cpx.Fe) && ...
        isfinite(MgNum_cpx_denominator) && MgNum_cpx_denominator > 0
    MgNum_cpx = cpx.Mg ./ MgNum_cpx_denominator;
else
    MgNum_cpx = NaN;
end

% Experimental Fe-Mg equilibrium coefficient from Putirka (2008), p. 95:
% KD = (Fe/Mg)cpx / (Fe/Mg)opx.
if isfinite(cpx.Fe) && isfinite(cpx.Mg) && ...
        isfinite(opx.Fe) && isfinite(opx.Mg) && ...
        cpx.Mg > 0 && opx.Mg > 0
    KD_FeMg_cpx_opx = (cpx.Fe ./ cpx.Mg) ./ (opx.Fe ./ opx.Mg);
else
    KD_FeMg_cpx_opx = NaN;
end

% Putirka (2008), Equation (37), p. 96.
denom_T = ...
    13.4 ...
    - 3.4 .* ln_XEnFs_ratio ...
    + 5.59 .* ln_XMg_cpx ...
    - 8.8 .* MgNum_cpx ...
    + 23.85 .* opx.Mn ...
    + 6.48 .* comp_opx.XFmAl2SiO6 ...
    - 2.38 .* comp_cpx.XDi ...
    - 0.044 .* P_kbar;

% Equation (37) gives T directly in degreeC. Require a finite positive
% denominator; invalid points are retained as NaN.
T_deg = nan(nP, 1);
validDenominator = isfinite(denom_T) & denom_T > 0;
T_deg(validDenominator) = 10000 ./ denom_T(validDenominator);
T_K = T_deg + 273.15;

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PrimaryEquation = repmat("Putirka2008_Eq37", nP, 1);

% Orthopyroxene cations.
row.Si_opx = repmat(opx.Si, nP, 1);
row.Al_opx = repmat(opx.Al, nP, 1);
row.Fe_opx = repmat(opx.Fe, nP, 1);
row.Mg_opx = repmat(opx.Mg, nP, 1);
row.Ca_opx = repmat(opx.Ca, nP, 1);
row.Na_opx = repmat(opx.Na, nP, 1);
row.K_opx = repmat(opx.K, nP, 1);
row.Mn_opx = repmat(opx.Mn, nP, 1);
row.Ti_opx = repmat(opx.Ti, nP, 1);
row.Cr_opx = repmat(opx.Cr, nP, 1);
row.cationSum_opx = repmat(opx.cationSum, nP, 1);

% Clinopyroxene cations.
row.Si_cpx = repmat(cpx.Si, nP, 1);
row.Al_cpx = repmat(cpx.Al, nP, 1);
row.Fe_cpx = repmat(cpx.Fe, nP, 1);
row.Mg_cpx = repmat(cpx.Mg, nP, 1);
row.Ca_cpx = repmat(cpx.Ca, nP, 1);
row.Na_cpx = repmat(cpx.Na, nP, 1);
row.K_cpx = repmat(cpx.K, nP, 1);
row.Mn_cpx = repmat(cpx.Mn, nP, 1);
row.Ti_cpx = repmat(cpx.Ti, nP, 1);
row.Cr_cpx = repmat(cpx.Cr, nP, 1);
row.cationSum_cpx = repmat(cpx.cationSum, nP, 1);

% Orthopyroxene components.
row.AlIV_opx = repmat(comp_opx.AlIV, nP, 1);
row.AlVI_opx = repmat(comp_opx.AlVI, nP, 1);
row.XNaAlSi2O6_opx = repmat(comp_opx.XNaAlSi2O6, nP, 1);
row.XFmTiAlSiO6_opx = repmat(comp_opx.XFmTiAlSiO6, nP, 1);
row.XCrAl2SiO6_opx = repmat(comp_opx.XCrAl2SiO6, nP, 1);
row.XFmAl2SiO6_opx = repmat(comp_opx.XFmAl2SiO6, nP, 1);
row.XCaFmSi2O6_opx = repmat(comp_opx.XCaFmSi2O6, nP, 1);
row.XEnFs_opx = repmat(comp_opx.XEnFs, nP, 1);

% Clinopyroxene components.
row.AlIV_cpx = repmat(comp_cpx.AlIV, nP, 1);
row.AlVI_cpx = repmat(comp_cpx.AlVI, nP, 1);
row.XFe3_cpx = repmat(comp_cpx.XFe3, nP, 1);
row.XJd_cpx = repmat(comp_cpx.XJd, nP, 1);
row.XCaTs_cpx = repmat(comp_cpx.XCaTs, nP, 1);
row.XCaTi_cpx = repmat(comp_cpx.XCaTi, nP, 1);
row.XCrCaTs_cpx = repmat(comp_cpx.XCrCaTs, nP, 1);
row.XDiHd_cpx = repmat(comp_cpx.XDiHd, nP, 1);
row.XEnFs_cpx = repmat(comp_cpx.XEnFs, nP, 1);
row.XDi_cpx = repmat(comp_cpx.XDi, nP, 1);

% Equilibrium, compositional-screening, equation, and temperature outputs.
row.XEnFs_ratio_cpx_opx = repmat(XEnFs_ratio, nP, 1);
row.ln_XEnFs_ratio = repmat(ln_XEnFs_ratio, nP, 1);
row.XMg_cpx = repmat(cpx.Mg, nP, 1);
row.ln_XMg_cpx = repmat(ln_XMg_cpx, nP, 1);
row.MgNum_cpx = repmat(MgNum_cpx, nP, 1);
row.KD_FeMg_cpx_opx = repmat(KD_FeMg_cpx_opx, nP, 1);
row.denom_T = denom_T;

row.T_K = T_K;
row.T_degreeC = T_deg;
row.T_deg = T_deg;

end

function px = preparePyroxeneRow(data_px, mineralLabel)
% preparePyroxeneRow
% Extract one pyroxene composition from normalized cation columns. Existing
% NaN values remain NaN. K is the only optional variable and is assigned
% zero only when its column is absent.

if height(data_px) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

px = struct();

px.Si = getVarRequired(data_px, 'Si_cation_apfu', mineralLabel);
px.Al = getVarRequired(data_px, 'Al_cation_apfu', mineralLabel);
px.Fe = getVarRequired(data_px, 'Fe_cation_apfu', mineralLabel);
px.Mg = getVarRequired(data_px, 'Mg_cation_apfu', mineralLabel);
px.Ca = getVarRequired(data_px, 'Ca_cation_apfu', mineralLabel);
px.Na = getVarRequired(data_px, 'Na_cation_apfu', mineralLabel);
px.Mn = getVarRequired(data_px, 'Mn_cation_apfu', mineralLabel);
px.Ti = getVarRequired(data_px, 'Ti_cation_apfu', mineralLabel);
px.Cr = getVarRequired(data_px, 'Cr_cation_apfu', mineralLabel);

px.K = getVarOptional(data_px, 'K_cation_apfu', 0, mineralLabel);

px.cationSum = ...
    px.Si + px.Al + px.Fe + px.Mg + px.Ca ...
    + px.Na + px.K + px.Mn + px.Ti + px.Cr;

end

function comp = calcOpxComponents(opx)
% calcOpxComponents
% Calculate Orthopyroxene components following Putirka (2008), Table 2.
% NaN values are preserved; only finite negative normative components are
% truncated to zero.

comp = struct();

comp.AlIV = truncateNegativePreserveNaN(2 - opx.Si);
comp.AlVI = truncateNegativePreserveNaN(opx.Al - comp.AlIV);

comp.XNaAlSi2O6 = minPreserveNaN(comp.AlVI, opx.Na);
comp.XFmTiAlSiO6 = opx.Ti;
comp.XCrAl2SiO6 = opx.Cr;

comp.XFmAl2SiO6 = truncateNegativePreserveNaN( ...
    comp.AlVI - comp.XNaAlSi2O6 - comp.XCrAl2SiO6);

comp.XCaFmSi2O6 = opx.Ca;

comp.XEnFs = truncateNegativePreserveNaN( ...
    ((opx.Fe + opx.Mn + opx.Mg) ...
    - comp.XFmTiAlSiO6 ...
    - comp.XFmAl2SiO6 ...
    - comp.XCaFmSi2O6) ./ 2);

end

function comp = calcCpxComponents(cpx)
% calcCpxComponents
% Calculate Clinopyroxene components following Putirka (2008), Table 3.
% NaN values are preserved; only finite negative normative components are
% truncated to zero.

comp = struct();

comp.AlIV = truncateNegativePreserveNaN(2 - cpx.Si);
comp.AlVI = truncateNegativePreserveNaN(cpx.Al - comp.AlIV);

comp.XFe3 = truncateNegativePreserveNaN( ...
    cpx.Na + comp.AlIV - comp.AlVI - 2 .* cpx.Ti - cpx.Cr);

comp.XJd = truncateNegativePreserveNaN( ...
    minPreserveNaN(comp.AlVI, cpx.Na));

comp.XCaTs = truncateNegativePreserveNaN(comp.AlVI - comp.XJd);

if isnan(comp.AlIV) || isnan(comp.XCaTs)
    comp.XCaTi = NaN;
elseif comp.AlIV > comp.XCaTs
    comp.XCaTi = truncateNegativePreserveNaN( ...
        (comp.AlIV - comp.XCaTs) ./ 2);
else
    comp.XCaTi = 0;
end

comp.XCrCaTs = truncateNegativePreserveNaN(cpx.Cr ./ 2);

comp.XDiHd = truncateNegativePreserveNaN( ...
    cpx.Ca - comp.XCaTi - comp.XCaTs - comp.XCrCaTs);

comp.XEnFs = truncateNegativePreserveNaN( ...
    (cpx.Fe + cpx.Mg - comp.XDiHd) ./ 2);

denom_Fm = cpx.Mg + cpx.Fe + cpx.Mn;
if isfinite(denom_Fm) && denom_Fm > 0 && isfinite(cpx.Mg)
    comp.XDi = (cpx.Mg ./ denom_Fm) .* comp.XDiHd;
else
    comp.XDi = NaN;
end

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify Equation (37) terms that prevent a finite temperature. Component
% terms are pressure-independent, whereas denom_T is checked over all rows.

maxTerms = 10;
termBuffer = strings(maxTerms, 1);
nTerms = 0;

if ~isfinite(row.XEnFs_opx(1)) || row.XEnFs_opx(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XEnFs_opx";
end
if ~isfinite(row.XEnFs_cpx(1)) || row.XEnFs_cpx(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XEnFs_cpx";
end
if ~isfinite(row.XEnFs_ratio_cpx_opx(1)) || ...
        row.XEnFs_ratio_cpx_opx(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XEnFs_cpx/XEnFs_opx";
end
if ~isfinite(row.XMg_cpx(1)) || row.XMg_cpx(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XMg_cpx";
end
if ~isfinite(row.ln_XMg_cpx(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "ln(XMg_cpx)";
end
if ~isfinite(row.MgNum_cpx(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "MgNum_cpx";
end
if ~isfinite(row.XFmAl2SiO6_opx(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XFmAl2SiO6_opx";
end
if ~isfinite(row.XDi_cpx(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XDi_cpx";
end
if ~isfinite(row.KD_FeMg_cpx_opx(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "KD(Fe-Mg)cpx-opx";
end
if any(~isfinite(row.denom_T) | row.denom_T <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "Equation (37) denominator";
end

invalidTerms = termBuffer(1:nTerms);

end

function validateRequiredColumns(tbl, mineralLabel)
% validateRequiredColumns
% Ensure that all normalized cation columns required for Equation (37) and
% its Table-2/Table-3 component calculations are present.

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    if ~ismember(variableName, tbl.Properties.VariableNames)
        error('%s table must contain variable: %s', ...
            mineralLabel, variableName);
    end
end

end

function value = getVarRequired(tbl, variableName, mineralLabel)
% getVarRequired
% Read a required numeric scalar while retaining NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);
validateScalarVariable(value, mineralLabel, variableName);

end

function value = getVarOptional( ...
        tbl, variableName, defaultValue, mineralLabel)
% getVarOptional
% Read an optional numeric scalar. An absent column receives defaultValue;
% an explicitly stored NaN remains NaN.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);
    validateScalarVariable(value, mineralLabel, variableName);
else
    value = defaultValue;
end

end

function validateScalarVariable(value, mineralLabel, variableName)
% validateScalarVariable
% Require one numeric scalar. NaN is allowed; negative finite values and Inf
% are handled by validateNonNegativeInputs.

if ~isnumeric(value) || ~isscalar(value)
    error('%s variable %s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end

function value = safeLogPositive(value)
% safeLogPositive
% Return ln(value) only when value is finite and strictly positive.

if isfinite(value) && value > 0
    value = log(value);
else
    value = NaN;
end

end

function value = truncateNegativePreserveNaN(value)
% truncateNegativePreserveNaN
% Preserve NaN; truncate only finite negative normative components to zero.

if isnan(value)
    return
end
if isfinite(value) && value < 0
    value = 0;
end

end

function value = minPreserveNaN(a, b)
% minPreserveNaN
% Return NaN if either argument is NaN; otherwise return min(a,b).

if isnan(a) || isnan(b)
    value = NaN;
else
    value = min(a, b);
end

end

function printTemperatureSummary( ...
        selectedCode_opx, selectedCode_cpx, temperatureValues)
% printTemperatureSummary
% Display either one temperature or the first-to-last range.

if isscalar(temperatureValues)
    disp([char(selectedCode_opx) ' & ' char(selectedCode_cpx) ': ' ...
        num2str(temperatureValues) ' degreeC']);
else
    disp([char(selectedCode_opx) ' & ' char(selectedCode_cpx) ': ' ...
        num2str(temperatureValues(1)) ' to ' ...
        num2str(temperatureValues(end)) ' degreeC']);
end

end

function printTemperatureRangeWarning( ...
        temperatureValues, minimumTemperature, maximumTemperature, ...
        selectedCode_opx, selectedCode_cpx)
% printTemperatureRangeWarning
% Warn when finite temperatures lie outside the approximate Figure-10
% comparison envelope. Results are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);
    fprintf(2, ...
        ['WARNING: Calculated Equation (37) temperature is outside the ' ...
         'approximate %.4g-%.4g degreeC comparison envelope shown in ' ...
         'Putirka (2008), Figure 10. This interval is a graphical ' ...
         'experimental envelope, not a strict universal calibration ' ...
         'boundary. %d of %d finite temperature point(s) are outside the ' ...
         'envelope; calculated finite range = %.4g-%.4g degreeC for %s & ' ...
         '%s.\n'], ...
        minimumTemperature, ...
        maximumTemperature, ...
        sum(outsideMask), ...
        sum(finiteMask), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(selectedCode_opx), ...
        char(selectedCode_cpx));
end

end

function printNonFiniteTemperatureWarning( ...
        temperatureValues, selectedCode_opx, selectedCode_cpx)
% printNonFiniteTemperatureWarning
% Report NaN and Inf temperatures without modifying or deleting them.

invalidMask = ~isfinite(temperatureValues);
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Equation (37) temperature values were ' ...
         'calculated for %s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the ' ...
         'calculation has not been stopped.\n'], ...
        char(selectedCode_opx), ...
        char(selectedCode_cpx), ...
        sum(invalidMask), ...
        numel(temperatureValues), ...
        sum(isnan(temperatureValues)), ...
        sum(isinf(temperatureValues)));
end

end
