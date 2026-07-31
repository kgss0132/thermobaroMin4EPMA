function results = RoederEmslie1970(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Olivine_Liquidus/RoederEmslie1970.m
% Tested with MATLAB R2024b
%
% Olivine-Liquid equilibrium thermometer / saturation relations
% Roeder, P.L., Emslie, R.F. (1970)
% Contributions to Mineralogy and Petrology, 29, 275–289
% DOI: https://doi.org/10.1007/BF00371276
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis selected from the
% input table with the first row of a Liquid dataset loaded by
% liquid.readLiquidExcel, and calculates temperatures from Eqs. (11), (12),
% and/or (13) of Roeder and Emslie (1970).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Olivine-Liquid pair, and appends the
% result blocks into a single output table.
%
% P_kbar may be a scalar or vector, allowing this function to be called from
% both startThermoCalc_fixedP and startThermoCalc_rangeP. The published
% equations contain no pressure term; therefore, the same calculated
% temperature is repeated for every supplied pressure value. Pressure is
% retained in the output for interface consistency and traceability.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Roeder and Emslie (1970) experimentally investigated olivine-liquid
% equilibria using three natural basaltic compositions over the following
% directly calibrated conditions:
%
%   Temperature : 1150–1300 degreeC
%   Pressure    : 1 atm total pressure (0.00101325 kbar)
%   Composition : natural olivine tholeiite, alkali olivine basalt, and
%                 ankaramite liquids
%   Redox range : approximately fO2 = 10^(-0.68) to 10^(-12) atm
%
% The temperature, pressure, starting-composition, and oxygen-fugacity
% ranges are summarized in the abstract on p. 275 and described in the
% Experimental Technique section on pp. 277–278. Eqs. (11)–(13) are given
% on p. 282.
%
% High-pressure experimental results at 4.5 and 9.0 kbar were plotted for
% comparison, but were not included in the regressions used to calculate
% Eqs. (11) and (12) (p. 282). Consequently, the equations are not directly
% pressure-calibrated beyond 1 atm and contain no pressure correction.
%
% Important application limitations stated or demonstrated in the paper:
%   - The selected olivine and liquid must represent crystal-liquid
%     equilibrium rather than later crystal-crystal re-equilibration
%     (p. 284).
%   - Liquid FeO must represent ferrous iron (Fe2+) expressed as FeO.
%     Electron-microprobe total iron should not be used without an
%     appropriate Fe2+/Fe3+ correction (pp. 280–281).
%   - The farther the temperature and liquid composition are extrapolated
%     from the experimental conditions, the greater the expected error
%     (p. 286).
%   - Very Fe-rich liquids and fayalitic olivines were not adequately tested
%     in the original calibration (pp. 283–285).
%   - For five Hawaiian natural liquids, calculated temperatures were
%     approximately 45–65 degreeC higher than observed temperatures
%     (p. 288). The relations should therefore be regarded as a first-order
%     estimate of the olivine saturation surface rather than a universally
%     precise thermometer.
%   - Eq. (13) has very weak temperature sensitivity. Experimental KD values
%     ranged from approximately 0.26 to 0.36, and the calculated change in KD
%     between 1000 and 1400 degreeC was only approximately 0.31 to 0.29
%     (p. 283). T13 should therefore be interpreted with particular caution.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) an input pressure differs from the directly calibrated pressure of
%      1 atm, or
%   2) a finite calculated temperature is outside 1150–1300 degreeC.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%
% The FIRST column of the Olivine table is treated as an identifier
% ("data code") displayed in the selection dialog.
%
% Olivine composition is retrieved using the first available format below,
% in this priority order:
%   1) Mg_cation_apfu and Fe_cation_apfu
%   2) XMg
%   3) Fo
%   4) MgO and FeO, or MgO and FeOt
%
% The Liquid dataset must contain MgO and FeO columns. Accepted column-name
% styles include MgO, MgOvalue, "MgO value", MgO_value, and corresponding
% FeO forms. Liquid FeO must be ferrous iron expressed as FeO.
%
% All finite composition values actually used by the calculation must be
% strictly greater than zero. XMg and converted Fo fractions must also be
% strictly less than one. NaN values are retained as missing values,
% propagated through the calculation, and reported by non-stopping fprintf
% warnings. NaN values are never replaced by zero.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Eq. (11):
%   log10(XMgO_ol / XMgO_liq) = 3740 / T(K) - 1.87
%
%   T11(K) = 3740 / [log10(XMgO_ol / XMgO_liq) + 1.87]
%
% Eq. (12):
%   log10(XFeO_ol / XFeO_liq) = 3911 / T(K) - 2.50
%
%   T12(K) = 3911 / [log10(XFeO_ol / XFeO_liq) + 2.50]
%
% Eq. (13):
%   log10(KD) = 171 / T(K) - 0.63
%
%   KD = (XFeO_ol / XMgO_ol) / (XFeO_liq / XMgO_liq)
%   T13(K) = 171 / [log10(KD) + 0.63]
%
% Temperatures are converted using:
%   T(degreeC) = T(K) - 273.15
%
% Notes:
% - Base-10 logarithms are used, following the published equations.
% - P_kbar is not used in Eqs. (11)–(13).
% - For backward compatibility with the original implementation, liquid MgO
%   and FeO proportions are normalized on an MgO + FeO basis before the
%   ratios are calculated.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = RoederEmslie1970(rawdata_struct, P_kbar)
%   results = RoederEmslie1970(rawdata_struct, P_kbar, ...
%       'Equation', equationOption)
%
% Inputs:
%   rawdata_struct : struct containing an Olivine table (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Name-value option:
%   'Equation' : 11, 12, 13, or 'both' (default: 'both')
%                'both' retains the behavior of the original script and
%                returns T11_C, T12_C, and T13_C.
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Olivine-Liquid pair.
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or an
% invalid pressure array.
if nargin < 2
    error('RoederEmslie1970 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% Options
% Validate the equation selector once before entering the interactive loop.
ip = inputParser;
ip.addParameter('Equation', 'both', ...
    @(x) isnumeric(x) || ischar(x) || isstring(x));
ip.parse(varargin{:});
equationOption = normalizeEquationOption(ip.Results.Equation);

%% 1) Retrieve datasets and constants
% Extract the Olivine table, load molecular-weight information, and read the
% Liquid dataset. The first Liquid row is used, matching the original script.
disp('=== Step 1: Preparing Olivine and Liquid datasets ===');

if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end

dataset_ol = rawdata_struct.Olivine;
MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();

if isempty(liqAll)
    error('Selected Liquid dataset is empty.');
end

if height(liqAll) > 1
    fprintf(2, ['WARNING: The selected Liquid dataset contains %d rows. ' ...
        'Only the first row will be used, following the original ' ...
        'RoederEmslie1970 behavior.\n'], height(liqAll));
end

selectedIdx_liq = 1;
selectedData_liq = liqAll(selectedIdx_liq, :);

disp('=== Preparing Olivine and Liquid datasets has been finished ===');

%% 2) Initialize output container and calibration checks
% Repeated table concatenation inside the interactive loop is avoided. Each
% result is stored as one table block in a preallocated cell buffer and all
% blocks are concatenated once after the loop finishes.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct experimental calibration limits from Roeder and Emslie (1970).
calibrationT_min_degC = 1150;
calibrationT_max_degC = 1300;
calibrationP_kbar = 0.00101325;  % 1 atm in kbar
pressureTolerance_kbar = 1.0e-6;

% Pressure is common to all selected Olivine rows in this function call.
% Therefore, the pressure warning is printed only once, after the first
% calculation result is displayed.
pressureOutsideCalibration = ...
    abs(P_kbar - calibrationP_kbar) > pressureTolerance_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–4) Interactive selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Olivine) ===');

while true
    % ----- Olivine selection -----
    % Assumption: the first column stores an identifier displayed to the user.
    dataCodes_ol = dataset_ol{:, 1};

    [selectedIdx_ol, ok] = listdlg( ...
        'PromptString', 'Please select the Olivine data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_ol, ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_ol)
        disp('Selection canceled');
        break;
    end

    selectedCode_ol = dataCodes_ol(selectedIdx_ol);
    selectedData_ol = dataset_ol(selectedIdx_ol, :);

    disp(['Olivine selected: ' char(string(selectedCode_ol))]);
    disp('Liquid selected: (auto) Row 1');

    % ----- Input checks -----
    % NaN values are recorded for later warning messages but do not stop the
    % calculation. Finite zero/negative values, fractions outside their
    % mathematical range, and infinite values stop the calculation.
    nanInputNames = findNaNInputs(selectedData_ol, selectedData_liq);
    validatePositiveInputs(selectedData_ol, selectedData_liq);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    row = calcTemp( ...
        selectedData_ol, selectedData_liq, P_kbar, MWinfo, equationOption);

    % Store identifiers with one value per pressure row.
    nRows = height(row);
    row.dataCode_ol = repmat(string(selectedCode_ol), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = localAttachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_ol', 'dataRow_liq'}, 'Before', 1);

    % Store this result as one table block. The buffer is doubled only when
    % its current capacity is exhausted, avoiding growth on every iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % ----- Result display -----
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    printTemperatureSummary(row, selectedCode_ol);

    % Warn once when any input pressure differs from the direct 1 atm
    % calibration pressure. The calculation is not stopped because pressure
    % is not included in the published equations.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the direct experimental calibration ' ...
             'pressure of Roeder and Emslie (1970): 1 atm ' ...
             '(0.00101325 kbar). %d of %d pressure point(s) differ from ' ...
             '1 atm; input range = %.6g–%.6g kbar. The published equations ' ...
             'contain no pressure correction, so pressure does not change ' ...
             'the calculated temperature.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn separately for each returned equation when finite temperatures
    % lie outside the direct experimental range of 1150–1300 degreeC.
    printTemperatureRangeWarnings( ...
        row, selectedCode_ol, calibrationT_min_degC, calibrationT_max_degC);

    % Report NaN inputs immediately after the temperature result. NaN values
    % remain NaN in the output and are never converted to zero.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         The calculation was continued, NaN values were retained, ' ...
             'and affected temperature results may be NaN.\n'], ...
            char(string(selectedCode_ol)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain a result-based check for NaN/Inf values caused by any numerical
    % issue, including an explicitly NaN input or an invalid equation result.
    printNonFiniteResultWarnings(row, selectedCode_ol);

    disp('--------------------------------------------------');

    % Ask whether to repeat using another Olivine analysis and the same
    % Liquid dataset.
    userAction = questdlg( ...
        'Continue with another data selection (same Liquid dataset)?', ...
        'RoederEmslie1970', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered result blocks only once after all selections have
% been completed. Return an empty table if no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

% Preserve metadata returned by liquid.readLiquidExcel.
results.Properties.UserData = struct('liquid', metaLiq);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function equationOption = normalizeEquationOption(rawOption)
% normalizeEquationOption
% Convert the Equation name-value input to one validated string value.

if isnumeric(rawOption)
    if ~isscalar(rawOption) || ~isfinite(rawOption) || ...
            ~ismember(rawOption, [11, 12, 13])
        error('Equation must be 11, 12, 13, or ''both''.');
    end
    equationOption = string(rawOption);
    return;
end

if ~(ischar(rawOption) || (isstring(rawOption) && isscalar(rawOption)))
    error('Equation must be 11, 12, 13, or ''both''.');
end

equationOption = lower(strtrim(string(rawOption)));
validOptions = ["11", "12", "13", "both"];

if ~any(equationOption == validOptions)
    error('Equation must be 11, 12, 13, or ''both''.');
end

end

function nanInputNames = findNaNInputs(data_olivine, data_liquid)
% findNaNInputs
% Return the names of required inputs that contain NaN. This function does
% not throw an error for NaN values.

[~, olivineNames, olivineValues] = localGetOlivineInputValues(data_olivine);
[liquidNames, liquidValues] = localGetLiquidInputValues(data_liquid);

allNames = [olivineNames; liquidNames];
allValues = [olivineValues; liquidValues];
nanMask = isnan(allValues);

nanInputNames = allNames(nanMask);

end

function validatePositiveInputs(data_olivine, data_liquid)
% validatePositiveInputs
% Stop the calculation when an input used by the thermometer is infinite,
% zero, negative, or outside a required fractional range. NaN is deliberately
% allowed so that it propagates and is reported by fprintf warnings.

[olivineMode, olivineNames, olivineValues] = ...
    localGetOlivineInputValues(data_olivine);
[liquidNames, liquidValues] = localGetLiquidInputValues(data_liquid);

allNames = [olivineNames; liquidNames];
allValues = [olivineValues; liquidValues];

infMask = isinf(allValues);
if any(infMask)
    error(['RoederEmslie1970: infinite composition value(s) were found in: ' ...
        char(strjoin(allNames(infMask), ', ')) '.']);
end

nonPositiveMask = isfinite(allValues) & allValues <= 0;
if any(nonPositiveMask)
    error(['RoederEmslie1970: composition values used by the calculation ' ...
        'must be > 0. Zero or negative finite value(s) were found in: ' ...
        char(strjoin(allNames(nonPositiveMask), ', ')) '.']);
end

% XMg and Fo are fractions after Fo conversion. A value of one would make
% the complementary Fe fraction zero, which is invalid for Eqs. (12)–(13).
if olivineMode == "XMg" || olivineMode == "Fo"
    fractionValue = olivineValues(1);
    if isfinite(fractionValue) && fractionValue >= 1
        error(['RoederEmslie1970: ' char(olivineNames(1)) ...
            ' must be > 0 and < 1 after conversion to a fraction.']);
    end
end

end

function row = calcTemp(data_olivine, data_liquid, P_kbar, MWinfo, equationOption)
% calcTemp
% Compute temperature estimates for one Olivine row and one Liquid row.
% One output row is returned for each supplied pressure value. Because the
% published equations contain no pressure term, composition-dependent values
% and temperatures are repeated across the pressure vector.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

% ----- Olivine composition -----
[XMg_ol_scalar, XFe_ol_scalar] = ...
    localGetOlMgFeFractions(data_olivine, MWinfo);

% ----- Liquid composition -----
% FeO must represent ferrous iron expressed as FeO.
MgO_liq_wt_scalar = localGetLiqOxRequired(data_liquid, 'MgO');
FeO_liq_wt_scalar = localGetLiqOxRequired(data_liquid, 'FeO');

% Mole proportions of MgO and FeO in the Liquid dataset. NaN values remain
% NaN and propagate through every subsequent operation.
nMg_liq_scalar = MgO_liq_wt_scalar .* ...
    MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
nFe_liq_scalar = FeO_liq_wt_scalar .* ...
    MWinfo.Cat.FeO ./ MWinfo.MW.FeO;

den_liq_scalar = nMg_liq_scalar + nFe_liq_scalar;
XMgO_liq_scalar = nMg_liq_scalar ./ den_liq_scalar;
XFeO_liq_scalar = nFe_liq_scalar ./ den_liq_scalar;

% For the Mg-Fe binary Olivine approximation used in the original script,
% XMg and XFe are treated as the MgO and FeO Olivine component fractions.
XMgO_ol_scalar = XMg_ol_scalar;
XFeO_ol_scalar = XFe_ol_scalar;

% ----- Intermediate ratios -----
ratio11_scalar = XMgO_ol_scalar ./ XMgO_liq_scalar;
ratio12_scalar = XFeO_ol_scalar ./ XFeO_liq_scalar;
KD_scalar = (XFeO_ol_scalar ./ XMgO_ol_scalar) ./ ...
    (XFeO_liq_scalar ./ XMgO_liq_scalar);

% ----- Temperature equations -----
% Direct evaluation intentionally allows NaN to propagate. localKToC changes
% non-positive or non-finite Kelvin results to NaN without stopping.
T11_K_scalar = 3740 ./ (log10(ratio11_scalar) + 1.87);
T12_K_scalar = 3911 ./ (log10(ratio12_scalar) + 2.50);
T13_K_scalar = 171 ./ (log10(KD_scalar) + 0.63);

T11_C_scalar = localKToC(T11_K_scalar);
T12_C_scalar = localKToC(T12_K_scalar);
T13_C_scalar = localKToC(T13_K_scalar);

% ----- Pack outputs -----
row = table();
row.P_kbar = P_kbar;

row.XMg_ol = repmat(XMg_ol_scalar, nP, 1);
row.XFe_ol = repmat(XFe_ol_scalar, nP, 1);
row.MgO_liq = repmat(MgO_liq_wt_scalar, nP, 1);
row.FeO_liq = repmat(FeO_liq_wt_scalar, nP, 1);
row.XMgO_liq = repmat(XMgO_liq_scalar, nP, 1);
row.XFeO_liq = repmat(XFeO_liq_scalar, nP, 1);
row.XMgO_ol = repmat(XMgO_ol_scalar, nP, 1);
row.XFeO_ol = repmat(XFeO_ol_scalar, nP, 1);
row.ratio11_XMg = repmat(ratio11_scalar, nP, 1);
row.ratio12_XFe = repmat(ratio12_scalar, nP, 1);
row.KD_FeMg_ol_liq = repmat(KD_scalar, nP, 1);
row.T11_K = repmat(T11_K_scalar, nP, 1);
row.T12_K = repmat(T12_K_scalar, nP, 1);
row.T13_K = repmat(T13_K_scalar, nP, 1);

switch equationOption
    case "11"
        row.T11_C = repmat(T11_C_scalar, nP, 1);
    case "12"
        row.T12_C = repmat(T12_C_scalar, nP, 1);
    case "13"
        row.T13_C = repmat(T13_C_scalar, nP, 1);
    otherwise
        row.T11_C = repmat(T11_C_scalar, nP, 1);
        row.T12_C = repmat(T12_C_scalar, nP, 1);
        row.T13_C = repmat(T13_C_scalar, nP, 1);
end

end

function T_C = localKToC(T_K)
% localKToC
% Convert Kelvin to degreeC. Non-positive or non-finite Kelvin results are
% returned as NaN so that invalid numerical results remain explicit.

if isnan(T_K)
    T_C = NaN;
elseif ~isfinite(T_K) || T_K <= 0
    T_C = NaN;
else
    T_C = T_K - 273.15;
end

end

function row = localAttachLiquidIDs(row, data_liquid)
% localAttachLiquidIDs
% Copy optional Liquid identifiers into every pressure row of the result.

nRows = height(row);
variableNames = data_liquid.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    rawValue = data_liquid.('Index');
    if isnumeric(rawValue) || islogical(rawValue)
        row.liq_Index = repmat(rawValue(1), nRows, 1);
    else
        row.liq_Index = repmat(string(rawValue(1)), nRows, 1);
    end
end

if any(strcmp(variableNames, 'Experiment'))
    rawExperiment = data_liquid.('Experiment');
    row.liq_Experiment = repmat(string(rawExperiment(1)), nRows, 1);
end

if any(strcmp(variableNames, 'Citation'))
    rawCitation = data_liquid.('Citation');
    row.liq_Citation = repmat(string(rawCitation(1)), nRows, 1);
end

end

function [XMg_ol, XFe_ol] = localGetOlMgFeFractions(data_olivine, MWinfo)
% localGetOlMgFeFractions
% Retrieve Olivine Mg and Fe fractions from the first supported input format.
% NaN values are retained and propagated.

[olivineMode, ~, olivineValues] = ...
    localGetOlivineInputValues(data_olivine);

switch olivineMode
    case "apfu"
        Mg = olivineValues(1);
        Fe = olivineValues(2);
        denominator = Mg + Fe;
        XMg_ol = Mg ./ denominator;
        XFe_ol = Fe ./ denominator;

    case "XMg"
        XMg_ol = olivineValues(1);
        XFe_ol = 1 - XMg_ol;

    case "Fo"
        XMg_ol = olivineValues(1);
        XFe_ol = 1 - XMg_ol;

    case "oxide"
        MgO = olivineValues(1);
        FeO = olivineValues(2);
        nMg = MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
        nFe = FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
        denominator = nMg + nFe;
        XMg_ol = nMg ./ denominator;
        XFe_ol = nFe ./ denominator;

    otherwise
        error('Unsupported Olivine input mode.');
end

end

function [olivineMode, inputNames, inputValues] = ...
        localGetOlivineInputValues(data_olivine)
% localGetOlivineInputValues
% Identify the preferred available Olivine input format and return its names
% and scalar numeric values. Fo percentages are converted to fractions.

variableNames = data_olivine.Properties.VariableNames;

hasMgApfu = any(strcmp(variableNames, 'Mg_cation_apfu'));
hasFeApfu = any(strcmp(variableNames, 'Fe_cation_apfu'));

if hasMgApfu && hasFeApfu
    olivineMode = "apfu";
    inputNames = ["Olivine.Mg_cation_apfu"; "Olivine.Fe_cation_apfu"];
    inputValues = [ ...
        localToScalarDouble(data_olivine.Mg_cation_apfu, NaN); ...
        localToScalarDouble(data_olivine.Fe_cation_apfu, NaN)];
    return;
end

if any(strcmp(variableNames, 'XMg'))
    olivineMode = "XMg";
    inputNames = "Olivine.XMg";
    inputValues = localToScalarDouble(data_olivine.XMg, NaN);
    return;
end

if any(strcmp(variableNames, 'Fo'))
    olivineMode = "Fo";
    inputNames = "Olivine.Fo";
    inputValues = localToScalarDouble(data_olivine.Fo, NaN);
    if isfinite(inputValues) && inputValues > 1.5
        inputValues = inputValues ./ 100;
    end
    return;
end

hasMgO = any(strcmp(variableNames, 'MgO'));
hasFeO = any(strcmp(variableNames, 'FeO'));
hasFeOt = any(strcmp(variableNames, 'FeOt'));

if hasMgO && (hasFeO || hasFeOt)
    olivineMode = "oxide";

    if hasFeO
        ironVariableName = 'FeO';
    else
        ironVariableName = 'FeOt';
    end

    inputNames = ["Olivine.MgO"; "Olivine." + string(ironVariableName)];
    inputValues = [ ...
        localToScalarDouble(data_olivine.MgO, NaN); ...
        localToScalarDouble(data_olivine.(ironVariableName), NaN)];
    return;
end

error(['Selected Olivine row lacks a supported composition format. ' ...
    'Required alternatives are Mg_cation_apfu + Fe_cation_apfu, XMg, Fo, ' ...
    'or MgO + FeO/FeOt.']);

end

function [inputNames, inputValues] = localGetLiquidInputValues(data_liquid)
% localGetLiquidInputValues
% Retrieve the Liquid MgO and ferrous-FeO values used by the equations.

mgColumnName = localFindOxideColumn( ...
    data_liquid.Properties.VariableNames, 'MgO');
feColumnName = localFindOxideColumn( ...
    data_liquid.Properties.VariableNames, 'FeO');

if isempty(mgColumnName)
    error('Selected Liquid row lacks required oxide column: MgO');
end
if isempty(feColumnName)
    error(['Selected Liquid row lacks required oxide column: FeO. ' ...
        'FeO must represent ferrous iron expressed as FeO.']);
end

inputNames = ["Liquid." + string(mgColumnName); ...
              "Liquid." + string(feColumnName)];
inputValues = [ ...
    localToScalarDouble(data_liquid.(mgColumnName), NaN); ...
    localToScalarDouble(data_liquid.(feColumnName), NaN)];

end

function value = localGetLiqOxRequired(data_liquid, oxide)
% localGetLiqOxRequired
% Return one required Liquid oxide value. Missing columns stop calculation;
% NaN values are returned unchanged and handled by non-stopping warnings.

columnName = localFindOxideColumn( ...
    data_liquid.Properties.VariableNames, oxide);

if isempty(columnName)
    error('Selected Liquid row lacks required oxide column: %s', oxide);
end

value = localToScalarDouble(data_liquid.(columnName), NaN);

end

function columnName = localFindOxideColumn(variableNames, oxide)
% localFindOxideColumn
% Match oxide columns after removing spaces, underscores, and hyphens.

canonicalNames = cell(size(variableNames));

for i = 1:numel(variableNames)
    currentName = lower(variableNames{i});
    currentName = strrep(currentName, ' ', '');
    currentName = strrep(currentName, '_', '');
    currentName = strrep(currentName, '-', '');
    canonicalNames{i} = currentName;
end

canonicalOxide = lower(oxide);
canonicalOxide = strrep(canonicalOxide, ' ', '');
canonicalOxide = strrep(canonicalOxide, '_', '');
canonicalOxide = strrep(canonicalOxide, '-', '');

targets = {[canonicalOxide 'value'], canonicalOxide};
columnName = '';

for i = 1:numel(targets)
    matchedIndex = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(matchedIndex)
        columnName = variableNames{matchedIndex};
        return;
    end
end

end

function value = localToScalarDouble(rawValue, defaultValue)
% localToScalarDouble
% Convert the first element of a table value to double. Empty, missing, or
% non-convertible content returns defaultValue, which is NaN in this script.

value = defaultValue;

if isempty(rawValue)
    return;
end

if isnumeric(rawValue) || islogical(rawValue)
    value = double(rawValue(1));
    return;
end

if isstring(rawValue)
    if ismissing(rawValue(1))
        return;
    end
    convertedValue = str2double(rawValue(1));
    if ~isnan(convertedValue)
        value = convertedValue;
    end
    return;
end

if ischar(rawValue)
    convertedValue = str2double(string(rawValue));
    if ~isnan(convertedValue)
        value = convertedValue;
    end
    return;
end

if iscell(rawValue)
    firstValue = rawValue{1};
    if isempty(firstValue)
        return;
    end

    if isnumeric(firstValue) || islogical(firstValue)
        value = double(firstValue(1));
        return;
    end

    if isstring(firstValue) || ischar(firstValue)
        convertedValue = str2double(string(firstValue));
        if ~isnan(convertedValue)
            value = convertedValue;
        end
    end
end

end

function printTemperatureSummary(row, selectedCode_ol)
% printTemperatureSummary
% Display the first-to-last value for each temperature equation returned.

variableNames = row.Properties.VariableNames;
temperatureVariables = {'T11_C', 'T12_C', 'T13_C'};
equationLabels = {'Eq11', 'Eq12', 'Eq13'};

for i = 1:numel(temperatureVariables)
    variableName = temperatureVariables{i};
    if any(strcmp(variableNames, variableName))
        temperatureValues = row.(variableName);
        if height(row) == 1
            disp([equationLabels{i} ' | ' char(string(selectedCode_ol)) ...
                ': ' num2str(temperatureValues(1)) ' degreeC']);
        else
            disp([equationLabels{i} ' | ' char(string(selectedCode_ol)) ...
                ': ' num2str(temperatureValues(1)) ' to ' ...
                num2str(temperatureValues(end)) ' degreeC']);
        end
    end
end

end

function printTemperatureRangeWarnings( ...
        row, selectedCode_ol, minimumTemperature, maximumTemperature)
% printTemperatureRangeWarnings
% Print one non-stopping calibration warning for each returned equation that
% contains a finite temperature outside the direct calibration range.

variableNames = row.Properties.VariableNames;
temperatureVariables = {'T11_C', 'T12_C', 'T13_C'};
equationLabels = {'Eq. (11)', 'Eq. (12)', 'Eq. (13)'};

for i = 1:numel(temperatureVariables)
    variableName = temperatureVariables{i};

    if ~any(strcmp(variableNames, variableName))
        continue;
    end

    temperatureValues = row.(variableName);
    finiteMask = isfinite(temperatureValues);
    outsideMask = finiteMask & ...
        (temperatureValues < minimumTemperature | ...
         temperatureValues > maximumTemperature);

    if any(outsideMask)
        finiteValues = temperatureValues(finiteMask);
        fprintf(2, ...
            ['WARNING: %s calculated temperature is outside the direct experimental ' ...
             'calibration range of Roeder and Emslie (1970): ' ...
             '1150–1300 degreeC. %d of %d finite temperature point(s) are ' ...
             'outside the range; calculated finite range = %.6g–%.6g degreeC ' ...
             'for %s.\n'], ...
            equationLabels{i}, ...
            sum(outsideMask), ...
            sum(finiteMask), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_ol)));
    end
end

end

function printNonFiniteResultWarnings(row, selectedCode_ol)
% printNonFiniteResultWarnings
% Print one non-stopping warning for each returned equation containing NaN or
% Inf. Values remain in the output table and calculation is not interrupted.

variableNames = row.Properties.VariableNames;
temperatureVariables = {'T11_C', 'T12_C', 'T13_C'};
equationLabels = {'Eq. (11)', 'Eq. (12)', 'Eq. (13)'};

for i = 1:numel(temperatureVariables)
    variableName = temperatureVariables{i};

    if ~any(strcmp(variableNames, variableName))
        continue;
    end

    temperatureValues = row.(variableName);
    invalidMask = ~isfinite(temperatureValues);

    if any(invalidMask)
        fprintf(2, ...
            ['WARNING: Non-finite %s temperature value(s) were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            equationLabels{i}, ...
            char(string(selectedCode_ol)), ...
            sum(invalidMask), ...
            numel(temperatureValues), ...
            sum(isnan(temperatureValues)), ...
            sum(isinf(temperatureValues)));
    end
end

end
