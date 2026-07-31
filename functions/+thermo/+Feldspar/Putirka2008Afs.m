function results = Putirka2008Afs(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Feldspar_Liquid/Putirka2008Afs.m
% Tested with MATLAB R2024b
%
% Alkali feldspar-liquid thermometer
% Putirka, K.D. (2008)
% Reviews in Mineralogy & Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one alkali feldspar analysis and pairs
% it with one liquid analysis, then calculates temperature using Equation
% (24b) of Putirka (2008).
%
% The function accepts pressure as either a scalar or vector. Therefore, it
% can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for each pressure value
% for every selected alkali feldspar-liquid pair.
%
% The function is designed for repeated calculations: after each run it asks
% whether another alkali feldspar analysis should be calculated and stores
% all result blocks in a preallocated cell buffer. The blocks are concatenated
% only once after the interactive loop finishes.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka (2008) calibrated the alkali feldspar-liquid thermometer in
% Equation (24b) using experiments on natural compositions conducted at:
%
%   Temperature : less than 1050 degreeC
%   Model error : +/-23 degreeC for the calibration dataset
%   Composition : natural-composition experimental systems
%   Volatiles   : nearly all calibration experiments were hydrous, although
%                 Equation (24b) does not contain an explicit H2O term
%
% Equation (24b) is presented on p. 79. The calibration-data description is
% given in the caption to Figure 5 on p. 80, and the calibration limit,
% uncertainty, high-temperature extrapolation problem, hydrous-data note,
% and equilibrium-test discussion are given on p. 82.
%
% Important application notes from Putirka (2008):
%   1) The thermometer was calibrated only at T < 1050 degreeC.
%      Application at T >= 1050 degreeC may produce strong systematic error.
%   2) The model does not perform well for the synthetic system of Holtz et
%      al. (2005). It should be applied most cautiously to compositions unlike
%      the natural-composition experiments used for calibration.
%   3) Alkali feldspar and liquid must have approached equilibrium. Putirka
%      (2008, p. 82) reports KD(An-Ab)afs-liq = 0.27 +/- 0.18 for 60
%      experiments as a possible, but imprecise, equilibrium test. Or-Ab
%      exchange is too dependent on T-P-H2O to serve as the same test.
%   4) All liquid and mineral components must be calculated in the same way
%      as in the original calibration. Component-calculation requirements are
%      described on pp. 68-71. Liquid cation fractions are calculated on an
%      anhydrous basis.
%   5) Putirka (2008) does not report a single numerical pressure calibration
%      range specifically for Equation (24b). Consequently, this implementation
%      does not invent pressure limits. Instead, it prints a non-stopping
%      fprintf warning stating that the supplied pressure cannot be screened
%      against a published numerical pressure range.
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) a finite calculated temperature is >= 1050 degreeC,
%   2) pressure applicability cannot be screened because no numerical range
%      is reported for Equation (24b),
%   3) an input used by the thermometer contains NaN, or
%   4) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Feldspar : table
%
% The FIRST column is treated as an identifier ("data code") displayed in
% the selection dialog. Oxide columns are located by names equivalent to the
% following, optionally with a "_value" suffix:
%
%   Alkali feldspar:
%     SiO2, Al2O3, CaO, Na2O, K2O
%
%   Liquid:
%     SiO2, TiO2, Al2O3, FeO or FeOt, MnO, MgO, CaO, Na2O, K2O,
%     V2O3, Cr2O3, NiO, P2O5, SO3, F, Cl, Fe2O3, H2O
%
% SiO2 and Al2O3 columns must exist in the feldspar table. Optional oxide
% columns that do not exist are treated as zero. However, when an existing
% input cell contains NaN, that NaN is retained rather than converted to
% zero. It therefore propagates through cation normalization and temperature
% calculation and is reported by a non-stopping fprintf warning.
%
% All finite oxide values used by the calculation must be greater than or
% equal to zero. A finite value below zero stops the calculation with an
% error. NaN is allowed and retained. Inf is not allowed.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Equation (24b):
%
%   10^4 / T(K) =
%       17.3
%     - 1.03 ln(KAb_afs_liq)
%     - 200 XCaO_liq
%     - 2.42 XNaO0.5_liq
%     - 29.8 XKO0.5_liq
%     + 13500 (XCaO_liq - 0.0037)^2
%     - 550 (XKO0.5_liq - 0.056)(XNaO0.5_liq - 0.089)
%     - 0.078 P(kbar)
%
% where:
%
%   XAb_afs = XNaO0.5_afs /
%             (XCaO_afs + XNaO0.5_afs + XKO0.5_afs)
%
%   KAb_afs_liq = XAb_afs /
%                 [XNaO0.5_liq * XAlO1.5_liq * (XSiO2_liq)^3]
%
%   T(K)       = 10^4 / denominator
%   T(degreeC) = T(K) - 273.15
%
% Mineral and liquid oxide abundances are converted to cation fractions
% using molecular weights and cation numbers supplied by
% liquid.getMolarWeights(). Liquid fractions are normalized on an anhydrous
% basis; H2O is retained in the output for reference but is not included in
% Equation (24b) or in the anhydrous cation-fraction denominator.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008Afs(rawdata_struct, P_kbar)
%   results = Putirka2008Afs(rawdata_struct, P_kbar, ...
%                               'LiquidRow', liquidRowNumber)
%
% Inputs:
%   rawdata_struct : struct containing table rawdata_struct.Feldspar
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Option:
%   'LiquidRow'    : positive integer scalar selecting a liquid table row.
%                    If omitted, the first liquid row is used.
%
% Output:
%   results : table containing one row per pressure value for every selected
%             alkali feldspar-liquid pair.
%

%% Input validation
% Basic argument checks prevent silent failures caused by missing inputs or
% invalid pressure values.
if nargin < 2
    error('Putirka2008Afs requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

if ~isfield(rawdata_struct, 'Feldspar') || ...
        ~istable(rawdata_struct.Feldspar)
    error('rawdata_struct must contain table: rawdata_struct.Feldspar');
end

ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && fix(x) == x));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve feldspar and liquid datasets
% The feldspar table is supplied in rawdata_struct. The liquid table and its
% metadata are read using the same project-level helper as the original file.
disp('=== Step 1: Preparing feldspar and liquid datasets ===');

dataset_afs = rawdata_struct.Feldspar;
MWinfo = liquid.getMolarWeights();

[liqAll, metaLiq] = liquid.readLiquidExcel();
if isempty(liqAll)
    error('Selected liquid dataset is empty.');
end

% Resolve the liquid row once because it is common to every feldspar
% selection in this function call.
if isempty(liquidRowOpt)
    if height(liqAll) > 1
        fprintf(2, ...
            ['WARNING: Liquid dataset has multiple rows (%d). ' ...
             'The first row will be used.\n'], ...
            height(liqAll));
    end
    idxLiq = 1;
else
    idxLiq = liquidRowOpt;
    if idxLiq > height(liqAll)
        error('Requested LiquidRow (%d) exceeds liquid rows (%d).', ...
            idxLiq, height(liqAll));
    end
end

selectedData_liq = liqAll(idxLiq, :);

disp('=== Preparing feldspar and liquid datasets has been finished ===');

%% 2) Initialize output container and calibration checks
% Repeated table concatenation inside the selection loop is avoided. Each
% result table is stored in a preallocated cell buffer and all blocks are
% concatenated once after the loop.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Putirka (2008, p. 82) states that Equation (24b) was calibrated at
% temperatures below 1050 degreeC. No explicit lower-temperature limit or
% numerical pressure range is reported specifically for Equation (24b).
calibrationT_max_degC = 1050;
pressureRangeWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop and calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Feldspar) ===');

while true
    % ----- Alkali feldspar selection -----
    % Assumption: the first table column stores the identifier shown to the
    % user in the selection dialog.
    dataCodes_afs = dataset_afs{:, 1};

    [selectedIdx_afs, ok] = listdlg( ...
        'PromptString', ...
        'Please select the Feldspar data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_afs, ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_afs)
        disp('Selection canceled');
        break;
    end

    selectedCode_afs = dataCodes_afs(selectedIdx_afs);
    selectedData_afs = dataset_afs(selectedIdx_afs, :);

    disp(['Feldspar selected: ' char(string(selectedCode_afs))]);
    disp(['Liquid selected: Row ' num2str(idxLiq)]);

    % ----- Input checks -----
    % Negative finite values and Inf are rejected. NaN is allowed so it can
    % propagate through the calculation and remain visible in the output.
    validateInputValues(selectedData_afs, selectedData_liq);
    nanInputNames = findNaNInputs(selectedData_afs, selectedData_liq);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    row = calcTemp(selectedData_afs, selectedData_liq, P_kbar, MWinfo);

    % Add identifiers to every pressure row for traceability.
    nRows = height(row);
    row.dataCode_afs = repmat(string(selectedCode_afs), nRows, 1);
    row.dataRow_liq = repmat(idxLiq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_afs', 'dataRow_liq'}, 'Before', 1);

    % Store this result block. The capacity is doubled only when necessary,
    % rather than changing the result table size on every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % ----- Immediate result display -----
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_afs)) ' & Liquid row ' ...
            num2str(idxLiq) ': ' num2str(row.T_Eq24b_C) ' degreeC']);
    else
        disp([char(string(selectedCode_afs)) ' & Liquid row ' ...
            num2str(idxLiq) ': ' num2str(row.T_Eq24b_C(1)) ' to ' ...
            num2str(row.T_Eq24b_C(end)) ' degreeC']);
    end

    % ----- Pressure applicability warning -----
    % Putirka (2008) gives no numerical pressure calibration interval for
    % Equation (24b), so an evidence-based outside-range test cannot be made.
    % Print this caution only once per function call.
    if ~pressureRangeWarningIssued
        fprintf(2, ...
            ['WARNING: Putirka (2008) does not report a numerical pressure ' ...
             'calibration range specifically for Equation (24b) (pp. 79-82). ' ...
             'The input pressure therefore cannot be screened against a ' ...
             'published pressure range; input range = %.4g-%.4g kbar.\n'], ...
            min(P_kbar), max(P_kbar));
        pressureRangeWarningIssued = true;
    end

    % ----- Temperature calibration warning -----
    % The published condition is T < 1050 degreeC. NaN and Inf are handled
    % separately by the non-finite-result warning below.
    finiteTemperature = isfinite(row.T_Eq24b_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        row.T_Eq24b_C >= calibrationT_max_degC;

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_Eq24b_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the published ' ...
             'calibration condition of Putirka (2008) Equation (24b): ' ...
             'T < 1050 degreeC (p. 82). %d of %d finite temperature ' ...
             'point(s) are outside the condition; calculated finite range ' ...
             '= %.4g-%.4g degreeC for %s & Liquid row %d. Strong systematic ' ...
             'error may occur at T >= 1050 degreeC.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_afs)), ...
            idxLiq);
    end

    % ----- NaN input warning -----
    % Existing NaN values are retained and never replaced by zero.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & ' ...
             'Liquid row %d: %s.\n' ...
             '         The NaN value(s) were retained and the calculation ' ...
             'continued; calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_afs)), ...
            idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % ----- Non-finite result warning -----
    % This also catches mathematically invalid zero-valued combinations, such
    % as a zero denominator in a component ratio, without discarding results.
    invalidTemperature = ~isfinite(row.T_Eq24b_C);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s ' ...
             '& Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_afs)), ...
            idxLiq, ...
            sum(invalidTemperature), ...
            numel(row.T_Eq24b_C), ...
            sum(isnan(row.T_Eq24b_C)), ...
            sum(isinf(row.T_Eq24b_C)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another feldspar analysis.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Putirka2008Afs', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered result blocks once. Return an empty table when no
% calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

% Preserve the liquid-file metadata used by the original implementation.
results.Properties.UserData = struct('liquid', metaLiq);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_afs, data_liq)
% findNaNInputs
% Return names of existing feldspar and liquid inputs that contain NaN.
% Missing optional columns are not listed because they are intentionally
% treated as zero. This function never throws an error for NaN.

feldsparOxides = {'SiO2', 'Al2O3', 'CaO', 'Na2O', 'K2O'};
liquidOxides = {'SiO2', 'TiO2', 'Al2O3', 'MnO', 'MgO', 'CaO', ...
    'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', ...
    'F', 'Cl', 'Fe2O3'};

maxNames = numel(feldsparOxides) + numel(liquidOxides) + 1;
nanNamesBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(feldsparOxides)
    oxide = feldsparOxides{i};
    name = findOxideColumn(data_afs.Properties.VariableNames, oxide);
    if ~isempty(name)
        value = toScalarDouble(data_afs.(name));
        if isnan(value)
            nNames = nNames + 1;
            nanNamesBuffer(nNames) = "Feldspar." + string(name);
        end
    end
end

% FeO and FeOt are alternatives. Use exactly the column that the calculation
% will use, preserving NaN when that existing column contains NaN.
feName = selectLiquidFeColumnName(data_liq.Properties.VariableNames);
if ~isempty(feName)
    feValue = toScalarDouble(data_liq.(feName));
    if isnan(feValue)
        nNames = nNames + 1;
        nanNamesBuffer(nNames) = "Liquid." + string(feName);
    end
end

for i = 1:numel(liquidOxides)
    oxide = liquidOxides{i};
    name = findOxideColumn(data_liq.Properties.VariableNames, oxide);
    if ~isempty(name)
        value = toScalarDouble(data_liq.(name));
        if isnan(value)
            nNames = nNames + 1;
            nanNamesBuffer(nNames) = "Liquid." + string(name);
        end
    end
end

nanInputNames = nanNamesBuffer(1:nNames);

end

function validateInputValues(data_afs, data_liq)
% validateInputValues
% Stop when an existing input used in cation calculations contains Inf or a
% finite value below zero. Zero is allowed. NaN is deliberately allowed and
% retained so that it propagates through the calculation.

feldsparOxides = {'SiO2', 'Al2O3', 'CaO', 'Na2O', 'K2O'};
liquidOxides = {'SiO2', 'TiO2', 'Al2O3', 'MnO', 'MgO', 'CaO', ...
    'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', ...
    'F', 'Cl', 'Fe2O3', 'H2O'};

maxNames = numel(feldsparOxides) + numel(liquidOxides) + 1;
invalidNamesBuffer = strings(maxNames, 1);
nInvalid = 0;

for i = 1:numel(feldsparOxides)
    oxide = feldsparOxides{i};
    name = findOxideColumn(data_afs.Properties.VariableNames, oxide);
    if ~isempty(name)
        value = toScalarDouble(data_afs.(name));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            invalidNamesBuffer(nInvalid) = "Feldspar." + string(name);
        end
    end
end

feName = selectLiquidFeColumnName(data_liq.Properties.VariableNames);
if ~isempty(feName)
    feValue = toScalarDouble(data_liq.(feName));
    if isinf(feValue) || (isfinite(feValue) && feValue < 0)
        nInvalid = nInvalid + 1;
        invalidNamesBuffer(nInvalid) = "Liquid." + string(feName);
    end
end

for i = 1:numel(liquidOxides)
    oxide = liquidOxides{i};
    name = findOxideColumn(data_liq.Properties.VariableNames, oxide);
    if ~isempty(name)
        value = toScalarDouble(data_liq.(name));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            invalidNamesBuffer(nInvalid) = "Liquid." + string(name);
        end
    end
end

if nInvalid > 0
    invalidInputNames = invalidNamesBuffer(1:nInvalid);
    error(['Putirka2008Afs: oxide inputs must be non-negative. ' ...
           'Negative finite or infinite value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_afs, data_liq, P_kbar, MWinfo)
% calcTemp
% Calculate Equation (24b) for one alkali feldspar row, one liquid row, and
% a scalar or vector of pressure values. Scalar composition values are
% expanded to nP rows before table construction, preventing dimension
% mismatches in range-pressure mode.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% Calculate feldspar and liquid cation fractions once for the selected pair.
afs = prepareFeldsparRow(data_afs, MWinfo);
liq = prepareLiquidRow(data_liq, MWinfo);

% Expand scalar composition terms to one value per pressure point.
row.XCaO_afs = repmat(afs.XCaO, nP, 1);
row.XNaO0_5_afs = repmat(afs.XNaO0_5, nP, 1);
row.XKO0_5_afs = repmat(afs.XKO0_5, nP, 1);
row.XAn_afs = repmat(afs.XAn, nP, 1);
row.XAb_afs = repmat(afs.XAb, nP, 1);
row.XOr_afs = repmat(afs.XOr, nP, 1);

row.XSiO2_liq = repmat(liq.XSiO2, nP, 1);
row.XAlO1_5_liq = repmat(liq.XAlO1_5, nP, 1);
row.XCaO_liq = repmat(liq.XCaO, nP, 1);
row.XNaO0_5_liq = repmat(liq.XNaO0_5, nP, 1);
row.XKO0_5_liq = repmat(liq.XKO0_5, nP, 1);
row.H2O_liq = repmat(liq.H2O, nP, 1);

KAb = afs.XAb ./ ...
    (liq.XNaO0_5 .* liq.XAlO1_5 .* (liq.XSiO2 .^ 3));
KAbVector = repmat(KAb, nP, 1);
row.KAb_afs_liq = KAbVector;

den24b = ...
    17.3 ...
    - 1.03 .* log(KAbVector) ...
    - 200 .* liq.XCaO ...
    - 2.42 .* liq.XNaO0_5 ...
    - 29.8 .* liq.XKO0_5 ...
    + 13500 .* ((liq.XCaO - 0.0037) .^ 2) ...
    - 550 .* (liq.XKO0_5 - 0.056) .* ...
              (liq.XNaO0_5 - 0.089) ...
    - 0.078 .* P_kbar;

row.den_Eq24b = den24b;
[row.T_Eq24b_K, row.T_Eq24b_C] = temperatureFrom10e4OverT(den24b);

end

function [T_K, T_C] = temperatureFrom10e4OverT(den)
% temperatureFrom10e4OverT
% Convert the vector 10^4/T denominator to Kelvin and degreeC. Invalid or
% non-positive denominator values are retained as NaN without stopping.

T_K = nan(size(den));
validDenominator = isfinite(den) & den > 0;
T_K(validDenominator) = 1e4 ./ den(validDenominator);
T_C = T_K - 273.15;

end

function fs = prepareFeldsparRow(data_fs, MWinfo)
% prepareFeldsparRow
% Convert one feldspar analysis to anhydrous cation fractions and An-Ab-Or
% fractions. Existing NaN values are retained. Missing optional CaO, Na2O,
% or K2O columns are treated as zero.

SiO2 = getMineralOxideRequired(data_fs, 'SiO2');
Al2O3 = getMineralOxideRequired(data_fs, 'Al2O3');
CaO = getMineralOxideOptional(data_fs, 'CaO', 0);
Na2O = getMineralOxideOptional(data_fs, 'Na2O', 0);
K2O = getMineralOxideOptional(data_fs, 'K2O', 0);

nSi = SiO2 ./ MWinfo.MW.SiO2 .* MWinfo.Cat.SiO2;
nAl = Al2O3 ./ MWinfo.MW.Al2O3 .* MWinfo.Cat.Al2O3;
nCa = CaO ./ MWinfo.MW.CaO .* MWinfo.Cat.CaO;
nNa = Na2O ./ MWinfo.MW.Na2O .* MWinfo.Cat.Na2O;
nK = K2O ./ MWinfo.MW.K2O .* MWinfo.Cat.K2O;

catSum = nSi + nAl + nCa + nNa + nK;
alkSum = nCa + nNa + nK;

% Division by a zero or NaN sum intentionally produces NaN. The calling
% function prints a non-stopping result warning rather than converting the
% value to zero or terminating the calculation.
fs.XSiO2 = nSi ./ catSum;
fs.XAlO1_5 = nAl ./ catSum;
fs.XCaO = nCa ./ catSum;
fs.XNaO0_5 = nNa ./ catSum;
fs.XKO0_5 = nK ./ catSum;

fs.XAn = nCa ./ alkSum;
fs.XAb = nNa ./ alkSum;
fs.XOr = nK ./ alkSum;

end

function liq = prepareLiquidRow(data_liq, MWinfo)
% prepareLiquidRow
% Convert one liquid analysis to anhydrous cation fractions. Missing optional
% columns are treated as zero, while NaN in an existing column is preserved.

SiO2 = getLiquidOxideOptional(data_liq, 'SiO2', 0);
TiO2 = getLiquidOxideOptional(data_liq, 'TiO2', 0);
Al2O3 = getLiquidOxideOptional(data_liq, 'Al2O3', 0);

% Prefer FeO when that column exists. Otherwise use FeOt. If the selected
% existing column contains NaN, retain NaN rather than falling back to zero.
feName = selectLiquidFeColumnName(data_liq.Properties.VariableNames);
if isempty(feName)
    FeO = 0;
else
    FeO = toScalarDouble(data_liq.(feName));
end

MnO = getLiquidOxideOptional(data_liq, 'MnO', 0);
MgO = getLiquidOxideOptional(data_liq, 'MgO', 0);
CaO = getLiquidOxideOptional(data_liq, 'CaO', 0);
Na2O = getLiquidOxideOptional(data_liq, 'Na2O', 0);
K2O = getLiquidOxideOptional(data_liq, 'K2O', 0);
V2O3 = getLiquidOxideOptional(data_liq, 'V2O3', 0);
Cr2O3 = getLiquidOxideOptional(data_liq, 'Cr2O3', 0);
NiO = getLiquidOxideOptional(data_liq, 'NiO', 0);
P2O5 = getLiquidOxideOptional(data_liq, 'P2O5', 0);
SO3 = getLiquidOxideOptional(data_liq, 'SO3', 0);
F = getLiquidOxideOptional(data_liq, 'F', 0);
Cl = getLiquidOxideOptional(data_liq, 'Cl', 0);
Fe2O3 = getLiquidOxideOptional(data_liq, 'Fe2O3', 0);
H2O = getLiquidOxideOptional(data_liq, 'H2O', 0);

nSi = SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
nTi = TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
nAl = Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
nFe = FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
nMn = MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
nMg = MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
nCa = CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
nNa = Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
nK = K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
nV = V2O3 .* MWinfo.Cat.V2O3 ./ MWinfo.MW.V2O3;
nCr = Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
nNi = NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;
nP = P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;
nS = SO3 .* MWinfo.Cat.SO3 ./ MWinfo.MW.SO3;
nF = F .* MWinfo.Cat.F ./ MWinfo.MW.F;
nCl = Cl .* MWinfo.Cat.Cl ./ MWinfo.MW.Cl;
nFe3 = Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

tot = nSi + nTi + nAl + nFe + nMn + nMg + nCa + nNa + nK + ...
      nV + nCr + nNi + nP + nS + nF + nCl + nFe3;

% As for feldspar, a zero or NaN total is allowed to propagate as NaN.
liq.XSiO2 = nSi ./ tot;
liq.XAlO1_5 = nAl ./ tot;
liq.XCaO = nCa ./ tot;
liq.XNaO0_5 = nNa ./ tot;
liq.XKO0_5 = nK ./ tot;
liq.H2O = H2O;

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Attach optional liquid identifiers to every pressure row.

variableNames = data_liq.Properties.VariableNames;
nRows = height(row);

if any(strcmp(variableNames, 'Index'))
    row.liq_Index = repeatScalarValue(data_liq.('Index'), nRows);
end
if any(strcmp(variableNames, 'Experiment'))
    row.liq_Experiment = repmat(string(data_liq.('Experiment')), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    row.liq_Citation = repmat(string(data_liq.('Citation')), nRows, 1);
end

end

function repeatedValue = repeatScalarValue(rawValue, nRows)
% repeatScalarValue
% Repeat the first value of a one-row table variable while preserving its
% basic MATLAB type where practical.

if iscell(rawValue)
    repeatedValue = repmat(rawValue(1), nRows, 1);
elseif ischar(rawValue)
    repeatedValue = repmat(string(rawValue), nRows, 1);
else
    repeatedValue = repmat(rawValue(1), nRows, 1);
end

end

function value = getMineralOxideRequired(data_tbl, oxide)
% getMineralOxideRequired
% Read a required mineral oxide column. The column must exist, but an
% existing NaN value is returned unchanged.

name = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(name)
    error('Selected mineral row must contain variable: %s', oxide);
end

value = toScalarDouble(data_tbl.(name));

end

function value = getMineralOxideOptional(data_tbl, oxide, defaultValue)
% getMineralOxideOptional
% Return defaultValue only when the column is absent. NaN in an existing
% column remains NaN.

name = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(name)
    value = defaultValue;
else
    value = toScalarDouble(data_tbl.(name));
end

end

function value = getLiquidOxideOptional(data_liq, oxide, defaultValue)
% getLiquidOxideOptional
% Return defaultValue only when the column is absent. NaN in an existing
% liquid column remains NaN.

name = findOxideColumn(data_liq.Properties.VariableNames, oxide);
if isempty(name)
    value = defaultValue;
else
    value = toScalarDouble(data_liq.(name));
end

end

function name = selectLiquidFeColumnName(variableNames)
% selectLiquidFeColumnName
% Select FeO when present; otherwise select FeOt. Return an empty character
% vector when neither column exists.

name = findOxideColumn(variableNames, 'FeO');
if isempty(name)
    name = findOxideColumn(variableNames, 'FeOt');
end

end

function name = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match oxide names after removing spaces, underscores, and hyphens. A
% column with the suffix "value" is preferred over the bare oxide name.

canonicalNames = cell(size(variableNames));
for i = 1:numel(variableNames)
    textValue = lower(variableNames{i});
    textValue = strrep(textValue, ' ', '');
    textValue = strrep(textValue, '_', '');
    textValue = strrep(textValue, '-', '');
    canonicalNames{i} = textValue;
end

canonicalOxide = lower(oxide);
canonicalOxide = strrep(canonicalOxide, ' ', '');
canonicalOxide = strrep(canonicalOxide, '_', '');
canonicalOxide = strrep(canonicalOxide, '-', '');

targets = {[canonicalOxide 'value'], canonicalOxide};

name = '';
for i = 1:numel(targets)
    index = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(index)
        name = variableNames{index};
        return;
    end
end

end

function value = toScalarDouble(rawValue)
% toScalarDouble
% Convert the first table value to double. Empty, missing, non-numeric, or
% explicitly NaN content is represented as NaN and is never replaced by zero.

value = NaN;

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
    value = str2double(rawValue(1));
    return;
end

if ischar(rawValue)
    value = str2double(string(rawValue));
    return;
end

if iscell(rawValue)
    if isempty(rawValue{1})
        return;
    end

    firstValue = rawValue{1};
    if isnumeric(firstValue) || islogical(firstValue)
        value = double(firstValue(1));
    elseif isstring(firstValue)
        if ~ismissing(firstValue(1))
            value = str2double(firstValue(1));
        end
    elseif ischar(firstValue)
        value = str2double(string(firstValue));
    end
end

end
