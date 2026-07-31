function results = Putirka2008ol(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Olivine_Liquidus/Putirka2008ol.m
% Tested structurally against the MATLAB R2024b coding style used by
% Ballhaus1991.m
%
% Olivine–Liquid thermometers: Equations 19 and 22
% Putirka, K.D. (2008)
% Reviews in Mineralogy and Geochemistry, 69, 61–120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis with the first row
% of the Liquid dataset selected through liquid.readLiquidExcel(), and
% calculates temperature using Putirka's (2008) Equations 19 and/or 22.
%
% The function accepts pressure as either a scalar or a vector. Therefore,
% it can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for each pressure value
% for every selected Olivine–Liquid pair.
%
% The function is designed for repeated calculations. After each run it asks
% whether another olivine analysis should be calculated using the same
% selected Liquid dataset, and stores all result blocks in one output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Putirka (2008) discusses the olivine–liquid thermometers and their tests on
% pp. 73–78. Equations 19, 21, and 22 are presented on pp. 76–77, and the
% comparison of the thermometers is shown in Figure 4 on p. 78.
%
% The global experimental database used for the tests shown in Figures 2
% and 4 spans approximately:
%
%   Temperature       : 729–2000 degreeC
%   Pressure          : 0.0001–14.4 GPa
%   Liquid SiO2       : 31.5–73.64 wt%
%   Liquid Na2O + K2O : 0–14.3 wt%
%   Liquid H2O        : 0–18.6 wt%
%
% These ranges are summarized for the experimental database on p. 73, with
% the olivine–liquid thermometer tests and data linkage described on
% pp. 76–78. They represent the broad test-data range and are not a strict,
% equation-specific guarantee of accuracy throughout the entire interval.
%
% Important equation-specific cautions stated by Putirka (2008) include:
%
%   Equation 19:
%     - Performs best for anhydrous systems.
%     - Putirka (2008, p. 77) states that no new calibration is required for
%       anhydrous systems at T < 1650 degreeC.
%     - It tends to overestimate temperature for hydrous compositions.
%     - Systematic error may occur at very high temperature and pressure;
%       Putirka (2008, pp. 76–78) discusses the Herzberg and O'Hara (2002)
%       pressure correction used in the Figure 4 evaluation.
%
%   Equation 22:
%     - Is preferred when H2O is present (Putirka, 2008, pp. 77–78).
%     - The test dataset includes H2O contents from 0 to 18.6 wt%.
%
%   Equilibrium requirement:
%     - Olivine and liquid must approach equilibrium; otherwise the
%       calculated temperature has no physical meaning (pp. 73–76).
%     - For common basaltic systems at P < approximately 2–3 GPa,
%       KD(Fe–Mg)ol-liq is commonly near 0.30, but it varies with pressure
%       and liquid composition. An appropriate equilibrium test should be
%       performed before interpreting calculated temperatures.
%
%   Special compositions:
%     - For peridotitic systems containing approximately 2–25 wt% CO2,
%       Putirka (2008, p. 77) recommends considering the Sisson and Grove
%       (1993b) expression because it outperforms Equations 19 and 22 for
%       those special compositions.
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) input pressure is outside 0.0001–14.4 GPa;
%   2) a finite Equation 19 result is outside 729–1650 degreeC;
%   3) a finite Equation 22 result is outside 729–2000 degreeC;
%   4) Equation 19 is used with a finite H2O value greater than zero;
%   5) Equation 22 is used with H2O outside 0–18.6 wt%;
%   6) required inputs contain NaN; or
%   7) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%
% The FIRST column is treated as an identifier ("data code") displayed in
% the selection dialog.
%
% Olivine composition is read in the following priority order:
%   1) Mg_cation_apfu and Fe_cation_apfu
%   2) XMg
%   3) Fo
%   4) MgO and FeO (or FeOt)
%
% Liquid composition is read through liquid.readLiquidExcel(). Component
% normalization uses the following oxide/element columns when available:
%   SiO2 TiO2 Al2O3 FeO MnO MgO CaO Na2O K2O V2O3 Cr2O3 NiO
%   P2O5 SO3 F Cl Fe2O3
%
% Missing optional Liquid columns are treated as zero, preserving the policy
% of the original script. In contrast, an existing cell containing NaN is
% retained as NaN and is never converted to zero. Because all listed
% components enter the normalization sum, a NaN among them propagates into
% the calculated components and temperature. H2O is also retained as NaN
% when the H2O cell itself contains NaN.
%
% Finite inputs that are mathematically required to be positive are checked:
%   - Olivine Mg and Fe values used to calculate XMg must be > 0.
%   - Liquid SiO2, FeO, and MgO must be > 0.
% Other finite additive Liquid components, including H2O, may be zero but
% must not be negative. This distinction preserves valid anhydrous and
% zero-concentration analyses while preventing mathematically invalid input.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Liquid components are calculated as cation fractions on an anhydrous basis
% using molecular weights and cation numbers loaded by
% liquid.getMolarWeights().
%
%   C_SiO2   = X_liq_SiO2
%   C_liq_NM = X_liq_FeO + X_liq_MnO + X_liq_MgO
%              + X_liq_CaO + X_liq_NiO
%   C_liq_NF = X_liq_SiO2 + X_liq_NaO0.5
%              + X_liq_KO0.5 + X_liq_TiO2
%   NF       = (7/2)ln(1 - X_liq_AlO1.5)
%              + 7ln(1 - X_liq_TiO2)
%
%   D_Mg(ol-liq) = XMg_ol / XMg_liq
%
% Equation 19:
%   T(degreeC) =
%     [13603 + 4.943e-7(P(Pa) - 1e-5)] /
%     [6.26 + 2lnD_Mg + 2ln(1.5C_liq_NM)
%      + 2ln(3C_SiO2) - NF] - 273.15
%
% Equation 22:
%   This file retains the iterative Equation 21–22 implementation of the
%   supplied Putirka2008ol.m file for reproducibility. The iteration is
%   performed independently for every pressure point and is never seeded by
%   replacing a NaN input with zero.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008ol(rawdata_struct, P_kbar)
%   results = Putirka2008ol(..., 'Equation', 19)
%   results = Putirka2008ol(..., 'Equation', 22)
%   results = Putirka2008ol(..., 'Equation', 'both')
%
% Inputs:
%   rawdata_struct : struct containing an Olivine table
%   P_kbar         : pressure in kbar; finite non-negative scalar or vector
%
% Options:
%   Equation : 19, 22, or 'both' (default: 'both')
%   MaxIter  : maximum iterations for Equation 22 (default: 100)
%   Tol      : Equation 22 convergence tolerance in degreeC (default: 1e-8)
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Olivine–Liquid pair
%

%% Input validation
% Basic argument checks prevent silent failures and permit either a fixed
% pressure scalar or a pressure vector supplied by startThermoCalc_rangeP.
if nargin < 2
    error('Putirka2008ol requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end
if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end

P_kbar = P_kbar(:);

%% Options
ip = inputParser;
ip.addParameter('Equation', 'both', ...
    @(x) isnumeric(x) || ischar(x) || isstring(x));
ip.addParameter('MaxIter', 100, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && fix(x) == x);
ip.addParameter('Tol', 1e-8, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
ip.parse(varargin{:});

eqOpt = normalizeEquationOption(ip.Results.Equation);
maxIter = double(ip.Results.MaxIter);
tol = double(ip.Results.Tol);

%% 1) Retrieve datasets and calculation constants
disp('=== Step 1: Preparing Olivine and Liquid datasets ===');

dataset_ol = rawdata_struct.Olivine;
MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();

if isempty(liqAll)
    error('Selected Liquid dataset is empty.');
end

% The original interface uses the first row from the selected Liquid dataset.
selectedIdx_liq = 1;
selectedData_liq = liqAll(selectedIdx_liq, :);

if height(liqAll) > 1
    fprintf(2, ...
        ['WARNING: The selected Liquid dataset contains %d rows. ' ...
         'Only the first row is used, following the supplied Putirka2008ol.m behavior.\n'], ...
        height(liqAll));
end

disp('=== Preparing Olivine and Liquid datasets has been finished ===');

%% 2) Initialize output container and warning limits
% Result blocks are buffered in a preallocated cell array and concatenated
% once after the interactive loop. This avoids resizing the results table on
% every iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Broad experimental test-data limits discussed by Putirka (2008).
testT_min_degC = 729;
testT_max_degC = 2000;
eq19RecommendedT_max_degC = 1650;
testP_min_GPa = 0.0001;
testP_max_GPa = 14.4;
testH2O_min_wt = 0;
testH2O_max_wt = 18.6;

P_GPa_input = P_kbar ./ 10;
pressureOutsideTestRange = ...
    P_GPa_input < testP_min_GPa | P_GPa_input > testP_max_GPa;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Olivine) ===');

while true
    % ----- Olivine selection -----
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

    % ----- Liquid selection -----
    % liquid.readLiquidExcel() already handles the dataset choice. This
    % thermometer follows the original script and uses row 1 automatically.
    disp('=== Step 4: Selecting Liquid data ===');
    disp('Liquid selected: (auto) Row 1');

    % ----- Input checks and calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    % NaN is reported but is intentionally allowed to propagate.
    nanInputNames = findNaNInputs(selectedData_ol, selectedData_liq, eqOpt);

    % Stop only for finite values that are mathematically invalid.
    validateCompositionInputs(selectedData_ol, selectedData_liq, eqOpt);

    row = calcTemp(selectedData_ol, selectedData_liq, P_kbar, ...
        MWinfo, eqOpt, maxIter, tol);

    % Replicate identifiers to match one row per pressure value.
    nRows = height(row);
    row.dataCode_ol = repmat(string(selectedCode_ol), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_ol', 'dataRow_liq'}, 'Before', 1);

    % Store this calculation as one buffered table block.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % ----- Immediate result display -----
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if equationRequested(eqOpt, 19)
        if height(row) == 1
            disp([char(string(selectedCode_ol)) ' Eq.19: ' ...
                num2str(row.T19_C) ' degreeC']);
        else
            disp([char(string(selectedCode_ol)) ' Eq.19: ' ...
                num2str(row.T19_C(1)) ' to ' ...
                num2str(row.T19_C(end)) ' degreeC']);
        end
    end

    if equationRequested(eqOpt, 22)
        if height(row) == 1
            disp([char(string(selectedCode_ol)) ' Eq.22: ' ...
                num2str(row.T22_C) ' degreeC']);
        else
            disp([char(string(selectedCode_ol)) ' Eq.22: ' ...
                num2str(row.T22_C(1)) ' to ' ...
                num2str(row.T22_C(end)) ' degreeC']);
        end
    end

    % ----- Pressure warning -----
    % Pressure is common to every selected olivine in this function call, so
    % print this warning only once.
    if any(pressureOutsideTestRange) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the broad experimental test-data range ' ...
             'reported for the Putirka (2008) olivine–liquid thermometer tests: ' ...
             '0.0001–14.4 GPa (0.001–144 kbar; pp. 73–78). ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g GPa.\n'], ...
            sum(pressureOutsideTestRange), numel(P_GPa_input), ...
            min(P_GPa_input), max(P_GPa_input));
        pressureWarningIssued = true;
    end

    % ----- Equation-specific temperature warnings -----
    if equationRequested(eqOpt, 19)
        finiteT19 = isfinite(row.T19_C);
        outsideT19 = finiteT19 & ...
            (row.T19_C < testT_min_degC | ...
             row.T19_C > eq19RecommendedT_max_degC);

        if any(outsideT19)
            finiteValues = row.T19_C(finiteT19);
            fprintf(2, ...
                ['WARNING: Equation 19 calculated temperature is outside the ' ...
                 'recommended/tested interval used here for screening: ' ...
                 '729–1650 degreeC. The lower limit represents the broad test-data ' ...
                 'range, whereas Putirka (2008, p. 77) specifically recommends the ' ...
                 'anhydrous model at T < 1650 degreeC. ' ...
                 '%d of %d finite point(s) are outside; finite range = ' ...
                 '%.4g–%.4g degreeC for %s.\n'], ...
                sum(outsideT19), sum(finiteT19), ...
                min(finiteValues), max(finiteValues), ...
                char(string(selectedCode_ol)));
        end
    end

    if equationRequested(eqOpt, 22)
        finiteT22 = isfinite(row.T22_C);
        outsideT22 = finiteT22 & ...
            (row.T22_C < testT_min_degC | row.T22_C > testT_max_degC);

        if any(outsideT22)
            finiteValues = row.T22_C(finiteT22);
            fprintf(2, ...
                ['WARNING: Equation 22 calculated temperature is outside the broad ' ...
                 'experimental test-data range of 729–2000 degreeC reported for the ' ...
                 'Putirka (2008) tests (pp. 73–78). ' ...
                 '%d of %d finite point(s) are outside; finite range = ' ...
                 '%.4g–%.4g degreeC for %s.\n'], ...
                sum(outsideT22), sum(finiteT22), ...
                min(finiteValues), max(finiteValues), ...
                char(string(selectedCode_ol)));
        end
    end

    % ----- Water-content cautions -----
    H2O_value = getLiquidOxide(selectedData_liq, 'H2O', 0);

    if equationRequested(eqOpt, 19) && isfinite(H2O_value) && H2O_value > 0
        fprintf(2, ...
            ['WARNING: Equation 19 is being applied to a hydrous Liquid composition ' ...
             '(H2O = %.4g wt%%). Putirka (2008, pp. 76–78) reports that Equation 19 ' ...
             'overestimates temperature for hydrous experiments and recommends ' ...
             'Equation 22 when H2O is present.\n'], ...
            H2O_value);
    end

    if equationRequested(eqOpt, 22) && isfinite(H2O_value) && ...
            (H2O_value < testH2O_min_wt || H2O_value > testH2O_max_wt)
        fprintf(2, ...
            ['WARNING: Liquid H2O is outside the 0–18.6 wt%% range represented in ' ...
             'the Putirka (2008) Equation 22 test dataset: H2O = %.4g wt%%.\n'], ...
            H2O_value);
    end

    % ----- NaN-input warning -----
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         NaN was retained as a missing value and was not replaced by zero.\n' ...
             '         The calculation continued, and affected outputs may remain NaN.\n'], ...
            char(string(selectedCode_ol)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % ----- Non-finite-result and convergence warnings -----
    if equationRequested(eqOpt, 19)
        invalidT19 = ~isfinite(row.T19_C);
        if any(invalidT19)
            fprintf(2, ...
                ['WARNING: Non-finite Equation 19 temperature values were calculated ' ...
                 'for %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
                 '         These values remain in the output table; calculation was not stopped.\n'], ...
                char(string(selectedCode_ol)), ...
                sum(invalidT19), numel(row.T19_C), ...
                sum(isnan(row.T19_C)), sum(isinf(row.T19_C)));
        end
    end

    if equationRequested(eqOpt, 22)
        invalidT22 = ~isfinite(row.T22_C);
        if any(invalidT22)
            fprintf(2, ...
                ['WARNING: Non-finite Equation 22 temperature values were calculated ' ...
                 'for %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
                 '         These values remain in the output table; calculation was not stopped.\n'], ...
                char(string(selectedCode_ol)), ...
                sum(invalidT22), numel(row.T22_C), ...
                sum(isnan(row.T22_C)), sum(isinf(row.T22_C)));
        end

        validButNotConverged = isfinite(row.T22_C) & ~row.Eq22_converged;
        if any(validButNotConverged)
            fprintf(2, ...
                ['WARNING: Equation 22 did not satisfy the convergence tolerance ' ...
                 '(Tol = %.4g degreeC) within MaxIter = %d for %d of %d finite point(s) ' ...
                 'for %s. The last finite iteration values were retained.\n'], ...
                tol, maxIter, sum(validButNotConverged), ...
                sum(isfinite(row.T22_C)), char(string(selectedCode_ol)));
        end
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another olivine and the same Liquid row.
    userAction = questdlg( ...
        'Continue with another data selection (same Liquid dataset)?', ...
        'Putirka2008ol', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all result blocks once after the loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'liquid', metaLiq, ...
    'equation', eqOpt, ...
    'sourceFunction', 'Putirka2008ol');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function eqOpt = normalizeEquationOption(value)
% normalizeEquationOption
% Convert the Equation option to one of the stable strings: "19", "22",
% or "both". Invalid values stop before the interactive calculation begins.

if isnumeric(value)
    if ~isscalar(value) || ~isfinite(value) || ~(value == 19 || value == 22)
        error('Equation must be 19, 22, or ''both''.');
    end
    eqOpt = string(value);
    return;
end

valueString = lower(strtrim(string(value)));
if numel(valueString) ~= 1 || ...
        ~(valueString == "19" || valueString == "22" || valueString == "both")
    error('Equation must be 19, 22, or ''both''.');
end

eqOpt = valueString;

end

function tf = equationRequested(eqOpt, equationNumber)
% equationRequested
% Return true when the normalized option requests the specified equation.

tf = eqOpt == "both" || eqOpt == string(equationNumber);

end

function nanInputNames = findNaNInputs(data_olivine, data_liquid, eqOpt)
% findNaNInputs
% Return names of input values that contain NaN. Existing NaN values are
% retained and never replaced with zero.

% Maximum possible entries: 2 olivine variables + 17 normalization inputs
% + H2O. Preallocation avoids growing the string array within the loops.
maxNames = 20;
nameBuffer = strings(maxNames, 1);
nNames = 0;

olivineNames = getOlivineSourceVariableNames(data_olivine);
for i = 1:numel(olivineNames)
    variableName = olivineNames{i};
    variableValue = toScalarDouble(data_olivine.(variableName), NaN);
    if isnan(variableValue)
        nNames = nNames + 1;
        nameBuffer(nNames) = "Olivine." + string(variableName);
    end
end

normalizationOxides = { ...
    'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', ...
    'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', ...
    'F', 'Cl', 'Fe2O3'};

for i = 1:numel(normalizationOxides)
    oxide = normalizationOxides{i};
    columnName = findOxideColumn(data_liquid.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDouble(data_liquid.(columnName), NaN);
        if isnan(value)
            nNames = nNames + 1;
            nameBuffer(nNames) = "Liquid." + string(columnName);
        end
    end
end

if equationRequested(eqOpt, 22)
    columnName = findOxideColumn(data_liquid.Properties.VariableNames, 'H2O');
    if ~isempty(columnName)
        value = toScalarDouble(data_liquid.(columnName), NaN);
        if isnan(value)
            nNames = nNames + 1;
            nameBuffer(nNames) = "Liquid." + string(columnName);
        end
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateCompositionInputs(data_olivine, data_liquid, eqOpt)
% validateCompositionInputs
% Stop for finite values that are mathematically invalid. NaN is excluded
% deliberately so that it propagates and is reported by fprintf warnings.

invalidBuffer = strings(24, 1);
nInvalid = 0;

% ----- Olivine -----
olivineNames = getOlivineSourceVariableNames(data_olivine);

if numel(olivineNames) == 1 && ...
        (strcmp(olivineNames{1}, 'XMg') || strcmp(olivineNames{1}, 'Fo'))
    value = toScalarDouble(data_olivine.(olivineNames{1}), NaN);
    if strcmp(olivineNames{1}, 'Fo') && isfinite(value) && value > 1.5
        value = value ./ 100;
    end
    if isfinite(value) && (value <= 0 || value > 1)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Olivine." + string(olivineNames{1}) + ...
            " (must resolve to 0 < XMg <= 1)";
    end
else
    for i = 1:numel(olivineNames)
        variableName = olivineNames{i};
        value = toScalarDouble(data_olivine.(variableName), NaN);
        if isfinite(value) && value <= 0
            nInvalid = nInvalid + 1;
            invalidBuffer(nInvalid) = "Olivine." + string(variableName) + ...
                " (must be > 0)";
        end
    end
end

% ----- Liquid values that must be strictly positive -----
positiveOxides = {'SiO2', 'FeO', 'MgO'};
for i = 1:numel(positiveOxides)
    oxide = positiveOxides{i};
    columnName = findOxideColumn(data_liquid.Properties.VariableNames, oxide);
    if isempty(columnName)
        error('Selected Liquid row lacks the required composition field: %s.', oxide);
    end
    value = toScalarDouble(data_liquid.(columnName), NaN);
    if isfinite(value) && value <= 0
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = "Liquid." + string(columnName) + ...
            " (must be > 0)";
    end
end

% ----- Other additive Liquid values may be zero but not negative -----
nonNegativeOxides = { ...
    'TiO2', 'Al2O3', 'MnO', 'CaO', 'Na2O', 'K2O', 'V2O3', ...
    'Cr2O3', 'NiO', 'P2O5', 'SO3', 'F', 'Cl', 'Fe2O3'};

if equationRequested(eqOpt, 22)
    nonNegativeOxides{end + 1} = 'H2O';
end

for i = 1:numel(nonNegativeOxides)
    oxide = nonNegativeOxides{i};
    columnName = findOxideColumn(data_liquid.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDouble(data_liquid.(columnName), NaN);
        if isfinite(value) && value < 0
            nInvalid = nInvalid + 1;
            invalidBuffer(nInvalid) = "Liquid." + string(columnName) + ...
                " (must be >= 0)";
        end
    end
end

if nInvalid > 0
    invalidNames = invalidBuffer(1:nInvalid);
    error(['Putirka2008ol: mathematically invalid finite composition ' ...
           'value(s) were found: ' char(strjoin(invalidNames, ', ')) '.']);
end

end

function names = getOlivineSourceVariableNames(data_olivine)
% getOlivineSourceVariableNames
% Return the exact olivine variables used by getOlivineXMg, following the
% same priority order.

variableNames = data_olivine.Properties.VariableNames;

if any(strcmp(variableNames, 'Mg_cation_apfu')) && ...
        any(strcmp(variableNames, 'Fe_cation_apfu'))
    names = {'Mg_cation_apfu', 'Fe_cation_apfu'};
    return;
end

if any(strcmp(variableNames, 'XMg'))
    names = {'XMg'};
    return;
end

if any(strcmp(variableNames, 'Fo'))
    names = {'Fo'};
    return;
end

hasMgO = any(strcmp(variableNames, 'MgO'));
hasFeO = any(strcmp(variableNames, 'FeO'));
hasFeOt = any(strcmp(variableNames, 'FeOt'));

if hasMgO && hasFeO
    names = {'MgO', 'FeO'};
    return;
end

if hasMgO && hasFeOt
    names = {'MgO', 'FeOt'};
    return;
end

error(['Selected Olivine row lacks required composition fields. Provide ' ...
       'Mg_cation_apfu + Fe_cation_apfu, XMg, Fo, or MgO + FeO/FeOt.']);

end

function row = calcTemp(data_olivine, data_liquid, P_kbar, ...
        MWinfo, eqOpt, maxIter, tol)
% calcTemp
% Calculate Equation 19 and/or Equation 22 for one Olivine–Liquid pair over
% a scalar or vector of pressures. Composition values are replicated to the
% pressure-vector length so that every table variable has a stable size.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% ----- Olivine composition -----
XMg_ol_scalar = getOlivineXMg(data_olivine, MWinfo);
XMg_ol = repmat(XMg_ol_scalar, nP, 1);
row.XMg_ol = XMg_ol;

% ----- Liquid oxide values -----
% Existing NaN values remain NaN. Only completely absent optional columns
% use the supplied default value of zero.
SiO2  = getLiquidOxide(data_liquid, 'SiO2',  NaN);
TiO2  = getLiquidOxide(data_liquid, 'TiO2',  0);
Al2O3 = getLiquidOxide(data_liquid, 'Al2O3', 0);
FeO   = getLiquidOxide(data_liquid, 'FeO',   NaN);
MnO   = getLiquidOxide(data_liquid, 'MnO',   0);
MgO   = getLiquidOxide(data_liquid, 'MgO',   NaN);
CaO   = getLiquidOxide(data_liquid, 'CaO',   0);
Na2O  = getLiquidOxide(data_liquid, 'Na2O',  0);
K2O   = getLiquidOxide(data_liquid, 'K2O',   0);
V2O3  = getLiquidOxide(data_liquid, 'V2O3',  0);
Cr2O3 = getLiquidOxide(data_liquid, 'Cr2O3', 0);
NiO   = getLiquidOxide(data_liquid, 'NiO',   0);
P2O5  = getLiquidOxide(data_liquid, 'P2O5',  0);
SO3   = getLiquidOxide(data_liquid, 'SO3',   0);
F     = getLiquidOxide(data_liquid, 'F',     0);
Cl    = getLiquidOxide(data_liquid, 'Cl',    0);
Fe2O3 = getLiquidOxide(data_liquid, 'Fe2O3', 0);
H2O   = getLiquidOxide(data_liquid, 'H2O',   0);

% Store liquid wt% inputs with one value per pressure point.
row.SiO2_liq = repmat(SiO2, nP, 1);
row.TiO2_liq = repmat(TiO2, nP, 1);
row.Al2O3_liq = repmat(Al2O3, nP, 1);
row.FeO_liq = repmat(FeO, nP, 1);
row.MnO_liq = repmat(MnO, nP, 1);
row.MgO_liq = repmat(MgO, nP, 1);
row.CaO_liq = repmat(CaO, nP, 1);
row.Na2O_liq = repmat(Na2O, nP, 1);
row.K2O_liq = repmat(K2O, nP, 1);
row.NiO_liq = repmat(NiO, nP, 1);
row.Cr2O3_liq = repmat(Cr2O3, nP, 1);
row.Fe2O3_liq = repmat(Fe2O3, nP, 1);
row.H2O_liq = repmat(H2O, nP, 1);

% ----- Cation proportions -----
n = struct();
n.SiO2  = SiO2  .* MWinfo.Cat.SiO2  ./ MWinfo.MW.SiO2;
n.TiO2  = TiO2  .* MWinfo.Cat.TiO2  ./ MWinfo.MW.TiO2;
n.Al2O3 = Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n.FeO   = FeO   .* MWinfo.Cat.FeO   ./ MWinfo.MW.FeO;
n.MnO   = MnO   .* MWinfo.Cat.MnO   ./ MWinfo.MW.MnO;
n.MgO   = MgO   .* MWinfo.Cat.MgO   ./ MWinfo.MW.MgO;
n.CaO   = CaO   .* MWinfo.Cat.CaO   ./ MWinfo.MW.CaO;
n.Na2O  = Na2O  .* MWinfo.Cat.Na2O  ./ MWinfo.MW.Na2O;
n.K2O   = K2O   .* MWinfo.Cat.K2O   ./ MWinfo.MW.K2O;
n.V2O3  = V2O3  .* MWinfo.Cat.V2O3  ./ MWinfo.MW.V2O3;
n.Cr2O3 = Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
n.NiO   = NiO   .* MWinfo.Cat.NiO   ./ MWinfo.MW.NiO;
n.P2O5  = P2O5  .* MWinfo.Cat.P2O5  ./ MWinfo.MW.P2O5;
n.SO3   = SO3   .* MWinfo.Cat.SO3   ./ MWinfo.MW.SO3;
n.F     = F     .* MWinfo.Cat.F     ./ MWinfo.MW.F;
n.Cl    = Cl    .* MWinfo.Cat.Cl    ./ MWinfo.MW.Cl;
n.Fe2O3 = Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

totalCations = n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + ...
    n.MgO + n.CaO + n.Na2O + n.K2O + n.V2O3 + n.Cr2O3 + ...
    n.NiO + n.P2O5 + n.SO3 + n.F + n.Cl + n.Fe2O3;

% Initialize all intermediate scalars as NaN. They remain NaN when any
% normalization input is NaN or the total is otherwise invalid.
X_SiO2 = NaN;
X_TiO2 = NaN;
X_AlO1_5 = NaN;
X_FeO = NaN;
X_MnO = NaN;
X_MgO = NaN;
X_CaO = NaN;
X_NiO = NaN;
X_NaO0_5 = NaN;
X_KO0_5 = NaN;
C_SiO2 = NaN;
C_liq_NM = NaN;
C_liq_NF = NaN;
NF = NaN;
XMg_liq_scalar = NaN;
DMg_scalar = NaN;

if isfinite(totalCations) && totalCations > 0
    X_SiO2   = n.SiO2  ./ totalCations;
    X_TiO2   = n.TiO2  ./ totalCations;
    X_AlO1_5 = n.Al2O3 ./ totalCations;
    X_FeO    = n.FeO   ./ totalCations;
    X_MnO    = n.MnO   ./ totalCations;
    X_MgO    = n.MgO   ./ totalCations;
    X_CaO    = n.CaO   ./ totalCations;
    X_NiO    = n.NiO   ./ totalCations;
    X_NaO0_5 = n.Na2O  ./ totalCations;
    X_KO0_5  = n.K2O   ./ totalCations;

    C_SiO2 = X_SiO2;
    C_liq_NM = X_FeO + X_MnO + X_MgO + X_CaO + X_NiO;
    C_liq_NF = X_SiO2 + X_NaO0_5 + X_KO0_5 + X_TiO2;

    if isfinite(1 - X_AlO1_5) && (1 - X_AlO1_5) > 0 && ...
            isfinite(1 - X_TiO2) && (1 - X_TiO2) > 0
        NF = (7 / 2) .* log(1 - X_AlO1_5) + ...
            7 .* log(1 - X_TiO2);
    end

    denominatorMg = n.MgO + n.FeO;
    if isfinite(denominatorMg) && denominatorMg > 0
        XMg_liq_scalar = n.MgO ./ denominatorMg;
        if isfinite(XMg_ol_scalar) && isfinite(XMg_liq_scalar) && ...
                XMg_liq_scalar > 0
            DMg_scalar = XMg_ol_scalar ./ XMg_liq_scalar;
        end
    end
end

% Replicate intermediate values so table height always equals numel(P_kbar).
row.C_SiO2 = repmat(C_SiO2, nP, 1);
row.C_liq_NM = repmat(C_liq_NM, nP, 1);
row.C_liq_NF = repmat(C_liq_NF, nP, 1);
row.NF = repmat(NF, nP, 1);
row.XMg_liq = repmat(XMg_liq_scalar, nP, 1);
row.DMg_ol_liq = repmat(DMg_scalar, nP, 1);

% ----- Equation 19 -----
commonLogInputsValid = ...
    isfinite(DMg_scalar) && DMg_scalar > 0 && ...
    isfinite(C_liq_NM) && C_liq_NM > 0 && ...
    isfinite(C_SiO2) && C_SiO2 > 0 && ...
    isfinite(NF);

T19_C = NaN(nP, 1);
if commonLogInputsValid
    denominator19 = 6.26 ...
        + 2 .* log(DMg_scalar) ...
        + 2 .* log(1.5 .* C_liq_NM) ...
        + 2 .* log(3 .* C_SiO2) ...
        - NF;

    T19_C = (13603 + 4.943e-7 .* (P_GPa .* 1e9 - 1e-5)) ./ ...
        denominator19 - 273.15;
end

% ----- Equation 22 iterative calculation -----
T22_C = NaN(nP, 1);
Eq22_converged = false(nP, 1);
Eq22_iterations = zeros(nP, 1);

validEq22Input = commonLogInputsValid && ...
    isfinite(H2O) && isfinite(Na2O) && isfinite(K2O);

if validEq22Input
    % A finite initial estimate is used only when all Equation 22 inputs are
    % finite. NaN inputs never receive a numerical replacement.
    Twork = repmat(1200, nP, 1);
    failed = false(nP, 1);

    for iteration = 1:maxIter
        active = ~Eq22_converged & ~failed;
        if ~any(active)
            break;
        end

        activeIndex = find(active);
        Tactive = Twork(activeIndex);
        Pactive = P_GPa(activeIndex);

        lnD = -2.158 ...
            + 55.09 .* (Pactive ./ Tactive) ...
            - 6.213e-2 .* H2O ...
            + 4430 ./ Tactive ...
            + 5.115e-2 .* (Na2O + K2O);

        denominator22 = 8.048 ...
            + 2.8352 .* lnD ...
            + 2.097 .* log(1.5 .* C_liq_NM) ...
            + 2.575 .* log(3 .* C_SiO2) ...
            - 1.41 .* NF ...
            + 0.222 .* H2O ...
            + 0.5 .* Pactive;

        Tnew = (15294.6 + 1318.8 .* Pactive ...
            + 2.4834 .* (Pactive .^ 2)) ./ denominator22;

        finiteNew = isfinite(Tnew);
        failedIndex = activeIndex(~finiteNew);
        failed(failedIndex) = true;
        T22_C(failedIndex) = NaN;
        Eq22_iterations(failedIndex) = iteration;

        goodIndex = activeIndex(finiteNew);
        if ~isempty(goodIndex)
            goodNew = Tnew(finiteNew);
            difference = abs(goodNew - Twork(goodIndex));

            Twork(goodIndex) = goodNew;
            T22_C(goodIndex) = goodNew;
            Eq22_iterations(goodIndex) = iteration;

            newlyConverged = difference < tol;
            Eq22_converged(goodIndex(newlyConverged)) = true;
        end
    end
end

% ----- Pack outputs according to the selected equation option -----
if equationRequested(eqOpt, 19)
    row.T19_C = T19_C;
end

if equationRequested(eqOpt, 22)
    row.T22_C = T22_C;
    row.Eq22_converged = Eq22_converged;
    row.Eq22_iterations = Eq22_iterations;
end

end

function row = attachLiquidIDs(row, data_liquid)
% attachLiquidIDs
% Replicate available Liquid identifiers so their lengths match the pressure
% vector and attach them to the result table.

nRows = height(row);
variableNames = data_liquid.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    indexValues = data_liquid.('Index');
    row.liq_Index = repmat(indexValues(1), nRows, 1);
end
if any(strcmp(variableNames, 'Experiment'))
    experimentValues = data_liquid.('Experiment');
    row.liq_Experiment = repmat(string(experimentValues(1)), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    citationValues = data_liquid.('Citation');
    row.liq_Citation = repmat(string(citationValues(1)), nRows, 1);
end

end

function XMg_ol = getOlivineXMg(data_olivine, MWinfo)
% getOlivineXMg
% Calculate olivine XMg using the documented priority order. Existing NaN
% values are retained and propagate through the calculation.

variableNames = data_olivine.Properties.VariableNames;

if any(strcmp(variableNames, 'Mg_cation_apfu')) && ...
        any(strcmp(variableNames, 'Fe_cation_apfu'))
    Mg = toScalarDouble(data_olivine.Mg_cation_apfu, NaN);
    Fe = toScalarDouble(data_olivine.Fe_cation_apfu, NaN);
    XMg_ol = Mg ./ (Mg + Fe);
    return;
end

if any(strcmp(variableNames, 'XMg'))
    XMg_ol = toScalarDouble(data_olivine.XMg, NaN);
    return;
end

if any(strcmp(variableNames, 'Fo'))
    Fo = toScalarDouble(data_olivine.Fo, NaN);
    if isfinite(Fo) && Fo > 1.5
        Fo = Fo ./ 100;
    end
    XMg_ol = Fo;
    return;
end

MgO = toScalarDouble(data_olivine.MgO, NaN);
if any(strcmp(variableNames, 'FeO'))
    FeO = toScalarDouble(data_olivine.FeO, NaN);
else
    FeO = toScalarDouble(data_olivine.FeOt, NaN);
end

nMg = MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
nFe = FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
XMg_ol = nMg ./ (nMg + nFe);

end

function value = getLiquidOxide(data_liquid, oxide, missingDefault)
% getLiquidOxide
% Read one scalar Liquid value. A completely absent optional column returns
% missingDefault. An existing NaN or missing cell returns NaN and is never
% converted to missingDefault.

columnName = findOxideColumn(data_liquid.Properties.VariableNames, oxide);

if isempty(columnName)
    value = missingDefault;
    return;
end

value = toScalarDouble(data_liquid.(columnName), NaN);

end

function name = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Robustly match names such as SiO2value, SiO2 value, SiO2_value, or SiO2.

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

function value = toScalarDouble(rawValue, missingDefault)
% toScalarDouble
% Convert the first table value to a scalar double. Existing NaN, missing,
% empty strings, and unconvertible values return missingDefault. Callers use
% NaN as missingDefault when an existing missing value must remain NaN.

value = missingDefault;

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
    converted = str2double(rawValue(1));
    if ~isnan(converted)
        value = converted;
    end
    return;
end

if ischar(rawValue)
    converted = str2double(string(rawValue));
    if ~isnan(converted)
        value = converted;
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
    if isstring(firstValue)
        if ismissing(firstValue(1))
            return;
        end
        converted = str2double(firstValue(1));
        if ~isnan(converted)
            value = converted;
        end
        return;
    end
    if ischar(firstValue)
        converted = str2double(string(firstValue));
        if ~isnan(converted)
            value = converted;
        end
        return;
    end
end

end
