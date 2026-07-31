function results = Putirka2005(rawdata_struct, P_kbar, varargin)
% Putirka2005.m
% Tested with MATLAB R2024b
%
% Plagioclase-Liquid thermometer, Table 2, Model A
% Putirka, K.D. (2005)
% American Mineralogist, 90, 336–346
% DOI: https://doi.org/10.2138/am.2005.1449
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Plagioclase analysis from
% rawdata_struct.Plagioclase and combines it with one row from a liquid
% dataset loaded by liquid.readLiquidExcel(). Temperature is calculated
% using Putirka (2005), Table 2, Model A.
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Plagioclase-Liquid pair, one output row is returned for each
% pressure value. It can therefore be called by both startThermoCalc_fixedP
% and startThermoCalc_rangeP.
%
% The function is designed for repeated calculations: after each run it asks
% whether another Plagioclase analysis should be calculated using the same
% loaded liquid dataset, and stores all result blocks in one output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Putirka (2005) does not list a strict, Model-A-only minimum and maximum
% calibration range separately from the experimental database. The new
% calibrations were selected from the experimental dataset used in the
% paper, which spans approximately:
%
%   Temperature : 850–1430 degreeC
%   Pressure    : 0.001–27 kbar
%   Liquid SiO2 : 42–73 wt%
%   Conditions  : anhydrous and hydrous, plagioclase-saturated experiments
%
% The full experimental-database ranges are stated on p. 337, and the
% individual experimental studies are listed in Table 1 on p. 338. The new
% calibration strategy and selection of calibration/test data are described
% on p. 341. These limits are used here as practical non-stopping warning
% limits, but they should not be interpreted as sharply defined universal
% boundaries for every composition represented by Model A.
%
% Important application cautions:
%   1) Use a texturally and compositionally plausible equilibrium pair of
%      plagioclase and melt. The calibration is based on coexisting
%      experimental plagioclase and glass compositions (p. 341).
%
%   2) Na loss during 1-atm experiments or microbeam analysis of hydrous
%      glasses can bias liquid compositions and thermometer results. This
%      issue is discussed on p. 337. Natural glass analyses should therefore
%      be screened for possible alkali loss.
%
%   3) H2O is entered independently in wt% and is excluded from the
%      anhydrous liquid cation-fraction calculation (Table 2, p. 342).
%      Uncertainty in H2O directly affects Model A temperatures. When an H2O
%      column is absent, this implementation assumes H2O = 0 wt%; when an H2O
%      value is present but NaN, the NaN is retained and propagated.
%
%   4) Putirka (2005) warns that Model A should not be solved simultaneously
%      with a plagioclase-liquid hygrometer because T and H2O are strongly
%      cross-correlated (discussion on p. 343; Fig. 6 on p. 345). Model A is
%      best used with independently measured or estimated liquid H2O.
%
%   5) Reported Model A standard errors of estimate are approximately:
%        Calibration data  : 24 degreeC (anhydrous), 19 degreeC (hydrous)
%        Independent tests : 24 degreeC (anhydrous), 33 degreeC (hydrous)
%      The overall Model A performance is shown in Fig. 4a on p. 344 and
%      discussed on p. 343. Analytical, pressure, H2O, and disequilibrium
%      uncertainties are additional to these model errors.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 0.001–27 kbar,
%   2) a finite calculated temperature is outside 850–1430 degreeC,
%   3) finite liquid SiO2 is outside 42–73 wt%,
%   4) a calculation input contains NaN, or
%   5) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Plagioclase : table
%
% The FIRST column of the Plagioclase table is treated as an identifier
% ("data code") displayed in the selection dialog.
%
% Required Plagioclase oxide columns:
%   SiO2, Al2O3, CaO, Na2O
%
% Optional Plagioclase oxide columns (0 is used only when the column is
% absent; an explicitly stored NaN remains NaN):
%   TiO2, FeO or FeOt, MnO, MgO, K2O, Cr2O3
%
% The liquid table is loaded by liquid.readLiquidExcel(). The following
% liquid oxides are used to calculate the anhydrous cation fractions exactly
% as specified by Putirka (2005, pp. 341–342):
%   SiO2, TiO2, Al2O3, FeO, MnO, MgO, CaO, Na2O, K2O, Cr2O3
%
% SiO2, Al2O3, CaO, and Na2O columns are required in the liquid table.
% TiO2, FeO/FeOt, MnO, MgO, K2O, and Cr2O3 are treated as 0 only when the
% corresponding column is absent. H2O/H2Ot is used separately in wt% and is
% treated as 0 only when its column is absent.
%
% All finite Plagioclase and liquid oxide values used by the calculation
% must be greater than or equal to zero. Negative finite values stop the
% calculation. NaN values are retained, propagated through the calculation,
% and reported by non-stopping fprintf warnings.
%
% -------------------------------------------------------------------------
% COMPONENT CALCULATIONS
%
% Putirka (2005) states that liquid oxide wt% values are not renormalized
% before calculation of cation fractions (p. 341). The liquid denominator in
% this implementation therefore contains only the published components:
%
%   SiO2, TiO2, AlO1.5, FeO, MnO, MgO, CaO,
%   NaO0.5, KO0.5, and CrO1.5
%
% H2O is excluded from this denominator and is added separately in wt%.
% Additional reported liquid components such as Fe2O3, P2O5, SO3, F, Cl,
% V2O3, or NiO are not included in the Model A cation-fraction denominator.
%
% Plagioclase components are calculated as:
%   An_pl = Ca / (Ca + Na + K)
%   Ab_pl = Na / (Ca + Na + K)
%   Or_pl = K  / (Ca + Na + K)
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Putirka, 2005, Table 2, Model A; p. 342)
%
%   T(K) = 10^4 / D
%
%   D = 6.12
%       + 0.257 * ln[An_pl / (Ca_liq * Al_liq^2 * Si_liq^2)]
%       - 3.166 * Ca_liq
%       + 0.2166 * H2O_liq
%       - 3.137 * [Al_liq / (Al_liq + Si_liq)]
%       + 1.216 * Ab_pl^2
%       - 2.475e-2 * P_kbar
%
% where Ca_liq, Al_liq, and Si_liq are anhydrous liquid cation
% fractions, H2O_liq is in wt%, pressure is in kbar, and T is in Kelvin.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2005(rawdata_struct, P_kbar)
%   results = Putirka2005(rawdata_struct, P_kbar, 'LiquidRow', n)
%
% Inputs:
%   rawdata_struct : struct containing a Plagioclase table
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Name-value option:
%   LiquidRow      : positive integer row number in the loaded liquid table.
%                    If empty, row 1 is used.
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Plagioclase-Liquid pair. T_K and T_degreeC are included for
%             compatibility with ThermoCalc launchers; the legacy Model A
%             variable names are also retained.
%

%% Input validation
% Basic argument checks prevent silent failures while allowing scalar and
% vector pressure inputs from both fixed-P and range-P launchers.
if nargin < 2
    error('Putirka2005 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

if ~isfield(rawdata_struct, 'Plagioclase') || ...
        ~istable(rawdata_struct.Plagioclase)
    error('rawdata_struct must contain table: rawdata_struct.Plagioclase');
end

%% Options
ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve Plagioclase and liquid datasets
% Molecular weights and cation numbers are loaded through the existing
% liquid utilities. The liquid file is selected once and reused during the
% interactive Plagioclase-selection loop.
disp('=== Step 1: Preparing Plagioclase and liquid datasets ===');

MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();

if isempty(liqAll) || ~istable(liqAll)
    error('Selected liquid dataset is empty or is not a table.');
end

dataset_pl = rawdata_struct.Plagioclase;

% Verify required oxide columns before the interactive loop. Values are not
% required to be finite here because NaN values must remain available for
% propagation and warning output.
validateRequiredColumns(dataset_pl, ...
    {{'SiO2'}, {'Al2O3'}, {'CaO'}, {'Na2O'}}, 'Plagioclase');
validateRequiredColumns(liqAll, ...
    {{'SiO2'}, {'Al2O3'}, {'CaO'}, {'Na2O'}}, 'Liquid');

disp('=== Preparing Plagioclase and liquid datasets has been finished ===');

%% 2) Select the liquid row and initialize the output container
% The selected liquid row is fixed for all Plagioclase selections in this
% function call, matching the behavior of the original implementation.
disp('=== Step 2: Preparing liquid selection and output container ===');

if isempty(liquidRowOpt)
    idxLiq = 1;
    if height(liqAll) > 1
        fprintf(2, ...
            ['WARNING: The selected liquid dataset contains %d rows. ' ...
             'Liquid row 1 is being used. Specify ''LiquidRow'' to select ' ...
             'a different row.\n'], ...
            height(liqAll));
    end
else
    idxLiq = liquidRowOpt;
    if idxLiq > height(liqAll)
        error(['Requested LiquidRow (%d) exceeds the number of rows in the ' ...
               'selected liquid dataset (%d).'], idxLiq, height(liqAll));
    end
end

selectedData_liq = liqAll(idxLiq, :);

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Practical warning limits based on the full experimental database used by
% Putirka (2005), rather than strict Model-A-only boundaries.
calibrationT_min_degC = 850;
calibrationT_max_degC = 1430;
calibrationP_min_kbar = 0.001;
calibrationP_max_kbar = 27;
calibrationSiO2_min_wt = 42;
calibrationSiO2_max_wt = 73;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

% The liquid composition is common to all selected Plagioclase analyses, so
% its compositional-range warning needs to be printed only once.
SiO2_liq_forRange = getLiqOxRequired(selectedData_liq, {'SiO2'});
liquidCompositionWarningIssued = false;

disp(['Liquid selected: Row ' num2str(idxLiq)]);
disp('=== Preparing liquid selection and output container has been finished ===');

%% 3–4) Interactive Plagioclase selection and temperature calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a completed calculation.
disp('=== Step 3: Selecting a data code from the list (Plagioclase) ===');

while true
    % ----- Plagioclase selection -----
    % Assumption: the first column contains an identifier shown to the user.
    dataCodes_pl = dataset_pl{:, 1};

    [selectedIdx_pl, ok] = listdlg( ...
        'PromptString', 'Please select the Plagioclase data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_pl, ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_pl)
        disp('Selection canceled');
        break;
    end

    selectedCode_pl = dataCodes_pl(selectedIdx_pl);
    selectedData_pl = dataset_pl(selectedIdx_pl, :);

    disp(['Plagioclase selected: ' char(string(selectedCode_pl))]);
    disp('=== Step 4: Calculating the temperature ===');

    % Check only actual calculation inputs. Missing optional columns are
    % interpreted as zero, but an explicitly stored NaN is retained.
    nanInputNames = findNaNInputs(selectedData_pl, selectedData_liq);

    % Negative finite oxide values are prohibited. NaN is intentionally
    % allowed and propagated; zero is allowed and may generate a non-finite
    % mathematical result that is retained and reported below.
    validateNonNegativeInputs(selectedData_pl, selectedData_liq);

    row = calcTemp(selectedData_pl, selectedData_liq, P_kbar, MWinfo);

    % Repeat identifiers and liquid metadata to match the number of pressure
    % points returned by calcTemp.
    nRows = height(row);
    row.dataCode_pl = repmat(string(selectedCode_pl), nRows, 1);
    row.dataRow_liq = repmat(idxLiq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_pl', 'dataRow_liq'}, 'Before', 1);

    % Store the table block in a preallocated cell buffer. Capacity is
    % doubled only when necessary, avoiding repeated growth of the complete
    % output table during every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature immediately.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_pl)) ' & Liquid row ' num2str(idxLiq) ...
            ': ' num2str(row.T_degreeC) ' degreeC']);
    else
        disp([char(string(selectedCode_pl)) ' & Liquid row ' num2str(idxLiq) ...
            ': ' num2str(row.T_degreeC(1)) ' to ' ...
            num2str(row.T_degreeC(end)) ' degreeC']);
    end

    % Warn once when any pressure input lies outside the experimental
    % database range used as the practical warning range.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental-database ' ...
             'range used by Putirka (2005): 0.001–27 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn once if the selected liquid SiO2 composition is outside the
    % experimental database range. NaN is handled by the NaN-input warning.
    if ~liquidCompositionWarningIssued && isfinite(SiO2_liq_forRange) && ...
            (SiO2_liq_forRange < calibrationSiO2_min_wt || ...
             SiO2_liq_forRange > calibrationSiO2_max_wt)
        fprintf(2, ...
            ['WARNING: Liquid SiO2 is outside the experimental-database ' ...
             'range reported by Putirka (2005): 42–73 wt%%. ' ...
             'Selected liquid SiO2 = %.4g wt%% (Liquid row %d).\n'], ...
            SiO2_liq_forRange, idxLiq);
        liquidCompositionWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the practical
    % temperature range. NaN and Inf are handled separately below.
    finiteTemperature = isfinite(row.T_degreeC);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_degreeC < calibrationT_min_degC | ...
         row.T_degreeC > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_degreeC(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the experimental-' ...
             'database range used by Putirka (2005): 850–1430 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s and ' ...
             'Liquid row %d.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_pl)), ...
            idxLiq);
    end

    % Print the names of all explicitly NaN calculation inputs. The values
    % are not replaced by zero and remain in the calculation.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s ' ...
             'and Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; the ' ...
             'calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_pl)), ...
            idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report NaN/Inf results instead of replacing or deleting them.
    invalidTemperature = ~isfinite(row.T_degreeC);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s ' ...
             'and Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_pl)), ...
            idxLiq, ...
            sum(invalidTemperature), ...
            numel(row.T_degreeC), ...
            sum(isnan(row.T_degreeC)), ...
            sum(isinf(row.T_degreeC)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Plagioclase selection (same Liquid dataset)?', ...
        'Putirka2005', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered result blocks once. Return an empty table when no
% calculation was completed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct('liquid', metaLiq);
disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_plagioclase, data_liquid)
% findNaNInputs
% Return names of explicitly stored calculation inputs whose selected value
% is NaN. Missing optional columns are not listed because they are treated as
% zero by definition in this implementation.

plagAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'Cr2O3'}};
plagLabels = [ ...
    "Plagioclase.SiO2", "Plagioclase.TiO2", ...
    "Plagioclase.Al2O3", "Plagioclase.FeO/FeOt", ...
    "Plagioclase.MnO", "Plagioclase.MgO", ...
    "Plagioclase.CaO", "Plagioclase.Na2O", ...
    "Plagioclase.K2O", "Plagioclase.Cr2O3"];

liqAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'Cr2O3'}, ...
    {'H2O', 'H2Ot'}};
liqLabels = [ ...
    "Liquid.SiO2", "Liquid.TiO2", "Liquid.Al2O3", ...
    "Liquid.FeO/FeOt", "Liquid.MnO", "Liquid.MgO", ...
    "Liquid.CaO", "Liquid.Na2O", "Liquid.K2O", ...
    "Liquid.Cr2O3", "Liquid.H2O/H2Ot"];

maxNames = numel(plagAliases) + numel(liqAliases);
nanNamesBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(plagAliases)
    [value, ~, columnFound] = getOxideValue( ...
        data_plagioclase, plagAliases{i}, 0, false, 'Plagioclase');
    if columnFound && isnan(value)
        nNames = nNames + 1;
        nanNamesBuffer(nNames) = plagLabels(i);
    end
end

for i = 1:numel(liqAliases)
    [value, ~, columnFound] = getOxideValue( ...
        data_liquid, liqAliases{i}, 0, false, 'Liquid');
    if columnFound && isnan(value)
        nNames = nNames + 1;
        nanNamesBuffer(nNames) = liqLabels(i);
    end
end

nanInputNames = nanNamesBuffer(1:nNames);

end

function validateNonNegativeInputs(data_plagioclase, data_liquid)
% validateNonNegativeInputs
% Stop when an explicitly stored calculation input is negative or infinite.
% Zero is allowed. NaN is intentionally allowed and propagated.

plagAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'Cr2O3'}};
plagLabels = [ ...
    "Plagioclase.SiO2", "Plagioclase.TiO2", ...
    "Plagioclase.Al2O3", "Plagioclase.FeO/FeOt", ...
    "Plagioclase.MnO", "Plagioclase.MgO", ...
    "Plagioclase.CaO", "Plagioclase.Na2O", ...
    "Plagioclase.K2O", "Plagioclase.Cr2O3"];

liqAliases = { ...
    {'SiO2'}, {'TiO2'}, {'Al2O3'}, {'FeO', 'FeOt'}, {'MnO'}, ...
    {'MgO'}, {'CaO'}, {'Na2O'}, {'K2O'}, {'Cr2O3'}, ...
    {'H2O', 'H2Ot'}};
liqLabels = [ ...
    "Liquid.SiO2", "Liquid.TiO2", "Liquid.Al2O3", ...
    "Liquid.FeO/FeOt", "Liquid.MnO", "Liquid.MgO", ...
    "Liquid.CaO", "Liquid.Na2O", "Liquid.K2O", ...
    "Liquid.Cr2O3", "Liquid.H2O/H2Ot"];

maxNames = numel(plagAliases) + numel(liqAliases);
invalidNamesBuffer = strings(maxNames, 1);
nInvalid = 0;

for i = 1:numel(plagAliases)
    [value, ~, columnFound] = getOxideValue( ...
        data_plagioclase, plagAliases{i}, 0, false, 'Plagioclase');
    if columnFound && (isinf(value) || (isfinite(value) && value < 0))
        nInvalid = nInvalid + 1;
        invalidNamesBuffer(nInvalid) = plagLabels(i);
    end
end

for i = 1:numel(liqAliases)
    [value, ~, columnFound] = getOxideValue( ...
        data_liquid, liqAliases{i}, 0, false, 'Liquid');
    if columnFound && (isinf(value) || (isfinite(value) && value < 0))
        nInvalid = nInvalid + 1;
        invalidNamesBuffer(nInvalid) = liqLabels(i);
    end
end

if nInvalid > 0
    invalidInputNames = invalidNamesBuffer(1:nInvalid);
    error(['Putirka2005: calculation inputs must not be negative or ' ...
           'infinite. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_plagioclase, data_liquid, P_kbar, MWinfo)
% calcTemp
% Calculate Putirka (2005) Table 2, Model A for one selected
% Plagioclase-Liquid pair and a scalar or vector of pressures.
%
% Inputs:
%   data_plagioclase : one-row Plagioclase table
%   data_liquid      : one-row liquid table
%   P_kbar           : finite non-negative pressure column vector
%   MWinfo           : molecular-weight and cation-number structure
%
% Output:
%   row : table with one row per pressure value, including intermediate
%         components and final temperatures.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

% --- Plagioclase components ---
pl = preparePlagRow(data_plagioclase, MWinfo);

% --- Liquid oxides (wt%) ---
% Missing optional columns are interpreted as zero. Existing NaN values are
% never replaced by zero.
SiO2  = getLiqOxRequired(data_liquid, {'SiO2'});
TiO2  = getLiqOxOptional(data_liquid, {'TiO2'}, 0);
Al2O3 = getLiqOxRequired(data_liquid, {'Al2O3'});
FeO   = getLiqOxOptional(data_liquid, {'FeO', 'FeOt'}, 0);
MnO   = getLiqOxOptional(data_liquid, {'MnO'}, 0);
MgO   = getLiqOxOptional(data_liquid, {'MgO'}, 0);
CaO   = getLiqOxRequired(data_liquid, {'CaO'});
Na2O  = getLiqOxRequired(data_liquid, {'Na2O'});
K2O   = getLiqOxOptional(data_liquid, {'K2O'}, 0);
Cr2O3 = getLiqOxOptional(data_liquid, {'Cr2O3'}, 0);
H2O   = getLiqOxOptional(data_liquid, {'H2O', 'H2Ot'}, 0);

% --- Liquid cation proportions on the published anhydrous basis ---
% Only the ten components explicitly listed by Putirka (2005, p. 341) are
% included in the cation-fraction denominator.
nSiO2  = SiO2  .* MWinfo.Cat.SiO2  ./ MWinfo.MW.SiO2;
nTiO2  = TiO2  .* MWinfo.Cat.TiO2  ./ MWinfo.MW.TiO2;
nAl2O3 = Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
nFeO   = FeO   .* MWinfo.Cat.FeO   ./ MWinfo.MW.FeO;
nMnO   = MnO   .* MWinfo.Cat.MnO   ./ MWinfo.MW.MnO;
nMgO   = MgO   .* MWinfo.Cat.MgO   ./ MWinfo.MW.MgO;
nCaO   = CaO   .* MWinfo.Cat.CaO   ./ MWinfo.MW.CaO;
nNa2O  = Na2O  .* MWinfo.Cat.Na2O  ./ MWinfo.MW.Na2O;
nK2O   = K2O   .* MWinfo.Cat.K2O   ./ MWinfo.MW.K2O;
nCr2O3 = Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;

cationTotal_liq = nSiO2 + nTiO2 + nAl2O3 + nFeO + nMnO + ...
    nMgO + nCaO + nNa2O + nK2O + nCr2O3;

XSiO2_liq   = nSiO2  ./ cationTotal_liq;
XTiO2_liq   = nTiO2  ./ cationTotal_liq;
XAlO1_5_liq = nAl2O3 ./ cationTotal_liq;
XFeO_liq    = nFeO   ./ cationTotal_liq;
XMnO_liq    = nMnO   ./ cationTotal_liq;
XMgO_liq    = nMgO   ./ cationTotal_liq;
XCaO_liq    = nCaO   ./ cationTotal_liq;
XNaO0_5_liq = nNa2O  ./ cationTotal_liq;
XKO0_5_liq  = nK2O   ./ cationTotal_liq;
XCrO1_5_liq = nCr2O3 ./ cationTotal_liq;

% --- Putirka (2005) Model A terms ---
An_pl = pl.An;
Ab_pl = pl.Ab;
Ca_liq = XCaO_liq;
Al_liq = XAlO1_5_liq;
Si_liq = XSiO2_liq;
H2O_liq = H2O;

lnK_ModelA = log(An_pl ./ ...
    (Ca_liq .* (Al_liq .^ 2) .* (Si_liq .^ 2)));
AlRatio_ModelA = Al_liq ./ (Al_liq + Si_liq);

denom_ModelA = ...
    6.12 ...
    + 0.257 .* lnK_ModelA ...
    - 3.166 .* Ca_liq ...
    + 0.2166 .* H2O_liq ...
    - 3.137 .* AlRatio_ModelA ...
    + 1.216 .* (Ab_pl .^ 2) ...
    - 2.475e-2 .* P_kbar;

% No NaN/Inf replacement is performed. MATLAB arithmetic propagates NaN and
% retains Inf values, which are reported by the calling function.
TModelA_K = 1e4 ./ denom_ModelA;
TModelA_C = TModelA_K - 273.15;

% --- Pack outputs ---
row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

row.SiO2_pl = repmat(pl.SiO2, nP, 1);
row.TiO2_pl = repmat(pl.TiO2, nP, 1);
row.Al2O3_pl = repmat(pl.Al2O3, nP, 1);
row.FeO_pl = repmat(pl.FeO, nP, 1);
row.MnO_pl = repmat(pl.MnO, nP, 1);
row.MgO_pl = repmat(pl.MgO, nP, 1);
row.CaO_pl = repmat(pl.CaO, nP, 1);
row.Na2O_pl = repmat(pl.Na2O, nP, 1);
row.K2O_pl = repmat(pl.K2O, nP, 1);
row.Cr2O3_pl = repmat(pl.Cr2O3, nP, 1);

row.XSi_pl = repmat(pl.XSi, nP, 1);
row.XTi_pl = repmat(pl.XTi, nP, 1);
row.XAl_pl = repmat(pl.XAl, nP, 1);
row.XFe_pl = repmat(pl.XFe, nP, 1);
row.XMn_pl = repmat(pl.XMn, nP, 1);
row.XMg_pl = repmat(pl.XMg, nP, 1);
row.XCa_pl = repmat(pl.XCa, nP, 1);
row.XNa_pl = repmat(pl.XNa, nP, 1);
row.XK_pl = repmat(pl.XK, nP, 1);
row.XCr_pl = repmat(pl.XCr, nP, 1);
row.cationSum_pl = repmat(pl.cationSum, nP, 1);
row.An_pl = repmat(pl.An, nP, 1);
row.Ab_pl = repmat(pl.Ab, nP, 1);
row.Or_pl = repmat(pl.Or, nP, 1);

row.SiO2_liq = repmat(SiO2, nP, 1);
row.TiO2_liq = repmat(TiO2, nP, 1);
row.Al2O3_liq = repmat(Al2O3, nP, 1);
row.FeO_liq = repmat(FeO, nP, 1);
row.MnO_liq = repmat(MnO, nP, 1);
row.MgO_liq = repmat(MgO, nP, 1);
row.CaO_liq = repmat(CaO, nP, 1);
row.Na2O_liq = repmat(Na2O, nP, 1);
row.K2O_liq = repmat(K2O, nP, 1);
row.Cr2O3_liq = repmat(Cr2O3, nP, 1);
row.H2O_liq_wt = repmat(H2O, nP, 1);

row.cationTotal_liq = repmat(cationTotal_liq, nP, 1);
row.XSiO2_liq = repmat(XSiO2_liq, nP, 1);
row.XTiO2_liq = repmat(XTiO2_liq, nP, 1);
row.XAlO1_5_liq = repmat(XAlO1_5_liq, nP, 1);
row.XFeO_liq = repmat(XFeO_liq, nP, 1);
row.XMnO_liq = repmat(XMnO_liq, nP, 1);
row.XMgO_liq = repmat(XMgO_liq, nP, 1);
row.XCaO_liq = repmat(XCaO_liq, nP, 1);
row.XNaO0_5_liq = repmat(XNaO0_5_liq, nP, 1);
row.XKO0_5_liq = repmat(XKO0_5_liq, nP, 1);
row.XCrO1_5_liq = repmat(XCrO1_5_liq, nP, 1);

row.term_An_pl = repmat(An_pl, nP, 1);
row.term_Ab_pl = repmat(Ab_pl, nP, 1);
row.term_Ca_liq = repmat(Ca_liq, nP, 1);
row.term_Al_liq = repmat(Al_liq, nP, 1);
row.term_Si_liq = repmat(Si_liq, nP, 1);
row.term_H2O_liq = repmat(H2O_liq, nP, 1);
row.lnK_ModelA = repmat(lnK_ModelA, nP, 1);
row.AlRatio_ModelA = repmat(AlRatio_ModelA, nP, 1);
row.denom_ModelA = denom_ModelA;

% Standardized names used by ThermoCalc launchers and plotting routines.
row.T_K = TModelA_K;
row.T_degreeC = TModelA_C;
row.T_deg = TModelA_C;

% Legacy names retained for backward compatibility with existing output
% processing scripts.
row.TModelA_K = TModelA_K;
row.TModelA_C = TModelA_C;

end

function pl = preparePlagRow(data_plagioclase, MWinfo)
% preparePlagRow
% Convert one Plagioclase oxide row to cations on an 8-oxygen basis and
% calculate An-Ab-Or components. Existing NaN values are retained.

SiO2  = getMineralOxRequired(data_plagioclase, {'SiO2'});
TiO2  = getMineralOxOptional(data_plagioclase, {'TiO2'}, 0);
Al2O3 = getMineralOxRequired(data_plagioclase, {'Al2O3'});
FeO   = getMineralOxOptional(data_plagioclase, {'FeO', 'FeOt'}, 0);
MnO   = getMineralOxOptional(data_plagioclase, {'MnO'}, 0);
MgO   = getMineralOxOptional(data_plagioclase, {'MgO'}, 0);
CaO   = getMineralOxRequired(data_plagioclase, {'CaO'});
Na2O  = getMineralOxRequired(data_plagioclase, {'Na2O'});
K2O   = getMineralOxOptional(data_plagioclase, {'K2O'}, 0);
Cr2O3 = getMineralOxOptional(data_plagioclase, {'Cr2O3'}, 0);

molSiO2  = SiO2  ./ MWinfo.MW.SiO2;
molTiO2  = TiO2  ./ MWinfo.MW.TiO2;
molAl2O3 = Al2O3 ./ MWinfo.MW.Al2O3;
molFeO   = FeO   ./ MWinfo.MW.FeO;
molMnO   = MnO   ./ MWinfo.MW.MnO;
molMgO   = MgO   ./ MWinfo.MW.MgO;
molCaO   = CaO   ./ MWinfo.MW.CaO;
molNa2O  = Na2O  ./ MWinfo.MW.Na2O;
molK2O   = K2O   ./ MWinfo.MW.K2O;
molCr2O3 = Cr2O3 ./ MWinfo.MW.Cr2O3;

oxygenSum = ...
    2 .* molSiO2 ...
    + 2 .* molTiO2 ...
    + 3 .* molAl2O3 ...
    + molFeO ...
    + molMnO ...
    + molMgO ...
    + molCaO ...
    + molNa2O ...
    + molK2O ...
    + 3 .* molCr2O3;

oxygenRenormalizationFactor = 8 ./ oxygenSum;

XSi = molSiO2 .* oxygenRenormalizationFactor;
XTi = molTiO2 .* oxygenRenormalizationFactor;
XAl = 2 .* molAl2O3 .* oxygenRenormalizationFactor;
XFe = molFeO .* oxygenRenormalizationFactor;
XMn = molMnO .* oxygenRenormalizationFactor;
XMg = molMgO .* oxygenRenormalizationFactor;
XCa = molCaO .* oxygenRenormalizationFactor;
XNa = 2 .* molNa2O .* oxygenRenormalizationFactor;
XK  = 2 .* molK2O .* oxygenRenormalizationFactor;
XCr = 2 .* molCr2O3 .* oxygenRenormalizationFactor;

cationSum = XSi + XTi + XAl + XFe + XMn + XMg + XCa + XNa + XK + XCr;

feldsparSiteSum = XCa + XNa + XK;
An = XCa ./ feldsparSiteSum;
Ab = XNa ./ feldsparSiteSum;
Or = XK ./ feldsparSiteSum;

pl = struct();
pl.SiO2 = SiO2;
pl.TiO2 = TiO2;
pl.Al2O3 = Al2O3;
pl.FeO = FeO;
pl.MnO = MnO;
pl.MgO = MgO;
pl.CaO = CaO;
pl.Na2O = Na2O;
pl.K2O = K2O;
pl.Cr2O3 = Cr2O3;
pl.XSi = XSi;
pl.XTi = XTi;
pl.XAl = XAl;
pl.XFe = XFe;
pl.XMn = XMn;
pl.XMg = XMg;
pl.XCa = XCa;
pl.XNa = XNa;
pl.XK = XK;
pl.XCr = XCr;
pl.cationSum = cationSum;
pl.An = An;
pl.Ab = Ab;
pl.Or = Or;

end

function row = attachLiquidIDs(row, data_liquid)
% attachLiquidIDs
% Copy commonly used liquid identifiers to every pressure row.

nRows = height(row);
variableNames = data_liquid.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    row.liq_Index = repmat(data_liquid.('Index'), nRows, 1);
end
if any(strcmp(variableNames, 'Experiment'))
    row.liq_Experiment = repmat(string(data_liquid.('Experiment')), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    row.liq_Citation = repmat(string(data_liquid.('Citation')), nRows, 1);
end

end

function validateRequiredColumns(data_table, aliasGroups, phaseName)
% validateRequiredColumns
% Confirm that each required oxide has at least one accepted column alias.

missingBuffer = strings(numel(aliasGroups), 1);
nMissing = 0;

for i = 1:numel(aliasGroups)
    aliases = aliasGroups{i};
    [~, matchedName] = findOxideColumnAliases( ...
        data_table.Properties.VariableNames, aliases);
    if isempty(matchedName)
        nMissing = nMissing + 1;
        missingBuffer(nMissing) = strjoin(string(aliases), '/');
    end
end

if nMissing > 0
    missingNames = missingBuffer(1:nMissing);
    error('%s table is missing required oxide column(s): %s.', ...
        phaseName, char(strjoin(missingNames, ', ')));
end

end

function value = getMineralOxRequired(data_table, aliases)
% getMineralOxRequired
% Read a required mineral oxide. Missing columns stop the calculation;
% existing NaN values are returned unchanged.

[value, ~, ~] = getOxideValue( ...
    data_table, aliases, NaN, true, 'Plagioclase');

end

function value = getMineralOxOptional(data_table, aliases, defaultValue)
% getMineralOxOptional
% Read an optional mineral oxide. defaultValue is used only when no accepted
% column alias exists; existing NaN values are returned unchanged.

[value, ~, ~] = getOxideValue( ...
    data_table, aliases, defaultValue, false, 'Plagioclase');

end

function value = getLiqOxRequired(data_table, aliases)
% getLiqOxRequired
% Read a required liquid oxide. Missing columns stop the calculation;
% existing NaN values are returned unchanged.

[value, ~, ~] = getOxideValue( ...
    data_table, aliases, NaN, true, 'Liquid');

end

function value = getLiqOxOptional(data_table, aliases, defaultValue)
% getLiqOxOptional
% Read an optional liquid oxide. defaultValue is used only when no accepted
% column alias exists; existing NaN values are returned unchanged.

[value, ~, ~] = getOxideValue( ...
    data_table, aliases, defaultValue, false, 'Liquid');

end

function [value, matchedName, columnFound] = getOxideValue( ...
        data_table, aliases, defaultValue, isRequired, phaseName)
% getOxideValue
% Resolve one of several oxide-column aliases and convert the selected
% one-row table value to a scalar double without replacing NaN.

if ischar(aliases) || isstring(aliases)
    aliases = cellstr(string(aliases));
end

[~, matchedName] = findOxideColumnAliases( ...
    data_table.Properties.VariableNames, aliases);
columnFound = ~isempty(matchedName);

if ~columnFound
    if isRequired
        error('%s table must contain oxide column: %s.', ...
            phaseName, char(strjoin(string(aliases), ' or ')));
    end
    value = defaultValue;
    return
end

value = toScalarDouble(data_table.(matchedName));

end

function [columnIndex, matchedName] = findOxideColumnAliases(varNames, aliases)
% findOxideColumnAliases
% Find the first table variable matching one of the supplied oxide aliases.
% Matching is case-insensitive and ignores spaces, underscores, and hyphens.

canonicalVarNames = strings(numel(varNames), 1);
for i = 1:numel(varNames)
    canonicalVarNames(i) = canonicalizeName(varNames{i});
end

columnIndex = [];
matchedName = '';

for i = 1:numel(aliases)
    canonicalAlias = canonicalizeName(aliases{i});
    targets = [canonicalAlias + "value", canonicalAlias];

    for j = 1:numel(targets)
        idx = find(canonicalVarNames == targets(j), 1, 'first');
        if ~isempty(idx)
            columnIndex = idx;
            matchedName = varNames{idx};
            return
        end
    end
end

end

function canonicalName = canonicalizeName(inputName)
% canonicalizeName
% Convert a variable or oxide name to the canonical form used for matching.

canonicalName = lower(string(inputName));
canonicalName = replace(canonicalName, " ", "");
canonicalName = replace(canonicalName, "_", "");
canonicalName = replace(canonicalName, "-", "");

end

function value = toScalarDouble(rawValue)
% toScalarDouble
% Convert the first value from a selected one-row table variable to double.
% Missing, empty, or non-convertible values become NaN. NaN is never replaced
% by a default value in this function.

value = NaN;

if isempty(rawValue)
    return
end

if isnumeric(rawValue) || islogical(rawValue)
    value = double(rawValue(1));
    return
end

if iscell(rawValue)
    if isempty(rawValue{1})
        return
    end
    rawValue = rawValue{1};
end

if iscategorical(rawValue)
    rawValue = string(rawValue(1));
elseif isstring(rawValue)
    if ismissing(rawValue(1))
        return
    end
    rawValue = rawValue(1);
elseif ischar(rawValue)
    rawValue = string(rawValue);
end

if isnumeric(rawValue) || islogical(rawValue)
    value = double(rawValue(1));
    return
end

convertedValue = str2double(string(rawValue));
if ~isempty(convertedValue)
    value = convertedValue(1);
end

end
