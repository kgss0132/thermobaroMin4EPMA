function results = Hoisch1990(rawdata_struct, T_degreeC)
% functions/+baro/+Mica/Hoisch1990.m
% Tested with MATLAB R2024b
%
% Empirical geobarometers R1-R6 for mica-bearing metapelites
% Hoisch, T.D. (1990)
% Contributions to Mineralogy and Petrology, 104, 225-234
% DOI: https://doi.org/10.1007/BF00306445
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Muscovite analysis, one Biotite
% analysis, one Garnet analysis, and one Plagioclase analysis, and calculates
% pressure using the six empirical geobarometers R1-R6 of Hoisch (1990).
%
% The function accepts temperature as either a scalar or a vector. It is
% therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected mineral combination, one output row
% is returned for every input temperature value.
%
% Six pressure estimates are retained separately as P_R1_kbar through
% P_R6_kbar. For compatibility with the standardized range-mode P-T plot,
% P_kbar is also supplied as the unweighted arithmetic mean of R1-R6. This
% mean is an implementation-level summary and is NOT a separate calibration
% published by Hoisch (1990). NaN is not omitted from this mean; therefore,
% if any R1-R6 pressure is NaN, P_kbar also remains NaN.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Hoisch (1990) empirically calibrated the six equilibria using 43 natural
% samples containing quartz + muscovite + biotite + garnet + plagioclase +
% Al2SiO5 (kyanite or sillimanite). The calibration data in Table 2 span:
%
%   Temperature : 468-752 degreeC
%   Pressure    : 2059-9585 bar (2.059-9.585 kbar)
%
% The calibration dataset and its P-T range are listed in Table 2 on p. 226.
% The abstract on p. 225 states that applications should be restricted to
% rocks whose equilibrium constants and compositional variables fall within
% the ranges used for calibration.
%
% Hoisch (1990) gives the following calibration-variable ranges in Table 9
% on p. 232:
%
%   ln(KR1) :  3.61 to  6.62
%   ln(KR2) :  1.85 to  4.89
%   ln(KR3) :  1.53 to  5.46
%   ln(KR4) : -0.92 to  1.62
%   ln(KR5) :  6.60 to 11.66
%   ln(KR6) :  2.02 to  6.31
%   XMg_M   : 0.00851 to 0.04330   (Muscovite octahedral Mg fraction)
%   XMg_B   : 0.31900 to 0.41700   (Biotite octahedral Mg fraction)
%   XFe_B   : 0.33900 to 0.44300   (Biotite octahedral Fe fraction)
%   XAl_B   : 0.08940 to 0.17000   (Biotite octahedral Al fraction)
%   XTi_B   : 0.02400 to 0.07320   (Biotite octahedral Ti fraction)
%
% Important limitations discussed on pp. 230-234:
%   1) The calibration is empirical and inherits systematic uncertainty from
%      the Garnet-Biotite thermometer and GASP barometer used to assign the
%      calibration P-T conditions (discussion of uncertainties, p. 231).
%   2) R1-R4 are especially sensitive to extrapolation of Muscovite and
%      Biotite compositions. The calibration micas are Al-rich and occupy
%      narrow compositional ranges (pp. 232-233).
%   3) R3-R4 use a Henry-law treatment for the MgAl-celadonite component in
%      Muscovite and should be restricted to compositions similar to the
%      calibration Muscovites (p. 233).
%   4) R1-R2 use a symmetric mixing model for Biotite. Extrapolation outside
%      the narrow calibration range is explicitly questioned (p. 233).
%   5) R5-R6 use ideal mica-component mixing, but the safest procedure for
%      all six geobarometers is still to remain within Table 9 (p. 233).
%   6) The six P-T lines and the Garnet-Biotite exchange thermometer should
%      converge. Strong divergence or non-overlapping prediction intervals
%      may indicate disequilibrium or inappropriate compositions (p. 233).
%   7) Applications to higher-variance assemblages with fewer saturating
%      phases may be adversely affected by the missing phases (p. 234).
%
% Reaction-specific phase requirements are summarized on p. 234:
%   R1-R2 : Quartz + Garnet + Biotite + Plagioclase
%   R3    : Quartz + Garnet + Muscovite + Plagioclase
%   R4    : Quartz + Garnet + Muscovite + Biotite + Plagioclase
%   R5-R6 : Garnet + Muscovite + Biotite + Plagioclase; quartz saturation
%           may still be important for successful application.
%
% Structural-formula assumptions from pp. 225-226 and pp. 231-232:
%   - Micas were normalized to 11 anhydrous oxygens.
%   - All Fe in Muscovite and Biotite was treated as Fe2+.
%   - Garnet Fe3+ was estimated by charge balance using 8 cations and
%     12 oxygens in the original calibration dataset.
%   - Octahedral vacancies in Biotite, a possible small trioctahedral
%     Muscovite component, mica polytype effects, and selected interlayer or
%     tetrahedral mixing terms were not explicitly represented.
%
% Regression residual standard deviations are 0.092-0.338 kbar, whereas
% propagated analytical standard deviations are approximately 0.208-0.331
% kbar (Table 7, p. 229). These values describe precision within the adopted
% empirical framework and should not be interpreted as the total absolute
% accuracy of every natural-rock application (discussion, pp. 231-232).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 468-752 degreeC,
%   2) finite calculated pressure is outside 2.059-9.585 kbar,
%   3) a Table 9 compositional variable or ln(KR1)-ln(KR6) is outside its
%      calibration range,
%   4) a required calculation input contains NaN, or
%   5) a calculated pressure is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Mica        : table used for both Muscovite and Biotite
%   rawdata_struct.Garnet      : table
%   rawdata_struct.Plagioclase : table
%
% The FIRST column of each table is treated as the identifier (data code)
% displayed in the selection dialog.
%
% Required Mica variables:
%   Mg_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   Fe_cation_apfu     % treated as Fe2+ following Hoisch (1990)
%   Ti_cation_apfu
%
% Required Garnet variables:
%   Mg_cation_apfu
%   Fe_cation_apfu     % Fe2+ used in the implemented activity model
%   Ca_cation_apfu
%   Mn_cation_apfu
%
% Required Plagioclase variables:
%   Ca_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%
% Optional Mica variables retained in the output when present:
%   Na_cation_apfu
%   K_cation_apfu
%
% Values used in the calculation must be real numeric scalars. Finite values
% must be non-negative. NaN is allowed, retained, propagated through the
% calculation, and reported by fprintf. Inf and finite negative values are
% prohibited. Zero is retained; division by zero or log(0) may consequently
% produce NaN or Inf, which remains in the output and is reported.
%
% No liquid composition is used by this geobarometer. Therefore, exclusion
% of F and Cl from cationTotal_liq and from liquid NaN warnings is not
% applicable to Hoisch1990.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION
% Activities, equilibrium constants, grossular-volume correction, and the
% six calibrated pressure expressions follow Tables 3, 4, and 8 on
% pp. 227-231 of Hoisch (1990). Pressure is calculated in bars and converted
% to kbar by division by 1000. Temperature is used in K.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Hoisch1990(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing Mica, Garnet, and Plagioclase tables
%   T_degreeC      : non-negative real numeric scalar or vector in degreeC.
%                    NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per input temperature value for every
%             user-selected mineral combination.
%

%% Input validation
% Accept scalar or vector temperature input so that fixed-temperature and
% temperature-range launchers use the same implementation.
if nargin < 2
    error('Hoisch1990 requires (rawdata_struct, T_degreeC).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(T_degreeC) || ~isreal(T_degreeC) || isempty(T_degreeC) || ...
        ~isvector(T_degreeC) || any(isinf(T_degreeC(:))) || ...
        any(isfinite(T_degreeC(:)) & T_degreeC(:) < 0)
    error(['T_degreeC must be a non-negative real numeric scalar or vector. ' ...
           'NaN is allowed and retained, but Inf and finite negative values are prohibited.']);
end

T_degreeC = T_degreeC(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. Source tables are not
% modified; selected rows are validated immediately before calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end
if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Plagioclase') || ~istable(rawdata_struct.Plagioclase)
    error('rawdata_struct must contain table: rawdata_struct.Plagioclase');
end

% The same Mica table is presented separately for Muscovite and Biotite
% selection, matching the original implementation.
dataset_ms = rawdata_struct.Mica;
dataset_bt = rawdata_struct.Mica;
dataset_grt = rawdata_struct.Garnet;
dataset_plg = rawdata_struct.Plagioclase;

validateRequiredVariables(dataset_ms, dataset_grt, dataset_plg);

dataCodes_ms = dataset_ms{:, 1};
dataCodes_bt = dataset_bt{:, 1};
dataCodes_grt = dataset_grt{:, 1};
dataCodes_plg = dataset_plg{:, 1};

dataCodeList_ms = cellstr(string(dataCodes_ms));
dataCodeList_bt = cellstr(string(dataCodes_bt));
dataCodeList_grt = cellstr(string(dataCodes_grt));
dataCodeList_plg = cellstr(string(dataCodes_plg));

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container and calibration limits
% Store each selected combination as one table block and concatenate only
% once after the interactive loop. The buffer grows geometrically rather
% than changing size on every iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degreeC = 468;
calibrationT_max_degreeC = 752;
calibrationP_min_kbar = 2.059;
calibrationP_max_kbar = 9.585;

lnK_min = [3.61, 1.85, 1.53, -0.92, 6.60, 2.02];
lnK_max = [6.62, 4.89, 5.46,  1.62, 11.66, 6.31];

compositionNames = {'XMg_M', 'XMg_B', 'XFe_B', 'XAl_B', 'XTi_B'};
compositionMin = [0.00851, 0.31900, 0.33900, 0.08940, 0.02400];
compositionMax = [0.04330, 0.41700, 0.44300, 0.17000, 0.07320];

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);
temperatureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-6) Interactive selection loop, calculation, and diagnostics
while true
    % ----- Muscovite selection -----
    disp('=== Step 3: Selecting a data code from the list (Muscovite) ===');

    [selectedIdx_ms, ok] = listdlg( ...
        'PromptString', 'Please select the Muscovite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodeList_ms, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_ms)
        disp('Selection canceled');
        break;
    end

    selectedCode_ms = dataCodes_ms(selectedIdx_ms);
    disp(['Muscovite selected: ' char(string(selectedCode_ms))]);

    % ----- Biotite selection -----
    disp('=== Step 4: Selecting a data code from the list (Biotite) ===');

    [selectedIdx_bt, ok] = listdlg( ...
        'PromptString', 'Please select the Biotite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodeList_bt, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_bt)
        disp('Selection canceled');
        break;
    end

    selectedCode_bt = dataCodes_bt(selectedIdx_bt);
    disp(['Biotite selected: ' char(string(selectedCode_bt))]);

    % ----- Garnet selection -----
    disp('=== Step 5: Selecting a data code from the list (Garnet) ===');

    [selectedIdx_grt, ok] = listdlg( ...
        'PromptString', 'Please select the Garnet data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodeList_grt, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_grt)
        disp('Selection canceled');
        break;
    end

    selectedCode_grt = dataCodes_grt(selectedIdx_grt);
    disp(['Garnet selected: ' char(string(selectedCode_grt))]);

    % ----- Plagioclase selection -----
    disp('=== Step 6: Selecting a data code from the list (Plagioclase) ===');

    [selectedIdx_plg, ok] = listdlg( ...
        'PromptString', 'Please select the Plagioclase data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodeList_plg, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_plg)
        disp('Selection canceled');
        break;
    end

    selectedCode_plg = dataCodes_plg(selectedIdx_plg);
    disp(['Plagioclase selected: ' char(string(selectedCode_plg))]);

    % ----- Prepare selected rows -----
    selectedData_ms = dataset_ms(selectedIdx_ms, :);
    selectedData_bt = dataset_bt(selectedIdx_bt, :);
    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_plg = dataset_plg(selectedIdx_plg, :);

    validateNonNegativeInputs(selectedData_ms, selectedData_bt, ...
        selectedData_grt, selectedData_plg);

    nanInputNames = findNaNInputs(selectedData_ms, selectedData_bt, ...
        selectedData_grt, selectedData_plg, T_degreeC);

    % ----- Calculation -----
    disp('=== Step 7: Calculating pressure using Hoisch (1990) R1-R6 ===');

    row = calcPressure(selectedData_ms, selectedData_bt, ...
        selectedData_grt, selectedData_plg, T_degreeC);

    % Store identifiers once for every temperature row.
    row.dataCode_ms = repmat(string(selectedCode_ms), height(row), 1);
    row.dataCode_bt = repmat(string(selectedCode_bt), height(row), 1);
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_plg = repmat(string(selectedCode_plg), height(row), 1);

    row = movevars(row, ...
        {'dataCode_ms', 'dataCode_bt', 'dataCode_grt', 'dataCode_plg'}, ...
        'Before', 1);

    % Store this calculation as one buffered table block.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % ----- Echo pressure results -----
    disp('--------------------------------------------------');
    disp('=== Pressures were calculated: ===');

    for reactionIndex = 1:6
        variableName = sprintf('P_R%d_kbar', reactionIndex);
        printResultRange(sprintf('R%d', reactionIndex), row.(variableName), 'kbar');
    end

    printResultRange('Mean P_kbar', row.P_kbar, 'kbar');

    % ----- Calibration-range and NaN diagnostics -----
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));

        fprintf(2, ...
            ['WARNING: Input temperature is outside the Hoisch (1990) ' ...
             'calibration-data range of 468-752 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'finite input range = %.4g-%.4g degreeC (Table 2, p. 226).\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));

        temperatureWarningIssued = true;
    end

    % Warn for each individual pressure equation outside the calibration
    % pressure envelope. Calculation is not stopped.
    for reactionIndex = 1:6
        variableName = sprintf('P_R%d_kbar', reactionIndex);
        pressureValues = row.(variableName);
        finitePressure = isfinite(pressureValues);
        outsidePressure = finitePressure & ...
            (pressureValues < calibrationP_min_kbar | ...
             pressureValues > calibrationP_max_kbar);

        if any(outsidePressure)
            finiteValues = pressureValues(finitePressure);

            fprintf(2, ...
                ['WARNING: Hoisch (1990) R%d calculated pressure is outside ' ...
                 'the calibration-data range of 2.059-9.585 kbar. ' ...
                 '%d of %d finite pressure point(s) are outside the range; ' ...
                 'finite calculated range = %.4g-%.4g kbar for %s, %s, %s, and %s ' ...
                 '(Table 2, p. 226).\n'], ...
                reactionIndex, ...
                sum(outsidePressure), ...
                sum(finitePressure), ...
                min(finiteValues), ...
                max(finiteValues), ...
                char(string(selectedCode_ms)), ...
                char(string(selectedCode_bt)), ...
                char(string(selectedCode_grt)), ...
                char(string(selectedCode_plg)));
        end
    end

    % Warn when selected mica compositional variables fall outside Table 9.
    for compositionIndex = 1:numel(compositionNames)
        variableName = compositionNames{compositionIndex};
        variableValues = row.(variableName);
        value = variableValues(1);

        if isfinite(value) && ...
                (value < compositionMin(compositionIndex) || ...
                 value > compositionMax(compositionIndex))
            fprintf(2, ...
                ['WARNING: %s = %.6g is outside the Hoisch (1990) calibration ' ...
                 'range %.6g-%.6g (Table 9, p. 232) for %s, %s, %s, and %s.\n'], ...
                variableName, ...
                value, ...
                compositionMin(compositionIndex), ...
                compositionMax(compositionIndex), ...
                char(string(selectedCode_ms)), ...
                char(string(selectedCode_bt)), ...
                char(string(selectedCode_grt)), ...
                char(string(selectedCode_plg)));
        end
    end

    % Warn when temperature-dependent equilibrium constants fall outside
    % their Table 9 calibration ranges.
    for reactionIndex = 1:6
        variableName = sprintf('lnKR%d', reactionIndex);
        values = row.(variableName);
        finiteValuesMask = isfinite(values);
        outsideValues = finiteValuesMask & ...
            (values < lnK_min(reactionIndex) | ...
             values > lnK_max(reactionIndex));

        if any(outsideValues)
            finiteValues = values(finiteValuesMask);

            fprintf(2, ...
                ['WARNING: ln(KR%d) is outside the Hoisch (1990) calibration ' ...
                 'range %.4g-%.4g. %d of %d finite value(s) are outside; ' ...
                 'finite calculated range = %.4g-%.4g for %s, %s, %s, and %s ' ...
                 '(Table 9, p. 232).\n'], ...
                reactionIndex, ...
                lnK_min(reactionIndex), ...
                lnK_max(reactionIndex), ...
                sum(outsideValues), ...
                sum(finiteValuesMask), ...
                min(finiteValues), ...
                max(finiteValues), ...
                char(string(selectedCode_ms)), ...
                char(string(selectedCode_bt)), ...
                char(string(selectedCode_grt)), ...
                char(string(selectedCode_plg)));
        end
    end

    % Report all required inputs that were NaN. NaN remains unchanged and
    % propagates through the equations.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Hoisch (1990) calculation input(s) ' ...
             'for %s, %s, %s, and %s: %s.\n' ...
             '         NaN was retained and propagated; calculated pressure may be NaN.\n'], ...
            char(string(selectedCode_ms)), ...
            char(string(selectedCode_bt)), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_plg)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report non-finite results separately for each reaction.
    for reactionIndex = 1:6
        variableName = sprintf('P_R%d_kbar', reactionIndex);
        pressureValues = row.(variableName);
        invalidPressure = ~isfinite(pressureValues);

        if any(invalidPressure)
            fprintf(2, ...
                ['WARNING: Non-finite pressure values were calculated by ' ...
                 'Hoisch (1990) R%d for %s, %s, %s, and %s ' ...
                 '(%d of %d points; NaN: %d, Inf: %d).\n' ...
                 '         These values remain in the output table, and the calculation was not stopped.\n'], ...
                reactionIndex, ...
                char(string(selectedCode_ms)), ...
                char(string(selectedCode_bt)), ...
                char(string(selectedCode_grt)), ...
                char(string(selectedCode_plg)), ...
                sum(invalidPressure), ...
                numel(pressureValues), ...
                sum(isnan(pressureValues)), ...
                sum(isinf(pressureValues)));
        end
    end

    if any(~isfinite(row.P_kbar))
        fprintf(2, ...
            ['WARNING: The launcher-compatible mean P_kbar contains non-finite ' ...
             'value(s) because at least one of R1-R6 is non-finite. ' ...
             'NaN was not omitted from the arithmetic mean.\n']);
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Hoisch1990', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks once after all selections are complete.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ------------------------------------------------------------------------
function row = calcPressure(data_ms, data_bt, data_grt, data_plg, T_degreeC)
% Calculate Hoisch (1990) R1-R6 for one selected mineral combination and a
% scalar or vector of temperatures.

R = 8.3144;
T_degreeC = T_degreeC(:);
T_K = T_degreeC + 273.15;
nT = numel(T_degreeC);

ms = prepareMicaRow(data_ms, 'Muscovite');
bt = prepareMicaRow(data_bt, 'Biotite');
grt = prepareGarnetRow(data_grt);
plg = preparePlagioclaseRow(data_plg);

%% Mica site fractions
XMg_M = ms.Mg ./ 2;
XAl_M = (ms.Al + ms.Si - 4) ./ 2;

XMg_B = bt.Mg ./ 3;
XFe_B = bt.Fe2 ./ 3;
XTi_B = bt.Ti ./ 3;
XAl_B = (bt.Al + bt.Si - 4) ./ 3;

Xcel = 4 .* XMg_M .* XAl_M;
Xms = XAl_M .^ 2;

Xphl = XMg_B .^ 3;
Xann = XFe_B .^ 3;
Xeas = 6.75 .* (XMg_B .^ 2) .* XAl_B;
Xsid = 6.75 .* (XFe_B .^ 2) .* XAl_B;

%% Garnet mole fractions
% Direct arithmetic is retained. A zero denominator gives NaN through 0/0,
% and a NaN input remains NaN rather than being replaced by zero.
denGrt = grt.Mg + grt.Fe2 + grt.Ca + grt.Mn;

Xpy = safeDivide(grt.Mg, denGrt);
Xal = safeDivide(grt.Fe2, denGrt);
Xgr = safeDivide(grt.Ca, denGrt);
Xsp = safeDivide(grt.Mn, denGrt);

%% Garnet activities
a_py = (Xpy .* exp(-((13809 - 6.28 .* T_K) ./ (R .* T_K)) .* ...
    (Xgr .^ 2 + Xal .* Xgr + Xgr .* Xsp))) .^ 3;

a_al = (Xal .* exp(((6.28 .* T_K - 13809) ./ (R .* T_K)) .* ...
    (Xpy .* Xgr))) .^ 3;

a_gr = (Xgr .* exp(-((13809 - 6.28 .* T_K) ./ (R .* T_K)) .* ...
    (Xpy .^ 2 + Xal .* Xpy + Xpy .* Xsp))) .^ 3;

%% Plagioclase anorthite activity
denPlg = plg.Ca + plg.Na + plg.K;
Xan = safeDivide(plg.Ca, denPlg);

a_an = (Xan .* (1 + Xan) .^ 2 ./ 4) .* ...
    exp(((1 - Xan) .^ 2 ./ (R .* T_K)) .* (8578 + 39300 .* Xan));

%% Grossular volume correction
DeltaVgr = calcDeltaVgr(Xal, Xpy, Xgr);

%% Equilibrium constants
KR1 = safeDivide(Xphl .* a_an .^ 2, ...
    a_py .^ (1/3) .* a_gr .^ (2/3) .* Xeas);
KR2 = safeDivide(Xann .* a_an .^ 2, ...
    a_al .^ (1/3) .* a_gr .^ (2/3) .* Xsid);
KR3 = safeDivide(Xcel .* a_an .^ 2, ...
    a_py .^ (1/3) .* a_gr .^ (2/3) .* Xms);
KR4 = safeDivide(Xcel .* a_an, ...
    Xphl .^ (1/3) .* a_gr .^ (1/3) .* Xms .^ (2/3));
KR5 = safeDivide(Xphl .* a_an .^ 3, a_py .* a_gr .* Xms);
KR6 = safeDivide(Xann .* a_an .^ 3, a_al .* a_gr .* Xms);

lnKR1 = safeLog(KR1);
lnKR2 = safeLog(KR2);
lnKR3 = safeLog(KR3);
lnKR4 = safeLog(KR4);
lnKR5 = safeLog(KR5);
lnKR6 = safeLog(KR6);

%% Hoisch (1990) calibrated pressure expressions, P in bars
P_R1_bar = safeDivide( ...
    -31830.6 + 79.0281 .* T_K - R .* T_K .* lnKR1 ...
    - 26968.7 .* (XAl_B - XMg_B) ...
    + 32604.5 .* XFe_B ...
    + 42855.4 .* XTi_B, ...
    3.8145 - (2 ./ 3) .* DeltaVgr);

P_R2_bar = safeDivide( ...
    -46707.2 + 85.5824 .* T_K - R .* T_K .* lnKR2 ...
    - 30960.2 .* (XAl_B - XFe_B) ...
    + 24289.6 .* XMg_B ...
    + 37265.6 .* XTi_B, ...
    3.8986 - (2 ./ 3) .* DeltaVgr);

P_R3_bar = safeDivide( ...
    -20681.4 + 69.8341 .* T_K - R .* T_K .* lnKR3 ...
    - 185443 .* (XMg_M .* (XMg_M - 2)), ...
    4.1740 - (2 ./ 3) .* DeltaVgr);

P_R4_bar = safeDivide( ...
    -21664.0 + 33.7500 .* T_K - R .* T_K .* lnKR4 ...
    - 172479 .* (XMg_M .* (XMg_M - 2)), ...
    2.19415 - (1 ./ 3) .* DeltaVgr);

P_R5_bar = safeDivide( ...
    -3546.01 + 121.347 .* T_K - R .* T_K .* lnKR5, ...
    6.37161 - DeltaVgr);

P_R6_bar = safeDivide( ...
    -55530.4 + 140.635 .* T_K - R .* T_K .* lnKR6, ...
    6.59940 - DeltaVgr);

P_R1_kbar = P_R1_bar ./ 1000;
P_R2_kbar = P_R2_bar ./ 1000;
P_R3_kbar = P_R3_bar ./ 1000;
P_R4_kbar = P_R4_bar ./ 1000;
P_R5_kbar = P_R5_bar ./ 1000;
P_R6_kbar = P_R6_bar ./ 1000;

% Standardized launcher-compatible output. NaN is deliberately not omitted.
pressureMatrix = [P_R1_kbar, P_R2_kbar, P_R3_kbar, ...
    P_R4_kbar, P_R5_kbar, P_R6_kbar];
P_kbar = mean(pressureMatrix, 2);
P_MPa = P_kbar .* 100;

%% Pack results
% Temperature-dependent values remain vectors. Composition-dependent values
% are explicitly expanded to nT rows to avoid implicit table-size changes.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.P_kbar = P_kbar;
row.P_MPa = P_MPa;

row.P_R1_kbar = P_R1_kbar;
row.P_R2_kbar = P_R2_kbar;
row.P_R3_kbar = P_R3_kbar;
row.P_R4_kbar = P_R4_kbar;
row.P_R5_kbar = P_R5_kbar;
row.P_R6_kbar = P_R6_kbar;

row.P_R1_bar = P_R1_bar;
row.P_R2_bar = P_R2_bar;
row.P_R3_bar = P_R3_bar;
row.P_R4_bar = P_R4_bar;
row.P_R5_bar = P_R5_bar;
row.P_R6_bar = P_R6_bar;

row.Mg_ms = expandScalar(ms.Mg, nT);
row.Al_ms = expandScalar(ms.Al, nT);
row.Si_ms = expandScalar(ms.Si, nT);
row.Fe2_ms = expandScalar(ms.Fe2, nT);
row.Ti_ms = expandScalar(ms.Ti, nT);
row.Na_ms = expandScalar(ms.Na, nT);
row.K_ms = expandScalar(ms.K, nT);

row.Mg_bt = expandScalar(bt.Mg, nT);
row.Al_bt = expandScalar(bt.Al, nT);
row.Si_bt = expandScalar(bt.Si, nT);
row.Fe2_bt = expandScalar(bt.Fe2, nT);
row.Ti_bt = expandScalar(bt.Ti, nT);
row.Na_bt = expandScalar(bt.Na, nT);
row.K_bt = expandScalar(bt.K, nT);

row.Mg_grt = expandScalar(grt.Mg, nT);
row.Fe2_grt = expandScalar(grt.Fe2, nT);
row.Ca_grt = expandScalar(grt.Ca, nT);
row.Mn_grt = expandScalar(grt.Mn, nT);

row.Ca_plg = expandScalar(plg.Ca, nT);
row.Na_plg = expandScalar(plg.Na, nT);
row.K_plg = expandScalar(plg.K, nT);

row.XMg_M = expandScalar(XMg_M, nT);
row.XAl_M = expandScalar(XAl_M, nT);
row.XMg_B = expandScalar(XMg_B, nT);
row.XFe_B = expandScalar(XFe_B, nT);
row.XTi_B = expandScalar(XTi_B, nT);
row.XAl_B = expandScalar(XAl_B, nT);

row.Xcel = expandScalar(Xcel, nT);
row.Xms = expandScalar(Xms, nT);
row.Xphl = expandScalar(Xphl, nT);
row.Xann = expandScalar(Xann, nT);
row.Xeas = expandScalar(Xeas, nT);
row.Xsid = expandScalar(Xsid, nT);

row.Xpy = expandScalar(Xpy, nT);
row.Xal = expandScalar(Xal, nT);
row.Xgr = expandScalar(Xgr, nT);
row.Xsp = expandScalar(Xsp, nT);
row.Xan = expandScalar(Xan, nT);

row.a_py = a_py;
row.a_al = a_al;
row.a_gr = a_gr;
row.a_an = a_an;
row.DeltaVgr_J_per_bar = expandScalar(DeltaVgr, nT);

row.KR1 = KR1;
row.KR2 = KR2;
row.KR3 = KR3;
row.KR4 = KR4;
row.KR5 = KR5;
row.KR6 = KR6;

row.lnKR1 = lnKR1;
row.lnKR2 = lnKR2;
row.lnKR3 = lnKR3;
row.lnKR4 = lnKR4;
row.lnKR5 = lnKR5;
row.lnKR6 = lnKR6;

end

%% ------------------------------------------------------------------------
function DeltaVgr = calcDeltaVgr(Xal, Xpy, Xgr)
% Calculate the grossular partial-molar-volume correction from Table 3.

z1 = (1 - Xgr - 0.914) ./ 0.066;
z2 = (1 - Xgr - 0.940) ./ 0.083;

Val = 125.24 + 1.482 .* (1 - Xgr) .^ 2 ...
    - 0.48 .* (1 + z1 .* ((1 - Xgr) ./ 0.066) .* exp(-(z1 .^ 2) ./ 2));

Vpy = 125.24 + 0.512 .* (1 - Xgr) .^ 2 ...
    - 0.418 .* (1 + z2 .* ((1 - Xgr) ./ 0.083) .* exp(-(z2 .^ 2) ./ 2));

den = Xal + Xpy;

Vgr = 0.1 .* (safeDivide(Val .* Xal, den) + ...
    safeDivide(Vpy .* Xpy, den));

DeltaVgr = Vgr - 12.53;

end

%% ------------------------------------------------------------------------
function mica = prepareMicaRow(data_tbl, mineralLabel)
% Extract one selected Mica row while retaining present NaN values.

if height(data_tbl) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mica = struct();
mica.Mg = getRequiredScalar(data_tbl, 'Mg_cation_apfu', mineralLabel);
mica.Al = getRequiredScalar(data_tbl, 'Al_cation_apfu', mineralLabel);
mica.Si = getRequiredScalar(data_tbl, 'Si_cation_apfu', mineralLabel);
mica.Fe2 = getRequiredScalar(data_tbl, 'Fe_cation_apfu', mineralLabel);
mica.Ti = getRequiredScalar(data_tbl, 'Ti_cation_apfu', mineralLabel);
mica.Na = getOptionalScalar(data_tbl, 'Na_cation_apfu', mineralLabel);
mica.K = getOptionalScalar(data_tbl, 'K_cation_apfu', mineralLabel);

end

%% ------------------------------------------------------------------------
function grt = prepareGarnetRow(data_tbl)
% Extract one selected Garnet row while retaining present NaN values.

if height(data_tbl) ~= 1
    error('Garnet input must be a 1-row table.');
end

grt = struct();
grt.Mg = getRequiredScalar(data_tbl, 'Mg_cation_apfu', 'Garnet');
grt.Fe2 = getRequiredScalar(data_tbl, 'Fe_cation_apfu', 'Garnet');
grt.Ca = getRequiredScalar(data_tbl, 'Ca_cation_apfu', 'Garnet');
grt.Mn = getRequiredScalar(data_tbl, 'Mn_cation_apfu', 'Garnet');

end

%% ------------------------------------------------------------------------
function plg = preparePlagioclaseRow(data_tbl)
% Extract one selected Plagioclase row while retaining present NaN values.

if height(data_tbl) ~= 1
    error('Plagioclase input must be a 1-row table.');
end

plg = struct();
plg.Ca = getRequiredScalar(data_tbl, 'Ca_cation_apfu', 'Plagioclase');
plg.Na = getRequiredScalar(data_tbl, 'Na_cation_apfu', 'Plagioclase');
plg.K = getRequiredScalar(data_tbl, 'K_cation_apfu', 'Plagioclase');

end

%% ------------------------------------------------------------------------
function validateRequiredVariables(dataset_mica, dataset_grt, dataset_plg)
% Confirm required calculation columns before opening selection dialogs.

requiredMica = { ...
    'Mg_cation_apfu', ...
    'Al_cation_apfu', ...
    'Si_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Ti_cation_apfu'};

requiredGarnet = { ...
    'Mg_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Mn_cation_apfu'};

requiredPlagioclase = { ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu'};

assertTableVariables(dataset_mica, requiredMica, 'Mica');
assertTableVariables(dataset_grt, requiredGarnet, 'Garnet');
assertTableVariables(dataset_plg, requiredPlagioclase, 'Plagioclase');

end

%% ------------------------------------------------------------------------
function assertTableVariables(tbl, requiredNames, label)
% Raise one clear error listing all absent required variables.

missingNames = requiredNames(~ismember(requiredNames, tbl.Properties.VariableNames));

if ~isempty(missingNames)
    error('%s table is missing required variable(s): %s', ...
        label, strjoin(missingNames, ', '));
end

end

%% ------------------------------------------------------------------------
function validateNonNegativeInputs(data_ms, data_bt, data_grt, data_plg)
% Reject Inf, complex, non-scalar, and finite negative calculation inputs.
% NaN and zero are retained.

validateTableValues(data_ms, ...
    {'Mg_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu'}, ...
    'Muscovite');

validateTableValues(data_bt, ...
    {'Mg_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu', ...
     'Fe_cation_apfu', 'Ti_cation_apfu'}, ...
    'Biotite');

validateTableValues(data_grt, ...
    {'Mg_cation_apfu', 'Fe_cation_apfu', 'Ca_cation_apfu', ...
     'Mn_cation_apfu'}, ...
    'Garnet');

validateTableValues(data_plg, ...
    {'Ca_cation_apfu', 'Na_cation_apfu', 'K_cation_apfu'}, ...
    'Plagioclase');

end

%% ------------------------------------------------------------------------
function validateTableValues(tbl, variableNames, label)
% Validate one-row values used numerically by the geobarometer.

for variableIndex = 1:numel(variableNames)
    variableName = variableNames{variableIndex};
    value = tbl.(variableName);

    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
        error('%s.%s must be a real numeric scalar.', label, variableName);
    end

    if isinf(value)
        error('%s.%s must not be Inf.', label, variableName);
    end

    if isfinite(value) && value < 0
        error('%s.%s must be non-negative.', label, variableName);
    end
end

end

%% ------------------------------------------------------------------------
function nanInputNames = findNaNInputs(data_ms, data_bt, data_grt, data_plg, T_degreeC)
% Return names of required calculation inputs whose selected value is NaN.
% The maximum possible list length is allocated before scanning the inputs.

msVariables = {'Mg_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu'};
btVariables = {'Mg_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu', ...
    'Fe_cation_apfu', 'Ti_cation_apfu'};
grtVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', 'Ca_cation_apfu', ...
    'Mn_cation_apfu'};
plgVariables = {'Ca_cation_apfu', 'Na_cation_apfu', 'K_cation_apfu'};

maximumNameCount = numel(msVariables) + numel(btVariables) + ...
    numel(grtVariables) + numel(plgVariables) + numel(T_degreeC);

nanInputNames = strings(maximumNameCount, 1);
nNames = 0;

[nanInputNames, nNames] = appendNaNTableNames( ...
    nanInputNames, nNames, data_ms, msVariables, 'Muscovite');
[nanInputNames, nNames] = appendNaNTableNames( ...
    nanInputNames, nNames, data_bt, btVariables, 'Biotite');
[nanInputNames, nNames] = appendNaNTableNames( ...
    nanInputNames, nNames, data_grt, grtVariables, 'Garnet');
[nanInputNames, nNames] = appendNaNTableNames( ...
    nanInputNames, nNames, data_plg, plgVariables, 'Plagioclase');

nanTemperatureIndices = find(isnan(T_degreeC));

for indexPosition = 1:numel(nanTemperatureIndices)
    nNames = nNames + 1;
    nanInputNames(nNames) = ...
        "T_degreeC(" + string(nanTemperatureIndices(indexPosition)) + ")";
end

nanInputNames = nanInputNames(1:nNames);

end

%% ------------------------------------------------------------------------
function [names, nNames] = appendNaNTableNames( ...
        names, nNames, tbl, variableNames, label)
% Add qualified names for selected NaN table values to a preallocated list.

for variableIndex = 1:numel(variableNames)
    variableName = variableNames{variableIndex};

    if isnan(tbl.(variableName))
        nNames = nNames + 1;
        names(nNames) = string(label) + "." + string(variableName);
    end
end

end

%% ------------------------------------------------------------------------
function value = getRequiredScalar(tbl, varName, mineralLabel)
% Retrieve a required scalar. NaN is retained; Inf and negatives are rejected.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);

if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('%s.%s must be a real numeric scalar.', mineralLabel, varName);
end
if isinf(value)
    error('%s.%s must not be Inf.', mineralLabel, varName);
end
if isfinite(value) && value < 0
    error('%s.%s must be non-negative.', mineralLabel, varName);
end

end

%% ------------------------------------------------------------------------
function value = getOptionalScalar(tbl, varName, mineralLabel)
% Retrieve an optional output-only scalar. Missing variables are stored as
% NaN rather than zero so that missing information is not misrepresented.

if ~ismember(varName, tbl.Properties.VariableNames)
    value = NaN;
    return;
end

value = tbl.(varName);

if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('%s.%s must be a real numeric scalar.', mineralLabel, varName);
end
if isinf(value)
    error('%s.%s must not be Inf.', mineralLabel, varName);
end
if isfinite(value) && value < 0
    error('%s.%s must be non-negative.', mineralLabel, varName);
end

end

%% ------------------------------------------------------------------------
function y = safeDivide(numerator, denominator)
% Perform element-wise division with MATLAB's native NaN/Inf propagation.
% Implicit expansion supports scalar composition terms and temperature vectors.

y = numerator ./ denominator;

end

%% ------------------------------------------------------------------------
function y = safeLog(x)
% Retain log(0) as -Inf, positive values as real logs, and NaN as NaN.
% Negative derived equilibrium constants are outside the real logarithm
% domain and are converted to NaN rather than complex values.

y = NaN(size(x));
nonNegative = x >= 0;
y(nonNegative) = log(x(nonNegative));
y(isnan(x)) = NaN;

end

%% ------------------------------------------------------------------------
function values = expandScalar(value, nRows)
% Expand one scalar to the row count required by the temperature vector.

values = repmat(value, nRows, 1);

end

%% ------------------------------------------------------------------------
function printResultRange(label, values, unitText)
% Print scalar or vector result summaries without changing stored values.

if isscalar(values)
    fprintf('%s = %.10g %s\n', label, values, unitText);
    return;
end

finiteValues = values(isfinite(values));

if isempty(finiteValues)
    fprintf('%s = no finite values (%d point(s))\n', label, numel(values));
else
    fprintf('%s = %.10g to %.10g %s (%d point(s))\n', ...
        label, min(finiteValues), max(finiteValues), unitText, numel(values));
end

end
