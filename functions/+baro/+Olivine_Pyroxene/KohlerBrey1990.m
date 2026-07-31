function results = KohlerBrey1990(rawdata_struct, T_degreeC)
% functions/+baro/+Olivine_Pyroxene/KohlerBrey1990.m
% Tested with MATLAB R2024b
%
% Ca-in-Olivine / Clinopyroxene geobarometer
% Köhler, T.P. and Brey, G.P. (1990)
% Geochimica et Cosmochimica Acta, 54, 2375-2388
% DOI: https://doi.org/10.1016/0016-7037(90)90226-B
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis and one
% Clinopyroxene analysis and calculates pressure using the Ca exchange
% between coexisting olivine and clinopyroxene calibrated by Köhler and
% Brey (1990).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Olivine-Clinopyroxene pair, one
% output row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Köhler and Brey (1990) calibrated the barometer using experiments in
% natural lherzolitic compositions over the following ranges:
%
%   Temperature : 900-1400 degreeC
%   Pressure    : 2-60 kbar
%   Olivine     : Mg# = 0.87-1.00, approximately Fo87-100
%
% The experimental pressure-temperature ranges are stated in the abstract
% and experimental discussion (pp. 2375, 2376, and 2380). The olivine Mg#
% restriction is discussed on p. 2383. Application below 900 degreeC is an
% extrapolation beyond the experimental temperature range (p. 2380).
%
% Equations (5) and (6) reproduce the experimental pressures to
% approximately +/-1.65 kbar (1 sigma) when temperature is independently
% known (p. 2383). When the barometer is combined with the Brey-Köhler
% two-pyroxene thermometer, the reported pressure reproducibility decreases
% to approximately +/-5.4 kbar (1 sigma; p. 2384).
%
% IMPORTANT APPLICATION LIMITATIONS:
%
% 1) The barometer is strongly temperature dependent. Accurate temperature
%    estimates are essential (p. 2384).
%
% 2) Olivine and clinopyroxene must represent an equilibrium pair. Because
%    Ca diffuses approximately 3-4 orders of magnitude faster in olivine
%    than in clinopyroxene, short-lived cooling can produce pressure
%    overestimates and short-lived heating can produce pressure
%    underestimates (pp. 2375 and 2385-2386).
%
% 3) Ca zoning in olivine must be checked. Ca-rich outer rims, commonly
%    observed within approximately 30-50 micrometres of grain boundaries,
%    may record heating during xenolith transport rather than mantle
%    equilibration (pp. 2386-2387).
%
% 4) Phase-boundary fluorescence can artificially increase measured Ca in
%    olivine when analyses are made close to Ca-rich phases such as
%    clinopyroxene. Köhler and Brey separated mineral grains and used WDS
%    trace-element analyses to minimize this effect (pp. 2379-2380).
%
% 5) Fertile harzburgites containing exsolved clinopyroxene, rocks with
%    multiple pyroxene generations, reaction textures, or evidence of
%    partial re-equilibration may yield geologically meaningless pressures
%    (pp. 2384-2386).
%
% 6) The transition between the two empirical equations is pressure
%    dependent rather than fixed at 1100 degreeC:
%
%      T_split(K) = 1275.25 + 2.827*P(kbar)
%
%    Selection of Eq. (5) or Eq. (6) must therefore be checked for
%    self-consistency. If both candidate equations are valid or neither is
%    valid, this implementation retains the candidate nearest to its branch
%    boundary and issues a non-stopping fprintf warning.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 900-1400 degreeC,
%   2) finite calculated pressure is outside 2-60 kbar,
%   3) finite olivine Mg# is outside 0.87-1.00,
%   4) required or available calculation inputs contain NaN,
%   5) a calculated pressure is NaN or Inf,
%   6) a calculated pressure is negative, or
%   7) equation-branch selection is ambiguous.
%
% These checks do not establish mineral equilibrium. Petrographic context,
% zoning, and analytical quality must be assessed independently.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%
% and either:
%   rawdata_struct.Cpx             : table
% or
%   rawdata_struct.Clinopyroxene   : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The selected rows must contain the
% following normalized-cation variables, or one of the accepted aliases:
%
%   Olivine:
%     Ca_cation_apfu, Ca_cation, or Ca
%
%   Clinopyroxene:
%     Ca_cation_apfu, Ca_cation, or Ca
%
% Optional variables used only to evaluate Mg# are:
%
%   Mg_cation_apfu, Mg_cation, or Mg
%   Fe2_cation_apfu, Fe_cation_apfu, Fe2_cation, Fe_cation,
%   Fe2, or Fe
%
% Ca_ol and Ca_cpx must be atomic proportions in structural formulae based
% on 4 oxygens for olivine and 6 oxygens for clinopyroxene, respectively
% (pp. 2375 and 2382-2383).
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if it makes a ratio or logarithm
% undefined, the resulting NaN or Inf is retained and reported.
%
% No liquid composition is used by this barometer. Therefore, Liq F and Cl,
% cationTotal_liq, and the exclusion of F and Cl from cationTotal_liq are not
% applicable to this function. F and Cl are consequently not included in
% any NaN-input warning.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION
%
% Define:
%
%   Dca = Ca_ol / Ca_cpx
%
% where Ca_ol and Ca_cpx are atomic proportions of Ca in olivine and
% clinopyroxene structural formulae based on 4 and 6 oxygens, respectively.
% Temperature T is in Kelvin and pressure is returned in kbar.
%
% Köhler and Brey (1990), Eq. (5), high-temperature branch (p. 2383):
%
%   P(kbar) = [-T*ln(Dca) - 11982 + 3.61*T] / 56.2
%
%   valid when T >= 1275.25 + 2.827*P
%
% Köhler and Brey (1990), Eq. (6), low-temperature branch (p. 2383):
%
%   P(kbar) = [-T*ln(Dca) - 5792 - 1.25*T] / 42.5
%
%   valid when T <= 1275.25 + 2.827*P
%
% -------------------------------------------------------------------------
% Syntax:
%   results = KohlerBrey1990(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing Olivine and Clinopyroxene tables
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Olivine-Clinopyroxene pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('KohlerBrey1990 requires (rawdata_struct, T_degreeC).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(T_degreeC) || isempty(T_degreeC) || ~isvector(T_degreeC) || ...
        any(isinf(T_degreeC(:))) || ...
        any(isfinite(T_degreeC(:)) & T_degreeC(:) < 0)
    error(['T_degreeC must be a non-negative numeric scalar or vector. ' ...
           'NaN is allowed and retained, but Inf and finite negative values are prohibited.']);
end

T_degreeC = T_degreeC(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The source tables are
% not modified; selected rows are read during each calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end

if isfield(rawdata_struct, 'Cpx') && istable(rawdata_struct.Cpx)
    dataset_cpx = rawdata_struct.Cpx;
elseif isfield(rawdata_struct, 'Clinopyroxene') && ...
        istable(rawdata_struct.Clinopyroxene)
    dataset_cpx = rawdata_struct.Clinopyroxene;
else
    error(['rawdata_struct must contain table: rawdata_struct.Cpx or ' ...
           'rawdata_struct.Clinopyroxene']);
end

dataset_ol = rawdata_struct.Olivine;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Experimental calibration limits from Köhler and Brey (1990).
experimentalT_min_degreeC = 900;
experimentalT_max_degreeC = 1400;
experimentalP_min_kbar = 2;
experimentalP_max_kbar = 60;
olivineMg_min = 0.87;
olivineMg_max = 1.00;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < experimentalT_min_degreeC | ...
     T_degreeC > experimentalT_max_degreeC);
temperatureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Olivine) ===');

while true
    % ----- Olivine selection -----
    % The first table column is used only as the displayed data identifier.
    dataCodes_ol = dataset_ol{:, 1};

    [selectedIdx_ol, ok] = listdlg( ...
        'PromptString', 'Please select the Olivine data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_ol)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_ol)
        disp('Selection canceled');
        break;
    end

    selectedCode_ol = dataCodes_ol(selectedIdx_ol);
    disp(['Olivine selected: ' char(string(selectedCode_ol))]);

    % ----- Clinopyroxene selection -----
    disp('=== Step 4: Selecting a data code from the list (Clinopyroxene) ===');

    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Clinopyroxene data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIdx_cpx);
    disp(['Clinopyroxene selected: ' char(string(selectedCode_cpx))]);

    % ----- Calculation -----
    % Olivine and clinopyroxene are selected independently; row numbers are
    % not assumed to correspond between the two tables.
    disp('=== Step 5: Calculating the pressure ===');

    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    ol = prepareOlivineRow(selectedData_ol);
    cpx = prepareClinopyroxeneRow(selectedData_cpx);

    % List NaN values without changing them. Optional Mg and Fe variables
    % are listed only when the corresponding input column exists.
    nanInputNames = findNaNInputs(ol, cpx, T_degreeC);

    % Reject only Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(ol, cpx);

    row = calcPressure(ol, cpx, T_degreeC);

    % Repeat identifiers for all temperatures in the current calculation.
    row.dataCode_olivine = repmat(string(selectedCode_ol), height(row), 1);
    row.dataCode_clinopyroxene = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, ...
        {'dataCode_olivine','dataCode_clinopyroxene'}, 'Before', 1);

    % Store one block per selected mineral pair. The buffer grows only when
    % its reserved capacity is exhausted, rather than on every iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure or finite pressure range.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    finitePressureValues = row.P_kbar(isfinite(row.P_kbar));
    if height(row) == 1
        disp([char(string(selectedCode_ol)) ' & ' ...
            char(string(selectedCode_cpx)) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    elseif ~isempty(finitePressureValues)
        disp([char(string(selectedCode_ol)) ' & ' ...
            char(string(selectedCode_cpx)) ': finite P range = ' ...
            num2str(min(finitePressureValues)) ' to ' ...
            num2str(max(finitePressureValues)) ' kbar']);
    else
        disp([char(string(selectedCode_ol)) ' & ' ...
            char(string(selectedCode_cpx)) ': all calculated P values are non-finite']);
    end

    % Input temperature is common to all selected pairs, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperatureValues = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the 900-1400 degreeC ' ...
             'experimental calibration range of Köhler and Brey (1990; ' ...
             'pp. 2375-2376 and 2380). %d of %d finite temperature point(s) ' ...
             'are outside the range; input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperatureValues), ...
            min(finiteTemperatureValues), ...
            max(finiteTemperatureValues));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures lie outside the experimental
    % pressure calibration range.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < experimentalP_min_kbar | ...
         row.P_kbar > experimentalP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the 2-60 kbar ' ...
             'experimental calibration range of Köhler and Brey (1990; ' ...
             'pp. 2375-2376 and 2380). %d of %d finite pressure point(s) ' ...
             'are outside the range; calculated finite range = %.4g-%.4g ' ...
             'kbar for %s & %s.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_cpx)));
    end

    % Warn when olivine Mg# is outside the compositional range of the
    % calibration, or when it cannot be evaluated.
    XMg_ol_scalar = row.XMg_ol(1);
    if isfinite(XMg_ol_scalar) && ...
            (XMg_ol_scalar < olivineMg_min || XMg_ol_scalar > olivineMg_max)
        fprintf(2, ...
            ['WARNING: Olivine Mg# = %.4g is outside the Mg# 0.87-1.00 ' ...
             'application range of Köhler and Brey (1990; p. 2383) for ' ...
             '%s & %s.\n'], ...
            XMg_ol_scalar, ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_cpx)));
    elseif ~isfinite(XMg_ol_scalar)
        fprintf(2, ...
            ['WARNING: Olivine Mg# could not be evaluated for %s & %s because ' ...
             'Mg and/or Fe is missing, NaN, or has a zero sum. The Mg# ' ...
             '0.87-1.00 application criterion (p. 2383) could not be checked.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_cpx)));
    end

    % List the exact available calculation inputs containing NaN. For vector
    % temperature input, the NaN element indices are included.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s & %s: %s.\n' ...
             '         NaN values were retained and propagated; calculated ' ...
             'pressure values may remain NaN.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Warn when the Ca ratio or logarithm is non-positive or non-finite.
    invalidDca = ~isfinite(row.Dca_ol_cpx) | row.Dca_ol_cpx <= 0;
    invalidLnDca = ~isfinite(row.lnDca);
    if any(invalidDca | invalidLnDca)
        fprintf(2, ...
            ['WARNING: Dca = Ca_ol/Ca_cpx or ln(Dca) is non-positive or ' ...
             'non-finite for %s & %s. Zero and NaN inputs were retained; ' ...
             'the corresponding pressure result may be NaN or Inf.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_cpx)));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, +Inf: %d, -Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_cpx)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar) & row.P_kbar > 0), ...
            sum(isinf(row.P_kbar) & row.P_kbar < 0));
    end

    % Negative finite pressures are retained for diagnosis but are outside
    % the calibration range and have already triggered the range warning.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & %s ' ...
             '(%d of %d points). The values were retained for diagnostic ' ...
             'purposes.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_cpx)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    % Report ambiguous branch selection. The selected candidate is retained
    % and identified explicitly in the output branch column.
    if any(row.isBranchAmbiguous)
        fprintf(2, ...
            ['WARNING: Equation-branch selection was ambiguous for %s & %s ' ...
             'at %d of %d temperature point(s). In those rows, both candidate ' ...
             'equations were valid or neither candidate was valid. The candidate ' ...
             'nearest to its pressure-dependent branch boundary was retained. ' ...
             'Inspect branch, P_highT_Eq5_kbar, P_lowT_Eq6_kbar, and the ' ...
             'boundary-temperature columns.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_cpx)), ...
            sum(row.isBranchAmbiguous), ...
            numel(row.isBranchAmbiguous));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'KohlerBrey1990', ...
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

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(ol, cpx, T_degreeC)
% findNaNInputs
% Return the names of available barometer inputs containing NaN. NaN values
% do not cause an error and are not replaced by zero.

maxNames = 7;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = "T_degreeC(indices=" + indexText + ")";
end

if isnan(ol.Ca)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Olivine." + ol.CaVariableName;
end
if isnan(cpx.Ca)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Clinopyroxene." + cpx.CaVariableName;
end

if ol.MgFound && isnan(ol.Mg)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Olivine." + ol.MgVariableName;
end
if ol.FeFound && isnan(ol.Fe)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Olivine." + ol.FeVariableName;
end
if cpx.MgFound && isnan(cpx.Mg)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Clinopyroxene." + cpx.MgVariableName;
end
if cpx.FeFound && isnan(cpx.Fe)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Clinopyroxene." + cpx.FeVariableName;
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(ol, cpx)
% validateNonNegativeInputs
% Reject Inf and finite negative values in extracted cation variables. Zero
% and NaN are intentionally allowed and retained.

maxNames = 6;
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

olivineFields = {'Ca', 'Mg', 'Fe'};
for i = 1:numel(olivineFields)
    fieldName = olivineFields{i};
    value = ol.(fieldName);
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Olivine." + getSourceVariableName(ol, fieldName);
    end
end

cpxFields = {'Ca', 'Mg', 'Fe'};
for i = 1:numel(cpxFields)
    fieldName = cpxFields{i};
    value = cpx.(fieldName);
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Clinopyroxene." + getSourceVariableName(cpx, fieldName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['KohlerBrey1990: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(ol, cpx, T_degreeC)
% calcPressure
% Compute pressure for one olivine row and one clinopyroxene row at one or
% more input temperatures. NaN and Inf generated by zero or NaN inputs are
% retained for diagnosis.
%
% Inputs:
%   ol          : prepared one-row Olivine structure
%   cpx         : prepared one-row Clinopyroxene structure
%   T_degreeC   : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Expand composition-dependent scalars to the temperature-vector length.
Ca_ol = repmat(ol.Ca, nT, 1);
Mg_ol = repmat(ol.Mg, nT, 1);
Fe2_ol = repmat(ol.Fe, nT, 1);

Ca_cpx = repmat(cpx.Ca, nT, 1);
Mg_cpx = repmat(cpx.Mg, nT, 1);
Fe2_cpx = repmat(cpx.Fe, nT, 1);

% Ca exchange coefficient. NaN and zero values propagate naturally. Because
% finite negative inputs are prohibited, log does not produce complex output.
Dca = Ca_ol ./ Ca_cpx;
lnDca = log(Dca);

% Köhler and Brey (1990), Eq. (5): high-temperature branch.
P_highT_kbar = (-T_K .* lnDca - 11982 + 3.61 .* T_K) ./ 56.2;

% Köhler and Brey (1990), Eq. (6): low-temperature branch.
P_lowT_kbar = (-T_K .* lnDca - 5792 - 1.25 .* T_K) ./ 42.5;

% Evaluate each candidate pressure against its own pressure-dependent branch
% boundary. Comparisons with NaN return false; +/-Inf generated from zero
% inputs is retained when a branch condition can still be evaluated.
boundary_high_K = 1275.25 + 2.827 .* P_highT_kbar;
boundary_low_K = 1275.25 + 2.827 .* P_lowT_kbar;

isHighTBranchValid = T_K >= boundary_high_K;
isLowTBranchValid = T_K <= boundary_low_K;

P_kbar = NaN(nT, 1);
branch = strings(nT, 1);
isBranchAmbiguous = false(nT, 1);

onlyHigh = isHighTBranchValid & ~isLowTBranchValid;
onlyLow = isLowTBranchValid & ~isHighTBranchValid;
bothValid = isHighTBranchValid & isLowTBranchValid;
neitherValid = ~isHighTBranchValid & ~isLowTBranchValid;

P_kbar(onlyHigh) = P_highT_kbar(onlyHigh);
branch(onlyHigh) = "highT_Eq5";

P_kbar(onlyLow) = P_lowT_kbar(onlyLow);
branch(onlyLow) = "lowT_Eq6";

boundaryDistanceHigh = abs(T_K - boundary_high_K);
boundaryDistanceLow = abs(T_K - boundary_low_K);

chooseHighBoth = bothValid & ...
    boundaryDistanceHigh <= boundaryDistanceLow;
chooseLowBoth = bothValid & ...
    boundaryDistanceHigh > boundaryDistanceLow;

P_kbar(chooseHighBoth) = P_highT_kbar(chooseHighBoth);
branch(chooseHighBoth) = "highT_Eq5_bothValid";

P_kbar(chooseLowBoth) = P_lowT_kbar(chooseLowBoth);
branch(chooseLowBoth) = "lowT_Eq6_bothValid";

finiteBoundaryDistances = ...
    isfinite(boundaryDistanceHigh) & isfinite(boundaryDistanceLow);
chooseHighNeither = neitherValid & finiteBoundaryDistances & ...
    boundaryDistanceHigh <= boundaryDistanceLow;
chooseLowNeither = neitherValid & finiteBoundaryDistances & ...
    boundaryDistanceHigh > boundaryDistanceLow;

P_kbar(chooseHighNeither) = P_highT_kbar(chooseHighNeither);
branch(chooseHighNeither) = "highT_Eq5_nearestBoundary";

P_kbar(chooseLowNeither) = P_lowT_kbar(chooseLowNeither);
branch(chooseLowNeither) = "lowT_Eq6_nearestBoundary";

unresolvedBranch = neitherValid & ~finiteBoundaryDistances;
branch(unresolvedBranch) = "unresolved_nonfinite";

isBranchAmbiguous(bothValid | chooseHighNeither | chooseLowNeither) = true;

% Mg# values are diagnostic composition checks and do not enter the pressure
% equation. Missing or NaN Mg/Fe values produce NaN Mg# rather than zero.
XMg_ol_scalar = calculateMgNumber(ol.Mg, ol.Fe);
XMg_cpx_scalar = calculateMgNumber(cpx.Mg, cpx.Fe);
XMg_ol = repmat(XMg_ol_scalar, nT, 1);
XMg_cpx = repmat(XMg_cpx_scalar, nT, 1);

% Applicability flags. These are diagnostic and do not stop calculation.
isWithinExperimentalTRange = isfinite(T_degreeC) & ...
    T_degreeC >= 900 & T_degreeC <= 1400;
isWithinExperimentalPRange = isfinite(P_kbar) & ...
    P_kbar >= 2 & P_kbar <= 60;
isWithinOlivineMgRange = isfinite(XMg_ol) & ...
    XMg_ol >= 0.87 & XMg_ol <= 1.00;
isRecommended = isWithinExperimentalTRange & ...
    isWithinExperimentalPRange & isWithinOlivineMgRange & ...
    ~isBranchAmbiguous & isfinite(P_kbar);

% Pack outputs using equal-height, pre-sized vectors.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.Ca_ol = Ca_ol;
row.Ca_cpx = Ca_cpx;
row.Dca_ol_cpx = Dca;
row.lnDca = lnDca;

row.P_highT_Eq5_kbar = P_highT_kbar;
row.P_lowT_Eq6_kbar = P_lowT_kbar;

row.boundaryT_highT_Eq5_K = boundary_high_K;
row.boundaryT_lowT_Eq6_K = boundary_low_K;
row.boundaryDistance_highT_Eq5_K = boundaryDistanceHigh;
row.boundaryDistance_lowT_Eq6_K = boundaryDistanceLow;

row.isHighTBranchValid = isHighTBranchValid;
row.isLowTBranchValid = isLowTBranchValid;
row.isBranchAmbiguous = isBranchAmbiguous;

row.P_kbar = P_kbar;
row.branch = branch;

row.Mg_ol = Mg_ol;
row.Fe2_ol = Fe2_ol;
row.XMg_ol = XMg_ol;

row.Mg_cpx = Mg_cpx;
row.Fe2_cpx = Fe2_cpx;
row.XMg_cpx = XMg_cpx;

row.P_uncertainty_equation_1sigma_kbar = repmat(1.65, nT, 1);
row.P_uncertainty_with_TBKN_1sigma_kbar = repmat(5.4, nT, 1);

row.isWithinExperimentalTRange_900_1400C = isWithinExperimentalTRange;
row.isWithinExperimentalPRange_2_60kbar = isWithinExperimentalPRange;
row.isWithinOlivineMgRange_0p87_1p00 = isWithinOlivineMgRange;
row.isRecommended = isRecommended;

% Backward-compatible aliases from the original implementation.
row.is_exp_temperature_range_900_1400C = isWithinExperimentalTRange;
row.is_exp_pressure_range_2_60kbar = isWithinExperimentalPRange;
row.is_olivine_Mg_range_Fo87_100 = isWithinOlivineMgRange;

end

function ol = prepareOlivineRow(data_olivine)
% prepareOlivineRow
% Extract one-row olivine cation data without replacing NaN by zero.

if height(data_olivine) ~= 1
    error('Olivine input must be a 1-row table.');
end

ol = struct();
[ol.Ca, ol.CaVariableName, ol.CaFound] = getVarByAliases( ...
    data_olivine, {'Ca_cation_apfu','Ca_cation','Ca'}, true, 'Olivine');
[ol.Mg, ol.MgVariableName, ol.MgFound] = getVarByAliases( ...
    data_olivine, {'Mg_cation_apfu','Mg_cation','Mg'}, false, 'Olivine');
[ol.Fe, ol.FeVariableName, ol.FeFound] = getVarByAliases( ...
    data_olivine, ...
    {'Fe2_cation_apfu','Fe_cation_apfu','Fe2_cation', ...
     'Fe_cation','Fe2','Fe'}, ...
    false, 'Olivine');

end

function cpx = prepareClinopyroxeneRow(data_cpx)
% prepareClinopyroxeneRow
% Extract one-row clinopyroxene cation data without replacing NaN by zero.

if height(data_cpx) ~= 1
    error('Clinopyroxene input must be a 1-row table.');
end

cpx = struct();
[cpx.Ca, cpx.CaVariableName, cpx.CaFound] = getVarByAliases( ...
    data_cpx, {'Ca_cation_apfu','Ca_cation','Ca'}, true, 'Clinopyroxene');
[cpx.Mg, cpx.MgVariableName, cpx.MgFound] = getVarByAliases( ...
    data_cpx, {'Mg_cation_apfu','Mg_cation','Mg'}, false, 'Clinopyroxene');
[cpx.Fe, cpx.FeVariableName, cpx.FeFound] = getVarByAliases( ...
    data_cpx, ...
    {'Fe2_cation_apfu','Fe_cation_apfu','Fe2_cation', ...
     'Fe_cation','Fe2','Fe'}, ...
    false, 'Clinopyroxene');

end

function [value, matchedName, wasFound] = getVarByAliases( ...
        tbl, aliases, isRequired, mineralLabel)
% getVarByAliases
% Retrieve the first matching scalar numeric table variable. NaN is retained.
% Missing optional variables are represented by NaN, never by zero.

value = NaN;
matchedName = "missing_optional_variable";
wasFound = false;

for i = 1:numel(aliases)
    variableName = aliases{i};
    if ismember(variableName, tbl.Properties.VariableNames)
        candidateValue = tbl.(variableName);
        if ~isnumeric(candidateValue) || ~isscalar(candidateValue)
            error('%s variable %s must be a numeric scalar in a 1-row table.', ...
                mineralLabel, variableName);
        end

        value = double(candidateValue);
        matchedName = string(variableName);
        wasFound = true;
        return
    end
end

if isRequired
    error('%s table must contain one of: %s', ...
        mineralLabel, strjoin(aliases, ', '));
end

end

function sourceName = getSourceVariableName(mineralStruct, fieldName)
% getSourceVariableName
% Return the source table variable associated with a prepared structure
% field, or a descriptive placeholder for a missing optional variable.

switch fieldName
    case 'Ca'
        sourceName = mineralStruct.CaVariableName;
    case 'Mg'
        sourceName = mineralStruct.MgVariableName;
    case 'Fe'
        sourceName = mineralStruct.FeVariableName;
    otherwise
        sourceName = string(fieldName);
end

end

function XMg = calculateMgNumber(Mg, Fe)
% calculateMgNumber
% Calculate Mg/(Mg+Fe). Missing, non-finite, or zero-sum inputs return NaN.

XMg = NaN;
if isfinite(Mg) && isfinite(Fe) && (Mg + Fe) > 0
    XMg = Mg ./ (Mg + Fe);
end

end
