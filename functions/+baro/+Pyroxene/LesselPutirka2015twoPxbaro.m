function results = LesselPutirka2015twoPxbaro(rawdata_struct, T_degreeC)
% functions/+baro/+Pyroxene/LesselPutirka2015twoPxbaro.m
% Compatibility target: MATLAB R2024b
%
% Two-pyroxene barometer
% Lessel, J. and Putirka, K. (2015)
% American Mineralogist, 100, 2163-2171
% DOI: https://doi.org/10.2138/am-2015-4732
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis and one
% Clinopyroxene analysis selected from rawdata_struct.Opx and
% rawdata_struct.Cpx, respectively, and calculates pressure using Equation
% (6) of Lessel and Putirka (2015; p. 2166).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Orthopyroxene-Clinopyroxene pair,
% one output row is returned for every input temperature value.
%
% IMPORTANT: Equation (6) contains no temperature term. Therefore, for one
% selected Orthopyroxene-Clinopyroxene pair, the same calculated pressure is
% repeated for every input temperature. T_degreeC and T_K are retained in
% the output so that the common fixed-T and range-T launchers can use the
% same function interface. Input temperature is also used only for the
% optional Cpx-Opx Fe-Mg equilibrium diagnostic of Equation (15).
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% IMPORTANT: This barometer was calibrated specifically for MARTIAN
% igneous compositions and martian meteorite analog experiments. The paper
% emphasizes that martian materials are characteristically richer in FeO
% and poorer in Al2O3 than common terrestrial basaltic systems
% (pp. 2163-2164). Application to terrestrial rocks, or to compositions
% unlike the martian experimental database, is compositional extrapolation.
%
% The combined experimental database used in the paper spans broadly
% (pp. 2164-2165; Table 1):
%
%   Temperature : 950-1540 degreeC
%   Pressure    : approximately 1 atm-2.3 GPa
%                 (approximately 0.001-23 kbar)
%
% These are ranges for the COMBINED experimental database, not strict
% Equation (6)-specific calibration limits. Lessel and Putirka (2015) state
% that the two-pyroxene thermometer-barometer dataset contains 39
% experiments, of which approximately 21% were reserved for independent
% testing because the total dataset was limited (p. 2165).
%
% Equation (6) reproduces calibration pressures with R2 = 0.93 and
% RMSE = +/-0.17 GPa (n = 31), and predicts the independent test data with
% R2 = 0.96 and RMSE = +/-0.16 GPa (n = 8; pp. 2164 and 2166, Table 2 and
% Fig. 7). Results outside the experimental pressure-temperature envelope
% or outside martian pyroxene compositions should be treated as
% extrapolations.
%
% Orthopyroxene and Clinopyroxene must represent an equilibrium pair from
% the same crystallization or re-equilibration event. Exsolution lamellae,
% strong zoning, inherited grains, mixed crystal populations, and analyses
% from unrelated growth zones may yield misleading pressures.
%
% Lessel and Putirka (2015) give the following temperature-dependent
% Cpx-Opx Fe-Mg exchange relation as an equilibrium check (Eq. 15, p. 2168):
%
%   KD(Fe-Mg)Cpx-Opx = 0.115 + 7.693e-4*T(degreeC)
%
% with SEE = +/-0.09 (n = 38). The authors explicitly do NOT recommend this
% relation as a weak thermometer; it is intended as a check on an
% independently obtained temperature. This implementation reports the
% observed and predicted KD values and issues a non-stopping warning when
% their finite difference exceeds +/-0.09.
%
% The combined database includes volatile-bearing experiments (p. 2165),
% but Equation (6) uses only Orthopyroxene and Clinopyroxene compositions.
% No Liquid composition is used. Therefore, the requested exclusion of
% Liquid F and Cl from cationTotal_liq and from Liquid NaN warnings is not
% applicable to this function.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 950-1540 degreeC,
%   2) finite calculated pressure is outside approximately 0.001-23 kbar,
%   3) the observed Cpx-Opx KD differs from Equation (15) by more than
%      its reported SEE of +/-0.09,
%   4) a calculation input contains NaN,
%   5) a calculated pressure is NaN or Inf, or
%   6) a component or ratio required by Equation (6) is non-finite.
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if it makes a ratio undefined, the
% resulting NaN or Inf is retained and reported.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%   rawdata_struct.Cpx : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The following normalized cation
% variables are used by the pressure calculation:
%
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu     % total Fe expressed as FeO-equivalent cations
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   Mn_cation_apfu
%   Ti_cation_apfu
%   Cr_cation_apfu
%
% Si, Al, Fe, Mg, and Ca columns are required. Na, Mn, Ti, and Cr columns
% may be absent for backward compatibility, but an absent column is
% represented by NaN, never by zero. A present NaN value is also retained.
% Consequently, missing or NaN calculation inputs propagate to the pressure
% result rather than being silently interpreted as zero.
%
% Fe3_cation_apfu is optional and retained only as a diagnostic output. The
% Equation (6) component calculation follows the existing project
% convention using Fe_cation_apfu as total Fe. Clinopyroxene Fe3+ used for
% diagnostic site allocation is calculated from charge balance following
% Papike et al. (1974), as described by Lessel and Putirka (2015; p. 2166).
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
% Lessel and Putirka (2015), Equation (6), p. 2166:
%
%   P(GPa) =
%     -3.764
%     + 0.6739*XCaO_Cpx
%     + 33.45*XJd_Cpx
%     + 4.033*XDiHd_Cpx
%     - 5.945*XDi_Opx
%     + 3.320*(XEnFs_Cpx/XFm2Si2O6_Opx)
%     + 40.34*[Al(total)_Cpx*Al(total)_Opx]
%     - 75.80*[(Al(total)_Cpx*Na_Cpx)
%              + (Al(total)_Opx*Na_Opx)]
%
% Pyroxene components are calculated from cations normalized to six oxygen
% atoms, following Putirka (2008), as specified by Lessel and Putirka
% (2015; pp. 2164 and 2166).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = LesselPutirka2015twoPxbaro(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing Opx and Cpx tables
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Orthopyroxene-Clinopyroxene pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('LesselPutirka2015twoPxbaro requires (rawdata_struct, T_degreeC).');
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

%% 1) Retrieve pyroxene datasets
disp('=== Step 1: Preparing Orthopyroxene and Clinopyroxene datasets ===');

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end
if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

if height(rawdata_struct.Opx) < 1
    error('rawdata_struct.Opx is empty.');
end
if height(rawdata_struct.Cpx) < 1
    error('rawdata_struct.Cpx is empty.');
end

dataset_opx = rawdata_struct.Opx;
dataset_cpx = rawdata_struct.Cpx;

disp('=== Preparing Orthopyroxene and Clinopyroxene datasets has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a fixed-size cell buffer and
% concatenate once after the interactive loop. The result table is not
% repeatedly enlarged within the loop.
disp('=== Step 2: Preparing output container ===');

maxResultBlocks = max(1024, 16 * max(height(dataset_opx), height(dataset_cpx)));
resultBlocks = cell(maxResultBlocks, 1);
nResultBlocks = 0;

% Broad ranges of the combined experimental database. These are warning
% envelopes only and are not strict Equation (6)-specific limits.
experimentalT_min_degreeC = 950;
experimentalT_max_degreeC = 1540;
experimentalP_min_kbar = 0.00101325;
experimentalP_max_kbar = 23;

temperatureOutsideExperimental = isfinite(T_degreeC) & ...
    (T_degreeC < experimentalT_min_degreeC | ...
     T_degreeC > experimentalT_max_degreeC);
temperatureWarningIssued = false;
modelCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Orthopyroxene) ===');

while true
    % ----- Orthopyroxene selection -----
    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Orthopyroxene data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Orthopyroxene selected: ' char(string(selectedCode_opx))]);

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

    % ----- Input checks and calculation -----
    disp('=== Step 5: Checking calculation inputs and calculating pressure ===');

    selectedData_opx = dataset_opx(selectedIdx_opx, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    % List NaN values without changing them. Temperature is reported even
    % though it does not enter Equation (6), because it is part of the
    % common launcher interface and the equilibrium diagnostic.
    nanInputNames = findNaNInputs( ...
        selectedData_opx, selectedData_cpx, T_degreeC);

    % Reject Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(selectedData_opx, selectedData_cpx);

    row = calcPressure( ...
        selectedData_opx, selectedData_cpx, T_degreeC);

    % Repeat identifiers for all temperatures in the current calculation.
    nRows = height(row);
    row.dataCode_opx = repmat(string(selectedCode_opx), nRows, 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), nRows, 1);
    row = movevars(row, {'dataCode_opx', 'dataCode_cpx'}, 'Before', 1);

    % Store one block per selected pair without resizing the buffer.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > maxResultBlocks
        error(['Maximum number of interactive selections (%d) exceeded. ' ...
               'Restart the function to continue.'], maxResultBlocks);
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_cpx)) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_cpx)) ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Print the model-specific limitation once per function call.
    if ~modelCautionIssued
        fprintf(2, ...
            ['CAUTION: Lessel and Putirka (2015) Equation (6) was calibrated ' ...
             'specifically for martian igneous compositions (pp. 2163-2166). ' ...
             'Application to terrestrial or compositionally dissimilar systems is ' ...
             'extrapolation. Orthopyroxene and Clinopyroxene must be an equilibrium ' ...
             'pair. Equation (6) contains no temperature term, so pressure is ' ...
             'repeated unchanged across the input temperature vector.\n']);
        modelCautionIssued = true;
    end

    % Input temperature is common to all selected pairs, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideExperimental) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the 950-1540 degreeC range ' ...
             'of the combined experimental database of Lessel and Putirka ' ...
             '(2015; pp. 2164-2165). This is not a strict Equation (6)-specific ' ...
             'range, and Equation (6) itself is temperature-independent. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideExperimental), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures fall outside the broad pressure
    % envelope represented by the combined experiments.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideExperimental = finitePressure & ...
        (row.P_kbar < experimentalP_min_kbar | ...
         row.P_kbar > experimentalP_max_kbar);

    if any(pressureOutsideExperimental)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the approximately ' ...
             '1-atm-2.3-GPa (0.001-23 kbar) range of the combined experimental ' ...
             'database of Lessel and Putirka (2015; pp. 2164-2165). This is not ' ...
             'a strict Equation (6)-specific range. %d of %d finite pressure ' ...
             'point(s) are outside the range; calculated finite range = ' ...
             '%.4g-%.4g kbar for %s & %s.\n'], ...
            sum(pressureOutsideExperimental), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)));
    end

    % Compare observed Fe-Mg exchange with Equation (15). This is a
    % diagnostic equilibrium check and is not used in the pressure equation.
    printKDEquilibriumWarning(row, selectedCode_opx, selectedCode_cpx);

    % List the exact calculation-input names containing NaN. Missing
    % optional cation columns are represented as NaN and included here.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s & %s: %s.\n' ...
             '         NaN values were retained and propagated where they enter ' ...
             'Equation (6); calculated pressure may remain NaN. Because Equation ' ...
             '(6) is temperature-independent, NaN in T_degreeC does not by itself ' ...
             'force pressure to NaN.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressure is retained for diagnostics but is outside
    % the physical/useful domain and also triggers the range warning above.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & %s ' ...
             '(%d of %d points). The values were retained for diagnostic purposes.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    % Report non-finite component terms explicitly. They are not replaced
    % with zero and are allowed to propagate to pressure.
    invalidComponents = ...
        ~isfinite(row.XJd_cpx) | ...
        ~isfinite(row.XDiHd_cpx) | ...
        ~isfinite(row.XDi_opx) | ...
        ~isfinite(row.XEnFs_cpx) | ...
        ~isfinite(row.XFm2Si2O6_opx) | ...
        ~isfinite(row.ratio_EnFsCpx_Fm2Si2O6Opx) | ...
        ~isfinite(row.AlTotProduct_Eq6) | ...
        ~isfinite(row.AlNaTerm_Eq6);

    if any(invalidComponents)
        fprintf(2, ...
            ['WARNING: One or more Equation (6) component terms are non-finite ' ...
             'for %s & %s (%d of %d points). NaN or Inf values were retained ' ...
             'and propagated.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)), ...
            sum(invalidComponents), ...
            numel(invalidComponents));
    end

    % A negative calculated Fe2+ after the Cpx charge-balance estimate is a
    % diagnostic inconsistency. It is retained and does not stop Eq. (6).
    negativeCalculatedFe2 = isfinite(row.Fe2_calc_cpx) & row.Fe2_calc_cpx < 0;
    if any(negativeCalculatedFe2)
        fprintf(2, ...
            ['WARNING: Charge-balance calculation produced negative Fe2+ in ' ...
             'Clinopyroxene for %s & %s. The diagnostic value was retained; ' ...
             'check the cation normalization and mineral analysis.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_cpx)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'LesselPutirka2015twoPxbaro', ...
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
function nanInputNames = findNaNInputs(data_opx, data_cpx, T_degreeC)
% findNaNInputs
% Return names of Equation (6) inputs containing NaN. Missing optional
% cation columns are represented as NaN and reported. Values are not
% replaced by zero.

maxNames = 32;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu'};

optionalCalculationVariables = { ...
    'Na_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

mineralTables = {data_opx, data_cpx};
mineralLabels = {'Opx', 'Cpx'};

for iMineral = 1:numel(mineralTables)
    tbl = mineralTables{iMineral};
    label = mineralLabels{iMineral};

    for i = 1:numel(requiredVariables)
        variableName = requiredVariables{i};
        if ismember(variableName, tbl.Properties.VariableNames)
            value = toScalarDouble(tbl.(variableName));
            if isnan(value)
                nNanInputs = nNanInputs + 1;
                nanInputBuffer(nNanInputs) = ...
                    string(label) + "." + string(variableName);
            end
        end
    end

    for i = 1:numel(optionalCalculationVariables)
        variableName = optionalCalculationVariables{i};
        if ismember(variableName, tbl.Properties.VariableNames)
            value = toScalarDouble(tbl.(variableName));
            if isnan(value)
                nNanInputs = nNanInputs + 1;
                nanInputBuffer(nNanInputs) = ...
                    string(label) + "." + string(variableName);
            end
        else
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = ...
                string(label) + "." + string(variableName) + "(missing->NaN)";
        end
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_opx, data_cpx)
% validateNonNegativeInputs
% Reject Inf and finite negative values in direct cation inputs used by
% Equation (6). Zero and NaN are intentionally allowed and retained.

maxNames = 32;
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

calculationVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

mineralTables = {data_opx, data_cpx};
mineralLabels = {'Opx', 'Cpx'};

for iMineral = 1:numel(mineralTables)
    tbl = mineralTables{iMineral};
    label = mineralLabels{iMineral};

    for i = 1:numel(calculationVariables)
        variableName = calculationVariables{i};
        if ~ismember(variableName, tbl.Properties.VariableNames)
            continue
        end

        value = toScalarDouble(tbl.(variableName));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = ...
                string(label) + "." + string(variableName);
        end
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['LesselPutirka2015twoPxbaro: calculation inputs must be ' ...
           'non-negative. NaN is allowed, but Inf and finite negative ' ...
           'value(s) are prohibited. Invalid input(s): ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_opx, data_cpx, T_degreeC)
% calcPressure
% Compute Equation (6) pressure for one Opx-Cpx pair at one or more input
% temperatures. Equation (6) is temperature-independent, so the same
% pressure is repeated for all temperature rows. NaN values propagate.

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Extract one-row pyroxene cation data. Missing optional calculation
% variables are represented as NaN, not zero.
opx = preparePyroxeneRow(data_opx, 'Opx');
cpx = preparePyroxeneRow(data_cpx, 'Cpx');

% Calculate component terms. Helper functions preserve NaN values while
% retaining the non-negative component allocation used in the original
% project implementation.
opxTerms = calcOpxTerms(opx);
cpxTerms = calcCpxTerms(cpx);

ratio_EnFsCpx_Fm2Si2O6Opx_scalar = ...
    cpxTerms.XEnFs ./ opxTerms.XFm2Si2O6;

AlTotProduct_Eq6_scalar = cpx.Al .* opx.Al;
AlNaTerm_Eq6_scalar = ...
    (cpx.Al .* cpx.Na) + ...
    (opx.Al .* opx.Na);

% Pressure calculation. No finite-value gate is used: NaN and Inf component
% terms remain non-finite and propagate to pressure.
PEq6_GPa_scalar = ...
    -3.764 ...
    + 0.6739 .* cpx.Ca ...
    + 33.45 .* cpxTerms.XJd ...
    + 4.033 .* cpxTerms.XDiHd ...
    - 5.945 .* opxTerms.XDi ...
    + 3.320 .* ratio_EnFsCpx_Fm2Si2O6Opx_scalar ...
    + 40.34 .* AlTotProduct_Eq6_scalar ...
    - 75.80 .* AlNaTerm_Eq6_scalar;

PEq6_kbar_scalar = PEq6_GPa_scalar .* 10;

% Fe-Mg exchange diagnostic of Equation (15). Fe_cation_apfu is used as
% FeO-equivalent total Fe, matching the input-table convention.
KD_FeMg_CpxOpx_observed_scalar = ...
    (cpx.FeT ./ cpx.Mg) ./ (opx.FeT ./ opx.Mg);
KD_FeMg_CpxOpx_predicted = 0.115 + 7.693e-4 .* T_degreeC;
KD_FeMg_CpxOpx_difference = ...
    KD_FeMg_CpxOpx_observed_scalar - KD_FeMg_CpxOpx_predicted;

% Expand composition-dependent scalar values to the temperature-vector
% length. All output columns therefore have a consistent number of rows.
Si_opx = repmat(opx.Si, nT, 1);
Al_opx = repmat(opx.Al, nT, 1);
FeT_opx = repmat(opx.FeT, nT, 1);
Fe2_input_opx = repmat(opx.Fe2_input, nT, 1);
Fe3_input_opx = repmat(opx.Fe3_input, nT, 1);
Mg_opx = repmat(opx.Mg, nT, 1);
Ca_opx = repmat(opx.Ca, nT, 1);
Na_opx = repmat(opx.Na, nT, 1);
Mn_opx = repmat(opx.Mn, nT, 1);
Ti_opx = repmat(opx.Ti, nT, 1);
Cr_opx = repmat(opx.Cr, nT, 1);

Si_cpx = repmat(cpx.Si, nT, 1);
Al_cpx = repmat(cpx.Al, nT, 1);
FeT_cpx = repmat(cpx.FeT, nT, 1);
Fe2_input_cpx = repmat(cpx.Fe2_input, nT, 1);
Fe3_input_cpx = repmat(cpx.Fe3_input, nT, 1);
Mg_cpx = repmat(cpx.Mg, nT, 1);
Ca_cpx = repmat(cpx.Ca, nT, 1);
Na_cpx = repmat(cpx.Na, nT, 1);
Mn_cpx = repmat(cpx.Mn, nT, 1);
Ti_cpx = repmat(cpx.Ti, nT, 1);
Cr_cpx = repmat(cpx.Cr, nT, 1);

AlIV_opx = repmat(opxTerms.AlIV, nT, 1);
AlVI_opx = repmat(opxTerms.AlVI, nT, 1);
XDi_opx = repmat(opxTerms.XDi, nT, 1);
XFm2Si2O6_opx = repmat(opxTerms.XFm2Si2O6, nT, 1);
XFmAl2SiO6_opx = repmat(opxTerms.XFmAl2SiO6, nT, 1);

AlIV_cpx = repmat(cpxTerms.AlIV, nT, 1);
AlVI_cpx = repmat(cpxTerms.AlVI, nT, 1);
Fe3_calc_cpx = repmat(cpxTerms.Fe3, nT, 1);
Fe2_calc_cpx = repmat(cpxTerms.Fe2, nT, 1);
XJd_cpx = repmat(cpxTerms.XJd, nT, 1);
XCaTs_cpx = repmat(cpxTerms.XCaTs, nT, 1);
XCaTi_cpx = repmat(cpxTerms.XCaTi, nT, 1);
XCrCaTs_cpx = repmat(cpxTerms.XCrCaTs, nT, 1);
XDiHd_cpx = repmat(cpxTerms.XDiHd, nT, 1);
XEnFs_cpx = repmat(cpxTerms.XEnFs, nT, 1);
XCaO_cpx = repmat(cpx.Ca, nT, 1);

ratio_EnFsCpx_Fm2Si2O6Opx = ...
    repmat(ratio_EnFsCpx_Fm2Si2O6Opx_scalar, nT, 1);
AlTotProduct_Eq6 = repmat(AlTotProduct_Eq6_scalar, nT, 1);
AlNaTerm_Eq6 = repmat(AlNaTerm_Eq6_scalar, nT, 1);

PEq6_GPa = repmat(PEq6_GPa_scalar, nT, 1);
PEq6_kbar = repmat(PEq6_kbar_scalar, nT, 1);
KD_FeMg_CpxOpx_observed = ...
    repmat(KD_FeMg_CpxOpx_observed_scalar, nT, 1);

% Applicability and diagnostic flags.
isWithinExperimentalTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 950 & T_degreeC <= 1540;
isWithinExperimentalPRange = ...
    isfinite(PEq6_kbar) & PEq6_kbar >= 0.00101325 & PEq6_kbar <= 23;
isWithinKDEquilibriumSEE = ...
    isfinite(KD_FeMg_CpxOpx_difference) & ...
    abs(KD_FeMg_CpxOpx_difference) <= 0.09;
isEquationFinite = ...
    isfinite(PEq6_GPa) & ...
    isfinite(ratio_EnFsCpx_Fm2Si2O6Opx);

% Pack outputs using pre-sized vectors of equal height.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.Si_opx = Si_opx;
row.Al_opx = Al_opx;
row.Fe2_opx = Fe2_input_opx;
row.Fe3_opx = Fe3_input_opx;
row.FeTot_opx = FeT_opx;
row.Fe3_input_opx = Fe3_input_opx;
row.Mg_opx = Mg_opx;
row.Ca_opx = Ca_opx;
row.Na_opx = Na_opx;
row.Mn_opx = Mn_opx;
row.Ti_opx = Ti_opx;
row.Cr_opx = Cr_opx;

row.Si_cpx = Si_cpx;
row.Al_cpx = Al_cpx;
row.Fe2_cpx = Fe2_input_cpx;
row.Fe3_cpx = Fe3_input_cpx;
row.FeTot_cpx = FeT_cpx;
row.Fe3_input_cpx = Fe3_input_cpx;
row.Mg_cpx = Mg_cpx;
row.Ca_cpx = Ca_cpx;
row.Na_cpx = Na_cpx;
row.Mn_cpx = Mn_cpx;
row.Ti_cpx = Ti_cpx;
row.Cr_cpx = Cr_cpx;

row.AlIV_opx = AlIV_opx;
row.AlVI_opx = AlVI_opx;
row.XDi_opx = XDi_opx;
row.XFm2Si2O6_opx = XFm2Si2O6_opx;
row.XFmAl2SiO6_opx = XFmAl2SiO6_opx;

row.AlIV_cpx = AlIV_cpx;
row.AlVI_cpx = AlVI_cpx;
row.Fe3_calc_cpx = Fe3_calc_cpx;
row.Fe2_calc_cpx = Fe2_calc_cpx;
row.XJd_cpx = XJd_cpx;
row.XCaTs_cpx = XCaTs_cpx;
row.XCaTi_cpx = XCaTi_cpx;
row.XCrCaTs_cpx = XCrCaTs_cpx;
row.XDiHd_cpx = XDiHd_cpx;
row.XEnFs_cpx = XEnFs_cpx;
row.XCaO_cpx = XCaO_cpx;

row.ratio_EnFsCpx_Fm2Si2O6Opx = ...
    ratio_EnFsCpx_Fm2Si2O6Opx;
row.AlTotProduct_Eq6 = AlTotProduct_Eq6;
row.AlNaTerm_Eq6 = AlNaTerm_Eq6;

row.KD_FeMg_CpxOpx_observed = KD_FeMg_CpxOpx_observed;
row.KD_FeMg_CpxOpx_predicted_Eq15 = KD_FeMg_CpxOpx_predicted;
row.KD_FeMg_CpxOpx_difference = KD_FeMg_CpxOpx_difference;
row.KD_FeMg_CpxOpx_SEE = repmat(0.09, nT, 1);

row.PEq6_GPa = PEq6_GPa;
row.PEq6_kbar = PEq6_kbar;

% Common launcher/output aliases.
row.P_GPa = PEq6_GPa;
row.P_kbar = PEq6_kbar;

row.P_RMSE_calibration_GPa = repmat(0.17, nT, 1);
row.P_RMSE_test_GPa = repmat(0.16, nT, 1);
row.P_RMSE_calibration_kbar = repmat(1.7, nT, 1);
row.P_RMSE_test_kbar = repmat(1.6, nT, 1);

row.isWithinExperimentalTRange = isWithinExperimentalTRange;
row.isWithinExperimentalPRange = isWithinExperimentalPRange;
row.isWithinKDEquilibriumSEE = isWithinKDEquilibriumSEE;
row.isEquationFinite = isEquationFinite;

end

function px = preparePyroxeneRow(data_px, mineralLabel)
% preparePyroxeneRow
% Extract one-row normalized pyroxene cation data. Required variables must
% exist. Missing optional calculation variables are represented by NaN and
% are never replaced by zero.

if height(data_px) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

px = struct();
px.Si = getVarOrError(data_px, 'Si_cation_apfu', mineralLabel);
px.Al = getVarOrError(data_px, 'Al_cation_apfu', mineralLabel);
px.FeT = getVarOrError(data_px, 'Fe_cation_apfu', mineralLabel);
px.Mg = getVarOrError(data_px, 'Mg_cation_apfu', mineralLabel);
px.Ca = getVarOrError(data_px, 'Ca_cation_apfu', mineralLabel);

px.Na = getVarOrNaN(data_px, 'Na_cation_apfu');
px.Mn = getVarOrNaN(data_px, 'Mn_cation_apfu');
px.Ti = getVarOrNaN(data_px, 'Ti_cation_apfu');
px.Cr = getVarOrNaN(data_px, 'Cr_cation_apfu');
px.Fe3_input = getVarOrNaN(data_px, 'Fe3_cation_apfu');

% Reject Inf and finite negative values in extracted raw cation values.
fieldNames = fieldnames(px);
for i = 1:numel(fieldNames)
    value = px.(fieldNames{i});
    if ~isscalar(value)
        error('%s variable %s must be scalar in a 1-row table.', ...
            mineralLabel, fieldNames{i});
    end
    if isinf(value) || (isfinite(value) && value < 0)
        error('%s contains an invalid negative or Inf value for %s.', ...
            mineralLabel, fieldNames{i});
    end
end

if isfinite(px.Fe3_input) && isfinite(px.FeT) && px.Fe3_input > px.FeT
    error('%s contains Fe3_cation_apfu > Fe_cation_apfu.', mineralLabel);
end

% Backward-compatible Fe2 output. When an explicit ferric-iron value is
% present, subtract it from total Fe; otherwise retain total Fe as the
% FeO-equivalent ferrous value commonly used for EPMA input. This value is
% diagnostic only and does not replace total Fe in Equation (6).
if isfinite(px.Fe3_input)
    px.Fe2_input = px.FeT - px.Fe3_input;
else
    px.Fe2_input = px.FeT;
end

end

function terms = calcOpxTerms(opx)
% calcOpxTerms
% Calculate the Orthopyroxene terms required by Equation (6). NaN values
% remain NaN. Finite negative component allocations are limited to zero,
% following the original project implementation.

terms = struct();
terms.AlIV = lowerBoundZero(2 - opx.Si);
terms.AlVI = lowerBoundZero(opx.Al - terms.AlIV);

% Di component in Opx is approximated by Ca cations.
terms.XDi = opx.Ca;

% Ferromagnesian quadrilateral component in Opx.
terms.XFm2Si2O6 = lowerBoundZero(opx.FeT + opx.Mg + opx.Mn);

% Retained for diagnostics and compatibility with the thermometer version.
terms.XFmAl2SiO6 = lowerBoundZero(min(terms.AlIV, terms.AlVI));

end

function terms = calcCpxTerms(cpx)
% calcCpxTerms
% Calculate Clinopyroxene components following the normative scheme used by
% the existing project and described by Lessel and Putirka (2015; p. 2166).
% NaN values remain NaN.

terms = struct();
terms.AlIV = lowerBoundZero(2 - cpx.Si);
terms.AlVI = lowerBoundZero(cpx.Al - terms.AlIV);

% Papike et al. (1974) charge-balance estimate of Cpx Fe3+.
terms.Fe3 = lowerBoundZero( ...
    terms.AlIV + cpx.Na - terms.AlVI - cpx.Cr - 2 .* cpx.Ti);
terms.Fe2 = cpx.FeT - terms.Fe3;

terms.XJd = lowerBoundZero(min(terms.AlVI, cpx.Na));
terms.XCaTs = lowerBoundZero(terms.AlVI - terms.XJd);
terms.XCaTi = lowerBoundZero((terms.AlIV - terms.XCaTs) ./ 2);
terms.XCrCaTs = lowerBoundZero(cpx.Cr ./ 2);
terms.XDiHd = lowerBoundZero( ...
    cpx.Ca - terms.XCaTi - terms.XCaTs - terms.XCrCaTs);

% En-Fs component in Cpx, following the existing project convention.
terms.XEnFs = lowerBoundZero( ...
    (cpx.FeT + cpx.Mg - terms.XDiHd) ./ 2);

end

function value = lowerBoundZero(value)
% lowerBoundZero
% Apply a zero lower bound to finite values while preserving NaN and Inf.

if isfinite(value)
    value = max(value, 0);
end

end

function printKDEquilibriumWarning(row, selectedCode_opx, selectedCode_cpx)
% printKDEquilibriumWarning
% Compare observed Cpx-Opx Fe-Mg exchange with Lessel and Putirka (2015)
% Equation (15). The relation is an equilibrium diagnostic, not a
% thermometer and not part of Equation (6).

finiteComparison = ...
    isfinite(row.KD_FeMg_CpxOpx_observed) & ...
    isfinite(row.KD_FeMg_CpxOpx_predicted_Eq15) & ...
    isfinite(row.KD_FeMg_CpxOpx_difference);

outsideSEE = finiteComparison & ...
    abs(row.KD_FeMg_CpxOpx_difference) > row.KD_FeMg_CpxOpx_SEE;

if any(outsideSEE)
    finiteObserved = row.KD_FeMg_CpxOpx_observed(finiteComparison);
    finitePredicted = row.KD_FeMg_CpxOpx_predicted_Eq15(finiteComparison);
    finiteDifference = row.KD_FeMg_CpxOpx_difference(finiteComparison);

    fprintf(2, ...
        ['WARNING: Cpx-Opx Fe-Mg exchange differs from the temperature-dependent ' ...
         'equilibrium relation of Lessel and Putirka (2015, Eq. 15, p. 2168) ' ...
         'by more than its SEE of +/-0.09 for %s & %s. %d of %d finite ' ...
         'temperature point(s) are outside the SEE; observed KD = %.4g; ' ...
         'predicted KD range = %.4g-%.4g; difference range = %.4g-%.4g. ' ...
         'Equation (15) is only an equilibrium check and is not recommended ' ...
         'as a thermometer.\n'], ...
        char(string(selectedCode_opx)), ...
        char(string(selectedCode_cpx)), ...
        sum(outsideSEE), ...
        sum(finiteComparison), ...
        finiteObserved(1), ...
        min(finitePredicted), ...
        max(finitePredicted), ...
        min(finiteDifference), ...
        max(finiteDifference));
end

if ~any(finiteComparison)
    fprintf(2, ...
        ['WARNING: The Cpx-Opx Fe-Mg equilibrium diagnostic could not be ' ...
         'evaluated for %s & %s because KD or temperature is non-finite. ' ...
         'Diagnostic values were retained in the output.\n'], ...
        char(string(selectedCode_opx)), ...
        char(string(selectedCode_cpx)));
end

end

function value = getVarOrError(tbl, variableName, mineralLabel)
% getVarOrError
% Retrieve a required scalar table variable without altering NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = toScalarDouble(tbl.(variableName));

end

function value = getVarOrNaN(tbl, variableName)
% getVarOrNaN
% Retrieve an optional scalar variable. Missing variables are represented
% by NaN, never by zero.

if ismember(variableName, tbl.Properties.VariableNames)
    value = toScalarDouble(tbl.(variableName));
else
    value = NaN;
end

end

function value = toScalarDouble(raw)
% toScalarDouble
% Convert a one-row table value to one double scalar while preserving NaN.

if isempty(raw)
    value = NaN;
    return
end

if istable(raw)
    raw = raw{1, 1};
elseif iscell(raw)
    raw = raw{1};
elseif numel(raw) > 1
    raw = raw(1);
end

if isnumeric(raw) || islogical(raw)
    value = double(raw);
elseif isstring(raw)
    if ismissing(raw)
        value = NaN;
    else
        value = str2double(raw);
    end
elseif ischar(raw)
    value = str2double(string(raw));
else
    error('Unsupported table value type: %s.', class(raw));
end

if ~isscalar(value)
    error('Table value could not be converted to a scalar double.');
end

end
