function results = JohnsonRutherford1988(rawdata_struct, T_degreeC)
% functions/+baro/+Amphibole/JohnsonRutherford1988.m
% Tested with MATLAB R2024b
%
% Experimentally calibrated aluminum-in-hornblende geobarometer
% Johnson, M.C. and Rutherford, M.J. (1989)
% Geology, 17, 837-841
% DOI: https://doi.org/10.1130/0091-7613(1989)017<0837:ECOTAI>2.3.CO;2
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Amphibole analysis and calculates
% pressure using the experimentally calibrated total-Al-in-hornblende
% geobarometer of Johnson and Rutherford (1989).
%
% The function name is retained as JohnsonRutherford1988 for compatibility
% with the existing thermobaroMin launcher structure and with the earlier
% EOS abstract citation. The pressure equation implemented here is the final
% regression published by Johnson and Rutherford (1989), not the preliminary
% coefficients commonly attributed to the 1988 EOS abstract.
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Amphibole analysis, one output row
% is returned for every input temperature value.
%
% Temperature does not enter the published pressure equation. T_degreeC is
% retained in the output and used only for experimental-range screening.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Johnson and Rutherford (1989) experimentally calibrated:
%
%   P(kbar) = -3.46 + 4.23 * AlT
%
% where AlT is total Al in hornblende calculated on a 23-oxygen basis.
% The equation, coefficient uncertainties, r^2 = 0.99, and regression
% standard deviation of 0.22 kbar are given as Eq. (3) on p. 838 and
% discussed on p. 839.
%
% Experimental calibration:
%   Pressure    : 2-8 kbar
%   Temperature : mainly 740-780 degreeC
%   Materials   : natural volcanic and plutonic calc-alkaline compositions
%
% The abstract on p. 837 states that the required assemblage was
% equilibrated over 2-8 kbar at 740-780 degreeC. The Experimental techniques
% section on pp. 837-838 describes the full experimental temperature span
% as approximately 720-780 degreeC, whereas the main Fish Canyon Tuff
% experiments used for the pressure calibration were run at 740-780 degreeC
% (Results, p. 838). This implementation uses 740-780 degreeC as the direct
% calibration-warning range and documents 720-780 degreeC as the broader
% experimental envelope.
%
% Required equilibrium assemblage:
%   melt + fluid + quartz + sanidine + plagioclase (approximately An30) +
%   hornblende + biotite + sphene + magnetite or ilmenite
%
% The required phase assemblage is stated in the abstract and Introduction
% on p. 837 and described in the Results on p. 838. Melt must be present to
% reduce the thermodynamic degrees of freedom and to permit hornblende
% recrystallization and equilibration on experimental time scales (p. 838).
%
% Hornblende selection:
%   - Use euhedral, unaltered hornblende RIMS in contact and equilibrium
%     with melt or glass.
%   - The calibration used only euhedral hornblende rims in contact with
%     melt (Analytical techniques, p. 838).
%   - Cores, altered rims, oxide-breakdown rims, inherited crystals,
%     subsolidus amphiboles, and analyses not demonstrably equilibrated with
%     the final melt should not be used.
%
% Quartz saturation is critical:
%   Hornblende not equilibrated with quartz contains anomalously high AlT
%   and yields pressure estimates that are too high. In a 5 kbar experiment
%   lacking quartz, the calibration returned approximately 6.2 kbar, an
%   overestimate of about 1.2 kbar (p. 839).
%
% Sanidine/K-feldspar and the other required phases are also important
% because their coexistence buffers the relevant component activities.
% Sphene absence was accepted only for specific Long Valley samples in
% which Ti activity was considered buffered by two Fe-Ti oxides and the
% hornblende TiO2 contents matched the calibration dataset (p. 840). This is
% a specific exception, not a general statement that sphene is unnecessary.
%
% Low-pressure caution:
%   - Below approximately 2 kbar, hornblende Al may become sensitive to
%     temperature as well as pressure (p. 841).
%   - Igneous hornblende is generally unstable below approximately 1.5 kbar
%     unless sufficiently hydrous melt conditions are maintained (p. 841).
%   - A 1.5 kbar, 720 degreeC Long Valley experiment supported application
%     to that specific composition, but it does not extend the general
%     experimental calibration below 2 kbar.
%
% High-pressure caution:
%   The published calibration is 2-8 kbar. A cited 10 kbar experiment was
%   not reversed and may not have reached equilibrium; its significance was
%   considered uncertain (p. 840). Results above 8 kbar are extrapolations.
%
% Uncertainty:
%   - Regression standard deviation: 0.22 kbar (p. 838).
%   - Approximate application uncertainty: +/-0.5 kbar at lower pressures
%     and up to approximately +/-0.6 kbar near 8 kbar (p. 839).
%
% Experimental hornblende-composition envelope:
%   Table 2 on p. 838 spans approximately:
%     AlT : 1.27-2.74 apfu
%     Si  : 6.07-6.90 apfu
%     Ca  : 1.75-1.93 apfu
%   These are descriptive experimental-data limits, not formal universal
%   compositional filters.
%
% The supplementary Si < 7.5 and Ca > 1.6 apfu screen is retained for
% backward compatibility with earlier total-Al hornblende implementations,
% but it is not sufficient to demonstrate applicability of this calibration.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 740-780 degreeC,
%   2) finite calculated pressure is outside 2-8 kbar,
%   3) calculated pressure is below 2 kbar or below 1.5 kbar,
%   4) calculated pressure is above 8 kbar,
%   5) AlT is outside the approximate 1.27-2.74 apfu experimental envelope,
%   6) the supplementary Si-Ca screen is not satisfied,
%   7) a calculation or screening input contains NaN,
%   8) calculated pressure is NaN or Inf, or
%   9) negative finite pressure is calculated.
%
% The program cannot verify hornblende rim position, alteration, melt
% contact, quartz saturation, sanidine presence, phase equilibrium, or the
% complete required assemblage. These conditions must be assessed
% petrographically by the user.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include:
%
%   Required Amphibole variables:
%     Si_cation_apfu
%     Ti_cation_apfu
%     Al_cation_apfu          % AlT used by the pressure equation
%     Fe_cation_apfu          % total Fe in the input cation table
%     Mg_cation_apfu
%     Ca_cation_apfu
%     Na_cation_apfu
%
%   Optional Amphibole variables retained in the output when present:
%     Fe3_cation_apfu
%     Mn_cation_apfu
%     K_cation_apfu
%     Cr_cation_apfu
%
% Finite cation values must be non-negative. NaN is allowed and retained as
% missing; it is never replaced by zero. Inf and finite negative values are
% prohibited. A NaN Al_cation_apfu value propagates directly to P_kbar.
%
% This barometer does not use a liquid (Liq) composition. Therefore, rules
% concerning exclusion of F and Cl from cationTotal_liq and their exclusion
% from Liq NaN warnings are not applicable to this function.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
% 1) Amphibole Al-site quantities retained for output:
%     AlT  = total Al cations per 23 oxygens
%     AlIV = max(0, 8 - Si)
%     AlVI = max(0, AlT - AlIV)
%
% 2) Pressure, Johnson and Rutherford (1989), Eq. (3):
%     P(kbar) = 4.23 * AlT - 3.46
%
% Notes:
% - This final 1989 regression replaces the preliminary coefficients
%   P = 4.28*AlT - 3.54 commonly associated with the 1988 EOS abstract.
% - AlT is the only mineral-composition variable in the pressure equation.
% - T_degreeC does not enter the pressure equation.
% - Negative calculated pressure is retained for diagnostic purposes and
%   reported by a non-stopping warning.
% - NaN is retained and propagated rather than interpreted as zero.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = JohnsonRutherford1988(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing an Amphibole table
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per input temperature value for every
%             user-selected Amphibole analysis
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('JohnsonRutherford1988 requires (rawdata_struct, T_degreeC).');
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

%% 1) Retrieve cation dataset
% Extract the required Amphibole table. The source table is not modified.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Amphibole') || ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end

dataset_amp = rawdata_struct.Amphibole;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-analysis result in a cell buffer and concatenate once
% after the interactive loop. This avoids repeated growth of the full output
% table during every calculation.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct experimental calibration ranges.
calibrationT_min_degreeC = 740;
calibrationT_max_degreeC = 780;
calibrationP_min_kbar = 2;
calibrationP_max_kbar = 8;

% Approximate experimental hornblende-composition envelope from Table 2.
experimentalAlT_min_apfu = 1.27;
experimentalAlT_max_apfu = 2.74;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);

temperatureWarningIssued = false;
applicationCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % ----- Amphibole selection -----
    % The first table column is used only as the displayed identifier.
    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, ok] = listdlg( ...
        'PromptString', 'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_amp)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the pressure ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);

    % Reject Inf and finite negative cation values. NaN is allowed and is
    % retained for propagation through the calculation.
    validateNonNegativeInputs(selectedData_amp);

    % Identify NaN values in the pressure equation, temperature input, and
    % numerical screening variables.
    nanInputNames = findNaNInputs(selectedData_amp, T_degreeC);

    row = calcPressure(selectedData_amp, T_degreeC);

    % Repeat the selected identifier for every temperature row.
    row.dataCode_amphibole = ...
        repmat(string(selectedCode_amp), height(row), 1);
    row = movevars(row, 'dataCode_amphibole', 'Before', 1);

    % Store one result block per selected analysis. Expand the buffer only
    % when its current capacity is exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated result for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    finitePressure = isfinite(row.P_kbar);
    finiteTemperature = isfinite(row.T_degreeC);

    if height(row) == 1
        if finitePressure
            fprintf('%s: P = %.6g kbar at T = %.6g degreeC\n', ...
                char(string(selectedCode_amp)), ...
                row.P_kbar, row.T_degreeC);
        else
            fprintf('%s: P = %s kbar at T = %.6g degreeC\n', ...
                char(string(selectedCode_amp)), ...
                char(string(row.P_kbar)), row.T_degreeC);
        end
    else
        if any(finitePressure)
            finitePressureValues = row.P_kbar(finitePressure);
            if any(finiteTemperature)
                finiteTemperatureValues = row.T_degreeC(finiteTemperature);
                fprintf('%s: P = %.6g to %.6g kbar for finite T = %.6g to %.6g degreeC\n', ...
                    char(string(selectedCode_amp)), ...
                    min(finitePressureValues), max(finitePressureValues), ...
                    min(finiteTemperatureValues), max(finiteTemperatureValues));
            else
                fprintf('%s: P = %.6g to %.6g kbar; all %d input temperature values are NaN\n', ...
                    char(string(selectedCode_amp)), ...
                    min(finitePressureValues), max(finitePressureValues), ...
                    height(row));
            end
        else
            if any(finiteTemperature)
                finiteTemperatureValues = row.T_degreeC(finiteTemperature);
                fprintf('%s: P = NaN/Inf for all %d temperature point(s), finite T = %.6g to %.6g degreeC\n', ...
                    char(string(selectedCode_amp)), height(row), ...
                    min(finiteTemperatureValues), max(finiteTemperatureValues));
            else
                fprintf('%s: P = NaN/Inf for all %d point(s); all input temperature values are NaN\n', ...
                    char(string(selectedCode_amp)), height(row));
            end
        end
    end

    % Print the major non-numerical application limitations once per call.
    if ~applicationCautionIssued
        fprintf(2, ...
            ['CAUTION: Johnson and Rutherford (1989) calibrated this equation using ' ...
             'euhedral hornblende rims in contact and equilibrium with melt in the ' ...
             'assemblage melt + fluid + quartz + sanidine + plagioclase + hornblende + ' ...
             'biotite + sphene + magnetite or ilmenite (pp. 837-839). The program ' ...
             'cannot verify rim position, alteration, melt contact, quartz saturation, ' ...
             'sanidine presence, phase coexistence, or equilibrium.\n']);
        applicationCautionIssued = true;
    end

    % Input temperature is common to all selected analyses in the current
    % function call, so this warning is printed only once.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperatureValues = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the direct Johnson and Rutherford ' ...
             '(1989) calibration range of 740-780 degreeC (abstract, p. 837; Results, ' ...
             'p. 838). %d of %d finite temperature point(s) are outside the range; ' ...
             'input finite range = %.6g-%.6g degreeC. The broader experimental ' ...
             'temperature envelope was approximately 720-780 degreeC (pp. 837-838).\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperatureValues), ...
            min(finiteTemperatureValues), ...
            max(finiteTemperatureValues));
        temperatureWarningIssued = true;
    end

    % Warn when finite pressure results fall outside the published 2-8 kbar
    % experimental calibration range.
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the Johnson and Rutherford ' ...
             '(1989) experimental calibration range of 2-8 kbar (pp. 837-838). ' ...
             '%d of %d finite pressure point(s) are outside the range; ' ...
             'calculated finite range = %.6g-%.6g kbar for %s. ' ...
             'The values are retained in the output table.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_amp)));
    end

    % Below 2 kbar, temperature effects may become important.
    pressureBelowTwo = finitePressure & row.P_kbar < 2;
    if any(pressureBelowTwo)
        fprintf(2, ...
            ['WARNING: Calculated pressure is below 2 kbar for %s (%d of %d finite ' ...
             'point(s)). Johnson and Rutherford (1989) note that hornblende Al may ' ...
             'become sensitive to temperature as well as pressure below approximately ' ...
             '2 kbar (p. 841).\n'], ...
            char(string(selectedCode_amp)), ...
            sum(pressureBelowTwo), ...
            sum(finitePressure));
    end

    % Below approximately 1.5 kbar, igneous hornblende stability is limited.
    pressureBelowOnePointFive = finitePressure & row.P_kbar < 1.5;
    if any(pressureBelowOnePointFive)
        fprintf(2, ...
            ['WARNING: Calculated pressure is below approximately 1.5 kbar for %s ' ...
             '(%d of %d finite point(s)). Igneous hornblende is generally unstable ' ...
             'below about 1.5 kbar unless sufficiently hydrous melt conditions are ' ...
             'maintained. The 1.5 kbar test in the paper was composition-specific ' ...
             '(p. 841).\n'], ...
            char(string(selectedCode_amp)), ...
            sum(pressureBelowOnePointFive), ...
            sum(finitePressure));
    end

    % Above 8 kbar, the result is an extrapolation.
    pressureAboveEight = finitePressure & row.P_kbar > 8;
    if any(pressureAboveEight)
        fprintf(2, ...
            ['WARNING: Calculated pressure exceeds 8 kbar for %s (%d of %d finite ' ...
             'point(s)). This is an extrapolation beyond the experimental calibration. ' ...
             'A cited 10 kbar experiment was not reversed and may not have reached ' ...
             'equilibrium (p. 840).\n'], ...
            char(string(selectedCode_amp)), ...
            sum(pressureAboveEight), ...
            sum(finitePressure));
    end

    % Warn when AlT is outside the approximate composition range represented
    % by the experimental hornblendes in Table 2.
    finiteAlT = isfinite(row.AlT_amp(1));
    if finiteAlT && ...
            (row.AlT_amp(1) < experimentalAlT_min_apfu || ...
             row.AlT_amp(1) > experimentalAlT_max_apfu)
        fprintf(2, ...
            ['WARNING: Amphibole AlT = %.6g apfu for %s is outside the approximate ' ...
             'experimental hornblende range of 1.27-2.74 apfu represented in Table 2 ' ...
             'of Johnson and Rutherford (1989; p. 838). This is a descriptive data ' ...
             'envelope, not a formal universal compositional limit.\n'], ...
            row.AlT_amp(1), ...
            char(string(selectedCode_amp)));
    end

    % Retain the original Si-Ca screen as a supplementary diagnostic.
    if ~row.isIgneousChemicalScreen(1)
        fprintf(2, ...
            ['WARNING: Amphibole %s does not satisfy the supplementary igneous ' ...
             'calcic-amphibole screen retained from the original implementation: ' ...
             'Si < 7.5 and Ca > 1.6 apfu. Si = %.6g, Ca = %.6g. The calculation ' ...
             'was continued. This screen does not substitute for the required ' ...
             'petrographic and phase-equilibrium conditions.\n'], ...
            char(string(selectedCode_amp)), ...
            row.Si_amp(1), row.Ca_amp(1));
    end

    % List exact required input names containing NaN. Temperature-vector NaN
    % positions are reported by index.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s: %s.\n' ...
             '         NaN values were retained and were not replaced by zero; ' ...
             'the pressure and/or applicability fields may remain NaN or false.\n'], ...
            char(string(selectedCode_amp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all NaN/Inf pressure results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_amp)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressure is retained for diagnostic purposes.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s (%d of %d ' ...
             'points). The value was retained for diagnostic purposes and is outside ' ...
             'the 2-8 kbar calibration range.\n'], ...
            char(string(selectedCode_amp)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another amphibole analysis.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'JohnsonRutherford1988', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered blocks once. Return an empty table when the user
% canceled before completing any calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_amphibole, T_degreeC)
% findNaNInputs
% Return names of NaN inputs used by the pressure equation, temperature
% screening, or supplementary numerical applicability checks. NaN values do
% not cause an error and are never replaced by zero.

maxNames = 4;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

variableNames = { ...
    'Al_cation_apfu', ...
    'Si_cation_apfu', ...
    'Ca_cation_apfu'};

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    variableValue = data_amphibole.(variableName);

    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = ...
            "Amphibole." + string(variableName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_amphibole)
% validateNonNegativeInputs
% Reject Inf and finite negative values in all amphibole cation variables
% read by this implementation. NaN and zero are allowed and retained.

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu'};

optionalVariables = { ...
    'Fe3_cation_apfu', ...
    'Mn_cation_apfu', ...
    'K_cation_apfu', ...
    'Cr_cation_apfu'};

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        error('Amphibole table must contain variable: %s', variableName);
    end
end

allVariables = [requiredVariables, optionalVariables];
invalidInputBuffer = strings(numel(allVariables), 1);
nInvalidInputs = 0;

for i = 1:numel(allVariables)
    variableName = allVariables{i};

    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        continue;
    end

    variableValue = data_amphibole.(variableName);

    if ~isscalar(variableValue)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end

    if any(isinf(variableValue(:)) | ...
            (isfinite(variableValue(:)) & variableValue(:) < 0))
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Amphibole." + string(variableName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['JohnsonRutherford1988: amphibole cation values must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_amphibole, T_degreeC)
% calcPressure
% Compute pressure for one amphibole row at one or more input temperatures.
% Temperature does not enter the published pressure equation. NaN values
% propagate through all calculations in which they are used.
%
% Inputs:
%   data_amphibole : 1-row Amphibole table
%   T_degreeC      : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

amp = prepareAmphiboleRow(data_amphibole);

% Preserve NaN explicitly in derived Al-site quantities.
AlT_scalar = amp.Al;

AlIV_scalar = 8 - amp.Si;
if isfinite(AlIV_scalar) && AlIV_scalar < 0
    AlIV_scalar = 0;
end

AlVI_scalar = AlT_scalar - AlIV_scalar;
if isfinite(AlVI_scalar) && AlVI_scalar < 0
    AlVI_scalar = 0;
end

% Final experimental regression of Johnson and Rutherford (1989), Eq. (3).
P_kbar_scalar = 4.23 .* AlT_scalar - 3.46;

% Numerical screening and range flags.
isIgneousChemicalScreen_scalar = ...
    isfinite(amp.Si) && isfinite(amp.Ca) && ...
    amp.Si < 7.5 && amp.Ca > 1.6;

isWithinCalibrationTRange = ...
    isfinite(T_degreeC) & ...
    T_degreeC >= 740 & T_degreeC <= 780;

isWithinCalibrationP_scalar = ...
    isfinite(P_kbar_scalar) && ...
    P_kbar_scalar >= 2 && P_kbar_scalar <= 8;

isWithinExperimentalAlTRange_scalar = ...
    isfinite(AlT_scalar) && ...
    AlT_scalar >= 1.27 && AlT_scalar <= 2.74;

% This is only a numerical diagnostic flag. Full application conditions must
% be assessed from petrography and phase-equilibrium relationships.
isApplicable_numeric = ...
    isWithinCalibrationTRange & ...
    repmat(isWithinCalibrationP_scalar, nT, 1) & ...
    repmat(isWithinExperimentalAlTRange_scalar, nT, 1) & ...
    repmat(isIgneousChemicalScreen_scalar, nT, 1);

% Expand scalar composition and pressure outputs to the temperature-vector
% length. All table variables are created at their final size.
row = table();

row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.Si_amp = repmat(amp.Si, nT, 1);
row.Ti_amp = repmat(amp.Ti, nT, 1);
row.AlT_amp = repmat(AlT_scalar, nT, 1);
row.Al_tot_amp = row.AlT_amp;
row.AlIV_amp = repmat(AlIV_scalar, nT, 1);
row.Al_IV_amp = row.AlIV_amp;
row.AlVI_amp = repmat(AlVI_scalar, nT, 1);
row.Al_VI_amp = row.AlVI_amp;

row.FeT_amp = repmat(amp.FeT, nT, 1);
row.Fe2_amp = repmat(amp.Fe2, nT, 1);
row.Fe3_amp = repmat(amp.Fe3, nT, 1);
row.Mg_amp = repmat(amp.Mg, nT, 1);
row.Ca_amp = repmat(amp.Ca, nT, 1);
row.Na_amp = repmat(amp.Na, nT, 1);
row.K_amp = repmat(amp.K, nT, 1);
row.Mn_amp = repmat(amp.Mn, nT, 1);
row.Cr_amp = repmat(amp.Cr, nT, 1);

row.P_kbar = repmat(P_kbar_scalar, nT, 1);
row.P_uncertainty_kbar = repmat(0.5, nT, 1);
row.P_uncertainty_at_8kbar_kbar = repmat(0.6, nT, 1);
row.P_regression_std_kbar = repmat(0.22, nT, 1);

row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinCalibrationP = ...
    repmat(isWithinCalibrationP_scalar, nT, 1);
row.isWithinExperimentalAlTRange = ...
    repmat(isWithinExperimentalAlTRange_scalar, nT, 1);
row.isIgneousChemicalScreen = ...
    repmat(isIgneousChemicalScreen_scalar, nT, 1);
row.isApplicable_numeric = isApplicable_numeric;

% Backward-compatible output name. This remains a numerical diagnostic and
% does not verify the full petrographic application conditions.
row.isApplicable = isApplicable_numeric;

end

function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Extract one-row amphibole cation data. Required variables must exist.
% Missing optional variables are represented as NaN, never as zero.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();

amp.Si = getVarOrError(data_amphibole, ...
    'Si_cation_apfu', 'Amphibole');
amp.Ti = getVarOrError(data_amphibole, ...
    'Ti_cation_apfu', 'Amphibole');
amp.Al = getVarOrError(data_amphibole, ...
    'Al_cation_apfu', 'Amphibole');
amp.FeT = getVarOrError(data_amphibole, ...
    'Fe_cation_apfu', 'Amphibole');
amp.Mg = getVarOrError(data_amphibole, ...
    'Mg_cation_apfu', 'Amphibole');
amp.Ca = getVarOrError(data_amphibole, ...
    'Ca_cation_apfu', 'Amphibole');
amp.Na = getVarOrError(data_amphibole, ...
    'Na_cation_apfu', 'Amphibole');

amp.Fe3 = getVarOrNaN(data_amphibole, 'Fe3_cation_apfu');
amp.Mn = getVarOrNaN(data_amphibole, 'Mn_cation_apfu');
amp.K = getVarOrNaN(data_amphibole, 'K_cation_apfu');
amp.Cr = getVarOrNaN(data_amphibole, 'Cr_cation_apfu');

% Fe2 remains NaN when total Fe or Fe3 is NaN. Missing Fe3 is not interpreted
% as zero.
if isfinite(amp.Fe3) && isfinite(amp.FeT) && amp.Fe3 > amp.FeT
    error('Amphibole contains Fe3_cation_apfu > Fe_cation_apfu.');
end

amp.Fe2 = amp.FeT - amp.Fe3;

% A negative derived Fe2 value is physically invalid. This can only occur
% when both total Fe and Fe3 are finite and inconsistent.
if isfinite(amp.Fe2) && amp.Fe2 < 0
    error('Amphibole contains a negative derived Fe2 value.');
end

end

function value = getVarOrError(tbl, variableName, mineralLabel)
% getVarOrError
% Retrieve one required scalar variable without altering NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);

if ~isscalar(value)
    error('Variable %s must be scalar in a 1-row table.', variableName);
end

if isinf(value) || (isfinite(value) && value < 0)
    error('Variable %s must be non-negative or NaN.', variableName);
end

end

function value = getVarOrNaN(tbl, variableName)
% getVarOrNaN
% Retrieve one optional scalar variable. Missing variables and explicit NaN
% values remain NaN and are never interpreted as zero.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);

    if ~isscalar(value)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end

    if isinf(value) || (isfinite(value) && value < 0)
        error('Variable %s must be non-negative or NaN.', variableName);
    end
else
    value = NaN;
end

end
