function results = Ai1994(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/Ai1994.m
% Tested with MATLAB R2024b
%
% Garnet-clinopyroxene Fe2+-Mg exchange thermometer
% Ai, Y. (1994)
% Contributions to Mineralogy and Petrology, 115, 467-473
% DOI: https://doi.org/10.1007/BF00320979
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis (selected by the user from tables) and calculates
% temperature using the Ai (1994) garnet-clinopyroxene Fe2+-Mg exchange
% thermometer.
%
% The function accepts pressure as either a scalar or a vector. Therefore,
% it can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. For each selected Grt-Cpx pair, the output table
% contains one row per pressure value.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Grt-Cpx pair, and appends results into
% a single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Ai (1994) assembled experimental garnet-clinopyroxene equilibrium data
% from ultramafic and mafic systems over the following overall range:
%
%   Temperature : 600-1500 degreeC
%   Pressure    : 10-60 kbar (1-6 GPa)
%   Composition : predominantly ultramafic systems, with additional mafic
%                 systems and experiments designed specifically to examine
%                 garnet-clinopyroxene equilibrium
%
% The dataset and its P-T range are described on p. 468. After analytical
% quality screening, 380 Grt-Cpx pairs were considered; 109 statistical
% outliers were subsequently excluded, and the final regression used 271
% pairs (pp. 468-469).
%
% Important application notes from Ai (1994):
%   1) Experimental data are much less abundant below 900 degreeC and above
%      1300 degreeC than in the main 1000-1200 degreeC interval (p. 468).
%   2) Low-temperature experiments may fail to attain equilibrium,
%      especially in anhydrous systems. At high temperature, rapid
%      nucleation may preserve metastable mineral compositions (p. 468).
%   3) The thermometer reproduced the retained experimental run
%      temperatures within approximately +/-100 degreeC (p. 469; Fig. 4 on
%      p. 472). This is calibration-data reproduction, not a universal
%      uncertainty for natural samples.
%   4) Unreliable estimates of Fe3+ in garnet and clinopyroxene and lack of
%      mutual Grt-Cpx equilibrium are identified as two major error sources
%      in natural rocks (pp. 469-470). Zoned minerals should therefore be
%      paired using texturally and chemically equivalent domains.
%   5) The thermometer was applied to lower-crustal amphibolites,
%      granulites and eclogites, and to upper-mantle eclogite/lherzolite
%      xenoliths and mineral inclusions in diamonds (pp. 469-472).
%   6) Ai (1994) did not specify strict numerical calibration limits for
%      garnet XCa or Mg number. The equation is empirical, so unusually
%      extreme compositions should be assessed independently.
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) input pressure is outside 10-60 kbar,
%   2) a finite calculated temperature is outside 600-1500 degreeC,
%   3) a finite calculated temperature is below 900 degreeC or above
%      1300 degreeC, where calibration data are relatively sparse,
%   4) a thermometer input contains NaN, or
%   5) a calculated temperature is NaN or Inf. In this case, the input
%      values and identified calculation-domain causes are also printed.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Cpx    : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must contain
% normalized cation data.
%
% Variables used directly by the thermometer:
%   Garnet table:
%     Fe_cation_apfu         % Fe2+ in garnet
%     Mg_cation_apfu
%     Ca_cation_apfu
%     Mn_cation_apfu         % optional; zero if the column is absent
%
%   Clinopyroxene table:
%     Fe_cation_apfu         % Fe2+ in clinopyroxene
%     Mg_cation_apfu
%
% Optional variables retained in the output when available:
%   Fe3_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%
% IMPORTANT Fe note:
% Ai (1994) defines KD, garnet XCa, and garnet Mg number using Fe2+, not
% total Fe. Therefore Fe3_cation_apfu is not added to Fe_cation_apfu in the
% calculation. If Fe_cation_apfu contains total Fe rather than Fe2+, Fe2+
% must be estimated before this function is used.
%
% Negative finite values in variables used by the thermometer are not
% permitted. Zero values are retained; if they make the equation
% mathematically undefined, the affected result is returned as NaN and a
% non-stopping warning is printed. NaN values are retained as missing
% values, propagated through the calculation, and reported by non-stopping
% warnings. An absent optional Mn column is treated as zero, whereas an
% explicitly stored NaN value remains NaN.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Ai (1994) used the Fe2+-Mg distribution between coexisting garnet and
% clinopyroxene:
%
%   KD = (Fe2+/Mg)_Grt / (Fe2+/Mg)_Cpx
%
% Garnet composition terms:
%
%   XCa_Grt  = Ca / (Fe2+ + Mg + Mn + Ca)
%   MgNo_Grt = 100 * Mg / (Fe2+ + Mg)
%
% Temperature equation (Eq. 3 on p. 469):
%
%   T(K) = [-1629*(XCa_Grt)^2 + 3648.55*XCa_Grt
%           - 6.59*MgNo_Grt + 1987.98 + 17.66*P_kbar]
%          / [ln(KD) + 1.076]
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Ai1994(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Cpx tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Grt-Cpx pair
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('Ai1994 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The first column of each
% table is used only as the data-code identifier in the selection dialog.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

dataset_grt = rawdata_struct.Garnet;
dataset_cpx = rawdata_struct.Cpx;

validateRequiredVariables(dataset_grt, dataset_cpx);

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated concatenation of the complete results table inside the loop is
% avoided because it repeatedly reallocates and copies the table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Overall experimental calibration limits reported by Ai (1994, p. 468).
calibrationT_min_degC = 600;
calibrationT_max_degC = 1500;
calibrationP_min_kbar = 10;
calibrationP_max_kbar = 60;

% Temperature intervals with relatively sparse experimental data
% (Ai, 1994, p. 468).
sparseT_lower_degC = 900;
sparseT_upper_degC = 1300;

% Pressure is common to all selected mineral pairs in this function call.
% The warning is therefore printed only once, after the first calculation.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
    dataCodes_grt = dataset_grt{:, 1};

    [selectedIdx_grt, ok] = listdlg( ...
        'PromptString', 'Please select the Garnet data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_grt)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_grt)
        disp('Selection canceled');
        break;
    end

    selectedCode_grt = dataCodes_grt(selectedIdx_grt);
    disp(['Garnet selected: ' char(string(selectedCode_grt))]);

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
    % Garnet and clinopyroxene are selected independently; their row indices
    % are not assumed to correspond.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    % NaN is reported but deliberately allowed to propagate.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_cpx);

    % Only negative finite thermometer inputs stop the calculation. Zero and
    % NaN are retained and handled by the result-domain checks and warnings.
    validateNonNegativeInputs(selectedData_grt, selectedData_cpx);

    row = calcTemp(selectedData_grt, selectedData_cpx, P_kbar);

    % Store identifiers as one value per pressure point.
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_cpx'}, 'Before', 1);

    % Store this result as one table block. The capacity is doubled only when
    % the preallocated buffer is exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_C(1)) ' to ' num2str(row.T_C(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside 10-60 kbar.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental calibration range ' ...
             'of Ai (1994): 10-60 kbar (1-6 GPa). ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when a finite calculated temperature is outside 600-1500 degreeC.
    finiteTemperature = isfinite(row.T_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_C < calibrationT_min_degC | ...
         row.T_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the overall experimental ' ...
             'calibration range of Ai (1994): 600-1500 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g-%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)));
    end

    % Within the overall calibration range, identify the low- and
    % high-temperature intervals for which experimental data are sparse.
    temperatureInSparseRegion = finiteTemperature & ...
        ~temperatureOutsideCalibration & ...
        (row.T_C < sparseT_lower_degC | row.T_C > sparseT_upper_degC);

    if any(temperatureInSparseRegion)
        sparseValues = row.T_C(temperatureInSparseRegion);
        fprintf(2, ...
            ['CAUTION: Ai (1994) reports relatively sparse experimental data below ' ...
             '900 degreeC and above 1300 degreeC (p. 468). ' ...
             '%d calculated point(s) fall in these intervals; ' ...
             'range = %.4g-%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureInSparseRegion), ...
            min(sparseValues), ...
            max(sparseValues), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)));
    end

    % Report explicitly stored NaN values without interrupting calculation.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, and NaN was not replaced by zero.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Preserve and report NaN/Inf output values caused by missing inputs or
    % by a mathematically undefined exchange coefficient or denominator.
    % The original implementation reported only the number of non-finite
    % values. Here, the actual thermometer inputs and the identified causes
    % are also printed so that a zero-derived NaN is not mistaken for an
    % explicitly stored NaN input.
    invalidTemperature = ~isfinite(row.T_C);
    if any(invalidTemperature)
        nonFiniteCauses = findNonFiniteCauses(row);

        fprintf(2, ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
                    '(%d of %d points; NaN: %d, Inf: %d).\n' ...
                    '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), ...
            numel(row.T_C), ...
            sum(isnan(row.T_C)), ...
            sum(isinf(row.T_C)));

        fprintf(2, ['         Thermometer inputs used: ' ...
                    'Garnet.Fe_cation_apfu=%s, Garnet.Mg_cation_apfu=%s, ' ...
                    'Garnet.Ca_cation_apfu=%s, Garnet.Mn_cation_apfu=%s, ' ...
                    'Cpx.Fe_cation_apfu=%s, Cpx.Mg_cation_apfu=%s.\n'], ...
            formatNumericValue(row.Fe2_grt(1)), ...
            formatNumericValue(row.Mg_grt(1)), ...
            formatNumericValue(row.Ca_grt(1)), ...
            formatNumericValue(row.Mn_grt(1)), ...
            formatNumericValue(row.Fe2_cpx(1)), ...
            formatNumericValue(row.Mg_cpx(1)));

        if isempty(nonFiniteCauses)
            fprintf(2, ['         No explicit NaN, zero, Inf, or invalid intermediate value ' ...
                        'was identified; inspect the stored intermediate variables.\n']);
        else
            fprintf(2, '         Identified cause(s): %s.\n', ...
                char(strjoin(nonFiniteCauses, '; ')));
        end
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Ai1994', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered result blocks once after all selections are complete.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(dataset_grt, dataset_cpx)
% validateRequiredVariables
% Verify the columns that must exist before opening the selection dialogs.
% Mn is optional and is treated as zero only when its column is absent.

requiredGrt = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu'};
requiredCpx = {'Fe_cation_apfu', 'Mg_cation_apfu'};

missingNames = strings(numel(requiredGrt) + numel(requiredCpx), 1);
nMissing = 0;

for i = 1:numel(requiredGrt)
    if ~ismember(requiredGrt{i}, dataset_grt.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Garnet." + string(requiredGrt{i});
    end
end

for i = 1:numel(requiredCpx)
    if ~ismember(requiredCpx{i}, dataset_cpx.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Cpx." + string(requiredCpx{i});
    end
end

if nMissing > 0
    missingNames = missingNames(1:nMissing);
    error(['Ai1994: required variable(s) are missing: ' ...
        char(strjoin(missingNames, ', ')) '.']);
end

end

function nanInputNames = findNaNInputs(data_grt, data_cpx)
% findNaNInputs
% Return names of thermometer inputs containing NaN. Missing optional Mn is
% not reported because its documented default is zero.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

nanInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nNaN = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    if ismember(variableName, data_grt.Properties.VariableNames)
        variableValue = data_grt.(variableName);
        if any(isnan(variableValue(:)))
            nNaN = nNaN + 1;
            nanInputNames(nNaN) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(cpxVariables)
    variableName = cpxVariables{i};
    variableValue = data_cpx.(variableName);
    if any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanInputNames(nNaN) = "Cpx." + string(variableName);
    end
end

nanInputNames = nanInputNames(1:nNaN);

end

function validateNonNegativeInputs(data_grt, data_cpx)
% validateNonNegativeInputs
% Stop when a finite thermometer input is negative. Zero and NaN are
% deliberately allowed so that undefined results can be retained as NaN and
% reported by non-stopping warnings.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

invalidInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nInvalid = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    if ismember(variableName, data_grt.Properties.VariableNames)
        variableValue = data_grt.(variableName);
        if any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nInvalid = nInvalid + 1;
            invalidInputNames(nInvalid) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(cpxVariables)
    variableName = cpxVariables{i};
    variableValue = data_cpx.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Cpx." + string(variableName);
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNames(1:nInvalid);
    error(['Ai1994: thermometer inputs must be >= 0. ' ...
        'Negative finite value(s) were found in: ' ...
        char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function nonFiniteCauses = findNonFiniteCauses(row)
% findNonFiniteCauses
% Identify raw-input and derived-variable conditions that can produce a NaN
% or Inf temperature. The returned strings are used only for non-stopping
% diagnostics and do not alter any stored value.

maximumCauses = 20;
nonFiniteCauses = strings(maximumCauses, 1);
nCauses = 0;

inputValues = {row.Fe2_grt, row.Mg_grt, row.Ca_grt, row.Mn_grt, ...
    row.Fe2_cpx, row.Mg_cpx};
inputNames = {"Garnet.Fe_cation_apfu", "Garnet.Mg_cation_apfu", ...
    "Garnet.Ca_cation_apfu", "Garnet.Mn_cation_apfu", ...
    "Cpx.Fe_cation_apfu", "Cpx.Mg_cation_apfu"};

for i = 1:numel(inputValues)
    value = inputValues{i};
    if any(isnan(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is NaN";
    elseif any(isinf(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is Inf";
    elseif any(value(:) == 0)
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is zero";
    end
end

if any(~isfinite(row.KD))
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "KD is non-finite";
elseif any(row.KD <= 0)
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "KD is zero or negative";
end

if any(~isfinite(row.XCa_grt))
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "XCa_grt is non-finite";
elseif any(row.XCa_grt < 0)
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "XCa_grt is negative";
end

if any(~isfinite(row.MgNo_grt))
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "MgNo_grt is non-finite";
elseif any(row.MgNo_grt < 0)
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "MgNo_grt is negative";
end

if any(~isfinite(row.numerator))
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "temperature numerator is non-finite";
end

if any(~isfinite(row.denominator))
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "temperature denominator is non-finite";
elseif any(row.denominator == 0)
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "temperature denominator is zero";
end

nonFiniteCauses = nonFiniteCauses(1:nCauses);

end

function textValue = formatNumericValue(value)
% formatNumericValue
% Format a scalar numeric input for a compact diagnostic message while
% preserving the explicit labels NaN, Inf, and -Inf.

if isnan(value)
    textValue = 'NaN';
elseif isinf(value) && value > 0
    textValue = 'Inf';
elseif isinf(value) && value < 0
    textValue = '-Inf';
else
    textValue = sprintf('%.8g', value);
end

end

function row = calcTemp(data_grt, data_cpx, P_kbar)
% calcTemp
% Compute Ai (1994) temperatures for one garnet row, one clinopyroxene row,
% and a scalar or vector of pressures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);
P_bar = P_kbar .* 1000;

row = table();

% --- Store pressure in both convenient units ---
row.P_kbar = P_kbar;
row.P_bar = P_bar;

% --- Extract and expand garnet cations ---
% Fe_cation_apfu is Fe2+; Fe3+ is retained separately and is not used in the
% Ai (1994) equation.
Fe2_grt = repmat(getRequiredValue(data_grt, 'Fe_cation_apfu', 'Garnet'), nP, 1);
Fe3_grt = repmat(getOptionalValue(data_grt, 'Fe3_cation_apfu', 0, 'Garnet'), nP, 1);
Mg_grt  = repmat(getRequiredValue(data_grt, 'Mg_cation_apfu', 'Garnet'), nP, 1);
Mn_grt  = repmat(getOptionalValue(data_grt, 'Mn_cation_apfu', 0, 'Garnet'), nP, 1);
Ca_grt  = repmat(getRequiredValue(data_grt, 'Ca_cation_apfu', 'Garnet'), nP, 1);
Al_grt  = repmat(getOptionalValue(data_grt, 'Al_cation_apfu', 0, 'Garnet'), nP, 1);
Si_grt  = repmat(getOptionalValue(data_grt, 'Si_cation_apfu', 0, 'Garnet'), nP, 1);

% --- Extract and expand clinopyroxene cations ---
Fe2_cpx = repmat(getRequiredValue(data_cpx, 'Fe_cation_apfu', 'Cpx'), nP, 1);
Fe3_cpx = repmat(getOptionalValue(data_cpx, 'Fe3_cation_apfu', 0, 'Cpx'), nP, 1);
Mg_cpx  = repmat(getRequiredValue(data_cpx, 'Mg_cation_apfu', 'Cpx'), nP, 1);
Mn_cpx  = repmat(getOptionalValue(data_cpx, 'Mn_cation_apfu', 0, 'Cpx'), nP, 1);
Ca_cpx  = repmat(getOptionalValue(data_cpx, 'Ca_cation_apfu', 0, 'Cpx'), nP, 1);
Ti_cpx  = repmat(getOptionalValue(data_cpx, 'Ti_cation_apfu', 0, 'Cpx'), nP, 1);
Al_cpx  = repmat(getOptionalValue(data_cpx, 'Al_cation_apfu', 0, 'Cpx'), nP, 1);
Si_cpx  = repmat(getOptionalValue(data_cpx, 'Si_cation_apfu', 0, 'Cpx'), nP, 1);
Na_cpx  = repmat(getOptionalValue(data_cpx, 'Na_cation_apfu', 0, 'Cpx'), nP, 1);
K_cpx   = repmat(getOptionalValue(data_cpx, 'K_cation_apfu', 0, 'Cpx'), nP, 1);

% --- Ai (1994) composition terms ---
FeMg_grt = Fe2_grt ./ Mg_grt;
FeMg_cpx = Fe2_cpx ./ Mg_cpx;
KD = FeMg_grt ./ FeMg_cpx;
lnKD = log(KD);

xSiteSum_grt = Fe2_grt + Mg_grt + Mn_grt + Ca_grt;
XCa_grt = Ca_grt ./ xSiteSum_grt;

mgFeSum_grt = Fe2_grt + Mg_grt;
MgNo_grt = 100 .* Mg_grt ./ mgFeSum_grt;

% --- Temperature equation (Ai, 1994, Eq. 3) ---
numerator = -1629 .* (XCa_grt .^ 2) ...
    + 3648.55 .* XCa_grt ...
    - 6.59 .* MgNo_grt ...
    + 1987.98 ...
    + 17.66 .* P_kbar;
denominator = lnKD + 1.076;

T_K = numerator ./ denominator;

% Explicitly return NaN when a zero input or another domain problem makes
% the exchange expression undefined. Existing input NaNs remain NaN.
validDomain = isfinite(KD) & KD > 0 ...
    & isfinite(XCa_grt) & XCa_grt >= 0 ...
    & isfinite(MgNo_grt) & MgNo_grt >= 0 ...
    & isfinite(numerator) & isfinite(denominator) ...
    & denominator ~= 0;
T_K(~validDomain) = NaN;
T_C = T_K - 273.15;

% --- Pack outputs ---
% FeUsed is retained for compatibility with the original Ai1994 output, but
% it now correctly equals Fe2+ rather than Fe2+ + Fe3+.
row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.FeUsed_grt = Fe2_grt;
row.Mg_grt = Mg_grt;
row.Mn_grt = Mn_grt;
row.Ca_grt = Ca_grt;
row.Al_grt = Al_grt;
row.Si_grt = Si_grt;

row.Fe2_cpx = Fe2_cpx;
row.Fe3_cpx = Fe3_cpx;
row.FeUsed_cpx = Fe2_cpx;
row.Mg_cpx = Mg_cpx;
row.Mn_cpx = Mn_cpx;
row.Ca_cpx = Ca_cpx;
row.Ti_cpx = Ti_cpx;
row.Al_cpx = Al_cpx;
row.Si_cpx = Si_cpx;
row.Na_cpx = Na_cpx;
row.K_cpx = K_cpx;

row.FeMg_grt = FeMg_grt;
row.FeMg_cpx = FeMg_cpx;
row.XCa_grt = XCa_grt;
row.MgNo_grt = MgNo_grt;
row.KD = KD;
row.lnKD = lnKD;
row.numerator = numerator;
row.denominator = denominator;
row.T_K = T_K;
row.T_C = T_C;

end

function value = getRequiredValue(data_tbl, variableName, mineralLabel)
% getRequiredValue
% Extract a required scalar numeric value. NaN is returned unchanged.

if ~ismember(variableName, data_tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = data_tbl.(variableName);
if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end

function value = getOptionalValue(data_tbl, variableName, defaultValue, mineralLabel)
% getOptionalValue
% Extract an optional scalar numeric value. Use defaultValue only when the
% column is absent; an explicitly stored NaN is returned unchanged.

if ismember(variableName, data_tbl.Properties.VariableNames)
    value = data_tbl.(variableName);
    if ~isnumeric(value) || ~isscalar(value)
        error('%s.%s must be a numeric scalar in the selected row.', ...
            mineralLabel, variableName);
    end
else
    value = defaultValue;
end

end
