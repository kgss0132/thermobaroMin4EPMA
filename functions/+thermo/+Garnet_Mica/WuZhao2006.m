function results = WuZhao2006(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/WuZhao2006.m
% Tested with MATLAB R2024b
%
% Empirical Fe-Mg exchange thermometer between Garnet and Muscovite
% Wu, C.-M. and Zhao, G. (2006)
% Journal of Petrology, 47, 2357-2368
% DOI: https://doi.org/10.1093/petrology/egl047
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% analysis (selected independently by the user from tables) and calculates
% temperature using the Wu and Zhao (2006) garnet-muscovite thermometer.
%
% Two temperatures are returned:
%   Model A : assuming no Fe3+ in muscovite
%   Model B : assuming 50% of total Fe in muscovite is Fe3+
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Garnet-Mica pair, the output contains one row per pressure
% value. This allows the same function to be called by both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Garnet-Mica pair and appends the new
% table block to the final output.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Wu and Zhao (2006) calibrated the GM thermometer together with the GMPQ
% barometer using 103 natural metapelitic samples over the following ranges:
%
%   Temperature               : 450-760 degreeC
%                               (actual calibrant range: 452-758 degreeC)
%   Pressure                  : 0.8-11.1 kbar
%   Garnet composition        : XFe = 0.53-0.81
%                               XMg = 0.05-0.24
%                               XCa = 0.03-0.23
%   Muscovite composition     : Fe   = 0.04-0.16 apfu
%                               Mg   = 0.04-0.13 apfu
%                               AlVI = 1.74-1.96 apfu
%                               (11-oxygen basis)
%   Intended application      : equilibrated metapelitic assemblages
%
% The rounded P-T and mineral-composition ranges are summarized in the
% abstract on p. 2357. The 103 calibration samples, actual P-T range, mineral
% ranges, and Fe3+ assumptions are described on p. 2358. The thermometer
% equations (11a and 11b) are given on p. 2359.
%
% Important application cautions stated or demonstrated in the article:
% - This is an empirical calibration based on natural metapelites, not a
%   direct experimental calibration. Reference P-T values were obtained from
%   the Holdaway (2000) Garnet-Biotite thermometer and Holdaway (2001) GASP
%   barometer (pp. 2357-2358).
% - Fe3+ in muscovite was not measured directly. Model A assumes no ferric
%   iron, whereas Model B assumes that 50% of total muscovite Fe is ferric
%   (pp. 2358-2359). The application tests favored Model A; Model B produced
%   implausible results for some metamorphic sequences (pp. 2361-2363).
% - The calibration calculations assumed garnet Fe3+ to be 3% of total
%   garnet Fe (p. 2358). This implementation follows the existing input
%   convention in which Fe_cation_apfu supplies the Fe2+ term; an optional
%   Fe3_cation_apfu value is retained in the output but is not added to XFe.
% - Test samples with mineral compositions outside the calibration range
%   were excluded. The authors explicitly do not advocate applying the GM
%   thermometer to metapelites with muscovite Mg > 0.13 apfu or
%   Fe < 0.04 apfu on an 11-oxygen basis (pp. 2361-2363 and p. 2366).
% - The estimated total random error is no more than approximately
%   +/-16 degreeC, but this is not an absolute-accuracy statement. No direct
%   experiments calibrated the GM thermometer; by analogy with the
%   Garnet-Biotite thermometer, the absolute error may be approximately
%   70 degreeC (p. 2366).
% - Pressure sensitivity is small. For the calibration samples, a pressure
%   error of +/-2 kbar produces less than approximately +/-1.1 degreeC error
%   in the calculated GM temperature (pp. 2361 and 2366).
% - When biotite and aluminosilicate are both available, the authors prefer
%   the experimentally calibrated Holdaway (2000) Garnet-Biotite thermometer
%   and Holdaway (2001) GASP barometer (p. 2366).
%
% This implementation uses 450-760 degreeC and 0.8-11.1 kbar as the
% non-stopping warning limits. It also reports finite mineral compositions
% outside the ranges represented by the calibrants. All warning messages are
% printed with fprintf and do not stop the calculation.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Mica   : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following variable names (normalized cations):
%
%   Garnet table variables used in the calculation:
%     Fe_cation_apfu       % Fe2+ in garnet (assumed)
%     Mg_cation_apfu
%
%   Mica table variables used in the calculation:
%     Fe_cation_apfu       % Fe2+ input before Model A/B reassignment
%     Mg_cation_apfu
%     Al_cation_apfu
%     Si_cation_apfu
%
% Optional variables used in the calculation when present:
%   Garnet : Mn_cation_apfu, Ca_cation_apfu
%   Mica   : Fe3_cation_apfu
%
% Optional variables retained in the output table when present:
%   Garnet : Fe3_cation_apfu, Al_cation_apfu, Si_cation_apfu
%   Mica   : Mn_cation_apfu, Ca_cation_apfu, Ti_cation_apfu,
%            K_cation_apfu, Na_cation_apfu
%
% Missing optional variables are represented by zero. An optional variable
% that exists and contains NaN remains NaN; it is never converted to zero.
%
% All finite mineral-composition values used by the thermometer must be
% non-negative. Negative and infinite values stop the calculation. Zero is
% allowed, although it may produce a non-finite exchange coefficient or
% temperature. NaN values are retained, propagated through the calculation,
% and reported by non-stopping fprintf messages.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Wu and Zhao, 2006, Eqs. 11a-11b, p. 2359)
%
% Wu & Zhao (2006) recalibrated the garnet-muscovite (GM) thermometer:
%
%   Model A (no Fe3+ in muscovite):
%
%   T_A(K) =
%   [2325.8 + P(kbar){-0.1 - 13.5(Fe_b - Mg_b)} - 0.0135(Fe_c - Mg_c)
%    - 6541.2(X_Fe_mus - X_Mg_mus) - 1127.7 X_Al_mus]
%   / [1 + 0.0135{R ln(Kideal) + (Fe_a - Mg_a)}]
%
%   Model B (50% of total Fe in muscovite as Fe3+):
%
%   T_B(K) =
%   [2064.7 + P(kbar){-0.7 - 9.8(Fe_b - Mg_b)} - 0.0098(Fe_c - Mg_c)
%    - 7077.9(X_Fe_mus - X_Mg_mus) - 941.7 X_Al_mus]
%   / [1 + 0.0098{R ln(Kideal) + (Fe_a - Mg_a)}]
%
% where:
%
%   Kideal = (X_Fe_grt^3 * X_Mg_mus^3) / (X_Mg_grt^3 * X_Fe_mus^3)
%
%   X_Fe_grt = Fe2+ / (Fe2+ + Mg + Ca + Mn)
%   X_Mg_grt = Mg   / (Fe2+ + Mg + Ca + Mn)
%   X_Ca_grt = Ca   / (Fe2+ + Mg + Ca + Mn)
%   X_Mn_grt = Mn   / (Fe2+ + Mg + Ca + Mn)
%
%   X_Fe_mus = Fe2+ / (Fe2+ + Mg + AlVI)
%   X_Mg_mus = Mg   / (Fe2+ + Mg + AlVI)
%   X_Al_mus = AlVI / (Fe2+ + Mg + AlVI)
%
% Garnet activity terms Fe_a, Fe_b, Fe_c, Mg_a, Mg_b, Mg_c are computed
% from the Holdaway (2000, 2001) average garnet activity model as written
% in the Appendix of Wu et al. (2004b), which is explicitly referenced by
% Wu & Zhao (2006).
%
% NOTES
% - For garnet, Fe_cation_apfu is treated as Fe2+. If Fe3_cation_apfu
%   exists, it is stored in the output but is not used in the Fe2+-based
%   X_Fe_grt required by the thermometer.
% - For muscovite:
%     Model A uses total Fe as Fe2+  (Fe2_A = Fe2 + Fe3_input)
%     Model B uses 50% of total Fe as Fe2+ (Fe2_B = 0.5 * total Fe)
% - AlVI in mica is calculated as:
%     AlIV = max(4 - Si, 0)
%     AlVI = Al_total - AlIV
%   Existing NaN values remain NaN; a finite negative AlVI is rejected.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = WuZhao2006(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Mica tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Mica pair. The output variable set is
%             intended to remain stable for downstream processing.

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('WuZhao2006 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        ~isreal(P_kbar) || any(~isfinite(P_kbar(:))) || ...
        any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative real numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The tables are not
% modified here; only the relevant columns are read during calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_grt = rawdata_struct.Garnet;
dataset_mica = rawdata_struct.Mica;

requiredVariables_grt = {'Fe_cation_apfu', 'Mg_cation_apfu'};
requiredVariables_mica = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Al_cation_apfu', 'Si_cation_apfu'};

validateRequiredVariables(dataset_grt, requiredVariables_grt, 'Garnet');
validateRequiredVariables(dataset_mica, requiredVariables_mica, 'Mica');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated table concatenation inside the loop is avoided because it forces
% MATLAB to repeatedly reallocate and copy the entire results table.
%
% The cell buffer is allocated in blocks and expanded only when necessary.
% After the interactive loop finishes, all blocks are concatenated once.
disp('=== Step 2: Preparing output container ===');

bufferBlockSize = 16;
resultBlocks = cell(bufferBlockSize, 1);
nResultBlocks = 0;

% Empirical calibration limits reported by Wu and Zhao (2006).
calibrationT_min_degC = 450;
calibrationT_max_degC = 760;
calibrationP_min_kbar = 0.8;
calibrationP_max_kbar = 11.1;

% Pressure is common to all selected mineral pairs in this function call.
% The pressure warning is therefore printed only once, after the first
% completed calculation.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
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

    % ----- Mica selection -----
    disp('=== Step 4: Selecting a data code from the list (Mica) ===');

    dataCodes_mica = dataset_mica{:, 1};

    [selectedIdx_mica, ok] = listdlg( ...
        'PromptString', 'Please select the Mica data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_mica)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_mica)
        disp('Selection canceled');
        break;
    end

    selectedCode_mica = dataCodes_mica(selectedIdx_mica);
    disp(['Mica selected: ' char(string(selectedCode_mica))]);

    % ----- Calculation -----
    % Garnet and mica are selected independently; do not assume row indices
    % correspond between the two datasets.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_mica = dataset_mica(selectedIdx_mica, :);

    % Check every raw variable actually used by the thermometer. The
    % calculation is intentionally allowed to continue when NaN is present.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_mica);

    % Negative, infinite, nonnumeric, nonreal, and nonscalar calculation
    % inputs stop the calculation. NaN and zero remain allowed.
    validateNonNegativeInputs(selectedData_grt, selectedData_mica);

    row = calcTemp(selectedData_grt, selectedData_mica, P_kbar);

    % Store the selected identifiers for every pressure point.
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_mica = repmat(string(selectedCode_mica), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_mica'}, 'Before', 1);

    % Store this result as one table block. Expand the buffer only in fixed
    % blocks, rather than growing the results table on every loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(bufferBlockSize, 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperatures for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': Model A = ' num2str(row.TA_C) ' degreeC, Model B = ' ...
            num2str(row.TB_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': Model A = ' num2str(row.TA_C(1)) ' to ' ...
            num2str(row.TA_C(end)) ' degreeC, Model B = ' ...
            num2str(row.TB_C(1)) ' to ' num2str(row.TB_C(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside the empirical
    % calibration range of 0.8-11.1 kbar. The calculation is not stopped.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the empirical calibration ' ...
             'range of Wu and Zhao (2006): 0.8-11.1 kbar. %d of %d ' ...
             'pressure point(s) are outside the range; input range = ' ...
             '%.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when finite Model A or Model B temperatures lie outside the
    % empirical calibration range. NaN and Inf are handled separately below.
    printTemperatureRangeMessage( ...
        row.TA_C, 'Model A', calibrationT_min_degC, ...
        calibrationT_max_degC, selectedCode_grt, selectedCode_mica);
    printTemperatureRangeMessage( ...
        row.TB_C, 'Model B', calibrationT_min_degC, ...
        calibrationT_max_degC, selectedCode_grt, selectedCode_mica);

    % Report mineral compositions outside the finite ranges represented by
    % the empirical calibration dataset. NaN inputs are handled separately.
    printCompositionRangeMessages(row, selectedCode_grt, selectedCode_mica);

    % Print a non-stopping warning when an input used by the thermometer
    % contains NaN. fprintf is used so the message remains visible even when
    % MATLAB warnings have been disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, and NaN values were not converted to zero.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report NaN/Inf results caused by missing data, zero
    % denominators, logarithms at zero, or other non-finite operations.
    invalidTA = ~isfinite(row.TA_C);
    invalidTB = ~isfinite(row.TB_C);

    if any(invalidTA) || any(invalidTB)
        fprintf(2, ...
            ['WARNING: Non-finite WuZhao2006 temperature values were calculated ' ...
             'for %s & %s. Model A: %d of %d points (NaN: %d, Inf: %d); ' ...
             'Model B: %d of %d points (NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            sum(invalidTA), numel(row.TA_C), ...
            sum(isnan(row.TA_C)), sum(isinf(row.TA_C)), ...
            sum(invalidTB), numel(row.TB_C), ...
            sum(isnan(row.TB_C)), sum(isinf(row.TB_C)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'WuZhao2006', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate the buffered table blocks only once after all selections have
% been completed. Return an empty table when no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(dataTable, requiredVariables, mineralLabel)
% validateRequiredVariables
% Confirm that every variable required by the thermometer exists before the
% user enters the interactive selection loop.

missingVariables = requiredVariables(~ismember( ...
    requiredVariables, dataTable.Properties.VariableNames));

if ~isempty(missingVariables)
    error('%s table must contain variable(s): %s', ...
        mineralLabel, char(strjoin(string(missingVariables), ', ')));
end

end

function nanInputNames = findNaNInputs(data_garnet, data_mica)
% findNaNInputs
% Return the names of raw thermometer inputs that contain NaN. The output
% buffer is preallocated so its size does not grow inside the checking loops.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Fe3_cation_apfu', ...
    'Mg_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu'};

nanInputNames = strings(numel(garnetVariables) + numel(micaVariables), 1);
nNaNInputs = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if isnumeric(variableValue) && any(isnan(variableValue(:)))
            nNaNInputs = nNaNInputs + 1;
            nanInputNames(nNaNInputs) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
    if ismember(variableName, data_mica.Properties.VariableNames)
        variableValue = data_mica.(variableName);
        if isnumeric(variableValue) && any(isnan(variableValue(:)))
            nNaNInputs = nNaNInputs + 1;
            nanInputNames(nNaNInputs) = "Mica." + string(variableName);
        end
    end
end

nanInputNames = nanInputNames(1:nNaNInputs);

end

function validateNonNegativeInputs(data_garnet, data_mica)
% validateNonNegativeInputs
% Stop when a raw input used by the thermometer is negative, infinite,
% nonnumeric, nonreal, or nonscalar. NaN and zero are intentionally allowed.
% The error-name buffer is preallocated and does not grow inside the loops.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Fe3_cation_apfu', ...
    'Mg_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu'};

invalidInputNames = strings(numel(garnetVariables) + numel(micaVariables), 1);
nInvalidInputs = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if isInvalidCationValue(variableValue)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputNames(nInvalidInputs) = ...
                "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
    if ismember(variableName, data_mica.Properties.VariableNames)
        variableValue = data_mica.(variableName);
        if isInvalidCationValue(variableValue)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputNames(nInvalidInputs) = ...
                "Mica." + string(variableName);
        end
    end
end

invalidInputNames = invalidInputNames(1:nInvalidInputs);

if ~isempty(invalidInputNames)
    error(['WuZhao2006: thermometer inputs must be real, scalar, ' ...
           'non-negative, and not Inf. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function invalid = isInvalidCationValue(value)
% isInvalidCationValue
% NaN and zero are valid at this stage. Negative, infinite, nonnumeric,
% nonreal, and nonscalar values are invalid.

invalid = ~isnumeric(value) || ~isscalar(value) || ~isreal(value);

if ~invalid
    invalid = isinf(value) || (isfinite(value) && value < 0);
end

end

function row = calcTemp(data_grt, data_mica, P_kbar)
% calcTemp
% Compute Wu and Zhao (2006) Eqs. (11a) and (11b) for one garnet row,
% one mica row, and every pressure supplied in P_kbar.
%
% Inputs:
%   data_grt  : 1-row table containing garnet apfu cations
%   data_mica : 1-row table containing mica apfu cations
%   P_kbar    : pressure in kbar (column vector after input validation)
%
% Output:
%   row : table with one row per pressure value, including mineral inputs,
%         compositional terms, activity terms, and both temperatures.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

row = table();

% --- Physical constant ---
R_J = 8.31446261815324;

% --- Prepare one-row mineral data ---
grt = prepareGarnetRow(data_grt);
mica = prepareMicaRow(data_mica);

% --- Replicate scalar mineral inputs for every pressure point ---
Fe2_grt = repmat(grt.Fe2, nP, 1);
Fe3_grt = repmat(grt.Fe3, nP, 1);
Mg_grt = repmat(grt.Mg, nP, 1);
Mn_grt = repmat(grt.Mn, nP, 1);
Ca_grt = repmat(grt.Ca, nP, 1);
Al_grt = repmat(grt.Al, nP, 1);
Si_grt = repmat(grt.Si, nP, 1);

Fe2_mica_input = repmat(mica.Fe2Input, nP, 1);
Fe3_mica_input = repmat(mica.Fe3Input, nP, 1);
FeTotal_mica = repmat(mica.FeTotal, nP, 1);
Mg_mica = repmat(mica.Mg, nP, 1);
Mn_mica = repmat(mica.Mn, nP, 1);
Ca_mica = repmat(mica.Ca, nP, 1);
Al_mica = repmat(mica.AlTotal, nP, 1);
Si_mica = repmat(mica.Si, nP, 1);
AlIV_mica = repmat(mica.AlIV, nP, 1);
AlVI_mica = repmat(mica.AlVI, nP, 1);
K_mica = repmat(mica.K, nP, 1);
Na_mica = repmat(mica.Na, nP, 1);
Ti_mica = repmat(mica.Ti, nP, 1);

% --- Garnet site fractions (Fe2+-Mg-Ca-Mn basis) ---
sumDiv_grt = Fe2_grt + Mg_grt + Ca_grt + Mn_grt;

XFe_grt = Fe2_grt ./ sumDiv_grt;
XMg_grt = Mg_grt ./ sumDiv_grt;
XCa_grt = Ca_grt ./ sumDiv_grt;
XMn_grt = Mn_grt ./ sumDiv_grt;

% --- Garnet activity polynomial terms ---
[Fea, Feb, Fec, Mga, Mgb, Mgc] = calcGarnetActivityTerms( ...
    XFe_grt, XMg_grt, XCa_grt, XMn_grt);

% --- Mica Model A: no Fe3+ in muscovite ---
Fe2_mica_A = FeTotal_mica;
[XFe_mica_A, XMg_mica_A, XAl_mica_A, sumDiv_mica_A] = ...
    calcMicaSiteFractions(Fe2_mica_A, Mg_mica, AlVI_mica);

Kideal_A = (XFe_grt.^3 .* XMg_mica_A.^3) ./ ...
    (XMg_grt.^3 .* XFe_mica_A.^3);
lnKideal_A = log(Kideal_A);

denA = 1 + 0.0135 .* (R_J .* lnKideal_A + (Fea - Mga));
numA = 2325.8 ...
    + P_kbar .* (-0.1 - 13.5 .* (Feb - Mgb)) ...
    - 0.0135 .* (Fec - Mgc) ...
    - 6541.2 .* (XFe_mica_A - XMg_mica_A) ...
    - 1127.7 .* XAl_mica_A;

% Direct calculation intentionally retains NaN and Inf results.
TA_K = numA ./ denA;
TA_C = TA_K - 273.15;

% --- Mica Model B: 50% of total Fe in muscovite as Fe3+ ---
Fe2_mica_B = 0.5 .* FeTotal_mica;
[XFe_mica_B, XMg_mica_B, XAl_mica_B, sumDiv_mica_B] = ...
    calcMicaSiteFractions(Fe2_mica_B, Mg_mica, AlVI_mica);

Kideal_B = (XFe_grt.^3 .* XMg_mica_B.^3) ./ ...
    (XMg_grt.^3 .* XFe_mica_B.^3);
lnKideal_B = log(Kideal_B);

denB = 1 + 0.0098 .* (R_J .* lnKideal_B + (Fea - Mga));
numB = 2064.7 ...
    + P_kbar .* (-0.7 - 9.8 .* (Feb - Mgb)) ...
    - 0.0098 .* (Fec - Mgc) ...
    - 7077.9 .* (XFe_mica_B - XMg_mica_B) ...
    - 941.7 .* XAl_mica_B;

% Direct calculation intentionally retains NaN and Inf results.
TB_K = numB ./ denB;
TB_C = TB_K - 273.15;

% --- Pack pressure and constant ---
row.P_kbar = P_kbar;
row.P_bar = P_bar;
row.R_J = repmat(R_J, nP, 1);

% --- Pack garnet inputs ---
row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.Mg_grt = Mg_grt;
row.Mn_grt = Mn_grt;
row.Ca_grt = Ca_grt;
row.Al_grt = Al_grt;
row.Si_grt = Si_grt;

% --- Pack mica inputs ---
row.Fe2_mica_input = Fe2_mica_input;
row.Fe3_mica_input = Fe3_mica_input;
row.FeTotal_mica = FeTotal_mica;
row.Mg_mica = Mg_mica;
row.Mn_mica = Mn_mica;
row.Ca_mica = Ca_mica;
row.Al_mica = Al_mica;
row.Si_mica = Si_mica;
row.AlIV_mica = AlIV_mica;
row.AlVI_mica = AlVI_mica;
row.K_mica = K_mica;
row.Na_mica = Na_mica;
row.Ti_mica = Ti_mica;

% --- Pack garnet compositional and activity terms ---
row.sumDiv_grt = sumDiv_grt;
row.XFe_grt = XFe_grt;
row.XMg_grt = XMg_grt;
row.XCa_grt = XCa_grt;
row.XMn_grt = XMn_grt;

row.Fea = Fea;
row.Feb = Feb;
row.Fec = Fec;
row.Mga = Mga;
row.Mgb = Mgb;
row.Mgc = Mgc;

% --- Pack Model A terms and result ---
row.Fe2_mica_A = Fe2_mica_A;
row.sumDiv_mica_A = sumDiv_mica_A;
row.XFe_mica_A = XFe_mica_A;
row.XMg_mica_A = XMg_mica_A;
row.XAl_mica_A = XAl_mica_A;
row.Kideal_A = Kideal_A;
row.lnKideal_A = lnKideal_A;
row.numA = numA;
row.denA = denA;
row.TA_K = TA_K;
row.TA_C = TA_C;

% --- Pack Model B terms and result ---
row.Fe2_mica_B = Fe2_mica_B;
row.sumDiv_mica_B = sumDiv_mica_B;
row.XFe_mica_B = XFe_mica_B;
row.XMg_mica_B = XMg_mica_B;
row.XAl_mica_B = XAl_mica_B;
row.Kideal_B = Kideal_B;
row.lnKideal_B = lnKideal_B;
row.numB = numB;
row.denB = denB;
row.TB_K = TB_K;
row.TB_C = TB_C;

end

function mineral = prepareGarnetRow(dataTable)
% prepareGarnetRow
% Extract one-row garnet cation data. Missing optional variables are set to
% zero; existing NaN values are retained and never converted to zero.

if height(dataTable) ~= 1
    error('Garnet input must be a 1-row table.');
end

mineral = struct();

mineral.Fe2 = getVarOrError(dataTable, 'Fe_cation_apfu', 'Garnet');
mineral.Mg = getVarOrError(dataTable, 'Mg_cation_apfu', 'Garnet');

mineral.Fe3 = getVarOrZero(dataTable, 'Fe3_cation_apfu', 'Garnet');
mineral.Mn = getVarOrZero(dataTable, 'Mn_cation_apfu', 'Garnet');
mineral.Ca = getVarOrZero(dataTable, 'Ca_cation_apfu', 'Garnet');
mineral.Al = getVarOrZero(dataTable, 'Al_cation_apfu', 'Garnet');
mineral.Si = getVarOrZero(dataTable, 'Si_cation_apfu', 'Garnet');

end

function mineral = prepareMicaRow(dataTable)
% prepareMicaRow
% Extract one-row mica cation data and calculate AlVI. Missing optional
% variables are set to zero; existing NaN values are retained.

if height(dataTable) ~= 1
    error('Mica input must be a 1-row table.');
end

mineral = struct();

mineral.Fe2Input = getVarOrError( ...
    dataTable, 'Fe_cation_apfu', 'Mica');
mineral.Mg = getVarOrError(dataTable, 'Mg_cation_apfu', 'Mica');
mineral.AlTotal = getVarOrError( ...
    dataTable, 'Al_cation_apfu', 'Mica');
mineral.Si = getVarOrError(dataTable, 'Si_cation_apfu', 'Mica');

mineral.Fe3Input = getVarOrZero( ...
    dataTable, 'Fe3_cation_apfu', 'Mica');
mineral.Mn = getVarOrZero(dataTable, 'Mn_cation_apfu', 'Mica');
mineral.Ca = getVarOrZero(dataTable, 'Ca_cation_apfu', 'Mica');
mineral.Ti = getVarOrZero(dataTable, 'Ti_cation_apfu', 'Mica');
mineral.K = getVarOrZero(dataTable, 'K_cation_apfu', 'Mica');
mineral.Na = getVarOrZero(dataTable, 'Na_cation_apfu', 'Mica');

mineral.FeTotal = mineral.Fe2Input + mineral.Fe3Input;

% Avoid max(NaN, 0), which can omit NaN in some MATLAB usages. This form
% explicitly preserves NaN from the Si analysis.
mineral.AlIV = 4 - mineral.Si;
if ~isnan(mineral.AlIV)
    mineral.AlIV = max(mineral.AlIV, 0);
end

mineral.AlVI = mineral.AlTotal - mineral.AlIV;

if isfinite(mineral.AlVI) && mineral.AlVI < 0
    error(['WuZhao2006: calculated Mica.AlVI is negative. Check ' ...
           'Al_cation_apfu and Si_cation_apfu (11-oxygen basis).']);
end

end

function [XFe, XMg, XAl, denom] = calcMicaSiteFractions(Fe2, Mg, AlVI)
% calcMicaSiteFractions
% Compute Fe2+-Mg-AlVI ternary muscovite fractions. Direct element-wise
% arithmetic preserves NaN and Inf for the caller to report.

denom = Fe2 + Mg + AlVI;
XFe = Fe2 ./ denom;
XMg = Mg ./ denom;
XAl = AlVI ./ denom;

end

function [Fea, Feb, Fec, Mga, Mgb, Mgc] = calcGarnetActivityTerms(XFe, XMg, XCa, XMn)
% Garnet activity polynomial terms from Wu et al. (2004b) Appendix A4-A9.

Fea = ...
      5.993 .* (XMg).^2 ...
    - 9.67 .* (XCa).^2 ...
    - 11.006 .* XFe .* XMg ...
    + 0.674 .* XFe .* XCa ...
    - 9.5725 .* XMg .* XCa ...
    - 22.765 .* XMg .* XMn ...
    - 4.6665 .* XCa .* XMn ...
    + 11.006 .* (XFe).^2 .* XMg ...
    - 0.674 .* (XFe).^2 .* XCa ...
    - 11.986 .* (XMg).^2 .* XFe ...
    + 37.966 .* (XMg).^2 .* XCa ...
    + 46.02 .* (XMg).^2 .* XMn ...
    + 19.34 .* (XCa).^2 .* XFe ...
    + 25.746 .* (XCa).^2 .* XMg ...
    + 46.02 .* (XMn).^2 .* XMg ...
    + 19.145 .* XFe .* XMg .* XCa ...
    + 45.53 .* XFe .* XMg .* XMn ...
    + 9.333 .* XFe .* XCa .* XMn ...
    + 77.876 .* XMg .* XCa .* XMn;

Feb = ...
    - 0.034 .* (XMg).^2 ...
    + 0.135 .* (XCa).^2 ...
    + 0.024 .* (XMn).^2 ...
    + 0.1 .* XFe .* XMg ...
    + 0.08 .* XFe .* XCa ...
    + 0.048 .* XFe .* XMn ...
    - 0.0515 .* XMg .* XCa ...
    + 0.07 .* XMg .* XMn ...
    + 0.1765 .* XCa .* XMn ...
    - 0.1 .* (XFe).^2 .* XMg ...
    - 0.08 .* (XFe).^2 .* XCa ...
    - 0.048 .* (XFe).^2 .* XMn ...
    + 0.068 .* (XMg).^2 .* XFe ...
    - 0.136 .* (XMg).^2 .* XCa ...
    - 0.076 .* (XMg).^2 .* XMn ...
    - 0.27 .* (XCa).^2 .* XFe ...
    - 0.28 .* (XCa).^2 .* XMg ...
    - 0.13 .* (XCa).^2 .* XMn ...
    - 0.048 .* (XMn).^2 .* XFe ...
    - 0.076 .* (XMn).^2 .* XMg ...
    - 0.13 .* (XMn).^2 .* XCa ...
    + 0.103 .* XFe .* XMg .* XCa ...
    - 0.14 .* XFe .* XMg .* XMn ...
    - 0.353 .* XFe .* XCa .* XMn ...
    - 0.414 .* XMg .* XCa .* XMn;

Fec = ...
    - 5672.0 .* (XMg).^2 ...
    + 19932.0 .* (XCa).^2 ...
    + 1617.0 .* (XMn).^2 ...
    + 23244.0 .* XFe .* XMg ...
    - 2608.0 .* XFe .* XCa ...
    + 3234.0 .* XFe .* XMn ...
    + 31326.5 .* XMg .* XCa ...
    + 45841.0 .* XMg .* XMn ...
    + 12356.0 .* XCa .* XMn ...
    - 23244.0 .* (XFe).^2 .* XMg ...
    + 2608.0 .* (XFe).^2 .* XCa ...
    - 3234.0 .* (XFe).^2 .* XMn ...
    + 11344.0 .* (XMg).^2 .* XFe ...
    - 132228.0 .* (XMg).^2 .* XCa ...
    - 82498.0 .* (XMg).^2 .* XMn ...
    - 39864.0 .* (XCa).^2 .* XFe ...
    - 51518.0 .* (XCa).^2 .* XMg ...
    - 2850.0 .* (XCa).^2 .* XMn ...
    - 3234.0 .* (XMn).^2 .* XFe ...
    - 82498.0 .* (XMn).^2 .* XMg ...
    - 2850.0 .* (XMn).^2 .* XCa ...
    - 62653.0 .* XFe .* XMg .* XCa ...
    - 91682.0 .* XFe .* XMg .* XMn ...
    - 24712.0 .* XFe .* XCa .* XMn ...
    - 177221.0 .* XMg .* XCa .* XMn;

Mga = ...
    - 5.503 .* (XFe).^2 ...
    - 12.873 .* (XCa).^2 ...
    - 23.01 .* (XMn).^2 ...
    + 11.986 .* XFe .* XMg ...
    - 9.5725 .* XFe .* XCa ...
    - 22.765 .* XFe .* XMn ...
    - 37.966 .* XMg .* XCa ...
    - 46.02 .* XMg .* XMn ...
    - 38.938 .* XCa .* XMn ...
    + 11.006 .* (XFe).^2 .* XMg ...
    - 0.674 .* (XFe).^2 .* XCa ...
    - 11.986 .* (XMg).^2 .* XFe ...
    + 37.966 .* (XMg).^2 .* XCa ...
    + 46.02 .* (XMg).^2 .* XMn ...
    + 19.34 .* (XCa).^2 .* XFe ...
    + 25.746 .* (XCa).^2 .* XMg ...
    + 46.02 .* (XMn).^2 .* XMg ...
    + 19.145 .* XFe .* XMg .* XCa ...
    + 45.53 .* XFe .* XMg .* XMn ...
    + 9.333 .* XFe .* XCa .* XMn ...
    + 77.876 .* XMg .* XCa .* XMn;

Mgb = ...
      0.05 .* (XFe).^2 ...
    + 0.14 .* (XCa).^2 ...
    + 0.038 .* (XMn).^2 ...
    - 0.068 .* XFe .* XMg ...
    - 0.0515 .* XFe .* XCa ...
    + 0.07 .* XFe .* XMn ...
    + 0.136 .* XMg .* XCa ...
    + 0.076 .* XMg .* XMn ...
    + 0.207 .* XCa .* XMn ...
    - 0.1 .* (XFe).^2 .* XMg ...
    - 0.08 .* (XFe).^2 .* XCa ...
    - 0.048 .* (XFe).^2 .* XMn ...
    + 0.068 .* (XMg).^2 .* XFe ...
    - 0.136 .* (XMg).^2 .* XCa ...
    - 0.076 .* (XMg).^2 .* XMn ...
    - 0.27 .* (XCa).^2 .* XFe ...
    - 0.28 .* (XCa).^2 .* XMg ...
    - 0.13 .* (XCa).^2 .* XMn ...
    - 0.048 .* (XMn).^2 .* XFe ...
    - 0.076 .* (XMn).^2 .* XMg ...
    - 0.13 .* (XMn).^2 .* XCa ...
    + 0.103 .* XFe .* XMg .* XCa ...
    - 0.14 .* XFe .* XMg .* XMn ...
    - 0.353 .* XFe .* XCa .* XMn ...
    - 0.414 .* XMg .* XCa .* XMn;

Mgc = ...
    11622.0 .* (XFe).^2 ...
    + 25759.0 .* (XCa).^2 ...
    + 41249.0 .* (XMn).^2 ...
    - 11344.0 .* XFe .* XMg ...
    + 31326.5 .* XFe .* XCa ...
    + 45841.0 .* XFe .* XMn ...
    + 132228.0 .* XMg .* XCa ...
    + 82498.0 .* XMg .* XMn ...
    + 88610.5 .* XCa .* XMn ...
    - 23244.0 .* (XFe).^2 .* XMg ...
    + 2608.0 .* (XFe).^2 .* XCa ...
    - 3234.0 .* (XFe).^2 .* XMn ...
    + 11344.0 .* (XMg).^2 .* XFe ...
    - 132228.0 .* (XMg).^2 .* XCa ...
    - 82498.0 .* (XMg).^2 .* XMn ...
    - 39864.0 .* (XCa).^2 .* XFe ...
    - 51518.0 .* (XCa).^2 .* XMg ...
    - 2850.0 .* (XCa).^2 .* XMn ...
    - 3234.0 .* (XMn).^2 .* XFe ...
    - 82498.0 .* (XMn).^2 .* XMg ...
    - 2850.0 .* (XMn).^2 .* XCa ...
    - 62653.0 .* XFe .* XMg .* XCa ...
    - 91682.0 .* XFe .* XMg .* XMn ...
    - 24712.0 .* XFe .* XCa .* XMn ...
    - 177221.0 .* XMg .* XCa .* XMn;

end

function printTemperatureRangeMessage(temperature_degC, modelLabel, ...
        calibrationT_min_degC, calibrationT_max_degC, ...
        selectedCode_grt, selectedCode_mica)
% printTemperatureRangeMessage
% Print a non-stopping warning when a finite calculated temperature lies
% outside the empirical calibration range of Wu and Zhao (2006).

finiteTemperature = isfinite(temperature_degC);
temperatureOutsideCalibration = finiteTemperature & ...
    (temperature_degC < calibrationT_min_degC | ...
     temperature_degC > calibrationT_max_degC);

if any(temperatureOutsideCalibration)
    finiteValues = temperature_degC(finiteTemperature);
    fprintf(2, ...
        ['WARNING: Calculated %s temperature is outside the empirical ' ...
         'calibration range of Wu and Zhao (2006): 450-760 degreeC. ' ...
         '%d of %d finite temperature point(s) are outside the range; ' ...
         'calculated finite range = %.4g-%.4g degreeC for %s & %s.\n'], ...
        modelLabel, ...
        sum(temperatureOutsideCalibration), ...
        sum(finiteTemperature), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

end

function printCompositionRangeMessages(row, selectedCode_grt, selectedCode_mica)
% printCompositionRangeMessages
% Print non-stopping warnings when finite mineral compositions fall outside
% the ranges represented by the Wu and Zhao (2006) calibration dataset.

XFe_grt = row.XFe_grt(1);
XMg_grt = row.XMg_grt(1);
XCa_grt = row.XCa_grt(1);
Fe_mica = row.FeTotal_mica(1);
Mg_mica = row.Mg_mica(1);
AlVI_mica = row.AlVI_mica(1);

if isfinite(XFe_grt) && (XFe_grt < 0.53 || XFe_grt > 0.81)
    fprintf(2, ...
        ['WARNING: Garnet XFe = %.4g is outside the Wu and Zhao (2006) ' ...
         'calibration-composition range 0.53-0.81 for %s & %s.\n'], ...
        XFe_grt, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(XMg_grt) && (XMg_grt < 0.05 || XMg_grt > 0.24)
    fprintf(2, ...
        ['WARNING: Garnet XMg = %.4g is outside the Wu and Zhao (2006) ' ...
         'calibration-composition range 0.05-0.24 for %s & %s.\n'], ...
        XMg_grt, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(XCa_grt) && (XCa_grt < 0.03 || XCa_grt > 0.23)
    fprintf(2, ...
        ['WARNING: Garnet XCa = %.4g is outside the Wu and Zhao (2006) ' ...
         'calibration-composition range 0.03-0.23 for %s & %s.\n'], ...
        XCa_grt, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(Fe_mica) && (Fe_mica < 0.04 || Fe_mica > 0.16)
    fprintf(2, ...
        ['WARNING: Muscovite total Fe = %.4g apfu is outside the ' ...
         'Wu and Zhao (2006) calibration-composition range ' ...
         '0.04-0.16 apfu (11-oxygen basis) for %s & %s. ' ...
         'The authors explicitly do not advocate application when ' ...
         'muscovite Fe < 0.04 apfu (p. 2366).\n'], ...
        Fe_mica, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(Mg_mica) && (Mg_mica < 0.04 || Mg_mica > 0.13)
    fprintf(2, ...
        ['WARNING: Muscovite Mg = %.4g apfu is outside the Wu and Zhao ' ...
         '(2006) calibration-composition range 0.04-0.13 apfu ' ...
         '(11-oxygen basis) for %s & %s. The authors explicitly do not ' ...
         'advocate application when muscovite Mg > 0.13 apfu ' ...
         '(p. 2366).\n'], ...
        Mg_mica, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(AlVI_mica) && (AlVI_mica < 1.74 || AlVI_mica > 1.96)
    fprintf(2, ...
        ['WARNING: Muscovite AlVI = %.4g apfu is outside the Wu and Zhao ' ...
         '(2006) calibration-composition range 1.74-1.96 apfu ' ...
         '(11-oxygen basis) for %s & %s.\n'], ...
        AlVI_mica, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

end

function value = getVarOrError(dataTable, variableName, mineralLabel)
% getVarOrError
% Read a required scalar cation value. NaN and zero are retained; negative,
% infinite, nonnumeric, nonreal, and nonscalar values are rejected.

if ~ismember(variableName, dataTable.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = dataTable.(variableName);

if isInvalidCationValue(value)
    error('%s contains an invalid value for %s.', mineralLabel, variableName);
end

end

function value = getVarOrZero(dataTable, variableName, mineralLabel)
% getVarOrZero
% Read an optional scalar cation value. A missing variable is set to zero,
% whereas an existing NaN value remains NaN.

if ismember(variableName, dataTable.Properties.VariableNames)
    value = dataTable.(variableName);

    if isInvalidCationValue(value)
        error('%s contains an invalid value for %s.', ...
            mineralLabel, variableName);
    end
else
    value = 0;
end

end
