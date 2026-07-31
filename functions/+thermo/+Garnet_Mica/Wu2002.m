function results = Wu2002(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/Wu2002.m
% Tested with MATLAB R2024b
%
% Empirical Fe-Mg exchange thermometer between Garnet and Muscovite
% Wu, C.-M., Wang, X.-S., Yang, C.-H., Geng, Y.-S. and Liu, F.-L. (2002)
% Lithos, 62, 1-13
% DOI: https://doi.org/10.1016/S0024-4937(02)00096-8
% Erratum: Lithos, 66 (2003), 291
% Erratum DOI: https://doi.org/10.1016/S0024-4937(02)00251-7
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% analysis (selected independently by the user from tables) and calculates
% temperature using the Wu et al. (2002) garnet-muscovite Fe-Mg exchange
% thermometer.
%
% Two temperatures are calculated:
%   - Model A : no ferric iron in muscovite
%   - Model B : 50% of total iron in muscovite is ferric iron
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
% Wu et al. (2002) report the following ranges:
%
%   Direct empirical calibration : 530-700 degreeC and 3.0-14.0 kbar
%   Recommended application      : 480-700 degreeC and 3.0-14.0 kbar
%   Garnet composition           : Xalm = 0.51-0.82
%                                  Xpyr = 0.04-0.22
%                                  Xgrs = 0.03-0.24
%   Muscovite composition        : Fetot = 0.03-0.17 apfu
%                                  Mg    = 0.04-0.14 apfu
%                                  (11-oxygen basis)
%   Intended application         : equilibrated Garnet-Muscovite pairs in
%                                  metapelites
%
% The calibration and recommended application ranges are stated in the
% abstract on p. 1. The calibration sample selection and equilibration
% criteria are described on p. 3, the regression and compositional ranges
% are discussed on p. 5, and the recommended 480-700 degreeC and
% 3.0-14.0 kbar application range is reiterated in the conclusions on
% p. 11. The thermometer equations (9a and 9b) are presented on p. 5.
%
% Important application cautions stated or demonstrated in the article:
% - This is an empirical calibration based on natural metapelites, not a
%   direct experimental calibration. Calibration P-T values were obtained
%   using the Holdaway (2000) Garnet-Biotite thermometer and Holdaway (2001)
%   GASP barometer, so their systematic uncertainties can be inherited by
%   this thermometer (pp. 3-5).
% - The direct calibration range begins at 530 degreeC. Temperatures from
%   480 to <530 degreeC lie within the authors' recommended application range
%   but outside the direct calibration dataset (pp. 1 and 11).
% - Most calibration samples cluster at 550-650 degreeC, and some samples
%   were used more than once in the regression to balance temperature bins.
%   The ends of the calibration range are therefore less well constrained
%   than its central part (p. 5).
% - Model A and Model B use deliberately different assumptions about ferric
%   iron in muscovite. The authors describe the 0% and 50% Fe3+ treatments
%   as somewhat arbitrary, although the separately fitted models generally
%   give similar results (p. 5). Model B produced an anomalous metamorphic-
%   grade ordering in one tested inverted metamorphic sequence (p. 11).
% - The reported approximately +/-5 degreeC value is a propagated random
%   uncertainty under assumed analytical and pressure errors; systematic
%   error could not be evaluated because suitable experimental metapelite
%   data were unavailable. Agreement with the reference Garnet-Biotite
%   temperatures is generally within approximately +/-50 degreeC
%   (pp. 5-6).
% - Use spatially and texturally equilibrated analyses. For zoned garnet,
%   the calibration used garnet rims and matrix-mineral rim compositions in
%   contact with garnet. Retrogressed samples were excluded (p. 3).
% - The calibration pressures depend partly on the GASP barometer. The paper
%   notes special caution at low anorthite and grossular contents and cites
%   the refined limits Xan > 17% and Xgrs > 3% (p. 4).
% - The 2003 erratum corrects the printed definitions of the muscovite molar
%   fractions to Fe2+/(Fe2+ + Mg + AlVI), Mg/(Fe2+ + Mg + AlVI), and
%   AlVI/(Fe2+ + Mg + AlVI). The corrected definitions are used here
%   (Erratum, p. 291).
%
% This implementation therefore issues non-stopping messages using fprintf
% when:
%   1) input pressure is outside 3.0-14.0 kbar;
%   2) a finite calculated temperature is outside the recommended
%      480-700 degreeC application range;
%   3) a finite temperature is 480 to <530 degreeC and is therefore within
%      the recommended application range but below the direct calibration;
%   4) a finite Garnet or Muscovite composition is outside the reported
%      calibration-composition range; or
%   5) an input or calculated temperature is NaN or Inf.
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
%     Fe_cation_apfu       % measured ferrous-iron term
%     Mg_cation_apfu
%     Al_cation_apfu
%     Si_cation_apfu
%
% Optional variables used in the calculation when present:
%   Garnet : Mn_cation_apfu, Ca_cation_apfu
%   Mica   : Fe3_cation_apfu
%
% Optional variables retained in the output table when present:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Missing optional variables are represented by zero. An optional variable
% that exists and contains NaN remains NaN; it is never converted to zero.
% AlVI in mica is calculated from Al_total and Si on an 11-oxygen basis as
% AlIV = max(0, 4 - Si) and AlVI = Al_total - AlIV.
%
% All finite mineral-composition values used by the thermometer must be
% non-negative. Negative and infinite values stop the calculation. Zero is
% allowed, although it may produce a non-finite exchange coefficient or
% temperature. NaN values are retained, propagated through the calculation,
% and reported by non-stopping fprintf messages.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Wu et al., 2002, Eqs. 9a-b, p. 5)
%
% Model A, assuming no ferric iron in muscovite:
%
%   T(A) (K) =
%   [ 969.9 + P(kbar)*(1.3 - 9.1*Gb) - 0.0091*Gc
%     - 4393.8*(XFe_mus_A - XMg_mus_A) + 200.4*XAl_mus_A ]
%   / [ 1 + 0.0091*(3*R*lnKd_A + Ga) ]
%
% Model B, assuming 50% ferric iron in muscovite:
%
%   T(B) (K) =
%   [ -1167.3 - P(kbar)*(0.2 + 8.8*Gb) - 0.0088*Gc
%     - 6878.1*(XFe_mus_B - XMg_mus_B) + 2469.0*XAl_mus_B ]
%   / [ 1 + 0.0088*(3*R*lnKd_B + Ga) ]
%
% where:
%   Kd = (Fe2+/Mg)_grt / (Fe2+/Mg)_mus
%   R  = 8.3144 J mol^-1 K^-1
%
% Garnet non-ideal terms Ga, Gb, and Gc are the Holdaway-type polynomial
% expressions adopted by Wu et al. (2002).
%
% For muscovite:
%   Fetot_mus = Fe_cation_apfu + Fe3_cation_apfu (when present)
%   Model A uses Fe2_mus = Fetot_mus
%   Model B uses Fe2_mus = 0.5 * Fetot_mus
%
% For garnet, Fe2_grt = Fe_cation_apfu is used in Kd and the garnet site
% fractions. Fe3_cation_apfu, when present, is retained in the output but is
% not added to the exchange Fe term.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Wu2002(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Mica tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Mica pair. The output variable set is
%             intended to remain stable for downstream processing.
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('Wu2002 requires (rawdata_struct, P_kbar).');
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

% Limits reported by Wu et al. (2002).
calibrationT_min_degC = 530;
applicationT_min_degC = 480;
applicationT_max_degC = 700;
applicationP_min_kbar = 3;
applicationP_max_kbar = 14;

% Pressure is common to all selected mineral pairs in this function call.
% The pressure message is therefore printed only once, after the first
% completed calculation.
pressureOutsideApplication = ...
    P_kbar < applicationP_min_kbar | ...
    P_kbar > applicationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
    % Assumption: the first column stores an identifier displayed to the user.
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

    % Check all raw variables actually used by the thermometer. The
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

    % Warn once when any pressure lies outside the recommended application
    % range of 3.0-14.0 kbar. The calculation is not stopped.
    if any(pressureOutsideApplication) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the calibration and recommended ' ...
             'application range of Wu et al. (2002): 3.0-14.0 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideApplication), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Print temperature-range messages independently for Models A and B.
    printTemperatureRangeMessages(row.TA_C, 'Model A', ...
        calibrationT_min_degC, applicationT_min_degC, ...
        applicationT_max_degC, selectedCode_grt, selectedCode_mica);
    printTemperatureRangeMessages(row.TB_C, 'Model B', ...
        calibrationT_min_degC, applicationT_min_degC, ...
        applicationT_max_degC, selectedCode_grt, selectedCode_mica);

    % Report compositions outside the finite ranges represented by the
    % empirical calibration dataset. NaN inputs are handled separately.
    printCompositionRangeMessages(row, selectedCode_grt, selectedCode_mica);

    % Print a non-stopping message when an input used by the thermometer
    % contains NaN. fprintf is used so that the message remains visible even
    % when MATLAB warnings have been disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, and NaN values were not converted to zero.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report NaN/Inf results caused by missing data, zero
    % denominators, logarithms at zero, or other non-finite operations.
    printNonFiniteTemperatureMessage(row.TA_C, 'Model A', ...
        selectedCode_grt, selectedCode_mica);
    printNonFiniteTemperatureMessage(row.TB_C, 'Model B', ...
        selectedCode_grt, selectedCode_mica);

    disp('--------------------------------------------------');

    % Ask whether to repeat.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Wu2002', ...
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
    error(['Wu2002: thermometer inputs must be real, scalar, ' ...
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

function row = calcTemp(data_garnet, data_mica, P_kbar)
% calcTemp
% Compute Wu et al. (2002) Models A and B for one garnet row, one mica row,
% and every pressure supplied in P_kbar.
%
% Inputs:
%   data_garnet : 1-row table containing garnet apfu cations
%   data_mica   : 1-row table containing mica apfu cations
%   P_kbar      : pressure in kbar (column vector after input validation)
%
% Output:
%   row : table with one row per pressure value, including mineral inputs,
%         compositional terms, exchange variables, and both temperatures.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

row = table();

% --- Physical constant ---
R_J = 8.3144;

% --- Prepare one-row mineral data ---
grt = prepareGarnetRow(data_garnet);
mica = prepareMicaRow(data_mica);

% --- Replicate scalar mineral inputs for every pressure point ---
Fe2_grt = repmat(grt.Fe2, nP, 1);
Fe3_grt = repmat(grt.Fe3, nP, 1);
Mg_grt = repmat(grt.Mg, nP, 1);
Mn_grt = repmat(grt.Mn, nP, 1);
Ca_grt = repmat(grt.Ca, nP, 1);
Al_grt = repmat(grt.Al, nP, 1);
Si_grt = repmat(grt.Si, nP, 1);

Fe2_mica_input = repmat(mica.Fe2, nP, 1);
Fe3_mica_input = repmat(mica.Fe3, nP, 1);
Mg_mica = repmat(mica.Mg, nP, 1);
Al_mica = repmat(mica.Al, nP, 1);
Si_mica = repmat(mica.Si, nP, 1);
AlIV_mica = repmat(mica.AlIV, nP, 1);
AlVI_mica = repmat(mica.AlVI, nP, 1);
Ti_mica = repmat(mica.Ti, nP, 1);
K_mica = repmat(mica.K, nP, 1);
Na_mica = repmat(mica.Na, nP, 1);

% --- Garnet site fractions ---
sumDivalent_grt = Fe2_grt + Mg_grt + Ca_grt + Mn_grt;

XFe_grt = Fe2_grt ./ sumDivalent_grt;
XMg_grt = Mg_grt ./ sumDivalent_grt;
XCa_grt = Ca_grt ./ sumDivalent_grt;
XMn_grt = Mn_grt ./ sumDivalent_grt;

% --- Garnet non-ideal terms Ga, Gb, and Gc ---
Ga = ...
      12.4   .* (XFe_grt.^2) ...
    + 22.09  .* (XMg_grt.^2) ...
    - 12.02  .* (XCa_grt.^2) ...
    + 23.01  .* (XMn_grt.^2) ...
    - 68.98  .* XFe_grt .* XMg_grt ...
    + 37.33  .* XFe_grt .* XCa_grt ...
    + 18.165 .* XFe_grt .* XMn_grt ...
    + 35.33  .* XMg_grt .* XCa_grt ...
    + 27.855 .* XMg_grt .* XMn_grt ...
    + 35.165 .* XCa_grt .* XMn_grt;

Gb = ...
     -0.05   .* (XFe_grt.^2) ...
    -0.034   .* (XMg_grt.^2) ...
    -0.005   .* (XCa_grt.^2) ...
    -0.014   .* (XMn_grt.^2) ...
    +0.168   .* XFe_grt .* XMg_grt ...
    +0.1565  .* XFe_grt .* XCa_grt ...
    -0.022   .* XFe_grt .* XMn_grt ...
    -0.2125  .* XMg_grt .* XCa_grt ...
    -0.006   .* XMg_grt .* XMn_grt ...
    -0.0305  .* XCa_grt .* XMn_grt;

Gc = ...
     -22265.0  .* (XFe_grt.^2) ...
    -24166.0   .* (XMg_grt.^2) ...
    + 3220.0   .* (XCa_grt.^2) ...
    -39632.0   .* (XMn_grt.^2) ...
    +92862.0   .* XFe_grt .* XMg_grt ...
    -67328.0   .* XFe_grt .* XCa_grt ...
    -38681.5   .* XFe_grt .* XMn_grt ...
    -99262.0   .* XMg_grt .* XCa_grt ...
    -40582.5   .* XMg_grt .* XMn_grt ...
    -79669.5   .* XCa_grt .* XMn_grt;

% --- Muscovite Fe assumptions for Model A and Model B ---
FeTotal_mica = Fe2_mica_input + Fe3_mica_input;

Fe2_mica_A = FeTotal_mica;
Fe2_mica_B = 0.5 .* FeTotal_mica;

denMus_A = Fe2_mica_A + Mg_mica + AlVI_mica;
denMus_B = Fe2_mica_B + Mg_mica + AlVI_mica;

% Corrected muscovite molar-fraction definitions from the 2003 erratum.
XFe_mus_A = Fe2_mica_A ./ denMus_A;
XMg_mus_A = Mg_mica ./ denMus_A;
XAl_mus_A = AlVI_mica ./ denMus_A;

XFe_mus_B = Fe2_mica_B ./ denMus_B;
XMg_mus_B = Mg_mica ./ denMus_B;
XAl_mus_B = AlVI_mica ./ denMus_B;

% --- Fe-Mg exchange coefficients ---
Kd_A = (Fe2_grt ./ Mg_grt) ./ (Fe2_mica_A ./ Mg_mica);
Kd_B = (Fe2_grt ./ Mg_grt) ./ (Fe2_mica_B ./ Mg_mica);

lnKd_A = log(Kd_A);
lnKd_B = log(Kd_B);

% --- Solve temperature using Eqs. (9a) and (9b) ---
denTA = 1 + 0.0091 .* (3 .* R_J .* lnKd_A + Ga);
numTA = 969.9 + P_kbar .* (1.3 - 9.1 .* Gb) - 0.0091 .* Gc ...
    - 4393.8 .* (XFe_mus_A - XMg_mus_A) + 200.4 .* XAl_mus_A;

denTB = 1 + 0.0088 .* (3 .* R_J .* lnKd_B + Ga);
numTB = -1167.3 - P_kbar .* (0.2 + 8.8 .* Gb) - 0.0088 .* Gc ...
    - 6878.1 .* (XFe_mus_B - XMg_mus_B) + 2469.0 .* XAl_mus_B;

TA_K = numTA ./ denTA;
TB_K = numTB ./ denTB;
TA_C = TA_K - 273.15;
TB_C = TB_K - 273.15;

% --- Pack pressure ---
row.P_kbar = P_kbar;
row.P_bar = P_bar;

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
row.Al_mica = Al_mica;
row.Si_mica = Si_mica;
row.AlIV_mica = AlIV_mica;
row.AlVI_mica = AlVI_mica;
row.Ti_mica = Ti_mica;
row.K_mica = K_mica;
row.Na_mica = Na_mica;

% --- Pack garnet fractions and non-ideal terms ---
row.XFe_grt = XFe_grt;
row.XMg_grt = XMg_grt;
row.XCa_grt = XCa_grt;
row.XMn_grt = XMn_grt;
row.Ga = Ga;
row.Gb = Gb;
row.Gc = Gc;

% --- Pack Model A outputs ---
row.Fe2_mica_A = Fe2_mica_A;
row.XFe_mus_A = XFe_mus_A;
row.XMg_mus_A = XMg_mus_A;
row.XAl_mus_A = XAl_mus_A;
row.Kd_A = Kd_A;
row.lnKd_A = lnKd_A;
row.TA_K = TA_K;
row.TA_C = TA_C;

% --- Pack Model B outputs ---
row.Fe2_mica_B = Fe2_mica_B;
row.XFe_mus_B = XFe_mus_B;
row.XMg_mus_B = XMg_mus_B;
row.XAl_mus_B = XAl_mus_B;
row.Kd_B = Kd_B;
row.lnKd_B = lnKd_B;
row.TB_K = TB_K;
row.TB_C = TB_C;

end

function grt = prepareGarnetRow(dataTable)
% prepareGarnetRow
% Extract one-row garnet cation data. Missing optional variables are set to
% zero; existing NaN values are retained and never converted to zero.

if height(dataTable) ~= 1
    error('Garnet input must be a 1-row table.');
end

grt = struct();

grt.Fe2 = getVarOrError(dataTable, 'Fe_cation_apfu', 'Garnet');
grt.Mg = getVarOrError(dataTable, 'Mg_cation_apfu', 'Garnet');

grt.Fe3 = getVarOrZero(dataTable, 'Fe3_cation_apfu', 'Garnet');
grt.Mn = getVarOrZero(dataTable, 'Mn_cation_apfu', 'Garnet');
grt.Ca = getVarOrZero(dataTable, 'Ca_cation_apfu', 'Garnet');
grt.Ti = getVarOrZero(dataTable, 'Ti_cation_apfu', 'Garnet');
grt.Al = getVarOrZero(dataTable, 'Al_cation_apfu', 'Garnet');
grt.Si = getVarOrZero(dataTable, 'Si_cation_apfu', 'Garnet');
grt.K = getVarOrZero(dataTable, 'K_cation_apfu', 'Garnet');
grt.Na = getVarOrZero(dataTable, 'Na_cation_apfu', 'Garnet');

end

function mica = prepareMicaRow(dataTable)
% prepareMicaRow
% Extract one-row mica cation data. Missing optional variables are set to
% zero; existing NaN values are retained and never converted to zero.

if height(dataTable) ~= 1
    error('Mica input must be a 1-row table.');
end

mica = struct();

mica.Fe2 = getVarOrError(dataTable, 'Fe_cation_apfu', 'Mica');
mica.Mg = getVarOrError(dataTable, 'Mg_cation_apfu', 'Mica');
mica.Al = getVarOrError(dataTable, 'Al_cation_apfu', 'Mica');
mica.Si = getVarOrError(dataTable, 'Si_cation_apfu', 'Mica');

mica.Fe3 = getVarOrZero(dataTable, 'Fe3_cation_apfu', 'Mica');
mica.Mn = getVarOrZero(dataTable, 'Mn_cation_apfu', 'Mica');
mica.Ca = getVarOrZero(dataTable, 'Ca_cation_apfu', 'Mica');
mica.Ti = getVarOrZero(dataTable, 'Ti_cation_apfu', 'Mica');
mica.K = getVarOrZero(dataTable, 'K_cation_apfu', 'Mica');
mica.Na = getVarOrZero(dataTable, 'Na_cation_apfu', 'Mica');

% Preserve NaN explicitly rather than allowing max to treat it as missing.
if isnan(mica.Si)
    mica.AlIV = NaN;
else
    mica.AlIV = max(0, 4 - mica.Si);
end

mica.AlVI = mica.Al - mica.AlIV;

if isfinite(mica.AlVI) && mica.AlVI < 0
    error(['Wu2002: calculated Mica.AlVI is negative. Check the ' ...
           'Al and Si cation normalization (the calculation assumes an ' ...
           '11-oxygen mica basis).']);
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

function printTemperatureRangeMessages(temperature_degC, modelLabel, ...
        calibrationT_min_degC, applicationT_min_degC, ...
        applicationT_max_degC, selectedCode_grt, selectedCode_mica)
% printTemperatureRangeMessages
% Print non-stopping temperature-range messages for one thermometer model.

finiteTemperature = isfinite(temperature_degC);
outsideApplication = finiteTemperature & ...
    (temperature_degC < applicationT_min_degC | ...
     temperature_degC > applicationT_max_degC);

if any(outsideApplication)
    finiteValues = temperature_degC(finiteTemperature);
    fprintf(2, ...
        ['WARNING: %s calculated temperature is outside the recommended ' ...
         'application range of Wu et al. (2002): 480-700 degreeC. ' ...
         '%d of %d finite temperature point(s) are outside the range; ' ...
         'calculated finite range = %.4g-%.4g degreeC for %s & %s.\n'], ...
        modelLabel, ...
        sum(outsideApplication), ...
        sum(finiteTemperature), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

belowDirectCalibration = finiteTemperature & ...
    temperature_degC >= applicationT_min_degC & ...
    temperature_degC < calibrationT_min_degC;

if any(belowDirectCalibration)
    relevantValues = temperature_degC(belowDirectCalibration);
    fprintf(2, ...
        ['CAUTION: %s calculated temperature is within the authors'' ' ...
         'recommended range but below the direct empirical calibration ' ...
         'range of 530-700 degreeC. %d point(s); range = %.4g-%.4g ' ...
         'degreeC for %s & %s.\n'], ...
        modelLabel, ...
        sum(belowDirectCalibration), ...
        min(relevantValues), ...
        max(relevantValues), ...
        char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

end

function printCompositionRangeMessages(row, selectedCode_grt, selectedCode_mica)
% printCompositionRangeMessages
% Print non-stopping messages when finite mineral compositions fall outside
% the ranges represented by the Wu et al. (2002) calibration dataset.

Xalm = row.XFe_grt(1);
Xpyr = row.XMg_grt(1);
Xgrs = row.XCa_grt(1);
FeTotal_mica = row.FeTotal_mica(1);
Mg_mica = row.Mg_mica(1);

if isfinite(Xalm) && (Xalm < 0.51 || Xalm > 0.82)
    fprintf(2, ...
        ['WARNING: Garnet Xalm = %.4g is outside the Wu et al. (2002) ' ...
         'calibration-composition range 0.51-0.82 for %s & %s.\n'], ...
        Xalm, char(string(selectedCode_grt)), char(string(selectedCode_mica)));
end

if isfinite(Xpyr) && (Xpyr < 0.04 || Xpyr > 0.22)
    fprintf(2, ...
        ['WARNING: Garnet Xpyr = %.4g is outside the Wu et al. (2002) ' ...
         'calibration-composition range 0.04-0.22 for %s & %s.\n'], ...
        Xpyr, char(string(selectedCode_grt)), char(string(selectedCode_mica)));
end

if isfinite(Xgrs) && (Xgrs < 0.03 || Xgrs > 0.24)
    fprintf(2, ...
        ['WARNING: Garnet Xgrs = %.4g is outside the Wu et al. (2002) ' ...
         'calibration-composition range 0.03-0.24 for %s & %s.\n'], ...
        Xgrs, char(string(selectedCode_grt)), char(string(selectedCode_mica)));
end

if isfinite(FeTotal_mica) && (FeTotal_mica < 0.03 || FeTotal_mica > 0.17)
    fprintf(2, ...
        ['WARNING: Muscovite Fetot = %.4g apfu is outside the Wu et al. ' ...
         '(2002) calibration-composition range 0.03-0.17 apfu ' ...
         '(11-oxygen basis) for %s & %s.\n'], ...
        FeTotal_mica, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

if isfinite(Mg_mica) && (Mg_mica < 0.04 || Mg_mica > 0.14)
    fprintf(2, ...
        ['WARNING: Muscovite Mg = %.4g apfu is outside the Wu et al. ' ...
         '(2002) calibration-composition range 0.04-0.14 apfu ' ...
         '(11-oxygen basis) for %s & %s.\n'], ...
        Mg_mica, char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

end

function printNonFiniteTemperatureMessage(temperature_degC, modelLabel, ...
        selectedCode_grt, selectedCode_mica)
% printNonFiniteTemperatureMessage
% Retain and report NaN/Inf temperatures without stopping the calculation.

invalidTemperature = ~isfinite(temperature_degC);

if any(invalidTemperature)
    fprintf(2, ...
        ['WARNING: Non-finite %s temperature values were calculated for ' ...
         '%s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
        modelLabel, ...
        char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)), ...
        sum(invalidTemperature), ...
        numel(temperature_degC), ...
        sum(isnan(temperature_degC)), ...
        sum(isinf(temperature_degC)));
end

end
