function results = KroghRaheim1978(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/KroghRfiheim1978.m
% Tested with MATLAB R2024b
%
% Empirical garnet-phengite Fe-Mg exchange thermometer
% Krogh, E.J. and Raheim, A. (1978; original spelling: Råheim)
% Contributions to Mineralogy and Petrology, 66, 75-80
% DOI: https://doi.org/10.1007/BF00376087
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% (phengite) analysis and calculates temperature using the pressure-
% dependent Fe-Mg exchange thermometer of Krogh & Raheim (1978).
%
% The function accepts pressure as either a scalar or a vector. Therefore,
% it can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for every pressure
% value for each user-selected Garnet-Mica pair.
%
% -------------------------------------------------------------------------
% CALIBRATION BASIS, RANGE, AND APPLICATION NOTES
%
% Krogh & Raheim (1978) combined a small experimental dataset with one
% natural-rock pressure anchor. The calibration basis is:
%
%   Experimental T : 700, 900, and 1000 degreeC
%   Experimental P : 30 kbar for all three usable experiments
%   Natural anchor : 670 +/- 20 degreeC and 11 +/- 1 kbar
%   Main rock type : low-Fe3+ basaltic eclogite
%   Experimental KD: approximately 3.3, 4.9, and 14.1
%
% The experimental temperatures and pressure are reported in Table 1 and
% the associated text on p. 76. A phengite-bearing 800 degreeC experiment
% was not used because no useful phengite analysis was obtained (p. 76).
% The lower-pressure constraint comes from one natural Tasmanian eclogite,
% not from a lower-pressure experiment (pp. 76-77).
%
% Consequently, 670-1000 degreeC and 11-30 kbar describe the envelope of
% the combined calibration anchors, but they are NOT a formally established
% rectangular experimental calibration range. In particular, pressure
% dependence below 30 kbar was inferred from the single natural anchor.
% Krogh & Raheim also applied the relation outside this envelope, including
% a non-basaltic garnet-phengite schist at 375 degreeC and 2.5 kbar; this
% was presented as a promising example, not as an experimental extension
% of the calibration (p. 79).
%
% Important application cautions stated in the original paper:
% - The calibration is preliminary and requires additional experimental
%   studies for refinement (p. 80).
% - The thermometer is particularly applicable to metamorphic rocks of
%   basaltic composition and low Fe3+ content (abstract, p. 75; conclusion,
%   p. 80).
% - All Fe was treated as Fe2+. If Fe3+ is present in phengite, the
%   calculated KD is too low and the calculated temperature at a given
%   pressure is too high (p. 77).
% - Estimating Fe3+ in mica from electron-microprobe stoichiometry is
%   unreliable because mica has vacant cation sites and variable H2O
%   content (p. 77).
% - Non-basaltic systems such as metapelites require much greater care.
%   More natural data and preferably separate experimental calibrations for
%   different bulk compositions are required before extended use can be
%   accepted with confidence (p. 79).
% - Experimental phengites were very fine grained and commonly affected by
%   neighbouring clinopyroxene during analysis; their Mg values were partly
%   estimated by graphical extrapolation (pp. 75-76).
% - Core and rim compositions may record different stages of a prograde
%   P-T history. Garnet and phengite compositions must represent the same
%   equilibrium stage and must be screened for later recrystallization or
%   homogenization (Table 2, p. 78; discussion, p. 80).
%
% This implementation therefore issues non-stopping fprintf messages when:
%   1) input pressure is outside the 11-30 kbar combined calibration-anchor
%      envelope,
%   2) a finite temperature is outside the 670-1000 degreeC combined
%      calibration-anchor envelope,
%   3) a required calculation input is NaN, or
%   4) a calculated temperature is NaN or Inf.
%
% A range message indicates extrapolation beyond the calibration anchors;
% it does not assert that the result is automatically invalid.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Mica   : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The following normalized cation
% variables are required because they enter the thermometer directly:
%
%   Garnet table variables:
%     Fe_cation_apfu
%     Mg_cation_apfu
%
%   Mica table variables:
%     Fe_cation_apfu
%     Mg_cation_apfu
%
% Optional calculation variable in both tables:
%   Fe3_cation_apfu
%
% In this implementation, Fe_used is calculated as:
%   Fe_used = Fe_cation_apfu + Fe3_cation_apfu
%
% when Fe3_cation_apfu is present. This preserves the original project's
% total-Fe-as-ferrous convention and reproduces the paper's simplifying
% treatment only when Fe_cation_apfu stores Fe2+ and Fe3_cation_apfu stores
% a separate, non-overlapping ferric component. If Fe_cation_apfu already
% represents total Fe, Fe3_cation_apfu must not be added again. If measured
% Fe3+ is known, users should assess whether treating it as exchangeable
% Fe2+ is appropriate in light of the paper's warning on p. 77.
%
% Optional output-only variables retained when present:
%   Mn_cation_apfu, Ca_cation_apfu, Ti_cation_apfu,
%   Al_cation_apfu, Si_cation_apfu, K_cation_apfu,
%   Na_cation_apfu
%
% All finite mineral-composition values used in the calculation must be
% greater than or equal to zero. Negative finite values and Inf stop the
% calculation with an error. Zero is allowed, although it may generate a
% non-finite mathematical result that is retained and reported. NaN is
% never replaced by zero: it is propagated and reported by non-stopping
% fprintf messages. A missing Fe3_cation_apfu column is interpreted as no
% separately stored ferric component and therefore contributes zero; a
% present NaN remains NaN.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Fe-Mg ratios:
%   FeMg_Grt = Fe_used_Grt / Mg_Grt
%   FeMg_Ph  = Fe_used_Ph  / Mg_Ph
%
% Distribution coefficient:
%   KD = FeMg_Grt / FeMg_Ph
%
% Krogh & Raheim (1978) thermometer:
%   T(K) = (3685 + 77.1 * P(kbar)) / (ln(KD) + 3.52)
%   T(degreeC) = T(K) - 273.15
%
% The denominator must be finite and greater than zero. Otherwise, the
% temperature is retained as NaN and reported without stopping the loop.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = KroghRfiheim1978(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Mica tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Mica pair
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('KroghRfiheim1978 requires (rawdata_struct, P_kbar).');
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
% Extract the required tables and confirm that calculation variables exist.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_grt = rawdata_struct.Garnet;
dataset_mica = rawdata_struct.Mica;

requiredVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
validateRequiredVariables(dataset_grt, requiredVariables, 'Garnet');
validateRequiredVariables(dataset_mica, requiredVariables, 'Mica');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated table concatenation inside the loop is avoided because it forces
% MATLAB to repeatedly reallocate and copy the complete results table.
%
% The cell buffer is preallocated and doubled only when necessary. After the
% interactive loop finishes, all blocks are concatenated once with vertcat.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Envelope of the experimental and natural calibration anchors. These are
% caution limits rather than a formal rectangular experimental range.
calibrationT_min_degC = 670;
calibrationT_max_degC = 1000;
calibrationP_min_kbar = 11;
calibrationP_max_kbar = 30;

% Pressure is common to every selected mineral pair in this function call,
% so its warning is printed only once, after the first calculation.
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
    % Assumption: the first table column stores the displayed identifier.
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
    % Garnet and mica are selected independently; row indices need not match.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_mica = dataset_mica(selectedIdx_mica, :);

    % NaN is reported but deliberately permitted to propagate.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_mica);

    % Negative finite values and Inf are prohibited only in variables that
    % enter the thermometer. Zero and NaN are permitted at this stage.
    validateNonNegativeInputs(selectedData_grt, selectedData_mica);

    row = calcTemp(selectedData_grt, selectedData_mica, P_kbar);

    % Store the selected identifiers in every pressure row.
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_mica = repmat(string(selectedCode_mica), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_mica'}, 'Before', 1);

    % Store this result as one block. Buffer growth is occasional rather than
    % occurring after every selected pair.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': Krogh-Raheim 1978 = ' num2str(row.T_KR1978_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': Krogh-Raheim 1978 = ' num2str(row.T_KR1978_C(1)) ' to ' ...
            num2str(row.T_KR1978_C(end)) ' degreeC']);
    end

    % Warn once when pressure lies outside the combined calibration-anchor
    % envelope. The message explicitly distinguishes it from a formal
    % rectangular experimental calibration range.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['CAUTION: Input pressure is outside the combined experimental-and-' ...
             'natural calibration-anchor envelope of Krogh & Raheim (1978): ' ...
             '11-30 kbar (pp. 76-77). This is not a formal rectangular ' ...
             'experimental range, and the lower-pressure constraint is based ' ...
             'on one natural sample. %d of %d pressure point(s) are outside; ' ...
             'input range = %.4g-%.4g kbar. Calculation was continued.\n'], ...
            sum(pressureOutsideCalibration), numel(P_kbar), ...
            min(P_kbar), max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite result lies outside the combined calibration-
    % anchor temperature envelope. NaN and Inf are handled separately.
    finiteTemperature = isfinite(row.T_KR1978_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_KR1978_C < calibrationT_min_degC | ...
         row.T_KR1978_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_KR1978_C(finiteTemperature);
        fprintf(2, ...
            ['CAUTION: Calculated temperature is outside the combined ' ...
             'calibration-anchor envelope of Krogh & Raheim (1978): ' ...
             '670-1000 degreeC (pp. 76-77). This is not a formal rectangular ' ...
             'experimental range; the paper also presents a 375 degreeC, ' ...
             '2.5 kbar application example (p. 79). %d of %d finite point(s) ' ...
             'are outside; calculated finite range = %.4g-%.4g degreeC for ' ...
             '%s & %s. Calculation was continued.\n'], ...
            sum(temperatureOutsideCalibration), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)));
    end

    % Report NaN inputs without stopping or replacing them with zero.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; it was not replaced by zero.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report NaN/Inf caused by missing values, zeros, or an
    % invalid logarithm/denominator.
    invalidTemperature = ~isfinite(row.T_KR1978_C);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)), ...
            sum(invalidTemperature), numel(row.T_KR1978_C), ...
            sum(isnan(row.T_KR1978_C)), sum(isinf(row.T_KR1978_C)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'KroghRfiheim1978', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered blocks only once. Return an empty table if the user
% canceled before completing any calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(dataset, requiredVariables, mineralLabel)
% validateRequiredVariables
% Confirm that required calculation columns exist and contain numeric data.

missingVariables = requiredVariables( ...
    ~ismember(requiredVariables, dataset.Properties.VariableNames));

if ~isempty(missingVariables)
    error('%s table is missing required variable(s): %s.', ...
        mineralLabel, char(strjoin(string(missingVariables), ', ')));
end

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    if ~isnumeric(dataset.(variableName))
        error('%s.%s must be numeric.', mineralLabel, variableName);
    end
end

if ismember('Fe3_cation_apfu', dataset.Properties.VariableNames) && ...
        ~isnumeric(dataset.Fe3_cation_apfu)
    error('%s.Fe3_cation_apfu must be numeric when present.', mineralLabel);
end

end

function nanInputNames = findNaNInputs(data_garnet, data_mica)
% findNaNInputs
% Return calculation input names containing NaN. The output buffer is
% preallocated and NaN is never changed to zero.

calculationVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu'};
maxNames = 2 * numel(calculationVariables);
nameBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(calculationVariables)
    variableName = calculationVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if any(isnan(variableValue(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(calculationVariables)
    variableName = calculationVariables{i};
    if ismember(variableName, data_mica.Properties.VariableNames)
        variableValue = data_mica.(variableName);
        if any(isnan(variableValue(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Mica." + string(variableName);
        end
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_garnet, data_mica)
% validateNonNegativeInputs
% Reject negative finite values and Inf in calculation variables. Zero and
% NaN are allowed so that resulting non-finite values can be retained and
% reported by the caller.

calculationVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu'};
maxNames = 2 * numel(calculationVariables);
negativeNameBuffer = strings(maxNames, 1);
infiniteNameBuffer = strings(maxNames, 1);
nNegative = 0;
nInfinite = 0;

for i = 1:numel(calculationVariables)
    variableName = calculationVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nNegative = nNegative + 1;
            negativeNameBuffer(nNegative) = "Garnet." + string(variableName);
        end
        if any(isinf(variableValue(:)))
            nInfinite = nInfinite + 1;
            infiniteNameBuffer(nInfinite) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(calculationVariables)
    variableName = calculationVariables{i};
    if ismember(variableName, data_mica.Properties.VariableNames)
        variableValue = data_mica.(variableName);
        if any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nNegative = nNegative + 1;
            negativeNameBuffer(nNegative) = "Mica." + string(variableName);
        end
        if any(isinf(variableValue(:)))
            nInfinite = nInfinite + 1;
            infiniteNameBuffer(nInfinite) = "Mica." + string(variableName);
        end
    end
end

if nNegative > 0
    invalidNames = negativeNameBuffer(1:nNegative);
    error(['KroghRfiheim1978: mineral-composition values used by the ' ...
           'thermometer must be >= 0. Negative finite value(s) were found ' ...
           'in: ' char(strjoin(invalidNames, ', ')) '.']);
end

if nInfinite > 0
    invalidNames = infiniteNameBuffer(1:nInfinite);
    error(['KroghRfiheim1978: Inf is not permitted in calculation ' ...
           'inputs. Inf was found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_garnet, data_mica, P_kbar)
% calcTemp
% Calculate the Krogh & Raheim (1978) thermometer for one Garnet-Mica pair
% and return one output row for each pressure value.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;

grt = prepareMineralRow(data_garnet, 'Garnet');
mica = prepareMineralRow(data_mica, 'Mica');

% Composition terms are scalar for the selected pair. A present NaN is
% retained through every following expression.
FeMg_grt_scalar = grt.FeUsed ./ grt.Mg;
FeMg_mica_scalar = mica.FeUsed ./ mica.Mg;
KD_FeMg_scalar = FeMg_grt_scalar ./ FeMg_mica_scalar;
lnKD_scalar = log(KD_FeMg_scalar);
denominator_scalar = lnKD_scalar + 3.52;

numerator = 3685 + 77.1 .* P_kbar;

% Preserve the original implementation's physical denominator guard while
% allowing all composition arithmetic, including NaN propagation, to occur.
if isfinite(denominator_scalar) && denominator_scalar > 0
    T_K = numerator ./ denominator_scalar;
    T_C = T_K - 273.15;
else
    T_K = NaN(nP, 1);
    T_C = NaN(nP, 1);
end

% Replicate pressure-independent composition outputs to vector length.
row.Fe2_grt = repmat(grt.Fe2, nP, 1);
row.Fe3_grt = repmat(grt.Fe3, nP, 1);
row.FeUsed_grt = repmat(grt.FeUsed, nP, 1);
row.Mg_grt = repmat(grt.Mg, nP, 1);
row.Mn_grt = repmat(grt.Mn, nP, 1);
row.Ca_grt = repmat(grt.Ca, nP, 1);
row.Al_grt = repmat(grt.Al, nP, 1);
row.Si_grt = repmat(grt.Si, nP, 1);

row.Fe2_mica = repmat(mica.Fe2, nP, 1);
row.Fe3_mica = repmat(mica.Fe3, nP, 1);
row.FeUsed_mica = repmat(mica.FeUsed, nP, 1);
row.Mg_mica = repmat(mica.Mg, nP, 1);
row.Mn_mica = repmat(mica.Mn, nP, 1);
row.Ca_mica = repmat(mica.Ca, nP, 1);
row.Ti_mica = repmat(mica.Ti, nP, 1);
row.Al_mica = repmat(mica.Al, nP, 1);
row.Si_mica = repmat(mica.Si, nP, 1);
row.K_mica = repmat(mica.K, nP, 1);
row.Na_mica = repmat(mica.Na, nP, 1);

row.FeMg_grt = repmat(FeMg_grt_scalar, nP, 1);
row.FeMg_mica = repmat(FeMg_mica_scalar, nP, 1);
row.KD_FeMg = repmat(KD_FeMg_scalar, nP, 1);
row.lnKD = repmat(lnKD_scalar, nP, 1);
row.numerator_K = numerator;
row.denominator = repmat(denominator_scalar, nP, 1);

row.T_KR1978_K = T_K;
row.T_KR1978_C = T_C;

end

function mineral = prepareMineralRow(dataTable, mineralLabel)
% prepareMineralRow
% Extract one-row mineral data. NaN in a present calculation variable is
% retained. Missing output-only variables are represented by NaN.

if height(dataTable) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();
mineral.Fe2 = getRequiredScalar(dataTable, 'Fe_cation_apfu', mineralLabel);
mineral.Mg = getRequiredScalar(dataTable, 'Mg_cation_apfu', mineralLabel);

% A missing Fe3 column means no separately stored ferric component. A
% present NaN is retained rather than converted to zero.
mineral.Fe3 = getOptionalScalar(dataTable, 'Fe3_cation_apfu', ...
    mineralLabel, 0);
mineral.FeUsed = mineral.Fe2 + mineral.Fe3;

% The remaining fields are output-only and do not enter the thermometer.
mineral.Mn = getOptionalScalar(dataTable, 'Mn_cation_apfu', ...
    mineralLabel, NaN);
mineral.Ca = getOptionalScalar(dataTable, 'Ca_cation_apfu', ...
    mineralLabel, NaN);
mineral.Ti = getOptionalScalar(dataTable, 'Ti_cation_apfu', ...
    mineralLabel, NaN);
mineral.Al = getOptionalScalar(dataTable, 'Al_cation_apfu', ...
    mineralLabel, NaN);
mineral.Si = getOptionalScalar(dataTable, 'Si_cation_apfu', ...
    mineralLabel, NaN);
mineral.K = getOptionalScalar(dataTable, 'K_cation_apfu', ...
    mineralLabel, NaN);
mineral.Na = getOptionalScalar(dataTable, 'Na_cation_apfu', ...
    mineralLabel, NaN);

end

function value = getRequiredScalar(dataTable, variableName, mineralLabel)
% getRequiredScalar
% Return a required numeric scalar without changing NaN to zero.

value = dataTable.(variableName);
if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end

function value = getOptionalScalar(dataTable, variableName, ...
        mineralLabel, missingValue)
% getOptionalScalar
% Return a present optional scalar unchanged. Use the supplied missingValue
% only when the column itself is absent; a present NaN remains NaN.

if ismember(variableName, dataTable.Properties.VariableNames)
    value = dataTable.(variableName);
    if ~isnumeric(value) || ~isscalar(value)
        error('%s.%s must be a numeric scalar in the selected row.', ...
            mineralLabel, variableName);
    end
else
    value = missingValue;
end

end
