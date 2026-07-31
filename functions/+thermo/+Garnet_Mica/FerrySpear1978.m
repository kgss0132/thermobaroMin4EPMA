function results = FerrySpear1978(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/FerrySpear1978.m
% Tested with MATLAB R2024b
%
% Garnet-biotite Fe-Mg exchange thermometer
% Ferry, J.M. and Spear, F.S. (1978)
% Contributions to Mineralogy and Petrology, 66, 113-117
% DOI: https://doi.org/10.1007/BF00372150
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% analysis (selected by the user from tables) and calculates temperature
% using the Ferry and Spear (1978) garnet-biotite Fe-Mg exchange
% thermometer.
%
% Two temperature expressions are returned:
%   1) T_Eq1 : the directly calibrated expression at 2.07 kbar
%   2) T_Eq7 : the pressure-corrected polybaric expression
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Garnet-Mica pair, and appends results
% into a single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Ferry and Spear (1978) calibrated Fe2+-Mg partitioning between synthetic
% garnet and synthetic biotite under the following conditions:
%
%   Temperature : 550-800 degreeC
%   Pressure    : 2.07 kbar (one experimental pressure)
%   Garnet      : principally Alm90Py10, with additional Alm80Py20 tests
%   Biotite     : annite-phlogopite Fe-Mg binary compositions
%   Redox       : oxygen fugacity controlled by a graphite-methane buffer
%
% The temperature and pressure conditions are summarized in the abstract
% on p. 113 and described in the Experimental Procedures on p. 114. The
% equilibrium brackets and fixed-pressure calibration are shown in Fig. 3
% and Eq. (1) on pp. 115-116. Equilibrium was demonstrated by reversal of
% the exchange reaction at every experimental temperature (pp. 113-115).
%
% IMPORTANT PRESSURE LIMITATION:
% All exchange experiments were performed at 2.07 kbar. Eq. (1) is therefore
% a direct calibration only at 2.07 kbar. The polybaric Eq. (7) was obtained
% by adding a reaction-volume term from thermodynamic data; it was not
% calibrated using experiments over a pressure interval (p. 116). An Eq. (7)
% calculation at any pressure other than 2.07 kbar is consequently a
% pressure extrapolation from the direct experimental condition.
%
% IMPORTANT COMPOSITION LIMITATIONS:
% The thermometer was developed for garnet and biotite close to binary
% Fe-Mg compositions (abstract, p. 113; Applications, p. 117). Ferry and
% Spear (1978) explicitly caution that Ca and Mn in garnet and Ti and
% octahedral Al in biotite affect the exchange coefficient (p. 117).
%
% Comparison with natural calibrations suggested that the thermometer may
% be useful without component corrections up to approximately:
%
%   (Ca + Mn)/(Ca + Mn + Fe + Mg) in garnet <= 0.20
%   (AlVI + Ti)/(AlVI + Ti + Fe + Mg) in biotite <= 0.15
%
% These values are empirical application guides from natural-sample
% comparisons, not direct multicomponent experimental calibration limits
% (p. 117). The additional experiments suggest ideal Fe-Mg mixing at least
% over Fe/(Fe + Mg) = 0.80-1.00 for the Fe-rich garnet compositions tested
% (p. 116).
%
% The experiments were extrapolated to lower temperature and agreed with
% the Thompson (1975) natural calibration over 400-600 degreeC, but this
% comparison does not extend the direct experimental calibration below
% 550 degreeC (p. 117).
%
% The maximum practical resolution reported by Ferry and Spear (1978) is
% approximately +/-50 degreeC (p. 117). The experimental temperature and
% pressure accuracies of +/-6 degreeC and +/-50 bar (p. 114) are apparatus
% uncertainties and must not be interpreted as the total thermometer error.
%
% The exchange reaction is formulated for Fe2+-Mg. Fe3_cation_apfu is
% retained in the output for inspection but is not added to
% Fe_cation_apfu in the mole ratios or K calculation. The selected Mica
% analysis must represent biotite/phlogopite-annite, not an arbitrary mica.
% Natural mineral pairs must represent the same equilibrium assemblage.
%
% This implementation issues non-stopping messages using fprintf when:
%   1) input pressure differs from the sole experimental pressure, 2.07 kbar,
%   2) a finite calculated T_Eq1 or T_Eq7 is outside 550-800 degreeC,
%   3) an input used by the thermometer contains NaN,
%   4) Fe2+ or Mg is zero in either mineral, or
%   5) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Mica   : table containing biotite analyses
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns should contain
% normalized cation data.
%
% Required variables in both Garnet and Mica tables:
%   Fe_cation_apfu          % Fe2+ used in the exchange calculation
%   Mg_cation_apfu
%
% Optional variables:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Existing NaN values are retained as missing values. NaN in Fe2+ or Mg is
% propagated through the calculation and produces NaN temperature. A missing
% optional variable is set to zero, but an existing optional variable whose
% value is NaN remains NaN in the output.
%
% All finite mineral-composition values must be non-negative. Negative finite
% values are prohibited. Fe2+ and Mg must additionally be greater than zero
% for the logarithmic exchange calculation. To preserve the requested
% non-negative-input policy, a zero Fe2+ or Mg value does not raise an error;
% the affected exchange and temperature results are set to NaN and reported
% by fprintf. This prevents an invalid 0 K (-273.15 degreeC) result.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Exchange coefficient:
%
%   K = (Mg/Fe2+)_garnet / (Mg/Fe2+)_biotite
%
% Eq. (1), directly calibrated at 2.07 kbar:
%
%   lnK = -2109/T(K) + 0.782
%   T_Eq1(K) = 2109/(0.782 - lnK)
%
% Eq. (7), polybaric expression:
%
%   12454 - 4.662T(K) + 0.057P(bar) + 3RT(K)lnK = 0
%   T_Eq7(K) = (12454 + 0.057P(bar))/(4.662 - 3RlnK)
%
% where:
%   R      : 1.987 cal mol^-1 K^-1
%   P_bar  : 1000 * P_kbar
%
% -------------------------------------------------------------------------
% Syntax:
%   results = FerrySpear1978(rawdata_struct, P_kbar)
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
% Pressure is converted to a column vector so that the fixed-pressure and
% pressure-range launchers use the same calculation path.
if nargin < 2
    error('FerrySpear1978 requires (rawdata_struct, P_kbar).');
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
% Extract the required tables without modifying their contents.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_grt = rawdata_struct.Garnet;
dataset_mica = rawdata_struct.Mica;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Calculation results are buffered as table blocks and concatenated once
% after the interactive loop, avoiding repeated results-table reallocation.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degC = 550;
calibrationT_max_degC = 800;
experimentalP_kbar = 2.07;
pressureComparisonTolerance_kbar = 1e-9;

% There is one direct experimental pressure rather than a calibrated
% pressure interval. Print this pressure-extrapolation warning only once.
pressureOutsideDirectCalibration = ...
    abs(P_kbar - experimentalP_kbar) > pressureComparisonTolerance_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
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
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_mica = dataset_mica(selectedIdx_mica, :);

    % NaN inputs remain NaN; their variable names are collected for the
    % non-stopping warning printed after the result.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_mica);

    % A zero Fe2+ or Mg input is mathematically invalid for log(K). It is not
    % treated as a negative-input error; it is converted to NaN in calcTemp.
    zeroExchangeInputNames = findZeroExchangeInputs( ...
        selectedData_grt, selectedData_mica);

    % Negative finite mineral compositions are prohibited. NaN and zero are
    % deliberately excluded from the error condition.
    validateNonnegativeInputs(selectedData_grt, selectedData_mica);

    row = calcTemp(selectedData_grt, selectedData_mica, P_kbar);

    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_mica = repmat(string(selectedCode_mica), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_mica'}, 'Before', 1);

    % Store the result as one block. The buffer grows geometrically only when
    % its capacity is exhausted, rather than on every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': Eq. (1) = ' num2str(row.TEq1_C) ...
            ' degreeC, Eq. (7) = ' num2str(row.TEq7_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': Eq. (1) = ' num2str(row.TEq1_C(1)) ' to ' ...
            num2str(row.TEq1_C(end)) ' degreeC, Eq. (7) = ' ...
            num2str(row.TEq7_C(1)) ' to ' num2str(row.TEq7_C(end)) ...
            ' degreeC']);
    end

    % Warn once if any input pressure differs from the sole directly tested
    % pressure. Eq. (7) calculations continue as pressure extrapolations.
    if any(pressureOutsideDirectCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Ferry and Spear (1978) directly calibrated the thermometer ' ...
             'at 2.07 kbar only. %d of %d pressure point(s) differ from 2.07 kbar; ' ...
             'input range = %.4g-%.4g kbar. Eq. (7) results at those pressures ' ...
             'use a thermodynamic pressure extrapolation.\n'], ...
            sum(pressureOutsideDirectCalibration), ...
            numel(P_kbar), min(P_kbar), max(P_kbar));
        pressureWarningIssued = true;
    end

    % Temperature-range warnings are evaluated separately for Eqs. (1) and
    % (7). NaN and Inf are handled by the non-finite warning below.
    printTemperatureCalibrationWarning( ...
        row.TEq1_C, 'Eq. (1)', calibrationT_min_degC, ...
        calibrationT_max_degC, selectedCode_grt, selectedCode_mica);
    printTemperatureCalibrationWarning( ...
        row.TEq7_C, 'Eq. (7)', calibrationT_min_degC, ...
        calibrationT_max_degC, selectedCode_grt, selectedCode_mica);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         The calculation was continued, and NaN was retained as a missing value.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    if ~isempty(zeroExchangeInputNames)
        fprintf(2, ...
            ['WARNING: Zero was found in Fe2+ or Mg used by the exchange calculation ' ...
             'for %s & %s: %s.\n' ...
             '         log(K) is undefined for this input; temperature was retained as NaN, not 0 K.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            char(strjoin(zeroExchangeInputNames, ', ')));
    end

    invalidEq1 = ~isfinite(row.TEq1_C);
    invalidEq7 = ~isfinite(row.TEq7_C);
    if any(invalidEq1) || any(invalidEq7)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(Eq. (1): %d of %d; Eq. (7): %d of %d; total NaN: %d; total Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            sum(invalidEq1), numel(row.TEq1_C), ...
            sum(invalidEq7), numel(row.TEq7_C), ...
            sum(isnan(row.TEq1_C)) + sum(isnan(row.TEq7_C)), ...
            sum(isinf(row.TEq1_C)) + sum(isinf(row.TEq7_C)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'FerrySpear1978', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered result blocks once after all selections.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_garnet, data_mica)
% findNaNInputs
% Return the names of Fe2+ and Mg calculation inputs that contain NaN.

variableNames = {'Fe_cation_apfu', 'Mg_cation_apfu'};
nameBuffer = strings(2 * numel(variableNames), 1);
nNames = 0;

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    variableValue = data_garnet.(variableName);
    if any(isnan(variableValue(:)))
        nNames = nNames + 1;
        nameBuffer(nNames) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    variableValue = data_mica.(variableName);
    if any(isnan(variableValue(:)))
        nNames = nNames + 1;
        nameBuffer(nNames) = "Mica." + string(variableName);
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function zeroInputNames = findZeroExchangeInputs(data_garnet, data_mica)
% findZeroExchangeInputs
% Return the names of finite Fe2+ and Mg calculation inputs equal to zero.

variableNames = {'Fe_cation_apfu', 'Mg_cation_apfu'};
nameBuffer = strings(2 * numel(variableNames), 1);
nNames = 0;

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    variableValue = data_garnet.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) == 0)
        nNames = nNames + 1;
        nameBuffer(nNames) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    variableValue = data_mica.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) == 0)
        nNames = nNames + 1;
        nameBuffer(nNames) = "Mica." + string(variableName);
    end
end

zeroInputNames = nameBuffer(1:nNames);

end

function validateNonnegativeInputs(data_garnet, data_mica)
% validateNonnegativeInputs
% Reject negative finite mineral compositions and Inf. NaN and zero are
% allowed; zero Fe2+ or Mg is handled as NaN by the exchange calculation.

variableNames = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu', 'Mn_cation_apfu', 'Ca_cation_apfu', ...
    'Ti_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu', ...
    'K_cation_apfu', 'Na_cation_apfu'};

nameBuffer = strings(2 * numel(variableNames), 1);
nNames = 0;

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if ~isnumeric(variableValue) || ~isscalar(variableValue) || ...
                any(isinf(variableValue(:))) || ...
                any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nNames = nNames + 1;
            nameBuffer(nNames) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(variableNames)
    variableName = variableNames{i};
    if ismember(variableName, data_mica.Properties.VariableNames)
        variableValue = data_mica.(variableName);
        if ~isnumeric(variableValue) || ~isscalar(variableValue) || ...
                any(isinf(variableValue(:))) || ...
                any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nNames = nNames + 1;
            nameBuffer(nNames) = "Mica." + string(variableName);
        end
    end
end

invalidInputNames = nameBuffer(1:nNames);
if ~isempty(invalidInputNames)
    error(['FerrySpear1978: mineral-composition values must be finite ' ...
           'non-negative scalars or NaN. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_garnet, data_mica, P_kbar)
% calcTemp
% Calculate Eqs. (1) and (7) for one garnet-biotite pair over a scalar or
% vector of pressures, returning one table row per pressure value.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);
R_cal = 1.987;

garnet = prepareMineralRow(data_garnet, 'Garnet');
mica = prepareMineralRow(data_mica, 'Mica');

% Repeat composition values to their final vector size before calculation.
Fe2_grt = repmat(garnet.Fe2, nP, 1);
Fe3_grt = repmat(garnet.Fe3, nP, 1);
FeUsed_grt = Fe2_grt;
Mg_grt = repmat(garnet.Mg, nP, 1);
Mn_grt = repmat(garnet.Mn, nP, 1);
Ca_grt = repmat(garnet.Ca, nP, 1);
Al_grt = repmat(garnet.Al, nP, 1);
Si_grt = repmat(garnet.Si, nP, 1);

Fe2_mica = repmat(mica.Fe2, nP, 1);
Fe3_mica = repmat(mica.Fe3, nP, 1);
FeUsed_mica = Fe2_mica;
Mg_mica = repmat(mica.Mg, nP, 1);
Mn_mica = repmat(mica.Mn, nP, 1);
Ca_mica = repmat(mica.Ca, nP, 1);
Ti_mica = repmat(mica.Ti, nP, 1);
Al_mica = repmat(mica.Al, nP, 1);
Si_mica = repmat(mica.Si, nP, 1);
K_mica = repmat(mica.K, nP, 1);
Na_mica = repmat(mica.Na, nP, 1);

% Calculate the Mg/Fe2+ exchange coefficient. NaN propagates naturally.
MgFe_grt = Mg_grt ./ FeUsed_grt;
MgFe_mica = Mg_mica ./ FeUsed_mica;
K_exchange = MgFe_grt ./ MgFe_mica;

% A finite zero in any exchange input makes log(K) mathematically invalid.
% Replace the exchange result with NaN before calculating temperature, so a
% denominator of +/-Inf cannot produce an invalid 0 K result.
zeroExchangeInput = ...
    (isfinite(FeUsed_grt) & FeUsed_grt == 0) | ...
    (isfinite(Mg_grt) & Mg_grt == 0) | ...
    (isfinite(FeUsed_mica) & FeUsed_mica == 0) | ...
    (isfinite(Mg_mica) & Mg_mica == 0);
K_exchange(zeroExchangeInput) = NaN;
lnK = log(K_exchange);

% Eq. (1): directly calibrated at 2.07 kbar. Its temperature is repeated for
% every pressure row because pressure does not appear in the equation.
denEq1 = 0.782 - lnK;
TEq1_K = 2109 ./ denEq1;
invalidEq1 = ~isfinite(denEq1) | denEq1 <= 0;
TEq1_K(invalidEq1) = NaN;
TEq1_C = TEq1_K - 273.15;

% Eq. (7): pressure-corrected polybaric form.
numEq7 = 12454 + 0.057 .* P_bar;
denEq7 = 4.662 - 3 .* R_cal .* lnK;
TEq7_K = numEq7 ./ denEq7;
invalidEq7 = ~isfinite(numEq7) | ~isfinite(denEq7) | denEq7 <= 0;
TEq7_K(invalidEq7) = NaN;
TEq7_C = TEq7_K - 273.15;

% Pack all variables at their final nP-row size.
row = table( ...
    P_kbar, P_bar, repmat(R_cal, nP, 1), ...
    Fe2_grt, Fe3_grt, FeUsed_grt, Mg_grt, Mn_grt, Ca_grt, Al_grt, Si_grt, ...
    Fe2_mica, Fe3_mica, FeUsed_mica, Mg_mica, Mn_mica, Ca_mica, ...
    Ti_mica, Al_mica, Si_mica, K_mica, Na_mica, ...
    MgFe_grt, MgFe_mica, K_exchange, lnK, zeroExchangeInput, ...
    denEq1, TEq1_K, TEq1_C, numEq7, denEq7, TEq7_K, TEq7_C, ...
    'VariableNames', { ...
    'P_kbar', 'P_bar', 'R_cal', ...
    'Fe2_grt', 'Fe3_grt', 'FeUsed_grt', 'Mg_grt', 'Mn_grt', 'Ca_grt', ...
    'Al_grt', 'Si_grt', ...
    'Fe2_mica', 'Fe3_mica', 'FeUsed_mica', 'Mg_mica', 'Mn_mica', ...
    'Ca_mica', 'Ti_mica', 'Al_mica', 'Si_mica', 'K_mica', 'Na_mica', ...
    'MgFe_grt', 'MgFe_mica', 'K_exchange', 'lnK', ...
    'zeroExchangeInput', 'denEq1', 'TEq1_K', 'TEq1_C', ...
    'numEq7', 'denEq7', 'TEq7_K', 'TEq7_C'});

end

function mineral = prepareMineralRow(data_table, mineralLabel)
% prepareMineralRow
% Read scalar cation values from a one-row table. Missing optional variables
% become zero; existing NaN values remain NaN.

if height(data_table) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();
mineral.Fe2 = getVarOrError( ...
    data_table, 'Fe_cation_apfu', mineralLabel);
mineral.Mg = getVarOrError( ...
    data_table, 'Mg_cation_apfu', mineralLabel);

mineral.Fe3 = getVarOrZero(data_table, 'Fe3_cation_apfu', mineralLabel);
mineral.Mn = getVarOrZero(data_table, 'Mn_cation_apfu', mineralLabel);
mineral.Ca = getVarOrZero(data_table, 'Ca_cation_apfu', mineralLabel);
mineral.Ti = getVarOrZero(data_table, 'Ti_cation_apfu', mineralLabel);
mineral.Al = getVarOrZero(data_table, 'Al_cation_apfu', mineralLabel);
mineral.Si = getVarOrZero(data_table, 'Si_cation_apfu', mineralLabel);
mineral.K = getVarOrZero(data_table, 'K_cation_apfu', mineralLabel);
mineral.Na = getVarOrZero(data_table, 'Na_cation_apfu', mineralLabel);

end

function value = getVarOrError(data_table, variableName, mineralLabel)
% getVarOrError
% Read a required scalar. NaN and zero are allowed; negative finite values
% and Inf are rejected.

if ~ismember(variableName, data_table.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = data_table.(variableName);
if ~isnumeric(value) || ~isscalar(value) || isinf(value) || ...
        (isfinite(value) && value < 0)
    error(['%s.%s must be a finite non-negative numeric scalar or NaN.'], ...
        mineralLabel, variableName);
end

end

function value = getVarOrZero(data_table, variableName, mineralLabel)
% getVarOrZero
% Read an optional scalar. Only a missing variable becomes zero; an existing
% NaN remains NaN. Negative finite values and Inf are rejected.

if ismember(variableName, data_table.Properties.VariableNames)
    value = data_table.(variableName);
    if ~isnumeric(value) || ~isscalar(value) || isinf(value) || ...
            (isfinite(value) && value < 0)
        error(['%s.%s must be a finite non-negative numeric scalar or NaN.'], ...
            mineralLabel, variableName);
    end
else
    value = 0;
end

end

function printTemperatureCalibrationWarning( ...
        temperature_degC, equationName, calibrationMin_degC, ...
        calibrationMax_degC, selectedCode_grt, selectedCode_mica)
% printTemperatureCalibrationWarning
% Print an fprintf warning for finite temperatures outside 550-800 degreeC.

finiteTemperature = isfinite(temperature_degC);
outsideCalibration = finiteTemperature & ...
    (temperature_degC < calibrationMin_degC | ...
     temperature_degC > calibrationMax_degC);

if any(outsideCalibration)
    finiteValues = temperature_degC(finiteTemperature);
    fprintf(2, ...
        ['WARNING: Calculated %s temperature is outside the direct experimental ' ...
         'calibration range of Ferry and Spear (1978): 550-800 degreeC. ' ...
         '%d of %d finite temperature point(s) are outside the range; ' ...
         'calculated finite range = %.4g-%.4g degreeC for %s & %s.\n'], ...
        equationName, ...
        sum(outsideCalibration), ...
        sum(finiteTemperature), ...
        min(finiteValues), max(finiteValues), ...
        char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

end
