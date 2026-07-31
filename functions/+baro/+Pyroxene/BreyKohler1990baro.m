function results = BreyKohler1990baro(rawdata_struct, T_degreeC)
% functions/+baro/+Pyroxene/BreyKohler1990baro.m
% Tested with MATLAB R2024b
%
% Al-in-Orthopyroxene geobarometer for Garnet peridotites
% Brey, G.P. and Kohler, T. (1990)
% Journal of Petrology, 31, 1353-1378
% DOI: https://doi.org/10.1093/petrology/31.6.1353
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis and one
% Garnet analysis and calculates pressure using the PBKN Al-in-Opx
% geobarometer of Brey and Kohler (1990).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Opx-Grt pair, one output row is
% returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% The PBKN formulation was developed for equilibrated garnet lherzolite /
% garnet peridotite assemblages containing coexisting orthopyroxene and
% garnet. The direct natural-system validation shown for PBKN spans broadly:
%
%   Temperature : approximately 900-1400 degreeC
%   Pressure    : approximately 28-60 kbar
%   Assemblage  : four-phase garnet lherzolite / garnet peridotite
%
% The 28-60 kbar test range is discussed on p. 1368, the PBKN formulation
% and site occupancies are given in Table 3 on p. 1372, and reproduction of
% the natural-system experiments is shown in Fig. 10 on p. 1373. PBKN alone
% reproduces those experiments to approximately +/-2.2 kbar (1 sigma;
% p. 1373). The preferred simultaneous TBKN + PBKN combination reproduces
% the experimental conditions to approximately +/-20 degreeC and +/-3 kbar
% (1 sigma; pp. 1374-1375).
%
% Orthopyroxene and garnet must represent the same equilibrium assemblage.
% Zoning, reaction rims, exsolution, retrograde re-equilibration, or pairing
% analyses from different generations can produce misleading pressures.
% Temperature errors propagate into pressure because PBKN is explicitly
% temperature dependent. Brey and Kohler (1990) further caution that
% correlated analytical or systematic errors can generate apparent P-T
% trends along conductive geotherms (pp. 1373-1376).
%
% Site occupancies must follow Table 3 (p. 1372), including the Opx
% corrections involving Na, Cr, Fe3+, Ti, and Mn. The input Opx cations
% should be normalized on the project Opx apfu basis and the Garnet cations
% on the project Garnet apfu basis consistently with the source tables.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 900-1400 degreeC,
%   2) finite calculated pressure is outside 28-60 kbar,
%   3) a required calculation input contains NaN,
%   4) a derived site occupancy, KD, or discriminant is outside its
%      mathematical domain, or
%   5) a calculated pressure is NaN, Inf, or negative.
%
% Calculations outside the ranges above are extrapolations. Passing the
% range checks does not by itself demonstrate mineral equilibrium.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx    : table
%   rawdata_struct.Garnet : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The following cation variables are
% required because they are used directly in Table 3:
%
%   Opx table variables:
%     Al_cation_apfu
%     Mg_cation_apfu
%     Fe_cation_apfu
%     Ca_cation_apfu
%     Na_cation_apfu
%     Mn_cation_apfu
%     Cr_cation_apfu
%     Ti_cation_apfu
%     Fe3_cation_apfu
%
%   Garnet table variables:
%     Al_cation_apfu
%     Mg_cation_apfu
%     Fe_cation_apfu
%     Ca_cation_apfu
%     Mn_cation_apfu
%     Cr_cation_apfu
%
%   Optional Garnet variable retained in the output when present:
%     Fe3_cation_apfu
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if it makes a ratio or logarithm
% undefined, the resulting NaN/Inf is retained and reported.
%
% No liquid composition is used by this barometer. Therefore, exclusion of
% Liq F and Cl from cationTotal_liq and from NaN-input warnings is not
% applicable to this function.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (worked-example-consistent implementation)
%
% Site occupancies follow Table 3 (p. 1372):
%
%   jadeite_opx = Na - Cr - Fe3+ - 2Ti
%
%   X_Al_M1_opx = (Al + Na - Cr - Fe3+ - 2Ti) / 2
%
%   if jadeite_opx < 0:
%       X_Al_TS_M1 = (Al + jadeite_opx) / 2
%   else:
%       X_Al_TS_M1 = (Al - jadeite_opx) / 2
%
%   X_MF_opx = Mg / (Mg + Fe)
%   X_MF_M1_capacity = 1 - X_Al_M1_opx - Cr - Fe3+ - Ti
%   X_MF_M2_capacity = 1 - Ca - Na - Mn
%
%   X_Mg_M1 = X_MF_M1_capacity * X_MF_opx
%   X_Fe_M1 = X_MF_M1_capacity * (1 - X_MF_opx)
%   X_Mg_M2 = X_MF_M2_capacity * X_MF_opx
%   X_Fe_M2 = X_MF_M2_capacity * (1 - X_MF_opx)
%
%   X_Al_grt = Al / (Al + Cr)
%   X_Cr_grt = Cr / (Cr + Al)
%   X_Ca_grt = Ca / (Ca + Mg + Fe + Mn)
%   X_Fe_grt = Fe / (Ca + Mg + Fe + Mn)
%   X_Mg_grt = Mg / (Ca + Mg + Fe + Mn)
%
%   KD = [(1 - X_Ca_grt)^3 * (X_Al_grt)^2] /
%        [X_MF_M1_capacity * (X_MF_M2_capacity)^2 * X_Al_TS_M1]
%
%   C1 = -R*T*ln(KD) - 5510 + 88.91*T - 19*sqrt(T)
%        + 3*(X_Ca_grt)^2*82458
%        + X_Mg_M1*X_Fe_M1*(80942 - 46.7*T)
%        - 3*X_Fe_grt*X_Ca_grt*17793
%        - X_Ca_grt*X_Cr_grt*(1.164e6 - 420.4*T)
%        - X_Fe_grt*X_Cr_grt*(-1.25e6 + 565*T)
%
%   C2 = -0.832 - 8.78e-5*(T - 298)
%        + 3*(X_Ca_grt)^2*3.305
%        - X_Ca_grt*X_Cr_grt*13.45
%        + X_Fe_grt*X_Cr_grt*10.5
%
%   C3 = 16.6e-4
%
%   P_BKN(kbar) =
%     [-C2 + sqrt(C2^2 + 4*C3*C1/1000)] / (20*C3)
%
% IMPORTANT IMPLEMENTATION NOTE:
% The literal quadratic-root typography printed in Table 3 does not
% reproduce the worked example in the Appendix (p. 1378). The positive root
% and factor-of-10 unit conversion used above reproduce the published worked
% example (approximately 59.1 kbar for the tabulated example). The printed
% Table 3 site occupancies and KD definition are otherwise followed here.
%
% R = 8.3143 J K^-1 mol^-1; T is in Kelvin; P is returned in kbar.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = BreyKohler1990baro(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing Opx and Garnet tables
%   T_degreeC      : non-negative numeric scalar or vector in degreeC.
%                    NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Opx-Grt pair.
%

%% Input validation
% Accept scalar or vector temperature input so that the fixed-temperature
% and temperature-range launchers use the same implementation.
if nargin < 2
    error('BreyKohler1990baro requires (rawdata_struct, T_degreeC).');
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

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end
if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end

dataset_opx = rawdata_struct.Opx;
dataset_grt = rawdata_struct.Garnet;

requiredOpxVariables = {'Al_cation_apfu', 'Mg_cation_apfu', ...
    'Fe_cation_apfu', 'Ca_cation_apfu', 'Na_cation_apfu', ...
    'Mn_cation_apfu', 'Cr_cation_apfu', 'Ti_cation_apfu', ...
    'Fe3_cation_apfu'};
requiredGarnetVariables = {'Al_cation_apfu', 'Mg_cation_apfu', ...
    'Fe_cation_apfu', 'Ca_cation_apfu', 'Mn_cation_apfu', ...
    'Cr_cation_apfu'};

checkRequiredVariables(dataset_opx, requiredOpxVariables, 'Opx');
checkRequiredVariables(dataset_grt, requiredGarnetVariables, 'Garnet');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the results table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct natural-system validation envelope used for non-stopping warnings.
calibrationT_min_degreeC = 900;
calibrationT_max_degreeC = 1400;
calibrationP_min_kbar = 28;
calibrationP_max_kbar = 60;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);
temperatureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a completed calculation.
disp('=== Step 3: Selecting a data code from the list (Opx) ===');

while true
    % ----- Orthopyroxene selection -----
    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Orthopyroxene (Opx) data:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Opx selected: ' char(string(selectedCode_opx))]);

    % ----- Garnet selection -----
    disp('=== Step 4: Selecting a data code from the list (Garnet) ===');

    dataCodes_grt = dataset_grt{:, 1};

    [selectedIdx_grt, ok] = listdlg( ...
        'PromptString', 'Please select the Garnet data:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_grt)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_grt)
        disp('Selection canceled');
        break;
    end

    selectedCode_grt = dataCodes_grt(selectedIdx_grt);
    disp(['Garnet selected: ' char(string(selectedCode_grt))]);

    % ----- Calculation -----
    % Opx and Garnet are selected independently; row numbers are not assumed
    % to correspond between the two source tables.
    disp('=== Step 5: Calculating the pressure ===');

    selectedData_opx = dataset_opx(selectedIdx_opx, :);
    selectedData_grt = dataset_grt(selectedIdx_grt, :);

    % List NaN only in variables directly used by PBKN. NaN values are not
    % replaced and do not stop the calculation.
    nanInputNames = findNaNInputs( ...
        selectedData_opx, selectedData_grt, T_degreeC);

    % Reject only Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(selectedData_opx, selectedData_grt);

    row = calcPressure(selectedData_opx, selectedData_grt, T_degreeC);

    % Repeat identifiers for all temperatures in the current calculation.
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row = movevars(row, {'dataCode_opx', 'dataCode_grt'}, 'Before', 1);

    % Store one block per selected mineral pair. Expand the cell buffer only
    % when its current capacity has been exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_grt)) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_grt)) ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Input temperature is common to all selected pairs, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperatureValues = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the approximate direct ' ...
             'natural-system validation range of Brey and Kohler (1990): ' ...
             '900-1400 degreeC (Figs. 6 and 10; pp. 1364 and 1373). ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperatureValues), ...
            min(finiteTemperatureValues), ...
            max(finiteTemperatureValues));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures fall outside the direct PBKN
    % natural-system validation envelope.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the approximate ' ...
             '28-60 kbar PBKN natural-system validation range of Brey and ' ...
             'Kohler (1990; pp. 1368 and 1373). %d of %d finite pressure ' ...
             'point(s) are outside the range; calculated finite range = ' ...
             '%.4g-%.4g kbar for %s & %s.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)));
    end

    % List the exact required input names containing NaN. For vector
    % temperature input, the NaN element indices are included.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the PBKN input(s) for %s & %s: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'pressure may remain NaN.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Report derived site quantities that fall outside the mathematical
    % domain without replacing the source cation values.
    invalidSiteNames = findInvalidSiteQuantities(row);
    if ~isempty(invalidSiteNames)
        fprintf(2, ...
            ['WARNING: Derived PBKN site quantity or equilibrium term is outside ' ...
             'its valid mathematical domain for %s & %s: %s.\n' ...
             '         The affected pressure values were retained as NaN/Inf ' ...
             'where applicable.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)), ...
            char(strjoin(invalidSiteNames, ', ')));
    end

    % Explicitly report negative finite discriminants before the result check.
    negativeDiscriminant = isfinite(row.discriminant) & row.discriminant < 0;
    if any(negativeDiscriminant)
        fprintf(2, ...
            ['WARNING: Negative quadratic discriminant was calculated for %s & %s ' ...
             '(%d of %d points). The corresponding pressure values were retained ' ...
             'as NaN.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)), ...
            sum(negativeDiscriminant), ...
            numel(row.discriminant));
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
            char(string(selectedCode_grt)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressure is retained for diagnosis but is physically
    % outside the application domain.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & %s ' ...
             '(%d of %d points). The values were retained for diagnostic purposes.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'BreyKohler1990baro', ...
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
function nanInputNames = findNaNInputs(data_opx, data_grt, T_degreeC)
% findNaNInputs
% Return names of PBKN calculation inputs containing NaN. NaN values do not
% cause an error and are never replaced by zero.

opxVariables = {'Al_cation_apfu', 'Mg_cation_apfu', ...
    'Fe_cation_apfu', 'Ca_cation_apfu', 'Na_cation_apfu', ...
    'Mn_cation_apfu', 'Cr_cation_apfu', 'Ti_cation_apfu', ...
    'Fe3_cation_apfu'};
garnetVariables = {'Al_cation_apfu', 'Mg_cation_apfu', ...
    'Fe_cation_apfu', 'Ca_cation_apfu', 'Mn_cation_apfu', ...
    'Cr_cation_apfu'};

maxNames = 1 + numel(opxVariables) + numel(garnetVariables);
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

for i = 1:numel(opxVariables)
    variableName = opxVariables{i};
    variableValue = data_opx.(variableName);
    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Opx." + string(variableName);
    end
end

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    variableValue = data_grt.(variableName);
    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Garnet." + string(variableName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_opx, data_grt)
% validateNonNegativeInputs
% Reject Inf and finite negative values in variables directly used by PBKN.
% Zero and NaN are intentionally allowed and retained.

opxVariables = {'Al_cation_apfu', 'Mg_cation_apfu', ...
    'Fe_cation_apfu', 'Ca_cation_apfu', 'Na_cation_apfu', ...
    'Mn_cation_apfu', 'Cr_cation_apfu', 'Ti_cation_apfu', ...
    'Fe3_cation_apfu'};
garnetVariables = {'Al_cation_apfu', 'Mg_cation_apfu', ...
    'Fe_cation_apfu', 'Ca_cation_apfu', 'Mn_cation_apfu', ...
    'Cr_cation_apfu'};

maxNames = numel(opxVariables) + numel(garnetVariables);
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

for i = 1:numel(opxVariables)
    variableName = opxVariables{i};
    variableValue = data_opx.(variableName);
    if any(isinf(variableValue(:)) | ...
            (isfinite(variableValue(:)) & variableValue(:) < 0))
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Opx." + string(variableName);
    end
end

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    variableValue = data_grt.(variableName);
    if any(isinf(variableValue(:)) | ...
            (isfinite(variableValue(:)) & variableValue(:) < 0))
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Garnet." + string(variableName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['BreyKohler1990baro: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_opx, data_grt, T_degreeC)
% calcPressure
% Compute PBKN pressure for one Opx row and one Garnet row at one or more
% input temperatures. NaN values propagate through the calculation.
%
% Inputs:
%   data_opx   : 1-row Opx table
%   data_grt   : 1-row Garnet table
%   T_degreeC  : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Physical constant and quadratic coefficient.
R = 8.3143;
C3_scalar = 16.6e-4;

% Extract one-row cation data. Required columns are checked explicitly;
% optional Garnet Fe3 is retained as NaN when absent.
opx = prepareOpxRow(data_opx);
grt = prepareGarnetRow(data_grt);

% ----- Opx site occupancies after Table 3 (p. 1372) -----
jadeite_opx_scalar = opx.Na - opx.Cr - opx.Fe3 - 2 .* opx.Ti;
X_Al_M1_opx_scalar = ...
    (opx.Al + opx.Na - opx.Cr - opx.Fe3 - 2 .* opx.Ti) ./ 2;

if isnan(opx.Al) || isnan(jadeite_opx_scalar)
    X_Al_TS_M1_scalar = NaN;
elseif jadeite_opx_scalar < 0
    X_Al_TS_M1_scalar = (opx.Al + jadeite_opx_scalar) ./ 2;
else
    X_Al_TS_M1_scalar = (opx.Al - jadeite_opx_scalar) ./ 2;
end

X_MF_opx_scalar = opx.Mg ./ (opx.Mg + opx.Fe);
X_MF_M1_capacity_scalar = ...
    1 - X_Al_M1_opx_scalar - opx.Cr - opx.Fe3 - opx.Ti;
X_MF_M2_capacity_scalar = 1 - opx.Ca - opx.Na - opx.Mn;

X_Mg_M1_scalar = X_MF_M1_capacity_scalar .* X_MF_opx_scalar;
X_Fe_M1_scalar = X_MF_M1_capacity_scalar .* (1 - X_MF_opx_scalar);
X_Mg_M2_scalar = X_MF_M2_capacity_scalar .* X_MF_opx_scalar;
X_Fe_M2_scalar = X_MF_M2_capacity_scalar .* (1 - X_MF_opx_scalar);

% ----- Garnet site fractions after Table 3 (p. 1372) -----
grt_AlCr_sum_scalar = grt.Al + grt.Cr;
grt_X_site_sum_scalar = grt.Ca + grt.Mg + grt.Fe + grt.Mn;

X_Al_grt_scalar = grt.Al ./ grt_AlCr_sum_scalar;
X_Cr_grt_scalar = grt.Cr ./ grt_AlCr_sum_scalar;
X_Ca_grt_scalar = grt.Ca ./ grt_X_site_sum_scalar;
X_Fe_grt_scalar = grt.Fe ./ grt_X_site_sum_scalar;
X_Mg_grt_scalar = grt.Mg ./ grt_X_site_sum_scalar;

% ----- Equilibrium constant -----
KD_scalar = ((1 - X_Ca_grt_scalar).^3 .* (X_Al_grt_scalar).^2) ./ ...
    (X_MF_M1_capacity_scalar .* (X_MF_M2_capacity_scalar).^2 .* ...
     X_Al_TS_M1_scalar);

% Avoid a complex logarithm for non-positive KD. NaN is retained.
lnKD_scalar = NaN;
if isfinite(KD_scalar) && KD_scalar > 0
    lnKD_scalar = log(KD_scalar);
elseif isinf(KD_scalar) && KD_scalar > 0
    lnKD_scalar = Inf;
end

% Expand composition-dependent scalars to the temperature-vector length.
R_output = repmat(R, nT, 1);
C3 = repmat(C3_scalar, nT, 1);

Al_opx = repmat(opx.Al, nT, 1);
Mg_opx = repmat(opx.Mg, nT, 1);
Fe_opx = repmat(opx.Fe, nT, 1);
Ca_opx = repmat(opx.Ca, nT, 1);
Na_opx = repmat(opx.Na, nT, 1);
Mn_opx = repmat(opx.Mn, nT, 1);
Cr_opx = repmat(opx.Cr, nT, 1);
Ti_opx = repmat(opx.Ti, nT, 1);
Fe3_opx = repmat(opx.Fe3, nT, 1);

Al_grt = repmat(grt.Al, nT, 1);
Mg_grt = repmat(grt.Mg, nT, 1);
Fe_grt = repmat(grt.Fe, nT, 1);
Ca_grt = repmat(grt.Ca, nT, 1);
Mn_grt = repmat(grt.Mn, nT, 1);
Cr_grt = repmat(grt.Cr, nT, 1);
Fe3_grt = repmat(grt.Fe3, nT, 1);

jadeite_opx = repmat(jadeite_opx_scalar, nT, 1);
X_Al_M1_opx = repmat(X_Al_M1_opx_scalar, nT, 1);
X_Al_TS_M1 = repmat(X_Al_TS_M1_scalar, nT, 1);
X_MF_opx = repmat(X_MF_opx_scalar, nT, 1);
X_MF_M1_capacity = repmat(X_MF_M1_capacity_scalar, nT, 1);
X_MF_M2_capacity = repmat(X_MF_M2_capacity_scalar, nT, 1);
X_Mg_M1 = repmat(X_Mg_M1_scalar, nT, 1);
X_Fe_M1 = repmat(X_Fe_M1_scalar, nT, 1);
X_Mg_M2 = repmat(X_Mg_M2_scalar, nT, 1);
X_Fe_M2 = repmat(X_Fe_M2_scalar, nT, 1);

grt_AlCr_sum = repmat(grt_AlCr_sum_scalar, nT, 1);
grt_X_site_sum = repmat(grt_X_site_sum_scalar, nT, 1);
X_Al_grt = repmat(X_Al_grt_scalar, nT, 1);
X_Cr_grt = repmat(X_Cr_grt_scalar, nT, 1);
X_Ca_grt = repmat(X_Ca_grt_scalar, nT, 1);
X_Fe_grt = repmat(X_Fe_grt_scalar, nT, 1);
X_Mg_grt = repmat(X_Mg_grt_scalar, nT, 1);
KD = repmat(KD_scalar, nT, 1);
lnKD = repmat(lnKD_scalar, nT, 1);

% ----- PBKN coefficients -----
C1 = -R .* T_K .* lnKD ...
    - 5510 ...
    + 88.91 .* T_K ...
    - 19 .* sqrt(T_K) ...
    + 3 .* (X_Ca_grt).^2 .* 82458 ...
    + X_Mg_M1 .* X_Fe_M1 .* (80942 - 46.7 .* T_K) ...
    - 3 .* X_Fe_grt .* X_Ca_grt .* 17793 ...
    - X_Ca_grt .* X_Cr_grt .* (1.164e6 - 420.4 .* T_K) ...
    - X_Fe_grt .* X_Cr_grt .* (-1.25e6 + 565 .* T_K);

C2 = -0.832 ...
    - 8.78e-5 .* (T_K - 298) ...
    + 3 .* (X_Ca_grt).^2 .* 3.305 ...
    - X_Ca_grt .* X_Cr_grt .* 13.45 ...
    + X_Fe_grt .* X_Cr_grt .* 10.5;

% The positive quadratic root and factor-of-10 pressure-unit conversion are
% required to reproduce the Appendix worked example on p. 1378.
discriminant = C2.^2 + 4 .* C3 .* C1 ./ 1000;
sqrtDiscriminant = NaN(nT, 1);
validDiscriminant = ~isnan(discriminant) & discriminant >= 0;
sqrtDiscriminant(validDiscriminant) = sqrt(discriminant(validDiscriminant));

P_kbar = (-C2 + sqrtDiscriminant) ./ (20 .* C3);

% Applicability and diagnostic flags. These flags do not prove equilibrium.
isWithinCalibrationTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 900 & T_degreeC <= 1400;
isWithinCalibrationPRange = ...
    isfinite(P_kbar) & P_kbar >= 28 & P_kbar <= 60;
isWithinSiteDomain = ...
    isfinite(X_Al_TS_M1) & X_Al_TS_M1 > 0 & ...
    isfinite(X_MF_M1_capacity) & X_MF_M1_capacity > 0 & ...
    isfinite(X_MF_M2_capacity) & X_MF_M2_capacity > 0 & ...
    isfinite(X_MF_opx) & X_MF_opx >= 0 & X_MF_opx <= 1 & ...
    isfinite(grt_AlCr_sum) & grt_AlCr_sum > 0 & ...
    isfinite(grt_X_site_sum) & grt_X_site_sum > 0 & ...
    isfinite(KD) & KD > 0 & ...
    isfinite(discriminant) & discriminant >= 0;

% Pack outputs using equal-length, pre-sized vectors.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.R = R_output;
row.P_kbar = P_kbar;
row.P_uncertainty_1sigma_kbar = repmat(2.2, nT, 1);

row.KD = KD;
row.lnKD = lnKD;
row.C1 = C1;
row.C2 = C2;
row.C3 = C3;
row.discriminant = discriminant;
row.sqrtDiscriminant = sqrtDiscriminant;

row.jadeite_opx = jadeite_opx;
row.X_Al_M1_opx = X_Al_M1_opx;
row.X_Al_TS_M1 = X_Al_TS_M1;
row.X_MF_opx = X_MF_opx;
row.X_MF_M1_capacity = X_MF_M1_capacity;
row.X_MF_M2_capacity = X_MF_M2_capacity;
row.X_Mg_M1 = X_Mg_M1;
row.X_Fe_M1 = X_Fe_M1;
row.X_Mg_M2 = X_Mg_M2;
row.X_Fe_M2 = X_Fe_M2;

row.grt_AlCr_sum = grt_AlCr_sum;
row.grt_X_site_sum = grt_X_site_sum;
row.X_Al_grt = X_Al_grt;
row.X_Cr_grt = X_Cr_grt;
row.X_Ca_grt = X_Ca_grt;
row.X_Fe_grt = X_Fe_grt;
row.X_Mg_grt = X_Mg_grt;

row.Al_opx = Al_opx;
row.Mg_opx = Mg_opx;
row.Fe_opx = Fe_opx;
row.Ca_opx = Ca_opx;
row.Na_opx = Na_opx;
row.Mn_opx = Mn_opx;
row.Cr_opx = Cr_opx;
row.Ti_opx = Ti_opx;
row.Fe3_opx = Fe3_opx;

row.Al_grt = Al_grt;
row.Mg_grt = Mg_grt;
row.Fe_grt = Fe_grt;
row.Ca_grt = Ca_grt;
row.Mn_grt = Mn_grt;
row.Cr_grt = Cr_grt;
row.Fe3_grt = Fe3_grt;

row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinCalibrationPRange = isWithinCalibrationPRange;
row.isWithinSiteDomain = isWithinSiteDomain;

end

function opx = prepareOpxRow(data_opx)
% prepareOpxRow
% Extract one-row Opx cation data without replacing NaN by zero.

if height(data_opx) ~= 1
    error('Opx input must be a 1-row table.');
end

opx = struct();
opx.Al = getVarOrError(data_opx, 'Al_cation_apfu', 'Opx');
opx.Mg = getVarOrError(data_opx, 'Mg_cation_apfu', 'Opx');
opx.Fe = getVarOrError(data_opx, 'Fe_cation_apfu', 'Opx');
opx.Ca = getVarOrError(data_opx, 'Ca_cation_apfu', 'Opx');
opx.Na = getVarOrError(data_opx, 'Na_cation_apfu', 'Opx');
opx.Mn = getVarOrError(data_opx, 'Mn_cation_apfu', 'Opx');
opx.Cr = getVarOrError(data_opx, 'Cr_cation_apfu', 'Opx');
opx.Ti = getVarOrError(data_opx, 'Ti_cation_apfu', 'Opx');
opx.Fe3 = getVarOrError(data_opx, 'Fe3_cation_apfu', 'Opx');

end

function grt = prepareGarnetRow(data_grt)
% prepareGarnetRow
% Extract one-row Garnet cation data without replacing NaN by zero.

if height(data_grt) ~= 1
    error('Garnet input must be a 1-row table.');
end

grt = struct();
grt.Al = getVarOrError(data_grt, 'Al_cation_apfu', 'Garnet');
grt.Mg = getVarOrError(data_grt, 'Mg_cation_apfu', 'Garnet');
grt.Fe = getVarOrError(data_grt, 'Fe_cation_apfu', 'Garnet');
grt.Ca = getVarOrError(data_grt, 'Ca_cation_apfu', 'Garnet');
grt.Mn = getVarOrError(data_grt, 'Mn_cation_apfu', 'Garnet');
grt.Cr = getVarOrError(data_grt, 'Cr_cation_apfu', 'Garnet');
grt.Fe3 = getVarOrNaN(data_grt, 'Fe3_cation_apfu');

end

function invalidSiteNames = findInvalidSiteQuantities(row)
% findInvalidSiteQuantities
% Return names of derived scalar site quantities outside their mathematical
% domain. Composition-dependent quantities are identical in every row, so
% the first row is sufficient for these checks.

maxNames = 8;
invalidSiteBuffer = strings(maxNames, 1);
nInvalidSite = 0;

if isempty(row)
    invalidSiteNames = strings(0, 1);
    return;
end

if ~isfinite(row.X_Al_M1_opx(1)) || row.X_Al_M1_opx(1) < 0
    nInvalidSite = nInvalidSite + 1;
    invalidSiteBuffer(nInvalidSite) = "X_Al_M1_opx";
end
if ~isfinite(row.X_Al_TS_M1(1)) || row.X_Al_TS_M1(1) <= 0
    nInvalidSite = nInvalidSite + 1;
    invalidSiteBuffer(nInvalidSite) = "X_Al_TS_M1";
end
if ~isfinite(row.X_MF_M1_capacity(1)) || row.X_MF_M1_capacity(1) <= 0
    nInvalidSite = nInvalidSite + 1;
    invalidSiteBuffer(nInvalidSite) = "X_MF_M1_capacity";
end
if ~isfinite(row.X_MF_M2_capacity(1)) || row.X_MF_M2_capacity(1) <= 0
    nInvalidSite = nInvalidSite + 1;
    invalidSiteBuffer(nInvalidSite) = "X_MF_M2_capacity";
end
if ~isfinite(row.X_MF_opx(1)) || row.X_MF_opx(1) < 0 || row.X_MF_opx(1) > 1
    nInvalidSite = nInvalidSite + 1;
    invalidSiteBuffer(nInvalidSite) = "X_MF_opx";
end
if ~isfinite(row.grt_AlCr_sum(1)) || row.grt_AlCr_sum(1) <= 0
    nInvalidSite = nInvalidSite + 1;
    invalidSiteBuffer(nInvalidSite) = "grt_AlCr_sum";
end
if ~isfinite(row.grt_X_site_sum(1)) || row.grt_X_site_sum(1) <= 0
    nInvalidSite = nInvalidSite + 1;
    invalidSiteBuffer(nInvalidSite) = "grt_X_site_sum";
end
if ~isfinite(row.KD(1)) || row.KD(1) <= 0
    nInvalidSite = nInvalidSite + 1;
    invalidSiteBuffer(nInvalidSite) = "KD";
end

invalidSiteNames = invalidSiteBuffer(1:nInvalidSite);

end

function checkRequiredVariables(tbl, requiredVariables, tableName)
% checkRequiredVariables
% Stop before the interactive loop when required calculation columns are
% absent. Existing NaN values within present columns are allowed.

missingVariables = ...
    requiredVariables(~ismember(requiredVariables, tbl.Properties.VariableNames));

if ~isempty(missingVariables)
    error('%s table must contain variable(s): %s', ...
        tableName, strjoin(missingVariables, ', '));
end

end

function value = getVarOrError(tbl, variableName, mineralLabel)
% getVarOrError
% Retrieve a required scalar table variable without altering NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);
if ~isscalar(value)
    error('Variable %s must be scalar in a 1-row table.', variableName);
end

end

function value = getVarOrNaN(tbl, variableName)
% getVarOrNaN
% Retrieve an optional scalar variable. Missing optional variables are
% represented by NaN, never by zero.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);
    if ~isscalar(value)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end
else
    value = NaN;
end

end
