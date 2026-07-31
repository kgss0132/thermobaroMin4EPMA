function results = Wells1979(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/Wells1979.m
% Tested with MATLAB R2024b
%
% Garnet-clinopyroxene Fe-Mg exchange thermometer
% Wells, P.R.A. (1979)
% Journal of Petrology, 20, 187-226
% DOI: https://doi.org/10.1093/petrology/20.2.187
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis and calculates temperature using equation (3) of
% Wells (1979).
%
% Pressure may be supplied as a scalar or vector. Therefore, this function
% can be called from both startThermoCalc_fixedP and startThermoCalc_rangeP.
% One output row is returned for every pressure value for each selected pair.
%
% -------------------------------------------------------------------------
% APPLICATION RANGE AND IMPORTANT CAUTIONS
%
% Wells (1979) did not report a single formal experimental calibration box
% for equation (3). The equation combines thermochemical data obtained by
% linear regression of experimental results from several earlier studies.
% It assumes that departures from ideal mixing in garnet and clinopyroxene
% mutually cancel (p. 197). Consequently, the following values are a
% demonstrated NATURAL-APPLICATION range, not formal calibration limits:
%
%   Temperature calculated by equation (3) : 750-830 degreeC at 10 kbar
%                                             (pp. 197-198)
%   Pressure of the investigated granulites : 8.5-11.5 kbar; most likely
%                                             9.5-11.0 kbar (p. 199)
%   Uncertainty printed with equation (3)    : approximately +/-40 K
%                                             (p. 197)
%
% Wells (1979) considered mineral pairs suitable for thermobarometry only
% when the phases were in chemical equilibrium. The listed criteria are:
%   1) absence of mineral zoning,
%   2) granoblastic polygonal textures,
%   3) regular Fe, Mg, Ni, Cr, and Mn distribution between coexisting
%      phases, and
%   4) internally consistent pressure and temperature estimates (p. 196).
%
% Equation (3) describes ferrous Fe-Mg exchange. Fe_cation_apfu is therefore
% used as Fe2+. Fe3_cation_apfu is retained in the output only for
% traceability and is NOT added to Fe2+. If Fe_cation_apfu already represents
% total Fe rather than Fe2+, the user must correct or explicitly justify the
% iron treatment before applying this thermometer.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside 8.5-11.5 kbar,
%   2) a finite calculated temperature is outside 750-830 degreeC,
%   3) a calculation input contains NaN, or
%   4) the exchange coefficient or calculated temperature is non-finite.
% The first two checks refer only to the demonstrated natural application in
% Wells (1979), not to formal experimental calibration limits.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Cpx    : table
%
% The FIRST column of each table is treated as an identifier ("data code").
%
% Required Garnet variables used in the equation:
%   Fe_cation_apfu         % Fe2+
%   Mg_cation_apfu
%
% Required Clinopyroxene variables used in the equation:
%   Fe_cation_apfu         % Fe2+
%   Mg_cation_apfu
%
% Optional variables retained for traceability when present:
%   Fe3_cation_apfu, Mn_cation_apfu, Ca_cation_apfu,
%   Ti_cation_apfu, Al_cation_apfu, Si_cation_apfu,
%   Na_cation_apfu, K_cation_apfu
%
% Existing NaN values are never converted to zero. They propagate through
% the calculation when they occur in Fe2+ or Mg and are reported by fprintf.
% A missing optional column is reported in the output as NaN. Negative or
% infinite numeric inputs are prohibited; zero and NaN are allowed.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% The original Wells (1979) exchange term is the reciprocal of the commonly
% used Fe-Mg distribution coefficient:
%
%   K_Wells = (XMg_grt * XFe_cpx) / (XFe_grt * XMg_cpx)
%
% This implementation retains the conventional output definition:
%
%   KD = (Fe2+/Mg)_grt / (Fe2+/Mg)_cpx = 1 / K_Wells
%
% Therefore equation (3) is evaluated in the algebraically equivalent form:
%
%   T(K) = [24440 + 0.06524*(P_bar - 1)]
%          / [13.41 + 3*R*ln(KD)]
%
% where:
%   R     = 1.987 cal mol^-1 K^-1
%   P_bar = 1000 * P_kbar
%
% The plus sign before 3*R*ln(KD) is required because KD, as defined above,
% is the reciprocal of the logarithmic exchange term printed by Wells
% (1979). Using a minus sign with this conventional KD definition produces
% erroneous temperatures.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Wells1979(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Cpx tables
%   P_kbar         : finite, non-negative numeric scalar or vector
%
% Output:
%   results : table containing one row per pressure value for each pair
%

%% Input validation
if nargin < 2
    error('Wells1979 requires (rawdata_struct, P_kbar).');
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
% Table blocks are buffered and concatenated only once after the loop.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% These are demonstrated natural-application values, not formal
% experimental calibration limits.
applicationT_min_degC = 750;
applicationT_max_degC = 830;
applicationP_min_kbar = 8.5;
applicationP_max_kbar = 11.5;

pressureOutsideApplication = P_kbar < applicationP_min_kbar | ...
    P_kbar > applicationP_max_kbar;
pressureMessageIssued = false;

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
    disp('=== Step 5: Calculating the temperature ===');
    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    validateNonnegativeInputs(selectedData_grt, selectedData_cpx);
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_cpx);
    row = calcTemp(selectedData_grt, selectedData_cpx, P_kbar);

    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_cpx'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_C(1)) ' to ' num2str(row.T_C(end)) ' degreeC']);
    end

    if any(pressureOutsideApplication) && ~pressureMessageIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the Wells (1979) demonstrated ' ...
             'natural-application range of 8.5-11.5 kbar. %d of %d pressure ' ...
             'point(s) are outside; input range = %.4g-%.4g kbar. This is not ' ...
             'a formal calibration limit. Calculation was continued.\n'], ...
            sum(pressureOutsideApplication), numel(P_kbar), ...
            min(P_kbar), max(P_kbar));
        pressureMessageIssued = true;
    end

    finiteTemperature = isfinite(row.T_C);
    temperatureOutsideApplication = finiteTemperature & ...
        (row.T_C < applicationT_min_degC | row.T_C > applicationT_max_degC);
    if any(temperatureOutsideApplication)
        finiteValues = row.T_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the Wells (1979) ' ...
             'demonstrated equation-(3) range of 750-830 degreeC at about ' ...
             '10 kbar. %d of %d finite point(s) are outside; finite range = ' ...
             '%.4g-%.4g degreeC for %s & %s. This is not a formal calibration ' ...
             'limit. Calculation was continued.\n'], ...
            sum(temperatureOutsideApplication), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_grt)), char(string(selectedCode_cpx)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in calculation input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; calculation was not stopped.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidExchange = ~isfinite(row.KD) | row.KD <= 0 | ...
        ~isfinite(row.lnKD) | ~isfinite(row.denominator) | ...
        row.denominator == 0;
    if any(invalidExchange)
        fprintf(2, ...
            ['WARNING: Non-finite, non-positive, or singular exchange term(s) ' ...
             'occurred for %s & %s (%d of %d pressure points). Values were ' ...
             'retained unchanged.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_cpx)), ...
            sum(invalidExchange), numel(invalidExchange));
    end

    invalidTemperature = ~isfinite(row.T_C);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d). Values remain unchanged ' ...
             'in the output table.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), numel(row.T_C), ...
            sum(isnan(row.T_C)), sum(isinf(row.T_C)));
    end

    disp('--------------------------------------------------');
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Wells1979', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

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
% Confirm that every mandatory equation input exists.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
missingNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nMissing = 0;

for i = 1:numel(grtVariables)
    if ~ismember(grtVariables{i}, dataset_grt.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Garnet." + string(grtVariables{i});
    end
end
for i = 1:numel(cpxVariables)
    if ~ismember(cpxVariables{i}, dataset_cpx.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Cpx." + string(cpxVariables{i});
    end
end

if nMissing > 0
    error('Wells1979: missing required variable(s): %s.', ...
        char(strjoin(missingNames(1:nMissing), ', ')));
end

end

function [inputNames, inputValues] = collectCalculationInputs(data_grt, data_cpx)
% collectCalculationInputs
% Collect the four Fe2+-Mg exchange inputs used by equation (3).

inputNames = ["Garnet.Fe_cation_apfu"; "Garnet.Mg_cation_apfu"; ...
    "Cpx.Fe_cation_apfu"; "Cpx.Mg_cation_apfu"];
inputValues = NaN(4, 1);

inputValues(1) = getRequiredVar(data_grt, 'Fe_cation_apfu', 'Garnet');
inputValues(2) = getRequiredVar(data_grt, 'Mg_cation_apfu', 'Garnet');
inputValues(3) = getRequiredVar(data_cpx, 'Fe_cation_apfu', 'Cpx');
inputValues(4) = getRequiredVar(data_cpx, 'Mg_cation_apfu', 'Cpx');

end

function nanInputNames = findNaNInputs(data_grt, data_cpx)
% findNaNInputs
% List equation inputs containing NaN without replacing them.

[inputNames, inputValues] = collectCalculationInputs(data_grt, data_cpx);
nanInputNames = inputNames(isnan(inputValues));

end

function validateNonnegativeInputs(data_grt, data_cpx)
% validateNonnegativeInputs
% Reject negative or infinite equation inputs. Zero and NaN are allowed.

[inputNames, inputValues] = collectCalculationInputs(data_grt, data_cpx);
invalidInput = inputValues < 0 | isinf(inputValues);
if any(invalidInput)
    error(['Wells1979: calculation inputs must be non-negative and not ' ...
           'Inf. Invalid value(s): %s.'], ...
        char(strjoin(inputNames(invalidInput), ', ')));
end

end

function row = calcTemp(data_grt, data_cpx, P_kbar)
% calcTemp
% Calculate one Garnet-Cpx pair over a scalar or vector of pressures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);
R_cal = 1.987;
P_bar = P_kbar .* 1000;

Fe2_grt = repmat(getRequiredVar(data_grt, 'Fe_cation_apfu', 'Garnet'), nP, 1);
Mg_grt = repmat(getRequiredVar(data_grt, 'Mg_cation_apfu', 'Garnet'), nP, 1);
Fe2_cpx = repmat(getRequiredVar(data_cpx, 'Fe_cation_apfu', 'Cpx'), nP, 1);
Mg_cpx = repmat(getRequiredVar(data_cpx, 'Mg_cation_apfu', 'Cpx'), nP, 1);

FeMg_grt = Fe2_grt ./ Mg_grt;
FeMg_cpx = Fe2_cpx ./ Mg_cpx;
KD = FeMg_grt ./ FeMg_cpx;
lnKD = log(KD);
K_Wells = 1 ./ KD;
lnK_Wells = log(K_Wells);

numerator = 24440 + 0.06524 .* (P_bar - 1);
denominator = 13.41 + 3 .* R_cal .* lnKD;
T_K = numerator ./ denominator;
T_C = T_K - 273.15;

row = table();
row.P_kbar = P_kbar;
row.P_bar = P_bar;

row.Fe2_grt = Fe2_grt;
row.Fe3_grt = repmat(getOptionalVar(data_grt, ...
    'Fe3_cation_apfu', 'Garnet', NaN), nP, 1);
% Retain the legacy output name, but do not add Fe3 to the exchange Fe.
row.FeUsed_grt = Fe2_grt;
row.Mg_grt = Mg_grt;
row.Mn_grt = repmat(getOptionalVar(data_grt, ...
    'Mn_cation_apfu', 'Garnet', NaN), nP, 1);
row.Ca_grt = repmat(getOptionalVar(data_grt, ...
    'Ca_cation_apfu', 'Garnet', NaN), nP, 1);
row.Ti_grt = repmat(getOptionalVar(data_grt, ...
    'Ti_cation_apfu', 'Garnet', NaN), nP, 1);
row.Al_grt = repmat(getOptionalVar(data_grt, ...
    'Al_cation_apfu', 'Garnet', NaN), nP, 1);
row.Si_grt = repmat(getOptionalVar(data_grt, ...
    'Si_cation_apfu', 'Garnet', NaN), nP, 1);
row.Na_grt = repmat(getOptionalVar(data_grt, ...
    'Na_cation_apfu', 'Garnet', NaN), nP, 1);
row.K_grt = repmat(getOptionalVar(data_grt, ...
    'K_cation_apfu', 'Garnet', NaN), nP, 1);

row.Fe2_cpx = Fe2_cpx;
row.Fe3_cpx = repmat(getOptionalVar(data_cpx, ...
    'Fe3_cation_apfu', 'Cpx', NaN), nP, 1);
% Retain the legacy output name, but do not add Fe3 to the exchange Fe.
row.FeUsed_cpx = Fe2_cpx;
row.Mg_cpx = Mg_cpx;
row.Mn_cpx = repmat(getOptionalVar(data_cpx, ...
    'Mn_cation_apfu', 'Cpx', NaN), nP, 1);
row.Ca_cpx = repmat(getOptionalVar(data_cpx, ...
    'Ca_cation_apfu', 'Cpx', NaN), nP, 1);
row.Ti_cpx = repmat(getOptionalVar(data_cpx, ...
    'Ti_cation_apfu', 'Cpx', NaN), nP, 1);
row.Al_cpx = repmat(getOptionalVar(data_cpx, ...
    'Al_cation_apfu', 'Cpx', NaN), nP, 1);
row.Si_cpx = repmat(getOptionalVar(data_cpx, ...
    'Si_cation_apfu', 'Cpx', NaN), nP, 1);
row.Na_cpx = repmat(getOptionalVar(data_cpx, ...
    'Na_cation_apfu', 'Cpx', NaN), nP, 1);
row.K_cpx = repmat(getOptionalVar(data_cpx, ...
    'K_cation_apfu', 'Cpx', NaN), nP, 1);

row.FeMg_grt = FeMg_grt;
row.FeMg_cpx = FeMg_cpx;
row.KD = KD;
row.lnKD = lnKD;
row.K_Wells = K_Wells;
row.lnK_Wells = lnK_Wells;
row.numerator = numerator;
row.denominator = denominator;
row.T_K = T_K;
row.T_C = T_C;

end

function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read a required scalar, preserving NaN and zero.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end
value = tbl.(varName);
validateScalarValue(value, varName, mineralLabel);

end

function value = getOptionalVar(tbl, varName, mineralLabel, missingDefault)
% getOptionalVar
% Read an optional scalar. Existing NaN is preserved; missingDefault is used
% only when the column itself does not exist.

if ~ismember(varName, tbl.Properties.VariableNames)
    value = missingDefault;
    return;
end
value = tbl.(varName);
validateScalarValue(value, varName, mineralLabel);

end

function validateScalarValue(value, varName, mineralLabel)
% validateScalarValue
% Require a numeric scalar and prohibit negative or infinite values.

if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in a one-row table.', ...
        mineralLabel, varName);
end
if isinf(value) || (~isnan(value) && value < 0)
    error('%s.%s must be non-negative or NaN, and must not be Inf.', ...
        mineralLabel, varName);
end

end
