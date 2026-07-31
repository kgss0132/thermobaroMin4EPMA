function results = Sudholz2022GrtCpx(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/Sudholz2022GrtCpx.m
% Tested with MATLAB R2024b
%
% Garnet-clinopyroxene Fe-Mg exchange thermometer
% Sudholz, Z.J., Green, D.H., Yaxley, G.M., Jaques, A.L. (2022)
% Contributions to Mineralogy and Petrology, 177, Article 77
% DOI: https://doi.org/10.1007/s00410-022-01944-3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis and calculates temperature using the Sudholz et al.
% (2022) garnet-clinopyroxene Fe-Mg exchange thermometer.
%
% Pressure may be supplied as a scalar or vector. Therefore, this function
% can be called from both startThermoCalc_fixedP and startThermoCalc_rangeP.
% One output row is returned for every pressure value for each selected pair.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% The thermometer was calibrated using 193 peridotite, pyroxenite, and
% eclogite synthesis experiments over the following formal range:
%
%   Temperature : 850-1750 degreeC
%   Pressure    : 15-70 kbar
%   XCa_grt     : 0.08-0.37
%   XMg_grt     : 0.25-0.82
%   lnKd        : 0.24-1.84
%   Jd_cpx      : 0-0.64
%
% These ranges and the regression dataset are described in the
% Garnet-clinopyroxene Fe-Mg exchange geothermometry section on pp. 10-11.
% The calibration uncertainty is interpreted as no more than +/-75 degreeC
% (p. 12). An independent test returned a standard deviation near 101
% degreeC, so +/-75 degreeC is not a guaranteed uncertainty for every
% natural sample (pp. 12-13).
%
% Application cautions reported by Sudholz et al. (2022):
%   1) Use equilibrated, coexisting garnet and clinopyroxene. The calibration
%      dataset was screened for analytical quality and equilibrium (pp. 3-4).
%   2) Analyses used in the study generally had oxide totals of 98-102 wt%,
%      garnet totals of 7.95-8.05 cpfu on 12 O, and pyroxene totals of
%      3.98-4.02 cpfu on 6 O. Hydrous eclogitic Cpx used 3.95-4.05 cpfu
%      (pp. 3-4).
%   3) All analysed Fe was treated as Fe2+ ("Fe2+ = total Fe"; pp. 4 and
%      11). If Fe_cation_apfu and Fe3_cation_apfu are separate Fe2+ and Fe3+
%      columns, they must be summed. If Fe_cation_apfu already stores total
%      Fe, Fe3_cation_apfu must not be supplied separately.
%   4) Exercise caution when Cpx total Fe exceeds 0.35 cpfu; high Cpx FeO
%      (>10 wt%) tended to yield overestimated temperatures (p. 13).
%   5) Independent tests performed reliably near XMg_grt = 0.15-0.80,
%      Jd_cpx = 0-0.65, and XCa_grt = 0.1-0.4 (p. 13). These test ranges do
%      not replace the formal regression ranges listed above.
%   6) Natural-xenolith comparisons showed systematic discrepancies for
%      garnet CaO <4.5 wt% and Cpx Al2O3 <1 wt% (pp. 16-17).
%   7) Effects of ferric iron, H2O, K, and Na require further evaluation
%      (Conclusion, p. 17).
%
% This implementation issues non-stopping fprintf messages when pressure,
% temperature, or calculated composition indices are outside the formal
% calibration range; when Cpx Fe exceeds 0.35 cpfu; and when NaN/Inf occurs.
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
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%
% Required Clinopyroxene variables used in the equation:
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Na_cation_apfu
%
% Optional calculation variables:
%   Fe3_cation_apfu       % Garnet and Cpx; missing column means zero
%   Cr_cation_apfu        % Cpx Jd correction; missing column means zero
%   Ti_cation_apfu        % Cpx Jd correction; missing column means zero
%
% Other optional variables retained for traceability when present:
%   Mn_cation_apfu, Al_cation_apfu, Si_cation_apfu, Ca_cation_apfu (Cpx)
%
% Existing NaN values are never converted to zero. They propagate through
% the calculation and are reported by fprintf. A missing optional column is
% distinct from an existing NaN: only a missing Fe3, Cr, or Ti column is
% assigned zero. Negative or infinite calculation inputs are prohibited;
% zero and NaN are allowed.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
%   FeTotal = Fe2 + Fe3, when Fe2 and Fe3 are stored separately
%
%   Kd = (FeTotal_grt * Mg_cpx) / (FeTotal_cpx * Mg_grt)
%   XCa_grt = Ca_grt / (Ca_grt + FeTotal_grt + Mg_grt)
%   XMg_grt = Mg_grt / (Ca_grt + FeTotal_grt + Mg_grt)
%   Jd_cpx = Na_cpx - Cr_cpx - 2*Ti_cpx
%
%   T(degreeC) = 3356.34 / [ln(Kd) - 0.008*P(kbar)
%                            + 0.259*XCa_grt + 0.914*XMg_grt
%                            - 0.159*Jd_cpx + 1.265] - 273
%
% Garnet cations are normalized to 12 O and Cpx cations to 6 O.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Sudholz2022GrtCpx(rawdata_struct, P_kbar)
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
    error('Sudholz2022GrtCpx requires (rawdata_struct, P_kbar).');
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

calibrationT_min_degC = 850;
calibrationT_max_degC = 1750;
calibrationP_min_kbar = 15;
calibrationP_max_kbar = 70;

pressureOutsideCalibration = P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
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
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ' degreeC']);
    end

    if any(pressureOutsideCalibration) && ~pressureMessageIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the Sudholz et al. (2022) ' ...
             'Grt-Cpx calibration range of 15-70 kbar. %d of %d pressure ' ...
             'point(s) are outside; input range = %.4g-%.4g kbar. ' ...
             'Calculation was continued.\n'], ...
            sum(pressureOutsideCalibration), numel(P_kbar), ...
            min(P_kbar), max(P_kbar));
        pressureMessageIssued = true;
    end

    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | row.T_deg > calibrationT_max_degC);
    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the Sudholz et al. ' ...
             '(2022) Grt-Cpx calibration range of 850-1750 degreeC. ' ...
             '%d of %d finite point(s) are outside; finite range = ' ...
             '%.4g-%.4g degreeC for %s & %s. Calculation was continued.\n'], ...
            sum(temperatureOutsideCalibration), sum(finiteTemperature), ...
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

    printCompositionMessages(row, selectedCode_grt, selectedCode_cpx);

    invalidExchange = ~isfinite(row.Kd_grt_cpx) | row.Kd_grt_cpx <= 0 | ...
        ~isfinite(row.lnKd_grt_cpx) | ~isfinite(row.denominator) | ...
        row.denominator == 0;
    if any(invalidExchange)
        fprintf(2, ...
            ['WARNING: Non-finite, non-positive, or singular exchange term(s) ' ...
             'occurred for %s & %s (%d of %d pressure points). Values were ' ...
             'retained unchanged.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_cpx)), ...
            sum(invalidExchange), numel(invalidExchange));
    end

    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d). Values remain unchanged ' ...
             'in the output table.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), numel(row.T_deg), ...
            sum(isnan(row.T_deg)), sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Sudholz2022GrtCpx', ...
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

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Na_cation_apfu'};
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
    error('Sudholz2022GrtCpx: missing required variable(s): %s.', ...
        char(strjoin(missingNames(1:nMissing), ', ')));
end

end

function [inputNames, inputValues, activeInputs] = collectCalculationInputs(data_grt, data_cpx)
% collectCalculationInputs
% Collect all mandatory inputs and optional inputs that actually exist.

inputNames = ["Garnet.Fe_cation_apfu"; "Garnet.Mg_cation_apfu"; ...
    "Garnet.Ca_cation_apfu"; "Garnet.Fe3_cation_apfu"; ...
    "Cpx.Fe_cation_apfu"; "Cpx.Mg_cation_apfu"; ...
    "Cpx.Na_cation_apfu"; "Cpx.Fe3_cation_apfu"; ...
    "Cpx.Cr_cation_apfu"; "Cpx.Ti_cation_apfu"];

inputValues = NaN(10, 1);
activeInputs = false(10, 1);

inputValues(1) = getRequiredVar(data_grt, 'Fe_cation_apfu', 'Garnet');
inputValues(2) = getRequiredVar(data_grt, 'Mg_cation_apfu', 'Garnet');
inputValues(3) = getRequiredVar(data_grt, 'Ca_cation_apfu', 'Garnet');
inputValues(5) = getRequiredVar(data_cpx, 'Fe_cation_apfu', 'Cpx');
inputValues(6) = getRequiredVar(data_cpx, 'Mg_cation_apfu', 'Cpx');
inputValues(7) = getRequiredVar(data_cpx, 'Na_cation_apfu', 'Cpx');
activeInputs([1 2 3 5 6 7]) = true;

optionalRows = [4 8 9 10];
optionalTables = {data_grt, data_cpx, data_cpx, data_cpx};
optionalLabels = {'Garnet', 'Cpx', 'Cpx', 'Cpx'};
optionalNames = {'Fe3_cation_apfu', 'Fe3_cation_apfu', ...
    'Cr_cation_apfu', 'Ti_cation_apfu'};

for i = 1:numel(optionalRows)
    if ismember(optionalNames{i}, optionalTables{i}.Properties.VariableNames)
        rowIndex = optionalRows(i);
        inputValues(rowIndex) = getOptionalVar(optionalTables{i}, ...
            optionalNames{i}, optionalLabels{i}, 0);
        activeInputs(rowIndex) = true;
    end
end

end

function nanInputNames = findNaNInputs(data_grt, data_cpx)
% findNaNInputs
% List existing equation inputs containing NaN without replacing them.

[inputNames, inputValues, activeInputs] = ...
    collectCalculationInputs(data_grt, data_cpx);
nanInputNames = inputNames(activeInputs & isnan(inputValues));

end

function validateNonnegativeInputs(data_grt, data_cpx)
% validateNonnegativeInputs
% Reject negative or infinite equation inputs. Zero and NaN are allowed.

[inputNames, inputValues, activeInputs] = ...
    collectCalculationInputs(data_grt, data_cpx);
invalidInput = activeInputs & (inputValues < 0 | isinf(inputValues));
if any(invalidInput)
    error(['Sudholz2022GrtCpx: calculation inputs must be non-negative ' ...
           'and not Inf. Invalid value(s): %s.'], ...
        char(strjoin(inputNames(invalidInput), ', ')));
end

end

function row = calcTemp(data_grt, data_cpx, P_kbar)
% calcTemp
% Calculate one Garnet-Cpx pair over a scalar or vector of pressures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

Fe2_grt = repmat(getRequiredVar(data_grt, 'Fe_cation_apfu', 'Garnet'), nP, 1);
Fe3_grt = repmat(getOptionalVar(data_grt, 'Fe3_cation_apfu', 'Garnet', 0), nP, 1);
Mg_grt = repmat(getRequiredVar(data_grt, 'Mg_cation_apfu', 'Garnet'), nP, 1);
Ca_grt = repmat(getRequiredVar(data_grt, 'Ca_cation_apfu', 'Garnet'), nP, 1);

Fe2_cpx = repmat(getRequiredVar(data_cpx, 'Fe_cation_apfu', 'Cpx'), nP, 1);
Fe3_cpx = repmat(getOptionalVar(data_cpx, 'Fe3_cation_apfu', 'Cpx', 0), nP, 1);
Mg_cpx = repmat(getRequiredVar(data_cpx, 'Mg_cation_apfu', 'Cpx'), nP, 1);
Na_cpx = repmat(getRequiredVar(data_cpx, 'Na_cation_apfu', 'Cpx'), nP, 1);
Cr_cpx = repmat(getOptionalVar(data_cpx, 'Cr_cation_apfu', 'Cpx', 0), nP, 1);
Ti_cpx = repmat(getOptionalVar(data_cpx, 'Ti_cation_apfu', 'Cpx', 0), nP, 1);

FeTotal_grt = Fe2_grt + Fe3_grt;
FeTotal_cpx = Fe2_cpx + Fe3_cpx;

denGrt = Ca_grt + FeTotal_grt + Mg_grt;
Kd_grt_cpx = (FeTotal_grt .* Mg_cpx) ./ (FeTotal_cpx .* Mg_grt);
lnKd_grt_cpx = log(Kd_grt_cpx);
XCa_grt = Ca_grt ./ denGrt;
XMg_grt = Mg_grt ./ denGrt;
Jd_cpx = Na_cpx - Cr_cpx - 2 .* Ti_cpx;

denominator = lnKd_grt_cpx - 0.008 .* P_kbar ...
    + 0.259 .* XCa_grt + 0.914 .* XMg_grt ...
    - 0.159 .* Jd_cpx + 1.265;

T_K = 3356.34 ./ denominator;
T_deg = T_K - 273;

row = table();
row.P_kbar = P_kbar;

row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.FeTotal_grt = FeTotal_grt;
row.Mg_grt = Mg_grt;
row.Ca_grt = Ca_grt;
row.Mn_grt = repmat(getOptionalVar(data_grt, 'Mn_cation_apfu', 'Garnet', NaN), nP, 1);
row.Al_grt = repmat(getOptionalVar(data_grt, 'Al_cation_apfu', 'Garnet', NaN), nP, 1);
row.Si_grt = repmat(getOptionalVar(data_grt, 'Si_cation_apfu', 'Garnet', NaN), nP, 1);

row.Fe2_cpx = Fe2_cpx;
row.Fe3_cpx = Fe3_cpx;
row.FeTotal_cpx = FeTotal_cpx;
row.Mg_cpx = Mg_cpx;
row.Ca_cpx = repmat(getOptionalVar(data_cpx, 'Ca_cation_apfu', 'Cpx', NaN), nP, 1);
row.Na_cpx = Na_cpx;
row.Cr_cpx = Cr_cpx;
row.Ti_cpx = Ti_cpx;
row.Mn_cpx = repmat(getOptionalVar(data_cpx, 'Mn_cation_apfu', 'Cpx', NaN), nP, 1);
row.Al_cpx = repmat(getOptionalVar(data_cpx, 'Al_cation_apfu', 'Cpx', NaN), nP, 1);
row.Si_cpx = repmat(getOptionalVar(data_cpx, 'Si_cation_apfu', 'Cpx', NaN), nP, 1);

row.XCa_grt = XCa_grt;
row.XMg_grt = XMg_grt;
row.Jd_cpx = Jd_cpx;
row.Kd_grt_cpx = Kd_grt_cpx;
row.lnKd_grt_cpx = lnKd_grt_cpx;
row.denominator = denominator;
row.T_K = T_K;
row.T_deg = T_deg;
row.T_Sudholz2022_K = T_K;
row.T_Sudholz2022_C = T_deg;

end

function printCompositionMessages(row, selectedCode_grt, selectedCode_cpx)
% printCompositionMessages
% Report formal compositional extrapolation and the high-Fe Cpx caution.

checks = [row.XCa_grt(1) < 0.08 || row.XCa_grt(1) > 0.37; ...
          row.XMg_grt(1) < 0.25 || row.XMg_grt(1) > 0.82; ...
          row.lnKd_grt_cpx(1) < 0.24 || row.lnKd_grt_cpx(1) > 1.84; ...
          row.Jd_cpx(1) < 0 || row.Jd_cpx(1) > 0.64];
names = ["XCa_grt"; "XMg_grt"; "lnKd_grt_cpx"; "Jd_cpx"];
values = [row.XCa_grt(1); row.XMg_grt(1); ...
    row.lnKd_grt_cpx(1); row.Jd_cpx(1)];

for i = 1:numel(checks)
    if isfinite(values(i)) && checks(i)
        fprintf(2, ...
            ['WARNING: %s = %.5g is outside the formal Sudholz et al. ' ...
             '(2022) Grt-Cpx calibration composition range for %s & %s.\n'], ...
            char(names(i)), values(i), char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)));
    end
end

if isfinite(row.FeTotal_cpx(1)) && row.FeTotal_cpx(1) > 0.35
    fprintf(2, ...
        ['WARNING: Cpx total Fe = %.5g cpfu exceeds 0.35 cpfu for %s & %s. ' ...
         'Sudholz et al. (2022) recommend caution because high-Fe Cpx may ' ...
         'yield overestimated temperature.\n'], ...
        row.FeTotal_cpx(1), char(string(selectedCode_grt)), ...
        char(string(selectedCode_cpx)));
end

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
