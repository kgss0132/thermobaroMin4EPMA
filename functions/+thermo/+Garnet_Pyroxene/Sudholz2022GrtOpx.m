function results = Sudholz2022GrtOpx(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/Sudholz2022GrtOpx.m
% Tested with MATLAB R2024b
%
% Garnet-orthopyroxene Fe-Mg exchange thermometer
% Sudholz, Z.J., Green, D.H., Yaxley, G.M., Jaques, A.L. (2022)
% Contributions to Mineralogy and Petrology, 177, Article 77
% DOI: https://doi.org/10.1007/s00410-022-01944-3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Orthopyroxene analysis and calculates temperature using the Sudholz et al.
% (2022) garnet-orthopyroxene Fe-Mg exchange thermometer.
%
% Pressure may be supplied as a scalar or vector. Therefore, this function
% can be called from both startThermoCalc_fixedP and startThermoCalc_rangeP.
% One output row is returned for every pressure value for each selected pair.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% The thermometer was calibrated using 166 peridotite and pyroxenite
% synthesis experiments over the following formal range:
%
%   Temperature : 850-1525 degreeC
%   Pressure    : 16-70 kbar
%   lnKd        : 0.52-1.29
%   Lithologies : peridotite and pyroxenite, including refractory
%                 harzburgitic compositions; eclogite was not calibrated
%
% These ranges are reported in the Garnet-orthopyroxene Fe-Mg exchange
% geothermometry section on pp. 14-15. The maximum interpreted uncertainty
% is +/-100 degreeC (p. 15). The calibration dataset has a standard
% deviation of approximately 82 degreeC, and internally consistent tests
% returned approximately 73 degreeC (pp. 15-16).
%
% Application cautions reported by Sudholz et al. (2022):
%   1) Use equilibrated, coexisting garnet and orthopyroxene. The calibration
%      dataset was screened for analytical quality and equilibrium (pp. 3-4).
%   2) Analyses used in the study generally had oxide totals of 98-102 wt%,
%      garnet totals of 7.95-8.05 cpfu on 12 O, and Opx totals of
%      3.98-4.02 cpfu on 6 O (pp. 3-4).
%   3) All analysed Fe was treated as Fe2+ ("Fe2+ = total Fe"; pp. 4 and
%      15). If Fe_cation_apfu and Fe3_cation_apfu are separate Fe2+ and Fe3+
%      columns, they must be summed. If Fe_cation_apfu already stores total
%      Fe, Fe3_cation_apfu must not be supplied separately.
%   4) Independent tests performed reliably at 30-70 kbar and
%      900-1400 degreeC, but this does not replace the formal range above
%      (pp. 15-16).
%   5) Very low-Ca Opx (<0.03 cpfu Ca) more commonly produced positive
%      Delta-T values, although deviations remained within +/-150 degreeC
%      in the experimental test (p. 15).
%   6) Natural comparisons showed increasing Opx SiO2 and NiO associated
%      with larger negative temperature differences (p. 17).
%   7) Ferric-iron effects remain incompletely evaluated and are especially
%      relevant to Grt-Opx exchange; effects of H2O, K, and Na also require
%      further work (pp. 4 and 17).
%
% This implementation issues non-stopping fprintf messages when pressure,
% temperature, or lnKd is outside the formal calibration range; when Opx Ca
% is below 0.03 cpfu; and when NaN/Inf occurs.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Opx    : table
%
% The FIRST column of each table is treated as an identifier ("data code").
%
% Required Garnet variables used in the equation:
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%
% Required Orthopyroxene variables used in the equation:
%   Fe_cation_apfu
%   Mg_cation_apfu
%
% Optional calculation variable:
%   Fe3_cation_apfu       % Garnet and Opx; missing column means zero
%
% Other optional variables retained for traceability when present:
%   Mn_cation_apfu, Ca_cation_apfu (Opx), Ti_cation_apfu,
%   Al_cation_apfu, Si_cation_apfu, Na_cation_apfu, Cr_cation_apfu
%
% Existing NaN values are never converted to zero. They propagate through
% the calculation and are reported by fprintf. A missing Fe3 column is
% distinct from an existing NaN and is assigned zero. Negative or infinite
% calculation inputs are prohibited; zero and NaN are allowed.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
%   FeTotal = Fe2 + Fe3, when Fe2 and Fe3 are stored separately
%
%   Kd = (FeTotal_grt * Mg_opx) / (FeTotal_opx * Mg_grt)
%   XCa_grt = Ca_grt / (Ca_grt + FeTotal_grt + Mg_grt)
%
%   T(degreeC) = 1851.85 / [ln(Kd) - 0.007*P(kbar)
%                            - 1.83*XCa_grt + 1.08] - 273
%
% Garnet cations are normalized to 12 O and Opx cations to 6 O.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Sudholz2022GrtOpx(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Opx tables
%   P_kbar         : finite, non-negative numeric scalar or vector
%
% Output:
%   results : table containing one row per pressure value for each pair
%

%% Input validation
if nargin < 2
    error('Sudholz2022GrtOpx requires (rawdata_struct, P_kbar).');
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
if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end

dataset_grt = rawdata_struct.Garnet;
dataset_opx = rawdata_struct.Opx;
validateRequiredVariables(dataset_grt, dataset_opx);

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Table blocks are buffered and concatenated only once after the loop.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degC = 850;
calibrationT_max_degC = 1525;
calibrationP_min_kbar = 16;
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

    % ----- Orthopyroxene selection -----
    disp('=== Step 4: Selecting a data code from the list (Orthopyroxene) ===');
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

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');
    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_opx = dataset_opx(selectedIdx_opx, :);

    validateNonnegativeInputs(selectedData_grt, selectedData_opx);
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_opx);
    row = calcTemp(selectedData_grt, selectedData_opx, P_kbar);

    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_opx'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_opx)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_opx)) ...
            ': ' num2str(row.T_deg(1)) ' to ' num2str(row.T_deg(end)) ' degreeC']);
    end

    if any(pressureOutsideCalibration) && ~pressureMessageIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the Sudholz et al. (2022) ' ...
             'Grt-Opx calibration range of 16-70 kbar. %d of %d pressure ' ...
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
             '(2022) Grt-Opx calibration range of 850-1525 degreeC. ' ...
             '%d of %d finite point(s) are outside; finite range = ' ...
             '%.4g-%.4g degreeC for %s & %s. Calculation was continued.\n'], ...
            sum(temperatureOutsideCalibration), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_grt)), char(string(selectedCode_opx)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in calculation input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; calculation was not stopped.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_opx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    printCompositionMessages(row, selectedCode_grt, selectedCode_opx);

    invalidExchange = ~isfinite(row.Kd_grt_opx) | row.Kd_grt_opx <= 0 | ...
        ~isfinite(row.lnKd_grt_opx) | ~isfinite(row.denominator) | ...
        row.denominator == 0;
    if any(invalidExchange)
        fprintf(2, ...
            ['WARNING: Non-finite, non-positive, or singular exchange term(s) ' ...
             'occurred for %s & %s (%d of %d pressure points). Values were ' ...
             'retained unchanged.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_opx)), ...
            sum(invalidExchange), numel(invalidExchange));
    end

    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d). Values remain unchanged ' ...
             'in the output table.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_opx)), ...
            sum(invalidTemperature), numel(row.T_deg), ...
            sum(isnan(row.T_deg)), sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Sudholz2022GrtOpx', ...
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
function validateRequiredVariables(dataset_grt, dataset_opx)
% validateRequiredVariables
% Confirm that every mandatory equation input exists.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu'};
opxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
missingNames = strings(numel(grtVariables) + numel(opxVariables), 1);
nMissing = 0;

for i = 1:numel(grtVariables)
    if ~ismember(grtVariables{i}, dataset_grt.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Garnet." + string(grtVariables{i});
    end
end
for i = 1:numel(opxVariables)
    if ~ismember(opxVariables{i}, dataset_opx.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Opx." + string(opxVariables{i});
    end
end

if nMissing > 0
    error('Sudholz2022GrtOpx: missing required variable(s): %s.', ...
        char(strjoin(missingNames(1:nMissing), ', ')));
end

end

function [inputNames, inputValues, activeInputs] = collectCalculationInputs(data_grt, data_opx)
% collectCalculationInputs
% Collect mandatory inputs and optional Fe3 inputs that actually exist.

inputNames = ["Garnet.Fe_cation_apfu"; "Garnet.Mg_cation_apfu"; ...
    "Garnet.Ca_cation_apfu"; "Garnet.Fe3_cation_apfu"; ...
    "Opx.Fe_cation_apfu"; "Opx.Mg_cation_apfu"; ...
    "Opx.Fe3_cation_apfu"];
inputValues = NaN(7, 1);
activeInputs = false(7, 1);

inputValues(1) = getRequiredVar(data_grt, 'Fe_cation_apfu', 'Garnet');
inputValues(2) = getRequiredVar(data_grt, 'Mg_cation_apfu', 'Garnet');
inputValues(3) = getRequiredVar(data_grt, 'Ca_cation_apfu', 'Garnet');
inputValues(5) = getRequiredVar(data_opx, 'Fe_cation_apfu', 'Opx');
inputValues(6) = getRequiredVar(data_opx, 'Mg_cation_apfu', 'Opx');
activeInputs([1 2 3 5 6]) = true;

if ismember('Fe3_cation_apfu', data_grt.Properties.VariableNames)
    inputValues(4) = getOptionalVar(data_grt, ...
        'Fe3_cation_apfu', 'Garnet', 0);
    activeInputs(4) = true;
end
if ismember('Fe3_cation_apfu', data_opx.Properties.VariableNames)
    inputValues(7) = getOptionalVar(data_opx, ...
        'Fe3_cation_apfu', 'Opx', 0);
    activeInputs(7) = true;
end

end

function nanInputNames = findNaNInputs(data_grt, data_opx)
% findNaNInputs
% List existing equation inputs containing NaN without replacing them.

[inputNames, inputValues, activeInputs] = ...
    collectCalculationInputs(data_grt, data_opx);
nanInputNames = inputNames(activeInputs & isnan(inputValues));

end

function validateNonnegativeInputs(data_grt, data_opx)
% validateNonnegativeInputs
% Reject negative or infinite equation inputs. Zero and NaN are allowed.

[inputNames, inputValues, activeInputs] = ...
    collectCalculationInputs(data_grt, data_opx);
invalidInput = activeInputs & (inputValues < 0 | isinf(inputValues));
if any(invalidInput)
    error(['Sudholz2022GrtOpx: calculation inputs must be non-negative ' ...
           'and not Inf. Invalid value(s): %s.'], ...
        char(strjoin(inputNames(invalidInput), ', ')));
end

end

function row = calcTemp(data_grt, data_opx, P_kbar)
% calcTemp
% Calculate one Garnet-Opx pair over a scalar or vector of pressures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

Fe2_grt = repmat(getRequiredVar(data_grt, 'Fe_cation_apfu', 'Garnet'), nP, 1);
Fe3_grt = repmat(getOptionalVar(data_grt, 'Fe3_cation_apfu', 'Garnet', 0), nP, 1);
Mg_grt = repmat(getRequiredVar(data_grt, 'Mg_cation_apfu', 'Garnet'), nP, 1);
Ca_grt = repmat(getRequiredVar(data_grt, 'Ca_cation_apfu', 'Garnet'), nP, 1);

Fe2_opx = repmat(getRequiredVar(data_opx, 'Fe_cation_apfu', 'Opx'), nP, 1);
Fe3_opx = repmat(getOptionalVar(data_opx, 'Fe3_cation_apfu', 'Opx', 0), nP, 1);
Mg_opx = repmat(getRequiredVar(data_opx, 'Mg_cation_apfu', 'Opx'), nP, 1);

FeTotal_grt = Fe2_grt + Fe3_grt;
FeTotal_opx = Fe2_opx + Fe3_opx;
denGrt = Ca_grt + FeTotal_grt + Mg_grt;
Kd_grt_opx = (FeTotal_grt .* Mg_opx) ./ (FeTotal_opx .* Mg_grt);
lnKd_grt_opx = log(Kd_grt_opx);
XCa_grt = Ca_grt ./ denGrt;

denominator = lnKd_grt_opx - 0.007 .* P_kbar ...
    - 1.83 .* XCa_grt + 1.08;
T_K = 1851.85 ./ denominator;
T_deg = T_K - 273;

row = table();
row.P_kbar = P_kbar;

row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.FeTotal_grt = FeTotal_grt;
row.Mg_grt = Mg_grt;
row.Ca_grt = Ca_grt;
row.Mn_grt = repmat(getOptionalVar(data_grt, 'Mn_cation_apfu', 'Garnet', NaN), nP, 1);
row.Ti_grt = repmat(getOptionalVar(data_grt, 'Ti_cation_apfu', 'Garnet', NaN), nP, 1);
row.Al_grt = repmat(getOptionalVar(data_grt, 'Al_cation_apfu', 'Garnet', NaN), nP, 1);
row.Si_grt = repmat(getOptionalVar(data_grt, 'Si_cation_apfu', 'Garnet', NaN), nP, 1);
row.Na_grt = repmat(getOptionalVar(data_grt, 'Na_cation_apfu', 'Garnet', NaN), nP, 1);
row.Cr_grt = repmat(getOptionalVar(data_grt, 'Cr_cation_apfu', 'Garnet', NaN), nP, 1);

row.Fe2_opx = Fe2_opx;
row.Fe3_opx = Fe3_opx;
row.FeTotal_opx = FeTotal_opx;
row.Mg_opx = Mg_opx;
row.Ca_opx = repmat(getOptionalVar(data_opx, 'Ca_cation_apfu', 'Opx', NaN), nP, 1);
row.Mn_opx = repmat(getOptionalVar(data_opx, 'Mn_cation_apfu', 'Opx', NaN), nP, 1);
row.Ti_opx = repmat(getOptionalVar(data_opx, 'Ti_cation_apfu', 'Opx', NaN), nP, 1);
row.Al_opx = repmat(getOptionalVar(data_opx, 'Al_cation_apfu', 'Opx', NaN), nP, 1);
row.Si_opx = repmat(getOptionalVar(data_opx, 'Si_cation_apfu', 'Opx', NaN), nP, 1);
row.Na_opx = repmat(getOptionalVar(data_opx, 'Na_cation_apfu', 'Opx', NaN), nP, 1);
row.Cr_opx = repmat(getOptionalVar(data_opx, 'Cr_cation_apfu', 'Opx', NaN), nP, 1);

row.XCa_grt = XCa_grt;
row.Kd_grt_opx = Kd_grt_opx;
row.lnKd_grt_opx = lnKd_grt_opx;
row.denominator = denominator;
row.T_K = T_K;
row.T_deg = T_deg;
row.T_Sudholz2022_K = T_K;
row.T_Sudholz2022_C = T_deg;

end

function printCompositionMessages(row, selectedCode_grt, selectedCode_opx)
% printCompositionMessages
% Report lnKd extrapolation and the low-Ca Opx caution.

lnKdValue = row.lnKd_grt_opx(1);
if isfinite(lnKdValue) && (lnKdValue < 0.52 || lnKdValue > 1.29)
    fprintf(2, ...
        ['WARNING: lnKd_grt_opx = %.5g is outside the formal Sudholz et al. ' ...
         '(2022) calibration range of 0.52-1.29 for %s & %s.\n'], ...
        lnKdValue, char(string(selectedCode_grt)), ...
        char(string(selectedCode_opx)));
end

opxCa = row.Ca_opx(1);
if isfinite(opxCa) && opxCa < 0.03
    fprintf(2, ...
        ['WARNING: Opx Ca = %.5g cpfu is below 0.03 cpfu for %s & %s. ' ...
         'Sudholz et al. (2022) observed more frequent positive Delta-T ' ...
         'values for very low-Ca Opx.\n'], ...
        opxCa, char(string(selectedCode_grt)), ...
        char(string(selectedCode_opx)));
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
