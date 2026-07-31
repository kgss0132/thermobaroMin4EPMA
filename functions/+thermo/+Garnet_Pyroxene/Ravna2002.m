function results = Ravna2002(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/Ravna2002.m
% Tested with MATLAB R2024b
%
% Garnet-clinopyroxene Fe2+-Mg exchange thermometer
% Ravna, E.K. (2002)
% Journal of Metamorphic Geology, 18, 211-219
% DOI: https://doi.org/10.1046/j.1525-1314.2000.00247.x
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis and calculates temperature using the Ravna (2000)
% garnet-clinopyroxene Fe2+-Mg exchange thermometer.
%
% The function accepts pressure as either a scalar or a vector. Therefore,
% it can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for every pressure
% value for each user-selected Garnet-Clinopyroxene pair.
%
% The function is designed for repeated calculations: after each run it asks
% whether another Garnet-Clinopyroxene pair should be calculated and stores
% all result blocks in a single output table.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Ravna (2000) derived the calibration by multiple regression of 311
% experimental Garnet-Clinopyroxene pairs and 49 natural Mn-rich pairs
% (n = 360; r2 = 0.98). Of 404 candidate experimental pairs, 93 points
% suspected to represent disequilibrium were excluded from the regression
% (Results and Table 3, p. 213).
%
% Calibration-data ranges:
%   Temperature, experimental data : 600-1740 degreeC
%   Pressure, experimental data    : 1.0-6.0 GPa (10-60 kbar)
%   Pressure, all regression data  : 0.7-6.0 GPa (7-60 kbar), because
%                                    natural Mn-rich pairs at 0.7 GPa were
%                                    also included
%   XCa_grt                        : 0.09-0.47
%   XMgNum_grt                     : 0.035-0.89
%   XMn_grt                        : 0-0.526
%   XNa_cpx                        : 0-0.51
%
% The experimental P-T ranges are reported in the Dataset section and
% Table 1 (p. 212). The compositional ranges are reported in Table 2
% (p. 213). The equation reproduces the experimental runs retained in the
% regression to within approximately +/-100 degreeC (p. 213). This is a
% calibration residual, not a complete uncertainty estimate for a natural
% sample.
%
% Important application cautions ("A cautionary note", pp. 217-218):
%   1) The thermometer assumes equilibrium between coexisting garnet and
%      clinopyroxene. Cooling may reset peak compositions by Fe-Mg exchange
%      and by combined net-transfer and exchange reactions.
%   2) KD is defined with Fe2+, not total Fe. Uncertainty in Fe2+/Fe3+ may
%      strongly affect the result, especially for Fe-poor garnet and Fe-poor
%      omphacite or diopside.
%   3) Stoichiometric Fe3+ estimates are sensitive to analytical error.
%      Overestimating Fe3+ in garnet lowers KD and increases calculated T.
%   4) No significant Na effect was recognized within XNa_cpx = 0-0.51
%      (abstract, p. 211; Results and Fig. 1, pp. 213-214), but this does not
%      remove the Fe2+/Fe3+ problem for Fe-poor omphacite.
%   5) The Pattison and Newton (1989) dataset was excluded because of
%      suspected systematic inconsistencies with the other datasets
%      (Dataset section, p. 212).
%
% This implementation issues non-stopping messages with fprintf when:
%   1) input pressure is outside 0.7-6.0 GPa (7-60 kbar),
%   2) a finite calculated temperature is outside 600-1740 degreeC,
%   3) any required thermometer input is NaN,
%   4) KD, ln(KD), the equation denominator, or temperature is non-finite.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Cpx    : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns contain
% normalized cation data.
%
% Required Garnet variables used by the thermometer:
%   Fe_cation_apfu         % Fe2+, not total Fe
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Mn_cation_apfu
%
% Required Clinopyroxene variables used by the thermometer:
%   Fe_cation_apfu         % Fe2+, not total Fe
%   Mg_cation_apfu
%
% Optional variables retained in the output when present:
%   Fe3_cation_apfu
%   Ti_cation_apfu
%   Cr_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%   Ca_cation_apfu         % Clinopyroxene only; not used in the equation
%   Mn_cation_apfu         % Clinopyroxene only; not used in the equation
%
% Required inputs may be NaN. NaN is retained and propagated through the
% calculation; it is never replaced by zero. A non-stopping message lists
% the exact NaN input variable(s). Infinite or negative required inputs are
% prohibited. Zero is accepted, but may produce NaN or Inf through the
% published algebra and will then trigger a non-stopping result message.
%
% Fe3_cation_apfu is not added to Fe_cation_apfu. Ravna (2000) explicitly
% defines KD and the garnet mole fractions using Fe2+.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
%   KD = (Fe2+/Mg)_grt / (Fe2+/Mg)_cpx
%
%   XCa_grt    = Ca / (Ca + Mn + Fe2+ + Mg)
%   XMn_grt    = Mn / (Ca + Mn + Fe2+ + Mg)
%   XMgNum_grt = Mg / (Mg + Fe2+)
%   P_GPa      = P_kbar / 10
%
%   T(K) = [1939.9 + 3270*XCa_grt - 1396*XCa_grt^2 ...
%           + 3319*XMn_grt - 3535*XMn_grt^2 ...
%           + 1105*XMgNum_grt - 3561*XMgNum_grt^2 ...
%           + 2324*XMgNum_grt^3 + 169.4*P_GPa] ...
%          / [ln(KD) + 1.223]
%
%   T(degreeC) = T(K) - 273.15
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Ravna2000(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Cpx tables
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Clinopyroxene pair
%

%% Input validation
% Basic argument checks prevent silent failure due to missing inputs or an
% invalid pressure array.
if nargin < 2
    error('Ravna2000 requires (rawdata_struct, P_kbar).');
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
% Extract the required tables without changing the input data.
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
% Result tables are buffered and concatenated only once after the selection
% loop. This avoids repeatedly reallocating and copying the full table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Regression-supported P-T limits from Ravna (2000).
calibrationT_min_degC = 600;
calibrationT_max_degC = 1740;
calibrationP_min_GPa = 0.7;
calibrationP_max_GPa = 6.0;

% Pressure is common to all selected pairs and is therefore reported only
% once per function call.
P_GPa_input = P_kbar ./ 10;
pressureOutsideCalibration = ...
    P_GPa_input < calibrationP_min_GPa | ...
    P_GPa_input > calibrationP_max_GPa;
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

    % NaN values are listed but deliberately allowed to propagate.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_cpx);

    % Negative and infinite calculation inputs stop the calculation. Zero
    % and NaN are allowed according to the requested handling rules.
    validateNonnegativeInputs(selectedData_grt, selectedData_cpx);

    row = calcTemp(selectedData_grt, selectedData_cpx, P_kbar);

    % Repeat identifiers so that vector pressure input creates a complete
    % and traceable output row for every pressure value.
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_cpx'}, 'Before', 1);

    % Store one result-table block. The buffer grows geometrically only when
    % its preallocated capacity is exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated result for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    % Report pressure extrapolation only once per function call.
    if any(pressureOutsideCalibration) && ~pressureMessageIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the regression-supported range ' ...
             'of Ravna (2000): 0.7-6.0 GPa (7-60 kbar). ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g GPa. Calculation was continued.\n'], ...
            sum(pressureOutsideCalibration), numel(P_GPa_input), ...
            min(P_GPa_input), max(P_GPa_input));
        pressureMessageIssued = true;
    end

    % Report finite temperatures outside the experimental temperature range.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the experimental range ' ...
             'of Ravna (2000): 600-1740 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g-%.4g degreeC for %s & %s. ' ...
             'Calculation was continued.\n'], ...
            sum(temperatureOutsideCalibration), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_grt)), char(string(selectedCode_cpx)));
    end

    % List the exact required input variables containing NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; calculation was not stopped.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Report undefined or singular exchange terms without altering them.
    invalidExchange = ~isfinite(row.KD) | row.KD <= 0 | ...
        ~isfinite(row.lnKD) | ~isfinite(row.denominator) | ...
        row.denominator == 0;
    if any(invalidExchange)
        fprintf(2, ...
            ['WARNING: Non-finite, non-positive, or singular exchange term(s) ' ...
             'occurred for %s & %s (%d of %d pressure points).\n' ...
             '         KD, lnKD, denominator, and calculated results were retained unchanged.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_cpx)), ...
            sum(invalidExchange), numel(invalidExchange));
    end

    % Retain and report every NaN or Inf result.
    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         Values remain unchanged in the output table; calculation was not stopped.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), numel(row.T_deg), ...
            sum(isnan(row.T_deg)), sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Ravna2000', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate only once after all interactive selections are complete.
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
% Confirm that all variables used by the published equation exist.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
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
    missingNames = missingNames(1:nMissing);
    error('Ravna2000: missing required variable(s): %s.', ...
        char(strjoin(missingNames, ', ')));
end

end

function nanInputNames = findNaNInputs(data_grt, data_cpx)
% findNaNInputs
% Return the exact required thermometer input names containing NaN.

inputNames = ["Garnet.Fe_cation_apfu"; ...
              "Garnet.Mg_cation_apfu"; ...
              "Garnet.Ca_cation_apfu"; ...
              "Garnet.Mn_cation_apfu"; ...
              "Cpx.Fe_cation_apfu"; ...
              "Cpx.Mg_cation_apfu"];

inputValues = [data_grt.Fe_cation_apfu; ...
               data_grt.Mg_cation_apfu; ...
               data_grt.Ca_cation_apfu; ...
               data_grt.Mn_cation_apfu; ...
               data_cpx.Fe_cation_apfu; ...
               data_cpx.Mg_cation_apfu];

nanInputNames = inputNames(isnan(inputValues));

end

function validateNonnegativeInputs(data_grt, data_cpx)
% validateNonnegativeInputs
% Reject negative or infinite required calculation inputs. NaN and zero are
% intentionally allowed so that the published algebra determines the result.

inputNames = ["Garnet.Fe_cation_apfu"; ...
              "Garnet.Mg_cation_apfu"; ...
              "Garnet.Ca_cation_apfu"; ...
              "Garnet.Mn_cation_apfu"; ...
              "Cpx.Fe_cation_apfu"; ...
              "Cpx.Mg_cation_apfu"];

inputValues = [data_grt.Fe_cation_apfu; ...
               data_grt.Mg_cation_apfu; ...
               data_grt.Ca_cation_apfu; ...
               data_grt.Mn_cation_apfu; ...
               data_cpx.Fe_cation_apfu; ...
               data_cpx.Mg_cation_apfu];

invalidInput = inputValues < 0 | isinf(inputValues);

if any(invalidInput)
    invalidNames = inputNames(invalidInput);
    error(['Ravna2000: required calculation inputs must be non-negative ' ...
           'and not Inf. Invalid value(s) were found in: %s.'], ...
        char(strjoin(invalidNames, ', ')));
end

end

function row = calcTemp(data_grt, data_cpx, P_kbar)
% calcTemp
% Compute Ravna (2000) temperatures for one garnet row, one clinopyroxene
% row, and a scalar or vector of pressure values.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

% Extract the cations used by the thermometer. Fe_cation_apfu is treated as
% Fe2+ and Fe3_cation_apfu is not added to it.
Fe2_grt = repmat(getRequiredVar(data_grt, 'Fe_cation_apfu', 'Garnet'), nP, 1);
Mg_grt  = repmat(getRequiredVar(data_grt, 'Mg_cation_apfu', 'Garnet'), nP, 1);
Ca_grt  = repmat(getRequiredVar(data_grt, 'Ca_cation_apfu', 'Garnet'), nP, 1);
Mn_grt  = repmat(getRequiredVar(data_grt, 'Mn_cation_apfu', 'Garnet'), nP, 1);

Fe2_cpx = repmat(getRequiredVar(data_cpx, 'Fe_cation_apfu', 'Clinopyroxene'), nP, 1);
Mg_cpx  = repmat(getRequiredVar(data_cpx, 'Mg_cation_apfu', 'Clinopyroxene'), nP, 1);

% Optional values are retained for traceability. Missing optional variables
% remain NaN rather than being silently replaced by zero.
Fe3_grt = repmat(getOptionalVar(data_grt, 'Fe3_cation_apfu', 'Garnet'), nP, 1);
Al_grt  = repmat(getOptionalVar(data_grt, 'Al_cation_apfu', 'Garnet'), nP, 1);
Si_grt  = repmat(getOptionalVar(data_grt, 'Si_cation_apfu', 'Garnet'), nP, 1);

Fe3_cpx = repmat(getOptionalVar(data_cpx, 'Fe3_cation_apfu', 'Clinopyroxene'), nP, 1);
Mn_cpx  = repmat(getOptionalVar(data_cpx, 'Mn_cation_apfu', 'Clinopyroxene'), nP, 1);
Ca_cpx  = repmat(getOptionalVar(data_cpx, 'Ca_cation_apfu', 'Clinopyroxene'), nP, 1);
Ti_cpx  = repmat(getOptionalVar(data_cpx, 'Ti_cation_apfu', 'Clinopyroxene'), nP, 1);
Cr_cpx  = repmat(getOptionalVar(data_cpx, 'Cr_cation_apfu', 'Clinopyroxene'), nP, 1);
Al_cpx  = repmat(getOptionalVar(data_cpx, 'Al_cation_apfu', 'Clinopyroxene'), nP, 1);
Si_cpx  = repmat(getOptionalVar(data_cpx, 'Si_cation_apfu', 'Clinopyroxene'), nP, 1);
Na_cpx  = repmat(getOptionalVar(data_cpx, 'Na_cation_apfu', 'Clinopyroxene'), nP, 1);
K_cpx   = repmat(getOptionalVar(data_cpx, 'K_cation_apfu', 'Clinopyroxene'), nP, 1);

% Apply the published algebra directly. NaN and Inf are not replaced.
FeMg_grt = Fe2_grt ./ Mg_grt;
FeMg_cpx = Fe2_cpx ./ Mg_cpx;
KD = FeMg_grt ./ FeMg_cpx;
lnKD = log(KD);

xSiteSum_grt = Fe2_grt + Mg_grt + Mn_grt + Ca_grt;
XCa_grt = Ca_grt ./ xSiteSum_grt;
XMn_grt = Mn_grt ./ xSiteSum_grt;
XMgNum_grt = Mg_grt ./ (Fe2_grt + Mg_grt);

numerator = 1939.9 ...
    + 3270 .* XCa_grt ...
    - 1396 .* (XCa_grt .^ 2) ...
    + 3319 .* XMn_grt ...
    - 3535 .* (XMn_grt .^ 2) ...
    + 1105 .* XMgNum_grt ...
    - 3561 .* (XMgNum_grt .^ 2) ...
    + 2324 .* (XMgNum_grt .^ 3) ...
    + 169.4 .* P_GPa;

denominator = lnKD + 1.223;
T_K = numerator ./ denominator;
T_deg = T_K - 273.15;

% Pack outputs. FeUsed is retained as a compatibility field, but it is
% explicitly equal to Fe2+ rather than Fe2+ + Fe3+.
row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

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
row.Cr_cpx = Cr_cpx;
row.Al_cpx = Al_cpx;
row.Si_cpx = Si_cpx;
row.Na_cpx = Na_cpx;
row.K_cpx = K_cpx;

row.FeMg_grt = FeMg_grt;
row.FeMg_cpx = FeMg_cpx;
row.XCa_grt = XCa_grt;
row.XMn_grt = XMn_grt;
row.XMgNum_grt = XMgNum_grt;
row.KD = KD;
row.lnKD = lnKD;
row.numerator = numerator;
row.denominator = denominator;
row.T_K = T_K;
row.T_deg = T_deg;
row.T_C = T_deg; % Legacy alias retained for compatibility.

end

function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read one required scalar value. NaN and zero are retained; negative and
% infinite values are prohibited.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);
if ~isscalar(value)
    error('%s.%s must be scalar in a one-row table.', mineralLabel, varName);
end
if isinf(value) || (~isnan(value) && value < 0)
    error('%s.%s must be non-negative or NaN, and must not be Inf.', ...
        mineralLabel, varName);
end

end

function value = getOptionalVar(tbl, varName, mineralLabel)
% getOptionalVar
% Read one optional scalar value. A missing variable is represented by NaN;
% an existing NaN is preserved and never replaced by zero.

if ~ismember(varName, tbl.Properties.VariableNames)
    value = NaN;
    return;
end

value = tbl.(varName);
if ~isscalar(value)
    error('%s.%s must be scalar in a one-row table.', mineralLabel, varName);
end
if isinf(value) || (~isnan(value) && value < 0)
    error('%s.%s must be non-negative or NaN, and must not be Inf.', ...
        mineralLabel, varName);
end

end
