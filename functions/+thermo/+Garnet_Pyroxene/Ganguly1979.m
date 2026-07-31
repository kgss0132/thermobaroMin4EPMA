function results = Ganguly1979(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/Ganguly1979.m
% Tested with MATLAB R2024b
%
% Garnet-clinopyroxene Fe2+-Mg exchange thermometer
% Ganguly, J. (1979)
% Geochimica et Cosmochimica Acta, 43, 1021-1029
% DOI: https://doi.org/10.1016/0016-7037(79)90091-7
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis (selected by the user from tables) and calculates
% temperature using the Ganguly (1979) garnet-clinopyroxene Fe2+-Mg
% exchange thermometer.
%
% Pressure may be supplied as either a scalar or a vector. Therefore, this
% function can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. For each selected Grt-Cpx pair, the output table
% contains one row per pressure value.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Grt-Cpx pair, and appends results into
% a single output table.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Ganguly (1979) is a thermodynamic and semi-empirical formulation rather
% than a single broad experimental calibration. Its principal high-P/T
% experimental anchor is the Wood (1976) Grt-Cpx dataset (p. 1024):
%
%   Temperature : 1100-1400 degreeC
%   Pressure    : 20-45 kbar
%
% The model was also compared with Na-poor natural and experimental Grt-Cpx
% pairs over a broader range. Ganguly states that, within the predominantly
% six-component system considered, temperatures in the following range
% should be reliable to approximately +/-35 degreeC (pp. 1027-1028):
%
%   Temperature : approximately 800-1400 degreeC
%   Pressure    : approximately 10-40.5 kbar in the comparison dataset
%
% The paper does not define a formal hard pressure-calibration interval.
% This implementation therefore uses 800-1400 degreeC and 10-45 kbar as a
% practical published application/comparison range for non-stopping
% warnings. The 10-45 kbar interval combines the comparison dataset with the
% direct 20-45 kbar experimental anchor and is not a formal rectangular
% calibration envelope.
%
% Important application notes from Ganguly (1979):
%   1) The thermometer must currently be restricted to Na-poor bulk
%      compositions because adequate mixing data for jadeite with diopside
%      and hedenbergite were unavailable (abstract, p. 1021; p. 1025).
%      The paper gives no numerical Na cutoff. Na-rich omphacitic Cpx should
%      therefore be treated as outside the demonstrated model.
%   2) The thermodynamic garnet model is based mainly on grossular and
%      spessartine contents up to approximately 30 mol% each. Application to
%      garnets significantly beyond those demonstrated ranges should be made
%      with caution (pp. 1023 and 1027).
%   3) In the Kelley Iron Formation comparison, Grt spans approximately
%      17-34 mol% Ca and 1-32 mol% Mn (p. 1026). For very Mn-rich
%      metasediments, the neglected influence of Mn in Cpx may become
%      important, but the required Cpx mixing data were unavailable
%      (p. 1023).
%   4) The natural and experimental comparison set in Table 3 contains Grt
%      with approximately 14-25 mol% Ca and Cpx with as much as 23 mol% CaTs
%      (pp. 1026-1027).
%   5) The model should not be applied without additional corrections when
%      Grt or Cpx contains significant components outside the system treated
%      in the paper. The possible influence of Cpx Ca was not explicitly
%      included (p. 1027).
%   6) KD is defined using Fe2+, not total Fe. Treating all probe Fe as FeO
%      introduces an additional temperature error when Fe3+ is significant
%      (p. 1028). Fe3_cation_apfu is therefore retained for reference in this
%      implementation but is not added to Fe_cation_apfu.
%   7) Approximately 5 mol% (Fe+Cr)3+ in the six-fold garnet site did not
%      appear to have a serious effect in the limited comparison dataset,
%      but this does not remove the requirement to use Fe2+ in KD (p. 1027).
%   8) The assumed weak dependence of the corrected equilibrium function on
%      Mg/(Mg+Fe2+) was supported over approximately 0.06-0.85, but Ganguly
%      notes that the underlying assumption was not adequately justified and
%      was tested mainly through its predictive consequences (p. 1024).
%   9) Uncertainty may exceed +/-35 degreeC below approximately 800 degreeC
%      because errors in the low-temperature extrapolation are magnified.
%      Better constrained low-temperature experiments were explicitly
%      identified as necessary (p. 1028).
%  10) Equilibrium mineral compositions are required. Wood's 1100 degreeC
%      data included a reversed determination, whereas the longer high-
%      temperature runs were judged adequate for equilibrium (p. 1024).
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside 10-45 kbar,
%   2) a finite calculated temperature is outside 800-1400 degreeC,
%   3) a required thermometer input contains NaN, or
%   4) a calculated temperature is NaN or Inf. In this case, the actual
%      thermometer inputs and identified calculation-domain causes are also
%      displayed.
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
% Variables used directly by the thermometer:
%   Garnet table:
%     Fe_cation_apfu         % Fe2+ in garnet
%     Mg_cation_apfu
%     Ca_cation_apfu
%     Mn_cation_apfu
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
%   Ca_cation_apfu and Mn_cation_apfu in Cpx
%
% IMPORTANT Fe note:
% The exchange reaction and Eqs. (10a-b) are defined using Fe2+. Therefore,
% Fe3_cation_apfu is retained for reference but is not added to
% Fe_cation_apfu. If Fe_cation_apfu contains total Fe rather than Fe2+, Fe2+
% should be estimated before this function is used.
%
% Negative finite values in variables used by the thermometer are not
% permitted. Zero values are retained; if they make the equation
% mathematically undefined, the affected result is returned as NaN and a
% non-stopping warning is printed. NaN values are never replaced by zero;
% they propagate through the calculation and remain in the output.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Fe2+-Mg distribution coefficient:
%
%   KD = (Fe2+/Mg)_Grt / (Fe2+/Mg)_Cpx
%
% Garnet divalent-site mole fractions:
%
%   XCa_Gt = Ca / (Fe2+ + Mg + Mn + Ca)
%   XMn_Gt = Mn / (Fe2+ + Mg + Mn + Ca)
%
% Ganguly (1979) gives two linear approximations (Eqs. 10a-b, pp. 1024-
% 1025). After neglecting the jadeite mixing term for Na-poor compositions
% and solving each expression for T:
%
% For T >= 1333 K (Eq. 10a):
%
%          4100 + 11.07*P_kbar + 1586*XCa_Gt + 1308*XMn_Gt
%   T(K) = ------------------------------------------------------
%                             ln(KD) + 2.400
%
% For T <= 1333 K (Eq. 10b):
%
%          4801 + 11.07*P_kbar + 1586*XCa_Gt + 1308*XMn_Gt
%   T(K) = ------------------------------------------------------
%                             ln(KD) + 2.930
%
% Both candidate temperatures are calculated. If only one is consistent
% with its stated temperature interval, that branch is used. If the choice
% is not evident, Ganguly instructs the user to accept the expression that
% yields the higher temperature (p. 1025); this rule is applied point by
% point. The selected branch and both candidates are retained in the output.
%
% The previous implementation used 4801 with ln(KD)+1.9034. The constant
% 1.9034 does not belong to Ganguly's Eqs. (10a-b) and is not used here.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Ganguly1979(rawdata_struct, P_kbar)
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
if nargin < 2
    error('Ganguly1979 requires (rawdata_struct, P_kbar).');
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
% Result blocks are buffered and concatenated once after the interactive
% loop, avoiding repeated reallocation of the complete output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Practical application/comparison limits from Ganguly (1979).
applicationT_min_degC = 800;
applicationT_max_degC = 1400;
applicationP_min_kbar = 10;
applicationP_max_kbar = 45;

pressureOutsideApplication = ...
    P_kbar < applicationP_min_kbar | P_kbar > applicationP_max_kbar;
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

    nanInputNames = findNaNInputs(selectedData_grt, selectedData_cpx);
    validateNonNegativeInputs(selectedData_grt, selectedData_cpx);

    row = calcTemp(selectedData_grt, selectedData_cpx, P_kbar);

    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_cpx'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo calculated temperature values.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_C(1)) ' to ' num2str(row.T_C(end)) ' degreeC']);
    end

    % Pressure warning is printed once because the same pressure vector is
    % used for every mineral pair in this function call.
    if any(pressureOutsideApplication) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the practical published ' ...
             'application/comparison range used here for Ganguly (1979): ' ...
             '10-45 kbar. %d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar. The direct Wood (1976) ' ...
             'experimental anchor is 20-45 kbar, and Ganguly did not define ' ...
             'a formal rectangular pressure-calibration limit (pp. 1024, 1027).\n'], ...
            sum(pressureOutsideApplication), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn for finite temperatures outside the stated reliability interval.
    finiteTemperature = isfinite(row.T_C);
    temperatureOutsideApplication = finiteTemperature & ...
        (row.T_C < applicationT_min_degC | ...
         row.T_C > applicationT_max_degC);

    if any(temperatureOutsideApplication)
        finiteValues = row.T_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximate ' ...
             '800-1400 degreeC reliability range discussed by Ganguly (1979). ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g-%.4g degreeC for %s & %s. ' ...
             'Uncertainty may be greater below 800 degreeC because of ' ...
             'low-temperature extrapolation (pp. 1027-1028).\n'], ...
            sum(temperatureOutsideApplication), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)));
    end

    % Report explicitly stored NaN thermometer inputs. Calculation continues
    % and NaN values are not replaced by zero.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, and NaN was not replaced by zero.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Preserve non-finite results and print both raw inputs and identified
    % calculation-domain causes.
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
                        'was identified; inspect the stored branch candidates and intermediate variables.\n']);
        else
            fprintf(2, '         Identified cause(s): %s.\n', ...
                char(strjoin(nonFiniteCauses, '; ')));
        end
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Ganguly1979', ...
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
% Verify all columns required by the Ganguly (1979) equation.

requiredGrt = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
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
    error(['Ganguly1979: required variable(s) are missing: ' ...
        char(strjoin(missingNames, ', ')) '.']);
end

end


function nanInputNames = findNaNInputs(data_grt, data_cpx)
% findNaNInputs
% Return names of required thermometer inputs containing NaN.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

nanInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nNaN = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    variableValue = data_grt.(variableName);
    if any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanInputNames(nNaN) = "Garnet." + string(variableName);
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
% Stop when a finite required thermometer input is negative. Zero and NaN
% remain available to the calculation and result diagnostics.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

invalidInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nInvalid = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    variableValue = data_grt.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Garnet." + string(variableName);
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
    error(['Ganguly1979: thermometer inputs must be >= 0. ' ...
        'Negative finite value(s) were found in: ' ...
        char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function row = calcTemp(data_grt, data_cpx, P_kbar)
% calcTemp
% Compute Ganguly (1979) Grt-Cpx temperatures for one mineral pair and a
% scalar or vector of pressures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;

% --- Extract and expand garnet cations ---
Fe2_grt = repmat(getRequiredValue(data_grt, 'Fe_cation_apfu', 'Garnet'), nP, 1);
Fe3_grt = repmat(getOptionalValue(data_grt, 'Fe3_cation_apfu', 0, 'Garnet'), nP, 1);
Mg_grt  = repmat(getRequiredValue(data_grt, 'Mg_cation_apfu', 'Garnet'), nP, 1);
Mn_grt  = repmat(getRequiredValue(data_grt, 'Mn_cation_apfu', 'Garnet'), nP, 1);
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

% --- Fe2+-Mg distribution coefficient ---
% Fe3+ is deliberately excluded from the exchange coefficient.
FeMg_grt = Fe2_grt ./ Mg_grt;
FeMg_cpx = Fe2_cpx ./ Mg_cpx;
KD = FeMg_grt ./ FeMg_cpx;
lnKD = log(KD);

% --- Garnet divalent-site mole fractions ---
xSiteSum_grt = Fe2_grt + Mg_grt + Mn_grt + Ca_grt;
XCa_grt = Ca_grt ./ xSiteSum_grt;
XMn_grt = Mn_grt ./ xSiteSum_grt;

% --- Ganguly (1979), Eqs. (10a-b), pp. 1024-1025 ---
highT_numerator = 4100 + 11.07 .* P_kbar ...
    + 1586 .* XCa_grt + 1308 .* XMn_grt;
highT_denominator = lnKD + 2.400;
T_highCandidate_K = highT_numerator ./ highT_denominator;

lowT_numerator = 4801 + 11.07 .* P_kbar ...
    + 1586 .* XCa_grt + 1308 .* XMn_grt;
lowT_denominator = lnKD + 2.930;
T_lowCandidate_K = lowT_numerator ./ lowT_denominator;

% Both branch candidates must first satisfy the shared mathematical domain.
mathValid = isfinite(KD) & KD > 0 ...
    & isfinite(XCa_grt) & XCa_grt >= 0 ...
    & isfinite(XMn_grt) & XMn_grt >= 0 ...
    & isfinite(highT_numerator) & isfinite(highT_denominator) ...
    & isfinite(lowT_numerator) & isfinite(lowT_denominator) ...
    & highT_denominator > 0 & lowT_denominator > 0 ...
    & isfinite(T_highCandidate_K) & isfinite(T_lowCandidate_K);

branchBoundary_K = 1333;
highBranchConsistent = mathValid & T_highCandidate_K >= branchBoundary_K;
lowBranchConsistent = mathValid & T_lowCandidate_K <= branchBoundary_K;

useHighOnly = highBranchConsistent & ~lowBranchConsistent;
useLowOnly = lowBranchConsistent & ~highBranchConsistent;
ambiguousBranch = mathValid & ...
    ((highBranchConsistent & lowBranchConsistent) | ...
     (~highBranchConsistent & ~lowBranchConsistent));
useHighWhenAmbiguous = ambiguousBranch & ...
    T_highCandidate_K >= T_lowCandidate_K;
useLowWhenAmbiguous = ambiguousBranch & ...
    T_lowCandidate_K > T_highCandidate_K;

T_K = NaN(nP, 1);
branchUsed = repmat("invalid", nP, 1);

T_K(useHighOnly) = T_highCandidate_K(useHighOnly);
branchUsed(useHighOnly) = "Eq10a_highT";

T_K(useLowOnly) = T_lowCandidate_K(useLowOnly);
branchUsed(useLowOnly) = "Eq10b_lowT";

% Ganguly's instruction for an unclear branch choice is to retain the
% expression that gives the higher temperature.
T_K(useHighWhenAmbiguous) = T_highCandidate_K(useHighWhenAmbiguous);
branchUsed(useHighWhenAmbiguous) = "Eq10a_higher_when_ambiguous";

T_K(useLowWhenAmbiguous) = T_lowCandidate_K(useLowWhenAmbiguous);
branchUsed(useLowWhenAmbiguous) = "Eq10b_higher_when_ambiguous";

T_C = T_K - 273.15;

% --- Pack outputs ---
% FeUsed is retained for compatibility with the original output variable
% set, but now correctly equals Fe2+ rather than Fe2+ + Fe3+.
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
row.XMn_grt = XMn_grt;
row.KD = KD;
row.lnKD = lnKD;
row.highT_numerator = highT_numerator;
row.highT_denominator = highT_denominator;
row.T_highCandidate_K = T_highCandidate_K;
row.lowT_numerator = lowT_numerator;
row.lowT_denominator = lowT_denominator;
row.T_lowCandidate_K = T_lowCandidate_K;
row.highBranchConsistent = highBranchConsistent;
row.lowBranchConsistent = lowBranchConsistent;
row.ambiguousBranch = ambiguousBranch;
row.branchUsed = branchUsed;
row.T_K = T_K;
row.T_C = T_C;

end


function nonFiniteCauses = findNonFiniteCauses(row)
% findNonFiniteCauses
% Identify raw-input and derived-variable conditions that can produce a NaN
% or Inf temperature. These diagnostics do not alter stored values.

maximumCauses = 24;
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

derivedValues = {row.KD, row.XCa_grt, row.XMn_grt, ...
    row.highT_denominator, row.lowT_denominator, ...
    row.T_highCandidate_K, row.T_lowCandidate_K};
derivedNames = {"KD", "XCa_grt", "XMn_grt", ...
    "Eq10a denominator", "Eq10b denominator", ...
    "Eq10a temperature candidate", "Eq10b temperature candidate"};

for i = 1:numel(derivedValues)
    value = derivedValues{i};
    if any(~isfinite(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is non-finite";
    elseif any(value(:) == 0) && ...
            (derivedNames{i} == "KD" || contains(derivedNames{i}, "denominator"))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is zero";
    elseif any(value(:) < 0) && ...
            (derivedNames{i} == "KD" || derivedNames{i} == "XCa_grt" || ...
             derivedNames{i} == "XMn_grt")
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is negative";
    end
end

if any(row.branchUsed == "invalid")
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "no mathematically valid Eq. (10a) or Eq. (10b) result";
end

nonFiniteCauses = nonFiniteCauses(1:nCauses);

end


function textValue = formatNumericValue(value)
% formatNumericValue
% Format a scalar numeric value for a compact diagnostic message.

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
% Use defaultValue only when an optional column is absent. An explicitly
% stored NaN is returned unchanged and is never converted to zero.

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
