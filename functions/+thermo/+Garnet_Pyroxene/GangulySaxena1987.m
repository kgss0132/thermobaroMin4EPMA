function results = GangulySaxena1987(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/GangulySaxena1987.m
% Tested with MATLAB R2024b
%
% Garnet-clinopyroxene Fe2+-Mg exchange thermometer
% Ganguly, J. and Saxena, S.K. (1987)
% Mixtures and Mineral Reactions, Minerals and Rocks, 19, 1-291
% DOI: https://doi.org/10.1007/978-3-642-46601-4
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis (selected by the user from tables) and calculates
% temperature using the Ganguly and Saxena (1987) garnet-clinopyroxene
% Fe2+-Mg exchange thermometer.
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
% The thermometer formulation is summarized in Appendix B.1.4,
% "Garnet-Clinopyroxene", of Ganguly and Saxena (1987), pp. 233-236.
%
% Ganguly and Saxena rejected the formerly suspected low-temperature change
% in slope and retained only the high-temperature Ganguly formulation. The
% experimental basis of that formulation is constrained principally by the
% following Wood dataset (p. 235; see also Ganguly, 1979, p. 1024):
%
%   Temperature : 1100-1400 degreeC
%   Pressure    : 20-45 kbar (2.0-4.5 GPa)
%
% The 1987 book does not define a new independent rectangular pressure
% calibration envelope for Eq. (B.27). This implementation uses 20-45 kbar
% as the practical pressure-warning interval because it is the experimental
% range underlying the retained high-temperature formulation.
%
% Important application notes from Ganguly and Saxena (1987):
%   1) The thermometer describes Fe2+-Mg exchange between coexisting garnet
%      and clinopyroxene. KD must be calculated using Fe2+, not total Fe
%      (p. 233). Fe3_cation_apfu is retained for reference but is not added
%      to Fe_cation_apfu in this implementation.
%   2) The general Eq. (B.28) assumes that mixing in garnet and
%      clinopyroxene can be approximated by the simple-mixture model
%      (p. 235).
%   3) Adequate mixing data for sodic clinopyroxene components were not
%      available. Diopside-jadeite mixing is expected to be approximately
%      ideal above 1000 degreeC but shows significant positive non-ideality
%      near 600 degreeC. Na may therefore have little effect at high
%      temperature but a large effect at lower temperature (p. 235).
%   4) Uncorrected Grt-Cpx Fe-Mg thermometry can be grossly erroneous for
%      low-temperature eclogites and blueschist-facies rocks when sodic Cpx
%      components are important (p. 235). The simplified implementation
%      below omits the explicit jadeite mixing term and should therefore be
%      treated as a Na-poor/high-temperature form of Eqs. (B.27-B.28).
%   5) Ganguly and Bhattacharya's correction to the Ganguly formulation was
%      developed for Cpx with XCa approximately 0.40-0.45 (p. 234). Use
%      outside the demonstrated compositional region requires caution.
%   6) The alternative modified Raheim-Green expression, Eq. (B.29), is
%      based on experiments at 800-1400 degreeC and 30-40 kbar, mostly at
%      30 kbar. It is restricted to XCa_Gt = 0.18 +/- 0.04 and 20-40 mol%
%      jadeite in Cpx (pp. 235-236). Eq. (B.29) is not implemented here and
%      its calibration range must not be treated as the calibration range
%      of the present Eq. (B.27-B.28) implementation.
%   7) Applying Eq. (B.29) to garnet with significantly lower XCa than its
%      experimental range overestimates temperature, and the reverse occurs
%      for higher XCa (p. 236).
%   8) Coexisting garnet and clinopyroxene must represent an equilibrated
%      mineral pair. Analyses from unrelated growth zones or different
%      generations should not be paired without petrographic justification.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside 20-45 kbar,
%   2) a finite calculated temperature is outside 1100-1400 degreeC,
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
% The exchange reaction is defined using Fe2+. Therefore,
% Fe3_cation_apfu is retained for reference but is not added to
% Fe_cation_apfu. If Fe_cation_apfu contains total Fe rather than Fe2+, Fe2+
% should be estimated before this function is used.
%
% Negative finite values in variables used by the thermometer are not
% permitted. Zero values are retained; if they make the equation
% mathematically undefined, the affected result is returned as NaN or Inf
% and a non-stopping warning is printed. NaN values are never replaced by
% zero; they propagate through the calculation and remain in the output.
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
% Ganguly and Saxena (1987), Eqs. (B.27-B.28), pp. 235-236, after
% neglecting the explicit jadeite mixing term and solving for temperature:
%
%          4100 + 11.07*P_kbar + 1510*(XCa_Gt + XMn_Gt)
%   T(K) = ----------------------------------------------------
%                          ln(KD) + 2.400
%
% The previous implementation divided by ln(KD) alone and added Fe3+ to
% Fe2+. Both behaviors are inconsistent with the formulation documented in
% Appendix B.1.4 and are not used here.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = GangulySaxena1987(rawdata_struct, P_kbar)
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
    error('GangulySaxena1987 requires (rawdata_struct, P_kbar).');
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

% Experimental limits underlying the retained high-temperature equation.
calibrationT_min_degC = 1100;
calibrationT_max_degC = 1400;
calibrationP_min_kbar = 20;
calibrationP_max_kbar = 45;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
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
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental range ' ...
             'underlying the retained Ganguly and Saxena (1987) high-' ...
             'temperature formulation: 20-45 kbar (2.0-4.5 GPa). ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar (p. 235; Ganguly, 1979, p. 1024).\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn for finite temperatures outside the experimental interval.
    finiteTemperature = isfinite(row.T_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_C < calibrationT_min_degC | ...
         row.T_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the experimental ' ...
             'range underlying the retained Ganguly and Saxena (1987) ' ...
             'high-temperature formulation: 1100-1400 degreeC. ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g-%.4g degreeC for %s & %s ' ...
             '(p. 235).\n'], ...
            sum(temperatureOutsideCalibration), ...
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
                        'was identified; inspect the stored intermediate variables.\n']);
        else
            fprintf(2, '         Identified cause(s): %s.\n', ...
                char(strjoin(nonFiniteCauses, '; ')));
        end
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'GangulySaxena1987', ...
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
% Verify all columns required by the Ganguly and Saxena (1987) equation.

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
    error(['GangulySaxena1987: required variable(s) are missing: ' ...
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
% Stop when a required thermometer input is negative or infinite. Zero and
% NaN remain available to the calculation and result diagnostics.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

invalidInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nInvalid = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    variableValue = data_grt.(variableName);
    if any(isinf(variableValue(:)) | variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(cpxVariables)
    variableName = cpxVariables{i};
    variableValue = data_cpx.(variableName);
    if any(isinf(variableValue(:)) | variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Cpx." + string(variableName);
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNames(1:nInvalid);
    error(['GangulySaxena1987: thermometer inputs must be finite or NaN ' ...
        'and must be >= 0. Negative or infinite value(s) were found in: ' ...
        char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function row = calcTemp(data_grt, data_cpx, P_kbar)
% calcTemp
% Compute Ganguly and Saxena (1987) Grt-Cpx temperatures for one mineral
% pair and a scalar or vector of pressures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_bar = P_kbar .* 1000;

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
% Fe3+ is also excluded from the divalent-site denominator.
xSiteSum_grt = Fe2_grt + Mg_grt + Mn_grt + Ca_grt;
XCa_grt = Ca_grt ./ xSiteSum_grt;
XMn_grt = Mn_grt ./ xSiteSum_grt;

% --- Ganguly and Saxena (1987), Eqs. (B.27-B.28), pp. 235-236 ---
temperatureNumerator = 4100 + 11.07 .* P_kbar ...
    + 1510 .* (XCa_grt + XMn_grt);
temperatureDenominator = lnKD + 2.400;

% No input or intermediate NaN is replaced. Mathematically undefined
% combinations (for example KD <= 0 or a non-positive denominator) are
% stored as NaN rather than being converted into an apparently finite
% temperature by IEEE arithmetic such as finite/(-Inf) = -0.
calculationDomainValid = isfinite(KD) & KD > 0 ...
    & isfinite(xSiteSum_grt) & xSiteSum_grt > 0 ...
    & isfinite(XCa_grt) & XCa_grt >= 0 ...
    & isfinite(XMn_grt) & XMn_grt >= 0 ...
    & isfinite(temperatureNumerator) ...
    & isfinite(temperatureDenominator) & temperatureDenominator > 0;

T_K = NaN(nP, 1);
T_K(calculationDomainValid) = ...
    temperatureNumerator(calculationDomainValid) ./ ...
    temperatureDenominator(calculationDomainValid);
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
row.xSiteSum_grt = xSiteSum_grt;
row.XCa_grt = XCa_grt;
row.XMn_grt = XMn_grt;
row.KD = KD;
row.lnKD = lnKD;
row.temperatureNumerator = temperatureNumerator;
row.temperatureDenominator = temperatureDenominator;
row.calculationDomainValid = calculationDomainValid;
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

derivedValues = {row.xSiteSum_grt, row.KD, row.lnKD, ...
    row.XCa_grt, row.XMn_grt, row.temperatureNumerator, ...
    row.temperatureDenominator};
derivedNames = {"garnet X-site sum", "KD", "lnKD", ...
    "XCa_grt", "XMn_grt", "temperature numerator", ...
    "temperature denominator"};

for i = 1:numel(derivedValues)
    value = derivedValues{i};
    if any(~isfinite(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is non-finite";
    elseif any(value(:) <= 0) && ...
            (derivedNames{i} == "garnet X-site sum" || ...
             derivedNames{i} == "KD" || ...
             derivedNames{i} == "temperature denominator")
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is zero or negative";
    elseif any(value(:) < 0) && ...
            (derivedNames{i} == "KD" || ...
             derivedNames{i} == "XCa_grt" || ...
             derivedNames{i} == "XMn_grt")
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is negative";
    end
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
