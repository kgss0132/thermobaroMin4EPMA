function results = Krogh1988(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/Krogh1988.m
% Tested structurally against MATLAB R2024b syntax
%
% Garnet-clinopyroxene Fe-Mg exchange thermometer
% Krogh, E.J. (1988)
% Contributions to Mineralogy and Petrology, 99, 44-48
% DOI: https://doi.org/10.1007/BF00399364
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis (selected independently from tables) and calculates
% temperature using the Krogh (1988) Grt-Cpx Fe2+-Mg thermometer.
%
% A scalar pressure from startThermoCalc_fixedP or a pressure vector from
% startThermoCalc_rangeP can be supplied. For each selected Grt-Cpx pair,
% the output contains one row for every pressure value.
%
% The function supports repeated mineral-pair selections. Each result is
% buffered as a table block and the blocks are concatenated only once after
% the interactive selection loop has finished.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Krogh (1988) reinterpreted experimental Grt-Cpx Fe2+-Mg partitioning data
% from Raheim and Green (1974a, b), Mori and Green (1978), and Ellis and
% Green (1979). The following limits and qualifications should be considered
% when applying Eq. (14):
%
%   Temperature : 600-1300 degreeC in the compiled calibration dataset
%                 (Table 1, pp. 45-46).
%
%   Pressure    : The calibration is centred on 30 kbar. Most data were at
%                 30 kbar; experiments at 20 and 40 kbar were normalised to
%                 30 kbar using C = -10 (Table 1 and its footnote, pp. 45-46).
%                 The pressure correction in Eq. (14) adopts dV/R = -10
%                 from Powell (1985), rather than independently fitting a
%                 broad pressure dataset (p. 46). Consequently, 20-40 kbar
%                 is used here as the experimental pressure-warning range,
%                 but the calibration should still be regarded as primarily
%                 anchored at 30 kbar.
%
%   Garnet XCa : The quadratic relationship is considered a satisfactory
%                 approximation for XCa_Grt = 0.10-0.50 (p. 47), where
%
%                   XCa_Grt = Ca/(Ca + Mn + Fe2+ + Mg).
%
% Important application cautions stated or demonstrated in the paper:
%
%   1) Experimental data at lower, crustal temperatures are scarce (p. 47).
%      Below 900 degreeC, Krogh (1988) gives progressively lower values than
%      Powell (1985), especially at lower XCa_Grt (abstract, p. 44). This
%      implementation therefore prints a separate low-temperature caution
%      for finite results below 900 degreeC, even when they remain within the
%      formal 600-1300 degreeC experimental-data interval.
%
%   2) The ln(KD)-XCa relationship is poorly constrained in the low-XCa
%      region. This is especially important for garnet peridotites,
%      orthopyroxene eclogites, and other rocks with low-Ca garnet
%      (pp. 47-48). The quadratic approximation should not be extrapolated
%      casually outside XCa_Grt = 0.10-0.50.
%
%   3) KD is defined using ferrous iron, Fe2+, not total iron (Eq. 2, p. 44).
%      Determining Fe2+/Fe3+ in clinopyroxene remains a major problem.
%      Stoichiometric recalculation is sensitive to analytical error,
%      nonstoichiometry, and submicroscopic intergrowth or exsolution; the
%      problem can be severe for low-total-Fe clinopyroxene from garnet
%      ultrabasites (p. 48). This implementation uses Fe_cation_apfu as Fe2+
%      and does not add Fe3_cation_apfu to the exchange ratio or XCa_Grt.
%
%   4) Possible clinopyroxene-composition effects are not incorporated in
%      the calibration. The paper reports little apparent jadeite effect for
%      XJd approximately 0.10-0.60, but warns that XJd > 0.60 and departures
%      of (Na + Ca)_Cpx from approximately 1 may affect Fe-Mg partitioning;
%      the problem may increase with temperature (pp. 47-48).
%
%   5) A possible dependence of KD on garnet Mg number remains unresolved.
%      It was not observed in the natural test suite over Mg/(Mg + Fe) =
%      0.17-0.55, but had been reported by experimental studies (p. 48).
%
%   6) The natural application suite comprised eclogites and omphacite-
%      bearing gneisses with XCa_Grt = 0.15-0.44, garnet Mg number =
%      0.17-0.54, and clinopyroxene Na content = 0.11-0.44; calculated
%      temperatures were 690-790 degreeC (pp. 46-47). These are application
%      examples, not independent experimental calibration limits.
%
%   7) Table 1 residuals extend approximately from -80 to +84 degreeC
%      (pp. 45-46). Krogh (1988) did not present this as a universal formal
%      precision applicable to every rock and composition.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside 20-40 kbar,
%   2) a finite calculated temperature is outside 600-1300 degreeC,
%   3) a finite calculated temperature is below 900 degreeC,
%   4) XCa_Grt is outside 0.10-0.50,
%   5) a thermometer input contains NaN, or
%   6) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Cpx    : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns contain
% normalised cation data.
%
% Required thermometer variables:
%   Garnet:
%     Fe_cation_apfu       % Fe2+ used in KD and XCa_Grt
%     Mg_cation_apfu
%     Ca_cation_apfu
%
%   Clinopyroxene:
%     Fe_cation_apfu       % Fe2+ used in KD
%     Mg_cation_apfu
%
% Optional garnet thermometer variable:
%     Mn_cation_apfu       % used in XCa_Grt; defaults to 0 only if absent
%
% Optional variables retained in the output when present:
%     Fe3_cation_apfu, Ti_cation_apfu, Al_cation_apfu,
%     Si_cation_apfu, Na_cation_apfu, K_cation_apfu
%
% The optional Cpx Ca and Mn variables are also retained in the output but
% are not used in the Krogh (1988) equation. An optional column that is
% absent defaults to zero for output compatibility. An explicitly stored
% NaN is retained and is never converted to zero.
%
% All finite values actually used by the calculation must be >= 0. Negative
% or infinite values stop the calculation with an error. NaN values remain
% NaN, propagate through the calculation, and trigger non-stopping fprintf
% diagnostics. Zero is permitted as an input value, but a zero Fe or Mg
% value makes KD mathematically invalid; the affected temperature is then
% returned as NaN and reported without stopping the calculation.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Krogh (1988), Eqs. (2) and (14), pp. 44 and 46:
%
%   KD = (Fe2+/Mg)_Grt / (Fe2+/Mg)_Cpx
%
%   XCa_Grt = Ca_Grt/(Ca_Grt + Mn_Grt + Fe2+_Grt + Mg_Grt)
%
%          -6173*(XCa_Grt)^2 + 6731*XCa_Grt + 1879 + 10*P_kbar
%   T(K) = -----------------------------------------------------------
%                            ln(KD) + 1.393
%
%   T(degreeC) = T(K) - 273.15
%
% The minus sign on the 6173*(XCa_Grt)^2 term is essential. The original
% unmodified Krogh1988.m used a positive sign and added Fe3+ to Fe2+; both
% behaviours are corrected here to follow the published definition.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Krogh1988(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Cpx tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Cpx pair
%

%% Input validation
if nargin < 2
    error('Krogh1988 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || ~isreal(P_kbar) || isempty(P_kbar) || ...
        ~isvector(P_kbar) || ...
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
% Store each calculation as one table block. The cell buffer is preallocated
% and doubled only when its capacity is exhausted; the full output table is
% concatenated once after the interactive loop.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degC = 600;
calibrationT_max_degC = 1300;
lowTemperatureCaution_degC = 900;
calibrationP_min_kbar = 20;
calibrationP_max_kbar = 40;
calibrationXCa_min = 0.10;
calibrationXCa_max = 0.50;

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

    validateNonNegativeInputs(selectedData_grt, selectedData_cpx);
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

    % Echo the calculated temperature values.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_C(1)) ' to ' num2str(row.T_C(end)) ' degreeC']);
    end

    % The same pressure vector is used for all mineral pairs, so print this
    % warning only once per function call.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental pressure ' ...
             'range represented in Krogh (1988): 20-40 kbar. %d of %d ' ...
             'pressure point(s) are outside the range; input range = ' ...
             '%.4g-%.4g kbar (Table 1 and footnote, pp. 45-46).\n' ...
             '         The calibration is primarily anchored at 30 kbar; ' ...
             'the pressure correction in Eq. (14) was adopted from ' ...
             'Powell (1985) (p. 46).\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn for finite temperatures outside the compiled experimental range.
    finiteTemperature = isfinite(row.T_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_C < calibrationT_min_degC | ...
         row.T_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the experimental ' ...
             'data range used by Krogh (1988): 600-1300 degreeC. %d of %d ' ...
             'finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g-%.4g degreeC for %s & %s ' ...
             '(Table 1, pp. 45-46).\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)));
    end

    % The paper specifically cautions that lower-temperature experiments are
    % sparse, even though the compiled table reaches 600 degreeC.
    lowTemperatureResult = finiteTemperature & ...
        row.T_C < lowTemperatureCaution_degC;
    if any(lowTemperatureResult)
        fprintf(2, ...
            ['CAUTION: %d of %d finite temperature point(s) are below ' ...
             '900 degreeC for %s & %s. Krogh (1988) states that lower, ' ...
             'crustal-temperature experimental data are scarce, and that ' ...
             'differences from Powell (1985) increase toward lower ' ...
             'temperature and lower XCa_Grt (pp. 44, 47).\n'], ...
            sum(lowTemperatureResult), ...
            sum(finiteTemperature), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)));
    end

    % The quadratic XCa correction is recommended only over 0.10-0.50.
    finiteXCa = isfinite(row.XCa_grt);
    xCaOutsideCalibration = finiteXCa & ...
        (row.XCa_grt < calibrationXCa_min | ...
         row.XCa_grt > calibrationXCa_max);
    if any(xCaOutsideCalibration)
        finiteXCaValues = row.XCa_grt(finiteXCa);
        fprintf(2, ...
            ['WARNING: Garnet XCa is outside the compositional interval ' ...
             'for which the quadratic approximation is considered ' ...
             'satisfactory in Krogh (1988): XCa_Grt = 0.10-0.50. ' ...
             'Input XCa_Grt = %.8g for %s & %s (p. 47).\n'], ...
            finiteXCaValues(1), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)));
    end

    % Report explicitly stored NaN values in variables used by the equation.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for ' ...
             '%s & %s: %s.\n' ...
             '         The calculation was continued, and NaN was not ' ...
             'replaced by zero.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Preserve non-finite results and display every raw value used by Eq. 14.
    invalidTemperature = ~isfinite(row.T_C);
    if any(invalidTemperature)
        nonFiniteCauses = findNonFiniteCauses(row);

        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for ' ...
             '%s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), ...
            numel(row.T_C), ...
            sum(isnan(row.T_C)), ...
            sum(isinf(row.T_C)));

        fprintf(2, ...
            ['         Thermometer inputs used: ' ...
             'Garnet.Fe_cation_apfu=%s, Garnet.Mg_cation_apfu=%s, ' ...
             'Garnet.Mn_cation_apfu=%s, Garnet.Ca_cation_apfu=%s, ' ...
             'Cpx.Fe_cation_apfu=%s, Cpx.Mg_cation_apfu=%s.\n'], ...
            formatNumericValue(row.Fe2_grt(1)), ...
            formatNumericValue(row.Mg_grt(1)), ...
            formatNumericValue(row.Mn_grt(1)), ...
            formatNumericValue(row.Ca_grt(1)), ...
            formatNumericValue(row.Fe2_cpx(1)), ...
            formatNumericValue(row.Mg_cpx(1)));

        if isempty(nonFiniteCauses)
            fprintf(2, ...
                ['         No explicit NaN, zero, Inf, or invalid ' ...
                 'intermediate value was identified; inspect the stored ' ...
                 'intermediate variables.\n']);
        else
            fprintf(2, '         Identified cause(s): %s.\n', ...
                char(strjoin(nonFiniteCauses, '; ')));
        end
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Krogh1988', ...
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
% Verify all table columns required by Krogh (1988), Eq. (14).

requiredGrt = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu'};
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
    error(['Krogh1988: required variable(s) are missing: ' ...
        char(strjoin(missingNames, ', ')) '.']);
end

end


function nanInputNames = findNaNInputs(data_grt, data_cpx)
% findNaNInputs
% Return names of equation inputs containing explicitly stored NaN values.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

nanInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nNaN = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    if ismember(variableName, data_grt.Properties.VariableNames)
        variableValue = data_grt.(variableName);
        if any(isnan(variableValue(:)))
            nNaN = nNaN + 1;
            nanInputNames(nNaN) = "Garnet." + string(variableName) + "=NaN";
        end
    end
end

for i = 1:numel(cpxVariables)
    variableName = cpxVariables{i};
    variableValue = data_cpx.(variableName);
    if any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanInputNames(nNaN) = "Cpx." + string(variableName) + "=NaN";
    end
end

nanInputNames = nanInputNames(1:nNaN);

end


function validateNonNegativeInputs(data_grt, data_cpx)
% validateNonNegativeInputs
% Reject negative or infinite values actually used by Eq. (14). NaN remains
% available to propagate through the equation and trigger diagnostics.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

invalidInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nInvalid = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    if ~ismember(variableName, data_grt.Properties.VariableNames)
        continue;
    end
    variableValue = data_grt.(variableName);
    if ~isnumeric(variableValue) || ~isreal(variableValue) || ...
            ~isscalar(variableValue)
        error('Garnet.%s must be a numeric scalar in the selected row.', ...
            variableName);
    end
    if isinf(variableValue) || variableValue < 0
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(cpxVariables)
    variableName = cpxVariables{i};
    variableValue = data_cpx.(variableName);
    if ~isnumeric(variableValue) || ~isreal(variableValue) || ...
            ~isscalar(variableValue)
        error('Cpx.%s must be a numeric scalar in the selected row.', ...
            variableName);
    end
    if isinf(variableValue) || variableValue < 0
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Cpx." + string(variableName);
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNames(1:nInvalid);
    error(['Krogh1988: thermometer inputs must be finite or NaN and ' ...
        'must be >= 0. Negative or infinite value(s) were found in: ' ...
        char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function row = calcTemp(data_grt, data_cpx, P_kbar)
% calcTemp
% Compute Krogh (1988) temperatures for one Grt-Cpx pair and a scalar or
% vector of pressures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;

% --- Extract and expand garnet cations ---
Fe2_grt = repmat(getRequiredValue(data_grt, ...
    'Fe_cation_apfu', 'Garnet'), nP, 1);
Fe3_grt = repmat(getOptionalValue(data_grt, ...
    'Fe3_cation_apfu', 0, 'Garnet'), nP, 1);
Mg_grt = repmat(getRequiredValue(data_grt, ...
    'Mg_cation_apfu', 'Garnet'), nP, 1);
Mn_grt = repmat(getOptionalValue(data_grt, ...
    'Mn_cation_apfu', 0, 'Garnet'), nP, 1);
Ca_grt = repmat(getRequiredValue(data_grt, ...
    'Ca_cation_apfu', 'Garnet'), nP, 1);
Ti_grt = repmat(getOptionalValue(data_grt, ...
    'Ti_cation_apfu', 0, 'Garnet'), nP, 1);
Al_grt = repmat(getOptionalValue(data_grt, ...
    'Al_cation_apfu', 0, 'Garnet'), nP, 1);
Si_grt = repmat(getOptionalValue(data_grt, ...
    'Si_cation_apfu', 0, 'Garnet'), nP, 1);
Na_grt = repmat(getOptionalValue(data_grt, ...
    'Na_cation_apfu', 0, 'Garnet'), nP, 1);
K_grt = repmat(getOptionalValue(data_grt, ...
    'K_cation_apfu', 0, 'Garnet'), nP, 1);

% --- Extract and expand clinopyroxene cations ---
Fe2_cpx = repmat(getRequiredValue(data_cpx, ...
    'Fe_cation_apfu', 'Cpx'), nP, 1);
Fe3_cpx = repmat(getOptionalValue(data_cpx, ...
    'Fe3_cation_apfu', 0, 'Cpx'), nP, 1);
Mg_cpx = repmat(getRequiredValue(data_cpx, ...
    'Mg_cation_apfu', 'Cpx'), nP, 1);
Mn_cpx = repmat(getOptionalValue(data_cpx, ...
    'Mn_cation_apfu', 0, 'Cpx'), nP, 1);
Ca_cpx = repmat(getOptionalValue(data_cpx, ...
    'Ca_cation_apfu', 0, 'Cpx'), nP, 1);
Ti_cpx = repmat(getOptionalValue(data_cpx, ...
    'Ti_cation_apfu', 0, 'Cpx'), nP, 1);
Al_cpx = repmat(getOptionalValue(data_cpx, ...
    'Al_cation_apfu', 0, 'Cpx'), nP, 1);
Si_cpx = repmat(getOptionalValue(data_cpx, ...
    'Si_cation_apfu', 0, 'Cpx'), nP, 1);
Na_cpx = repmat(getOptionalValue(data_cpx, ...
    'Na_cation_apfu', 0, 'Cpx'), nP, 1);
K_cpx = repmat(getOptionalValue(data_cpx, ...
    'K_cation_apfu', 0, 'Cpx'), nP, 1);

% --- Krogh (1988), Eqs. (2) and (14) ---
% Fe3+ is deliberately excluded from KD and the garnet X-site denominator.
sum_grt_FeMgMnCa = Fe2_grt + Mg_grt + Mn_grt + Ca_grt;
FeMg_grt = Fe2_grt ./ Mg_grt;
FeMg_cpx = Fe2_cpx ./ Mg_cpx;
KD = FeMg_grt ./ FeMg_cpx;
lnKD = log(KD);
XCa_grt = Ca_grt ./ sum_grt_FeMgMnCa;

temperatureNumerator = ...
    -6173 .* (XCa_grt .^ 2) + 6731 .* XCa_grt + 1879 ...
    + 10 .* P_kbar;
temperatureDenominator = lnKD + 1.393;

% Preserve NaN and return NaN for invalid mathematical/physical domains.
calculationDomainValid = ...
    isfinite(Fe2_grt) & Fe2_grt > 0 ...
    & isfinite(Mg_grt) & Mg_grt > 0 ...
    & isfinite(Mn_grt) & Mn_grt >= 0 ...
    & isfinite(Ca_grt) & Ca_grt >= 0 ...
    & isfinite(Fe2_cpx) & Fe2_cpx > 0 ...
    & isfinite(Mg_cpx) & Mg_cpx > 0 ...
    & isfinite(sum_grt_FeMgMnCa) & sum_grt_FeMgMnCa > 0 ...
    & isfinite(KD) & KD > 0 ...
    & isfinite(lnKD) ...
    & isfinite(XCa_grt) & XCa_grt >= 0 & XCa_grt <= 1 ...
    & isfinite(temperatureNumerator) ...
    & isfinite(temperatureDenominator) & temperatureDenominator > 0;

T_K = NaN(nP, 1);
T_K(calculationDomainValid) = ...
    temperatureNumerator(calculationDomainValid) ./ ...
    temperatureDenominator(calculationDomainValid);
T_C = T_K - 273.15;

% --- Pack outputs ---
% FeUsed is retained for compatibility with the original output schema, but
% correctly equals Fe2+ rather than Fe2+ + Fe3+.
row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.FeUsed_grt = Fe2_grt;
row.Mg_grt = Mg_grt;
row.Mn_grt = Mn_grt;
row.Ca_grt = Ca_grt;
row.Ti_grt = Ti_grt;
row.Al_grt = Al_grt;
row.Si_grt = Si_grt;
row.Na_grt = Na_grt;
row.K_grt = K_grt;

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

row.sum_grt_FeMgMnCa = sum_grt_FeMgMnCa;
row.FeMg_grt = FeMg_grt;
row.FeMg_cpx = FeMg_cpx;
row.XCa_grt = XCa_grt;
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
% Identify raw-input and derived-variable conditions that can produce NaN or
% Inf without altering the stored output values.

maximumCauses = 24;
nonFiniteCauses = strings(maximumCauses, 1);
nCauses = 0;

inputValues = {row.Fe2_grt, row.Mg_grt, row.Mn_grt, row.Ca_grt, ...
    row.Fe2_cpx, row.Mg_cpx};
inputNames = {"Garnet.Fe_cation_apfu", "Garnet.Mg_cation_apfu", ...
    "Garnet.Mn_cation_apfu", "Garnet.Ca_cation_apfu", ...
    "Cpx.Fe_cation_apfu", "Cpx.Mg_cation_apfu"};
zeroIsInvalid = [true, true, false, false, true, true];

for i = 1:numel(inputValues)
    value = inputValues{i};
    if any(isnan(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is NaN";
    elseif any(isinf(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is Inf";
    elseif zeroIsInvalid(i) && any(value(:) == 0)
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is zero";
    end
end

derivedValues = {row.sum_grt_FeMgMnCa, row.FeMg_grt, row.FeMg_cpx, ...
    row.KD, row.lnKD, row.XCa_grt, row.temperatureNumerator, ...
    row.temperatureDenominator};
derivedNames = {"garnet Fe-Mg-Mn-Ca sum", "garnet Fe/Mg", ...
    "clinopyroxene Fe/Mg", "KD", "lnKD", "XCa_grt", ...
    "temperature numerator", "temperature denominator"};

for i = 1:numel(derivedValues)
    value = derivedValues{i};
    if any(~isfinite(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is non-finite";
    elseif any(value(:) <= 0) && ...
            (derivedNames{i} == "garnet Fe-Mg-Mn-Ca sum" || ...
             derivedNames{i} == "garnet Fe/Mg" || ...
             derivedNames{i} == "clinopyroxene Fe/Mg" || ...
             derivedNames{i} == "KD" || ...
             derivedNames{i} == "temperature denominator")
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is zero or negative";
    elseif derivedNames{i} == "XCa_grt" && ...
            any(value(:) < 0 | value(:) > 1)
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = "XCa_grt is outside 0-1";
    end
end

if any(~row.calculationDomainValid)
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "Krogh (1988) calculation domain is invalid";
end

nonFiniteCauses = nonFiniteCauses(1:nCauses);

end


function textValue = formatNumericValue(value)
% formatNumericValue
% Format one scalar value for compact fprintf diagnostics.

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
% Extract a required scalar numeric value; an explicit NaN is unchanged.

if ~ismember(variableName, data_tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = data_tbl.(variableName);
if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end


function value = getOptionalValue(data_tbl, variableName, defaultValue, mineralLabel)
% getOptionalValue
% Use defaultValue only when a column is absent. Explicit NaN is unchanged.

if ismember(variableName, data_tbl.Properties.VariableNames)
    value = data_tbl.(variableName);
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
        error('%s.%s must be a numeric scalar in the selected row.', ...
            mineralLabel, variableName);
    end
else
    value = defaultValue;
end

end
