function results = MysenHeier1972(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/MysenHeier1972.m
% Tested structurally against MATLAB R2024b syntax
%
% Garnet-clinopyroxene Fe-Mg exchange thermometer
% Mysen, B.O. and Heier, K.S. (1972)
% Contributions to Mineralogy and Petrology, 36, 73-94
% DOI: https://doi.org/10.1007/BF00372836
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis (selected independently from tables) and calculates
% temperature using the Mysen and Heier (1972) Grt-Cpx Fe2+-Mg thermometer.
%
% A scalar pressure from startThermoCalc_fixedP or a pressure vector from
% startThermoCalc_rangeP can be supplied. For each selected Grt-Cpx pair,
% the output contains one row for every pressure value. The thermometer is
% pressure-independent, so the calculated temperature is repeated for all
% pressure rows while P_kbar is retained for downstream processing.
%
% The function supports repeated mineral-pair selections. Each result is
% buffered as a table block and the blocks are concatenated only once after
% the interactive selection loop has finished.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Mysen and Heier (1972) calibrated the Grt-Cpx Fe2+-Mg exchange relation
% from published natural garnet-clinopyroxene pairs with independently
% estimated temperatures. It is an empirical natural-rock calibration, not
% a new reversed experimental calibration. The following limits and
% qualifications should be considered when applying the thermometer:
%
%   Temperature : 600-1100 degreeC. The authors define the calibration as
%                 a straight line over this interval (pp. 88-89; Fig. 9).
%
%   Pressure    : 6-36 kbar. The rocks used to establish the relationship
%                 formed within this load-pressure interval (pp. 86-89).
%                 Pressure is not present in the thermometer equation
%                 because its effect on KD was assumed negligible (p. 86).
%                 This range is therefore the pressure domain represented
%                 by the calibration data, not a fitted pressure correction.
%
%   Assemblage  : Coexisting, equilibrated garnet and clinopyroxene from
%                 eclogites and related high-pressure rocks (pp. 86-89).
%
% Important application cautions stated or demonstrated in the paper:
%
%   1) Eclogite garnets formed below approximately 500 degreeC are often
%      chemically zoned, so element distribution between garnet and
%      clinopyroxene must be treated cautiously. Low-temperature deviation
%      from ideality may also be important (pp. 88-89). This implementation
%      prints a separate caution for finite results below 500 degreeC.
%
%   2) KD is defined using ferrous iron, Fe2+, not total iron (Eq. iii,
%      p. 86). Determination of the Cpx Fe3+/Fe2+ ratio is important
%      (p. 77), and wet-chemical Fe3+ may be overestimated (p. 89). This
%      implementation uses Fe_cation_apfu as Fe2+ and does not add
%      Fe3_cation_apfu to the exchange ratio.
%
%   3) The selected Grt-Cpx analyses must represent an equilibrium pair.
%      Garnet zoning, mismatched mineral generations, and retrograde
%      re-equilibration can invalidate KD. Omphacite in the Hareidland
%      eclogite locally breaks down to pyroxene-plagioclase symplectite, and
%      equilibrium was not attained during symplectitisation (p. 90).
%      Secondary symplectitic Cpx should therefore not be paired with peak
%      garnet to estimate peak temperature.
%
%   4) Variable composition had a negligible effect within the Hareidland
%      dataset (p. 86), but the paper does not establish universal XCa, XJd,
%      Mn, or other compositional limits. Extrapolation to markedly different
%      mineral compositions remains uncertain.
%
%   5) The 625 +/- 30 degreeC result reported for the Hareidland eclogite
%      is a result for that sample suite, not a universal thermometer
%      precision (p. 89).
%
%   6) The approximately 14 kbar estimate for the Hareidland eclogite was
%      obtained separately from phase-equilibrium constraints and is not a
%      fixed pressure or pressure calibration of this thermometer (p. 89).
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside 6-36 kbar,
%   2) a finite calculated temperature is outside 600-1100 degreeC,
%   3) a finite calculated temperature is below 500 degreeC,
%   4) a thermometer input contains NaN, or
%   5) a calculated temperature is NaN or Inf.
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
% Required thermometer variables in both mineral tables:
%     Fe_cation_apfu       % Fe2+ used in KD
%     Mg_cation_apfu
%
% Optional variables retained in the output when present:
%     Fe3_cation_apfu, Mn_cation_apfu, Ca_cation_apfu,
%     Ti_cation_apfu, Al_cation_apfu, Si_cation_apfu,
%     Na_cation_apfu, K_cation_apfu
%
% An optional column that is absent defaults to zero for output
% compatibility. An explicitly stored NaN is retained and is never
% converted to zero.
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
% Mysen and Heier (1972), Eq. (iii) and Fig. 9, pp. 86-89:
%
%   KD = (Fe2+/Mg)_Grt / (Fe2+/Mg)_Cpx
%
%   T(K) = 2475 / (ln(KD) + 0.781)
%
%   T(degreeC) = T(K) - 273.15
%
% Fe3+ is deliberately excluded from KD. Pressure is stored in the output
% but is not used in the published temperature equation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = MysenHeier1972(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Cpx tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Cpx pair

%% Input validation
if nargin < 2
    error('MysenHeier1972 requires (rawdata_struct, P_kbar).');
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
calibrationT_max_degC = 1100;
lowTemperatureCaution_degC = 500;
calibrationP_min_kbar = 6;
calibrationP_max_kbar = 36;

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
        'PromptString', ...
        'Please select the Clinopyroxene data you would like to use:', ...
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
        disp([char(string(selectedCode_grt)) ' & ' ...
            char(string(selectedCode_cpx)) ': ' ...
            num2str(row.T_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' ...
            char(string(selectedCode_cpx)) ': ' ...
            num2str(row.T_C(1)) ' to ' ...
            num2str(row.T_C(end)) ' degreeC']);
    end

    % The same pressure vector is used for all mineral pairs, so print this
    % warning only once per function call.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the empirical calibration ' ...
             'data range of Mysen and Heier (1972): 6-36 kbar. %d of %d ' ...
             'pressure point(s) are outside the range; input range = ' ...
             '%.4g-%.4g kbar (pp. 86-89).\n' ...
             '         Pressure is retained in the output but is not used ' ...
             'in this pressure-independent thermometer (p. 86).\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn for finite temperatures outside the interval over which the
    % authors define the straight-line calibration.
    finiteTemperature = isfinite(row.T_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_C < calibrationT_min_degC | ...
         row.T_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the empirical ' ...
             'calibration range of Mysen and Heier (1972): ' ...
             '600-1100 degreeC. %d of %d finite temperature point(s) are ' ...
             'outside the range; calculated finite range = ' ...
             '%.4g-%.4g degreeC for %s & %s (pp. 88-89; Fig. 9).\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)));
    end

    % Print the paper's specific low-temperature zoning caution.
    lowTemperatureResult = finiteTemperature & ...
        row.T_C < lowTemperatureCaution_degC;
    if any(lowTemperatureResult)
        fprintf(2, ...
            ['CAUTION: %d of %d finite temperature point(s) are below ' ...
             '500 degreeC for %s & %s. Mysen and Heier (1972) state that ' ...
             'eclogite garnets formed below approximately 500 degreeC are ' ...
             'often chemically zoned; low-temperature non-ideality may ' ...
             'also be important (pp. 88-89).\n'], ...
            sum(lowTemperatureResult), ...
            sum(finiteTemperature), ...
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

    % Preserve non-finite results and display every raw value used by the
    % Mysen and Heier (1972) equation.
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
             'Cpx.Fe_cation_apfu=%s, Cpx.Mg_cation_apfu=%s.\n'], ...
            formatNumericValue(row.Fe2_grt(1)), ...
            formatNumericValue(row.Mg_grt(1)), ...
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
        'MysenHeier1972', ...
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
% Verify the table columns required by the Mysen-Heier thermometer.

requiredGrt = {'Fe_cation_apfu', 'Mg_cation_apfu'};
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
    error(['MysenHeier1972: required variable(s) are missing: ' ...
        char(strjoin(missingNames, ', ')) '.']);
end

end


function nanInputNames = findNaNInputs(data_grt, data_cpx)
% findNaNInputs
% Return names and values of equation inputs containing explicitly stored NaN.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

nanInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nNaN = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    variableValue = data_grt.(variableName);
    if any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanInputNames(nNaN) = "Garnet." + string(variableName) + "=NaN";
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
% Reject negative or infinite values used by the equation. NaN is retained.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

invalidInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nInvalid = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
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
    error(['MysenHeier1972: thermometer inputs must be finite or NaN ' ...
        'and must be >= 0. Negative or infinite value(s) were found in: ' ...
        char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function row = calcTemp(data_grt, data_cpx, P_kbar)
% calcTemp
% Compute Mysen-Heier temperatures for one Grt-Cpx pair and a scalar or
% vector of pressures. Pressure is retained but not used in the equation.

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
Ca_grt = repmat(getOptionalValue(data_grt, ...
    'Ca_cation_apfu', 0, 'Garnet'), nP, 1);
Al_grt = repmat(getOptionalValue(data_grt, ...
    'Al_cation_apfu', 0, 'Garnet'), nP, 1);
Si_grt = repmat(getOptionalValue(data_grt, ...
    'Si_cation_apfu', 0, 'Garnet'), nP, 1);

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

% --- Mysen and Heier (1972), Eq. (iii) and Fig. 9 ---
% Fe3+ is deliberately excluded from both Fe-Mg exchange ratios.
FeMg_grt = Fe2_grt ./ Mg_grt;
FeMg_cpx = Fe2_cpx ./ Mg_cpx;
KD = FeMg_grt ./ FeMg_cpx;
lnKD = log(KD);
temperatureDenominator = lnKD + 0.781;

% Preserve NaN and return NaN for invalid mathematical/physical domains.
calculationDomainValid = ...
    isfinite(Fe2_grt) & Fe2_grt > 0 ...
    & isfinite(Mg_grt) & Mg_grt > 0 ...
    & isfinite(Fe2_cpx) & Fe2_cpx > 0 ...
    & isfinite(Mg_cpx) & Mg_cpx > 0 ...
    & isfinite(FeMg_grt) & FeMg_grt > 0 ...
    & isfinite(FeMg_cpx) & FeMg_cpx > 0 ...
    & isfinite(KD) & KD > 0 ...
    & isfinite(lnKD) ...
    & isfinite(temperatureDenominator) & temperatureDenominator > 0;

T_K = NaN(nP, 1);
T_K(calculationDomainValid) = ...
    2475 ./ temperatureDenominator(calculationDomainValid);
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
row.KD = KD;
row.lnKD = lnKD;
row.temperatureDenominator = temperatureDenominator;
row.calculationDomainValid = calculationDomainValid;
row.T_K = T_K;
row.T_C = T_C;

end


function nonFiniteCauses = findNonFiniteCauses(row)
% findNonFiniteCauses
% Identify input and derived-variable conditions that can produce NaN or Inf.

maximumCauses = 16;
nonFiniteCauses = strings(maximumCauses, 1);
nCauses = 0;

inputValues = {row.Fe2_grt, row.Mg_grt, row.Fe2_cpx, row.Mg_cpx};
inputNames = {"Garnet.Fe_cation_apfu", "Garnet.Mg_cation_apfu", ...
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

derivedValues = {row.FeMg_grt, row.FeMg_cpx, row.KD, row.lnKD, ...
    row.temperatureDenominator};
derivedNames = {"garnet Fe/Mg", "clinopyroxene Fe/Mg", "KD", "lnKD", ...
    "temperature denominator"};

for i = 1:numel(derivedValues)
    value = derivedValues{i};
    if any(~isfinite(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is non-finite";
    elseif any(value(:) <= 0) && ...
            (derivedNames{i} == "garnet Fe/Mg" || ...
             derivedNames{i} == "clinopyroxene Fe/Mg" || ...
             derivedNames{i} == "KD" || ...
             derivedNames{i} == "temperature denominator")
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is zero or negative";
    end
end

if any(~row.calculationDomainValid)
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = ...
        "Mysen and Heier (1972) calculation domain is invalid";
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
