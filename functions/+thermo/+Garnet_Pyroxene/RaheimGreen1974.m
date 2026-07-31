function results = RaheimGreen1974(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/RaheimGreen1974.m
% Tested structurally against MATLAB R2024b syntax
%
% Garnet-clinopyroxene Fe-Mg exchange thermometer
% Raheim, A. and Green, D.H. (1974)
% Contributions to Mineralogy and Petrology, 48, 179-203
% DOI: https://doi.org/10.1007/BF00383355
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis (selected independently from tables) and calculates
% temperature using the Raheim and Green (1974) Grt-Cpx Fe2+-Mg thermometer.
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
% Raheim and Green (1974) calibrated Grt-Cpx Fe2+-Mg partitioning using a
% mineral mix, a glass of typical tholeiite composition, and nine tholeiite
% glasses spanning a wide bulk-rock Mg value. The directly documented
% experimental and calibration conditions are summarised in the abstract
% on p. 179:
%
%   Experimental temperature range : 600-1500 degreeC
%   Temperature calibration        : 600-1400 degreeC at 30 kbar
%   Experimental pressure range    : 20-40 kbar
%   Pressure calibration           : 20-40 kbar at 1100 degreeC
%   Starting compositions          : tholeiitic basalt glasses
%   Experimental bulk Mg value     : 6.2-93, where
%                                    Mg value = 100*Mg/(Mg + Fe2+)
%   Composition-independent KD     : demonstrated for bulk Mg value 6.2-85
%
% The 600-1400 degreeC and 20-40 kbar intervals are therefore used for
% non-stopping range warnings. The wider 600-1500 degreeC interval is the
% range over which experimental mineral analyses of comparable accuracy
% were obtained, but the straight ln(KD)-1/T calibration at 30 kbar is
% described from 600 to 1400 degreeC (abstract, p. 179).
%
% Important application cautions stated or implied by the experimental
% design and scope of the paper:
%
%   1) The temperature dependence was determined at 30 kbar, whereas the
%      pressure dependence was determined at 1100 degreeC between 20 and
%      40 kbar. The whole P-T rectangle was not sampled uniformly. Results
%      near combined temperature and pressure extremes involve greater
%      extrapolation of the separate trends (abstract, p. 179).
%
%   2) Pressure has a larger effect on KD than previously predicted. A
%      pressure estimate is required to obtain a unique temperature from
%      Eq. (1). Uncertain pressure propagates directly into temperature
%      through the 28.35*P_kbar term (abstract, p. 179).
%
%   3) The calibration is intended for natural eclogites of basaltic bulk
%      composition. The experiments were based mainly on tholeiitic basalt
%      glasses, and the demonstrated absence of a bulk-composition effect
%      applies to bulk Mg value 6.2-85. Application to garnet peridotites,
%      pyroxenites, felsic granulites, calc-silicate rocks, extremely
%      Mg-rich bulk compositions, or other grossly different compositions
%      is an extrapolation (abstract, p. 179).
%
%   4) The mineral-mix starting material showed incomplete equilibration
%      and was considered unsuitable by itself. Agreement between its
%      minimum KD and the tholeiite-glass KD within experimental uncertainty
%      was used as a reversal of the partition relationship and to justify
%      the glass starting materials (abstract, p. 179). Natural applications
%      likewise require texturally coexisting, equilibrated Grt-Cpx pairs.
%
%   5) Experiments were performed in iron capsules to maintain iron as Fe2+
%      and eliminate uncertainty caused by variation in Fe3+/Fe2+ among
%      runs. KD is therefore defined with Fe2+, not total iron. This
%      implementation uses Fe_cation_apfu as Fe2+ and does not add
%      Fe3_cation_apfu to KD (abstract, p. 179).
%
%   6) Zoning, exsolution, symplectitisation, and mismatched core-rim or
%      mineral-generation analyses can violate the equilibrium assumption.
%      The selected Grt and Cpx analyses should represent the same textural
%      and metamorphic generation.
%
% The separate Raheim and Green (1974) paper on talc-garnet-kyanite-quartz
% schist (Contributions to Mineralogy and Petrology, 43, 223-231) reports an
% approximately 600 degreeC and 10 kbar natural-rock estimate. Those values
% are not calibration limits of this Grt-Cpx thermometer and are not used
% for the warnings implemented here.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside 20-40 kbar,
%   2) a finite calculated temperature is outside 600-1400 degreeC,
%   3) a thermometer input contains NaN, or
%   4) a calculated temperature is NaN or Inf.
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
% An optional output-only column that is absent defaults to zero. An
% explicitly stored NaN is retained and is never converted to zero.
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
% Raheim and Green (1974), abstract and thermometer expression, p. 179:
%
%   KD = (Fe2+/Mg)_Grt / (Fe2+/Mg)_Cpx
%
%          3686 + 28.35*P_kbar
%   T(K) = ----------------------
%              ln(KD) + 2.33
%
%   T(degreeC) = T(K) - 273.15
%
% Fe3+ is deliberately excluded from KD.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = RaheimGreen1974(rawdata_struct, P_kbar)
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
    error('RaheimGreen1974 requires (rawdata_struct, P_kbar).');
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
calibrationT_max_degC = 1400;
calibrationP_min_kbar = 20;
calibrationP_max_kbar = 40;

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
            ['WARNING: Input pressure is outside the experimental ' ...
             'calibration range of Raheim and Green (1974): 20-40 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; input ' ...
             'range = %.4g-%.4g kbar (abstract, p. 179).\n' ...
             '         Pressure dependence was determined at 1100 ' ...
             'degreeC, whereas temperature dependence was determined at ' ...
             '30 kbar (abstract, p. 179).\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    finiteTemperature = isfinite(row.T_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_C < calibrationT_min_degC | ...
         row.T_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the direct ' ...
             'straight-line calibration interval of Raheim and Green ' ...
             '(1974): 600-1400 degreeC at 30 kbar. %d of %d finite ' ...
             'temperature point(s) are outside the range; calculated ' ...
             'finite range = %.4g-%.4g degreeC for %s & %s ' ...
             '(abstract, p. 179).\n' ...
             '         The broader 600-1500 degreeC interval describes ' ...
             'the experimental analysis range, not the direct ' ...
             'ln(KD)-1/T calibration interval.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
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

    % Preserve non-finite results and display all four equation inputs.
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
        'RaheimGreen1974', ...
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
% Verify the table columns required by the Raheim-Green thermometer.

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
    error(['RaheimGreen1974: required variable(s) are missing: ' ...
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
    error(['RaheimGreen1974: thermometer inputs must be finite or NaN ' ...
        'and must be >= 0. Negative or infinite value(s) were found in: ' ...
        char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function row = calcTemp(data_grt, data_cpx, P_kbar)
% calcTemp
% Compute Raheim-Green temperatures for one Grt-Cpx pair and a scalar or
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

% --- Raheim and Green (1974) ---
% Fe3+ is deliberately excluded from both Fe-Mg exchange ratios.
FeMg_grt = Fe2_grt ./ Mg_grt;
FeMg_cpx = Fe2_cpx ./ Mg_cpx;
KD = FeMg_grt ./ FeMg_cpx;
lnKD = log(KD);
temperatureNumerator = 3686 + 28.35 .* P_kbar;
temperatureDenominator = lnKD + 2.33;

calculationDomainValid = ...
    isfinite(Fe2_grt) & Fe2_grt > 0 ...
    & isfinite(Mg_grt) & Mg_grt > 0 ...
    & isfinite(Fe2_cpx) & Fe2_cpx > 0 ...
    & isfinite(Mg_cpx) & Mg_cpx > 0 ...
    & isfinite(FeMg_grt) & FeMg_grt > 0 ...
    & isfinite(FeMg_cpx) & FeMg_cpx > 0 ...
    & isfinite(KD) & KD > 0 ...
    & isfinite(lnKD) ...
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
row.temperatureNumerator = temperatureNumerator;
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
    row.temperatureNumerator, row.temperatureDenominator};
derivedNames = {"garnet Fe/Mg", "clinopyroxene Fe/Mg", "KD", "lnKD", ...
    "temperature numerator", "temperature denominator"};

for i = 1:numel(derivedValues)
    value = derivedValues{i};
    if any(~isfinite(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is non-finite";
    elseif any(value(:) <= 0) && ...
            (derivedNames{i} == "garnet Fe/Mg" || ...
             derivedNames{i} == "clinopyroxene Fe/Mg" || ...
             derivedNames{i} == "KD" || ...
             derivedNames{i} == "temperature numerator" || ...
             derivedNames{i} == "temperature denominator")
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is zero or negative";
    end
end

if any(~row.calculationDomainValid)
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = ...
        "Raheim and Green (1974) calculation domain is invalid";
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
