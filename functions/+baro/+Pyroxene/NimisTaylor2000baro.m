function results = NimisTaylor2000baro(rawdata_struct, T_degreeC)
% functions/+baro/+Pyroxene/NimisTaylor2000baro.m
% Tested with MATLAB R2024b
%
% Single-clinopyroxene Cr-in-Cpx geobarometer
% Nimis, P. and Taylor, W.R. (2000)
% Contributions to Mineralogy and Petrology, 139, 541-554
% DOI: https://doi.org/10.1007/s004100000156
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Clinopyroxene analysis and
% calculates pressure using the Cr-in-Cpx barometer of Nimis and Taylor
% (2000), Equations (9) and (10).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Cpx analysis, one output row is
% returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Nimis and Taylor (2000) calibrated the Cr-in-Cpx barometer primarily for
% Cr-diopside that equilibrated with garnet in garnet-bearing lherzolitic
% systems. The pressure calibration dataset spans approximately:
%
%   Pressure    : 20-60 kbar
%   Temperature : 900-1400 degreeC
%   Composition : fertile pyrolite through depleted, high-Cr lherzolite
%
% These ranges and the adopted regression are described on p. 544 and in
% Equations (9)-(10). The broader 0-60 kbar and 850-1500 degreeC ranges in
% the abstract refer to the combined barometer-thermometer experimental
% database and should not be used as the direct Cr-in-Cpx pressure-
% calibration range.
%
% Calibration statistics are summarized in Table 2 on p. 547. For all 76
% pressure experiments, the standard deviation of calculated versus
% experimental pressure is approximately 2.3 kbar. The estimated 1-sigma
% uncertainty is approximately 2 kbar at P <= 40 kbar and slightly greater
% than 3 kbar at P > 40 kbar because relatively few calibration experiments
% extend above 40 kbar (p. 547).
%
% IMPORTANT APPLICATION NOTES:
%
% 1) This is a single-mineral barometer, but it is not universally valid for
%    every clinopyroxene. The analyzed Cr-diopside must have equilibrated
%    with garnet under garnet-peridotite conditions. Application to Cpx from
%    spinel peridotites, low-pressure metasomatic assemblages, pyroxenites,
%    magmatic phenocrysts, or grains of uncertain origin may yield
%    geologically meaningless apparent pressures (pp. 548 and 550-551).
%
% 2) The pressure-sensitive activity parameter must satisfy:
%
%      aCaCrTs_cpx >= 0.003
%
%    Nimis and Taylor (2000) treated 0.003 as a safe lower limit and
%    discarded analyses below it because analytical uncertainty becomes
%    strongly amplified at very low Cr contents (p. 548). Values <= 0 are
%    outside the logarithm domain and produce NaN in this implementation.
%
% 3) The natural-sample test was restricted to Cpx with Cr2O3 < 5 wt%.
%    Very Cr-rich, very Al-poor, kosmochlor-rich compositions may not be
%    represented adequately by the simplified CaCrTs activity model
%    (pp. 543 and 548).
%
% 4) Analytical quality and pyroxene stoichiometry must be checked. Nimis
%    and Taylor (2000) accepted analyses having both the tetrahedral-site
%    cation sum and the M1+M2-site cation sum greater than 1.990 on a
%    six-oxygen basis (p. 548). Analyses failing either criterion were
%    rejected in their natural-sample evaluation.
%
% 5) Cr, Al, and Na analytical precision is especially important. The
%    pressure expression contains logarithmic terms and becomes very
%    sensitive to small analytical errors when Cr and aCaCrTs are low
%    (pp. 548-549). Existing NaN values must not be interpreted as zero.
%
% 6) Equation (10) is printed using Na, but for natural samples the authors
%    added measured K to Na when K was reported (p. 548). This implementation
%    therefore uses (Na + K), as in the original script supplied for this
%    project.
%
% 7) Fe3+ is not explicitly corrected in the adopted barometer. If Fe3+
%    behaves similarly to Cr3+, ignoring substantial Fe3+ may overestimate
%    pressure. A test using independently measured Fe3+ lowered pressures
%    by an average of approximately 2.6 kbar and by less than 6 kbar for the
%    examined samples; however, the authors did not recommend a general Fe3+
%    correction (p. 549).
%
% 8) Pressure depends on the externally supplied temperature. The reported
%    temperature sensitivity is approximately 1.2-2.4 kbar per 50 degreeC,
%    depending on composition, with an illustrative value of approximately
%    1.5 kbar per 50 degreeC discussed on p. 549. Temperature and Cpx must
%    therefore represent the same equilibrium stage.
%
% 9) Isolated Cpx grains, diamond inclusions, and polymineralic inclusions
%    may preserve different stages of equilibration. Textural context and
%    coexisting phases must be used to determine whether the calculated
%    pressure represents mantle equilibration, entrapment, or later
%    re-equilibration (p. 548).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 900-1400 degreeC,
%   2) finite calculated pressure is outside 20-60 kbar,
%   3) finite calculated pressure is above 40 kbar, where calibration
%      uncertainty is larger,
%   4) a calculation input contains NaN,
%   5) aCaCrTs_cpx is below 0.003 or outside its mathematical domain,
%   6) the six-oxygen site-sum criteria are not satisfied,
%   7) Cr2O3 is >= 5 wt% when an oxide concentration is available, or
%   8) calculated pressure is NaN, Inf, or negative.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The Cpx table must contain the
% following six-oxygen cation variables:
%
%   Si_cation_apfu
%   Al_cation_apfu        % required by Cr# and site-sum checks
%   Fe_cation_apfu        % total Fe
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu        % required by aCaCrTs
%   K_cation_apfu         % added to Na following the natural-sample use
%   Cr_cation_apfu        % required by Cr# and aCaCrTs
%
% Optional variables retained for output and diagnostic calculations:
%   Mn_cation_apfu
%   Ti_cation_apfu
%   Fe3_cation_apfu
%
% Missing optional variables are represented by NaN, never by zero. A
% present NaN value is retained, propagated through calculations that use
% it, and reported by fprintf. Inf and finite negative values are rejected.
% Zero is retained; if it makes a ratio or logarithm undefined, the
% resulting pressure remains NaN and is reported without stopping.
%
% If a Cr2O3 or Cr2O3_value oxide column is present, it is used only for the
% <5 wt% applicability warning and is not used in the pressure equation.
%
% No liquid composition is used by this barometer. Therefore, treatment of
% Liq F and Cl and cationTotal_liq is not applicable to this function.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION
%
% Equations (9) and (10) of Nimis and Taylor (2000; p. 544):
%
%   P(kbar) =
%     -(T(K)/126.9) * ln(aCaCrTs_cpx)
%     + 15.483 * ln(CrNumber_cpx/T(K))
%     + T(K)/71.38
%     + 107.8
%
% where:
%
%   CrNumber_cpx = Cr / (Cr + Al)
%
%   aCaCrTs_cpx = Cr - 0.81*CrNumber_cpx*(Na + K)
%
% All cations are expressed per formula unit on a six-oxygen basis.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = NimisTaylor2000baro(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing a Cpx table
%   T_degreeC      : non-negative numeric scalar or vector in degreeC.
%                    NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Cpx analysis.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same implementation.
if nargin < 2
    error('NimisTaylor2000baro requires (rawdata_struct, T_degreeC).');
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
if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

T_degreeC = T_degreeC(:);

%% 1) Retrieve cation dataset
% Extract the Cpx table from the input struct. The source table is not
% modified; one selected row is read during each calculation.
disp('=== Step 1: Preparing cation dataset ===');

dataset_cpx = rawdata_struct.Cpx;

if isempty(dataset_cpx)
    error('rawdata_struct.Cpx must be a non-empty table.');
end

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-Cpx result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max(16, height(dataset_cpx));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degreeC = 900;
calibrationT_max_degreeC = 1400;
calibrationP_min_kbar = 20;
calibrationP_max_kbar = 60;
highPressureThreshold_kbar = 40;
minimumSafeActivity = 0.003;
minimumSiteSum = 1.990;
maximumCr2O3_wtpercent = 5;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);
temperatureWarningIssued = false;
originCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a completed calculation.
disp('=== Step 3: Selecting a data code from the list (Cpx) ===');

while true
    % ----- Clinopyroxene selection -----
    % The first table column is used only as the displayed data identifier.
    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Cpx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIdx_cpx);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);
    disp(['Cpx selected: ' char(string(selectedCode_cpx))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the pressure ===');

    % List NaN values without replacing them by zero.
    nanInputNames = findNaNInputs(selectedData_cpx, T_degreeC);

    % Reject Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(selectedData_cpx);

    row = calcPressure(selectedData_cpx, T_degreeC);

    % Repeat the selected identifier for every temperature row.
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_cpx'}, 'Before', 1);

    % Store one result block per selected Cpx. The buffer is expanded only
    % when its current capacity has been exhausted, not every iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo calculated pressure or pressure range.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_cpx)) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_cpx)) ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Print the mineralogical applicability condition once per function call.
    if ~originCautionIssued
        fprintf(2, ...
            ['CAUTION: Nimis and Taylor (2000) calibrated this single-Cpx ' ...
             'barometer for Cr-diopside that equilibrated with garnet in ' ...
             'garnet-peridotite systems. Application to Cpx from spinel ' ...
             'peridotites, metasomatic assemblages, pyroxenites, magmatic ' ...
             'phenocrysts, or grains of uncertain origin may yield ' ...
             'geologically meaningless apparent pressures (pp. 548, 550-551).\n']);
        originCautionIssued = true;
    end

    % Input temperature is common to all selected Cpx analyses, so this
    % warning is printed only once after the first completed calculation.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperatureValues = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the Nimis and Taylor ' ...
             '(2000) Cr-in-Cpx pressure-calibration range of 900-1400 ' ...
             'degreeC (p. 544). %d of %d finite temperature point(s) are ' ...
             'outside the range; input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperatureValues), ...
            min(finiteTemperatureValues), ...
            max(finiteTemperatureValues));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures are outside the direct
    % experimental calibration range.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the Nimis and Taylor ' ...
             '(2000) Cr-in-Cpx calibration range of 20-60 kbar (p. 544). ' ...
             '%d of %d finite pressure point(s) are outside the range; ' ...
             'calculated finite range = %.4g-%.4g kbar for %s.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_cpx)));
    end

    % Above 40 kbar, the paper reports a larger uncertainty because few
    % calibration experiments extend to this range.
    pressureAbove40 = finitePressure & row.P_kbar > highPressureThreshold_kbar;
    if any(pressureAbove40)
        fprintf(2, ...
            ['CAUTION: %d of %d finite pressure point(s) for %s exceed ' ...
             '40 kbar. Nimis and Taylor (2000) estimate approximately ' ...
             '2 kbar 1-sigma uncertainty at P <= 40 kbar and slightly ' ...
             'greater than 3 kbar at P > 40 kbar (p. 547).\n'], ...
            sum(pressureAbove40), ...
            sum(finitePressure), ...
            char(string(selectedCode_cpx)));
    end

    % Report exact NaN input names and vector-temperature indices.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Nimis-Taylor calculation ' ...
             'input(s) for %s: %s.\n' ...
             '         NaN values were retained and propagated; calculated ' ...
             'pressure may remain NaN.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Safe activity limit and mathematical-domain checks.
    activityScalar = row.aCaCrTs_cpx(1);
    if isfinite(activityScalar) && activityScalar > 0 && ...
            activityScalar < minimumSafeActivity
        fprintf(2, ...
            ['WARNING: aCaCrTs_cpx = %.6g is below the safe lower limit ' ...
             'of 0.003 adopted by Nimis and Taylor (2000; p. 548) for %s. ' ...
             'Analytical errors may be strongly amplified; the calculated ' ...
             'pressure has been retained.\n'], ...
            activityScalar, char(string(selectedCode_cpx)));
    end

    invalidEquationTerms = findInvalidEquationTerms(row);
    if ~isempty(invalidEquationTerms)
        fprintf(2, ...
            ['WARNING: Nimis-Taylor equation term(s) are outside their ' ...
             'mathematical domain for %s: %s.\n' ...
             '         The affected pressure values were retained as NaN/Inf ' ...
             'where applicable.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(invalidEquationTerms, ', ')));
    end

    % Six-oxygen crystal-chemical quality criteria from p. 548.
    TsiteScalar = row.T_site_sum_cpx(1);
    M1M2Scalar = row.M1M2_site_sum_cpx(1);

    if isfinite(TsiteScalar) && TsiteScalar <= minimumSiteSum
        fprintf(2, ...
            ['WARNING: Cpx tetrahedral-site cation sum = %.6g does not ' ...
             'satisfy the >1.990 analytical-quality criterion of Nimis ' ...
             'and Taylor (2000; p. 548) for %s.\n'], ...
            TsiteScalar, char(string(selectedCode_cpx)));
    end

    if isfinite(M1M2Scalar) && M1M2Scalar <= minimumSiteSum
        fprintf(2, ...
            ['WARNING: Cpx M1+M2-site cation sum = %.6g does not satisfy ' ...
             'the >1.990 analytical-quality criterion of Nimis and Taylor ' ...
             '(2000; p. 548) for %s.\n'], ...
            M1M2Scalar, char(string(selectedCode_cpx)));
    end

    % Cr2O3 applicability screening is possible only when an oxide column
    % is retained in the selected Cpx table.
    Cr2O3Scalar = row.Cr2O3_wtpercent(1);
    if isfinite(Cr2O3Scalar) && Cr2O3Scalar >= maximumCr2O3_wtpercent
        fprintf(2, ...
            ['WARNING: Cpx Cr2O3 = %.6g wt%% is outside the <5 wt%% ' ...
             'natural-sample screening range used by Nimis and Taylor ' ...
             '(2000; p. 548) for %s. The simplified CaCrTs activity model ' ...
             'may be inaccurate.\n'], ...
            Cr2O3Scalar, char(string(selectedCode_cpx)));
    end

    % Retain and report all non-finite calculated pressures.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_cpx)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressure is retained for diagnostic purposes.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s ' ...
             '(%d of %d points). Values were retained for diagnostic ' ...
             'purposes.\n'], ...
            char(string(selectedCode_cpx)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'NimisTaylor2000baro', ...
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
function nanInputNames = findNaNInputs(data_cpx, T_degreeC)
% findNaNInputs
% Return names of selected input variables containing NaN. Existing NaN
% values are not changed and do not stop the calculation.

cationVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Fe3_cation_apfu'};

maxNames = 1 + numel(cationVariables);
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

for i = 1:numel(cationVariables)
    variableName = cationVariables{i};
    if ismember(variableName, data_cpx.Properties.VariableNames)
        variableValue = data_cpx.(variableName);
        if any(isnan(variableValue(:)))
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = "Cpx." + string(variableName);
        end
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_cpx)
% validateNonNegativeInputs
% Reject Inf and finite negative values in all Cpx cation variables used by
% the equation or its diagnostics. Zero and NaN are intentionally allowed.

cationVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Fe3_cation_apfu'};

invalidInputBuffer = strings(numel(cationVariables), 1);
nInvalidInputs = 0;

for i = 1:numel(cationVariables)
    variableName = cationVariables{i};
    if ismember(variableName, data_cpx.Properties.VariableNames)
        variableValue = data_cpx.(variableName);
        if any(isinf(variableValue(:)) | ...
                (isfinite(variableValue(:)) & variableValue(:) < 0))
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = ...
                "Cpx." + string(variableName);
        end
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['NimisTaylor2000baro: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative values are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_cpx, T_degreeC)
% calcPressure
% Compute Nimis and Taylor (2000) pressure for one Cpx row at one or more
% input temperatures. NaN values propagate through the calculation.
%
% Inputs:
%   data_cpx    : 1-row Clinopyroxene table
%   T_degreeC   : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

cpx = prepareCpxRow(data_cpx);

CrNumber_scalar = calcCrNumber(cpx);
aCaCrTs_scalar = calcACaCrTsCpx(cpx, CrNumber_scalar);
site = calcCpxSiteSums(cpx);
Cr2O3_wtpercent_scalar = getOptionalOxideWt(data_cpx, 'Cr2O3');

% Expand composition-dependent scalars to the temperature-vector length.
Si_cpx = repmat(cpx.Si, nT, 1);
Al_cpx = repmat(cpx.Al, nT, 1);
Fe_cpx = repmat(cpx.Fe, nT, 1);
Fe3_cpx = repmat(cpx.Fe3, nT, 1);
Mg_cpx = repmat(cpx.Mg, nT, 1);
Ca_cpx = repmat(cpx.Ca, nT, 1);
Na_cpx = repmat(cpx.Na, nT, 1);
K_cpx = repmat(cpx.K, nT, 1);
Mn_cpx = repmat(cpx.Mn, nT, 1);
Ti_cpx = repmat(cpx.Ti, nT, 1);
Cr_cpx = repmat(cpx.Cr, nT, 1);

CrNumber_cpx = repmat(CrNumber_scalar, nT, 1);
aCaCrTs_cpx = repmat(aCaCrTs_scalar, nT, 1);

cationTotal_cpx = repmat(site.cationTotal, nT, 1);
AlIV_cpx = repmat(site.AlIV, nT, 1);
T_site_sum_cpx = repmat(site.TsiteSum, nT, 1);
M1M2_site_sum_cpx = repmat(site.M1M2siteSum, nT, 1);
Cr2O3_wtpercent = repmat(Cr2O3_wtpercent_scalar, nT, 1);

% Calculate logarithm terms only where their mathematical domains are valid.
% Invalid values are represented by NaN rather than causing an error or a
% complex result.
ln_aCaCrTs_cpx = NaN(nT, 1);
validActivity = isfinite(aCaCrTs_cpx) & aCaCrTs_cpx > 0;
ln_aCaCrTs_cpx(validActivity) = log(aCaCrTs_cpx(validActivity));

CrNumber_over_T = CrNumber_cpx ./ T_K;
ln_CrNumber_over_T = NaN(nT, 1);
validCrNumberTerm = isfinite(CrNumber_over_T) & CrNumber_over_T > 0;
ln_CrNumber_over_T(validCrNumberTerm) = ...
    log(CrNumber_over_T(validCrNumberTerm));

% Equations (9) and (10). NaN inputs remain NaN and propagate naturally.
P_kbar = ...
    -(T_K ./ 126.9) .* ln_aCaCrTs_cpx ...
    + 15.483 .* ln_CrNumber_over_T ...
    + T_K ./ 71.38 ...
    + 107.8;

P_GPa = P_kbar ./ 10;

% Diagnostic applicability flags. These do not replace geological judgment
% about whether the Cpx equilibrated with garnet.
isWithinCalibrationTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 900 & T_degreeC <= 1400;
isWithinCalibrationPRange = ...
    isfinite(P_kbar) & P_kbar >= 20 & P_kbar <= 60;
isHighPressureUncertainty = ...
    isfinite(P_kbar) & P_kbar > 40;
isAboveSafeActivityLimit = ...
    isfinite(aCaCrTs_cpx) & aCaCrTs_cpx >= 0.003;
isWithinEquationDomain = ...
    isfinite(CrNumber_cpx) & CrNumber_cpx > 0 & ...
    isfinite(aCaCrTs_cpx) & aCaCrTs_cpx > 0 & ...
    isfinite(T_K) & T_K > 0;
passesSiteSumCriteria = ...
    isfinite(T_site_sum_cpx) & T_site_sum_cpx > 1.990 & ...
    isfinite(M1M2_site_sum_cpx) & M1M2_site_sum_cpx > 1.990;
passesCr2O3Screen = ...
    isnan(Cr2O3_wtpercent) | Cr2O3_wtpercent < 5;

% Pack outputs using equal-length, pre-sized vectors.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.P_uncertainty_overall_1sigma_kbar = repmat(2.3, nT, 1);

row.Si_cpx = Si_cpx;
row.Al_cpx = Al_cpx;
row.Fe_cpx = Fe_cpx;
row.Fe3_cpx = Fe3_cpx;
row.Mg_cpx = Mg_cpx;
row.Ca_cpx = Ca_cpx;
row.Na_cpx = Na_cpx;
row.K_cpx = K_cpx;
row.Mn_cpx = Mn_cpx;
row.Ti_cpx = Ti_cpx;
row.Cr_cpx = Cr_cpx;

row.CrNumber_cpx = CrNumber_cpx;
row.aCaCrTs_cpx = aCaCrTs_cpx;
row.CrNumber_over_T = CrNumber_over_T;
row.ln_aCaCrTs_cpx = ln_aCaCrTs_cpx;
row.ln_CrNumber_over_T = ln_CrNumber_over_T;

row.cationTotal_cpx = cationTotal_cpx;
row.AlIV_cpx = AlIV_cpx;
row.T_site_sum_cpx = T_site_sum_cpx;
row.M1M2_site_sum_cpx = M1M2_site_sum_cpx;
row.Cr2O3_wtpercent = Cr2O3_wtpercent;

row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinCalibrationPRange = isWithinCalibrationPRange;
row.isHighPressureUncertainty = isHighPressureUncertainty;
row.isAboveSafeActivityLimit = isAboveSafeActivityLimit;
row.isWithinEquationDomain = isWithinEquationDomain;
row.passesSiteSumCriteria = passesSiteSumCriteria;
row.passesCr2O3Screen = passesCr2O3Screen;

end

function cpx = prepareCpxRow(data_cpx)
% prepareCpxRow
% Extract one-row six-oxygen Cpx cation data. Required columns must exist,
% but present NaN values are retained. Missing optional variables are NaN.

if height(data_cpx) ~= 1
    error('Cpx input must be a 1-row table.');
end

cpx = struct();

cpx.Si = getVarOrError(data_cpx, 'Si_cation_apfu', 'Cpx');
cpx.Al = getVarOrError(data_cpx, 'Al_cation_apfu', 'Cpx');
cpx.Fe = getVarOrError(data_cpx, 'Fe_cation_apfu', 'Cpx');
cpx.Mg = getVarOrError(data_cpx, 'Mg_cation_apfu', 'Cpx');
cpx.Ca = getVarOrError(data_cpx, 'Ca_cation_apfu', 'Cpx');
cpx.Na = getVarOrError(data_cpx, 'Na_cation_apfu', 'Cpx');
cpx.K = getVarOrError(data_cpx, 'K_cation_apfu', 'Cpx');
cpx.Cr = getVarOrError(data_cpx, 'Cr_cation_apfu', 'Cpx');

cpx.Mn = getVarOrNaN(data_cpx, 'Mn_cation_apfu');
cpx.Ti = getVarOrNaN(data_cpx, 'Ti_cation_apfu');
cpx.Fe3 = getVarOrNaN(data_cpx, 'Fe3_cation_apfu');

% Reject Inf and finite negative values in all extracted cation variables.
fieldNames = fieldnames(cpx);
for i = 1:numel(fieldNames)
    value = cpx.(fieldNames{i});
    if ~isscalar(value)
        error('Cpx variable %s must be scalar in a 1-row table.', fieldNames{i});
    end
    if isinf(value) || (isfinite(value) && value < 0)
        error('Cpx contains an invalid negative or Inf value for %s.', fieldNames{i});
    end
end

% Fe_cation_apfu is treated as total Fe. Check only when Fe3 is available.
if isfinite(cpx.Fe3) && isfinite(cpx.Fe) && cpx.Fe3 > cpx.Fe
    error('Cpx contains Fe3_cation_apfu > Fe_cation_apfu.');
end

end

function CrNumber_cpx = calcCrNumber(cpx)
% calcCrNumber
% Calculate Cr# = Cr/(Cr+Al). Invalid or zero denominators return NaN.

denominator = cpx.Cr + cpx.Al;

if isfinite(denominator) && denominator > 0
    CrNumber_cpx = cpx.Cr ./ denominator;
else
    CrNumber_cpx = NaN;
end

end

function aCaCrTs_cpx = calcACaCrTsCpx(cpx, CrNumber_cpx)
% calcACaCrTsCpx
% Calculate the corrected CaCrTs activity parameter. K is added to Na as in
% the natural-sample application described on p. 548.

aCaCrTs_cpx = ...
    cpx.Cr - 0.81 .* CrNumber_cpx .* (cpx.Na + cpx.K);

end

function site = calcCpxSiteSums(cpx)
% calcCpxSiteSums
% Calculate approximate six-oxygen tetrahedral and M1+M2 cation sums for the
% analytical-quality criteria on p. 548. NaN inputs remain NaN.

site = struct();

if isnan(cpx.Si) || isnan(cpx.Al)
    site.AlIV = NaN;
else
    tetrahedralDeficit = max(0, 2 - cpx.Si);
    site.AlIV = min(cpx.Al, tetrahedralDeficit);
end

% Fe_cation_apfu is total Fe; Fe3 is not added separately.
site.cationTotal = ...
    cpx.Si + cpx.Al + cpx.Fe + cpx.Mg + cpx.Ca + ...
    cpx.Na + cpx.K + cpx.Mn + cpx.Ti + cpx.Cr;

site.TsiteSum = cpx.Si + site.AlIV;
site.M1M2siteSum = site.cationTotal - site.TsiteSum;

end

function invalidEquationTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Return composition-dependent terms outside the logarithm domains. The
% first row is sufficient because composition is constant across T input.

maxNames = 4;
invalidBuffer = strings(maxNames, 1);
nInvalid = 0;

if isempty(row)
    invalidEquationTerms = strings(0, 1);
    return;
end

checks = {
    'Cr + Al denominator', row.Cr_cpx(1) + row.Al_cpx(1), 'positive'
    'CrNumber_cpx', row.CrNumber_cpx(1), 'positive'
    'aCaCrTs_cpx', row.aCaCrTs_cpx(1), 'positive'
    'CrNumber_cpx/T_K', row.CrNumber_over_T(1), 'positive'
    };

for i = 1:size(checks, 1)
    name = checks{i, 1};
    value = checks{i, 2};
    domain = checks{i, 3};

    if strcmp(domain, 'positive')
        invalid = ~isfinite(value) || value <= 0;
    else
        invalid = ~isfinite(value) || value < 0;
    end

    if invalid
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = string(name);
    end
end

invalidEquationTerms = invalidBuffer(1:nInvalid);

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
% represented by NaN and are never replaced by zero.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);
    if ~isscalar(value)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end
else
    value = NaN;
end

end

function value = getOptionalOxideWt(tbl, oxide)
% getOptionalOxideWt
% Retrieve an optional oxide concentration used only for an applicability
% warning. Missing or nonnumeric values are represented by NaN.

columnName = findOxideColumn(tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    value = NaN;
    return;
end

raw = tbl.(columnName);
value = toScalarDoubleOrNaN(raw);

end

function columnName = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match oxide columns after removing spaces, underscores, and hyphens. Both
% "Cr2O3" and "Cr2O3_value" styles are supported.

canonicalNames = cell(size(variableNames));
for i = 1:numel(variableNames)
    canonicalNames{i} = canonicalizeName(variableNames{i});
end

oxideCanonical = canonicalizeName(oxide);
targets = {[oxideCanonical 'value'], oxideCanonical};

columnName = '';
for i = 1:numel(targets)
    index = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(index)
        columnName = variableNames{index};
        return;
    end
end

end

function canonicalName = canonicalizeName(name)
% canonicalizeName
% Convert a table-variable name to the canonical form used for matching.

canonicalName = lower(char(string(name)));
canonicalName = strrep(canonicalName, ' ', '');
canonicalName = strrep(canonicalName, '_', '');
canonicalName = strrep(canonicalName, '-', '');

end

function value = toScalarDoubleOrNaN(raw)
% toScalarDoubleOrNaN
% Convert one table value to a scalar double. Missing, empty, and
% non-convertible values become NaN and are never converted to zero.

value = NaN;

if isempty(raw)
    return;
end

if isnumeric(raw) || islogical(raw)
    value = double(raw(1));
    return;
end

if isstring(raw)
    if ismissing(raw(1))
        return;
    end
    value = str2double(raw(1));
    return;
end

if ischar(raw)
    value = str2double(string(raw));
    return;
end

if iscell(raw)
    if isempty(raw{1})
        return;
    end
    value = toScalarDoubleOrNaN(raw{1});
    return;
end

end
