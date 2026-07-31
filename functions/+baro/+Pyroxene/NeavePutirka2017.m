function results = NeavePutirka2017(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Pyroxene/NeavePutirka2017.m
% Tested with MATLAB R2024b
%
% Clinopyroxene-Liquid jadeite-in-clinopyroxene geobarometer
% Neave, D.A. and Putirka, K.D. (2017)
% American Mineralogist, 102, 777-794
% DOI: https://doi.org/10.2138/am-2017-5968
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Clinopyroxene analysis and pairs
% it with one row from a liquid-composition dataset loaded through
% liquid.readLiquidExcel. Pressure is calculated using Equation (1) of
% Neave and Putirka (2017).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Cpx-Liq pair, one output row is
% returned for every input temperature value.
%
% The optional name-value argument 'LiquidRow' may be used to specify the
% liquid row. If it is omitted, the first row of the selected liquid dataset
% is used, preserving the behavior of the original implementation.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Neave and Putirka (2017) calibrated Equation (1) using 113
% clinopyroxene-saturated experiments covering:
%
%   Pressure    : 0.001-20 kbar
%   Temperature : 950-1400 degreeC
%   Composition : ultramafic to intermediate liquids
%   H2O         : anhydrous and hydrous experiments
%
% The calibration range, experimental sources, and Equation (1) are given
% on p. 781 and in Table 1. The calibration data are reproduced with a
% standard error of estimate of approximately 1.3 kbar when experimental
% temperatures are imposed. Iterative calculations with Equation (33) of
% Putirka (2008) give an SEE of approximately 1.4 kbar (pp. 781-782).
%
% IMPORTANT APPLICATION NOTES:
%
% 1) Calculations outside 950-1400 degreeC or 0.001-20 kbar are
%    extrapolations and are reported by non-stopping fprintf warnings.
%
% 2) The barometer may be used for hydrous and anhydrous samples, but the
%    calibrated compositional domain is ultramafic to intermediate. Highly
%    alkaline, strongly evolved, or otherwise compositionally dissimilar
%    liquids require caution (pp. 781-783).
%
% 3) Caution is required at fO2 >= QFM+1. The formulation does not calculate
%    Fe3+ or aegirine (Aeg) explicitly. Na assigned to Jd may actually belong
%    to Aeg, causing pressure overestimation. High-fO2 phonolitic experiments
%    produced errors of about 10 kbar in a global test (pp. 781-783).
%
% 4) When this barometer is solved iteratively with Putirka (2008) Equation
%    (33), temperatures are systematically overestimated at <=1100 degreeC.
%    Neave and Putirka (2017) therefore advise caution for this combined
%    barometer-thermometer application below or at 1100 degreeC (p. 781).
%
% 5) Clinopyroxene and liquid must represent an equilibrium pair. Fe-Mg
%    equilibrium alone is not sufficient in mixed or disequilibrium systems
%    because NaAl-CaMg exchange may equilibrate at a different rate.
%    Multiple-component equilibrium tests involving Fe-Mg, DiHd, CaTs, and
%    Ti are recommended for natural samples (pp. 778-779 and 786-787).
%
% 6) The carrier glass or host whole-rock composition must not be assumed to
%    be the equilibrium liquid for zoned, antecrystic, xenocrystic, or mixed
%    clinopyroxene populations. Equilibrium liquid selection is a major
%    source of uncertainty (pp. 786-787).
%
% 7) Cpx stoichiometry and analytical precision must be checked. Calibration
%    clinopyroxenes have six-oxygen cation sums of 3.97-4.05 (p. 781). For
%    the Icelandic application, the authors used a stricter 3.99-4.02 filter
%    and excluded XJd < 0.01 because Na2O was near the EPMA detection limit
%    (p. 784).
%
% 8) Sector zoning and rapid-growth disequilibrium can strongly affect Na,
%    Al, Ti, Ca, Mg, and Fe. The Al(6O) >= 0.11 filter used for Icelandic
%    tholeiites is application-specific and must not be treated as a
%    universal calibration limit (pp. 784-786).
%
% 9) Na2O analytical precision is critical because Na2O is commonly only
%    0.1-0.5 wt% in clinopyroxene and may migrate under the electron beam.
%    The approximately 1.4 kbar model precision may be close to the practical
%    limit of the Jd-liquid approach (p. 783; pp. 790-791).
%
% 10) A pressure distribution with a width comparable to approximately
%     1.4 kbar may reflect model, analytical, or liquid-matching uncertainty
%     rather than true polybaric crystallization (pp. 787 and 790-791).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 950-1400 degreeC,
%   2) finite calculated pressure is outside 0.001-20 kbar,
%   3) finite input temperature is <=1100 degreeC,
%   4) a required calculation input contains NaN,
%   5) Cpx cation sum or XJd is outside the recommended diagnostic domain,
%   6) an equation term is outside its mathematical domain, or
%   7) calculated pressure is NaN, Inf, or negative.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST column of the Cpx table is treated as an identifier ("data
% code") displayed in the selection dialog.
%
% The Cpx table and liquid table must contain oxide columns recognizable as:
%
%   Cpx calculation oxides:
%     SiO2, TiO2, Al2O3, FeO or FeOt, MnO, MgO, CaO, Na2O, K2O, Cr2O3
%
%   Liquid calculation oxides:
%     SiO2, TiO2, Al2O3, FeO or FeOt, MnO, MgO, CaO, Na2O, K2O,
%     V2O3, Cr2O3, NiO, P2O5, SO3, Fe2O3
%
% For auxiliary liquid oxides, an absent column is interpreted as zero.
% However, when a column exists and its selected value is NaN, that NaN is
% retained and propagated. Required major oxides must exist; NaN values are
% allowed, retained, and reported.
%
% Finite calculation inputs must be non-negative. NaN is allowed and is not
% replaced by zero. Inf and finite negative values are prohibited. Zero is
% retained; if it makes a logarithm or ratio undefined, the corresponding
% pressure remains NaN and is reported.
%
% Liquid F and Cl are explicitly excluded from cationTotal_liq and are also
% excluded from NaN-input warnings, as requested.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION
%
% Equation (1) of Neave and Putirka (2017; p. 781):
%
%   P(kbar) =
%     -26.27
%     + 39.16*T(K)/10^4 *
%       ln[XJd_cpx /
%          {XNaO0.5_liq*XAlO1.5_liq*(XSiO2_liq)^2}]
%     - 4.22*ln(XDiHd_cpx)
%     + 78.43*XAlO1.5_liq
%     + 393.81*(XNaO0.5_liq*XKO0.5_liq)^2
%
% Clinopyroxene and liquid components follow the calculation procedures
% summarized in Putirka (2008), as adopted by Neave and Putirka (2017).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = NeavePutirka2017(rawdata_struct, T_degreeC)
%   results = NeavePutirka2017(..., 'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing a Cpx table
%   T_degreeC      : non-negative numeric scalar or vector in degreeC.
%                    NaN is allowed and retained; Inf is prohibited.
%
% Name-value option:
%   LiquidRow      : positive integer scalar or [].
%
% Output:
%   results : table containing one row per temperature value for every
%             selected Cpx-Liq pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same implementation.
if nargin < 2
    error('NeavePutirka2017 requires (rawdata_struct, T_degreeC).');
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

%% Options
ip = inputParser;
ip.FunctionName = 'NeavePutirka2017';
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve datasets and constants
% Load molecular weights/cation numbers and the liquid dataset. The source
% tables are not modified.
disp('=== Step 1: Preparing Cpx and Liquid datasets ===');

MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();

if isempty(liqAll) || ~istable(liqAll)
    error('Selected liquid dataset must be a non-empty table.');
end

dataset_cpx = rawdata_struct.Cpx;

disp('=== Preparing Cpx and Liquid datasets has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max(16, height(dataset_cpx));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degreeC = 950;
calibrationT_max_degreeC = 1400;
calibrationP_min_kbar = 0.001;
calibrationP_max_kbar = 20;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);
lowTemperatureCaution = isfinite(T_degreeC) & T_degreeC <= 1100;
temperatureRangeWarningIssued = false;
lowTemperatureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Cpx) ===');

while true
    % ----- Clinopyroxene selection -----
    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Clinopyroxene (Cpx) data:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIdx_cpx);
    disp(['Cpx selected: ' char(string(selectedCode_cpx))]);

    % ----- Liquid-row selection -----
    disp('=== Step 4: Selecting a Liquid row ===');

    if isempty(liquidRowOpt)
        selectedIdx_liq = 1;
        if height(liqAll) > 1
            fprintf(2, ...
                ['CAUTION: The selected liquid dataset contains %d rows. ' ...
                 'LiquidRow was not specified, so row 1 is used, preserving ' ...
                 'the original NeavePutirka2017 behavior.\n'], ...
                height(liqAll));
        end
    else
        selectedIdx_liq = liquidRowOpt;
        if selectedIdx_liq > height(liqAll)
            error(['Requested LiquidRow (%d) exceeds the number of rows in ' ...
                   'the selected liquid dataset (%d).'], ...
                selectedIdx_liq, height(liqAll));
        end
    end

    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);
    selectedData_liq = liqAll(selectedIdx_liq, :);
    selectedCode_liq = getLiquidLabel(selectedData_liq, selectedIdx_liq);

    disp(['Liquid selected: ' char(selectedCode_liq)]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the pressure ===');

    % List NaN in calculation inputs. F and Cl are intentionally excluded.
    nanInputNames = findNaNInputs( ...
        selectedData_cpx, selectedData_liq, T_degreeC);

    % Reject Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(selectedData_cpx, selectedData_liq);

    row = calcPressure( ...
        selectedData_cpx, selectedData_liq, T_degreeC, MWinfo);

    % Repeat identifiers for all temperature rows.
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row.dataRow_liq = repmat(selectedIdx_liq, height(row), 1);
    row.dataCode_liq = repmat(string(selectedCode_liq), height(row), 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, ...
        {'dataCode_cpx', 'dataRow_liq', 'dataCode_liq'}, 'Before', 1);

    % Buffer one result block per selected Cpx-Liq pair.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo calculated pressure.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_cpx)) ' & ' char(selectedCode_liq) ...
            ': P = ' num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_cpx)) ' & ' char(selectedCode_liq) ...
            ': P = ' num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Input-temperature warnings are common to all selected pairs and are
    % therefore issued only once per function call.
    if any(temperatureOutsideCalibration) && ~temperatureRangeWarningIssued
        finiteTemperatureValues = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the Neave and Putirka ' ...
             '(2017) calibration range of 950-1400 degreeC (p. 781). ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperatureValues), ...
            min(finiteTemperatureValues), ...
            max(finiteTemperatureValues));
        temperatureRangeWarningIssued = true;
    end

    if any(lowTemperatureCaution) && ~lowTemperatureWarningIssued
        finiteTemperatureValues = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['CAUTION: %d of %d finite input temperature point(s) are at or ' ...
             'below 1100 degreeC. Neave and Putirka (2017) advise caution ' ...
             'when this barometer is combined iteratively with Putirka ' ...
             '(2008) Equation 33 at <=1100 degreeC because temperature is ' ...
             'systematically overestimated (p. 781). Input finite range = ' ...
             '%.4g-%.4g degreeC.\n'], ...
            sum(lowTemperatureCaution), ...
            numel(finiteTemperatureValues), ...
            min(finiteTemperatureValues), ...
            max(finiteTemperatureValues));
        lowTemperatureWarningIssued = true;
    end

    % Calculated-pressure calibration warning.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the Neave and Putirka ' ...
             '(2017) calibration range of 0.001-20 kbar (p. 781). ' ...
             '%d of %d finite pressure point(s) are outside the range; ' ...
             'calculated finite range = %.4g-%.4g kbar for %s & %s.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_cpx)), ...
            char(selectedCode_liq));
    end

    % Report exact NaN input names. F and Cl are not checked here.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Neave-Putirka calculation ' ...
             'input(s) for %s & %s: %s.\n' ...
             '         NaN values were retained and propagated; calculated ' ...
             'pressure may remain NaN.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(selectedCode_liq), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Cpx analytical and stoichiometric diagnostics.
    cationSumScalar = row.cationSum_cpx(1);
    if isfinite(cationSumScalar) && ...
            (cationSumScalar < 3.97 || cationSumScalar > 4.05)
        fprintf(2, ...
            ['WARNING: Cpx six-oxygen cation sum = %.4g is outside the ' ...
             '3.97-4.05 range represented by the calibration Cpx analyses ' ...
             '(Neave and Putirka, 2017; p. 781) for %s & %s.\n'], ...
            cationSumScalar, ...
            char(string(selectedCode_cpx)), ...
            char(selectedCode_liq));
    elseif isfinite(cationSumScalar) && ...
            (cationSumScalar < 3.99 || cationSumScalar > 4.02)
        fprintf(2, ...
            ['CAUTION: Cpx six-oxygen cation sum = %.4g is outside the ' ...
             'stricter 3.99-4.02 analytical filter used for the Icelandic ' ...
             'natural-sample application (p. 784), although it remains ' ...
             'within the broader calibration-data range, for %s & %s.\n'], ...
            cationSumScalar, ...
            char(string(selectedCode_cpx)), ...
            char(selectedCode_liq));
    end

    XJdScalar = row.XJd_cpx(1);
    if isfinite(XJdScalar) && XJdScalar < 0.01
        fprintf(2, ...
            ['CAUTION: XJd_cpx = %.4g is below the XJd >= 0.01 analytical ' ...
             'filter used by Neave and Putirka (2017; p. 784). Na2O may be ' ...
             'near the EPMA detection limit for %s & %s.\n'], ...
            XJdScalar, ...
            char(string(selectedCode_cpx)), ...
            char(selectedCode_liq));
    end

    % Report invalid mathematical terms without stopping the calculation.
    invalidEquationTerms = findInvalidEquationTerms(row);
    if ~isempty(invalidEquationTerms)
        fprintf(2, ...
            ['WARNING: Neave-Putirka equation term(s) are outside their ' ...
             'mathematical domain for %s & %s: %s.\n' ...
             '         The affected pressure values were retained as NaN/Inf ' ...
             'where applicable.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(selectedCode_liq), ...
            char(strjoin(invalidEquationTerms, ', ')));
    end

    % Retain and report all non-finite results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(selectedCode_liq), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & %s ' ...
             '(%d of %d points). Values were retained for diagnostic purposes.\n'], ...
            char(string(selectedCode_cpx)), ...
            char(selectedCode_liq), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Cpx selection using the same Liquid dataset ?', ...
        'NeavePutirka2017', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered blocks once.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct('liquid', metaLiq);
disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_cpx, data_liq, T_degreeC)
% findNaNInputs
% Return names of calculation inputs containing NaN. F and Cl are
% intentionally excluded. Existing NaN values are never replaced by zero.

cpxOxides = {'SiO2', 'TiO2', 'Al2O3', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'Cr2O3'};
liqOxides = {'SiO2', 'TiO2', 'Al2O3', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', ...
    'P2O5', 'SO3', 'Fe2O3'};

maxNames = 3 + numel(cpxOxides) + numel(liqOxides);
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

for i = 1:numel(cpxOxides)
    oxide = cpxOxides{i};
    columnName = findOxideColumn(data_cpx.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDoubleOrNaN(data_cpx.(columnName));
        if isnan(value)
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = "Cpx." + string(columnName);
        end
    end
end

[FeValueCpx, FeLabelCpx] = getFeOValue(data_cpx, 'Cpx');
if isnan(FeValueCpx)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = string(FeLabelCpx);
end

for i = 1:numel(liqOxides)
    oxide = liqOxides{i};
    columnName = findOxideColumn(data_liq.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDoubleOrNaN(data_liq.(columnName));
        if isnan(value)
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = "Liq." + string(columnName);
        end
    end
end

[FeValueLiq, FeLabelLiq] = getFeOValue(data_liq, 'Liq');
if isnan(FeValueLiq)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = string(FeLabelLiq);
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_cpx, data_liq)
% validateNonNegativeInputs
% Reject Inf and finite negative values in inputs used to calculate Cpx
% components or liquid cation fractions. F and Cl are excluded because they
% are not used in cationTotal_liq.

cpxOxides = {'SiO2', 'TiO2', 'Al2O3', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'Cr2O3'};
liqOxides = {'SiO2', 'TiO2', 'Al2O3', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', ...
    'P2O5', 'SO3', 'Fe2O3'};

maxNames = 2 + numel(cpxOxides) + numel(liqOxides);
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

for i = 1:numel(cpxOxides)
    oxide = cpxOxides{i};
    columnName = findOxideColumn(data_cpx.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDoubleOrNaN(data_cpx.(columnName));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = ...
                "Cpx." + string(columnName);
        end
    end
end

[FeValueCpx, FeLabelCpx] = getFeOValue(data_cpx, 'Cpx');
if isinf(FeValueCpx) || (isfinite(FeValueCpx) && FeValueCpx < 0)
    nInvalidInputs = nInvalidInputs + 1;
    invalidInputBuffer(nInvalidInputs) = string(FeLabelCpx);
end

for i = 1:numel(liqOxides)
    oxide = liqOxides{i};
    columnName = findOxideColumn(data_liq.Properties.VariableNames, oxide);
    if ~isempty(columnName)
        value = toScalarDoubleOrNaN(data_liq.(columnName));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = ...
                "Liq." + string(columnName);
        end
    end
end

[FeValueLiq, FeLabelLiq] = getFeOValue(data_liq, 'Liq');
if isinf(FeValueLiq) || (isfinite(FeValueLiq) && FeValueLiq < 0)
    nInvalidInputs = nInvalidInputs + 1;
    invalidInputBuffer(nInvalidInputs) = string(FeLabelLiq);
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['NeavePutirka2017: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative values are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_cpx, data_liq, T_degreeC, MWinfo)
% calcPressure
% Compute Neave and Putirka (2017) Equation (1) pressure for one Cpx row and
% one liquid row at one or more input temperatures. NaN values propagate.
%
% Inputs:
%   data_cpx    : 1-row Clinopyroxene table
%   data_liq    : 1-row Liquid table
%   T_degreeC   : scalar or vector temperature in degreeC
%   MWinfo      : molecular-weight and cation-number structure
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

cpx = prepareCpxRow(data_cpx, MWinfo);
liq = prepareLiquidRow(data_liq, MWinfo);

% Equation terms. Invalid logarithm arguments are represented by NaN rather
% than converted to zero or complex values.
logArgJd_scalar = cpx.XJd ./ ...
    (liq.XNaO0_5 .* liq.XAlO1_5 .* (liq.XSiO2 .^ 2));

logTermJd_scalar = NaN;
if isfinite(logArgJd_scalar) && logArgJd_scalar > 0
    logTermJd_scalar = log(logArgJd_scalar);
elseif isinf(logArgJd_scalar) && logArgJd_scalar > 0
    logTermJd_scalar = Inf;
end

logTermDiHd_scalar = NaN;
if isfinite(cpx.XDiHd) && cpx.XDiHd > 0
    logTermDiHd_scalar = log(cpx.XDiHd);
elseif isinf(cpx.XDiHd) && cpx.XDiHd > 0
    logTermDiHd_scalar = Inf;
end

% Expand composition-dependent scalars to the temperature-vector length.
XSi_cpx = repmat(cpx.XSi, nT, 1);
XTi_cpx = repmat(cpx.XTi, nT, 1);
XAl_cpx = repmat(cpx.XAl, nT, 1);
XFe_cpx = repmat(cpx.XFe, nT, 1);
XMn_cpx = repmat(cpx.XMn, nT, 1);
XMg_cpx = repmat(cpx.XMg, nT, 1);
XCa_cpx = repmat(cpx.XCa, nT, 1);
XNa_cpx = repmat(cpx.XNa, nT, 1);
XK_cpx = repmat(cpx.XK, nT, 1);
XCr_cpx = repmat(cpx.XCr, nT, 1);
cationSum_cpx = repmat(cpx.cationSum, nT, 1);

XAlIV_cpx = repmat(cpx.XAlIV, nT, 1);
XAlVI_cpx = repmat(cpx.XAlVI, nT, 1);
XFe3_cpx = repmat(cpx.XFe3, nT, 1);
XJd_cpx = repmat(cpx.XJd, nT, 1);
XCaTs_cpx = repmat(cpx.XCaTs, nT, 1);
XCaTi_cpx = repmat(cpx.XCaTi, nT, 1);
XCrCaTs_cpx = repmat(cpx.XCrCaTs, nT, 1);
XDiHd_cpx = repmat(cpx.XDiHd, nT, 1);
XEnFs_cpx = repmat(cpx.XEnFs, nT, 1);

SiO2_liq = repmat(liq.SiO2, nT, 1);
TiO2_liq = repmat(liq.TiO2, nT, 1);
Al2O3_liq = repmat(liq.Al2O3, nT, 1);
FeO_liq = repmat(liq.FeO, nT, 1);
MnO_liq = repmat(liq.MnO, nT, 1);
MgO_liq = repmat(liq.MgO, nT, 1);
CaO_liq = repmat(liq.CaO, nT, 1);
Na2O_liq = repmat(liq.Na2O, nT, 1);
K2O_liq = repmat(liq.K2O, nT, 1);
V2O3_liq = repmat(liq.V2O3, nT, 1);
Cr2O3_liq = repmat(liq.Cr2O3, nT, 1);
NiO_liq = repmat(liq.NiO, nT, 1);
P2O5_liq = repmat(liq.P2O5, nT, 1);
SO3_liq = repmat(liq.SO3, nT, 1);
Fe2O3_liq = repmat(liq.Fe2O3, nT, 1);

cationTotal_liq = repmat(liq.cationTotal, nT, 1);
XSiO2_liq = repmat(liq.XSiO2, nT, 1);
XTiO2_liq = repmat(liq.XTiO2, nT, 1);
XAlO1_5_liq = repmat(liq.XAlO1_5, nT, 1);
XFeO_liq = repmat(liq.XFeO, nT, 1);
XMnO_liq = repmat(liq.XMnO, nT, 1);
XMgO_liq = repmat(liq.XMgO, nT, 1);
XCaO_liq = repmat(liq.XCaO, nT, 1);
XNaO0_5_liq = repmat(liq.XNaO0_5, nT, 1);
XKO0_5_liq = repmat(liq.XKO0_5, nT, 1);

logArg_Eq1 = repmat(logArgJd_scalar, nT, 1);
logTerm_Jd_Eq1 = repmat(logTermJd_scalar, nT, 1);
logTerm_DiHd_Eq1 = repmat(logTermDiHd_scalar, nT, 1);

% Equation (1). NaN inputs remain NaN and propagate naturally.
P_kbar = ...
    -26.27 ...
    + (39.16 .* T_K ./ 1e4) .* logTerm_Jd_Eq1 ...
    - 4.22 .* logTerm_DiHd_Eq1 ...
    + 78.43 .* XAlO1_5_liq ...
    + 393.81 .* ((XNaO0_5_liq .* XKO0_5_liq) .^ 2);

P_GPa = P_kbar ./ 10;

isWithinCalibrationTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 950 & T_degreeC <= 1400;
isWithinCalibrationPRange = ...
    isfinite(P_kbar) & P_kbar >= 0.001 & P_kbar <= 20;
isLowTemperatureCaution = ...
    isfinite(T_degreeC) & T_degreeC <= 1100;
isWithinCpxCalibrationStoichiometry = ...
    isfinite(cationSum_cpx) & cationSum_cpx >= 3.97 & cationSum_cpx <= 4.05;
passesStrictNaturalCpxFilter = ...
    isfinite(cationSum_cpx) & cationSum_cpx >= 3.99 & cationSum_cpx <= 4.02 & ...
    isfinite(XJd_cpx) & XJd_cpx >= 0.01;
isWithinEquationDomain = ...
    isfinite(XJd_cpx) & XJd_cpx > 0 & ...
    isfinite(XDiHd_cpx) & XDiHd_cpx > 0 & ...
    isfinite(XSiO2_liq) & XSiO2_liq > 0 & ...
    isfinite(XAlO1_5_liq) & XAlO1_5_liq > 0 & ...
    isfinite(XNaO0_5_liq) & XNaO0_5_liq > 0 & ...
    isfinite(XKO0_5_liq) & XKO0_5_liq >= 0 & ...
    isfinite(logArg_Eq1) & logArg_Eq1 > 0;

% Pack outputs using equal-length vectors.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% Backward-compatible names from the original implementation.
row.PEq1_kbar = P_kbar;
row.PEq1_GPa = P_GPa;
row.P_uncertainty_SEE_kbar = repmat(1.4, nT, 1);

row.XSi_cpx = XSi_cpx;
row.XTi_cpx = XTi_cpx;
row.XAl_cpx = XAl_cpx;
row.XFe_cpx = XFe_cpx;
row.XMn_cpx = XMn_cpx;
row.XMg_cpx = XMg_cpx;
row.XCa_cpx = XCa_cpx;
row.XNa_cpx = XNa_cpx;
row.XK_cpx = XK_cpx;
row.XCr_cpx = XCr_cpx;
row.cationSum_cpx = cationSum_cpx;

row.XAlIV_cpx = XAlIV_cpx;
row.XAlVI_cpx = XAlVI_cpx;
row.XFe3_cpx = XFe3_cpx;
row.XJd_cpx = XJd_cpx;
row.XCaTs_cpx = XCaTs_cpx;
row.XCaTi_cpx = XCaTi_cpx;
row.XCrCaTs_cpx = XCrCaTs_cpx;
row.XDiHd_cpx = XDiHd_cpx;
row.XEnFs_cpx = XEnFs_cpx;

row.SiO2_liq = SiO2_liq;
row.TiO2_liq = TiO2_liq;
row.Al2O3_liq = Al2O3_liq;
row.FeO_liq = FeO_liq;
row.MnO_liq = MnO_liq;
row.MgO_liq = MgO_liq;
row.CaO_liq = CaO_liq;
row.Na2O_liq = Na2O_liq;
row.K2O_liq = K2O_liq;
row.V2O3_liq = V2O3_liq;
row.Cr2O3_liq = Cr2O3_liq;
row.NiO_liq = NiO_liq;
row.P2O5_liq = P2O5_liq;
row.SO3_liq = SO3_liq;
row.Fe2O3_liq = Fe2O3_liq;

row.cationTotal_liq = cationTotal_liq;
row.XSiO2_liq = XSiO2_liq;
row.XTiO2_liq = XTiO2_liq;
row.XAlO1_5_liq = XAlO1_5_liq;
row.XFeO_liq = XFeO_liq;
row.XMnO_liq = XMnO_liq;
row.XMgO_liq = XMgO_liq;
row.XCaO_liq = XCaO_liq;
row.XNaO0_5_liq = XNaO0_5_liq;
row.XKO0_5_liq = XKO0_5_liq;

row.logArg_Eq1 = logArg_Eq1;
row.logTerm_Jd_Eq1 = logTerm_Jd_Eq1;
row.logTerm_DiHd_Eq1 = logTerm_DiHd_Eq1;

row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinCalibrationPRange = isWithinCalibrationPRange;
row.isLowTemperatureCaution = isLowTemperatureCaution;
row.isWithinCpxCalibrationStoichiometry = ...
    isWithinCpxCalibrationStoichiometry;
row.passesStrictNaturalCpxFilter = passesStrictNaturalCpxFilter;
row.isWithinEquationDomain = isWithinEquationDomain;

end

function cpx = prepareCpxRow(data_cpx, MWinfo)
% prepareCpxRow
% Calculate six-oxygen clinopyroxene cations and components. Existing NaN
% values remain NaN. Required columns must exist; FeO may be supplied as
% FeO or FeOt.

if height(data_cpx) ~= 1
    error('Cpx input must be a 1-row table.');
end

SiO2 = getOxideRequired(data_cpx, 'SiO2', 'Cpx');
TiO2 = getOxideRequired(data_cpx, 'TiO2', 'Cpx');
Al2O3 = getOxideRequired(data_cpx, 'Al2O3', 'Cpx');
[FeO, ~] = getFeOValue(data_cpx, 'Cpx');
MnO = getOxideRequired(data_cpx, 'MnO', 'Cpx');
MgO = getOxideRequired(data_cpx, 'MgO', 'Cpx');
CaO = getOxideRequired(data_cpx, 'CaO', 'Cpx');
Na2O = getOxideRequired(data_cpx, 'Na2O', 'Cpx');
K2O = getOxideRequired(data_cpx, 'K2O', 'Cpx');
Cr2O3 = getOxideRequired(data_cpx, 'Cr2O3', 'Cpx');

molProp.SiO2 = SiO2 ./ MWinfo.MW.SiO2;
molProp.TiO2 = TiO2 ./ MWinfo.MW.TiO2;
molProp.Al2O3 = Al2O3 ./ MWinfo.MW.Al2O3;
molProp.FeO = FeO ./ MWinfo.MW.FeO;
molProp.MnO = MnO ./ MWinfo.MW.MnO;
molProp.MgO = MgO ./ MWinfo.MW.MgO;
molProp.CaO = CaO ./ MWinfo.MW.CaO;
molProp.Na2O = Na2O ./ MWinfo.MW.Na2O;
molProp.K2O = K2O ./ MWinfo.MW.K2O;
molProp.Cr2O3 = Cr2O3 ./ MWinfo.MW.Cr2O3;

oxygenSum = ...
    2 .* molProp.SiO2 + ...
    2 .* molProp.TiO2 + ...
    3 .* molProp.Al2O3 + ...
    molProp.FeO + ...
    molProp.MnO + ...
    molProp.MgO + ...
    molProp.CaO + ...
    molProp.Na2O + ...
    molProp.K2O + ...
    3 .* molProp.Cr2O3;

oxygenRenormalizationFactor = 6 ./ oxygenSum;

XSi = molProp.SiO2 .* oxygenRenormalizationFactor;
XTi = molProp.TiO2 .* oxygenRenormalizationFactor;
XAl = 2 .* molProp.Al2O3 .* oxygenRenormalizationFactor;
XFe = molProp.FeO .* oxygenRenormalizationFactor;
XMn = molProp.MnO .* oxygenRenormalizationFactor;
XMg = molProp.MgO .* oxygenRenormalizationFactor;
XCa = molProp.CaO .* oxygenRenormalizationFactor;
XNa = 2 .* molProp.Na2O .* oxygenRenormalizationFactor;
XK = 2 .* molProp.K2O .* oxygenRenormalizationFactor;
XCr = 2 .* molProp.Cr2O3 .* oxygenRenormalizationFactor;

cationSum = XSi + XTi + XAl + XFe + XMn + XMg + XCa + XNa + XK + XCr;

% Component calculation follows the original implementation. max/min
% operations are part of the component-allocation scheme; NaN inputs remain
% NaN in MATLAB element-wise max/min operations.
XAlIV = max(2 - XSi, 0);
XAlVI = max(XAl - XAlIV, 0);

XFe3 = XNa + XAlIV - XAlVI - 2 .* XTi - XCr;
XFe3 = max(XFe3, 0);

XJd = min(XAlVI, XNa);
XJd = max(XJd, 0);

XCaTs = max(XAlVI - XJd, 0);

if isnan(XAlIV) || isnan(XCaTs)
    XCaTi = NaN;
elseif XAlIV > XCaTs
    XCaTi = (XAlIV - XCaTs) ./ 2;
else
    XCaTi = 0;
end
XCaTi = max(XCaTi, 0);

XCrCaTs = max(XCr ./ 2, 0);
XDiHd = XCa - XCaTi - XCaTs - XCrCaTs;
XDiHd = max(XDiHd, 0);

XEnFs = (XFe + XMg - XDiHd) ./ 2;
XEnFs = max(XEnFs, 0);

cpx = struct();
cpx.XSi = XSi;
cpx.XTi = XTi;
cpx.XAl = XAl;
cpx.XFe = XFe;
cpx.XMn = XMn;
cpx.XMg = XMg;
cpx.XCa = XCa;
cpx.XNa = XNa;
cpx.XK = XK;
cpx.XCr = XCr;
cpx.cationSum = cationSum;
cpx.XAlIV = XAlIV;
cpx.XAlVI = XAlVI;
cpx.XFe3 = XFe3;
cpx.XJd = XJd;
cpx.XCaTs = XCaTs;
cpx.XCaTi = XCaTi;
cpx.XCrCaTs = XCrCaTs;
cpx.XDiHd = XDiHd;
cpx.XEnFs = XEnFs;

end

function liq = prepareLiquidRow(data_liq, MWinfo)
% prepareLiquidRow
% Calculate liquid cation fractions. F and Cl are not read and are excluded
% from cationTotal_liq. Existing NaN values in present columns propagate.
% Absent auxiliary-oxide columns are represented by zero.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liq = struct();
liq.SiO2 = getOxideRequired(data_liq, 'SiO2', 'Liq');
liq.TiO2 = getOxideRequired(data_liq, 'TiO2', 'Liq');
liq.Al2O3 = getOxideRequired(data_liq, 'Al2O3', 'Liq');
[liq.FeO, ~] = getFeOValue(data_liq, 'Liq');
liq.MnO = getOxideRequired(data_liq, 'MnO', 'Liq');
liq.MgO = getOxideRequired(data_liq, 'MgO', 'Liq');
liq.CaO = getOxideRequired(data_liq, 'CaO', 'Liq');
liq.Na2O = getOxideRequired(data_liq, 'Na2O', 'Liq');
liq.K2O = getOxideRequired(data_liq, 'K2O', 'Liq');

liq.V2O3 = getOxideOptionalZeroIfAbsent(data_liq, 'V2O3');
liq.Cr2O3 = getOxideOptionalZeroIfAbsent(data_liq, 'Cr2O3');
liq.NiO = getOxideOptionalZeroIfAbsent(data_liq, 'NiO');
liq.P2O5 = getOxideOptionalZeroIfAbsent(data_liq, 'P2O5');
liq.SO3 = getOxideOptionalZeroIfAbsent(data_liq, 'SO3');
liq.Fe2O3 = getOxideOptionalZeroIfAbsent(data_liq, 'Fe2O3');

n.SiO2 = liq.SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
n.TiO2 = liq.TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
n.Al2O3 = liq.Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n.FeO = liq.FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
n.MnO = liq.MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
n.MgO = liq.MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
n.CaO = liq.CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
n.Na2O = liq.Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
n.K2O = liq.K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
n.V2O3 = liq.V2O3 .* MWinfo.Cat.V2O3 ./ MWinfo.MW.V2O3;
n.Cr2O3 = liq.Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
n.NiO = liq.NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;
n.P2O5 = liq.P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;
n.SO3 = liq.SO3 .* MWinfo.Cat.SO3 ./ MWinfo.MW.SO3;
n.Fe2O3 = liq.Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

% F and Cl are explicitly excluded.
liq.cationTotal = ...
    n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + ...
    n.MgO + n.CaO + n.Na2O + n.K2O + n.V2O3 + ...
    n.Cr2O3 + n.NiO + n.P2O5 + n.SO3 + n.Fe2O3;

liq.XSiO2 = n.SiO2 ./ liq.cationTotal;
liq.XTiO2 = n.TiO2 ./ liq.cationTotal;
liq.XAlO1_5 = n.Al2O3 ./ liq.cationTotal;
liq.XFeO = n.FeO ./ liq.cationTotal;
liq.XMnO = n.MnO ./ liq.cationTotal;
liq.XMgO = n.MgO ./ liq.cationTotal;
liq.XCaO = n.CaO ./ liq.cationTotal;
liq.XNaO0_5 = n.Na2O ./ liq.cationTotal;
liq.XKO0_5 = n.K2O ./ liq.cationTotal;

end

function invalidEquationTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Return names of composition-dependent terms outside the logarithm or
% normalization domains. The first row is sufficient because these terms
% are identical for all temperatures.

maxNames = 8;
invalidBuffer = strings(maxNames, 1);
nInvalid = 0;

if isempty(row)
    invalidEquationTerms = strings(0, 1);
    return;
end

checks = {
    'cationTotal_liq', row.cationTotal_liq(1), 'positive'
    'XJd_cpx', row.XJd_cpx(1), 'positive'
    'XDiHd_cpx', row.XDiHd_cpx(1), 'positive'
    'XSiO2_liq', row.XSiO2_liq(1), 'positive'
    'XAlO1_5_liq', row.XAlO1_5_liq(1), 'positive'
    'XNaO0_5_liq', row.XNaO0_5_liq(1), 'positive'
    'XKO0_5_liq', row.XKO0_5_liq(1), 'nonnegative'
    'logArg_Eq1', row.logArg_Eq1(1), 'positive'
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

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Repeat commonly used liquid identifiers to match the result-table height.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    row.liq_Index = repmat(data_liq.('Index'), nRows, 1);
end
if any(strcmp(variableNames, 'Experiment'))
    row.liq_Experiment = repmat(string(data_liq.('Experiment')), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    row.liq_Citation = repmat(string(data_liq.('Citation')), nRows, 1);
end

end

function label = getLiquidLabel(data_liq, rowIndex)
% getLiquidLabel
% Construct a readable liquid-row label without requiring a specific ID
% column.

variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Experiment'))
    label = "Liquid " + string(data_liq.('Experiment'));
elseif any(strcmp(variableNames, 'Index'))
    label = "Liquid Index " + string(data_liq.('Index'));
else
    label = "Liquid Row " + string(rowIndex);
end

end

function value = getOxideRequired(data_tbl, oxide, datasetLabel)
% getOxideRequired
% Retrieve an oxide from a required column. A present NaN value is retained.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    error('%s table must contain variable: %s', datasetLabel, oxide);
end

value = toScalarDoubleOrNaN(data_tbl.(columnName));

end

function value = getOxideOptionalZeroIfAbsent(data_tbl, oxide)
% getOxideOptionalZeroIfAbsent
% Return zero only when an auxiliary oxide column is absent. If the column
% exists and contains NaN, retain NaN.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    value = 0;
else
    value = toScalarDoubleOrNaN(data_tbl.(columnName));
end

end

function [value, label] = getFeOValue(data_tbl, datasetLabel)
% getFeOValue
% Retrieve FeO, falling back to FeOt. If both columns are present, a finite
% FeO value is preferred; otherwise FeOt is used. If both exist but are NaN,
% NaN is retained. At least one column must exist.

FeOColumn = findOxideColumn(data_tbl.Properties.VariableNames, 'FeO');
FeOtColumn = findOxideColumn(data_tbl.Properties.VariableNames, 'FeOt');

if isempty(FeOColumn) && isempty(FeOtColumn)
    error('%s table must contain FeO or FeOt.', datasetLabel);
end

FeOValue = NaN;
FeOtValue = NaN;

if ~isempty(FeOColumn)
    FeOValue = toScalarDoubleOrNaN(data_tbl.(FeOColumn));
end
if ~isempty(FeOtColumn)
    FeOtValue = toScalarDoubleOrNaN(data_tbl.(FeOtColumn));
end

if isfinite(FeOValue) || isinf(FeOValue)
    value = FeOValue;
    label = string(datasetLabel) + "." + string(FeOColumn);
elseif isfinite(FeOtValue) || isinf(FeOtValue)
    value = FeOtValue;
    label = string(datasetLabel) + "." + string(FeOtColumn);
else
    value = NaN;
    if ~isempty(FeOColumn) && ~isempty(FeOtColumn)
        label = string(datasetLabel) + "." + ...
            string(FeOColumn) + "/" + string(FeOtColumn);
    elseif ~isempty(FeOColumn)
        label = string(datasetLabel) + "." + string(FeOColumn);
    else
        label = string(datasetLabel) + "." + string(FeOtColumn);
    end
end

end

function columnName = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match oxide columns after removing spaces, underscores, and hyphens.
% Both "SiO2" and "SiO2_value" styles are supported.

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
% non-convertible values become NaN; they are never converted to zero.

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
