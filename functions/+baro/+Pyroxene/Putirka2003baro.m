function results = Putirka2003baro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Clinopyroxene_Liquid/Putirka2003baro.m
% Tested with MATLAB R2024b
%
% Clinopyroxene-Liquid barometer
% Putirka, K.D., Mikaelian, H., Ryerson, F. and Shaw, H. (2003)
% American Mineralogist, 88, 1542-1554
% DOI: https://doi.org/10.2138/am-2003-1017
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Clinopyroxene analysis selected from
% rawdata_struct.Cpx with one Liquid analysis loaded by
% liquid.readLiquidExcel and calculates pressure using Model A of Putirka et
% al. (2003; Table 4, p. 1546).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Clinopyroxene-Liquid pair, one
% output row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Putirka et al. (2003) calibrated Model A using 77 experiments. For the
% regression dataset, R^2 = 0.97 and the standard error of estimate is
% 1.7 kbar (Fig. 2a and Table 4, p. 1546).
%
% Independent test data span approximately 1 bar to 110 kbar. Model A
% recovers pressures from approximately 1 bar to 40 kbar without a
% systematic offset, with an overall test-data SEE of approximately 4 kbar.
% At pressures above approximately 50 kbar, Model A systematically
% underestimates pressure. The authors therefore retained the formulation
% optimized for low to moderate pressures (pp. 1548-1549). In this
% implementation:
%
%   Recommended pressure envelope : 0-40 kbar
%   Strong high-pressure caution  : P > 50 kbar
%
% The new high-SiO2 and volatile-bearing experiments reported directly in
% Putirka et al. (2003) cover approximately:
%
%   Temperature : 850-1300 degreeC
%   Pressure    : 10-35 kbar
%   Liquid SiO2 : up to 71.3 wt%
%
% These values are summarized in Tables 1-3 and the experimental discussion
% on pp. 1542-1544 and 1547-1548. The 850-1300 degreeC interval is used here
% as a non-stopping warning envelope for the new experiments; it is not
% presented by the authors as a strict universal temperature limit for the
% complete multi-study regression dataset.
%
% IMPORTANT HYDROUS-SYSTEM CAUTION: The calibration includes nominally
% hydrous and volatile-bearing compositions, but Model A substantially
% underestimates pressure for some hydrous datasets. It gives negative mean
% pressures for the 1-2 kbar Sisson and Grove experiments and approximately
% 7 kbar for the 12 kbar Muentener et al. experiments. The cause was not
% resolved, and including those datasets in the regression degraded the fit
% elsewhere (p. 1550). Hydrous-system results therefore require particular
% caution.
%
% Clinopyroxene and Liquid must represent an equilibrium pair. A whole-rock
% composition is an acceptable liquid proxy only when Clinopyroxene formed
% near the liquidus, its modal abundance is small, and later magma mixing,
% wall-rock assimilation, or extensive fractionation did not modify the
% liquid. Otherwise, a matrix/glass composition or a reconstructed liquid
% should be paired with a Clinopyroxene rim in contact with that liquid
% (pp. 1550-1552).
%
% Putirka et al. (2003) recommend checking Clinopyroxene-Liquid Fe-Mg
% equilibrium using:
%
%   KD(Fe-Mg) = (Mg_liq * Fe_cpx) / (Mg_cpx * Fe_liq)
%
% with mean KD = 0.275 +/- 0.067. The approximate 3-sigma experimental
% interval is 0.105-0.488 (p. 1550). Values outside this interval are strong
% evidence that the selected Clinopyroxene and Liquid may not be an
% equilibrium pair.
%
% Pressure estimates from individual Clinopyroxene analyses may scatter
% because of phase heterogeneity and slow diffusion. Averaging several
% demonstrably equilibrium Clinopyroxene analyses can reduce the uncertainty
% of the mean pressure to less than approximately 1.5 kbar (p. 1549).
%
% Liquid components are calculated as cation fractions. H2O is excluded from
% cation-fraction normalization, and oxide wt% values are not renormalized to
% 100 before conversion (pp. 1545-1546). In this implementation, F and Cl
% are also excluded from cationTotal_liq and from NaN-input warnings because
% neither is used by Model A.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 850-1300 degreeC,
%   2) finite calculated pressure is outside 0-40 kbar,
%   3) finite calculated pressure exceeds 50 kbar,
%   4) Liquid SiO2 exceeds 71.3 wt%,
%   5) KD(Fe-Mg) is outside the experimental equilibrium interval,
%   6) a required calculation input contains NaN, or
%   7) a calculated pressure is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Cpx : table
%
% The FIRST column of the Cpx table is treated as an identifier ("data
% code") displayed in the selection dialog. Oxide columns may be named as
% either Oxide or Oxide_value, with spaces, underscores, and hyphens ignored
% during matching.
%
% Required Clinopyroxene oxide variables:
%   SiO2, Al2O3, MgO, CaO, Na2O, and either FeO or FeOt
%
% Optional Clinopyroxene oxide variables:
%   TiO2, MnO, K2O, Cr2O3
% Missing optional variables are treated as zero. If an existing variable
% contains NaN, NaN is retained and propagated; it is never replaced by
% zero.
%
% Liquid data are loaded by liquid.readLiquidExcel. Required Liquid oxide
% variables are:
%   SiO2, Al2O3, MgO, CaO, Na2O, and either FeO or FeOt
%
% Optional Liquid oxide variables contributing to cationTotal_liq are:
%   TiO2, MnO, K2O, V2O3, Cr2O3, NiO, P2O5, SO3, Fe2O3
%
% Missing optional Liquid variables are treated as zero. Existing NaN values
% are retained. H2O, F, and Cl do not contribute to cationTotal_liq. F and
% Cl are additionally excluded from NaN-input warnings.
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if it makes a logarithm or ratio
% undefined, the resulting NaN or Inf is retained and reported.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (Model A, Table 4, p. 1546)
%
%   P(kbar) =
%     -88.3
%     + 2.82e-3*T(K)*ln[Jd_cpx/(Na_liq*Al_liq*Si_liq^2)]
%     + 2.19e-2*T(K)
%     - 25.1*ln[Ca_liq*Si_liq]
%     + 7.03*MgPrime_liq
%     + 12.4*ln[Ca_liq]
%
% where:
%   Jd_cpx       = min(Na_cpx, AlVI_cpx)
%   MgPrime_liq  = Mg_liq/(Mg_liq + Fe_liq)
%   T            = temperature in K
%   P            = pressure in kbar
%
% Clinopyroxene cations are calculated on a six-oxygen basis. Component
% allocation follows pp. 1545-1547 of Putirka et al. (2003).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2003baro(rawdata_struct, T_degreeC)
%   results = Putirka2003baro(rawdata_struct, T_degreeC, ...
%       'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct containing the Cpx table
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Name-value option:
%   'LiquidRow'    : positive integer row number in the loaded Liquid table.
%                    Default is row 1.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Clinopyroxene-Liquid pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('Putirka2003baro requires (rawdata_struct, T_degreeC).');
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

T_degreeC = T_degreeC(:);

%% Options
ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOption = ip.Results.LiquidRow;

%% 1) Retrieve Cpx and Liquid datasets
% Extract the Cpx table, molar-weight information, and the external Liquid
% table. Source tables are not modified.
disp('=== Step 1: Preparing Cpx and Liquid datasets ===');

if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

dataset_cpx = rawdata_struct.Cpx;
MWinfo = liquid.getMolarWeights();
[dataset_liq, metaLiq] = liquid.readLiquidExcel();

if isempty(dataset_liq)
    error('Selected Liquid dataset is empty.');
end

if isempty(liquidRowOption)
    selectedIdx_liq = 1;
    if height(dataset_liq) > 1
        fprintf(2, ...
            ['WARNING: Liquid dataset contains %d rows. LiquidRow was not ' ...
             'specified, so row 1 will be used for all calculations.\n'], ...
            height(dataset_liq));
    end
else
    selectedIdx_liq = liquidRowOption;
    if selectedIdx_liq > height(dataset_liq)
        error('Requested LiquidRow (%d) exceeds Liquid rows (%d).', ...
            selectedIdx_liq, height(dataset_liq));
    end
end

selectedData_liq = dataset_liq(selectedIdx_liq, :);
disp(['Liquid selected: Row ' num2str(selectedIdx_liq)]);
disp('=== Preparing Cpx and Liquid datasets has been finished ===');

%% 2) Initialize output container and applicability limits
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

newExperimentT_min_degreeC = 850;
newExperimentT_max_degreeC = 1300;
recommendedP_min_kbar = 0;
recommendedP_max_kbar = 40;
highPressureCaution_kbar = 50;
maximumNewExperimentSiO2_wt = 71.3;
KD_mean = 0.275;
KD_sigma = 0.067;
KD_3sigma_min = 0.105;
KD_3sigma_max = 0.488;

temperatureOutsideNewExperimentRange = isfinite(T_degreeC) & ...
    (T_degreeC < newExperimentT_min_degreeC | ...
     T_degreeC > newExperimentT_max_degreeC);
temperatureWarningIssued = false;
applicationCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive Cpx selection + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Cpx) ===');

while true
    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Cpx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIdx_cpx);
    disp(['Cpx selected: ' char(string(selectedCode_cpx))]);

    disp('=== Step 4: Preparing selected Cpx-Liquid pair ===');

    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);
    cpx = prepareCpxRow(selectedData_cpx);
    liq = prepareLiquidRow(selectedData_liq);

    % List NaN values only for variables that contribute to Model A or to
    % cationTotal_liq. F, Cl, and H2O are intentionally excluded.
    nanInputNames = findNaNInputs(cpx, liq, T_degreeC);

    % Reject only Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(cpx, liq);

    disp('=== Step 5: Calculating the pressure ===');

    row = calcPressure(cpx, liq, T_degreeC, MWinfo, ...
        newExperimentT_min_degreeC, newExperimentT_max_degreeC, ...
        recommendedP_min_kbar, recommendedP_max_kbar, ...
        highPressureCaution_kbar, maximumNewExperimentSiO2_wt, ...
        KD_3sigma_min, KD_3sigma_max);

    % Repeat identifiers for all temperatures in the current calculation.
    nRows = height(row);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_cpx', 'dataRow_liq'}, 'Before', 1);

    % Store one block per selected mineral-liquid pair. Expand the cell
    % buffer only when its current capacity has been exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    finitePressureValues = row.P_kbar(isfinite(row.P_kbar));
    if height(row) == 1
        disp([char(string(selectedCode_cpx)) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    elseif isempty(finitePressureValues)
        disp([char(string(selectedCode_cpx)) ...
            ': all calculated pressures are NaN or Inf']);
    else
        disp([char(string(selectedCode_cpx)) ': finite P range = ' ...
            num2str(min(finitePressureValues)) ' to ' ...
            num2str(max(finitePressureValues)) ' kbar']);
    end

    % Print the principal application cautions once per function call.
    if ~applicationCautionIssued
        fprintf(2, ...
            ['CAUTION: Putirka et al. (2003) Model A requires an equilibrium ' ...
             'Clinopyroxene-Liquid pair (pp. 1550-1552). Some hydrous ' ...
             'experimental datasets yield substantial pressure ' ...
             'underestimates, including negative calculated pressures at ' ...
             '1-2 kbar (p. 1550). Check petrography, liquid selection, and ' ...
             'KD(Fe-Mg) before interpreting the result.\n']);
        applicationCautionIssued = true;
    end

    % Input temperature is common to all selected pairs, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideNewExperimentRange) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the 850-1300 degreeC ' ...
             'range of the new high-SiO2 and volatile-bearing experiments ' ...
             'reported by Putirka et al. (2003; pp. 1543-1544, 1547-1548). ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'input finite range = %.4g-%.4g degreeC. This is a warning ' ...
             'envelope, not a strict universal limit for the complete ' ...
             'multi-study regression dataset.\n'], ...
            sum(temperatureOutsideNewExperimentRange), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures fall outside the range for which
    % the independent tests show no systematic pressure offset.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideRecommended = finitePressure & ...
        (row.P_kbar < recommendedP_min_kbar | ...
         row.P_kbar > recommendedP_max_kbar);

    if any(pressureOutsideRecommended)
        pressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the recommended ' ...
             '0-40 kbar envelope of Putirka et al. (2003; pp. 1548-1549). ' ...
             '%d of %d finite pressure point(s) are outside the envelope; ' ...
             'calculated finite range = %.4g-%.4g kbar for Cpx %s and ' ...
             'Liquid row %d. Values were retained.\n'], ...
            sum(pressureOutsideRecommended), ...
            sum(finitePressure), ...
            min(pressureValues), ...
            max(pressureValues), ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq);
    end

    aboveHighPressureCaution = finitePressure & ...
        row.P_kbar > highPressureCaution_kbar;
    if any(aboveHighPressureCaution)
        fprintf(2, ...
            ['WARNING: %d calculated pressure point(s) exceed 50 kbar for ' ...
             'Cpx %s and Liquid row %d. Putirka et al. (2003) report ' ...
             'systematic pressure underestimation above approximately ' ...
             '50 kbar (pp. 1548-1549). Values were retained.\n'], ...
            sum(aboveHighPressureCaution), ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq);
    end

    % Warn when the selected liquid exceeds the maximum SiO2 content of the
    % new experimental glasses reported in this paper.
    liquidSiO2 = row.SiO2_liq(1);
    if isfinite(liquidSiO2) && liquidSiO2 > maximumNewExperimentSiO2_wt
        fprintf(2, ...
            ['WARNING: Liquid SiO2 = %.4g wt%% exceeds the maximum ' ...
             '71.3 wt%% SiO2 of the new experiments reported by Putirka ' ...
             'et al. (2003; pp. 1542-1544) for Cpx %s and Liquid row %d.\n'], ...
            liquidSiO2, char(string(selectedCode_cpx)), selectedIdx_liq);
    end

    % Report Fe-Mg equilibrium diagnostics. KD is composition dependent and
    % therefore identical for all temperatures of the selected pair.
    KD_value = row.KD_FeMg_cpx_liq(1);
    if isfinite(KD_value)
        if KD_value < KD_3sigma_min || KD_value > KD_3sigma_max
            fprintf(2, ...
                ['WARNING: KD(Fe-Mg) = %.4g is outside the approximate ' ...
                 '3-sigma experimental interval 0.105-0.488 of Putirka ' ...
                 'et al. (2003; p. 1550) for Cpx %s and Liquid row %d. ' ...
                 'The selected pair may not represent equilibrium.\n'], ...
                KD_value, char(string(selectedCode_cpx)), selectedIdx_liq);
        elseif abs(KD_value - KD_mean) > KD_sigma
            fprintf(2, ...
                ['CAUTION: KD(Fe-Mg) = %.4g lies outside the mean +/-1 ' ...
                 'sigma interval %.3f-%.3f but remains within the broad ' ...
                 'experimental interval 0.105-0.488 (Putirka et al., ' ...
                 '2003; p. 1550) for Cpx %s and Liquid row %d.\n'], ...
                KD_value, KD_mean - KD_sigma, KD_mean + KD_sigma, ...
                char(string(selectedCode_cpx)), selectedIdx_liq);
        end
    end

    % Warn when equation terms are zero or otherwise outside their strict
    % logarithmic/ratio domains. Values are not altered.
    reportEquationDomainWarnings(row, selectedCode_cpx, selectedIdx_liq);

    % List exact calculation inputs containing NaN. Temperature indices are
    % included for vector input. F and Cl are deliberately absent.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for Cpx %s ' ...
             'and Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; calculated ' ...
             'pressure values may remain NaN.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for Cpx ' ...
             '%s and Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative calculated pressure is retained as a diagnostic result. This
    % behavior is important because negative pressures are also reported for
    % some hydrous test datasets in the source paper.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for Cpx %s ' ...
             'and Liquid row %d (%d of %d points). The values were retained ' ...
             'for diagnostic purposes.\n'], ...
            char(string(selectedCode_cpx)), ...
            selectedIdx_liq, ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Cpx selection (same Liquid dataset)?', ...
        'Putirka2003baro', ...
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

results.Properties.UserData = struct('liquid', metaLiq);
disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(cpx, liq, T_degreeC)
% findNaNInputs
% Return the names of pressure-equation and cationTotal_liq inputs containing
% NaN. F, Cl, and H2O are intentionally excluded. NaN values are never
% replaced by zero.

maxNames = 32;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(transpose(nanTemperatureIndices)), ',');
    nanInputBuffer(nNanInputs) = "T_degreeC(indices=" + indexText + ")";
end

cpxFieldNames = { ...
    'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};
for i = 1:numel(cpxFieldNames)
    fieldName = cpxFieldNames{i};
    if isnan(cpx.(fieldName))
        nNanInputs = nNanInputs + 1;
        if strcmp(fieldName, 'FeO')
            nanInputBuffer(nNanInputs) = "Cpx.FeO/FeOt";
        else
            nanInputBuffer(nNanInputs) = "Cpx." + string(fieldName);
        end
    end
end

% These are exactly the liquid oxides included in cationTotal_liq. F, Cl,
% and H2O are excluded from both the total and this warning list.
liqFieldNames = { ...
    'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', ...
    'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', 'Fe2O3'};
for i = 1:numel(liqFieldNames)
    fieldName = liqFieldNames{i};
    if isnan(liq.(fieldName))
        nNanInputs = nNanInputs + 1;
        if strcmp(fieldName, 'FeO')
            nanInputBuffer(nNanInputs) = "Liquid.FeO/FeOt";
        else
            nanInputBuffer(nNanInputs) = "Liquid." + string(fieldName);
        end
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(cpx, liq)
% validateNonNegativeInputs
% Reject Inf and finite negative values in all variables used by Model A or
% cationTotal_liq. Zero and NaN are intentionally allowed and retained. F,
% Cl, and H2O are not calculation inputs and are not checked here.

maxNames = 32;
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

cpxFieldNames = { ...
    'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};
for i = 1:numel(cpxFieldNames)
    fieldName = cpxFieldNames{i};
    value = cpx.(fieldName);
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        if strcmp(fieldName, 'FeO')
            invalidInputBuffer(nInvalidInputs) = "Cpx.FeO/FeOt";
        else
            invalidInputBuffer(nInvalidInputs) = "Cpx." + string(fieldName);
        end
    end
end

liqFieldNames = { ...
    'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', ...
    'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', 'Fe2O3'};
for i = 1:numel(liqFieldNames)
    fieldName = liqFieldNames{i};
    value = liq.(fieldName);
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        if strcmp(fieldName, 'FeO')
            invalidInputBuffer(nInvalidInputs) = "Liquid.FeO/FeOt";
        else
            invalidInputBuffer(nInvalidInputs) = "Liquid." + string(fieldName);
        end
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['Putirka2003baro: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(cpx, liq, T_degreeC, MWinfo, ...
        T_min_degreeC, T_max_degreeC, P_min_kbar, P_max_kbar, ...
        highPressureCaution_kbar, maximumSiO2_wt, KD_min, KD_max)
% calcPressure
% Compute pressure for one Clinopyroxene row and one Liquid row at one or
% more input temperatures. NaN values propagate through the calculation.
%
% Inputs:
%   cpx       : scalar structure of Clinopyroxene oxide wt% values
%   liq       : scalar structure of Liquid oxide wt% values
%   T_degreeC : scalar or vector temperature in degreeC
%   MWinfo    : molar-weight and cation-number structure
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% ----- Clinopyroxene formula and components on a six-oxygen basis -----
molProp_cpx = struct();
molProp_cpx.SiO2 = cpx.SiO2 ./ MWinfo.MW.SiO2;
molProp_cpx.TiO2 = cpx.TiO2 ./ MWinfo.MW.TiO2;
molProp_cpx.Al2O3 = cpx.Al2O3 ./ MWinfo.MW.Al2O3;
molProp_cpx.FeO = cpx.FeO ./ MWinfo.MW.FeO;
molProp_cpx.MnO = cpx.MnO ./ MWinfo.MW.MnO;
molProp_cpx.MgO = cpx.MgO ./ MWinfo.MW.MgO;
molProp_cpx.CaO = cpx.CaO ./ MWinfo.MW.CaO;
molProp_cpx.Na2O = cpx.Na2O ./ MWinfo.MW.Na2O;
molProp_cpx.K2O = cpx.K2O ./ MWinfo.MW.K2O;
molProp_cpx.Cr2O3 = cpx.Cr2O3 ./ MWinfo.MW.Cr2O3;

oxygenSum_cpx = ...
    2 .* molProp_cpx.SiO2 + ...
    2 .* molProp_cpx.TiO2 + ...
    3 .* molProp_cpx.Al2O3 + ...
    molProp_cpx.FeO + ...
    molProp_cpx.MnO + ...
    molProp_cpx.MgO + ...
    molProp_cpx.CaO + ...
    molProp_cpx.Na2O + ...
    molProp_cpx.K2O + ...
    3 .* molProp_cpx.Cr2O3;

ORF_cpx = 6 ./ oxygenSum_cpx;

XSi_cpx_scalar = molProp_cpx.SiO2 .* ORF_cpx;
XTi_cpx_scalar = molProp_cpx.TiO2 .* ORF_cpx;
XAl_cpx_scalar = 2 .* molProp_cpx.Al2O3 .* ORF_cpx;
XFe_cpx_scalar = molProp_cpx.FeO .* ORF_cpx;
XMn_cpx_scalar = molProp_cpx.MnO .* ORF_cpx;
XMg_cpx_scalar = molProp_cpx.MgO .* ORF_cpx;
XCa_cpx_scalar = molProp_cpx.CaO .* ORF_cpx;
XNa_cpx_scalar = 2 .* molProp_cpx.Na2O .* ORF_cpx;
XK_cpx_scalar = 2 .* molProp_cpx.K2O .* ORF_cpx;
XCr_cpx_scalar = 2 .* molProp_cpx.Cr2O3 .* ORF_cpx;

cationSum_cpx_scalar = XSi_cpx_scalar + XTi_cpx_scalar + ...
    XAl_cpx_scalar + XFe_cpx_scalar + XMn_cpx_scalar + ...
    XMg_cpx_scalar + XCa_cpx_scalar + XNa_cpx_scalar + ...
    XK_cpx_scalar + XCr_cpx_scalar;

XAlIV_cpx_scalar = nonNegativeOrNaN(2 - XSi_cpx_scalar);
XAlVI_cpx_scalar = nonNegativeOrNaN(XAl_cpx_scalar - XAlIV_cpx_scalar);
XFe3_cpx_scalar = nonNegativeOrNaN(XNa_cpx_scalar + XAlIV_cpx_scalar - ...
    XAlVI_cpx_scalar - 2 .* XTi_cpx_scalar - XCr_cpx_scalar);
XJd_cpx_scalar = minimumOrNaN(XAlVI_cpx_scalar, XNa_cpx_scalar);
XCaTs_cpx_scalar = nonNegativeOrNaN(XAlVI_cpx_scalar - XJd_cpx_scalar);

if isnan(XAlIV_cpx_scalar) || isnan(XCaTs_cpx_scalar)
    XCaTi_cpx_scalar = NaN;
elseif XAlIV_cpx_scalar > XCaTs_cpx_scalar
    XCaTi_cpx_scalar = nonNegativeOrNaN( ...
        (XAlIV_cpx_scalar - XCaTs_cpx_scalar) ./ 2);
else
    XCaTi_cpx_scalar = 0;
end

XCrCaTs_cpx_scalar = nonNegativeOrNaN(XCr_cpx_scalar ./ 2);
XFm_cpx_scalar = XFe_cpx_scalar + XMg_cpx_scalar;

componentInputs = [XCa_cpx_scalar, XCaTs_cpx_scalar, XCaTi_cpx_scalar, ...
    XCrCaTs_cpx_scalar, XFm_cpx_scalar, XAlVI_cpx_scalar, XAlIV_cpx_scalar];
if any(isnan(componentInputs))
    XDiHd_cpx_scalar = NaN;
    XEnFs_cpx_scalar = NaN;
    XFmCaTs_cpx_scalar = NaN;
    XFmTi_cpx_scalar = NaN;
elseif XCa_cpx_scalar >= (XCaTs_cpx_scalar + XCaTi_cpx_scalar)
    XDiHd_cpx_scalar = nonNegativeOrNaN(XCa_cpx_scalar - ...
        XCaTi_cpx_scalar - XCaTs_cpx_scalar - XCrCaTs_cpx_scalar);
    XFmCaTs_cpx_scalar = 0;
    XFmTi_cpx_scalar = 0;
    XEnFs_cpx_scalar = nonNegativeOrNaN( ...
        (XFm_cpx_scalar - XDiHd_cpx_scalar) ./ 2);
else
    XDiHd_cpx_scalar = 0;
    XCaTs_cpx_scalar = XCa_cpx_scalar;
    VIAlex = nonNegativeOrNaN(XAlVI_cpx_scalar - XCaTs_cpx_scalar);
    XFmCaTs_cpx_scalar = nonNegativeOrNaN(VIAlex - XCaTs_cpx_scalar);
    XFmTi_cpx_scalar = nonNegativeOrNaN((XAlIV_cpx_scalar - ...
        XCaTs_cpx_scalar - XFmCaTs_cpx_scalar) ./ 2);
    XEnFs_cpx_scalar = nonNegativeOrNaN((XFm_cpx_scalar - ...
        XFmCaTs_cpx_scalar - XFmTi_cpx_scalar) ./ 2);
end

% ----- Liquid cation proportions -----
% F, Cl, and H2O are deliberately excluded from cationTotal_liq.
n_liq = struct();
n_liq.SiO2 = liq.SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
n_liq.TiO2 = liq.TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
n_liq.Al2O3 = liq.Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n_liq.FeO = liq.FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
n_liq.MnO = liq.MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
n_liq.MgO = liq.MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
n_liq.CaO = liq.CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
n_liq.Na2O = liq.Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
n_liq.K2O = liq.K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
n_liq.V2O3 = liq.V2O3 .* MWinfo.Cat.V2O3 ./ MWinfo.MW.V2O3;
n_liq.Cr2O3 = liq.Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
n_liq.NiO = liq.NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;
n_liq.P2O5 = liq.P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;
n_liq.SO3 = liq.SO3 .* MWinfo.Cat.SO3 ./ MWinfo.MW.SO3;
n_liq.Fe2O3 = liq.Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

% F and Cl proportions are calculated only for diagnostic output. They do
% not contribute to cationTotal_liq or to any Model A term.
n_liq.F = liq.F .* MWinfo.Cat.F ./ MWinfo.MW.F;
n_liq.Cl = liq.Cl .* MWinfo.Cat.Cl ./ MWinfo.MW.Cl;

cationTotal_liq_scalar = ...
    n_liq.SiO2 + n_liq.TiO2 + n_liq.Al2O3 + n_liq.FeO + ...
    n_liq.MnO + n_liq.MgO + n_liq.CaO + n_liq.Na2O + ...
    n_liq.K2O + n_liq.V2O3 + n_liq.Cr2O3 + n_liq.NiO + ...
    n_liq.P2O5 + n_liq.SO3 + n_liq.Fe2O3;

XSiO2_liq_scalar = n_liq.SiO2 ./ cationTotal_liq_scalar;
XTiO2_liq_scalar = n_liq.TiO2 ./ cationTotal_liq_scalar;
XAlO1_5_liq_scalar = n_liq.Al2O3 ./ cationTotal_liq_scalar;
XFeO_liq_scalar = n_liq.FeO ./ cationTotal_liq_scalar;
XMnO_liq_scalar = n_liq.MnO ./ cationTotal_liq_scalar;
XMgO_liq_scalar = n_liq.MgO ./ cationTotal_liq_scalar;
XCaO_liq_scalar = n_liq.CaO ./ cationTotal_liq_scalar;
XNaO0_5_liq_scalar = n_liq.Na2O ./ cationTotal_liq_scalar;
XKO0_5_liq_scalar = n_liq.K2O ./ cationTotal_liq_scalar;
XVO1_5_liq_scalar = n_liq.V2O3 ./ cationTotal_liq_scalar;
XCrO1_5_liq_scalar = n_liq.Cr2O3 ./ cationTotal_liq_scalar;
XNiO_liq_scalar = n_liq.NiO ./ cationTotal_liq_scalar;
XPO2_5_liq_scalar = n_liq.P2O5 ./ cationTotal_liq_scalar;
XSO3_liq_scalar = n_liq.SO3 ./ cationTotal_liq_scalar;
XFeO1_5_liq_scalar = n_liq.Fe2O3 ./ cationTotal_liq_scalar;
XF_liq_scalar = n_liq.F ./ cationTotal_liq_scalar;
XCl_liq_scalar = n_liq.Cl ./ cationTotal_liq_scalar;

XNa_liq_scalar = XNaO0_5_liq_scalar;
XAl_liq_scalar = XAlO1_5_liq_scalar;
XCa_liq_scalar = XCaO_liq_scalar;
XSi_liq_scalar = XSiO2_liq_scalar;

MgPrimeDenominator_liq_scalar = XMgO_liq_scalar + XFeO_liq_scalar;
MgPrime_liq_scalar = XMgO_liq_scalar ./ MgPrimeDenominator_liq_scalar;

% Fe-Mg exchange coefficient used as an equilibrium diagnostic.
KD_FeMg_cpx_liq_scalar = ...
    (XMgO_liq_scalar .* XFe_cpx_scalar) ./ ...
    (XMg_cpx_scalar .* XFeO_liq_scalar);

% ----- Model A terms and pressure -----
JdTerm_scalar = XJd_cpx_scalar ./ ...
    (XNa_liq_scalar .* XAl_liq_scalar .* XSi_liq_scalar.^2);
CaSiTerm_scalar = XCa_liq_scalar .* XSi_liq_scalar;

lnJdTerm_scalar = log(JdTerm_scalar);
lnCaSiTerm_scalar = log(CaSiTerm_scalar);
lnXCa_liq_scalar = log(XCa_liq_scalar);

% Expand all composition-dependent scalars to the temperature-vector length.
XSi_cpx = repmat(XSi_cpx_scalar, nT, 1);
XTi_cpx = repmat(XTi_cpx_scalar, nT, 1);
XAl_cpx = repmat(XAl_cpx_scalar, nT, 1);
XFe_cpx = repmat(XFe_cpx_scalar, nT, 1);
XMn_cpx = repmat(XMn_cpx_scalar, nT, 1);
XMg_cpx = repmat(XMg_cpx_scalar, nT, 1);
XCa_cpx = repmat(XCa_cpx_scalar, nT, 1);
XNa_cpx = repmat(XNa_cpx_scalar, nT, 1);
XK_cpx = repmat(XK_cpx_scalar, nT, 1);
XCr_cpx = repmat(XCr_cpx_scalar, nT, 1);
cationSum_cpx = repmat(cationSum_cpx_scalar, nT, 1);
oxygenSum_cpx_output = repmat(oxygenSum_cpx, nT, 1);

XAlIV_cpx = repmat(XAlIV_cpx_scalar, nT, 1);
XAlVI_cpx = repmat(XAlVI_cpx_scalar, nT, 1);
XFe3_cpx = repmat(XFe3_cpx_scalar, nT, 1);
XJd_cpx = repmat(XJd_cpx_scalar, nT, 1);
XCaTs_cpx = repmat(XCaTs_cpx_scalar, nT, 1);
XCaTi_cpx = repmat(XCaTi_cpx_scalar, nT, 1);
XCrCaTs_cpx = repmat(XCrCaTs_cpx_scalar, nT, 1);
XDiHd_cpx = repmat(XDiHd_cpx_scalar, nT, 1);
XEnFs_cpx = repmat(XEnFs_cpx_scalar, nT, 1);
XFmCaTs_cpx = repmat(XFmCaTs_cpx_scalar, nT, 1);
XFmTi_cpx = repmat(XFmTi_cpx_scalar, nT, 1);

cationTotal_liq = repmat(cationTotal_liq_scalar, nT, 1);
XSiO2_liq = repmat(XSiO2_liq_scalar, nT, 1);
XTiO2_liq = repmat(XTiO2_liq_scalar, nT, 1);
XAlO1_5_liq = repmat(XAlO1_5_liq_scalar, nT, 1);
XFeO_liq = repmat(XFeO_liq_scalar, nT, 1);
XMnO_liq = repmat(XMnO_liq_scalar, nT, 1);
XMgO_liq = repmat(XMgO_liq_scalar, nT, 1);
XCaO_liq = repmat(XCaO_liq_scalar, nT, 1);
XNaO0_5_liq = repmat(XNaO0_5_liq_scalar, nT, 1);
XKO0_5_liq = repmat(XKO0_5_liq_scalar, nT, 1);
XVO1_5_liq = repmat(XVO1_5_liq_scalar, nT, 1);
XCrO1_5_liq = repmat(XCrO1_5_liq_scalar, nT, 1);
XNiO_liq = repmat(XNiO_liq_scalar, nT, 1);
XPO2_5_liq = repmat(XPO2_5_liq_scalar, nT, 1);
XSO3_liq = repmat(XSO3_liq_scalar, nT, 1);
XFeO1_5_liq = repmat(XFeO1_5_liq_scalar, nT, 1);
XF_liq = repmat(XF_liq_scalar, nT, 1);
XCl_liq = repmat(XCl_liq_scalar, nT, 1);

XNa_liq = repmat(XNa_liq_scalar, nT, 1);
XAl_liq = repmat(XAl_liq_scalar, nT, 1);
XCa_liq = repmat(XCa_liq_scalar, nT, 1);
XSi_liq = repmat(XSi_liq_scalar, nT, 1);
MgPrimeDenominator_liq = repmat(MgPrimeDenominator_liq_scalar, nT, 1);
MgPrime_liq = repmat(MgPrime_liq_scalar, nT, 1);
KD_FeMg_cpx_liq = repmat(KD_FeMg_cpx_liq_scalar, nT, 1);

JdTerm = repmat(JdTerm_scalar, nT, 1);
CaSiTerm = repmat(CaSiTerm_scalar, nT, 1);
lnJdTerm = repmat(lnJdTerm_scalar, nT, 1);
lnCaSiTerm = repmat(lnCaSiTerm_scalar, nT, 1);
lnXCa_liq = repmat(lnXCa_liq_scalar, nT, 1);

P_kbar = ...
    -88.3 ...
    + 2.82e-3 .* T_K .* lnJdTerm ...
    + 2.19e-2 .* T_K ...
    - 25.1 .* lnCaSiTerm ...
    + 7.03 .* MgPrime_liq ...
    + 12.4 .* lnXCa_liq;
P_GPa = P_kbar ./ 10;

% Diagnostic applicability flags. They do not replace petrographic and
% equilibrium assessment.
isWithinEquationDomain = ...
    isfinite(T_K) & T_K > 0 & ...
    isfinite(XJd_cpx) & XJd_cpx > 0 & ...
    isfinite(XNa_liq) & XNa_liq > 0 & ...
    isfinite(XAl_liq) & XAl_liq > 0 & ...
    isfinite(XSi_liq) & XSi_liq > 0 & ...
    isfinite(XCa_liq) & XCa_liq > 0 & ...
    isfinite(MgPrimeDenominator_liq) & MgPrimeDenominator_liq > 0;

isWithinNewExperimentTRange = isfinite(T_degreeC) & ...
    T_degreeC >= T_min_degreeC & T_degreeC <= T_max_degreeC;
isWithinRecommendedPRange = isfinite(P_kbar) & ...
    P_kbar >= P_min_kbar & P_kbar <= P_max_kbar;
isAboveHighPressureCaution = isfinite(P_kbar) & ...
    P_kbar > highPressureCaution_kbar;
isWithinNewExperimentSiO2Maximum = isfinite(liq.SiO2) & ...
    liq.SiO2 <= maximumSiO2_wt;
isWithinKDExperimentalRange = isfinite(KD_FeMg_cpx_liq) & ...
    KD_FeMg_cpx_liq >= KD_min & KD_FeMg_cpx_liq <= KD_max;

% ----- Pack outputs using vectors of equal height -----
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.SiO2_cpx = repmat(cpx.SiO2, nT, 1);
row.TiO2_cpx = repmat(cpx.TiO2, nT, 1);
row.Al2O3_cpx = repmat(cpx.Al2O3, nT, 1);
row.FeO_cpx = repmat(cpx.FeO, nT, 1);
row.MnO_cpx = repmat(cpx.MnO, nT, 1);
row.MgO_cpx = repmat(cpx.MgO, nT, 1);
row.CaO_cpx = repmat(cpx.CaO, nT, 1);
row.Na2O_cpx = repmat(cpx.Na2O, nT, 1);
row.K2O_cpx = repmat(cpx.K2O, nT, 1);
row.Cr2O3_cpx = repmat(cpx.Cr2O3, nT, 1);

row.oxygenSum_cpx = oxygenSum_cpx_output;
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
row.XFmCaTs_cpx = XFmCaTs_cpx;
row.XFmTi_cpx = XFmTi_cpx;

row.SiO2_liq = repmat(liq.SiO2, nT, 1);
row.TiO2_liq = repmat(liq.TiO2, nT, 1);
row.Al2O3_liq = repmat(liq.Al2O3, nT, 1);
row.FeO_liq = repmat(liq.FeO, nT, 1);
row.MnO_liq = repmat(liq.MnO, nT, 1);
row.MgO_liq = repmat(liq.MgO, nT, 1);
row.CaO_liq = repmat(liq.CaO, nT, 1);
row.Na2O_liq = repmat(liq.Na2O, nT, 1);
row.K2O_liq = repmat(liq.K2O, nT, 1);
row.V2O3_liq = repmat(liq.V2O3, nT, 1);
row.Cr2O3_liq = repmat(liq.Cr2O3, nT, 1);
row.NiO_liq = repmat(liq.NiO, nT, 1);
row.P2O5_liq = repmat(liq.P2O5, nT, 1);
row.SO3_liq = repmat(liq.SO3, nT, 1);
row.Fe2O3_liq = repmat(liq.Fe2O3, nT, 1);
row.H2O_liq = repmat(liq.H2O, nT, 1);
row.F_liq = repmat(liq.F, nT, 1);
row.Cl_liq = repmat(liq.Cl, nT, 1);

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
row.XVO1_5_liq = XVO1_5_liq;
row.XCrO1_5_liq = XCrO1_5_liq;
row.XNiO_liq = XNiO_liq;
row.XPO2_5_liq = XPO2_5_liq;
row.XSO3_liq = XSO3_liq;
row.XFeO1_5_liq = XFeO1_5_liq;
row.XF_liq = XF_liq;
row.XCl_liq = XCl_liq;

row.XNa_liq = XNa_liq;
row.XAl_liq = XAl_liq;
row.XCa_liq = XCa_liq;
row.XSi_liq = XSi_liq;
row.MgPrimeDenominator_liq = MgPrimeDenominator_liq;
row.MgPrime_liq = MgPrime_liq;
row.KD_FeMg_cpx_liq = KD_FeMg_cpx_liq;

row.JdTerm_2003 = JdTerm;
row.CaSiTerm_2003 = CaSiTerm;
row.lnJdTerm_2003 = lnJdTerm;
row.lnCaSiTerm_2003 = lnCaSiTerm;
row.lnXCa_liq_2003 = lnXCa_liq;

% Generic barometer output names required by the launchers.
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% Backward-compatible aliases retained from the original implementation.
row.P2003_kbar = P_kbar;
row.P2003_GPa = P_GPa;

row.P_regression_SEE_kbar = repmat(1.7, nT, 1);
row.P_testData_SEE_kbar = repmat(4.0, nT, 1);

row.isWithinEquationDomain = isWithinEquationDomain;
row.isWithinNewExperimentTRange = isWithinNewExperimentTRange;
row.isWithinRecommendedPRange = isWithinRecommendedPRange;
row.isAboveHighPressureCaution = isAboveHighPressureCaution;
row.isWithinNewExperimentSiO2Maximum = ...
    repmat(isWithinNewExperimentSiO2Maximum, nT, 1);
row.isWithinKDExperimentalRange = isWithinKDExperimentalRange;

end

function cpx = prepareCpxRow(data_cpx)
% prepareCpxRow
% Extract one-row Clinopyroxene oxide data. Required variables must exist.
% Missing optional variables are represented by zero. Existing NaN values
% are preserved and are not converted to zero.

if height(data_cpx) ~= 1
    error('Cpx input must be a 1-row table.');
end

cpx = struct();
cpx.SiO2 = getOxideRequired(data_cpx, 'SiO2', 'Cpx');
cpx.TiO2 = getOxideOptional(data_cpx, 'TiO2', 0);
cpx.Al2O3 = getOxideRequired(data_cpx, 'Al2O3', 'Cpx');
cpx.FeO = getFeOWithFallback(data_cpx, 'Cpx');
cpx.MnO = getOxideOptional(data_cpx, 'MnO', 0);
cpx.MgO = getOxideRequired(data_cpx, 'MgO', 'Cpx');
cpx.CaO = getOxideRequired(data_cpx, 'CaO', 'Cpx');
cpx.Na2O = getOxideRequired(data_cpx, 'Na2O', 'Cpx');
cpx.K2O = getOxideOptional(data_cpx, 'K2O', 0);
cpx.Cr2O3 = getOxideOptional(data_cpx, 'Cr2O3', 0);

end

function liq = prepareLiquidRow(data_liq)
% prepareLiquidRow
% Extract one-row Liquid oxide data. F, Cl, and H2O are retained for output
% but are not included in cationTotal_liq. Existing NaN values are preserved.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liq = struct();
liq.SiO2 = getOxideRequired(data_liq, 'SiO2', 'Liquid');
liq.TiO2 = getOxideOptional(data_liq, 'TiO2', 0);
liq.Al2O3 = getOxideRequired(data_liq, 'Al2O3', 'Liquid');
liq.FeO = getFeOWithFallback(data_liq, 'Liquid');
liq.MnO = getOxideOptional(data_liq, 'MnO', 0);
liq.MgO = getOxideRequired(data_liq, 'MgO', 'Liquid');
liq.CaO = getOxideRequired(data_liq, 'CaO', 'Liquid');
liq.Na2O = getOxideRequired(data_liq, 'Na2O', 'Liquid');
liq.K2O = getOxideOptional(data_liq, 'K2O', 0);
liq.V2O3 = getOxideOptional(data_liq, 'V2O3', 0);
liq.Cr2O3 = getOxideOptional(data_liq, 'Cr2O3', 0);
liq.NiO = getOxideOptional(data_liq, 'NiO', 0);
liq.P2O5 = getOxideOptional(data_liq, 'P2O5', 0);
liq.SO3 = getOxideOptional(data_liq, 'SO3', 0);
liq.Fe2O3 = getOxideOptional(data_liq, 'Fe2O3', 0);
liq.H2O = getOxideOptional(data_liq, 'H2O', 0);
liq.F = getOxideOptional(data_liq, 'F', 0);
liq.Cl = getOxideOptional(data_liq, 'Cl', 0);

end

function reportEquationDomainWarnings(row, selectedCode_cpx, selectedIdx_liq)
% reportEquationDomainWarnings
% Print non-stopping warnings for zero or non-finite equation terms. Values
% are retained exactly as calculated.

labels = { ...
    'XJd_cpx', row.XJd_cpx(1); ...
    'XNa_liq', row.XNa_liq(1); ...
    'XAl_liq', row.XAl_liq(1); ...
    'XSi_liq', row.XSi_liq(1); ...
    'XCa_liq', row.XCa_liq(1); ...
    'MgPrimeDenominator_liq', row.MgPrimeDenominator_liq(1)};

for i = 1:size(labels, 1)
    label = labels{i, 1};
    value = labels{i, 2};
    if isfinite(value) && value <= 0
        fprintf(2, ...
            ['WARNING: %s = %.4g is not strictly positive for Cpx %s and ' ...
             'Liquid row %d. A logarithm or ratio in Putirka et al. (2003) ' ...
             'Model A may be undefined; NaN or Inf results are retained.\n'], ...
            label, value, char(string(selectedCode_cpx)), selectedIdx_liq);
    end
end

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Repeat optional Liquid identifiers to match the number of temperature rows.

variableNames = data_liq.Properties.VariableNames;
nRows = height(row);

if any(strcmp(variableNames, 'Index'))
    indexValue = toScalarDoublePreserveNaN(data_liq.('Index'));
    row.liq_Index = repmat(indexValue, nRows, 1);
end
if any(strcmp(variableNames, 'Experiment'))
    experimentRaw = data_liq.('Experiment');
    experimentValue = string(experimentRaw(1));
    row.liq_Experiment = repmat(experimentValue, nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    citationRaw = data_liq.('Citation');
    citationValue = string(citationRaw(1));
    row.liq_Citation = repmat(citationValue, nRows, 1);
end

end

function value = getOxideRequired(data_tbl, oxide, datasetLabel)
% getOxideRequired
% Retrieve a required oxide variable. A present NaN is returned unchanged.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    error('%s table must contain variable: %s', datasetLabel, oxide);
end

value = toScalarDoublePreserveNaN(data_tbl.(columnName));

end

function value = getOxideOptional(data_tbl, oxide, missingDefault)
% getOxideOptional
% Retrieve an optional oxide variable. A missing column uses missingDefault;
% a present NaN remains NaN.

columnName = findOxideColumn(data_tbl.Properties.VariableNames, oxide);
if isempty(columnName)
    value = missingDefault;
else
    value = toScalarDoublePreserveNaN(data_tbl.(columnName));
end

end

function value = getFeOWithFallback(data_tbl, datasetLabel)
% getFeOWithFallback
% Use finite FeO when available; otherwise use FeOt. If both are present but
% non-finite, retain NaN. At least one of the two columns must exist.

FeO_name = findOxideColumn(data_tbl.Properties.VariableNames, 'FeO');
FeOt_name = findOxideColumn(data_tbl.Properties.VariableNames, 'FeOt');

if isempty(FeO_name) && isempty(FeOt_name)
    error('%s table must contain either FeO or FeOt.', datasetLabel);
end

FeO_value = NaN;
FeOt_value = NaN;
if ~isempty(FeO_name)
    FeO_value = toScalarDoublePreserveNaN(data_tbl.(FeO_name));
end
if ~isempty(FeOt_name)
    FeOt_value = toScalarDoublePreserveNaN(data_tbl.(FeOt_name));
end

if isfinite(FeO_value)
    value = FeO_value;
elseif isfinite(FeOt_value)
    value = FeOt_value;
elseif isinf(FeO_value)
    value = FeO_value;
elseif isinf(FeOt_value)
    value = FeOt_value;
else
    value = NaN;
end

end

function columnName = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match Oxide or Oxide_value while ignoring spaces, underscores, and hyphens.

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
        return
    end
end

end

function canonical = canonicalizeName(inputName)
% canonicalizeName
% Convert a variable name to the comparison form used by findOxideColumn.

canonical = lower(char(string(inputName)));
canonical = strrep(canonical, ' ', '');
canonical = strrep(canonical, '_', '');
canonical = strrep(canonical, '-', '');

end

function value = toScalarDoublePreserveNaN(raw)
% toScalarDoublePreserveNaN
% Convert a one-row table value to scalar double. Missing or unparsable
% existing values become NaN and are not replaced by zero.

value = NaN;

if isempty(raw)
    return
end

if isnumeric(raw) || islogical(raw)
    value = double(raw(1));
    return
end

if isstring(raw)
    if ismissing(raw(1))
        return
    end
    value = str2double(raw(1));
    return
end

if ischar(raw)
    value = str2double(string(raw));
    return
end

if iscell(raw)
    if isempty(raw{1})
        return
    end
    firstValue = raw{1};
    if isnumeric(firstValue) || islogical(firstValue)
        value = double(firstValue(1));
    elseif isstring(firstValue) || ischar(firstValue)
        value = str2double(string(firstValue));
    end
end

end

function value = nonNegativeOrNaN(value)
% nonNegativeOrNaN
% Clamp a finite calculated component to zero while preserving NaN.

if isfinite(value)
    value = max(value, 0);
elseif isnan(value)
    value = NaN;
end

end

function value = minimumOrNaN(value1, value2)
% minimumOrNaN
% Return min(value1,value2) when both are finite; otherwise preserve NaN.

if isfinite(value1) && isfinite(value2)
    value = min(value1, value2);
else
    value = NaN;
end

end
