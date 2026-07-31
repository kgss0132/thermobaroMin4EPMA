function results = Putirka2016baro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Amphibole/Putirka2016baro.m
% Tested with MATLAB R2024b
%
% Amphibole-Liquid empirical barometers
% Putirka, K.D. (2016)
% American Mineralogist, 101, 841-858
% DOI: https://doi.org/10.2138/am-2016-5506
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Amphibole analysis with one Liquid
% analysis and calculates the tentative Amphibole-Liquid barometers of
% Putirka (2016), Equations 7a, 7b, and 7c.
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. Equations 7a-7c do not contain a temperature term;
% consequently, the same composition-dependent pressure is repeated for
% every input temperature, while temperature-dependent applicability flags
% are evaluated separately for each temperature value.
%
% Equation 7a is retained as the primary launcher-compatible pressure
% output, P_kbar. Equation 7b is retained for comparison and may perform
% better at low temperature and above approximately 10 kbar. Equation 7c is
% an anhydrous formulation and is retained only as a crude check.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% IMPORTANT: Putirka (2016) describes Equations 7a-7c as TENTATIVE
% Amphibole-Liquid barometers. Individual pressure estimates have large
% uncertainties and should not be treated as precise crystallization
% pressures. Multiple equilibrium Amphibole-Liquid pairs should be used,
% and multiple pressure estimates should be averaged. Independent checks
% using Clinopyroxene or other phase-equilibrium constraints are strongly
% recommended (abstract, p. 841; discussion, pp. 851-855).
%
% Practical pressure range:
%
%   Best-supported interval : approximately 1-8 kbar
%   8-10 kbar              : use cautiously
%   >10 kbar               : Equation 7b or Ridolfi and Renzulli (2011)
%                            model 1d may be preferable
%
% Putirka (2016) shows that the 1-8 kbar interval can be distinguished to
% approximately +/-1 kbar only when multiple estimates are averaged
% (pp. 852-854; Fig. 7). Individual estimates have approximately:
%
%   Calibration data : +/-2 kbar
%   Independent tests: +/-3.6 to 4.3 kbar, approximately +/-4 kbar
%
% The broad experimental database extends beyond 8 kbar, but it is not a
% universal high-pressure validation range. This implementation therefore
% uses 1-8 kbar as a non-stopping practical-range warning and issues a
% stronger warning above 10 kbar.
%
% Temperature:
% Equations 7a-7c are temperature independent, but Putirka (2016) confirms
% that Amphibole barometers are more precise and accurate under the
% Anderson and Smith (1995) restrictions:
%
%   T <= 800 degreeC
%   Fe#Amp = FeT_amp / (FeT_amp + Mg_amp) < 0.65
%
% These restrictions are discussed on pp. 841-842 and 851-852. The
% experimental Amphibole database spans broadly approximately 650-1175
% degreeC (Figs. 1-2, pp. 842-844). In this implementation, 650-1175
% degreeC is treated only as a broad experimental-temperature envelope,
% whereas T <= 800 degreeC is treated as the preferred reliability filter.
%
% Amphibole composition:
% The barometers should probably be applied only to CALCIC IGNEOUS
% AMPHIBOLES. Equations 7a and 7c return extreme false pressures for
% Mg-Fe-Mn-Li-group Amphiboles, and sodic-calcic Amphiboles were generally
% excluded from calibration and testing (Methods, pp. 844-846; discussion,
% p. 852). This function calculates a numerical calcic-group screening flag,
% but the user must independently confirm Amphibole classification and
% petrographic origin.
%
% Amphibole-Liquid equilibrium:
% The selected Amphibole and Liquid must represent the same magmatic
% equilibrium. Putirka (2016) proposes the Fe-Mg exchange coefficient:
%
%   KD_FeMg_Amp_Liq =
%     (FeT_amp / Mg_amp) / (XFeOt_liq / XMgO_liq)
%
% with a mean and standard deviation of 0.28 +/- 0.11 and a broad 10th-90th
% percentile interval of 0.13-0.41 (Eq. 2, pp. 848-850). This is a useful
% but imperfect test because Fe-Mg equilibrium does not guarantee Al or Na
% equilibrium. Zoning, inherited crystals, magma mixing, alteration, and
% mismatched crystal-liquid pairs may produce misleading pressures.
%
% Liquid H2O:
% Equations 7a and 7b use HYDROUS oxide mole fractions. Pressure estimates
% increase by approximately 0.4 kbar per 1 wt% increase in assumed H2O
% (pp. 850-855). Missing H2O is therefore retained as NaN by default and is
% never silently replaced by zero. The optional DefaultH2O_wt parameter may
% be supplied explicitly, but use of that substituted value is reported by
% fprintf.
%
% Equation choice:
%   Eq. 7a : preferred general formulation; hydrous mole fractions
%   Eq. 7b : comparison formulation; may perform better at low T and >10 kbar
%   Eq. 7c : anhydrous formulation; NOT recommended except as a crude check
%
% Equation 7c is explicitly described as less precise and not recommended
% except as a crude check (p. 852). Equations 7a and 7b should be compared
% with one another and with independent pressure estimates (pp. 852-855).
%
% Liquid and Amphibole normalization:
% Amphibole cations are used on a 23-oxygen basis. Total Fe is treated as
% FeOt for thermometry and barometry. Equations 7a and 7b use hydrous oxide
% mole fractions, whereas Equation 7c uses anhydrous oxide mole fractions
% (Tables 4b-4c, pp. 846-847; Equations 7a-7c, pp. 850-852).
%
% F and Cl are NOT included in the Liquid mole-fraction totals
% (cationTotal_liq_hyd and cationTotal_liq_anhyd) and are NOT included in
% the NaN-input warning list. This follows the oxide-component definition
% used in Table 4c and Equations 7a-7c.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input T is outside the broad 650-1175 degreeC envelope,
%   2) finite input T is above the preferred 800 degreeC filter,
%   3) finite calculated pressure is outside the practical 1-8 kbar range,
%   4) finite calculated pressure is above 10 kbar,
%   5) Amphibole Fe# is not below 0.65,
%   6) the Amphibole does not pass the numerical calcic-group screen,
%   7) Amphibole cation sum is outside 15-16 apfu,
%   8) KD_FeMg_Amp_Liq is outside 0.13-0.41,
%   9) an equation input contains NaN,
%  10) an equation returns NaN or Inf, or
%  11) a negative finite pressure is calculated.
%
% All warnings are diagnostic and do not stop the calculation. Passing all
% numerical checks does not establish complete equilibrium or applicability.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole : table
%
% The FIRST Amphibole-table column is treated as the data identifier shown
% in the selection dialog. The function searches for the following 23-O
% cation variables or compatible aliases:
%
%   Si_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu           % total Fe
%   Mn_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%   Cr_cation_apfu
%   Ni_cation_apfu
%
% Missing Amphibole variables are represented by NaN. Equations 7a and 7c
% require Al, Na, and K; Equation 7b requires Al. FeT, Mg, Ca, and Na are
% used by numerical applicability and equilibrium checks.
%
% Liquid analyses are read using liquid.readLiquidExcel(). The selected
% Liquid table should provide the following oxide wt% values:
%
%   SiO2, TiO2, Al2O3, FeOt or FeO/Fe2O3, MnO, MgO, CaO,
%   Na2O, K2O, P2O5, and H2O
%
% Missing or explicit NaN values remain NaN. No calculation input is
% silently replaced by zero. F and Cl, when present, are ignored in Liquid
% totals and NaN diagnostics.
%
% Finite calculation inputs must be non-negative. NaN is allowed, retained,
% propagated, and reported. Inf and finite negative values are rejected.
% Zero is retained; if it makes a ratio or logarithm invalid, the associated
% equation result remains NaN and a warning is printed.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATIONS
%
% Equation 7a, hydrous mole fractions:
%
%   P7a(kbar) =
%     -30.93 - 42.74*ln(DAl_hyd) - 42.16*ln(XAl2O3_hyd)
%     + 633*XP2O5_hyd + 12.64*XH2O_hyd
%     + 24.57*Al_amp + 18.6*K_amp + 4.0*ln(DNa_hyd)
%
% Equation 7b, hydrous mole fractions:
%
%   P7b(kbar) =
%     -64.79 - 6.064*ln(DAl_hyd) + 61.75*XSiO2_hyd
%     + 682*XP2O5_hyd - 101.9*XCaO_hyd + 7.85*Al_amp
%     - 46.46*ln(XSiO2_hyd)
%     - 4.81*ln(XNa2O_hyd + XK2O_hyd)
%
% Equation 7c, anhydrous mole fractions:
%
%   P7c(kbar) =
%     -45.5 - 46.3*ln(DAl_anhyd) - 51.1*ln(XAl2O3_anhyd)
%     + 439*XP2O5_anhyd + 26.6*Al_amp + 22.5*K_amp
%     + 5.23*ln(DNa_anhyd)
%
%   DAl = Al_amp / XAl2O3_liq
%   DNa = Na_amp / XNa2O_liq
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2016baro(rawdata_struct, T_degreeC)
%   results = Putirka2016baro(..., 'UseFirstLiquidRow', false)
%   results = Putirka2016baro(..., 'DefaultH2O_wt', value)
%
% Inputs:
%   rawdata_struct   : struct containing an Amphibole table
%   T_degreeC        : non-negative numeric scalar or vector; NaN allowed
%   UseFirstLiquidRow: logical scalar, default true
%   DefaultH2O_wt    : explicit fallback H2O wt%; default NaN
%
% Output:
%   results : table containing one row per temperature value for each
%             selected Amphibole-Liquid pair
%

%% Input validation
if nargin < 2
    error('Putirka2016baro requires (rawdata_struct, T_degreeC).');
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
if ~isfield(rawdata_struct, 'Amphibole') || ...
        ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end

T_degreeC = T_degreeC(:);

%% Options
ip = inputParser;
ip.addParameter('UseFirstLiquidRow', true, ...
    @(x) islogical(x) && isscalar(x));
ip.addParameter('DefaultH2O_wt', NaN, ...
    @(x) isnumeric(x) && isscalar(x) && ...
    (isnan(x) || (isfinite(x) && x >= 0)));
ip.parse(varargin{:});

useFirstLiquidRow = ip.Results.UseFirstLiquidRow;
defaultH2O_wt = ip.Results.DefaultH2O_wt;

%% 1) Retrieve Amphibole and Liquid datasets
disp('=== Step 1: Preparing cation and liquid datasets ===');

dataset_amp = rawdata_struct.Amphibole;

MWinfo = liquid.getMolarWeights();
[dataset_liq, metaLiq] = liquid.readLiquidExcel();

if isempty(dataset_liq) || ~istable(dataset_liq)
    error('Selected Liquid dataset is empty or is not a table.');
end

disp('=== Preparing cation and liquid datasets has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result block in a cell buffer and concatenate only
% once after the interactive loop.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

experimentalT_min_degreeC = 650;
experimentalT_max_degreeC = 1175;
preferredT_max_degreeC = 800;
preferredFeNumber_max = 0.65;
practicalP_min_kbar = 1;
practicalP_max_kbar = 8;
highP_threshold_kbar = 10;

temperatureOutsideExperimental = isfinite(T_degreeC) & ...
    (T_degreeC < experimentalT_min_degreeC | ...
     T_degreeC > experimentalT_max_degreeC);
temperatureAbovePreferred = isfinite(T_degreeC) & ...
    T_degreeC > preferredT_max_degreeC;

temperatureEnvelopeWarningIssued = false;
temperatureFilterWarningIssued = false;
applicationCautionIssued = false;
equationCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % ----- Amphibole selection -----
    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, ok] = listdlg( ...
        'PromptString', ...
        'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_amp)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

    % ----- Liquid selection -----
    disp('=== Step 4: Selecting a data code from the list (Liquid) ===');

    if useFirstLiquidRow
        selectedIdx_liq = 1;
        if height(dataset_liq) > 1
            fprintf(2, ['WARNING: The Liquid dataset contains %d rows. ' ...
                'UseFirstLiquidRow=true, so Liquid row 1 is used.\n'], ...
                height(dataset_liq));
        end
        disp('Liquid selected: row 1');
    else
        liquidList = buildLiquidList(dataset_liq);

        [selectedIdx_liq, okLiq] = listdlg( ...
            'PromptString', ...
            'Please select the Liquid data you would like to use:', ...
            'SelectionMode', 'single', ...
            'ListString', liquidList, ...
            'ListSize', [520 360]);

        if ~okLiq || isempty(selectedIdx_liq)
            disp('Liquid selection canceled');
            break;
        end

        disp(['Liquid selected: row ' num2str(selectedIdx_liq)]);
    end

    % ----- Calculation -----
    disp('=== Step 5: Calculating the pressure ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);
    selectedData_liq = dataset_liq(selectedIdx_liq, :);

    amp = prepareAmphiboleRow(selectedData_amp);
    liq = prepareLiquidMoleFractions( ...
        selectedData_liq, MWinfo, defaultH2O_wt);

    nanInputNames = findNaNInputs(amp, liq, T_degreeC);

    row = calcPressure(amp, liq, T_degreeC);

    % Repeat identifiers for all temperature rows.
    row.dataCode_amp = ...
        repmat(string(selectedCode_amp), height(row), 1);
    row.dataRow_liq = ...
        repmat(selectedIdx_liq, height(row), 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_amp','dataRow_liq'}, 'Before', 1);

    % Store one result block per selected pair.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = ...
            [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo pressures. Equations 7a-7c are temperature independent, so their
    % values are repeated across all temperature rows.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    fprintf('%s & Liquid row %d: Eq.7a = %.6g kbar, Eq.7b = %.6g kbar, Eq.7c = %.6g kbar', ...
        char(string(selectedCode_amp)), selectedIdx_liq, ...
        row.P_Putirka2016_7a_kbar(1), ...
        row.P_Putirka2016_7b_kbar(1), ...
        row.P_Putirka2016_7c_kbar(1));

    finiteTemperatureValues = row.T_degreeC(isfinite(row.T_degreeC));
    if isempty(finiteTemperatureValues)
        fprintf('; all %d input temperature value(s) are NaN\n', height(row));
    elseif height(row) == 1
        fprintf(' at T = %.6g degreeC\n', row.T_degreeC);
    else
        fprintf(' for finite T = %.6g to %.6g degreeC (%d rows)\n', ...
            min(finiteTemperatureValues), ...
            max(finiteTemperatureValues), ...
            height(row));
    end

    % Print limitations that cannot be established automatically once.
    if ~applicationCautionIssued
        fprintf(2, ['CAUTION: Putirka (2016) treats Equations 7a-7c as ' ...
            'tentative Amphibole-Liquid barometers. Multiple equilibrium ' ...
            'pairs should be averaged, and independent pressure checks are ' ...
            'required (p. 841; pp. 852-855). Numerical flags do not verify ' ...
            'petrographic equilibrium, zoning, inheritance, alteration, or ' ...
            'whether the selected Amphibole and Liquid record the same event.\n']);
        applicationCautionIssued = true;
    end

    if ~equationCautionIssued
        fprintf(2, ['CAUTION: Equation 7a is used as P_kbar. Equation 7b ' ...
            'should be compared with Equation 7a and may be preferable at ' ...
            'P > 10 kbar. Equation 7c is not recommended except as a crude ' ...
            'check (pp. 850-852, 854).\n']);
        equationCautionIssued = true;
    end

    % Temperature warnings are common to all pairs and printed once.
    if any(temperatureOutsideExperimental) && ...
            ~temperatureEnvelopeWarningIssued
        finiteT = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ['WARNING: Input temperature is outside the broad ' ...
            'approximately 650-1175 degreeC experimental Amphibole-data ' ...
            'envelope in Putirka (2016; Figs. 1-2, pp. 842-844). ' ...
            '%d of %d finite temperature point(s) are outside; finite ' ...
            'input range = %.6g-%.6g degreeC.\n'], ...
            sum(temperatureOutsideExperimental), ...
            numel(finiteT), min(finiteT), max(finiteT));
        temperatureEnvelopeWarningIssued = true;
    end

    if any(temperatureAbovePreferred) && ~temperatureFilterWarningIssued
        finiteT = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ['WARNING: Input temperature exceeds the preferred ' ...
            'T <= 800 degreeC reliability filter confirmed by Putirka ' ...
            '(2016; pp. 841-842, 851-852). %d of %d finite temperature ' ...
            'point(s) exceed 800 degreeC; finite input range = ' ...
            '%.6g-%.6g degreeC. Equations 7a-7c themselves do not contain ' ...
            'a temperature term.\n'], ...
            sum(temperatureAbovePreferred), ...
            numel(finiteT), min(finiteT), max(finiteT));
        temperatureFilterWarningIssued = true;
    end

    % Amphibole compositional and equilibrium diagnostics.
    if isfinite(row.FeNumber_amp(1)) && ...
            row.FeNumber_amp(1) >= preferredFeNumber_max
        fprintf(2, ['WARNING: Amphibole Fe# = %.6g is not below the ' ...
            'preferred Fe# < 0.65 filter of Putirka (2016; pp. 841-842, ' ...
            '851-852) for %s & Liquid row %d.\n'], ...
            row.FeNumber_amp(1), ...
            char(string(selectedCode_amp)), selectedIdx_liq);
    end

    if ~row.isCalcicAmphibole_numeric(1)
        fprintf(2, ['WARNING: The selected Amphibole does not pass the ' ...
            'numerical calcic-group screen for %s & Liquid row %d. ' ...
            'Putirka (2016) states that the barometers should probably be ' ...
            'applied only to calcic Amphiboles (pp. 844-846, 852). ' ...
            'Ca_B = %.6g, Na_B = %.6g, (Ca+Na)_B = %.6g.\n'], ...
            char(string(selectedCode_amp)), selectedIdx_liq, ...
            row.Ca_B_amp(1), row.Na_B_amp(1), row.CaNa_B_amp(1));
    end

    if isfinite(row.SumCat_amp_23O(1)) && ...
            ~row.isWithinCationSumRange(1)
        fprintf(2, ['WARNING: Amphibole 23-O cation sum = %.6g is outside ' ...
            'the expected approximately 15-16 apfu interval shown by ' ...
            'Putirka (2016; Table 4b, p. 847) for %s & Liquid row %d.\n'], ...
            row.SumCat_amp_23O(1), ...
            char(string(selectedCode_amp)), selectedIdx_liq);
    end

    if isfinite(row.KD_FeMg_Amp_Liq(1)) && ...
            ~row.isWithinBroadKDRange(1)
        fprintf(2, ['WARNING: KD_FeMg_Amp_Liq = %.6g is outside the broad ' ...
            '0.13-0.41 equilibrium-screening interval of Putirka (2016; ' ...
            'Eqs. 1-2, pp. 848-850) for %s & Liquid row %d. Fe-Mg ' ...
            'agreement alone does not guarantee Al or Na equilibrium.\n'], ...
            row.KD_FeMg_Amp_Liq(1), ...
            char(string(selectedCode_amp)), selectedIdx_liq);
    end

    if liq.usedDefaultH2O
        fprintf(2, ['WARNING: Liquid H2O was missing or NaN and was replaced ' ...
            'by the explicitly supplied DefaultH2O_wt = %.6g wt%% for %s ' ...
            '& Liquid row %d. Putirka (2016) reports approximately ' ...
            '0.4 kbar pressure change per 1 wt%% H2O (pp. 850-855).\n'], ...
            liq.H2O_wt, ...
            char(string(selectedCode_amp)), selectedIdx_liq);
    end

    % Pressure-range and non-finite diagnostics for each equation.
    printPressureDiagnostics( ...
        row.P_Putirka2016_7a_kbar, 'Equation 7a', ...
        practicalP_min_kbar, practicalP_max_kbar, highP_threshold_kbar, ...
        selectedCode_amp, selectedIdx_liq);

    printPressureDiagnostics( ...
        row.P_Putirka2016_7b_kbar, 'Equation 7b', ...
        practicalP_min_kbar, practicalP_max_kbar, highP_threshold_kbar, ...
        selectedCode_amp, selectedIdx_liq);

    printPressureDiagnostics( ...
        row.P_Putirka2016_7c_kbar, 'Equation 7c', ...
        practicalP_min_kbar, practicalP_max_kbar, highP_threshold_kbar, ...
        selectedCode_amp, selectedIdx_liq);

    % List exact inputs containing NaN. F and Cl are intentionally absent.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the Putirka (2016) ' ...
            'barometer input(s) for %s & Liquid row %d: %s.\n' ...
            '         NaN values were retained and were not replaced by ' ...
            'zero; affected derived values and pressures remain NaN.\n'], ...
            char(string(selectedCode_amp)), selectedIdx_liq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Report invalid equation domains separately.
    if ~row.isWithinEquation7aDomain(1)
        fprintf(2, ['WARNING: Equation 7a has an invalid ratio/logarithm ' ...
            'domain for %s & Liquid row %d. DAl_hyd = %.6g, ' ...
            'XAl2O3_hyd = %.6g, DNa_hyd = %.6g, XNa2O_hyd = %.6g. ' ...
            'The Equation 7a pressure remains NaN.\n'], ...
            char(string(selectedCode_amp)), selectedIdx_liq, ...
            row.DAl_hyd(1), row.XAl2O3_liq_hyd(1), ...
            row.DNa_hyd(1), row.XNa2O_liq_hyd(1));
    end

    if ~row.isWithinEquation7bDomain(1)
        fprintf(2, ['WARNING: Equation 7b has an invalid ratio/logarithm ' ...
            'domain for %s & Liquid row %d. DAl_hyd = %.6g, ' ...
            'XSiO2_hyd = %.6g, XNa2O_hyd + XK2O_hyd = %.6g. ' ...
            'The Equation 7b pressure remains NaN.\n'], ...
            char(string(selectedCode_amp)), selectedIdx_liq, ...
            row.DAl_hyd(1), row.XSiO2_liq_hyd(1), ...
            row.XNaK_liq_hyd(1));
    end

    if ~row.isWithinEquation7cDomain(1)
        fprintf(2, ['WARNING: Equation 7c has an invalid ratio/logarithm ' ...
            'domain for %s & Liquid row %d. DAl_anhyd = %.6g, ' ...
            'XAl2O3_anhyd = %.6g, DNa_anhyd = %.6g, ' ...
            'XNa2O_anhyd = %.6g. The Equation 7c pressure remains NaN.\n'], ...
            char(string(selectedCode_amp)), selectedIdx_liq, ...
            row.DAl_anhyd(1), row.XAl2O3_liq_anhyd(1), ...
            row.DNa_anhyd(1), row.XNa2O_liq_anhyd(1));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection (same Liquid dataset)?', ...
        'Putirka2016baro', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered blocks once.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct('liquid', metaLiq);
disp('=== Putirka2016baro finished ===');

end

%% ---- local functions ----
function amp = prepareAmphiboleRow(data_amp)
% prepareAmphiboleRow
% Extract one Amphibole row using 23-O cation variables or compatible
% aliases. Missing values remain NaN. Inf and finite negative values are
% rejected.

if height(data_amp) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();

amp.Si = getFirstVariableOrNaN(data_amp, ...
    {'Si_cation_apfu','Si_23O','Si23O','Si_cation'}, 'Amphibole.Si');
amp.Ti = getFirstVariableOrNaN(data_amp, ...
    {'Ti_cation_apfu','Ti_23O','Ti23O','Ti_cation'}, 'Amphibole.Ti');
amp.Al = getFirstVariableOrNaN(data_amp, ...
    {'Al_cation_apfu','Al_23O','Al23O','Al_cation'}, 'Amphibole.Al');
amp.FeT = getFirstVariableOrNaN(data_amp, ...
    {'Fe_cation_apfu','FeT_23O','Fe_23O','Fe23O','Fe_cation'}, ...
    'Amphibole.FeT');
amp.Mn = getFirstVariableOrNaN(data_amp, ...
    {'Mn_cation_apfu','Mn_23O','Mn23O','Mn_cation'}, 'Amphibole.Mn');
amp.Mg = getFirstVariableOrNaN(data_amp, ...
    {'Mg_cation_apfu','Mg_23O','Mg23O','Mg_cation'}, 'Amphibole.Mg');
amp.Ca = getFirstVariableOrNaN(data_amp, ...
    {'Ca_cation_apfu','Ca_23O','Ca23O','Ca_cation'}, 'Amphibole.Ca');
amp.Na = getFirstVariableOrNaN(data_amp, ...
    {'Na_cation_apfu','Na_23O','Na23O','Na_cation'}, 'Amphibole.Na');
amp.K = getFirstVariableOrNaN(data_amp, ...
    {'K_cation_apfu','K_23O','K23O','K_cation'}, 'Amphibole.K');
amp.Cr = getFirstVariableOrNaN(data_amp, ...
    {'Cr_cation_apfu','Cr_23O','Cr23O','Cr_cation'}, 'Amphibole.Cr');
amp.Ni = getFirstVariableOrNaN(data_amp, ...
    {'Ni_cation_apfu','Ni_23O','Ni23O','Ni_cation'}, 'Amphibole.Ni');

cationVector = [amp.Si, amp.Ti, amp.Al, amp.FeT, amp.Mn, ...
    amp.Mg, amp.Ca, amp.Na, amp.K, amp.Cr, amp.Ni];

if all(isfinite(cationVector))
    amp.SumCat = sum(cationVector);
else
    amp.SumCat = NaN;
end

denFeMg = amp.FeT + amp.Mg;
if isfinite(denFeMg) && denFeMg > 0
    amp.FeNumber = amp.FeT ./ denFeMg;
else
    amp.FeNumber = NaN;
end

% Numerical Leake-style group screening used only as a diagnostic.
if isfinite(amp.Ca)
    amp.Ca_B = min(amp.Ca, 2);
else
    amp.Ca_B = NaN;
end

if isfinite(amp.Na) && isfinite(amp.Ca_B)
    amp.Na_B = min(amp.Na, max(0, 2 - amp.Ca_B));
else
    amp.Na_B = NaN;
end

amp.CaNa_B = amp.Ca_B + amp.Na_B;

amp.isCalcic_numeric = ...
    isfinite(amp.CaNa_B) && amp.CaNa_B >= 1.0 && ...
    isfinite(amp.Na_B) && amp.Na_B <= 0.5;

end

function liq = prepareLiquidMoleFractions( ...
        data_liq, MWinfo, defaultH2O_wt)
% prepareLiquidMoleFractions
% Extract Liquid oxide wt% values and calculate hydrous and anhydrous oxide
% mole fractions. Missing values remain NaN. F and Cl are intentionally
% excluded from all totals and NaN diagnostics.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liq = struct();

liq.SiO2_wt = getOxideOrNaN(data_liq, {'SiO2'}, 'Liq.SiO2');
liq.TiO2_wt = getOxideOrNaN(data_liq, {'TiO2'}, 'Liq.TiO2');
liq.Al2O3_wt = getOxideOrNaN(data_liq, {'Al2O3'}, 'Liq.Al2O3');
liq.MnO_wt = getOxideOrNaN(data_liq, {'MnO'}, 'Liq.MnO');
liq.MgO_wt = getOxideOrNaN(data_liq, {'MgO'}, 'Liq.MgO');
liq.CaO_wt = getOxideOrNaN(data_liq, {'CaO'}, 'Liq.CaO');
liq.Na2O_wt = getOxideOrNaN(data_liq, {'Na2O'}, 'Liq.Na2O');
liq.K2O_wt = getOxideOrNaN(data_liq, {'K2O'}, 'Liq.K2O');
liq.P2O5_wt = getOxideOrNaN(data_liq, {'P2O5'}, 'Liq.P2O5');

% Total Fe is treated as FeOt. Prefer a direct total-Fe-as-FeO column.
[FeOt_direct, hasFeOt] = getOxideWithPresence(data_liq, ...
    {'FeOt','FeOT','FeOtotal','FeO_Total','FeOtot'}, 'Liq.FeOt');
[FeO_value, hasFeO] = getOxideWithPresence(data_liq, ...
    {'FeO'}, 'Liq.FeO');
[Fe2O3_value, hasFe2O3] = getOxideWithPresence(data_liq, ...
    {'Fe2O3'}, 'Liq.Fe2O3');

if hasFeOt
    liq.FeOt_wt = FeOt_direct;
elseif hasFeO && hasFe2O3
    if isnan(FeO_value) || isnan(Fe2O3_value)
        liq.FeOt_wt = NaN;
    else
        liq.FeOt_wt = FeO_value + ...
            Fe2O3_value .* (2 .* MWinfo.MW.FeO ./ MWinfo.MW.Fe2O3);
    end
elseif hasFeO
    liq.FeOt_wt = FeO_value;
elseif hasFe2O3
    if isnan(Fe2O3_value)
        liq.FeOt_wt = NaN;
    else
        liq.FeOt_wt = ...
            Fe2O3_value .* (2 .* MWinfo.MW.FeO ./ MWinfo.MW.Fe2O3);
    end
else
    liq.FeOt_wt = NaN;
end

[H2O_value, hasH2O] = getOxideWithPresence(data_liq, ...
    {'H2O','H2Ot','H2Ototal','H2O_Total'}, 'Liq.H2O');

liq.usedDefaultH2O = false;
if hasH2O && isfinite(H2O_value)
    liq.H2O_wt = H2O_value;
elseif isfinite(defaultH2O_wt)
    liq.H2O_wt = defaultH2O_wt;
    liq.usedDefaultH2O = true;
else
    liq.H2O_wt = NaN;
end

% Oxide mole proportions. F and Cl are not read and therefore cannot enter
% these totals.
n = struct();
n.SiO2 = liq.SiO2_wt ./ MWinfo.MW.SiO2;
n.TiO2 = liq.TiO2_wt ./ MWinfo.MW.TiO2;
n.Al2O3 = liq.Al2O3_wt ./ MWinfo.MW.Al2O3;
n.FeOt = liq.FeOt_wt ./ MWinfo.MW.FeO;
n.MnO = liq.MnO_wt ./ MWinfo.MW.MnO;
n.MgO = liq.MgO_wt ./ MWinfo.MW.MgO;
n.CaO = liq.CaO_wt ./ MWinfo.MW.CaO;
n.Na2O = liq.Na2O_wt ./ MWinfo.MW.Na2O;
n.K2O = liq.K2O_wt ./ MWinfo.MW.K2O;
n.P2O5 = liq.P2O5_wt ./ MWinfo.MW.P2O5;
n.H2O = liq.H2O_wt ./ 18.01528;

anhydVector = [n.SiO2, n.TiO2, n.Al2O3, n.FeOt, n.MnO, ...
    n.MgO, n.CaO, n.Na2O, n.K2O, n.P2O5];

if all(isfinite(anhydVector))
    liq.cationTotal_anhyd = sum(anhydVector);
else
    liq.cationTotal_anhyd = NaN;
end

if isfinite(liq.cationTotal_anhyd) && isfinite(n.H2O)
    liq.cationTotal_hyd = liq.cationTotal_anhyd + n.H2O;
else
    liq.cationTotal_hyd = NaN;
end

liq.anhyd = makeMoleFractionStruct(n, liq.cationTotal_anhyd);
liq.hyd = makeMoleFractionStruct(n, liq.cationTotal_hyd);

if isfinite(n.FeOt) && isfinite(n.MgO) && n.MgO > 0
    liq.FeMgMolarRatio = n.FeOt ./ n.MgO;
else
    liq.FeMgMolarRatio = NaN;
end

end

function x = makeMoleFractionStruct(n, denominator)
% makeMoleFractionStruct
% Calculate oxide mole fractions using a supplied positive finite total.

fieldNames = {'SiO2','TiO2','Al2O3','FeOt','MnO','MgO', ...
    'CaO','Na2O','K2O','P2O5','H2O'};

x = struct();
for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    if isfinite(denominator) && denominator > 0
        x.(fieldName) = n.(fieldName) ./ denominator;
    else
        x.(fieldName) = NaN;
    end
end

end

function row = calcPressure(amp, liq, T_degreeC)
% calcPressure
% Calculate Putirka (2016) Equations 7a-7c for one Amphibole-Liquid pair.
% The equations do not use T; scalar pressures are repeated to match the
% input temperature-vector length.

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

DAl_hyd_scalar = amp.Al ./ liq.hyd.Al2O3;
DNa_hyd_scalar = amp.Na ./ liq.hyd.Na2O;
DAl_anhyd_scalar = amp.Al ./ liq.anhyd.Al2O3;
DNa_anhyd_scalar = amp.Na ./ liq.anhyd.Na2O;

lnDAl_hyd_scalar = safeLogPositive(DAl_hyd_scalar);
lnDNa_hyd_scalar = safeLogPositive(DNa_hyd_scalar);
lnDAl_anhyd_scalar = safeLogPositive(DAl_anhyd_scalar);
lnDNa_anhyd_scalar = safeLogPositive(DNa_anhyd_scalar);
lnXAl2O3_hyd_scalar = safeLogPositive(liq.hyd.Al2O3);
lnXSiO2_hyd_scalar = safeLogPositive(liq.hyd.SiO2);
lnXNaK_hyd_scalar = safeLogPositive(liq.hyd.Na2O + liq.hyd.K2O);
lnXAl2O3_anhyd_scalar = safeLogPositive(liq.anhyd.Al2O3);

isEquation7aDomain_scalar = allPositiveFinite([ ...
    DAl_hyd_scalar, liq.hyd.Al2O3, ...
    DNa_hyd_scalar, liq.hyd.Na2O]) && ...
    all(isfinite([amp.Al, amp.Na, amp.K, ...
    liq.hyd.P2O5, liq.hyd.H2O]));

isEquation7bDomain_scalar = allPositiveFinite([ ...
    DAl_hyd_scalar, liq.hyd.SiO2, ...
    liq.hyd.Na2O + liq.hyd.K2O]) && ...
    all(isfinite([amp.Al, liq.hyd.P2O5, liq.hyd.CaO]));

isEquation7cDomain_scalar = allPositiveFinite([ ...
    DAl_anhyd_scalar, liq.anhyd.Al2O3, ...
    DNa_anhyd_scalar, liq.anhyd.Na2O]) && ...
    all(isfinite([amp.Al, amp.Na, amp.K, liq.anhyd.P2O5]));

P7a_scalar = NaN;
if isEquation7aDomain_scalar
    P7a_scalar = -30.93 ...
        - 42.74 .* lnDAl_hyd_scalar ...
        - 42.16 .* lnXAl2O3_hyd_scalar ...
        + 633 .* liq.hyd.P2O5 ...
        + 12.64 .* liq.hyd.H2O ...
        + 24.57 .* amp.Al ...
        + 18.6 .* amp.K ...
        + 4.0 .* lnDNa_hyd_scalar;
end

P7b_scalar = NaN;
if isEquation7bDomain_scalar
    P7b_scalar = -64.79 ...
        - 6.064 .* lnDAl_hyd_scalar ...
        + 61.75 .* liq.hyd.SiO2 ...
        + 682 .* liq.hyd.P2O5 ...
        - 101.9 .* liq.hyd.CaO ...
        + 7.85 .* amp.Al ...
        - 46.46 .* lnXSiO2_hyd_scalar ...
        - 4.81 .* lnXNaK_hyd_scalar;
end

P7c_scalar = NaN;
if isEquation7cDomain_scalar
    P7c_scalar = -45.5 ...
        - 46.3 .* lnDAl_anhyd_scalar ...
        - 51.1 .* lnXAl2O3_anhyd_scalar ...
        + 439 .* liq.anhyd.P2O5 ...
        + 26.6 .* amp.Al ...
        + 22.5 .* amp.K ...
        + 5.23 .* lnDNa_anhyd_scalar;
end

% Fe-Mg equilibrium-screening coefficient.
liquidFeMgRatio = liq.FeMgMolarRatio;
amphiboleFeMgRatio = amp.FeT ./ amp.Mg;
KD_FeMg_scalar = amphiboleFeMgRatio ./ liquidFeMgRatio;

isWithinBroadKDRange_scalar = ...
    isfinite(KD_FeMg_scalar) && ...
    KD_FeMg_scalar >= 0.13 && KD_FeMg_scalar <= 0.41;

isWithinMean1SigmaKDRange_scalar = ...
    isfinite(KD_FeMg_scalar) && ...
    KD_FeMg_scalar >= 0.17 && KD_FeMg_scalar <= 0.39;

% Expand all composition-dependent scalar values.
P7a = repmat(P7a_scalar, nT, 1);
P7b = repmat(P7b_scalar, nT, 1);
P7c = repmat(P7c_scalar, nT, 1);

isWithinExperimentalTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 650 & T_degreeC <= 1175;
isWithinPreferredTFilter = ...
    isfinite(T_degreeC) & T_degreeC <= 800;
isWithinPracticalPRange_7a = ...
    isfinite(P7a) & P7a >= 1 & P7a <= 8;
isWithinPracticalPRange_7b = ...
    isfinite(P7b) & P7b >= 1 & P7b <= 8;
isWithinPracticalPRange_7c = ...
    isfinite(P7c) & P7c >= 1 & P7c <= 8;

isWithinCationSumRange_scalar = ...
    isfinite(amp.SumCat) && amp.SumCat >= 15 && amp.SumCat <= 16;

isRecommended_AS1995_filter = ...
    isWithinPreferredTFilter & ...
    repmat(isfinite(amp.FeNumber) && amp.FeNumber < 0.65, nT, 1);

isApplicable_numeric_7a = ...
    repmat(isEquation7aDomain_scalar, nT, 1) & ...
    repmat(amp.isCalcic_numeric, nT, 1) & ...
    repmat(isWithinBroadKDRange_scalar, nT, 1) & ...
    isRecommended_AS1995_filter & ...
    isWithinPracticalPRange_7a;

row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.temperatureUsedInPressureEquation = false(nT, 1);

row.Si_amp_23O = repmat(amp.Si, nT, 1);
row.Ti_amp_23O = repmat(amp.Ti, nT, 1);
row.Al_amp_23O = repmat(amp.Al, nT, 1);
row.FeT_amp_23O = repmat(amp.FeT, nT, 1);
row.Mn_amp_23O = repmat(amp.Mn, nT, 1);
row.Mg_amp_23O = repmat(amp.Mg, nT, 1);
row.Ca_amp_23O = repmat(amp.Ca, nT, 1);
row.Na_amp_23O = repmat(amp.Na, nT, 1);
row.K_amp_23O = repmat(amp.K, nT, 1);
row.Cr_amp_23O = repmat(amp.Cr, nT, 1);
row.Ni_amp_23O = repmat(amp.Ni, nT, 1);
row.SumCat_amp_23O = repmat(amp.SumCat, nT, 1);
row.FeNumber_amp = repmat(amp.FeNumber, nT, 1);
row.Ca_B_amp = repmat(amp.Ca_B, nT, 1);
row.Na_B_amp = repmat(amp.Na_B, nT, 1);
row.CaNa_B_amp = repmat(amp.CaNa_B, nT, 1);

row.SiO2_liq = repmat(liq.SiO2_wt, nT, 1);
row.TiO2_liq = repmat(liq.TiO2_wt, nT, 1);
row.Al2O3_liq = repmat(liq.Al2O3_wt, nT, 1);
row.FeOt_liq = repmat(liq.FeOt_wt, nT, 1);
row.MnO_liq = repmat(liq.MnO_wt, nT, 1);
row.MgO_liq = repmat(liq.MgO_wt, nT, 1);
row.CaO_liq = repmat(liq.CaO_wt, nT, 1);
row.Na2O_liq = repmat(liq.Na2O_wt, nT, 1);
row.K2O_liq = repmat(liq.K2O_wt, nT, 1);
row.P2O5_liq = repmat(liq.P2O5_wt, nT, 1);
row.H2O_liq = repmat(liq.H2O_wt, nT, 1);

row.cationTotal_liq_anhyd = ...
    repmat(liq.cationTotal_anhyd, nT, 1);
row.cationTotal_liq_hyd = ...
    repmat(liq.cationTotal_hyd, nT, 1);
row.F_Cl_excluded_from_cationTotal_liq = true(nT, 1);

row.XSiO2_liq_hyd = repmat(liq.hyd.SiO2, nT, 1);
row.XTiO2_liq_hyd = repmat(liq.hyd.TiO2, nT, 1);
row.XAl2O3_liq_hyd = repmat(liq.hyd.Al2O3, nT, 1);
row.XFeOt_liq_hyd = repmat(liq.hyd.FeOt, nT, 1);
row.XMnO_liq_hyd = repmat(liq.hyd.MnO, nT, 1);
row.XMgO_liq_hyd = repmat(liq.hyd.MgO, nT, 1);
row.XCaO_liq_hyd = repmat(liq.hyd.CaO, nT, 1);
row.XNa2O_liq_hyd = repmat(liq.hyd.Na2O, nT, 1);
row.XK2O_liq_hyd = repmat(liq.hyd.K2O, nT, 1);
row.XP2O5_liq_hyd = repmat(liq.hyd.P2O5, nT, 1);
row.XH2O_liq_hyd = repmat(liq.hyd.H2O, nT, 1);
row.XNaK_liq_hyd = ...
    repmat(liq.hyd.Na2O + liq.hyd.K2O, nT, 1);

row.XSiO2_liq_anhyd = repmat(liq.anhyd.SiO2, nT, 1);
row.XTiO2_liq_anhyd = repmat(liq.anhyd.TiO2, nT, 1);
row.XAl2O3_liq_anhyd = repmat(liq.anhyd.Al2O3, nT, 1);
row.XFeOt_liq_anhyd = repmat(liq.anhyd.FeOt, nT, 1);
row.XMnO_liq_anhyd = repmat(liq.anhyd.MnO, nT, 1);
row.XMgO_liq_anhyd = repmat(liq.anhyd.MgO, nT, 1);
row.XCaO_liq_anhyd = repmat(liq.anhyd.CaO, nT, 1);
row.XNa2O_liq_anhyd = repmat(liq.anhyd.Na2O, nT, 1);
row.XK2O_liq_anhyd = repmat(liq.anhyd.K2O, nT, 1);
row.XP2O5_liq_anhyd = repmat(liq.anhyd.P2O5, nT, 1);

row.DAl_hyd = repmat(DAl_hyd_scalar, nT, 1);
row.DNa_hyd = repmat(DNa_hyd_scalar, nT, 1);
row.DAl_anhyd = repmat(DAl_anhyd_scalar, nT, 1);
row.DNa_anhyd = repmat(DNa_anhyd_scalar, nT, 1);
row.lnDAl_hyd = repmat(lnDAl_hyd_scalar, nT, 1);
row.lnDNa_hyd = repmat(lnDNa_hyd_scalar, nT, 1);
row.lnDAl_anhyd = repmat(lnDAl_anhyd_scalar, nT, 1);
row.lnDNa_anhyd = repmat(lnDNa_anhyd_scalar, nT, 1);

row.KD_FeMg_Amp_Liq = repmat(KD_FeMg_scalar, nT, 1);
row.isWithinBroadKDRange = ...
    repmat(isWithinBroadKDRange_scalar, nT, 1);
row.isWithinMean1SigmaKDRange = ...
    repmat(isWithinMean1SigmaKDRange_scalar, nT, 1);

% Primary launcher-compatible pressure and equation-specific outputs.
row.P_kbar = P7a;
row.P_Putirka2016_7a_kbar = P7a;
row.P_Putirka2016_7b_kbar = P7b;
row.P_Putirka2016_7c_kbar = P7c;
row.P7a_minus_P7b_kbar = P7a - P7b;

row.P_uncertainty_individual_calibration_kbar = repmat(2.0, nT, 1);
row.P_uncertainty_individual_test_approx_kbar = repmat(4.0, nT, 1);
row.P_uncertainty_averaged_1to8kbar_approx_kbar = repmat(1.0, nT, 1);
row.P_H2O_sensitivity_kbar_per_wtpercent = repmat(0.4, nT, 1);

row.isWithinEquation7aDomain = ...
    repmat(isEquation7aDomain_scalar, nT, 1);
row.isWithinEquation7bDomain = ...
    repmat(isEquation7bDomain_scalar, nT, 1);
row.isWithinEquation7cDomain = ...
    repmat(isEquation7cDomain_scalar, nT, 1);
row.isWithinExperimentalTRange = isWithinExperimentalTRange;
row.isWithinPreferredTFilter = isWithinPreferredTFilter;
row.isWithinPracticalPRange_7a = isWithinPracticalPRange_7a;
row.isWithinPracticalPRange_7b = isWithinPracticalPRange_7b;
row.isWithinPracticalPRange_7c = isWithinPracticalPRange_7c;
row.isFeNumberBelow065 = ...
    repmat(isfinite(amp.FeNumber) && amp.FeNumber < 0.65, nT, 1);
row.isCalcicAmphibole_numeric = ...
    repmat(amp.isCalcic_numeric, nT, 1);
row.isWithinCationSumRange = ...
    repmat(isWithinCationSumRange_scalar, nT, 1);
row.isRecommended_AS1995_filter = isRecommended_AS1995_filter;
row.isApplicable_numeric_7a = isApplicable_numeric_7a;

% Backward-compatible general diagnostic flag.
row.isApplicable = isApplicable_numeric_7a;

end

function nanInputNames = findNaNInputs(amp, liq, T_degreeC)
% findNaNInputs
% Return canonical names of NaN calculation or numerical-screening inputs.
% Liquid F and Cl are intentionally excluded.

maxNames = 20;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

ampNames = {'Si','Al','FeT','Mg','Ca','Na','K'};
for i = 1:numel(ampNames)
    fieldName = ampNames{i};
    if isnan(amp.(fieldName))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = ...
            "Amphibole." + string(fieldName) + "_cation_apfu";
    end
end

liqNames = {'SiO2','TiO2','Al2O3','FeOt','MnO','MgO', ...
    'CaO','Na2O','K2O','P2O5','H2O'};
liqValues = [liq.SiO2_wt, liq.TiO2_wt, liq.Al2O3_wt, ...
    liq.FeOt_wt, liq.MnO_wt, liq.MgO_wt, liq.CaO_wt, ...
    liq.Na2O_wt, liq.K2O_wt, liq.P2O5_wt, liq.H2O_wt];

for i = 1:numel(liqNames)
    if isnan(liqValues(i))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = ...
            "Liq." + string(liqNames{i});
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function printPressureDiagnostics( ...
        pressureVector, equationLabel, practicalP_min_kbar, ...
        practicalP_max_kbar, highP_threshold_kbar, ...
        selectedCode_amp, selectedIdx_liq)
% printPressureDiagnostics
% Print non-stopping range and finite-value warnings for one equation.

finitePressure = isfinite(pressureVector);

if any(finitePressure)
    finiteValues = pressureVector(finitePressure);
    outsidePractical = finitePressure & ...
        (pressureVector < practicalP_min_kbar | ...
         pressureVector > practicalP_max_kbar);

    if any(outsidePractical)
        fprintf(2, ['WARNING: %s pressure is outside the practical ' ...
            'approximately 1-8 kbar interval emphasized by Putirka ' ...
            '(2016; p. 841; pp. 852-854) for %s & Liquid row %d. ' ...
            '%d of %d finite row(s) are outside; finite pressure range = ' ...
            '%.6g-%.6g kbar. Values are retained.\n'], ...
            equationLabel, ...
            char(string(selectedCode_amp)), selectedIdx_liq, ...
            sum(outsidePractical), sum(finitePressure), ...
            min(finiteValues), max(finiteValues));
    end

    aboveHighP = finitePressure & pressureVector > highP_threshold_kbar;
    if any(aboveHighP)
        fprintf(2, ['WARNING: %s gives pressure above 10 kbar for %s & ' ...
            'Liquid row %d (%d of %d finite row(s)). Putirka (2016; ' ...
            'pp. 850-854) recommends comparing Equations 7a and 7b and ' ...
            'preferring Equation 7b or Ridolfi and Renzulli (2011) ' ...
            'model 1d at P > 10 kbar.\n'], ...
            equationLabel, ...
            char(string(selectedCode_amp)), selectedIdx_liq, ...
            sum(aboveHighP), sum(finitePressure));
    end

    negativePressure = finitePressure & pressureVector < 0;
    if any(negativePressure)
        fprintf(2, ['WARNING: %s calculated negative finite pressure for ' ...
            '%s & Liquid row %d (%d of %d row(s)). Negative values are ' ...
            'retained for diagnostic purposes and must not be interpreted ' ...
            'as physical negative pressure.\n'], ...
            equationLabel, ...
            char(string(selectedCode_amp)), selectedIdx_liq, ...
            sum(negativePressure), numel(pressureVector));
    end
end

invalidPressure = ~isfinite(pressureVector);
if any(invalidPressure)
    fprintf(2, ['WARNING: %s calculated non-finite pressure for %s & ' ...
        'Liquid row %d (%d of %d row(s); NaN: %d, Inf: %d). ' ...
        'The values remain in the output table.\n'], ...
        equationLabel, ...
        char(string(selectedCode_amp)), selectedIdx_liq, ...
        sum(invalidPressure), numel(pressureVector), ...
        sum(isnan(pressureVector)), sum(isinf(pressureVector)));
end

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Attach commonly used Liquid identifiers to every output row.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if ismember('Index', variableNames)
    row.liq_Index = repmat(string(data_liq.Index), nRows, 1);
end
if ismember('Experiment', variableNames)
    row.liq_Experiment = repmat(string(data_liq.Experiment), nRows, 1);
end
if ismember('Citation', variableNames)
    row.liq_Citation = repmat(string(data_liq.Citation), nRows, 1);
end

end

function items = buildLiquidList(dataset_liq)
% buildLiquidList
% Build a pre-sized display list for interactive Liquid selection.

nRows = height(dataset_liq);
items = cell(nRows, 1);

hasIndex = ismember('Index', dataset_liq.Properties.VariableNames);
hasExperiment = ...
    ismember('Experiment', dataset_liq.Properties.VariableNames);
hasCitation = ...
    ismember('Citation', dataset_liq.Properties.VariableNames);

for i = 1:nRows
    indexText = "";
    experimentText = "";
    citationText = "";

    if hasIndex
        indexText = " | Index=" + string(dataset_liq.Index(i));
    end
    if hasExperiment
        experimentText = " | " + string(dataset_liq.Experiment(i));
    end
    if hasCitation
        citationText = " | " + string(dataset_liq.Citation(i));
    end

    items{i} = char("Row " + string(i) + ...
        indexText + experimentText + citationText);
end

end

function value = getFirstVariableOrNaN( ...
        tableRow, candidateNames, displayName)
% getFirstVariableOrNaN
% Retrieve the first matching scalar variable. Missing variables and explicit
% NaN remain NaN. Inf and finite negative values are rejected.

value = NaN;
variableNames = tableRow.Properties.VariableNames;

for i = 1:numel(candidateNames)
    candidateName = candidateNames{i};
    if ismember(candidateName, variableNames)
        value = convertToScalarDouble(tableRow.(candidateName));
        validateNonNegativeScalar(value, displayName);
        return;
    end
end

end

function value = getOxideOrNaN(tableRow, candidates, displayName)
% getOxideOrNaN
% Retrieve a Liquid oxide scalar using flexible column-name matching.

[value, isPresent] = ...
    getOxideWithPresence(tableRow, candidates, displayName);
if ~isPresent
    value = NaN;
end

end

function [value, isPresent] = getOxideWithPresence( ...
        tableRow, candidates, displayName)
% getOxideWithPresence
% Retrieve a Liquid oxide scalar and report whether a matching column exists.

columnName = findColumnByCandidates( ...
    tableRow.Properties.VariableNames, candidates);

if isempty(columnName)
    value = NaN;
    isPresent = false;
    return;
end

value = convertToScalarDouble(tableRow.(columnName));
validateNonNegativeScalar(value, displayName);
isPresent = true;

end

function columnName = findColumnByCandidates(variableNames, candidates)
% findColumnByCandidates
% Find a table column after removing spaces, underscores, and hyphens.

canonicalVariables = strings(numel(variableNames), 1);
for i = 1:numel(variableNames)
    canonicalVariables(i) = canonicalizeName(variableNames{i});
end

columnName = '';
for i = 1:numel(candidates)
    canonicalCandidate = canonicalizeName(candidates{i});
    possibleNames = [canonicalCandidate; canonicalCandidate + "value"];

    for j = 1:numel(possibleNames)
        matchIndex = find( ...
            canonicalVariables == possibleNames(j), 1, 'first');
        if ~isempty(matchIndex)
            columnName = variableNames{matchIndex};
            return;
        end
    end
end

end

function canonicalName = canonicalizeName(inputName)
% canonicalizeName
% Convert a variable name to a lower-case comparison key.

canonicalName = lower(string(inputName));
canonicalName = replace(canonicalName, " ", "");
canonicalName = replace(canonicalName, "_", "");
canonicalName = replace(canonicalName, "-", "");
canonicalName = replace(canonicalName, ".", "");

end

function value = convertToScalarDouble(rawValue)
% convertToScalarDouble
% Convert one scalar numeric, string, char, or cell value to double. Missing
% or non-convertible values are returned as NaN.

if isempty(rawValue)
    value = NaN;
    return;
end

if numel(rawValue) ~= 1
    error('Selected table variable must contain one scalar value.');
end

if isnumeric(rawValue)
    value = double(rawValue);
    return;
end

if islogical(rawValue)
    value = double(rawValue);
    return;
end

if isstring(rawValue)
    if ismissing(rawValue)
        value = NaN;
    else
        value = str2double(rawValue);
    end
    return;
end

if ischar(rawValue)
    value = str2double(string(rawValue));
    return;
end

if iscell(rawValue)
    cellValue = rawValue{1};
    if isempty(cellValue)
        value = NaN;
    elseif isnumeric(cellValue) || islogical(cellValue)
        if numel(cellValue) ~= 1
            error('Selected cell value must contain one scalar value.');
        end
        value = double(cellValue);
    else
        value = str2double(string(cellValue));
    end
    return;
end

value = NaN;

end

function validateNonNegativeScalar(value, displayName)
% validateNonNegativeScalar
% Allow zero and NaN; reject Inf and finite negative values.

if isinf(value) || (isfinite(value) && value < 0)
    error('%s must be non-negative or NaN; Inf is prohibited.', displayName);
end

end

function logarithm = safeLogPositive(value)
% safeLogPositive
% Return log(value) only for a finite strictly positive scalar.

if isfinite(value) && value > 0
    logarithm = log(value);
else
    logarithm = NaN;
end

end

function tf = allPositiveFinite(values)
% allPositiveFinite
% Return true only when all values are finite and strictly positive.

tf = all(isfinite(values)) && all(values > 0);

end
