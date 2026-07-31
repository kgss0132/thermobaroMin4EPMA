function results = DSouza2020(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Olivine_Pyroxene/DSouza2020.m
% Tested with MATLAB R2024b
%
% Spinel-peridotite thermobarometer using Al and Ca in olivine
% D'Souza, R.J., Canil, D. & Coogan, L.A. (2020)
% Contributions to Mineralogy and Petrology, 175, 5
% DOI: https://doi.org/10.1007/s00410-019-1647-6
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function combines:
%
%   1) the Coogan et al. (2014) Al-in-olivine thermometer, using one
%      Olivine-Spinel pair, and
%   2) the Kohler & Brey (1990) Ca exchange barometer, using one
%      Olivine-Clinopyroxene pair.
%
% The pressure calculation therefore uses a temperature calculated
% internally from Olivine-Spinel Al partitioning. T_degreeC is accepted only
% for compatibility with startBaroCalc_fixedT and startBaroCalc_rangeT; it is
% not used in the D'Souza et al. (2020) calculation. The output variable
% T_degreeC is the internally calculated Al-in-olivine temperature. The
% launcher-supplied values are retained as T_launcher_degreeC.
%
% For a vector T_degreeC input, the internally calculated P-T result is
% repeated to preserve the common launcher row structure. The repeated rows
% differ only in T_launcher_degreeC and T_launcher_K.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% D'Souza et al. (2020) tested the pressure dependence of the Al-in-olivine
% thermometer at:
%
%   Temperature : 950-1250 degreeC
%   Pressure    : 1.5-2.4 GPa (15-24 kbar)
%
% No analytically resolvable pressure effect was observed over this range
% for the successfully equilibrated low-Cr spinel experiments. Experiments
% using Cr-rich spinel (Cr# approximately 0.72) did not demonstrably
% equilibrate. Application to high-Cr spinel therefore requires particular
% caution. Coogan et al. (2014) calibrated the compositional correction up
% to spinel Cr# approximately 0.69.
%
% The Kohler & Brey (1990) Ca-in-olivine barometer was calibrated at:
%
%   Temperature : 900-1400 degreeC
%   Pressure    : 2-60 kbar
%   Olivine Mg# : approximately 0.9
%   Precision   : +/-1.7 kbar (1 sigma)
%
% IMPORTANT APPLICATION REQUIREMENTS:
%   - The selected Olivine, Spinel, and Clinopyroxene must represent a
%     coexisting, equilibrated spinel-peridotite assemblage.
%   - Olivine Al and Ca must be measured with sufficient precision. Ca in
%     olivine is especially vulnerable to secondary fluorescence from nearby
%     clinopyroxene during EPMA analysis.
%   - Ca diffuses relatively rapidly in olivine and may be reset during
%     xenolith transport, transient heating, cooling, or exhumation.
%   - Al and Ca zoning in olivine should be assessed before interpreting the
%     calculated P-T conditions as a single equilibrium state.
%   - This function does not test textural equilibrium, phase coexistence,
%     zoning, or analytical secondary fluorescence.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT
% rawdata_struct must contain one-row-selectable tables:
%
%   rawdata_struct.Olivine
%   rawdata_struct.Cpx
%   rawdata_struct.Spinel   (rawdata_struct.Sp is also accepted)
%
% Required analytical information:
%
%   Olivine:
%     Ca_cation_apfu normalized to 4 O, or oxide wt% sufficient to calculate
%     the 4-O formula; and Al2O3 wt% or elemental Al in ppm.
%
%   Clinopyroxene:
%     Ca_cation_apfu normalized to 6 O, or oxide wt% sufficient to calculate
%     the 6-O formula.
%
%   Spinel:
%     Al2O3 wt%; Cr, Al, and Fe3+ cations normalized to 4 O. If Fe3+ is not
%     supplied, it is estimated by spinel charge balance when possible.
%
% Missing or non-numeric calculation inputs remain NaN and propagate through
% the calculation. Finite negative calculation inputs and Inf are rejected.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
% Coogan et al. (2014):
%
%   T_AlOl(K) = 10000 / [
%       0.575
%       + 0.884*Cr#_sp
%       - 0.897*ln(D_Al2O3_ol_sp)
%   ]
%
% where D'Souza et al. (2020, Table 4) define:
%
%   D_Al2O3_ol_sp = Al2O3_ol / Al2O3_sp
%   Cr#_sp = Cr / (Cr + Al + Fe3+)
%
% Al2O3 in Olivine and Spinel must use the same concentration unit. This
% implementation uses wt%. If elemental Al in olivine is supplied in ppm, it
% is converted to equivalent Al2O3 wt% before calculating the ratio.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION
% Kohler & Brey (1990):
%
%   D_Ca_ol_cpx = Ca_ol(4 O) / Ca_cpx(6 O)
%
% High-temperature branch:
%
%   P_high(kbar) =
%       [-T*ln(D_Ca) - 11982 + 3.61*T] / 56.2
%
%   valid when T >= 1275.25 + 2.827*P_high
%
% Low-temperature branch:
%
%   P_low(kbar) =
%       [-T*ln(D_Ca) - 5792 - 1.25*T] / 42.5
%
%   valid when T <= 1275.25 + 2.827*P_low
%
% T is in Kelvin. Each candidate pressure is tested against its own branch
% criterion. If both or neither criterion is satisfied because of numerical
% or compositional inconsistency, the candidate closest to the branch
% boundary is retained and the ambiguity is reported.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = DSouza2020(rawdata_struct, T_degreeC)
%   results = DSouza2020(rawdata_struct, T_degreeC, ...
%       'OlivineRow', olRow, 'CpxRow', cpxRow, 'SpinelRow', spRow)
%
% Inputs:
%   rawdata_struct : struct containing Olivine, Cpx, and Spinel/Sp tables
%   T_degreeC      : numeric scalar or vector supplied by the common
%                    barometer launcher. It is retained but not used.
%
% Optional name-value inputs:
%   OlivineRow : positive integer row number, or [] for interactive selection
%   CpxRow     : positive integer row number, or [] for interactive selection
%   SpinelRow  : positive integer row number, or [] for interactive selection
%
% Output:
%   results : table containing one repeated row per launcher temperature for
%             every selected Olivine-Cpx-Spinel triplet.
%

%% Input validation
if nargin < 2
    error('DSouza2020 requires (rawdata_struct, T_degreeC).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(T_degreeC) || isempty(T_degreeC) || ~isvector(T_degreeC) || ...
        any(isinf(T_degreeC(:)))
    error(['T_degreeC must be a numeric scalar or vector. NaN is allowed ' ...
           'and retained, but Inf is prohibited.']);
end

T_degreeC = T_degreeC(:);

[dataset_ol, olFieldName] = resolveMineralTable( ...
    rawdata_struct, 'Olivine', {'Ol', 'olivine'});
[dataset_cpx, cpxFieldName] = resolveMineralTable( ...
    rawdata_struct, 'Cpx', {'Clinopyroxene', 'cpx'});
[dataset_sp, spFieldName] = resolveMineralTable( ...
    rawdata_struct, 'Spinel', {'Sp', 'spinel'});

ip = inputParser;
rowValidator = @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == fix(x));
ip.addParameter('OlivineRow', [], rowValidator);
ip.addParameter('CpxRow', [], rowValidator);
ip.addParameter('SpinelRow', [], rowValidator);
ip.parse(varargin{:});

olivineRowOption = ip.Results.OlivineRow;
cpxRowOption = ip.Results.CpxRow;
spinelRowOption = ip.Results.SpinelRow;

validateRequestedRow(olivineRowOption, height(dataset_ol), 'OlivineRow');
validateRequestedRow(cpxRowOption, height(dataset_cpx), 'CpxRow');
validateRequestedRow(spinelRowOption, height(dataset_sp), 'SpinelRow');

%% 1) Prepare datasets
disp('=== Step 1: Preparing Olivine, Cpx, and Spinel datasets ===');

olivineItems = buildMineralList(dataset_ol, 'Olivine');
cpxItems = buildMineralList(dataset_cpx, 'Cpx');
spinelItems = buildMineralList(dataset_sp, 'Spinel');

disp('=== Preparing mineral datasets has been finished ===');

%% 2) Initialize output container
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max([16, height(dataset_ol), ...
    height(dataset_cpx), height(dataset_sp)]);
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

allRowsFixed = ~isempty(olivineRowOption) && ...
    ~isempty(cpxRowOption) && ~isempty(spinelRowOption);
fixedTripletProcessed = false;

inputTemperatureCautionIssued = false;
analyticalCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-6) Interactive selection and calculation
while true
    if allRowsFixed && fixedTripletProcessed
        break;
    end

    % ----- Olivine selection -----
    disp('=== Step 3: Selecting a data code from the list (Olivine) ===');
    [selectedIdx_ol, ok] = selectMineralRow( ...
        olivineItems, olivineRowOption, 'Olivine');
    if ~ok
        disp('Selection canceled');
        break;
    end

    % ----- Clinopyroxene selection -----
    disp('=== Step 4: Selecting a data code from the list (Cpx) ===');
    [selectedIdx_cpx, ok] = selectMineralRow( ...
        cpxItems, cpxRowOption, 'Clinopyroxene');
    if ~ok
        disp('Selection canceled');
        break;
    end

    % ----- Spinel selection -----
    disp('=== Step 5: Selecting a data code from the list (Spinel) ===');
    [selectedIdx_sp, ok] = selectMineralRow( ...
        spinelItems, spinelRowOption, 'Spinel');
    if ~ok
        disp('Selection canceled');
        break;
    end

    fixedTripletProcessed = true;

    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);
    selectedData_sp = dataset_sp(selectedIdx_sp, :);

    selectedCode_ol = getDataCode(selectedData_ol, selectedIdx_ol);
    selectedCode_cpx = getDataCode(selectedData_cpx, selectedIdx_cpx);
    selectedCode_sp = getDataCode(selectedData_sp, selectedIdx_sp);

    disp(['Olivine selected: ' selectedCode_ol]);
    disp(['Cpx selected: ' selectedCode_cpx]);
    disp(['Spinel selected: ' selectedCode_sp]);

    % ----- Calculation -----
    disp('=== Step 6: Checking inputs and calculating P-T ===');

    inputs = prepareCalculationInputs( ...
        selectedData_ol, selectedData_cpx, selectedData_sp);

    validateNonNegativeInputs(inputs);
    nanInputNames = findNaNInputs(inputs);

    row = calcPressure(inputs, T_degreeC);

    nRows = height(row);
    row.dataRow_ol = repmat(selectedIdx_ol, nRows, 1);
    row.dataRow_cpx = repmat(selectedIdx_cpx, nRows, 1);
    row.dataRow_sp = repmat(selectedIdx_sp, nRows, 1);
    row.dataCode_ol = repmat(string(selectedCode_ol), nRows, 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), nRows, 1);
    row.dataCode_sp = repmat(string(selectedCode_sp), nRows, 1);

    row = movevars(row, ...
        {'dataRow_ol', 'dataRow_cpx', 'dataRow_sp', ...
         'dataCode_ol', 'dataCode_cpx', 'dataCode_sp'}, ...
        'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % ----- Immediate output -----
    disp('--------------------------------------------------');
    disp('=== P-T conditions were calculated: ===');
    disp([selectedCode_ol ' & ' selectedCode_cpx ' & ' selectedCode_sp ...
        ': T_AlOl = ' num2str(row.T_AlOl_degreeC(1)) ...
        ' degreeC; P = ' num2str(row.P_kbar(1)) ' kbar']);

    if ~inputTemperatureCautionIssued
        fprintf(2, ...
            ['CAUTION: DSouza2020 calculates temperature internally from ' ...
             'Olivine-Spinel Al partitioning. The launcher-supplied ' ...
             'T_degreeC values are retained as T_launcher_degreeC but are ' ...
             'not used in pressure calculation.\n']);
        if numel(T_degreeC) > 1
            fprintf(2, ...
                ['         A vector input therefore produces repeated P-T ' ...
                 'rows that differ only in launcher-temperature metadata.\n']);
        end
        inputTemperatureCautionIssued = true;
    end

    if ~analyticalCautionIssued
        fprintf(2, ...
            ['CAUTION: High-precision Ca and Al analyses of olivine are ' ...
             'required. Check Ca secondary fluorescence, Ol-Cpx-Sp ' ...
             'equilibrium, and Ca/Al zoning before interpretation.\n']);
        analyticalCautionIssued = true;
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the calculation input(s) for %s, ' ...
             '%s, and %s: %s. NaN values were retained and propagated.\n'], ...
            selectedCode_ol, selectedCode_cpx, selectedCode_sp, ...
            char(strjoin(nanInputNames, ', ')));
    end

    Tcalc = row.T_AlOl_degreeC(1);
    Pcalc = row.P_kbar(1);
    CrNumber = row.CrNumber_sp_DSouza(1);

    if isfinite(Tcalc) && (Tcalc < 950 || Tcalc > 1250)
        fprintf(2, ...
            ['WARNING: Calculated T_AlOl = %.4g degreeC is outside the ' ...
             '950-1250 degreeC experimental range tested by D''Souza et al. ' ...
             '(2020).\n'], Tcalc);
    end

    if isfinite(Pcalc) && (Pcalc < 15 || Pcalc > 24)
        fprintf(2, ...
            ['CAUTION: Calculated P = %.4g kbar is outside the 15-24 kbar ' ...
             'range over which D''Souza et al. (2020) directly tested the ' ...
             'pressure independence of the Al-in-olivine thermometer.\n'], ...
            Pcalc);
    end

    if isfinite(Tcalc) && (Tcalc < 900 || Tcalc > 1400)
        fprintf(2, ...
            ['WARNING: Calculated T_AlOl = %.4g degreeC is outside the ' ...
             '900-1400 degreeC calibration range of the Kohler & Brey ' ...
             '(1990) Ca-in-olivine barometer.\n'], Tcalc);
    end

    if isfinite(Pcalc) && (Pcalc < 2 || Pcalc > 60)
        fprintf(2, ...
            ['WARNING: Calculated P = %.4g kbar is outside the 2-60 kbar ' ...
             'calibration range of the Kohler & Brey (1990) barometer.\n'], ...
            Pcalc);
    end

    if isfinite(CrNumber) && CrNumber > 0.69
        fprintf(2, ...
            ['CAUTION: Spinel Cr# = %.4g exceeds approximately 0.69. ' ...
             'D''Souza et al. (2020) found that their Cr-rich spinel ' ...
             'experiments near Cr# 0.72 did not demonstrably equilibrate.\n'], ...
            CrNumber);
    end

    if ~strcmp(row.KB90_branchStatus(1), "unique")
        fprintf(2, ...
            ['WARNING: Kohler & Brey branch selection was not unique ' ...
             '(status: %s). The candidate closest to the branch boundary ' ...
             'was retained. Inspect P_KB90_high_kbar, P_KB90_low_kbar, and ' ...
             'the branch-residual columns.\n'], ...
            char(row.KB90_branchStatus(1)));
    end

    if any(~isfinite(row.P_kbar)) || any(~isfinite(row.T_AlOl_K))
        fprintf(2, ...
            ['WARNING: A non-finite P or T result was calculated for the ' ...
             'selected mineral triplet. The result was retained.\n']);
    end

    if isfinite(Pcalc) && Pcalc < 0
        fprintf(2, ...
            ['WARNING: Negative pressure was calculated (%.4g kbar). The ' ...
             'value was retained for diagnostic purposes.\n'], Pcalc);
    end

    disp('--------------------------------------------------');

    if allRowsFixed
        break;
    end

    userAction = questdlg( ...
        'Continue with another Olivine-Cpx-Spinel selection?', ...
        'DSouza2020', ...
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

results.Properties.UserData = struct( ...
    'OlivineField', olFieldName, ...
    'CpxField', cpxFieldName, ...
    'SpinelField', spFieldName, ...
    'LauncherTemperatureUsed', false, ...
    'PrimaryTemperature', 'Coogan et al. (2014) Al-in-olivine', ...
    'PrimaryPressure', 'Kohler and Brey (1990) Ca-in-olivine');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function inputs = prepareCalculationInputs(data_ol, data_cpx, data_sp)
% Extract the variables required by the combined thermobarometer.

inputs = struct();

[inputs.Ca_ol_4O, inputs.Ca_ol_source] = ...
    getCationWithOxideFallback(data_ol, 'Ca', 4, 'Olivine');
[inputs.Ca_cpx_6O, inputs.Ca_cpx_source] = ...
    getCationWithOxideFallback(data_cpx, 'Ca', 6, 'Cpx');

[inputs.Al2O3_ol_wtpercent, inputs.Al2O3_ol_source] = ...
    getAl2O3WtPercent(data_ol, 'Olivine');
[inputs.Al2O3_sp_wtpercent, inputs.Al2O3_sp_source] = ...
    getAl2O3WtPercent(data_sp, 'Spinel');

sp = prepareSpinelCations(data_sp);
inputs.Cr_sp_4O = sp.Cr;
inputs.Al_sp_4O = sp.Al;
inputs.Fe3_sp_4O = sp.Fe3;
inputs.FeTotal_sp_4O = sp.FeTotal;
inputs.Fe3_sp_source = sp.Fe3Source;
inputs.spinelCationSource = sp.CationSource;

end

function row = calcPressure(inputs, T_launcher_degreeC)
% Calculate Coogan et al. (2014) T and Kohler & Brey (1990) P.

T_launcher_degreeC = T_launcher_degreeC(:);
nT = numel(T_launcher_degreeC);
T_launcher_K = T_launcher_degreeC + 273.15;

% --- Olivine-Spinel Al thermometer ---
D_Al2O3_scalar = inputs.Al2O3_ol_wtpercent ./ ...
    inputs.Al2O3_sp_wtpercent;

if isfinite(D_Al2O3_scalar) && D_Al2O3_scalar > 0
    lnD_Al2O3_scalar = log(D_Al2O3_scalar);
else
    lnD_Al2O3_scalar = NaN;
end

CrDenominator_DSouza = ...
    inputs.Cr_sp_4O + inputs.Al_sp_4O + inputs.Fe3_sp_4O;
if isfinite(CrDenominator_DSouza) && CrDenominator_DSouza > 0
    CrNumber_DSouza_scalar = ...
        inputs.Cr_sp_4O ./ CrDenominator_DSouza;
else
    CrNumber_DSouza_scalar = NaN;
end

CrDenominator_Coogan = inputs.Cr_sp_4O + inputs.Al_sp_4O;
if isfinite(CrDenominator_Coogan) && CrDenominator_Coogan > 0
    CrNumber_Coogan_scalar = inputs.Cr_sp_4O ./ CrDenominator_Coogan;
else
    CrNumber_Coogan_scalar = NaN;
end

T_denominator_scalar = ...
    0.575 + 0.884 .* CrNumber_DSouza_scalar ...
    - 0.897 .* lnD_Al2O3_scalar;

if isfinite(T_denominator_scalar) && T_denominator_scalar > 0
    T_AlOl_K_scalar = 10000 ./ T_denominator_scalar;
    T_AlOl_degreeC_scalar = T_AlOl_K_scalar - 273.15;
else
    T_AlOl_K_scalar = NaN;
    T_AlOl_degreeC_scalar = NaN;
end

% --- Olivine-Clinopyroxene Ca barometer ---
D_Ca_scalar = inputs.Ca_ol_4O ./ inputs.Ca_cpx_6O;
if isfinite(D_Ca_scalar) && D_Ca_scalar > 0
    lnD_Ca_scalar = log(D_Ca_scalar);
else
    lnD_Ca_scalar = NaN;
end

if isfinite(T_AlOl_K_scalar) && isfinite(lnD_Ca_scalar)
    P_high_scalar = ...
        (-T_AlOl_K_scalar .* lnD_Ca_scalar ...
         - 11982 + 3.61 .* T_AlOl_K_scalar) ./ 56.2;

    P_low_scalar = ...
        (-T_AlOl_K_scalar .* lnD_Ca_scalar ...
         - 5792 - 1.25 .* T_AlOl_K_scalar) ./ 42.5;
else
    P_high_scalar = NaN;
    P_low_scalar = NaN;
end

boundary_high_K_scalar = 1275.25 + 2.827 .* P_high_scalar;
boundary_low_K_scalar = 1275.25 + 2.827 .* P_low_scalar;
branchResidual_high_K_scalar = ...
    T_AlOl_K_scalar - boundary_high_K_scalar;
branchResidual_low_K_scalar = ...
    T_AlOl_K_scalar - boundary_low_K_scalar;

highCondition_scalar = isfinite(branchResidual_high_K_scalar) && ...
    branchResidual_high_K_scalar >= 0;
lowCondition_scalar = isfinite(branchResidual_low_K_scalar) && ...
    branchResidual_low_K_scalar <= 0;

if highCondition_scalar && ~lowCondition_scalar
    P_kbar_scalar = P_high_scalar;
    branch_scalar = "high-T";
    branchStatus_scalar = "unique";
elseif lowCondition_scalar && ~highCondition_scalar
    P_kbar_scalar = P_low_scalar;
    branch_scalar = "low-T";
    branchStatus_scalar = "unique";
elseif highCondition_scalar && lowCondition_scalar
    if abs(branchResidual_high_K_scalar) <= ...
            abs(branchResidual_low_K_scalar)
        P_kbar_scalar = P_high_scalar;
        branch_scalar = "high-T";
    else
        P_kbar_scalar = P_low_scalar;
        branch_scalar = "low-T";
    end
    branchStatus_scalar = "both-valid";
else
    if isfinite(branchResidual_high_K_scalar) && ...
            (~isfinite(branchResidual_low_K_scalar) || ...
             abs(branchResidual_high_K_scalar) <= ...
             abs(branchResidual_low_K_scalar))
        P_kbar_scalar = P_high_scalar;
        branch_scalar = "high-T";
    elseif isfinite(branchResidual_low_K_scalar)
        P_kbar_scalar = P_low_scalar;
        branch_scalar = "low-T";
    else
        P_kbar_scalar = NaN;
        branch_scalar = "undetermined";
    end
    branchStatus_scalar = "neither-valid";
end

P_GPa_scalar = P_kbar_scalar ./ 10;

% --- Expand composition-dependent values to launcher-vector length ---
row = table();
row.T_launcher_degreeC = T_launcher_degreeC;
row.T_launcher_K = T_launcher_K;
row.T_degreeC = repmat(T_AlOl_degreeC_scalar, nT, 1);
row.T_K = repmat(T_AlOl_K_scalar, nT, 1);
row.T_AlOl_degreeC = repmat(T_AlOl_degreeC_scalar, nT, 1);
row.T_AlOl_K = repmat(T_AlOl_K_scalar, nT, 1);
row.isLauncherTemperatureUsed = false(nT, 1);

row.Al2O3_ol_wtpercent = ...
    repmat(inputs.Al2O3_ol_wtpercent, nT, 1);
row.Al2O3_sp_wtpercent = ...
    repmat(inputs.Al2O3_sp_wtpercent, nT, 1);
row.Al2O3_ol_source = repmat(string(inputs.Al2O3_ol_source), nT, 1);
row.Al2O3_sp_source = repmat(string(inputs.Al2O3_sp_source), nT, 1);
row.D_Al2O3_ol_sp = repmat(D_Al2O3_scalar, nT, 1);
row.lnD_Al2O3_ol_sp = repmat(lnD_Al2O3_scalar, nT, 1);

row.Cr_sp_4O = repmat(inputs.Cr_sp_4O, nT, 1);
row.Al_sp_4O = repmat(inputs.Al_sp_4O, nT, 1);
row.Fe3_sp_4O = repmat(inputs.Fe3_sp_4O, nT, 1);
row.FeTotal_sp_4O = repmat(inputs.FeTotal_sp_4O, nT, 1);
row.Fe3_sp_source = repmat(string(inputs.Fe3_sp_source), nT, 1);
row.spinelCationSource = ...
    repmat(string(inputs.spinelCationSource), nT, 1);
row.CrNumber_sp_DSouza = ...
    repmat(CrNumber_DSouza_scalar, nT, 1);
row.CrNumber_sp_Coogan = ...
    repmat(CrNumber_Coogan_scalar, nT, 1);
row.T_AlOl_denominator = repmat(T_denominator_scalar, nT, 1);

row.Ca_ol_4O = repmat(inputs.Ca_ol_4O, nT, 1);
row.Ca_cpx_6O = repmat(inputs.Ca_cpx_6O, nT, 1);
row.Ca_ol_source = repmat(string(inputs.Ca_ol_source), nT, 1);
row.Ca_cpx_source = repmat(string(inputs.Ca_cpx_source), nT, 1);
row.D_Ca_ol_cpx = repmat(D_Ca_scalar, nT, 1);
row.lnD_Ca_ol_cpx = repmat(lnD_Ca_scalar, nT, 1);

row.P_KB90_high_kbar = repmat(P_high_scalar, nT, 1);
row.P_KB90_low_kbar = repmat(P_low_scalar, nT, 1);
row.KB90_boundary_high_K = repmat(boundary_high_K_scalar, nT, 1);
row.KB90_boundary_low_K = repmat(boundary_low_K_scalar, nT, 1);
row.KB90_branchResidual_high_K = ...
    repmat(branchResidual_high_K_scalar, nT, 1);
row.KB90_branchResidual_low_K = ...
    repmat(branchResidual_low_K_scalar, nT, 1);
row.KB90_highConditionSatisfied = ...
    repmat(highCondition_scalar, nT, 1);
row.KB90_lowConditionSatisfied = ...
    repmat(lowCondition_scalar, nT, 1);
row.KB90_branch = repmat(branch_scalar, nT, 1);
row.KB90_branchStatus = repmat(branchStatus_scalar, nT, 1);

row.P_kbar = repmat(P_kbar_scalar, nT, 1);
row.P_GPa = repmat(P_GPa_scalar, nT, 1);
row.P_DSouza2020_kbar = row.P_kbar;
row.P_DSouza2020_GPa = row.P_GPa;
row.P_uncertainty_1sigma_kbar = repmat(1.7, nT, 1);

row.isWithinDSouzaExperimentTRange = ...
    isfinite(row.T_AlOl_degreeC) & ...
    row.T_AlOl_degreeC >= 950 & row.T_AlOl_degreeC <= 1250;
row.isWithinDSouzaExperimentPRange = ...
    isfinite(row.P_kbar) & row.P_kbar >= 15 & row.P_kbar <= 24;
row.isWithinKB90TRange = ...
    isfinite(row.T_AlOl_degreeC) & ...
    row.T_AlOl_degreeC >= 900 & row.T_AlOl_degreeC <= 1400;
row.isWithinKB90PRange = ...
    isfinite(row.P_kbar) & row.P_kbar >= 2 & row.P_kbar <= 60;
row.isWithinCooganCrNumberUpperRange = ...
    isfinite(row.CrNumber_sp_DSouza) & ...
    row.CrNumber_sp_DSouza <= 0.69;
row.isAlPartitionFinitePositive = ...
    isfinite(row.D_Al2O3_ol_sp) & row.D_Al2O3_ol_sp > 0;
row.isCaPartitionFinitePositive = ...
    isfinite(row.D_Ca_ol_cpx) & row.D_Ca_ol_cpx > 0;
row.isKB90BranchUnique = row.KB90_branchStatus == "unique";

end

function sp = prepareSpinelCations(data_sp)
% Obtain spinel Cr, Al, Fe total, and Fe3+ on a 4-O basis.

if height(data_sp) ~= 1
    error('Spinel input must be a 1-row table.');
end

fallback = calculateCationsFromOxides(data_sp, 4);

[Cr_direct, Cr_found] = getCationColumn(data_sp, 'Cr');
[Al_direct, Al_found] = getCationColumn(data_sp, 'Al');
[Fe3_direct, Fe3_found] = getFe3CationColumn(data_sp);
[FeTotal_direct, FeTotal_found] = getTotalFeCationColumn(data_sp, Fe3_direct, Fe3_found);

if Cr_found
    Cr = Cr_direct;
else
    Cr = fallback.Cr;
end
if Al_found
    Al = Al_direct;
else
    Al = fallback.Al;
end
if FeTotal_found
    FeTotal = FeTotal_direct;
else
    FeTotal = fallback.FeTotal;
end

if Cr_found && Al_found && FeTotal_found
    cationSource = "cation_apfu columns";
else
    cationSource = "mixed cation_apfu and oxide fallback";
end

if Fe3_found
    Fe3 = Fe3_direct;
    Fe3Source = "Fe3 cation column";
elseif fallback.hasExplicitFe2O3
    Fe3 = fallback.Fe3;
    Fe3Source = "Fe2O3 oxide normalization";
else
    c = mergeSpinelCationsWithFallback(data_sp, fallback);
    nonFeCharge = ...
        4 .* c.Si + 4 .* c.Ti + 3 .* c.Al + 3 .* c.Cr ...
        + 3 .* c.V + 2 .* c.Mn + 2 .* c.Mg + 2 .* c.Ca ...
        + 2 .* c.Ni + 2 .* c.Zn + c.Na + c.K + 5 .* c.P;

    Fe3_charge = 8 - nonFeCharge - 2 .* FeTotal;
    tolerance = 1e-6;

    if isfinite(Fe3_charge) && isfinite(FeTotal) && ...
            Fe3_charge >= -tolerance && ...
            Fe3_charge <= FeTotal + tolerance
        Fe3 = min(max(Fe3_charge, 0), FeTotal);
        Fe3Source = "stoichiometric charge balance";
    else
        Fe3 = NaN;
        Fe3Source = "invalid stoichiometric charge balance";
    end
end

sp = struct();
sp.Cr = Cr;
sp.Al = Al;
sp.FeTotal = FeTotal;
sp.Fe3 = Fe3;
sp.Fe3Source = Fe3Source;
sp.CationSource = cationSource;

end

function c = mergeSpinelCationsWithFallback(data_sp, fallback)
% Prefer supplied cation columns and otherwise use oxide-normalized values.

elements = {'Si', 'Ti', 'Al', 'Cr', 'V', 'Mn', 'Mg', ...
    'Ca', 'Ni', 'Zn', 'Na', 'K', 'P'};

c = struct();
for i = 1:numel(elements)
    element = elements{i};
    [directValue, found] = getCationColumn(data_sp, element);
    if found
        c.(element) = directValue;
    else
        c.(element) = fallback.(element);
    end
end

end

function [value, source] = getCationWithOxideFallback( ...
        data_m, element, oxygenBasis, mineralLabel)
% Prefer phase-normalized cation columns; otherwise calculate from oxides.

[value, found] = getCationColumn(data_m, element);
if found
    source = [element '_cation_apfu column'];
    return;
end

cations = calculateCationsFromOxides(data_m, oxygenBasis);
value = cations.(element);
source = sprintf('oxide normalization to %.0f O (%s)', ...
    oxygenBasis, mineralLabel);

end

function [value, source] = getAl2O3WtPercent(data_m, mineralLabel)
% Retrieve Al2O3 wt%, with optional conversion from elemental Al.

[value, found] = getTableScalarByAliases(data_m, ...
    {'Al2O3', 'Al2O3_wtpercent', 'Al2O3_wt', 'Al2O3_percent'});
if found
    source = [mineralLabel '.Al2O3 wt%'];
    return;
end

[Al_ppm, foundPpm] = getTableScalarByAliases(data_m, ...
    {'Al_ppm', 'AlElement_ppm', 'Al_ppm_element'});
if foundPpm
    atomicWeightAl = 26.9815385;
    molarWeightAl2O3 = 101.961276;
    value = (Al_ppm ./ 10000) .* ...
        (molarWeightAl2O3 ./ (2 .* atomicWeightAl));
    source = [mineralLabel '.Al ppm converted to Al2O3 wt%'];
    return;
end

[Al_wtpercent, foundWt] = getTableScalarByAliases(data_m, ...
    {'Al_wtpercent', 'Al_wt', 'Al_percent'});
if foundWt
    atomicWeightAl = 26.9815385;
    molarWeightAl2O3 = 101.961276;
    value = Al_wtpercent .* ...
        (molarWeightAl2O3 ./ (2 .* atomicWeightAl));
    source = [mineralLabel '.Al wt% converted to Al2O3 wt%'];
    return;
end

value = NaN;
source = [mineralLabel '.Al2O3 unavailable'];

end

function c = calculateCationsFromOxides(data_m, oxygenBasis)
% Calculate selected mineral cations from oxide wt% on a specified O basis.
% Missing oxide columns are treated as zero only inside this fallback.

if height(data_m) ~= 1
    error('Mineral oxide input must be a 1-row table.');
end

[FeO, hasFeO] = getTableScalarByAliases(data_m, {'FeO'});
[FeOt, hasFeOt] = getTableScalarByAliases(data_m, {'FeOt', 'FeOtotal'});
[Fe2O3, hasFe2O3] = getTableScalarByAliases(data_m, {'Fe2O3'});

if ~hasFeO && hasFeOt
    FeO = FeOt;
    hasFeO = true;
end
if ~hasFeO
    FeO = 0;
end
if ~hasFe2O3
    Fe2O3 = 0;
end

oxideNames = {'SiO2', 'TiO2', 'Al2O3', 'MnO', 'MgO', 'CaO', ...
    'Na2O', 'K2O', 'Cr2O3', 'V2O3', 'NiO', 'ZnO', 'P2O5'};
molarWeights = [60.0843, 79.866, 101.961276, 70.93745, ...
    40.3044, 56.0774, 61.97894, 94.196, 151.9904, ...
    149.881, 74.6928, 81.379, 141.9445];
cationNumbers = [1, 1, 2, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2];
oxygenNumbers = [2, 2, 3, 1, 1, 1, 1, 1, 3, 3, 1, 1, 5];
elementNames = {'Si', 'Ti', 'Al', 'Mn', 'Mg', 'Ca', ...
    'Na', 'K', 'Cr', 'V', 'Ni', 'Zn', 'P'};

cationProportions = zeros(size(oxideNames));
oxygenProportions = zeros(size(oxideNames));

for i = 1:numel(oxideNames)
    [oxideValue, found] = getTableScalarByAliases(data_m, oxideNames(i));
    if ~found
        oxideValue = 0;
    end
    cationProportions(i) = ...
        oxideValue .* cationNumbers(i) ./ molarWeights(i);
    oxygenProportions(i) = ...
        oxideValue .* oxygenNumbers(i) ./ molarWeights(i);
end

Fe2Proportion = FeO ./ 71.844;
Fe2OxygenProportion = FeO ./ 71.844;
Fe3Proportion = Fe2O3 .* 2 ./ 159.688;
Fe3OxygenProportion = Fe2O3 .* 3 ./ 159.688;

oxygenTotal = sum(oxygenProportions) + ...
    Fe2OxygenProportion + Fe3OxygenProportion;

c = struct();
if ~isfinite(oxygenTotal) || oxygenTotal <= 0
    normalizationFactor = NaN;
else
    normalizationFactor = oxygenBasis ./ oxygenTotal;
end

for i = 1:numel(elementNames)
    c.(elementNames{i}) = ...
        cationProportions(i) .* normalizationFactor;
end

c.Fe2 = Fe2Proportion .* normalizationFactor;
c.Fe3 = Fe3Proportion .* normalizationFactor;
c.FeTotal = c.Fe2 + c.Fe3;
c.normalizationFactor = normalizationFactor;
c.hasExplicitFe2O3 = hasFe2O3 && isfinite(Fe2O3);

end

function [value, found] = getCationColumn(data_m, element)
% Retrieve a phase-normalized cation value without changing NaN.

aliases = { ...
    [element '_cation_apfu'], ...
    [element '_cation'], ...
    [element '_apfu']};
[value, found] = getTableScalarByAliases(data_m, aliases);

end

function [value, found] = getTotalFeCationColumn(data_m, Fe3, foundFe3)
% Retrieve total Fe cations when explicitly available. Explicit total-Fe
% columns are preferred. Otherwise Fe2 (including the project-standard
% Fe_cation_apfu alias) is combined with Fe3 when Fe3 is available.

[value, found] = getTableScalarByAliases(data_m, ...
    {'Fet_cation_apfu', 'FeTotal_cation_apfu', ...
     'Fet_cation', 'FeTotal_cation'});

if found
    return;
end

[Fe2, foundFe2] = getTableScalarByAliases(data_m, ...
    {'Fe2_cation_apfu', 'Fe_cation_apfu', ...
     'Fe2_cation', 'Fe_cation'});

if foundFe2
    if foundFe3
        value = Fe2 + Fe3;
    else
        value = Fe2;
    end
    found = true;
end

end

function [value, found] = getFe3CationColumn(data_m)
% Retrieve Fe3+ cations when explicitly available.

[value, found] = getTableScalarByAliases(data_m, ...
    {'Fe3_cation_apfu', 'Fe3plus_cation_apfu', ...
     'FeIII_cation_apfu', 'Fe3_cation'});

end

function nanInputNames = findNaNInputs(inputs)
% Return exact calculation-input names containing NaN.

fields = { ...
    'Ca_ol_4O', 'Ca_cpx_6O', ...
    'Al2O3_ol_wtpercent', 'Al2O3_sp_wtpercent', ...
    'Cr_sp_4O', 'Al_sp_4O', 'Fe3_sp_4O'};
labels = { ...
    'Olivine.Ca(4O)', 'Cpx.Ca(6O)', ...
    'Olivine.Al2O3', 'Spinel.Al2O3', ...
    'Spinel.Cr(4O)', 'Spinel.Al(4O)', 'Spinel.Fe3(4O)'};

buffer = strings(numel(fields), 1);
nFound = 0;
for i = 1:numel(fields)
    if isnan(inputs.(fields{i}))
        nFound = nFound + 1;
        buffer(nFound) = labels{i};
    end
end
nanInputNames = buffer(1:nFound);

end

function validateNonNegativeInputs(inputs)
% Reject Inf and finite negative calculation inputs. NaN and zero remain.

fields = { ...
    'Ca_ol_4O', 'Ca_cpx_6O', ...
    'Al2O3_ol_wtpercent', 'Al2O3_sp_wtpercent', ...
    'Cr_sp_4O', 'Al_sp_4O', 'Fe3_sp_4O', 'FeTotal_sp_4O'};
labels = { ...
    'Olivine.Ca(4O)', 'Cpx.Ca(6O)', ...
    'Olivine.Al2O3', 'Spinel.Al2O3', ...
    'Spinel.Cr(4O)', 'Spinel.Al(4O)', ...
    'Spinel.Fe3(4O)', 'Spinel.FeTotal(4O)'};

invalidBuffer = strings(numel(fields), 1);
nInvalid = 0;
for i = 1:numel(fields)
    value = inputs.(fields{i});
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalid = nInvalid + 1;
        invalidBuffer(nInvalid) = labels{i};
    end
end

if nInvalid > 0
    invalidNames = invalidBuffer(1:nInvalid);
    error(['DSouza2020 calculation inputs must be non-negative. NaN is ' ...
           'allowed, but Inf and finite negative values are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidNames, ', ')) '.']);
end

end

function [dataset, resolvedField] = resolveMineralTable( ...
        rawdata_struct, primaryField, aliases)
% Resolve a mineral table from accepted field names.

candidateFields = [{primaryField}, aliases];
structFields = fieldnames(rawdata_struct);
canonicalStructFields = cellfun(@canonicalizeName, structFields, ...
    'UniformOutput', false);

for i = 1:numel(candidateFields)
    target = canonicalizeName(candidateFields{i});
    idx = find(strcmp(canonicalStructFields, target), 1, 'first');
    if ~isempty(idx)
        resolvedField = structFields{idx};
        dataset = rawdata_struct.(resolvedField);
        if ~istable(dataset)
            error('rawdata_struct.%s must be a table.', resolvedField);
        end
        if isempty(dataset)
            error('rawdata_struct.%s is empty.', resolvedField);
        end
        return;
    end
end

error('rawdata_struct must contain a %s table.', primaryField);

end

function validateRequestedRow(rowNumber, numberOfRows, optionName)
% Confirm that a requested row lies inside its mineral table.

if ~isempty(rowNumber) && rowNumber > numberOfRows
    error('%s (%d) exceeds the number of available rows (%d).', ...
        optionName, rowNumber, numberOfRows);
end

end

function [selectedIdx, ok] = selectMineralRow(items, rowOption, label)
% Select a mineral row interactively or use a supplied fixed row.

if isempty(rowOption)
    [selectedIdx, ok] = listdlg( ...
        'PromptString', ['Please select the ' label ' data you would like to use:'], ...
        'SelectionMode', 'single', ...
        'ListString', items, ...
        'ListSize', [520 360]);
else
    selectedIdx = rowOption;
    ok = true;
end

end

function items = buildMineralList(dataset, mineralLabel)
% Build display labels for a mineral selection dialog.

nRows = height(dataset);
items = cell(nRows, 1);
for i = 1:nRows
    code = getDataCode(dataset(i, :), i);
    items{i} = [mineralLabel ' Row ' num2str(i) ' | ' code];
end

end

function code = getDataCode(data_m, rowNumber)
% Return the first-column value as a compact data-code label.

if width(data_m) < 1
    code = ['Row ' num2str(rowNumber)];
    return;
end

rawValue = data_m{1, 1};
if iscell(rawValue)
    rawValue = rawValue{1};
end

if isempty(rawValue)
    code = ['Row ' num2str(rowNumber)];
else
    code = char(string(rawValue));
end

end

function [value, found] = getTableScalarByAliases(data_m, aliases)
% Retrieve the first scalar value from the first matching variable name.
% Matching ignores spaces, underscores, and hyphens and accepts a trailing
% "Value" suffix. If a column exists but contains missing/non-numeric data,
% found remains true and value is NaN so the missing value is preserved.

variableNames = data_m.Properties.VariableNames;
canonicalVariables = cellfun(@canonicalizeName, variableNames, ...
    'UniformOutput', false);

value = NaN;
found = false;

for i = 1:numel(aliases)
    target = canonicalizeName(aliases{i});
    targets = {target, [target 'value']};
    idx = [];
    for j = 1:numel(targets)
        idx = find(strcmp(canonicalVariables, targets{j}), 1, 'first');
        if ~isempty(idx)
            break;
        end
    end

    if ~isempty(idx)
        found = true;
        rawValue = data_m.(variableNames{idx});
        value = toScalarDouble(rawValue);
        return;
    end
end

end

function name = canonicalizeName(name)
% Convert a name to a comparison-safe lowercase character vector.

name = lower(char(string(name)));
name = strrep(name, ' ', '');
name = strrep(name, '_', '');
name = strrep(name, '-', '');
name = strrep(name, '+', 'plus');

end

function value = toScalarDouble(rawValue)
% Convert the first table value to double; missing values remain NaN.

value = NaN;
if isempty(rawValue)
    return;
end

if istable(rawValue)
    rawValue = rawValue{1, 1};
end
if iscell(rawValue)
    if isempty(rawValue{1})
        return;
    end
    rawValue = rawValue{1};
end

if isnumeric(rawValue) || islogical(rawValue)
    value = double(rawValue(1));
elseif isstring(rawValue)
    if ~ismissing(rawValue(1))
        value = str2double(rawValue(1));
    end
elseif ischar(rawValue)
    value = str2double(string(rawValue));
elseif iscategorical(rawValue)
    value = str2double(string(rawValue(1)));
end

end
