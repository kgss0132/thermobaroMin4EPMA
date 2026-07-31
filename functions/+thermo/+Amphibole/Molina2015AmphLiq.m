function results = Molina2015AmphLiq(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Amphibole_Liquid/Molina2015AmphLiq.m
% Tested with MATLAB R2024b
%
% Calcic amphibole-liquid Mg-partition thermometer
% Molina, J.F., Moreno, J.A., Castro, A., Rodriguez, C., Fershtater, G.B. (2015)
% Lithos, 232, 286-305
% DOI: https://doi.org/10.1016/j.lithos.2015.06.027
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Amphibole analysis and one Liquid
% analysis and calculates temperature using Mg partitioning between calcic
% amphibole and silicate liquid.
%
% The function accepts pressure as either a scalar or a vector so that it can
% be called from both startThermoCalc_fixedP and startThermoCalc_rangeP.
% Pressure is retained in the output table for interface consistency and for
% checking the pressure range represented by the experimental database, but
% pressure is not used in the published thermometer equation.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Amphibole-Liquid pair, and appends the
% results into a single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Molina et al. (2015) calibrated this thermometer using selected
% amphibole-glass pairs compiled from experiments in the literature.
%
%   Temperature : 800-1100 degreeC
%   Pressure    : 0.8-20 kbar represented in the selected database
%   Amphibole   : calcic amphibole with
%                 #CaM4 = CaM4 / (CaM4 + NaM4) > 0.75
%   Formula     : amphibole normalized by the 13-CNK method on a
%                 23-oxygen basis
%   Liquid      : anhydrous mole fractions of cation components
%   Rock types  : alkaline and subalkaline igneous rocks
%
% The experimental database and data-selection procedure are described on
% pp. 287-289 (Section 2; Tables 2-3). The amphibole-liquid thermometer is
% calibrated and tested on pp. 296-298 (Section 6.2; Tables 6-7). Explicit
% application restrictions are stated on p. 298 (Section 7.3), and are
% repeated in the Conclusions on p. 304. Compositional screening against the
% experimental data fields is illustrated in Fig. 14 on p. 302.
%
% Molina et al. (2015) report a precision of approximately +/-35 degreeC for
% the calibration dataset and +/-45 degreeC for the independent test dataset
% (pp. 296-298; Table 7). Natural samples may have additional uncertainty
% caused by analytical error, disequilibrium, zoning, melt evolution, or an
% inappropriate pairing of amphibole and liquid compositions.
%
% IMPORTANT APPLICATION NOTES:
%   1) The selected amphibole and liquid must represent an equilibrium pair.
%      Pairing an amphibole core with a later matrix glass, a xenocryst with
%      an unrelated host melt, or a cumulate amphibole with an uncorrected
%      whole-rock composition can produce misleading temperatures
%      (discussion on pp. 301-304).
%   2) More reliable estimates are expected when amphibole and liquid
%      compositions lie within the experimental compositional fields shown
%      in Fig. 14 (p. 302).
%   3) The 0.8-20 kbar interval is the pressure range represented by the
%      selected amphibole-glass database (Table 3, p. 289). Pressure was not
%      statistically significant in the final regression and does not occur
%      in the thermometer equation (p. 296). A pressure-range warning is
%      therefore an extrapolation warning, not a pressure correction.
%   4) The equation is intended for calcic amphibole. Results for
%      #CaM4 <= 0.75 are outside the stated compositional applicability.
%   5) Kaersutite and other unusual high-Ti amphiboles require particular
%      care because the 13-CNK normalization may be less reliable for some
%      amphibole compositions (normalization discussion on p. 292).
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) input pressure is outside 0.8-20 kbar,
%   2) a finite calculated temperature is outside 800-1100 degreeC,
%   3) finite #CaM4 is <= 0.75,
%   4) a required input contains NaN, or
%   5) the calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole : table
%
% The FIRST column of the Amphibole table is treated as an identifier
% ("data code") displayed in the selection dialog.
%
% The liquid dataset is loaded with:
%   liquid.readLiquidExcel()
%
% Amphibole compositions are preferably read from oxide wt% columns and are
% normalized internally by the 13-CNK method. Common oxide column variants
% such as 'SiO2', 'SiO2value', 'SiO2 value', and 'SiO2_value' are accepted.
% If no SiO2 oxide column is available, the function attempts to use existing
% 23-oxygen or apfu cation columns.
%
% Liquid compositions are read as oxide wt% and converted to anhydrous mole
% fractions of cation components. The calculation uses:
%   SiO2, TiO2, Al2O3, FeO, Fe2O3, MnO, MgO, CaO, Na2O, K2O,
%   Cr2O3, NiO, and P2O5
%
% Missing optional oxide columns are treated as absent components (0), but
% an explicit NaN stored in an existing input column is retained as NaN. NaN
% is never replaced by zero. NaN values propagate through the calculation,
% remain in the output table, and trigger non-stopping fprintf warnings.
%
% All finite composition values used by this thermometer must be
% non-negative. Negative finite values stop the calculation with an error.
% Zero is permitted as an input value, but it may make a logarithm or ratio
% mathematically invalid; in that case the calculated temperature is kept as
% NaN and a warning is printed.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% 1) Amphibole Mg mole fraction
%     XMg_amp = Mg23O_amp / SumCat_amp
%     SumCat_amp = 16
%
% 2) Liquid cation mole fractions
%     XN_liq = N_liq / SumCat_liq
%
%     SumCat_liq = Si + Ti + Al + Cr + Ni + Fe + Mn + Mg + Ca
%                  + Na + K + P
%
% 3) Amphibole-liquid Mg partition coefficient
%     DMg_amp_liq = XMg_amp / XMg_liq
%
% 4) Temperature solution
%     T(degreeC) =
%       [71975 - 11896*ln(XCa_liq / (XCa_liq + XAl_liq))]
%       / [8.3144*ln(DMg_amp_liq) + 58] - 273
%
% Notes:
% - Natural logarithms are used.
% - The published equation subtracts 273, and that value is retained here
%   exactly for reproducibility.
% - Pressure is not used in the temperature equation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Molina2015AmphLiq(rawdata_struct, P_kbar)
%   results = Molina2015AmphLiq(rawdata_struct, P_kbar, ...
%       'UseFirstLiquidRow', false)
%
% Inputs:
%   rawdata_struct : struct containing an Amphibole table (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Name-value option:
%   'UseFirstLiquidRow' : logical scalar, default true
%       true  -> use row 1 of the selected liquid dataset automatically
%       false -> show a list dialog for liquid-row selection
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Amphibole-Liquid pair
%

%% Input validation
% Basic argument checks prevent silent failures caused by missing inputs or
% invalid pressure values.
if nargin < 2
    error('Molina2015AmphLiq requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end
if ~isfield(rawdata_struct, 'Amphibole') || ...
        ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end

P_kbar = P_kbar(:);

%% Options
% Keep the existing liquid-row selection option while validating its type.
ip = inputParser;
ip.addParameter('UseFirstLiquidRow', true, ...
    @(x) islogical(x) && isscalar(x));
ip.parse(varargin{:});
useFirstLiquidRow = ip.Results.UseFirstLiquidRow;

%% 1) Retrieve datasets and constants
% Read the Amphibole table, molecular weights, and the external liquid file.
disp('=== Step 1: Preparing Amphibole and Liquid datasets ===');

dataset_amp = rawdata_struct.Amphibole;
MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();

if isempty(liqAll) || ~istable(liqAll)
    error('Selected liquid dataset is empty or is not a table.');
end

disp('=== Preparing Amphibole and Liquid datasets has been finished ===');

%% 2) Initialize output container and calibration limits
% Each calculation result is stored temporarily as one table block.
% Repeated table concatenation inside the interactive loop is avoided because
% it forces MATLAB to repeatedly reallocate and copy the entire results table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Application and experimental-data limits reported by Molina et al. (2015).
calibrationT_min_degC = 800;
calibrationT_max_degC = 1100;
experimentalP_min_kbar = 0.8;
experimentalP_max_kbar = 20;
calcicCriterion_min = 0.75;

% Pressure is common to all selected pairs in this function call, so the
% pressure warning is printed only once after the first calculation.
pressureOutsideCalibration = ...
    P_kbar < experimentalP_min_kbar | ...
    P_kbar > experimentalP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop and calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % ----- Amphibole selection -----
    % Assumption: the first column stores an identifier displayed to the user.
    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, okAmp] = listdlg( ...
        'PromptString', ...
        'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_amp, ...
        'ListSize', [420 360]);

    if ~okAmp || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    selectedData_amp = dataset_amp(selectedIdx_amp, :);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

    % ----- Liquid selection -----
    disp('=== Step 4: Selecting a data row from the list (Liquid) ===');

    if useFirstLiquidRow
        if height(liqAll) > 1
            fprintf(2, ...
                ['WARNING: Liquid dataset has multiple rows (%d). ' ...
                 'Row 1 is being used automatically.\n'], ...
                height(liqAll));
        end
        selectedIdx_liq = 1;
        disp('Liquid selected: (auto) Row 1');
    else
        liquidItems = buildLiquidList(liqAll);
        [selectedIdx_liq, okLiq] = listdlg( ...
            'PromptString', ...
            'Please select the Liquid data you would like to use:', ...
            'SelectionMode', 'single', ...
            'ListString', liquidItems, ...
            'ListSize', [520 360]);

        if ~okLiq || isempty(selectedIdx_liq)
            disp('Selection canceled');
            break;
        end

        disp(['Liquid selected: Row ' num2str(selectedIdx_liq)]);
    end

    selectedData_liq = liqAll(selectedIdx_liq, :);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    % Check only values actually recognized and used by this implementation.
    % Explicit NaN values are recorded for later warning but do not stop the
    % calculation. Negative finite values are not physically meaningful and
    % stop the calculation before normalization or logarithms are evaluated.
    nanInputNames = findNaNInputs(selectedData_amp, selectedData_liq);
    validateNonNegativeInputs(selectedData_amp, selectedData_liq);

    row = calcTemp(selectedData_amp, selectedData_liq, P_kbar, MWinfo);

    % Store user-selected identifiers for traceability. Identifiers and liquid
    % metadata are replicated once for each supplied pressure value.
    nRows = height(row);
    row.dataCode_amp = repmat(string(selectedCode_amp), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_amp','dataRow_liq'}, 'Before', 1);

    % Store this result as one table block. If the preallocated cell buffer is
    % full, double its capacity; this happens only occasionally rather than on
    % every loop iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature immediately after each selection.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    liquidLabel = makeLiquidLabel(selectedData_liq, selectedIdx_liq);

    if height(row) == 1
        disp([char(string(selectedCode_amp)) ' & ' liquidLabel ': ' ...
            num2str(row.T_Molina2015_C) ' degreeC']);
    else
        disp([char(string(selectedCode_amp)) ' & ' liquidLabel ': ' ...
            num2str(row.T_Molina2015_C(1)) ' to ' ...
            num2str(row.T_Molina2015_C(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside the 0.8-20 kbar range
    % represented by the selected experimental amphibole-glass database.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the pressure range represented ' ...
             'by the selected experimental database of Molina et al. (2015): ' ...
             '0.8-20 kbar. %d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar. Pressure is not used in the ' ...
             'thermometer equation, so this is an extrapolation warning only.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside 800-1100 C.
    finiteTemperature = isfinite(row.T_Molina2015_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_Molina2015_C < calibrationT_min_degC | ...
         row.T_Molina2015_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_Molina2015_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the calibration range ' ...
             'of Molina et al. (2015): 800-1100 degreeC. %d of %d finite ' ...
             'temperature point(s) are outside the range; calculated finite ' ...
             'range = %.4g-%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_amp)), ...
            liquidLabel);
    end

    % Warn when the selected amphibole is outside the stated calcic-amphibole
    % criterion. NaN is handled by the separate NaN and non-finite warnings.
    finiteCaM4 = isfinite(row.CaM4_ratio);
    outsideCalcicCriterion = finiteCaM4 & ...
        row.CaM4_ratio <= calcicCriterion_min;

    if any(outsideCalcicCriterion)
        fprintf(2, ...
            ['WARNING: Amphibole composition is outside the stated applicability ' ...
             'of Molina et al. (2015): #CaM4 must be > 0.75. ' ...
             'Calculated #CaM4 = %.6g for %s.\n'], ...
            row.CaM4_ratio(find(finiteCaM4, 1, 'first')), ...
            char(string(selectedCode_amp)));
    end

    % Print a non-stopping warning when any recognized input contains NaN.
    % fprintf is used instead of warning so the message remains in the log even
    % when MATLAB warnings have been disabled globally.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and the calculation was continued; ' ...
             'the calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_amp)), ...
            liquidLabel, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain result-based checks for NaN/Inf generated by missing data, zero
    % denominator, invalid logarithm arguments, or other numerical problems.
    invalidTemperature = ~isfinite(row.T_Molina2015_C);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_amp)), ...
            liquidLabel, ...
            sum(invalidTemperature), ...
            numel(row.T_Molina2015_C), ...
            sum(isnan(row.T_Molina2015_C)), ...
            sum(isinf(row.T_Molina2015_C)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat using the same loaded liquid dataset.
    userAction = questdlg( ...
        'Continue with another data selection (same Liquid dataset)?', ...
        'Molina2015AmphLiq', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate the buffered table blocks once after all selections have been
% completed. Return an empty table when no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

% Retain liquid-file metadata for downstream inspection.
results.Properties.UserData = struct('liquid', metaLiq);

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function row = calcTemp(data_amp, data_liq, P_kbar, MWinfo)
% calcTemp
% Calculate temperature for one Amphibole row and one Liquid row.
% One output row is returned for each input pressure value.
%
% Pressure does not occur in the published Molina et al. (2015)
% thermometer equation. Therefore, the calculated temperature is identical
% for every pressure point, while P_kbar is retained for traceability and
% calibration-range checks.
%
% Inputs:
%   data_amp : one-row Amphibole table
%   data_liq : one-row Liquid table
%   P_kbar   : finite non-negative pressure scalar or vector in kbar
%   MWinfo   : molecular-weight and cation-number structure
%
% Output:
%   row : table containing one row for each supplied pressure value

P_kbar = P_kbar(:);
nPressure = numel(P_kbar);

row = table();

% Store pressure first for compatibility with fixed-P and range-P workflows.
row.P_kbar = P_kbar;

% -------------------------------------------------------------------------
% Amphibole: 13-CNK normalization
% -------------------------------------------------------------------------
amp = prepareAmphibole13CNK(data_amp, MWinfo);

row.Si_amp_23O        = repmat(amp.Si,          nPressure, 1);
row.Ti_amp_23O        = repmat(amp.Ti,          nPressure, 1);
row.Al_amp_23O        = repmat(amp.Al,          nPressure, 1);
row.AlIV_amp_23O      = repmat(amp.AlIV,        nPressure, 1);
row.AlVI_amp_23O      = repmat(amp.AlVI,        nPressure, 1);
row.Fe3_amp_23O       = repmat(amp.Fe3,         nPressure, 1);
row.Fe2_amp_23O       = repmat(amp.Fe2,         nPressure, 1);
row.Mn_amp_23O        = repmat(amp.Mn,          nPressure, 1);
row.Mg_amp_23O        = repmat(amp.Mg,          nPressure, 1);
row.Ca_amp_23O        = repmat(amp.Ca,          nPressure, 1);
row.Na_amp_23O        = repmat(amp.Na,          nPressure, 1);
row.K_amp_23O         = repmat(amp.K,           nPressure, 1);
row.Cr_amp_23O        = repmat(amp.Cr,          nPressure, 1);
row.Ni_amp_23O        = repmat(amp.Ni,          nPressure, 1);
row.NaA_amp_23O       = repmat(amp.NaA,         nPressure, 1);
row.NaM4_amp_23O      = repmat(amp.NaM4,        nPressure, 1);
row.CaM4_amp_23O      = repmat(amp.CaM4,        nPressure, 1);
row.Vc_amp            = repmat(amp.Vc,          nPressure, 1);
row.SumCat_amp        = repmat(amp.SumCat,      nPressure, 1);
row.XMg_amp           = repmat(amp.XMg,         nPressure, 1);
row.CaM4_ratio        = repmat(amp.CaM4_ratio,  nPressure, 1);
row.isCalcicCriterion = row.CaM4_ratio > 0.75;

% -------------------------------------------------------------------------
% Liquid oxide compositions
% -------------------------------------------------------------------------
% Missing optional columns are treated as absent components and assigned 0.
% Existing columns containing NaN retain NaN.
%
% SiO2, Al2O3, MgO, and CaO are required by the thermometer calculation.
% If these columns are absent, NaN is assigned and propagated.

SiO2  = getLiqOxOptional(data_liq, 'SiO2',  NaN);
TiO2  = getLiqOxOptional(data_liq, 'TiO2',  0);
Al2O3 = getLiqOxOptional(data_liq, 'Al2O3', NaN);
FeO   = getLiqOxOptional(data_liq, 'FeO',   0);
Fe2O3 = getLiqOxOptional(data_liq, 'Fe2O3', 0);
MnO   = getLiqOxOptional(data_liq, 'MnO',   0);
MgO   = getLiqOxOptional(data_liq, 'MgO',   NaN);
CaO   = getLiqOxOptional(data_liq, 'CaO',   NaN);
Na2O  = getLiqOxOptional(data_liq, 'Na2O',  0);
K2O   = getLiqOxOptional(data_liq, 'K2O',   0);
Cr2O3 = getLiqOxOptional(data_liq, 'Cr2O3', 0);
NiO   = getLiqOxOptional(data_liq, 'NiO',   0);
P2O5  = getLiqOxOptional(data_liq, 'P2O5',  0);

row.SiO2_liq  = repmat(SiO2,  nPressure, 1);
row.TiO2_liq  = repmat(TiO2,  nPressure, 1);
row.Al2O3_liq = repmat(Al2O3, nPressure, 1);
row.FeO_liq   = repmat(FeO,   nPressure, 1);
row.Fe2O3_liq = repmat(Fe2O3, nPressure, 1);
row.MnO_liq   = repmat(MnO,   nPressure, 1);
row.MgO_liq   = repmat(MgO,   nPressure, 1);
row.CaO_liq   = repmat(CaO,   nPressure, 1);
row.Na2O_liq  = repmat(Na2O,  nPressure, 1);
row.K2O_liq   = repmat(K2O,   nPressure, 1);
row.Cr2O3_liq = repmat(Cr2O3, nPressure, 1);
row.NiO_liq   = repmat(NiO,   nPressure, 1);
row.P2O5_liq  = repmat(P2O5,  nPressure, 1);

% -------------------------------------------------------------------------
% Liquid cation moles on an anhydrous basis
% -------------------------------------------------------------------------
nSi = SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
nTi = TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
nAl = Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;

nFe = ...
    FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO + ...
    Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

nMn = MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
nMg = MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
nCa = CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
nNa = Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
nK  = K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
nCr = Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
nNi = NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;

% Do not use nP here because nP may be confused with the number of
% pressure points. Use a descriptive variable name for phosphorus.
nPhosphorus = ...
    P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;

sumCatLiq = ...
    nSi + nTi + nAl + nCr + nNi + nFe + ...
    nMn + nMg + nCa + nNa + nK + nPhosphorus;

% -------------------------------------------------------------------------
% Initialize derived values
% -------------------------------------------------------------------------
% NaN is retained until all required denominators and logarithm arguments
% have been confirmed to be finite and mathematically valid.

XMg_liq = NaN;
XCa_liq = NaN;
XAl_liq = NaN;

DMg = NaN;
lnDMg = NaN;
lnCaCaAl = NaN;
T_C = NaN;

% -------------------------------------------------------------------------
% Calculate liquid mole fractions and temperature
% -------------------------------------------------------------------------
if isfinite(sumCatLiq) && sumCatLiq > 0

    XMg_liq = nMg ./ sumCatLiq;
    XCa_liq = nCa ./ sumCatLiq;
    XAl_liq = nAl ./ sumCatLiq;

    DMg = amp.XMg ./ XMg_liq;

    validLogArguments = ...
        isfinite(DMg) && DMg > 0 && ...
        isfinite(XCa_liq) && XCa_liq > 0 && ...
        isfinite(XAl_liq) && XAl_liq >= 0 && ...
        isfinite(XCa_liq + XAl_liq) && ...
        (XCa_liq + XAl_liq) > 0;

    if validLogArguments

        lnDMg = log(DMg);

        lnCaCaAl = log( ...
            XCa_liq ./ (XCa_liq + XAl_liq));

        denominator = 8.3144 .* lnDMg + 58;

        if isfinite(denominator) && denominator ~= 0
            T_C = ...
                (71975 - 11896 .* lnCaCaAl) ./ denominator - 273;
        end
    end
end

% -------------------------------------------------------------------------
% Replicate pressure-independent results for each pressure point
% -------------------------------------------------------------------------
row.SumCat_liq = repmat(sumCatLiq,   nPressure, 1);
row.XMg_liq = repmat(XMg_liq,        nPressure, 1);
row.XCa_liq = repmat(XCa_liq,        nPressure, 1);
row.XAl_liq = repmat(XAl_liq,        nPressure, 1);
row.DMg_amp_liq = repmat(DMg,        nPressure, 1);
row.ln_DMg_amp_liq = repmat(lnDMg,   nPressure, 1);

row.ln_XCa_over_CaPlusAl = ...
    repmat(lnCaCaAl, nPressure, 1);

row.T_Molina2015_C = ...
    repmat(T_C, nPressure, 1);

% -------------------------------------------------------------------------
% Applicability flags
% -------------------------------------------------------------------------
row.isTemperatureWithinCalibRange = ...
    isfinite(row.T_Molina2015_C) & ...
    row.T_Molina2015_C >= 800 & ...
    row.T_Molina2015_C <= 1100;

row.isPressureWithinDataRange = ...
    row.P_kbar >= 0.8 & ...
    row.P_kbar <= 20;

row.isWithinCalibRange = ...
    row.isTemperatureWithinCalibRange & ...
    row.isPressureWithinDataRange & ...
    row.isCalcicCriterion;

end

function amp = prepareAmphibole13CNK(data_amp, MWinfo)
% prepareAmphibole13CNK
% Prefer oxide-based 13-CNK recalculation. If no SiO2 oxide column exists,
% use recognized cation columns. Explicit NaN values are preserved.

amp = emptyAmpStruct();
variableNames = data_amp.Properties.VariableNames;

hasOxideData = ~isempty(findOxideColumn(variableNames, 'SiO2'));

if hasOxideData
    amp = prepareAmphibole13CNKFromOxides(data_amp, MWinfo);
else
    amp = prepareAmphibole13CNKFromCations(data_amp);
end

end

function amp = prepareAmphibole13CNKFromOxides(data_amp, MWinfo)
% prepareAmphibole13CNKFromOxides
% Recalculate amphibole cations using the 13-CNK scheme. Missing optional
% oxide columns are treated as zero; explicit NaN values remain NaN.

SiO2  = getAmpOxOptional(data_amp, 'SiO2',  NaN);
TiO2  = getAmpOxOptional(data_amp, 'TiO2',  0);
Al2O3 = getAmpOxOptional(data_amp, 'Al2O3', NaN);
FeO   = getAmpOxOptional(data_amp, 'FeO',   0);
Fe2O3 = getAmpOxOptional(data_amp, 'Fe2O3', 0);
MnO   = getAmpOxOptional(data_amp, 'MnO',   0);
MgO   = getAmpOxOptional(data_amp, 'MgO',   NaN);
CaO   = getAmpOxOptional(data_amp, 'CaO',   NaN);
Na2O  = getAmpOxOptional(data_amp, 'Na2O',  NaN);
K2O   = getAmpOxOptional(data_amp, 'K2O',   0);
Cr2O3 = getAmpOxOptional(data_amp, 'Cr2O3', 0);
NiO   = getAmpOxOptional(data_amp, 'NiO',   0);
F     = getAmpOxOptional(data_amp, 'F',     0);
Cl    = getAmpOxOptional(data_amp, 'Cl',    0);

% Raw cation moles.
catSi  = SiO2  .* MWinfo.Cat.SiO2  ./ MWinfo.MW.SiO2;
catTi  = TiO2  .* MWinfo.Cat.TiO2  ./ MWinfo.MW.TiO2;
catAl  = Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
catFe2 = FeO   .* MWinfo.Cat.FeO   ./ MWinfo.MW.FeO;
catFe3 = Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;
catMn  = MnO   .* MWinfo.Cat.MnO   ./ MWinfo.MW.MnO;
catMg  = MgO   .* MWinfo.Cat.MgO   ./ MWinfo.MW.MgO;
catCa  = CaO   .* MWinfo.Cat.CaO   ./ MWinfo.MW.CaO;
catNa  = Na2O  .* MWinfo.Cat.Na2O  ./ MWinfo.MW.Na2O;
catK   = K2O   .* MWinfo.Cat.K2O   ./ MWinfo.MW.K2O;
catCr  = Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
catNi  = NiO   .* MWinfo.Cat.NiO   ./ MWinfo.MW.NiO;

catFeT = catFe2 + catFe3;
sumNonCNK = catSi + catTi + catAl + catFeT + catMn + catMg + ...
    catCr + catNi;

if ~(isfinite(sumNonCNK) && sumNonCNK > 0)
    amp = emptyAmpStruct();
    return
end

% 13-CNK scaling.
scaleFactor = 13 ./ sumNonCNK;

Si = catSi .* scaleFactor;
Ti = catTi .* scaleFactor;
Al = catAl .* scaleFactor;
FeT = catFeT .* scaleFactor;
Mn = catMn .* scaleFactor;
Mg = catMg .* scaleFactor;
Ca = catCa .* scaleFactor;
Na = catNa .* scaleFactor;
K  = catK  .* scaleFactor;
Cr = catCr .* scaleFactor;
Ni = catNi .* scaleFactor;

% Oxygen sum on the 13-CNK basis using the Fe split supplied by the oxide
% columns. F and Cl are included in the oxygen-equivalent correction.
O_Si  = (SiO2  ./ MWinfo.MW.SiO2)  .* 2 .* scaleFactor;
O_Ti  = (TiO2  ./ MWinfo.MW.TiO2)  .* 2 .* scaleFactor;
O_Al  = (Al2O3 ./ MWinfo.MW.Al2O3) .* 3 .* scaleFactor;
O_Fe2 = (FeO   ./ MWinfo.MW.FeO)   .* 1 .* scaleFactor;
O_Fe3 = (Fe2O3 ./ MWinfo.MW.Fe2O3) .* 3 .* scaleFactor;
O_Mn  = (MnO   ./ MWinfo.MW.MnO)   .* 1 .* scaleFactor;
O_Mg  = (MgO   ./ MWinfo.MW.MgO)   .* 1 .* scaleFactor;
O_Ca  = (CaO   ./ MWinfo.MW.CaO)   .* 1 .* scaleFactor;
O_Na  = (Na2O  ./ MWinfo.MW.Na2O)  .* 1 .* scaleFactor;
O_K   = (K2O   ./ MWinfo.MW.K2O)   .* 1 .* scaleFactor;
O_Cr  = (Cr2O3 ./ MWinfo.MW.Cr2O3) .* 3 .* scaleFactor;
O_Ni  = (NiO   ./ MWinfo.MW.NiO)   .* 1 .* scaleFactor;
O_F   = (F     ./ MWinfo.MW.F)     .* 1 .* scaleFactor;
O_Cl  = (Cl    ./ MWinfo.MW.Cl)    .* 1 .* scaleFactor;

Ototal = O_Si + O_Ti + O_Al + O_Fe2 + O_Fe3 + O_Mn + O_Mg + ...
    O_Ca + O_Na + O_K + O_Cr + O_Ni - 0.5 .* (O_F + O_Cl);

if isfinite(Ototal)
    Fe3_est = max(0, 2 .* (23 - Ototal));
else
    Fe3_est = NaN;
end

if isfinite(Fe3_est) && isfinite(FeT)
    Fe3_est = min(Fe3_est, FeT);
end
Fe2_est = FeT - Fe3_est;

amp = finalizeAmpStruct( ...
    Si, Ti, Al, Fe3_est, Fe2_est, Mn, Mg, Ca, Na, K, Cr, Ni);

end

function amp = prepareAmphibole13CNKFromCations(data_amp)
% prepareAmphibole13CNKFromCations
% Use recognized 23-oxygen/apfu cation columns when oxide data are absent.
% Existing NaN values remain NaN and are not replaced by fallback values.

Si = getVarOptional(data_amp, ...
    {'Si_23O','Si23O','Si_cation_apfu','Si_cation'}, NaN);
Ti = getVarOptional(data_amp, ...
    {'Ti_23O','Ti23O','Ti_cation_apfu','Ti_cation'}, 0);
Al = getVarOptional(data_amp, ...
    {'Al_23O','Al23O','Al_cation_apfu','Al_cation'}, NaN);

[Fe3, hasFe3] = getVarOptional(data_amp, ...
    {'Fe3_23O','Fe3+_23O','Fe3_cation_apfu','Fe3_cation'}, 0);
[Fe2, hasFe2] = getVarOptional(data_amp, ...
    {'Fe2_23O','Fe2+_23O','Fe2_cation_apfu','Fe2_cation'}, NaN);

Mn = getVarOptional(data_amp, ...
    {'Mn_23O','Mn23O','Mn_cation_apfu','Mn_cation'}, 0);
Mg = getVarOptional(data_amp, ...
    {'Mg_23O','Mg23O','Mg_cation_apfu','Mg_cation'}, NaN);
Ca = getVarOptional(data_amp, ...
    {'Ca_23O','Ca23O','Ca_cation_apfu','Ca_cation'}, NaN);
Na = getVarOptional(data_amp, ...
    {'Na_23O','Na23O','Na_cation_apfu','Na_cation'}, NaN);
K = getVarOptional(data_amp, ...
    {'K_23O','K23O','K_cation_apfu','K_cation'}, 0);
Cr = getVarOptional(data_amp, ...
    {'Cr_23O','Cr23O','Cr_cation_apfu','Cr_cation'}, 0);
Ni = getVarOptional(data_amp, ...
    {'Ni_23O','Ni23O','Ni_cation_apfu','Ni_cation'}, 0);

% Use total Fe only when an Fe2+ column is genuinely absent. An explicit NaN
% in an existing Fe2+ column is preserved and is not replaced.
if ~hasFe2
    [FeT, hasFeT] = getVarOptional(data_amp, ...
        {'Fe_23O','Fe23O','Fe_cation_apfu','Fe_cation'}, 0);

    if hasFeT
        if ~hasFe3
            Fe3 = 0;
        end
        Fe2 = FeT - Fe3;
    else
        Fe2 = 0;
    end
end

% Approximate 13-CNK scaling if the required sum is finite and positive.
sumNonCNK = Si + Ti + Al + Fe3 + Fe2 + Mn + Mg + Cr + Ni;

if isfinite(sumNonCNK) && sumNonCNK > 0
    scaleFactor = 13 ./ sumNonCNK;
    Si  = Si  .* scaleFactor;
    Ti  = Ti  .* scaleFactor;
    Al  = Al  .* scaleFactor;
    Fe3 = Fe3 .* scaleFactor;
    Fe2 = Fe2 .* scaleFactor;
    Mn  = Mn  .* scaleFactor;
    Mg  = Mg  .* scaleFactor;
    Ca  = Ca  .* scaleFactor;
    Na  = Na  .* scaleFactor;
    K   = K   .* scaleFactor;
    Cr  = Cr  .* scaleFactor;
    Ni  = Ni  .* scaleFactor;
else
    amp = emptyAmpStruct();
    return
end

amp = finalizeAmpStruct( ...
    Si, Ti, Al, Fe3, Fe2, Mn, Mg, Ca, Na, K, Cr, Ni);

end

function amp = finalizeAmpStruct(Si, Ti, Al, Fe3, Fe2, Mn, Mg, Ca, Na, K, Cr, Ni)
% finalizeAmpStruct
% Allocate tetrahedral Al, M4 Ca/Na, A-site Na, vacancy, and mole fractions
% from normalized amphibole cations.

AlIV = max(8 - Si, 0);
AlVI = max(Al - AlIV, 0);

CaM4 = min(max(Ca, 0), 2);
NaM4 = min(max(2 - CaM4, 0), max(Na, 0));
NaA = max(Na - NaM4, 0);

sumNoVac = Si + Ti + Al + Fe3 + Fe2 + Mn + Mg + Ca + Na + K + Cr + Ni;
Vc = max(16 - sumNoVac, 0);

denCaM4 = CaM4 + NaM4;
if isfinite(denCaM4) && denCaM4 > 0
    CaM4_ratio = CaM4 ./ denCaM4;
else
    CaM4_ratio = NaN;
end

amp = struct( ...
    'Si', Si, 'Ti', Ti, 'Al', Al, 'Fe3', Fe3, 'Fe2', Fe2, ...
    'Mn', Mn, 'Mg', Mg, 'Ca', Ca, 'Na', Na, 'K', K, ...
    'Cr', Cr, 'Ni', Ni, ...
    'AlIV', AlIV, 'AlVI', AlVI, ...
    'CaM4', CaM4, 'NaM4', NaM4, 'NaA', NaA, ...
    'Vc', Vc, 'SumCat', 16, 'XMg', Mg ./ 16, ...
    'CaM4_ratio', CaM4_ratio);

end

function amp = emptyAmpStruct()
% emptyAmpStruct
% Return a stable amphibole structure filled with NaN for failed or missing
% normalization while retaining the fixed total-cation denominator of 16.

amp = struct( ...
    'Si', NaN, 'Ti', NaN, 'Al', NaN, 'Fe3', NaN, 'Fe2', NaN, ...
    'Mn', NaN, 'Mg', NaN, 'Ca', NaN, 'Na', NaN, 'K', NaN, ...
    'Cr', NaN, 'Ni', NaN, 'AlIV', NaN, 'AlVI', NaN, ...
    'CaM4', NaN, 'NaM4', NaN, 'NaA', NaN, ...
    'Vc', NaN, 'SumCat', 16, 'XMg', NaN, 'CaM4_ratio', NaN);

end

function nanInputNames = findNaNInputs(data_amp, data_liq)
% findNaNInputs
% Return recognized thermometer input names that explicitly contain NaN.
% The output buffer is preallocated and trimmed once after scanning.

maxNames = 32;
nameBuffer = strings(maxNames, 1);
nNames = 0;

ampVariableNames = data_amp.Properties.VariableNames;
usesAmpOxides = ~isempty(findOxideColumn(ampVariableNames, 'SiO2'));

if usesAmpOxides
    ampOxides = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO', ...
        'CaO','Na2O','K2O','Cr2O3','NiO','F','Cl'};

    for i = 1:numel(ampOxides)
        columnName = findOxideColumn(ampVariableNames, ampOxides{i});
        if ~isempty(columnName)
            value = toScalarDouble(data_amp.(columnName), NaN);
            if isnan(value)
                nNames = nNames + 1;
                nameBuffer(nNames) = "Amphibole." + string(columnName);
            end
        end
    end
else
    ampGroups = { ...
        {'Si_23O','Si23O','Si_cation_apfu','Si_cation'}, ...
        {'Ti_23O','Ti23O','Ti_cation_apfu','Ti_cation'}, ...
        {'Al_23O','Al23O','Al_cation_apfu','Al_cation'}, ...
        {'Fe3_23O','Fe3+_23O','Fe3_cation_apfu','Fe3_cation'}, ...
        {'Fe2_23O','Fe2+_23O','Fe2_cation_apfu','Fe2_cation'}, ...
        {'Fe_23O','Fe23O','Fe_cation_apfu','Fe_cation'}, ...
        {'Mn_23O','Mn23O','Mn_cation_apfu','Mn_cation'}, ...
        {'Mg_23O','Mg23O','Mg_cation_apfu','Mg_cation'}, ...
        {'Ca_23O','Ca23O','Ca_cation_apfu','Ca_cation'}, ...
        {'Na_23O','Na23O','Na_cation_apfu','Na_cation'}, ...
        {'K_23O','K23O','K_cation_apfu','K_cation'}, ...
        {'Cr_23O','Cr23O','Cr_cation_apfu','Cr_cation'}, ...
        {'Ni_23O','Ni23O','Ni_cation_apfu','Ni_cation'}};

    for i = 1:numel(ampGroups)
        [value, found, matchedName] = getVarOptional( ...
            data_amp, ampGroups{i}, NaN);
        if found && isnan(value)
            nNames = nNames + 1;
            nameBuffer(nNames) = "Amphibole." + string(matchedName);
        end
    end
end

liqVariableNames = data_liq.Properties.VariableNames;
liqOxides = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO', ...
    'CaO','Na2O','K2O','Cr2O3','NiO','P2O5'};

for i = 1:numel(liqOxides)
    columnName = findOxideColumn(liqVariableNames, liqOxides{i});
    if ~isempty(columnName)
        value = toScalarDouble(data_liq.(columnName), NaN);
        if isnan(value)
            nNames = nNames + 1;
            nameBuffer(nNames) = "Liquid." + string(columnName);
        end
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_amp, data_liq)
% validateNonNegativeInputs
% Stop when a recognized finite composition input is negative. Zero and NaN
% are intentionally allowed: zero may later produce NaN through an invalid
% ratio/logarithm, while NaN is retained and reported non-fatally.

maxNames = 32;
nameBuffer = strings(maxNames, 1);
nNames = 0;

ampVariableNames = data_amp.Properties.VariableNames;
usesAmpOxides = ~isempty(findOxideColumn(ampVariableNames, 'SiO2'));

if usesAmpOxides
    ampOxides = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO', ...
        'CaO','Na2O','K2O','Cr2O3','NiO','F','Cl'};

    for i = 1:numel(ampOxides)
        columnName = findOxideColumn(ampVariableNames, ampOxides{i});
        if ~isempty(columnName)
            value = toScalarDouble(data_amp.(columnName), NaN);
            if isfinite(value) && value < 0
                nNames = nNames + 1;
                nameBuffer(nNames) = "Amphibole." + string(columnName);
            end
        end
    end
else
    ampGroups = { ...
        {'Si_23O','Si23O','Si_cation_apfu','Si_cation'}, ...
        {'Ti_23O','Ti23O','Ti_cation_apfu','Ti_cation'}, ...
        {'Al_23O','Al23O','Al_cation_apfu','Al_cation'}, ...
        {'Fe3_23O','Fe3+_23O','Fe3_cation_apfu','Fe3_cation'}, ...
        {'Fe2_23O','Fe2+_23O','Fe2_cation_apfu','Fe2_cation'}, ...
        {'Fe_23O','Fe23O','Fe_cation_apfu','Fe_cation'}, ...
        {'Mn_23O','Mn23O','Mn_cation_apfu','Mn_cation'}, ...
        {'Mg_23O','Mg23O','Mg_cation_apfu','Mg_cation'}, ...
        {'Ca_23O','Ca23O','Ca_cation_apfu','Ca_cation'}, ...
        {'Na_23O','Na23O','Na_cation_apfu','Na_cation'}, ...
        {'K_23O','K23O','K_cation_apfu','K_cation'}, ...
        {'Cr_23O','Cr23O','Cr_cation_apfu','Cr_cation'}, ...
        {'Ni_23O','Ni23O','Ni_cation_apfu','Ni_cation'}};

    for i = 1:numel(ampGroups)
        [value, found, matchedName] = getVarOptional( ...
            data_amp, ampGroups{i}, NaN);
        if found && isfinite(value) && value < 0
            nNames = nNames + 1;
            nameBuffer(nNames) = "Amphibole." + string(matchedName);
        end
    end
end

liqVariableNames = data_liq.Properties.VariableNames;
liqOxides = {'SiO2','TiO2','Al2O3','FeO','Fe2O3','MnO','MgO', ...
    'CaO','Na2O','K2O','Cr2O3','NiO','P2O5'};

for i = 1:numel(liqOxides)
    columnName = findOxideColumn(liqVariableNames, liqOxides{i});
    if ~isempty(columnName)
        value = toScalarDouble(data_liq.(columnName), NaN);
        if isfinite(value) && value < 0
            nNames = nNames + 1;
            nameBuffer(nNames) = "Liquid." + string(columnName);
        end
    end
end

if nNames > 0
    invalidNames = nameBuffer(1:nNames);
    error(['Molina2015AmphLiq: composition values must be >= 0. ' ...
           'Negative finite value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Replicate available liquid identifiers for every pressure output row.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    indexValues = data_liq.('Index');
    row.liq_Index = repmat(indexValues(1), nRows, 1);
end
if any(strcmp(variableNames, 'Experiment'))
    experimentValues = data_liq.('Experiment');
    row.liq_Experiment = repmat(string(experimentValues(1)), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    citationValues = data_liq.('Citation');
    row.liq_Citation = repmat(string(citationValues(1)), nRows, 1);
end

end

function items = buildLiquidList(liq)
% buildLiquidList
% Build a preallocated list of liquid rows for the selection dialog.

nRows = height(liq);
items = cell(nRows, 1);

hasIndex = any(strcmp(liq.Properties.VariableNames, 'Index'));
hasExperiment = any(strcmp(liq.Properties.VariableNames, 'Experiment'));
hasCitation = any(strcmp(liq.Properties.VariableNames, 'Citation'));

indexValues = [];
experimentValues = [];
citationValues = [];

if hasIndex
    indexValues = liq.('Index');
end
if hasExperiment
    experimentValues = liq.('Experiment');
end
if hasCitation
    citationValues = liq.('Citation');
end

for i = 1:nRows
    partBuffer = strings(4, 1);
    nParts = 1;
    partBuffer(nParts) = "Row " + string(i);

    if hasIndex
        nParts = nParts + 1;
        partBuffer(nParts) = "Index=" + string(indexValues(i));
    end
    if hasExperiment
        nParts = nParts + 1;
        partBuffer(nParts) = string(experimentValues(i));
    end
    if hasCitation
        nParts = nParts + 1;
        partBuffer(nParts) = string(citationValues(i));
    end

    items{i} = char(strjoin(partBuffer(1:nParts), ' | '));
end

end

function label = makeLiquidLabel(data_liq, rowNumber)
% makeLiquidLabel
% Create a compact text identifier for console output and warning messages.

label = ['Liquid row ' num2str(rowNumber)];
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Experiment'))
    experimentValues = data_liq.('Experiment');
    experimentText = char(string(experimentValues(1)));
    if ~isempty(experimentText)
        label = experimentText;
        return
    end
end

if any(strcmp(variableNames, 'Index'))
    indexValues = data_liq.('Index');
    indexText = char(string(indexValues(1)));
    if ~isempty(indexText)
        label = ['Liquid Index ' indexText];
    end
end

end

function value = getAmpOxOptional(data_amp, oxide, defaultValue)
% getAmpOxOptional
% Return an amphibole oxide value. Missing columns use defaultValue; explicit
% NaN in an existing column is preserved.

columnName = findOxideColumn(data_amp.Properties.VariableNames, oxide);
if isempty(columnName)
    value = defaultValue;
else
    value = toScalarDouble(data_amp.(columnName), defaultValue);
end

end

function value = getLiqOxOptional(data_liq, oxide, defaultValue)
% getLiqOxOptional
% Return a liquid oxide value. Missing columns use defaultValue; explicit NaN
% in an existing column is preserved.

columnName = findOxideColumn(data_liq.Properties.VariableNames, oxide);
if isempty(columnName)
    value = defaultValue;
else
    value = toScalarDouble(data_liq.(columnName), defaultValue);
end

end

function [value, found, matchedName] = getVarOptional(dataRow, candidateNames, defaultValue)
% getVarOptional
% Return the first matching variable from candidateNames. The found flag
% distinguishes an absent column from an existing column that contains NaN.

value = defaultValue;
found = false;
matchedName = '';
variableNames = dataRow.Properties.VariableNames;

for i = 1:numel(candidateNames)
    if any(strcmp(variableNames, candidateNames{i}))
        matchedName = candidateNames{i};
        value = toScalarDouble(dataRow.(matchedName), defaultValue);
        found = true;
        return
    end
end

end

function name = findOxideColumn(varNames, oxide)
% findOxideColumn
% Match common oxide-column variants after removing spaces, underscores, and
% hyphens and converting names to lower case.

canonicalNames = cell(size(varNames));

for i = 1:numel(varNames)
    textValue = lower(varNames{i});
    textValue = strrep(textValue, ' ', '');
    textValue = strrep(textValue, '_', '');
    textValue = strrep(textValue, '-', '');
    canonicalNames{i} = textValue;
end

canonicalOxide = lower(oxide);
canonicalOxide = strrep(canonicalOxide, ' ', '');
canonicalOxide = strrep(canonicalOxide, '_', '');
canonicalOxide = strrep(canonicalOxide, '-', '');

targets = {[canonicalOxide 'value'], canonicalOxide};
name = '';

for i = 1:numel(targets)
    idx = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(idx)
        name = varNames{idx};
        return
    end
end

end

function value = toScalarDouble(raw, defaultValue)
% toScalarDouble
% Convert the first table value to double. Explicit numeric NaN and the text
% "NaN" are retained as NaN rather than replaced by defaultValue.

value = defaultValue;

if isempty(raw)
    return
end

if isnumeric(raw)
    value = double(raw(1));
    return
end

if islogical(raw)
    value = double(raw(1));
    return
end

if isstring(raw)
    if ismissing(raw(1))
        value = NaN;
        return
    end

    textValue = strtrim(raw(1));
    if strcmpi(textValue, "NaN")
        value = NaN;
        return
    end

    convertedValue = str2double(textValue);
    if ~isnan(convertedValue)
        value = convertedValue;
    end
    return
end

if ischar(raw)
    textValue = strtrim(string(raw));
    if strcmpi(textValue, "NaN")
        value = NaN;
        return
    end

    convertedValue = str2double(textValue);
    if ~isnan(convertedValue)
        value = convertedValue;
    end
    return
end

if iscell(raw)
    firstValue = raw{1};

    if isempty(firstValue)
        return
    end

    if isnumeric(firstValue) || islogical(firstValue)
        value = double(firstValue(1));
        return
    end

    if isstring(firstValue) || ischar(firstValue)
        textValue = strtrim(string(firstValue));
        if strcmpi(textValue, "NaN")
            value = NaN;
            return
        end

        convertedValue = str2double(textValue);
        if ~isnan(convertedValue)
            value = convertedValue;
        end
    end
end

end
