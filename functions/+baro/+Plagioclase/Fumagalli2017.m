function results = Fumagalli2017(rawdata_struct, T_degreeC)
% functions/+baro/+Plagioclase/Fumagalli2017.m
% Tested with MATLAB R2024b
%
% Forsterite-Anorthite-Ca-Tschermak-Enstatite (FACE) geobarometer
% Fumagalli, P., Borghini, G., Rampone, E. and Poli, S. (2017)
% Contributions to Mineralogy and Petrology, 172, Article 38
% DOI: https://doi.org/10.1007/s00410-017-1352-2
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Olivine analysis, one Plagioclase
% analysis, one Clinopyroxene analysis and one Orthopyroxene analysis and
% calculates pressure using the FACE geobarometer of Fumagalli et al.
% (2017).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Ol-Pl-Cpx-Opx assemblage, one
% output row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Fumagalli et al. (2017) calibrated the pressure-sensitive equilibrium:
%
%   Forsterite + Anorthite = Ca-Tschermak + Enstatite
%
%   Mg2SiO4(Ol) + CaAl2Si2O8(Pl)
%       = CaAl2SiO6(Cpx) + Mg2Si2O6(Opx)
%
% The new experiments of Fumagalli et al. (2017) were performed at:
%
%   Temperature : 1050-1150 degreeC
%   Pressure    : 5-10 kbar
%   Conditions  : nominally anhydrous, subsolidus
%
% These run conditions are reported on p. 4 and in Table 2. The combined
% calibration database, including Borghini et al. (2010, 2011), is reported
% in the abstract as 3-9 kbar and 1000-1150 degreeC (p. 1). Table 4 on
% p. 13 contains calibration experiments from 2.7 to 9 kbar. This
% implementation uses 3-9 kbar and 1000-1150 degreeC as the strict
% non-stopping warning range.
%
% The experimental bulk compositions cover approximately:
%
%   Na2O/CaO              : 0.08-0.13
%   XCr = Cr/(Cr + Al)    : 0.07-0.10
%   XMg                   : approximately 0.89-0.90
%
% These compositional ranges are summarized on pp. 1 and 3 and in Table 1.
% The selected mineral analyses alone do not provide whole-rock Na2O/CaO or
% XCr; therefore, this function cannot automatically test those bulk-rock
% calibration limits.
%
% IMPORTANT APPLICATION REQUIREMENT:
% A correct application requires accurate coupling of mineral composition
% and textural occurrence. Olivine, plagioclase, clinopyroxene and
% orthopyroxene must belong to the same equilibrated, texturally associated
% plagioclase-bearing assemblage. Core and rim analyses, different mineral
% generations, or unrelated microstructural domains must not be mixed.
% Fumagalli et al. (2017) emphasize texturally associated neoblasts and
% detailed microstructural control on pp. 15-16. They explicitly state that
% the lack of texturally controlled analyses prevents application to some
% Horoman and Uenzaru samples (p. 16).
%
% Experiments below 1000 degreeC were prevented by slow reaction rates under
% nominally anhydrous conditions (p. 4). Natural applications at lower
% temperature, including the examples on pp. 15-16, are therefore
% extrapolations beyond the experimental temperature calibration.
%
% The plagioclase-out boundary is composition dependent. At near-solidus
% conditions it occurs at approximately 7-10 kbar in variably depleted
% lherzolites and is controlled mainly by bulk Na2O/CaO and XCr (pp. 10 and
% 16). A pressure above the plagioclase stability field should trigger a
% re-examination of phase equilibrium, mineral pairing and bulk composition.
%
% The standard error of the regression is 0.42 kbar, and propagation of
% analytical and experimental uncertainties by Monte Carlo simulation gives
% an uncertainty of approximately +/-0.5 kbar (p. 14). When uncertainty in
% the fitted coefficients is also considered, the practical accuracy is
% discussed as approximately 2.5-3 kbar (p. 15).
%
% All Fe is treated as Fe2+ in the activity models, following the reducing
% experimental conditions. The specified activity-composition relations for
% plagioclase, clinopyroxene, orthopyroxene and olivine are described on
% pp. 12-13.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 1000-1150 degreeC,
%   2) finite calculated pressure is outside 3-9 kbar,
%   3) a calculation input contains NaN,
%   4) a required activity or equilibrium constant is non-positive or
%      non-finite, or
%   5) calculated pressure is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Olivine     : table
%   rawdata_struct.Plagioclase : table
%   rawdata_struct.Cpx         : table
%   rawdata_struct.Opx         : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following normalized-cation variables. The standard *_cation_apfu names
% are preferred; compatible aliases from the original implementation are
% also accepted.
%
%   Olivine calculation variables:
%     Mg_cation_apfu
%     Fe_cation_apfu or Fe2_cation_apfu
%
%   Plagioclase calculation variables:
%     Ca_cation_apfu
%     Na_cation_apfu
%
%   Clinopyroxene calculation variables:
%     Si_cation_apfu
%     Al_cation_apfu
%     Ca_cation_apfu
%     Na_cation_apfu
%     Ti_cation_apfu
%     Cr_cation_apfu
%
%   Orthopyroxene calculation variables:
%     Al_cation_apfu
%     Mg_cation_apfu
%     Fe_cation_apfu or Fe2_cation_apfu
%     Na_cation_apfu
%     Ti_cation_apfu
%     Cr_cation_apfu
%     Ca_cation_apfu
%
% Mn is not used in the Fumagalli et al. (2017) activity equations. If Mn
% columns are present, their selected values are retained in the output only.
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if it makes an activity, ratio or
% logarithm undefined, the resulting NaN or Inf is retained and reported.
%
% No liquid composition is used by this barometer. Therefore, exclusion of
% Liq F and Cl from cationTotal_liq and from Liq NaN warnings is not
% applicable to this function.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
%   P(kbar) = 7.2 + 0.0078*T(K) + 0.0022*T(K)*ln(K)
%
%   K = (aCaTs_cpx * aEn_opx) / (aAn_plag * aFo_ol)
%
% The pressure equation is Eq. (3) on p. 14.
%
% Olivine (p. 13):
%   XMg_ol = Mg / (Mg + Fe)
%   aFo_ol = XMg_ol^2
%
% Plagioclase (p. 12; Holland and Powell, 1992):
%   XAn = Ca / (Ca + Na)
%   Xb = 0.12 + 0.00038*T(K)
%   XAn_c = 0.25*XAn*(1 + XAn)^2
%   IAn = -R*T*ln(XAn_c/XAn) - (Wc - Wi)*(1 - Xb)^2
%   aAn = XAn_c*exp([Wc*(1 - XAn)^2 + IAn]/[R*T])
%   Wc = 1070 J/mol; Wi = 9790 J/mol
%
% Clinopyroxene (pp. 12-13):
%   AlIV = (Al + Cr + 2Ti - Na)/2
%   AlVI = Al - AlIV
%   aCaTs = 4*XCa_M2*XAl_M1*XAl_T*XSi_T
%
% Orthopyroxene (p. 13):
%   AlVI = (Al + Cr + 2Ti - Na)/2
%   aEn = XMg_M1*XMg_M2
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Fumagalli2017(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing Olivine, Plagioclase, Cpx and Opx
%                    tables
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Ol-Pl-Cpx-Opx assemblage.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('Fumagalli2017 requires (rawdata_struct, T_degreeC).');
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

%% 1) Retrieve cation datasets
% Extract required tables from the input struct. The source tables are not
% modified; selected rows are read during each calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end
if ~isfield(rawdata_struct, 'Plagioclase') || ~istable(rawdata_struct.Plagioclase)
    error('rawdata_struct must contain table: rawdata_struct.Plagioclase');
end
if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end
if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end

dataset_ol = rawdata_struct.Olivine;
dataset_plag = rawdata_struct.Plagioclase;
dataset_cpx = rawdata_struct.Cpx;
dataset_opx = rawdata_struct.Opx;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-assemblage result in a cell buffer and concatenate once
% after the interactive loop. This avoids repeated growth of the output
% table on every loop iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degreeC = 1000;
calibrationT_max_degreeC = 1150;
calibrationP_min_kbar = 3;
calibrationP_max_kbar = 9;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);
temperatureWarningIssued = false;
applicationCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-7) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
while true
    % ----- Olivine selection -----
    disp('=== Step 3: Selecting a data code from the list (Olivine) ===');
    dataCodes_ol = dataset_ol{:, 1};

    [selectedIdx_ol, ok] = listdlg( ...
        'PromptString', 'Please select the Olivine data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_ol)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_ol)
        disp('Selection canceled');
        break;
    end

    selectedCode_ol = dataCodes_ol(selectedIdx_ol);
    disp(['Olivine selected: ' char(string(selectedCode_ol))]);

    % ----- Plagioclase selection -----
    disp('=== Step 4: Selecting a data code from the list (Plagioclase) ===');
    dataCodes_plag = dataset_plag{:, 1};

    [selectedIdx_plag, ok] = listdlg( ...
        'PromptString', 'Please select the Plagioclase data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_plag)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_plag)
        disp('Selection canceled');
        break;
    end

    selectedCode_plag = dataCodes_plag(selectedIdx_plag);
    disp(['Plagioclase selected: ' char(string(selectedCode_plag))]);

    % ----- Clinopyroxene selection -----
    disp('=== Step 5: Selecting a data code from the list (Clinopyroxene) ===');
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

    % ----- Orthopyroxene selection -----
    disp('=== Step 6: Selecting a data code from the list (Orthopyroxene) ===');
    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Orthopyroxene data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Orthopyroxene selected: ' char(string(selectedCode_opx))]);

    % ----- Calculation -----
    disp('=== Step 7: Calculating the pressure ===');

    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_plag = dataset_plag(selectedIdx_plag, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);
    selectedData_opx = dataset_opx(selectedIdx_opx, :);

    ol = prepareOlivineRow(selectedData_ol);
    plag = preparePlagioclaseRow(selectedData_plag);
    cpx = prepareClinopyroxeneRow(selectedData_cpx);
    opx = prepareOrthopyroxeneRow(selectedData_opx);

    % NaN values are identified and reported, but they are not replaced.
    nanInputNames = findNaNInputs(ol, plag, cpx, opx, T_degreeC);

    % Reject Inf and finite negative values only. Zero and NaN are retained.
    validateNonNegativeInputs(ol, plag, cpx, opx);

    row = calcPressure(ol, plag, cpx, opx, T_degreeC);

    % Repeat identifiers for all temperatures in the current calculation.
    row.dataCode_ol = repmat(string(selectedCode_ol), height(row), 1);
    row.dataCode_pl = repmat(string(selectedCode_plag), height(row), 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row = movevars(row, ...
        {'dataCode_ol','dataCode_pl','dataCode_cpx','dataCode_opx'}, ...
        'Before', 1);

    % Store one block per selected assemblage. The buffer is expanded by
    % doubling only when its capacity is exhausted, not on every iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure result for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    assemblageLabel = [char(string(selectedCode_ol)) ' & ' ...
        char(string(selectedCode_plag)) ' & ' ...
        char(string(selectedCode_cpx)) ' & ' ...
        char(string(selectedCode_opx))];

    if height(row) == 1
        disp([assemblageLabel ': P = ' num2str(row.P_kbar) ' kbar']);
    else
        disp([assemblageLabel ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Print the textural-equilibrium limitation once per function call.
    if ~applicationCautionIssued
        fprintf(2, ...
            ['CAUTION: Correct use of the Fumagalli et al. (2017) FACE ' ...
             'geobarometer requires Ol-Pl-Cpx-Opx analyses from the same ' ...
             'texturally associated, equilibrated plagioclase-bearing ' ...
             'assemblage. Core-rim pairs and different mineral generations ' ...
             'must not be mixed (pp. 15-16).\n']);
        applicationCautionIssued = true;
    end

    % Input temperature is common to all selected assemblages, so this
    % warning is printed only once after the first completed calculation.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the Fumagalli et al. ' ...
             '(2017) calibration range of 1000-1150 degreeC (p. 1; Table 4, ' ...
             'p. 13). %d of %d finite temperature point(s) are outside the ' ...
             'range; input finite range = %.4g-%.4g degreeC. Values are ' ...
             'retained and calculated as extrapolations.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures fall outside the strict
    % calibration range used by this implementation.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the Fumagalli et al. ' ...
             '(2017) strict calibration range of 3-9 kbar (p. 1; Table 4, ' ...
             'p. 13). %d of %d finite pressure point(s) are outside the ' ...
             'range; calculated finite range = %.4g-%.4g kbar for %s. ' ...
             'Values are retained for diagnostic purposes.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            assemblageLabel);
    end

    % List the exact calculation inputs containing NaN. For vector
    % temperature input, the NaN element indices are included.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the FACE geobarometer input(s) for ' ...
             '%s: %s.\n' ...
             '         NaN values were retained and propagated; calculated ' ...
             'pressure may remain NaN.\n'], ...
            assemblageLabel, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Report derived quantities outside their mathematical/physical domain.
    reportDerivedDomainWarnings(row, assemblageLabel);

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, +Inf: %d, -Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            assemblageLabel, ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar) & row.P_kbar > 0), ...
            sum(isinf(row.P_kbar) & row.P_kbar < 0));
    end

    % Negative calculated pressure is retained but is outside the physical
    % and calibration domain.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s ' ...
             '(%d of %d points). Values were retained for diagnostic ' ...
             'purposes.\n'], ...
            assemblageLabel, ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Fumagalli2017', ...
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
function nanInputNames = findNaNInputs(ol, plag, cpx, opx, T_degreeC)
% findNaNInputs
% Return names of calculation inputs containing NaN. NaN values do not cause
% an error and are not replaced by zero.

[inputNames, inputValues] = buildCalculationInputList(ol, plag, cpx, opx);

maxNames = numel(inputNames) + 1;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

for i = 1:numel(inputValues)
    if isnan(inputValues(i))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = inputNames(i);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(ol, plag, cpx, opx)
% validateNonNegativeInputs
% Reject Inf and finite negative values in variables directly used by the
% pressure equation. Zero and NaN are intentionally allowed and retained.

[inputNames, inputValues] = buildCalculationInputList(ol, plag, cpx, opx);
invalidMask = isinf(inputValues) | ...
    (isfinite(inputValues) & inputValues < 0);

if any(invalidMask)
    invalidInputNames = inputNames(invalidMask);
    error(['Fumagalli2017: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are ' ...
           'prohibited. Invalid input(s): ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function [inputNames, inputValues] = buildCalculationInputList(ol, plag, cpx, opx)
% buildCalculationInputList
% Assemble names and scalar values used directly or indirectly by the FACE
% pressure calculation.

inputNames = [ ...
    "Olivine." + ol.MgName; ...
    "Olivine." + ol.FeName; ...
    "Plagioclase." + plag.CaName; ...
    "Plagioclase." + plag.NaName; ...
    "Cpx." + cpx.SiName; ...
    "Cpx." + cpx.AlName; ...
    "Cpx." + cpx.CaName; ...
    "Cpx." + cpx.NaName; ...
    "Cpx." + cpx.TiName; ...
    "Cpx." + cpx.CrName; ...
    "Opx." + opx.AlName; ...
    "Opx." + opx.MgName; ...
    "Opx." + opx.FeName; ...
    "Opx." + opx.NaName; ...
    "Opx." + opx.TiName; ...
    "Opx." + opx.CrName; ...
    "Opx." + opx.CaName];

inputValues = [ ...
    ol.Mg; ol.Fe; ...
    plag.Ca; plag.Na; ...
    cpx.Si; cpx.Al; cpx.Ca; cpx.Na; cpx.Ti; cpx.Cr; ...
    opx.Al; opx.Mg; opx.Fe; opx.Na; opx.Ti; opx.Cr; opx.Ca];

end

function row = calcPressure(ol, plag, cpx, opx, T_degreeC)
% calcPressure
% Compute FACE pressure for one selected Ol-Pl-Cpx-Opx assemblage at one or
% more input temperatures. NaN values propagate through the calculation.
%
% Inputs:
%   ol, plag, cpx, opx : prepared one-analysis composition structs
%   T_degreeC          : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

R_J = 8.31446261815324;
Wc_J = 1070.0;
Wi_J = 9790.0;

% -------------------------------------------------------------------------
% Olivine: forsterite activity
% -------------------------------------------------------------------------
% Fumagalli et al. (2017) define XMg using Mg and Fe only. Mn, when
% present, is retained in the output but is not included in this equation.
XMg_ol_scalar = ol.Mg ./ (ol.Mg + ol.Fe);
aFo_ol_scalar = XMg_ol_scalar.^2;

% -------------------------------------------------------------------------
% Plagioclase: anorthite activity
% -------------------------------------------------------------------------
XAn_pl_scalar = plag.Ca ./ (plag.Ca + plag.Na);
Xb_pl = 0.12 + 0.00038 .* T_K;
XAn_c_pl = 0.25 .* XAn_pl_scalar .* (1 + XAn_pl_scalar).^2 .* ones(nT, 1);

IAn_pl = -R_J .* T_K .* log(XAn_c_pl ./ XAn_pl_scalar) ...
    - (Wc_J - Wi_J) .* (1 - Xb_pl).^2;

aAn_pl = XAn_c_pl .* exp((Wc_J .* (1 - XAn_pl_scalar).^2 + IAn_pl) ...
    ./ (R_J .* T_K));

% -------------------------------------------------------------------------
% Clinopyroxene: Ca-Tschermak activity
% -------------------------------------------------------------------------
AlIV_cpx_scalar = (cpx.Al + cpx.Cr + 2 .* cpx.Ti - cpx.Na) ./ 2;
AlVI_cpx_scalar = cpx.Al - AlIV_cpx_scalar;

XAl_M1_cpx_scalar = AlVI_cpx_scalar;
XAl_T_cpx_scalar = AlIV_cpx_scalar ./ (cpx.Si + AlIV_cpx_scalar);
XSi_T_cpx_scalar = cpx.Si ./ (cpx.Si + AlIV_cpx_scalar);
XCa_M2_cpx_scalar = cpx.Ca;

aCaTs_cpx_scalar = 4 .* XCa_M2_cpx_scalar .* XAl_M1_cpx_scalar .* ...
    XAl_T_cpx_scalar .* XSi_T_cpx_scalar;

% -------------------------------------------------------------------------
% Orthopyroxene: enstatite activity
% -------------------------------------------------------------------------
[opxActivity, opxSites] = calcOpxEnstatiteActivity(opx);
aEn_opx_scalar = opxActivity;

% -------------------------------------------------------------------------
% Equilibrium constant and pressure
% -------------------------------------------------------------------------
aFo_ol = repmat(aFo_ol_scalar, nT, 1);
aCaTs_cpx = repmat(aCaTs_cpx_scalar, nT, 1);
aEn_opx = repmat(aEn_opx_scalar, nT, 1);

K_eq = (aCaTs_cpx .* aEn_opx) ./ (aFo_ol .* aAn_pl);

% Avoid complex logarithms. Positive Inf is retained as Inf; zero,
% negative and NaN K values are represented by NaN in lnK.
lnK = NaN(nT, 1);
validLogK = K_eq > 0 & ~isnan(K_eq);
lnK(validLogK) = log(K_eq(validLogK));

P_kbar = 7.2 + 0.0078 .* T_K + 0.0022 .* T_K .* lnK;

% -------------------------------------------------------------------------
% Expand selected composition data to the temperature-vector length
% -------------------------------------------------------------------------
Mg_ol = repmat(ol.Mg, nT, 1);
Fe_ol = repmat(ol.Fe, nT, 1);
Mn_ol = repmat(ol.Mn, nT, 1);
XMg_ol = repmat(XMg_ol_scalar, nT, 1);

Ca_pl = repmat(plag.Ca, nT, 1);
Na_pl = repmat(plag.Na, nT, 1);
XAn_pl = repmat(XAn_pl_scalar, nT, 1);

Si_cpx = repmat(cpx.Si, nT, 1);
Al_cpx = repmat(cpx.Al, nT, 1);
Ca_cpx = repmat(cpx.Ca, nT, 1);
Na_cpx = repmat(cpx.Na, nT, 1);
Ti_cpx = repmat(cpx.Ti, nT, 1);
Cr_cpx = repmat(cpx.Cr, nT, 1);
AlIV_cpx = repmat(AlIV_cpx_scalar, nT, 1);
AlVI_cpx = repmat(AlVI_cpx_scalar, nT, 1);
XAl_M1_cpx = repmat(XAl_M1_cpx_scalar, nT, 1);
XAl_T_cpx = repmat(XAl_T_cpx_scalar, nT, 1);
XSi_T_cpx = repmat(XSi_T_cpx_scalar, nT, 1);
XCa_M2_cpx = repmat(XCa_M2_cpx_scalar, nT, 1);

Al_opx = repmat(opx.Al, nT, 1);
Mg_opx = repmat(opx.Mg, nT, 1);
Fe_opx = repmat(opx.Fe, nT, 1);
Na_opx = repmat(opx.Na, nT, 1);
Ti_opx = repmat(opx.Ti, nT, 1);
Cr_opx = repmat(opx.Cr, nT, 1);
Ca_opx = repmat(opx.Ca, nT, 1);
Mn_opx = repmat(opx.Mn, nT, 1);
AlVI_opx = repmat(opxSites.AlVI, nT, 1);
M1_den_opx = repmat(opxSites.M1_den, nT, 1);
M2_den_opx = repmat(opxSites.M2_den, nT, 1);
XMg_M1_opx = repmat(opxSites.XMg_M1, nT, 1);
XMg_M2_opx = repmat(opxSites.XMg_M2, nT, 1);

R_J_output = repmat(R_J, nT, 1);
Wc_J_output = repmat(Wc_J, nT, 1);
Wi_J_output = repmat(Wi_J, nT, 1);

% -------------------------------------------------------------------------
% Applicability and diagnostic flags
% -------------------------------------------------------------------------
isWithinCalibrationTRange = isfinite(T_degreeC) & ...
    T_degreeC >= 1000 & T_degreeC <= 1150;

isWithinCalibrationPRange = isfinite(P_kbar) & ...
    P_kbar >= 3 & P_kbar <= 9;

isWithinBroadExperimentalEnvelope = isfinite(T_degreeC) & ...
    T_degreeC >= 1000 & T_degreeC <= 1150 & ...
    isfinite(P_kbar) & P_kbar >= 2.5 & P_kbar <= 10;

isWithinEquationDomain = ...
    isfinite(aFo_ol) & aFo_ol > 0 & ...
    isfinite(aAn_pl) & aAn_pl > 0 & ...
    isfinite(aCaTs_cpx) & aCaTs_cpx > 0 & ...
    isfinite(aEn_opx) & aEn_opx > 0 & ...
    isfinite(K_eq) & K_eq > 0;

isRecommended = isWithinCalibrationTRange & ...
    isWithinCalibrationPRange & isWithinEquationDomain;

% -------------------------------------------------------------------------
% Pack outputs using pre-sized vectors of equal height
% -------------------------------------------------------------------------
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.R_J = R_J_output;
row.Wc_J = Wc_J_output;
row.Wi_J = Wi_J_output;

row.Mg_ol = Mg_ol;
row.Fe_ol = Fe_ol;
row.Mn_ol = Mn_ol;
row.XMg_ol = XMg_ol;
row.aFo_ol = aFo_ol;

row.Ca_pl = Ca_pl;
row.Na_pl = Na_pl;
row.XAn_pl = XAn_pl;
row.Xb_pl = Xb_pl;
row.XAn_c_pl = XAn_c_pl;
row.IAn_pl = IAn_pl;
row.aAn_pl = aAn_pl;

row.Si_cpx = Si_cpx;
row.Al_cpx = Al_cpx;
row.Ca_cpx = Ca_cpx;
row.Na_cpx = Na_cpx;
row.Ti_cpx = Ti_cpx;
row.Cr_cpx = Cr_cpx;
row.AlIV_cpx = AlIV_cpx;
row.AlVI_cpx = AlVI_cpx;
row.XAl_M1_cpx = XAl_M1_cpx;
row.XAl_T_cpx = XAl_T_cpx;
row.XSi_T_cpx = XSi_T_cpx;
row.XCa_M2_cpx = XCa_M2_cpx;
row.aCaTs_cpx = aCaTs_cpx;

row.Al_opx = Al_opx;
row.Mg_opx = Mg_opx;
row.Fe_opx = Fe_opx;
row.Na_opx = Na_opx;
row.Ti_opx = Ti_opx;
row.Cr_opx = Cr_opx;
row.Ca_opx = Ca_opx;
row.Mn_opx = Mn_opx;
row.AlVI_opx = AlVI_opx;
row.M1_den_opx = M1_den_opx;
row.M2_den_opx = M2_den_opx;
row.XMg_M1_opx = XMg_M1_opx;
row.XMg_M2_opx = XMg_M2_opx;
row.aEn_opx = aEn_opx;

row.K_eq = K_eq;
row.lnK = lnK;
row.P_kbar = P_kbar;

row.P_standard_error_kbar = repmat(0.42, nT, 1);
row.P_MC_uncertainty_approx_kbar = repmat(0.5, nT, 1);
row.P_accuracy_with_fit_coefficients_kbar = repmat(3.0, nT, 1);

row.isWithinEquationDomain = isWithinEquationDomain;
row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinCalibrationPRange = isWithinCalibrationPRange;
row.isWithinBroadExperimentalEnvelope = isWithinBroadExperimentalEnvelope;
row.isRecommended = isRecommended;

% Backward-compatible aliases retained from the original implementation.
row.is_calibration_T_range = isWithinCalibrationTRange;
row.is_calibration_P_range = isWithinCalibrationPRange;
row.is_recommended = isRecommended;

end

function [aEn, site] = calcOpxEnstatiteActivity(opx)
% calcOpxEnstatiteActivity
% Calculate orthopyroxene enstatite activity exactly from the ideal
% two-site expressions reported by Fumagalli et al. (2017, p. 13):
%
%   XMg_M1 = Mg / (Fe + AlVI + Ti + Cr + Mg)
%   XMg_M2 = Mg / (Ca + Mg + Fe + Na)
%   aEn = XMg_M1 * XMg_M2
%
% NaN inputs remain NaN. Non-positive denominators or non-positive activity
% values are retained as NaN for subsequent non-stopping warnings.

site = struct();
site.AlVI = (opx.Al + opx.Cr + 2 .* opx.Ti - opx.Na) ./ 2;
site.M1_den = opx.Fe + site.AlVI + opx.Ti + opx.Cr + opx.Mg;
site.M2_den = opx.Ca + opx.Mg + opx.Fe + opx.Na;

site.XMg_M1 = opx.Mg ./ site.M1_den;
site.XMg_M2 = opx.Mg ./ site.M2_den;
aEn = site.XMg_M1 .* site.XMg_M2;

invalidDerived = ...
    (isfinite(site.AlVI) && site.AlVI < 0) || ...
    (isfinite(site.M1_den) && site.M1_den <= 0) || ...
    (isfinite(site.M2_den) && site.M2_den <= 0) || ...
    (isfinite(aEn) && aEn <= 0);

if invalidDerived
    aEn = NaN;
end

end

function reportDerivedDomainWarnings(row, assemblageLabel)
% reportDerivedDomainWarnings
% Print non-stopping warnings for activity-model results outside the domain
% required to calculate a real finite pressure.

if any(isfinite(row.XAn_pl) & ...
        (row.XAn_pl <= 0 | row.XAn_pl > 1))
    fprintf(2, ...
        ['WARNING: Plagioclase XAn is outside the mathematical interval ' ...
         '0 < XAn <= 1 for %s. The affected activity and pressure values ' ...
         'were retained as NaN or Inf.\n'], assemblageLabel);
end

activityNames = ["aFo_ol", "aAn_pl", "aCaTs_cpx", "aEn_opx"];
activityValues = {row.aFo_ol, row.aAn_pl, row.aCaTs_cpx, row.aEn_opx};

invalidActivityBuffer = strings(numel(activityNames), 1);
nInvalidActivities = 0;

for i = 1:numel(activityNames)
    value = activityValues{i};
    if any(~isfinite(value) | value <= 0)
        nInvalidActivities = nInvalidActivities + 1;
        invalidActivityBuffer(nInvalidActivities) = activityNames(i);
    end
end

if nInvalidActivities > 0
    invalidActivityNames = invalidActivityBuffer(1:nInvalidActivities);
    fprintf(2, ...
        ['WARNING: Non-positive or non-finite activity value(s) were ' ...
         'calculated for %s: %s. No value was replaced by zero; NaN/Inf ' ...
         'was propagated to K, lnK and pressure.\n'], ...
        assemblageLabel, ...
        char(strjoin(invalidActivityNames, ', ')));
end

if any(~isnan(row.K_eq) & row.K_eq <= 0)
    fprintf(2, ...
        ['WARNING: Non-positive equilibrium constant K was calculated for ' ...
         '%s. lnK and pressure were retained as NaN.\n'], assemblageLabel);
end

end

function ol = prepareOlivineRow(data_olivine)
% prepareOlivineRow
% Extract one-row olivine cation data without replacing NaN by zero.

if height(data_olivine) ~= 1
    error('Olivine input must be a 1-row table.');
end

ol = struct();
[ol.Mg, ol.MgName] = getVarOrError(data_olivine, ...
    {'Mg_cation_apfu','Mg_cation','Mg'}, 'Olivine');
[ol.Fe, ol.FeName] = getVarOrError(data_olivine, ...
    {'Fe_cation_apfu','Fe2_cation_apfu','Fe_cation','Fe2_cation','Fe','Fe2'}, ...
    'Olivine');
ol.Mn = getVarOrNaN(data_olivine, ...
    {'Mn_cation_apfu','Mn_cation','Mn'});

end

function plag = preparePlagioclaseRow(data_plagioclase)
% preparePlagioclaseRow
% Extract one-row plagioclase cation data without replacing NaN by zero.

if height(data_plagioclase) ~= 1
    error('Plagioclase input must be a 1-row table.');
end

plag = struct();
[plag.Ca, plag.CaName] = getVarOrError(data_plagioclase, ...
    {'Ca_cation_apfu','Ca_cation','Ca'}, 'Plagioclase');
[plag.Na, plag.NaName] = getVarOrError(data_plagioclase, ...
    {'Na_cation_apfu','Na_cation','Na'}, 'Plagioclase');

end

function cpx = prepareClinopyroxeneRow(data_cpx)
% prepareClinopyroxeneRow
% Extract one-row clinopyroxene cation data without replacing NaN by zero.

if height(data_cpx) ~= 1
    error('Cpx input must be a 1-row table.');
end

cpx = struct();
[cpx.Si, cpx.SiName] = getVarOrError(data_cpx, ...
    {'Si_cation_apfu','Si_cation','Si'}, 'Cpx');
[cpx.Al, cpx.AlName] = getVarOrError(data_cpx, ...
    {'Al_cation_apfu','Al_cation','Al'}, 'Cpx');
[cpx.Ca, cpx.CaName] = getVarOrError(data_cpx, ...
    {'Ca_cation_apfu','Ca_cation','Ca'}, 'Cpx');
[cpx.Na, cpx.NaName] = getVarOrError(data_cpx, ...
    {'Na_cation_apfu','Na_cation','Na'}, 'Cpx');
[cpx.Ti, cpx.TiName] = getVarOrError(data_cpx, ...
    {'Ti_cation_apfu','Ti_cation','Ti'}, 'Cpx');
[cpx.Cr, cpx.CrName] = getVarOrError(data_cpx, ...
    {'Cr_cation_apfu','Cr_cation','Cr'}, 'Cpx');

end

function opx = prepareOrthopyroxeneRow(data_opx)
% prepareOrthopyroxeneRow
% Extract one-row orthopyroxene cation data without replacing NaN by zero.

if height(data_opx) ~= 1
    error('Opx input must be a 1-row table.');
end

opx = struct();
[opx.Al, opx.AlName] = getVarOrError(data_opx, ...
    {'Al_cation_apfu','Al_cation','Al'}, 'Opx');
[opx.Mg, opx.MgName] = getVarOrError(data_opx, ...
    {'Mg_cation_apfu','Mg_cation','Mg'}, 'Opx');
[opx.Fe, opx.FeName] = getVarOrError(data_opx, ...
    {'Fe_cation_apfu','Fe2_cation_apfu','Fe_cation','Fe2_cation','Fe','Fe2'}, ...
    'Opx');
[opx.Na, opx.NaName] = getVarOrError(data_opx, ...
    {'Na_cation_apfu','Na_cation','Na'}, 'Opx');
[opx.Ti, opx.TiName] = getVarOrError(data_opx, ...
    {'Ti_cation_apfu','Ti_cation','Ti'}, 'Opx');
[opx.Cr, opx.CrName] = getVarOrError(data_opx, ...
    {'Cr_cation_apfu','Cr_cation','Cr'}, 'Opx');
[opx.Ca, opx.CaName] = getVarOrError(data_opx, ...
    {'Ca_cation_apfu','Ca_cation','Ca'}, 'Opx');
opx.Mn = getVarOrNaN(data_opx, ...
    {'Mn_cation_apfu','Mn_cation','Mn'});

end

function [value, variableName] = getVarOrError(tbl, acceptedNames, mineralLabel)
% getVarOrError
% Retrieve a required scalar table variable using the first accepted name
% found. Existing NaN is retained unchanged.

value = [];
variableName = "";

for i = 1:numel(acceptedNames)
    candidateName = acceptedNames{i};
    if ismember(candidateName, tbl.Properties.VariableNames)
        value = tbl.(candidateName);
        variableName = string(candidateName);
        break;
    end
end

if isempty(value)
    error('%s table must contain one of the variables: %s', ...
        mineralLabel, strjoin(acceptedNames, ', '));
end
if ~isscalar(value) || ~isnumeric(value)
    error('%s variable %s must be a numeric scalar in a 1-row table.', ...
        mineralLabel, char(variableName));
end

end

function value = getVarOrNaN(tbl, acceptedNames)
% getVarOrNaN
% Retrieve an optional scalar variable. A missing optional variable is
% represented by NaN, never by zero.

value = NaN;

for i = 1:numel(acceptedNames)
    candidateName = acceptedNames{i};
    if ismember(candidateName, tbl.Properties.VariableNames)
        value = tbl.(candidateName);
        if ~isscalar(value) || ~isnumeric(value)
            error('Optional variable %s must be a numeric scalar in a 1-row table.', ...
                candidateName);
        end
        return;
    end
end

end
