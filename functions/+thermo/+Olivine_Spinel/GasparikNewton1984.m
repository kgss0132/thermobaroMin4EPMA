function results = GasparikNewton1984(rawdata_struct, P_kbar)
% functions/+thermo/+Olivine_Spinel/GasparikNewton1984.m
% Tested with MATLAB R2024b
%
% Orthopyroxene-Spinel-Olivine geothermometer
% Gasparik, T., Newton, R.C. (1984)
% Contributions to Mineralogy and Petrology, 85, 186-196
% DOI: https://doi.org/10.1007/BF00371708
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Olivine analysis, one
% Orthopyroxene analysis, and one Spinel analysis from tables and calculates
% temperature using the natural spinel-lherzolite formulation derived from
% the reaction:
%
%   En + Sp = MgTs + Fo
%
% The function is designed for repeated calculations: after each run it asks
% whether another Ol-Opx-Sp combination should be calculated and stores all
% result blocks in a single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Gasparik and Newton (1984) reversed the alumina contents of orthopyroxene
% coexisting with spinel and forsterite in the MgO-Al2O3-SiO2 (MAS) system
% at 15 pressure-temperature conditions over the following direct
% experimental range:
%
%   Temperature : 1030-1600 degreeC
%   Pressure    : 10-28 kbar (1.0-2.8 GPa)
%   System      : MgO-Al2O3-SiO2
%   Assemblage  : orthopyroxene + spinel + forsterite
%
% These experimental conditions are stated in the abstract on p. 186 and
% listed in Table 3 on p. 189. When the related reversal data of Danckwerth
% and Newton (1978) and Perkins et al. (1981) are also included, the broader
% thermodynamic dataset covers approximately 900-1600 degreeC and
% 10-40 kbar (Conclusion, p. 195). The direct Opx-Sp-Fo calibration range
% of 1030-1600 degreeC and 10-28 kbar is used for the calibration warnings
% in this implementation.
%
% For natural spinel lherzolites, Gasparik and Newton (1984) used:
%
%   T(degreeC) = 3857*K + 443 + P(kbar)
%
% and noted that the pressure effect is approximately 1 degreeC/kbar
% (Application to geothermometry, p. 193). They assumed 10 kbar for their
% natural examples because this lies near the middle of the approximate
% 6-15 kbar spinel-lherzolite stability interval at about 1000 degreeC
% (p. 193). They later described an approximate possible recrystallization
% pressure range of 6-18 kbar for natural spinel lherzolites, bounded by
% plagioclase-lherzolite stability at lower pressure and garnet-lherzolite
% stability at higher pressure (p. 195).
%
% The thermometer should be applied only to natural systems that closely
% approach the MAS system (Conclusion, p. 195). Application to strongly
% Fe-, Cr-, Ti-, Na-, or Ca-rich compositions requires caution because the
% natural-system calculation depends on activity-composition models.
%
% Orthopyroxene MgTs and En components are calculated with corrections for
% Na, Ti, Cr, Ca, and Mn following Eqs. (23)-(25) on p. 193. Spinel
% MgAl2O4 activity is calculated with the simplified Sack (1982) model
% presented in Eqs. (26)-(27) on p. 193. Gasparik and Newton (1984) noted
% that MgAl2O4, FeCr2O4, and Fe3O4 commonly account for at least about 95%
% of lherzolite spinels. Compositions with substantial additional spinel
% components may therefore be poorly represented by this simplification.
%
% Texturally and chemically equilibrated minerals should be paired. The
% natural examples on pp. 193-195 show that orthopyroxene porphyroclasts
% and matrix/neoblast assemblages may record distinctly different
% temperatures. Mixing minerals from different generations can therefore
% produce geologically misleading results.
%
% Natural application results reported in the abstract are approximately
% 800-1350 degreeC for alpine peridotites and 850-1130 degreeC for volcanic
% upper-mantle inclusions (p. 186). These are application results, not the
% direct experimental calibration limits.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside the direct experimental range of
%      10-28 kbar;
%   2) input pressure is outside the approximate 6-18 kbar natural
%      spinel-lherzolite stability range;
%   3) a finite calculated temperature is outside the direct experimental
%      range of 1030-1600 degreeC;
%   4) thermometer inputs contain NaN; or
%   5) calculated temperatures contain NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%   rawdata_struct.Opx     : table
%   rawdata_struct.Spinel  : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% Required normalized-cation variable names:
%
%   Olivine table variables:
%     Mg_cation_apfu
%     Fe_cation_apfu        % Fe2+ in olivine (assumed)
%     Mn_cation_apfu
%
%   Orthopyroxene table variables:
%     Al_cation_apfu
%     Mg_cation_apfu
%     Na_cation_apfu
%     Ti_cation_apfu
%     Cr_cation_apfu
%     Ca_cation_apfu
%     Mn_cation_apfu
%
%   Spinel table variables:
%     Al_cation_apfu
%     Cr_cation_apfu
%     Mg_cation_apfu
%     Fe3_cation_apfu
%
% Variable aliases such as Mg_cation, Mg, Fe2_cation_apfu, Fe2, and
% corresponding element-name variants are also accepted by the local
% data-retrieval functions. All variables listed above are required; absent
% variables are not replaced by zero.
%
% All finite values in the required mineral-composition variables above must
% be strictly greater than zero. NaN values are retained as missing values,
% propagated through the calculation, and reported by non-stopping fprintf
% warnings after the calculated temperature is displayed.
%
% Expected cation bases:
%   Olivine       : cations per formula unit, commonly normalized to 4 O
%   Orthopyroxene : cations per 6 O
%   Spinel        : cations per 4 O
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% 1) Olivine activity
%     XMg_ol = Mg_ol / (Mg_ol + Fe2_ol + Mn_ol)
%     aFo_ol = XMg_ol^2
%
% 2) Orthopyroxene components and ideal two-site activity ratio
%     XMgTs = (Al - Na - 2*Ti - Cr) / 2
%     XEn   = Mg - Ca - Cr/2 - Mn - XMgTs
%
%     XAl_M1 = XMgTs + Na
%     XMg_M1 = XEn + Ca + Mn
%
%     aMgTs/aEn = XAl_M1^2 / XMg_M1
%
% 3) Simplified spinel components
%     X2 = MgAl2O4
%     X3 = FeCr2O4
%     X5 = Fe3O4
%
%     aSp is calculated from the simplified Sack (1982) expression quoted
%     by Gasparik and Newton (1984).
%
% 4) Equilibrium constant and temperature
%     K = aFo_ol * (aMgTs/aEn) / aSp
%
%     T(degreeC) = 3857*K + 443 + P(kbar)
%
% The spinel activity depends on temperature, so temperature and spinel
% activity are updated iteratively. This implementation retains the three
% update cycles and final activity/temperature update of the original code.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = GasparikNewton1984(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Olivine, Opx, and Spinel tables
%   P_kbar         : pressure in kbar (finite, non-negative numeric scalar
%                    or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Ol-Opx-Sp combination
%

%% Input validation
% Basic argument checks prevent silent failures caused by missing arguments
% or invalid pressure values.
if nargin < 2
    error('GasparikNewton1984 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The tables are not
% modified; selected rows are read during each calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end
if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end
if ~isfield(rawdata_struct, 'Spinel') || ~istable(rawdata_struct.Spinel)
    error('rawdata_struct must contain table: rawdata_struct.Spinel');
end

dataset_ol = rawdata_struct.Olivine;
dataset_opx = rawdata_struct.Opx;
dataset_sp = rawdata_struct.Spinel;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated table concatenation inside the loop is avoided because it forces
% MATLAB to repeatedly reallocate and copy the entire results table.
%
% The cell buffer is preallocated and doubled only when necessary. After the
% interactive loop finishes, all blocks are concatenated once with vertcat.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct Opx-Sp-Fo experimental limits reported by Gasparik and Newton
% (1984), pp. 186 and 189.
calibrationT_min_degC = 1030;
calibrationT_max_degC = 1600;
calibrationP_min_kbar = 10;
calibrationP_max_kbar = 28;

% Approximate natural spinel-lherzolite pressure interval discussed on
% p. 195. This is an application/phase-stability range, not the direct
% experimental calibration range.
naturalP_min_kbar = 6;
naturalP_max_kbar = 18;

% Pressure is common to all mineral combinations in this function call.
% Pressure messages are therefore printed only once, after the first result.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureOutsideNaturalRange = ...
    P_kbar < naturalP_min_kbar | P_kbar > naturalP_max_kbar;

calibrationPressureWarningIssued = false;
naturalPressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-6) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or selects
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Olivine) ===');

while true
    % ----- Olivine selection -----
    % Assumption: the first column stores an identifier displayed to the user.
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

    % ----- Orthopyroxene selection -----
    disp('=== Step 4: Selecting a data code from the list (Orthopyroxene) ===');

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

    % ----- Spinel selection -----
    disp('=== Step 5: Selecting a data code from the list (Spinel) ===');

    dataCodes_sp = dataset_sp{:, 1};

    [selectedIdx_sp, ok] = listdlg( ...
        'PromptString', 'Please select the Spinel data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_sp)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_sp)
        disp('Selection canceled');
        break;
    end

    selectedCode_sp = dataCodes_sp(selectedIdx_sp);
    disp(['Spinel selected: ' char(string(selectedCode_sp))]);

    % ----- Calculation -----
    % Mineral phases are selected independently. Do not assume that row
    % indices correspond among the three datasets.
    disp('=== Step 6: Calculating the temperature ===');

    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_opx = dataset_opx(selectedIdx_opx, :);
    selectedData_sp = dataset_sp(selectedIdx_sp, :);

    % NaN does not stop calculation. Names of affected inputs are collected
    % before calculation and reported immediately after the result.
    nanInputNames = findNaNInputs( ...
        selectedData_ol, selectedData_opx, selectedData_sp);

    % All finite compositional values used by the thermometer must be
    % strictly positive. NaN is intentionally allowed so that it propagates
    % through the calculation and is reported by non-stopping warnings.
    validatePositiveInputs( ...
        selectedData_ol, selectedData_opx, selectedData_sp);

    row = calcTemp( ...
        selectedData_ol, selectedData_opx, selectedData_sp, P_kbar);

    % Repeat identifiers for every pressure point in the result block.
    row.dataCode_ol = repmat(string(selectedCode_ol), height(row), 1);
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row.dataCode_sp = repmat(string(selectedCode_sp), height(row), 1);
    row = movevars(row, ...
        {'dataCode_ol','dataCode_opx','dataCode_sp'}, 'Before', 1);

    % Store this result as one table block. If the preallocated buffer is
    % full, double its capacity only occasionally rather than reallocating
    % the complete results table on every calculation.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        oldCapacity = numel(resultBlocks);
        resultBlocks(oldCapacity + 1 : 2 * oldCapacity, 1) = ...
            cell(oldCapacity, 1);
    end

    resultBlocks{nResultBlocks} = row;

    % Echo calculated temperatures for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    selectedCombinationText = [ ...
        char(string(selectedCode_ol)) ' & ' ...
        char(string(selectedCode_opx)) ' & ' ...
        char(string(selectedCode_sp))];

    if height(row) == 1
        disp([selectedCombinationText ': ' num2str(row.T_deg) ' degreeC']);
    else
        disp([selectedCombinationText ': ' ...
            num2str(row.T_deg(1)) ' to ' ...
            num2str(row.T_deg(end)) ' degreeC']);
    end

    % Warn once when any pressure lies outside the direct experimental
    % calibration range. Calculation is not stopped.
    if any(pressureOutsideCalibration) && ...
            ~calibrationPressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the direct experimental ' ...
             'calibration range of Gasparik and Newton (1984): ' ...
             '10-28 kbar (1.0-2.8 GPa). %d of %d pressure point(s) ' ...
             'are outside the range; input range = %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        calibrationPressureWarningIssued = true;
    end

    % Print a separate application note for pressures outside the
    % approximate natural spinel-lherzolite stability range.
    if any(pressureOutsideNaturalRange) && ...
            ~naturalPressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate ' ...
             '6-18 kbar natural spinel-lherzolite stability range ' ...
             'discussed by Gasparik and Newton (1984). %d of %d ' ...
             'pressure point(s) are outside this range; input range = ' ...
             '%.4g-%.4g kbar. This is an application/phase-stability ' ...
             'range, not the direct experimental calibration range.\n'], ...
            sum(pressureOutsideNaturalRange), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        naturalPressureWarningIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the direct
    % experimental range. NaN and Inf are handled separately below.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_deg < calibrationT_min_degC | ...
         row.T_deg > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the direct ' ...
             'experimental calibration range of Gasparik and Newton ' ...
             '(1984): 1030-1600 degreeC. %d of %d finite temperature ' ...
             'point(s) are outside the range; calculated finite range = ' ...
             '%.4g-%.4g degreeC for %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            selectedCombinationText);
    end

    % Print a non-stopping warning when a present thermometer input contains
    % NaN. fprintf is used so the message remains in the command-window log.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for ' ...
             '%s: %s.\n' ...
             '         The calculation was continued, but the calculated ' ...
             'temperature may be NaN.\n'], ...
            selectedCombinationText, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain a result-based check for NaN/Inf produced by any cause,
    % including invalid activity terms or division by zero.
    invalidTemperature = ~isfinite(row.T_deg);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for ' ...
             '%s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            selectedCombinationText, ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    % Ask whether another mineral combination should be calculated.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'GasparikNewton1984', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate result blocks once after the interactive loop. Return an empty
% table when no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs( ...
    data_olivine, data_opx, data_spinel)
% findNaNInputs
% Return names of required thermometer variables containing NaN.
% NaN values are reported but do not stop calculation. Missing required
% variables are handled separately by getVar. The output buffer is
% preallocated so it does not grow during the loop.

phaseTables = {data_olivine, data_opx, data_spinel};
phaseNames = { ...
    'Olivine', 'Olivine', 'Olivine', ...
    'Opx', 'Opx', 'Opx', 'Opx', 'Opx', 'Opx', 'Opx', ...
    'Spinel', 'Spinel', 'Spinel', 'Spinel'};
tableIndices = [ ...
    1, 1, 1, ...
    2, 2, 2, 2, 2, 2, 2, ...
    3, 3, 3, 3];
canonicalNames = { ...
    'Mg_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Al_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Al_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Fe3_cation_apfu'};
aliasGroups = { ...
    {'Mg_cation_apfu','Mg_cation','Mg'}, ...
    {'Fe2_cation_apfu','Fe_cation_apfu','Fe2_cation', ...
        'Fe_cation','Fe2','Fe'}, ...
    {'Mn_cation_apfu','Mn_cation','Mn'}, ...
    {'Al_cation_apfu','Al_cation','Al'}, ...
    {'Mg_cation_apfu','Mg_cation','Mg'}, ...
    {'Na_cation_apfu','Na_cation','Na'}, ...
    {'Ti_cation_apfu','Ti_cation','Ti'}, ...
    {'Cr_cation_apfu','Cr_cation','Cr'}, ...
    {'Ca_cation_apfu','Ca_cation','Ca'}, ...
    {'Mn_cation_apfu','Mn_cation','Mn'}, ...
    {'Al_cation_apfu','Al_cation','Al'}, ...
    {'Cr_cation_apfu','Cr_cation','Cr'}, ...
    {'Mg_cation_apfu','Mg_cation','Mg'}, ...
    {'Fe3_cation_apfu','Fe3_cation','Fe3'}};

nVariables = numel(aliasGroups);
nanInputBuffer = strings(nVariables, 1);
nNaNInputs = 0;

for i = 1:nVariables
    tbl = phaseTables{tableIndices(i)};
    [isPresent, variableValue] = getVarIfPresent(tbl, aliasGroups{i});

    if isPresent && any(isnan(variableValue(:)))
        nNaNInputs = nNaNInputs + 1;
        nanInputBuffer(nNaNInputs) = ...
            string(phaseNames{i}) + "." + string(canonicalNames{i});
    end
end

nanInputNames = nanInputBuffer(1:nNaNInputs);

end


function validatePositiveInputs(data_olivine, data_opx, data_spinel)
% validatePositiveInputs
% Stop the calculation when any finite required mineral-composition value
% is zero or negative. NaN is intentionally excluded from this check so that
% it remains NaN and propagates through the thermometer calculation.

phaseTables = {data_olivine, data_opx, data_spinel};
phaseNames = { ...
    'Olivine', 'Olivine', 'Olivine', ...
    'Opx', 'Opx', 'Opx', 'Opx', 'Opx', 'Opx', 'Opx', ...
    'Spinel', 'Spinel', 'Spinel', 'Spinel'};
tableIndices = [ ...
    1, 1, 1, ...
    2, 2, 2, 2, 2, 2, 2, ...
    3, 3, 3, 3];
canonicalNames = { ...
    'Mg_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Al_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Al_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Fe3_cation_apfu'};
aliasGroups = { ...
    {'Mg_cation_apfu','Mg_cation','Mg'}, ...
    {'Fe2_cation_apfu','Fe_cation_apfu','Fe2_cation', ...
        'Fe_cation','Fe2','Fe'}, ...
    {'Mn_cation_apfu','Mn_cation','Mn'}, ...
    {'Al_cation_apfu','Al_cation','Al'}, ...
    {'Mg_cation_apfu','Mg_cation','Mg'}, ...
    {'Na_cation_apfu','Na_cation','Na'}, ...
    {'Ti_cation_apfu','Ti_cation','Ti'}, ...
    {'Cr_cation_apfu','Cr_cation','Cr'}, ...
    {'Ca_cation_apfu','Ca_cation','Ca'}, ...
    {'Mn_cation_apfu','Mn_cation','Mn'}, ...
    {'Al_cation_apfu','Al_cation','Al'}, ...
    {'Cr_cation_apfu','Cr_cation','Cr'}, ...
    {'Mg_cation_apfu','Mg_cation','Mg'}, ...
    {'Fe3_cation_apfu','Fe3_cation','Fe3'}};

nVariables = numel(aliasGroups);
invalidInputBuffer = strings(nVariables, 1);
nInvalidInputs = 0;

for i = 1:nVariables
    tbl = phaseTables{tableIndices(i)};
    variableValue = getVar(tbl, aliasGroups{i});

    if any(isinf(variableValue(:)))
        error(['GasparikNewton1984: required mineral-composition ' ...
            'values must not contain Inf. Inf was found in %s.%s.'], ...
            phaseNames{i}, canonicalNames{i});
    end

    if any(isfinite(variableValue(:)) & variableValue(:) <= 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            string(phaseNames{i}) + "." + string(canonicalNames{i});
    end
end

invalidInputNames = invalidInputBuffer(1:nInvalidInputs);

if ~isempty(invalidInputNames)
    error(['GasparikNewton1984: all finite mineral-composition ' ...
        'values used by the thermometer must be > 0. Zero or negative ' ...
        'value(s) were found in: ' ...
        char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_olivine, data_opx, data_spinel, P_kbar)
% calcTemp
% Calculate one temperature per input pressure for a single selected
% Olivine-Orthopyroxene-Spinel combination.
%
% Inputs:
%   data_olivine : one-row table of olivine cations
%   data_opx     : one-row table of orthopyroxene cations
%   data_spinel  : one-row table of spinel cations
%   P_kbar       : pressure in kbar; scalar or column vector
%
% Output:
%   row : table with one row per pressure value, including cations,
%         intermediate activity terms, applicability flags, and temperature

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% -------------------------------------------------------------------------
% Extract selected cations
% -------------------------------------------------------------------------
Mg_ol_scalar = getVar( ...
    data_olivine, {'Mg_cation_apfu','Mg_cation','Mg'});
Fe2_ol_scalar = getVar( ...
    data_olivine, {'Fe2_cation_apfu','Fe_cation_apfu', ...
    'Fe2_cation','Fe_cation','Fe2','Fe'});
Mn_ol_scalar = getVar( ...
    data_olivine, {'Mn_cation_apfu','Mn_cation','Mn'});

Al_opx_scalar = getVar( ...
    data_opx, {'Al_cation_apfu','Al_cation','Al'});
Mg_opx_scalar = getVar( ...
    data_opx, {'Mg_cation_apfu','Mg_cation','Mg'});
Na_opx_scalar = getVar( ...
    data_opx, {'Na_cation_apfu','Na_cation','Na'});
Ti_opx_scalar = getVar( ...
    data_opx, {'Ti_cation_apfu','Ti_cation','Ti'});
Cr_opx_scalar = getVar( ...
    data_opx, {'Cr_cation_apfu','Cr_cation','Cr'});
Ca_opx_scalar = getVar( ...
    data_opx, {'Ca_cation_apfu','Ca_cation','Ca'});
Mn_opx_scalar = getVar( ...
    data_opx, {'Mn_cation_apfu','Mn_cation','Mn'});

Al_sp_scalar = getVar( ...
    data_spinel, {'Al_cation_apfu','Al_cation','Al'});
Cr_sp_scalar = getVar( ...
    data_spinel, {'Cr_cation_apfu','Cr_cation','Cr'});
Mg_sp_scalar = getVar( ...
    data_spinel, {'Mg_cation_apfu','Mg_cation','Mg'});
Fe3_sp_scalar = getVar( ...
    data_spinel, {'Fe3_cation_apfu','Fe3_cation','Fe3'});

% Each selected table entry must contain one scalar value per variable.
selectedValues = [ ...
    Mg_ol_scalar, Fe2_ol_scalar, Mn_ol_scalar, ...
    Al_opx_scalar, Mg_opx_scalar, Na_opx_scalar, Ti_opx_scalar, ...
    Cr_opx_scalar, Ca_opx_scalar, Mn_opx_scalar, ...
    Al_sp_scalar, Cr_sp_scalar, Mg_sp_scalar, Fe3_sp_scalar];

if numel(selectedValues) ~= 14
    error(['Each selected mineral table entry must contain one scalar ' ...
        'value per compositional variable.']);
end

% Infinite inputs are rejected. NaN inputs are retained so a non-stopping
% warning can be printed after the calculation result.
if any(isinf(selectedValues))
    error('Thermometer compositional inputs must not contain Inf.');
end

if any(isfinite(selectedValues) & selectedValues <= 0)
    error(['GasparikNewton1984: all finite mineral-composition ' ...
        'values used by the thermometer must be positive.']);
end

% Repeat the selected composition for every input pressure.
Mg_ol = repmat(Mg_ol_scalar, nP, 1);
Fe2_ol = repmat(Fe2_ol_scalar, nP, 1);
Mn_ol = repmat(Mn_ol_scalar, nP, 1);

Al_opx = repmat(Al_opx_scalar, nP, 1);
Mg_opx = repmat(Mg_opx_scalar, nP, 1);
Na_opx = repmat(Na_opx_scalar, nP, 1);
Ti_opx = repmat(Ti_opx_scalar, nP, 1);
Cr_opx = repmat(Cr_opx_scalar, nP, 1);
Ca_opx = repmat(Ca_opx_scalar, nP, 1);
Mn_opx = repmat(Mn_opx_scalar, nP, 1);

Al_sp = repmat(Al_sp_scalar, nP, 1);
Cr_sp = repmat(Cr_sp_scalar, nP, 1);
Mg_sp = repmat(Mg_sp_scalar, nP, 1);
Fe3_sp = repmat(Fe3_sp_scalar, nP, 1);

% -------------------------------------------------------------------------
% Olivine activity
% -------------------------------------------------------------------------
XMg_ol = Mg_ol ./ (Mg_ol + Fe2_ol + Mn_ol);
aFo_ol = XMg_ol.^2;

% -------------------------------------------------------------------------
% Orthopyroxene components and ideal two-site activity ratio
% Gasparik and Newton (1984), Eqs. (23)-(25), p. 193
% -------------------------------------------------------------------------
XMgTs_opx = (Al_opx - Na_opx - 2 .* Ti_opx - Cr_opx) ./ 2;
XEn_opx = Mg_opx - Ca_opx - Cr_opx ./ 2 - Mn_opx - XMgTs_opx;

XAl_M1_opx = XMgTs_opx + Na_opx;
XMg_M1_opx = XEn_opx + Ca_opx + Mn_opx;

invalidFiniteOpx = ...
    (isfinite(XMgTs_opx) & XMgTs_opx <= 0) | ...
    (isfinite(XEn_opx) & XEn_opx <= 0) | ...
    (isfinite(XAl_M1_opx) & XAl_M1_opx <= 0) | ...
    (isfinite(XMg_M1_opx) & XMg_M1_opx <= 0);

if any(invalidFiniteOpx)
    error(['Calculated finite orthopyroxene component or M1-site terms ' ...
        'are non-positive. Check Opx cation data.']);
end

aMgTs_over_aEn = XAl_M1_opx.^2 ./ XMg_M1_opx;

% -------------------------------------------------------------------------
% Simplified spinel components
% -------------------------------------------------------------------------
X2_raw = Al_sp ./ 2;
X3_raw = Cr_sp ./ 2;
X5_raw = Fe3_sp ./ 2;

sumX_sp = X2_raw + X3_raw + X5_raw;

if any(isfinite(sumX_sp) & sumX_sp <= 0)
    error('Finite spinel component sums must be positive.');
end

X2_sp = X2_raw ./ sumX_sp;
X3_sp = X3_raw ./ sumX_sp;
X5_sp = X5_raw ./ sumX_sp;

invalidFiniteSpinel = ...
    (isfinite(X2_sp) & X2_sp <= 0) | ...
    (isfinite(X3_sp) & X3_sp <= 0) | ...
    (isfinite(X5_sp) & X5_sp <= 0) | ...
    (isfinite(1 - X5_sp) & (1 - X5_sp) <= 0) | ...
    (isfinite(1 - X3_sp - X5_sp) & ...
        (1 - X3_sp - X5_sp) <= 0);

if any(invalidFiniteSpinel)
    error(['Calculated finite spinel component fractions are outside ' ...
        'the simplified Sack activity model.']);
end

% -------------------------------------------------------------------------
% Iterative spinel activity and temperature solution
% -------------------------------------------------------------------------
T_deg_iter = repmat(1000, nP, 1);

for iter = 1:3
    T_K_iter = T_deg_iter + 273.15;

    aSp_sp = calcSpinelActivitySack( ...
        X2_sp, X3_sp, X5_sp, T_K_iter);

    if any(isfinite(aSp_sp) & aSp_sp <= 0)
        error('Calculated finite spinel activities must be positive.');
    end

    K_eq = aFo_ol .* aMgTs_over_aEn ./ aSp_sp;

    if any(isfinite(K_eq) & K_eq <= 0)
        error('Calculated finite equilibrium constants K must be positive.');
    end

    T_deg_iter = 3857 .* K_eq + 443 + P_kbar;
end

% Final activity and temperature update.
T_K_iter = T_deg_iter + 273.15;
aSp_sp = calcSpinelActivitySack(X2_sp, X3_sp, X5_sp, T_K_iter);
K_eq = aFo_ol .* aMgTs_over_aEn ./ aSp_sp;
T_deg = 3857 .* K_eq + 443 + P_kbar;
T_K = T_deg + 273.15;

% -------------------------------------------------------------------------
% Applicability flags
% -------------------------------------------------------------------------
% Direct Opx-Sp-Fo experimental range.
is_within_experimental_pressure_range = ...
    P_kbar >= 10 & P_kbar <= 28;
is_within_experimental_temperature_range = ...
    T_deg >= 1030 & T_deg <= 1600;

% Existing natural-application flags are retained for compatibility.
% The 6-18 kbar interval is the approximate natural spinel-lherzolite
% pressure range discussed on p. 195. The 800-1350 degreeC interval spans
% the natural alpine-peridotite results reported in the abstract and is not
% a direct experimental calibration range.
is_spinel_lherzolite_range = P_kbar >= 6 & P_kbar <= 18;
is_T_reasonable = T_deg >= 800 & T_deg <= 1350;
is_recommended = is_spinel_lherzolite_range & is_T_reasonable;

% -------------------------------------------------------------------------
% Pack outputs
% -------------------------------------------------------------------------
row.Mg_ol = Mg_ol;
row.Fe2_ol = Fe2_ol;
row.Mn_ol = Mn_ol;
row.XMg_ol = XMg_ol;
row.aFo_ol = aFo_ol;

row.Al_opx = Al_opx;
row.Mg_opx = Mg_opx;
row.Na_opx = Na_opx;
row.Ti_opx = Ti_opx;
row.Cr_opx = Cr_opx;
row.Ca_opx = Ca_opx;
row.Mn_opx = Mn_opx;

row.XMgTs_opx = XMgTs_opx;
row.XEn_opx = XEn_opx;
row.XAl_M1_opx = XAl_M1_opx;
row.XMg_M1_opx = XMg_M1_opx;
row.aMgTs_over_aEn = aMgTs_over_aEn;

row.Al_sp = Al_sp;
row.Cr_sp = Cr_sp;
row.Mg_sp = Mg_sp;
row.Fe3_sp = Fe3_sp;

row.X2_MgAl2O4_sp = X2_sp;
row.X3_FeCr2O4_sp = X3_sp;
row.X5_Fe3O4_sp = X5_sp;
row.aSp_sp = aSp_sp;

row.K_eq = K_eq;

row.is_within_experimental_pressure_range = ...
    is_within_experimental_pressure_range;
row.is_within_experimental_temperature_range = ...
    is_within_experimental_temperature_range;
row.is_spinel_lherzolite_pressure_range = ...
    is_spinel_lherzolite_range;
row.is_T_reasonable = is_T_reasonable;
row.is_recommended = is_recommended;

row.T_K = T_K;
row.T_deg = T_deg;

end

function aSp = calcSpinelActivitySack(X2, X3, X5, T_K)
% calcSpinelActivitySack
% Calculate MgAl2O4 spinel activity using the simplified Sack (1982)
% formulation quoted by Gasparik and Newton (1984), pp. 193-194.

R = 8.31446261815324;

% Parameters are given in kJ/mol, except the explicitly coded 14713 J/mol
% term retained from the original implementation.
dmu23 = (38.87 - 0.0164 .* T_K) .* 1000;
dmu25 = (30.12 - 0.0054 .* T_K) .* 1000;

W13 = 28.66 .* 1000;
W31 = 19.25 .* 1000;
W15 = 55.23 .* 1000;
W53 = 37.66 .* 1000;
W51 = 58.58 .* 1000;

idealTerm = X2 .* (1 - X5) .* (1 - X3 - X5).^2;

if any(isfinite(idealTerm) & idealTerm <= 0)
    error('Finite ideal terms in the spinel activity model must be positive.');
end

RTln_a2 = R .* T_K .* log(idealTerm) ...
    + dmu23 .* (1 - X2) .* X3 ...
    + dmu25 .* (1 - X2) .* X5 ...
    + (W31 + 2 .* (14713 - W31) .* (1 - X3)) .* X3.^2 ...
    + (W51 + 2 .* (W15 - W51) .* (1 - X5)) .* X5.^2 ...
    - (W53 - W13 - W51) .* X3 .* X5;

aSp = exp(RTln_a2 ./ (R .* T_K));

end

function val = getVar(tbl, names)
% getVar
% Return the first matching required variable from a one-row table.

[isPresent, val] = getVarIfPresent(tbl, names);

if ~isPresent
    error('Required variable is missing. Accepted names include: %s', ...
        strjoin(names, ', '));
end

end


function [isPresent, val] = getVarIfPresent(tbl, names)
% getVarIfPresent
% Find the first accepted variable name in a table.

isPresent = false;
val = [];

for i = 1:numel(names)
    if ismember(names{i}, tbl.Properties.VariableNames)
        val = tbl.(names{i});
        isPresent = true;
        return
    end
end

end
