function results = EvansFrost1975(rawdata_struct, P_kbar)
% functions/+thermo/+Olivine_Spinel/EvansFrost1975.m
% Tested with MATLAB R2024b
%
% Semi-quantitative Fe-Mg exchange thermometer between Olivine and Spinel
% Evans, B.W., Frost, B.R. (1975)
% Geochimica et Cosmochimica Acta, 39, 959-972
% DOI: https://doi.org/10.1016/0016-7037(75)90041-1
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis and one Spinel
% analysis (selected by the user from tables) and estimates temperature
% using the graphical Fe-Mg exchange calibration of Evans and Frost (1975).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Ol-Sp pair, and appends results into a
% single output table.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Evans and Frost (1975) did not present a single experimentally calibrated
% thermometer equation. Instead, they proposed a tentative, revised,
% graphical calibration of olivine-chromite Fe-Mg exchange using natural
% mineral pairs (abstract, p. 959; discussion and Fig. 8, pp. 969-971).
%
% The two principal temperature reference lines used in their graphical
% calibration are approximately:
%
%   Temperature : 700-1200 degreeC
%   Lower anchor: approximately 700 degreeC for high-grade metamorphic
%                 chlorite-enstatite-olivine-spinel assemblages
%   Upper anchor: approximately 1200 degreeC for olivine-chromite pairs from
%                 basaltic pumice and oceanic basalt
%   Pressure    : no explicit numerical pressure calibration range reported
%
% The approximately 700 degreeC and 1200 degreeC reference isotherms and the
% interpolation concept are described on pp. 969-970 and shown in Fig. 8.
% Evans and Frost (1975) state that intermediate isotherms would be spaced
% approximately in proportion to 1/T only if reaction enthalpy were
% independent of temperature (p. 969). The calibration should therefore be
% regarded as semi-quantitative rather than as a rigorously calibrated
% absolute thermometer.
%
% IMPORTANT APPLICATION NOTES:
% - Use equilibrium olivine-chromite pairs. The study emphasizes thoroughly
%   equilibrated assemblages and avoids retrograde or metastable chrome-
%   spinel growth (pp. 959-960).
% - The approximately 700 degreeC metamorphic trend is not strictly
%   isothermal because spinel composition also depends on pressure,
%   temperature, and fluid conditions (p. 967).
% - Spinel solid-solution non-ideality produces curvature, particularly for
%   Al-poor magnetite-chromite compositions at lower temperatures
%   (pp. 967-969; Fig. 7).
% - More aluminous spinels from peridotite nodules show greater scatter than
%   more chromiferous spinels (p. 969).
% - Fe2+ and Fe3+ in spinel were calculated by ideal spinel stoichiometry;
%   analytical or ferric-iron estimation errors therefore propagate into the
%   exchange coefficient and temperature estimate (p. 960).
% - No explicit pressure range was calibrated. Several high-pressure natural
%   pairs were plotted without pressure correction, and the paper illustrates
%   a pressure correction for a 60 kbar diamond inclusion using 298 K data
%   (pp. 969-970; Fig. 8). Pressure effects are therefore not quantified by
%   the implemented interpolation.
% - Evans and Frost (1975) explicitly state that the empirical graph is not a
%   substitute for experimental determination of the temperature and
%   compositional dependence of the exchange coefficient (pp. 970-971).
%
% This implementation therefore issues non-stopping fprintf messages when:
%   1) a finite calculated temperature is outside 700-1200 degreeC, or
%   2) pressure is supplied. Because Evans and Frost (1975) did not define a
%      numerical pressure calibration range, the pressure message states that
%      pressure applicability cannot be classified and that pressure is not
%      used in the calculation.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%   rawdata_struct.Spinel  : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following variable names (normalized cations):
%
%   Olivine table variables:
%     Mg_cation_apfu
%     Fe_cation_apfu         % Fe2+ in olivine (assumed)
%
%   Spinel table variables:
%     Mg_cation_apfu
%     Fe_cation_apfu         % Fe2+ in spinel
%     Fe3_cation_apfu        % Fe3+ in spinel
%     Cr_cation_apfu
%     Al_cation_apfu
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (graphical interpolation implemented here)
%
% 1) Mole fractions / indices
%   Olivine:
%     XMg_ol = Mg_ol / (Mg_ol + Fe2_ol)
%     XFe_ol = Fe2_ol / (Mg_ol + Fe2_ol)
%
%   Spinel:
%     XMg_sp  = Mg_sp / (Mg_sp + Fe2_sp)
%     XFe2_sp = Fe2_sp / (Mg_sp + Fe2_sp)
%     YCr_sp  = Cr_sp / (Cr_sp + Al_sp + Fe3_sp)
%     YFe3_sp = Fe3_sp / (Cr_sp + Al_sp + Fe3_sp)
%
% 2) Fe-Mg exchange coefficient and ferric-iron normalization
%     KD = (XMg_ol * XFe2_sp) / (XFe_ol * XMg_sp)
%
%     ln(KD*) = ln(KD) + (0.05 - YFe3_sp) * 4
%
% 3) Approximate reference isotherms used by this implementation
%     700 degreeC : ln(KD*) =  0.85 + 2.94 * YCr_sp
%     1200 degreeC: ln(KD*) = -0.05 + 2.94 * YCr_sp
%
% 4) Temperature interpolation
%     T(degreeC) = 700 + 500 *
%                  [ln(KD*)_700 - ln(KD*)] /
%                  [ln(KD*)_700 - ln(KD*)_1200]
%
% Notes:
% - The reference-line equations and interpolation are a numerical
%   representation of the semi-quantitative graphical calibration, not a
%   published closed-form thermometer equation.
% - Pressure is not used in the calculation.
% - P_kbar is retained for compatibility with the common thermometer
%   interface and may be a scalar or vector. One output row is returned for
%   every supplied pressure value.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = EvansFrost1975(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Olivine and Spinel tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Ol-Sp pair. Because pressure is not used in the
%             thermometer, repeated rows for one pair have identical
%             calculated temperatures but retain their individual P_kbar.
%

%% Input validation
% Basic argument checks to prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('EvansFrost1975 requires (rawdata_struct, P_kbar).');
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
% Extract the required tables from the input struct. We do not modify the
% tables here; we only read the relevant columns during calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end
if ~isfield(rawdata_struct, 'Spinel') || ~istable(rawdata_struct.Spinel)
    error('rawdata_struct must contain table: rawdata_struct.Spinel');
end

dataset_ol = rawdata_struct.Olivine;
dataset_sp = rawdata_struct.Spinel;

requiredVars_ol = {'Mg_cation_apfu', 'Fe_cation_apfu'};
requiredVars_sp = {'Mg_cation_apfu', 'Fe_cation_apfu', ...
    'Fe3_cation_apfu', 'Cr_cation_apfu', 'Al_cation_apfu'};

for i = 1:numel(requiredVars_ol)
    if ~ismember(requiredVars_ol{i}, dataset_ol.Properties.VariableNames)
        error('Olivine table must contain variable: %s', requiredVars_ol{i});
    end
end

for i = 1:numel(requiredVars_sp)
    if ~ismember(requiredVars_sp{i}, dataset_sp.Properties.VariableNames)
        error('Spinel table must contain variable: %s', requiredVars_sp{i});
    end
end

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

% Semi-quantitative graphical temperature interval represented by the two
% reference isotherms in Evans and Frost (1975), Fig. 8.
calibrationT_min_degC = 700;
calibrationT_max_degC = 1200;

% Evans and Frost (1975) did not define a numerical pressure calibration
% range, and pressure is not used in this implementation. Print one pressure
% applicability message after the first calculation.
pressureMessageIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% "Finish" after a calculation.
disp('=== Step 3: Selecting a data code from the list (Olivine) ===');

while true
    % ----- Olivine selection -----
    % Assumption: the first column stores an identifier to display to the user.
    dataCodes_ol = dataset_ol{:, 1};

    [selectedIdx_ol, ok] = listdlg( ...
        'PromptString', 'Please select the Olivine data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_ol, ...
        'ListSize', [320 320]);

    % If the user cancels, exit the loop gracefully.
    if ~ok || isempty(selectedIdx_ol)
        disp('Selection canceled');
        break;
    end

    selectedCode_ol = dataCodes_ol(selectedIdx_ol);
    disp(['Olivine selected: ' char(string(selectedCode_ol))]);

    % ----- Spinel selection -----
    disp('=== Step 4: Selecting a data code from the list (Spinel) ===');

    dataCodes_sp = dataset_sp{:, 1};

    [selectedIdx_sp, ok] = listdlg( ...
        'PromptString', 'Please select the Spinel data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_sp, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_sp)
        disp('Selection canceled');
        break;
    end

    selectedCode_sp = dataCodes_sp(selectedIdx_sp);
    disp(['Spinel selected: ' char(string(selectedCode_sp))]);

    % ----- Calculation -----
    % IMPORTANT:
    % Olivine and spinel are selected independently; do not assume row indices
    % correspond between the two datasets.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_sp = dataset_sp(selectedIdx_sp, :);

    % Check only the variables that are actually used in the thermometer.
    % The calculation is intentionally allowed to continue when NaN occurs;
    % a non-stopping fprintf message is printed after the result.
    nanInputNames = findNaNInputs(selectedData_ol, selectedData_sp);

    row = calcTemp(selectedData_ol, selectedData_sp, P_kbar);

    % Store the user-selected identifiers for traceability.
    row.dataCode_ol = repmat(string(selectedCode_ol), height(row), 1);
    row.dataCode_sp = repmat(string(selectedCode_sp), height(row), 1);

    % Move identifiers to the front for readability.
    row = movevars(row, {'dataCode_ol','dataCode_sp'}, 'Before', 1);

    % Store this result as one table block. If the preallocated buffer is
    % full, double its capacity. This avoids reallocating the full results
    % table after every calculation.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
            ': ' num2str(row.T_degC) ' degreeC']);
    else
        finiteDisplayValues = row.T_degC(isfinite(row.T_degC));
        if isempty(finiteDisplayValues)
            disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
                ': NaN to NaN degreeC']);
        else
            disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
                ': ' num2str(finiteDisplayValues(1)) ' to ' ...
                num2str(finiteDisplayValues(end)) ' degreeC']);
        end
    end

    % Evans and Frost (1975) do not report a numerical pressure calibration
    % interval. Print this applicability message once per function call.
    if ~pressureMessageIssued
        fprintf(2, ...
            ['WARNING: Evans and Frost (1975) did not define a numerical pressure ' ...
             'calibration range for this graphical thermometer. Input pressure ' ...
             'range = %.4g-%.4g kbar (%d point(s)).\n' ...
             '         Pressure is retained for interface compatibility but is not ' ...
             'used in the temperature calculation; pressure-related applicability ' ...
             'cannot be evaluated from the original calibration.\n'], ...
            min(P_kbar), ...
            max(P_kbar), ...
            numel(P_kbar));
        pressureMessageIssued = true;
    end

    % Warn when any finite calculated temperature lies outside the graphical
    % interval of 700-1200 degreeC. NaN and Inf are handled separately below.
    finiteTemperature = isfinite(row.T_degC);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_degC < calibrationT_min_degC | ...
         row.T_degC > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_degC(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the semi-quantitative ' ...
             'graphical calibration interval represented by Evans and Frost ' ...
             '(1975): 700-1200 degreeC. %d of %d finite temperature point(s) ' ...
             'are outside the interval; calculated finite range = %.4g-%.4g ' ...
             'degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)));
    end

    % Print a non-stopping message when required thermometer input contains
    % NaN. fprintf is used instead of warning so that the message remains in
    % the command-window log even when MATLAB warnings are disabled.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, but the calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.T_degC);

    % Retain a result-based check for NaN/Inf values caused by invalid ratios,
    % division by zero, or other compositional problems.
    if any(invalidTemperature)
        fprintf(2, ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
                    '(%d of %d points; NaN: %d, Inf: %d).\n' ...
                    '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            sum(invalidTemperature), ...
            numel(row.T_degC), ...
            sum(isnan(row.T_degC)), ...
            sum(isinf(row.T_degC)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'EvansFrost1975', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate the buffered table blocks only once after all selections have
% been completed. Return an empty table when no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_olivine, data_spinel)
% findNaNInputs
% Return the names of required thermometer input variables that contain NaN.
% This function does not throw an error for NaN values; it only prepares a
% warning message for the calling function.

olivineVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
spinelVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', ...
    'Fe3_cation_apfu', 'Cr_cation_apfu', 'Al_cation_apfu'};

nanInputNames = strings(0, 1);

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isnan(variableValue(:)))
        nanInputNames(end + 1, 1) = "Olivine." + string(variableName); %#ok<AGROW>
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isnan(variableValue(:)))
        nanInputNames(end + 1, 1) = "Spinel." + string(variableName); %#ok<AGROW>
    end
end

end

function row = calcTemp(data_olivine, data_spinel, P_kbar)
% calcTemp
% Compute an Evans and Frost (1975) semi-quantitative temperature estimate
% for one olivine row and one spinel row. One output row is returned for each
% supplied pressure value, although pressure is not used in the calculation.
%
% Inputs:
%   data_olivine : 1-row table containing olivine apfu cations
%   data_spinel  : 1-row table containing spinel apfu cations
%   P_kbar       : pressure in kbar; retained but not used
%
% Output:
%   row : table containing one row per input pressure value

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();

% Store pressure for traceability. PressureUsed is false for every row because
% the implemented graphical interpolation has no pressure term.
row.P_kbar = P_kbar;
row.PressureUsed = false(nP, 1);

% --- Extract cation data ---
% Repeat the selected mineral compositions for every supplied pressure point.
Mg_ol  = repmat(data_olivine.Mg_cation_apfu, nP, 1);
Fe2_ol = repmat(data_olivine.Fe_cation_apfu, nP, 1);

Mg_sp  = repmat(data_spinel.Mg_cation_apfu, nP, 1);
Fe2_sp = repmat(data_spinel.Fe_cation_apfu, nP, 1);
Fe3_sp = repmat(data_spinel.Fe3_cation_apfu, nP, 1);
Cr_sp  = repmat(data_spinel.Cr_cation_apfu, nP, 1);
Al_sp  = repmat(data_spinel.Al_cation_apfu, nP, 1);

% --- Composition terms ---
XMg_ol = Mg_ol ./ (Mg_ol + Fe2_ol);
XFe_ol = Fe2_ol ./ (Mg_ol + Fe2_ol);

XMg_sp  = Mg_sp ./ (Mg_sp + Fe2_sp);
XFe2_sp = Fe2_sp ./ (Mg_sp + Fe2_sp);

trivalentSum_sp = Cr_sp + Al_sp + Fe3_sp;
YCr_sp  = Cr_sp ./ trivalentSum_sp;
YFe3_sp = Fe3_sp ./ trivalentSum_sp;

% --- KD and normalized lnKD* ---
KD = (XMg_ol .* XFe2_sp) ./ (XFe_ol .* XMg_sp);
lnKD = log(KD);
lnKD_star = lnKD + (0.05 - YFe3_sp) .* 4.0;

% --- Evans and Frost (1975) graphical calibration ---
lnKD_700 = 0.85 + 2.94 .* YCr_sp;
lnKD_1200 = -0.05 + 2.94 .* YCr_sp;

T_degC = 700 + (1200 - 700) .* ...
    (lnKD_700 - lnKD_star) ./ (lnKD_700 - lnKD_1200);
T_K = T_degC + 273.15;

% --- Quality flags ---
% Start with OK, then overwrite individual rows where composition or the
% graphical interval indicates a problem.
QualityFlag = repmat("OK", nP, 1);

invalidComposition = ...
    ~isfinite(KD) | KD <= 0 | ...
    ~isfinite(XMg_ol) | ~isfinite(XFe_ol) | ...
    ~isfinite(XMg_sp) | ~isfinite(XFe2_sp) | ...
    ~isfinite(YCr_sp) | ~isfinite(YFe3_sp) | ...
    trivalentSum_sp <= 0;

belowCalibration = ~invalidComposition & lnKD_star > lnKD_700;
aboveCalibration = ~invalidComposition & lnKD_star < lnKD_1200;

QualityFlag(invalidComposition) = "Invalid_composition_or_missing_data";
QualityFlag(belowCalibration) = "Below_700C_or_outside_calibration";
QualityFlag(aboveCalibration) = "Above_1200C_or_outside_calibration";

T_degC(invalidComposition) = NaN;
T_K(invalidComposition) = NaN;

% --- Pack outputs ---
row.XMg_ol = XMg_ol;
row.XFe_ol = XFe_ol;

row.XMg_sp = XMg_sp;
row.XFe2_sp = XFe2_sp;

row.YCr_sp = YCr_sp;
row.YFe3_sp = YFe3_sp;

row.KD = KD;
row.lnKD = lnKD;
row.lnKD_star = lnKD_star;

row.lnKD_700 = lnKD_700;
row.lnKD_1200 = lnKD_1200;

row.T_K = T_K;
row.T_degC = T_degC;
row.QualityFlag = QualityFlag;

end
