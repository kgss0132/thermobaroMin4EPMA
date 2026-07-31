function results = Fujii1978(rawdata_struct, P_kbar)
% functions/+thermo/+Olivine_Spinel/Fujii1978.m
% Target environment: MATLAB R2024b
%
% Empirical Fe-Mg exchange thermometer between Olivine and Spinel
% Fujii, T. (1978)
% Fe-Mg partitioning between olivine and spinel.
% Carnegie Institution of Washington, Geophysical Laboratory,
% Year Book 76, 563-569.
% DOI: None
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis and one Spinel
% analysis and calculates temperature from Fe2+-Mg partitioning between
% the two minerals using the Fujii (1978) calibration.
%
% The user-selection and output workflow follows Ballhaus1991.m so that the
% function can be called by the same fixed-pressure and pressure-range
% launchers. Fujii's temperature equation contains no explicit pressure
% term. P_kbar is therefore retained for interface compatibility and output
% traceability, but changing pressure does not change the calculated
% temperature for a given Olivine-Spinel pair.
%
% -------------------------------------------------------------------------
% CALIBRATION AND APPLICATION NOTES
%
% Fujii (1978) investigated Olivine-Spinel Fe-Mg partitioning using
% experiments at 1 atm and 15 kbar, principally over 1200-1350 degreeC.
% Natural low-temperature Olivine-Spinel pairs interpreted at approximately
% 550 and 700 degreeC were also used to examine the low-temperature trend.
% The final relation was presented as an empirical geothermometer,
% particularly for spinel peridotites.
%
% Important cautions:
%   1) The equation has no explicit pressure correction. The pressure values
%      passed to this function are not used in the temperature solution.
%   2) Fe3+ in Spinel must be estimated separately from Fe2+. Treating total
%      Spinel Fe as Fe2+ changes both KD and the ferric correction.
%   3) Primary, equilibrated Olivine-Spinel pairs should be used. Magnetite,
%      ferritchromite, altered Spinel rims, or mismatched core-rim analyses
%      may not preserve the relevant equilibrium.
%   4) Fe-Mg exchange may continue during cooling, so calculated values can
%      record subsolidus re-equilibration rather than peak temperature.
%   5) Fujii (1978) discussed this relation as empirical; ideal-solution
%      behavior or cancellation of non-idealities should not be assumed for
%      every natural composition.
%
% This implementation issues a non-stopping warning when a finite result is
% outside 550-1350 degreeC. This interval represents the approximate range
% illustrated by the natural low-temperature and experimental high-
% temperature reference data; it is not treated as a strict universal
% applicability limit.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%   rawdata_struct.Spinel  : table
%
% The FIRST column of each table is treated as an identifier displayed in
% the selection dialog. The following normalized-cation variables are used:
%
%   Olivine:
%     Mg_cation_apfu
%     Fe_cation_apfu         % Fe2+ in Olivine
%
%   Spinel:
%     Mg_cation_apfu
%     Fe_cation_apfu         % Fe2+ in Spinel
%     Fe3_cation_apfu        % Fe3+ in Spinel
%     Cr_cation_apfu
%     Al_cation_apfu
%
% Finite Mg and Fe2+ values must be strictly positive because they occur in
% ratios and logarithms. Finite Fe3+, Cr, and Al values may be zero but must
% not be negative. The trivalent-cation sum Al + Cr + Fe3+ must be positive.
% NaN values are allowed to propagate and are reported by non-stopping
% warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% 1) Divalent-cation fractions
%   XMg_ol  = Mg_ol  / (Mg_ol  + Fe2_ol)
%   XFe2_ol = Fe2_ol / (Mg_ol  + Fe2_ol)
%   XMg_sp  = Mg_sp  / (Mg_sp  + Fe2_sp)
%   XFe2_sp = Fe2_sp / (Mg_sp  + Fe2_sp)
%
% 2) Trivalent-cation fractions in Spinel
%   YCr_sp  = Cr_sp  / (Al_sp + Cr_sp + Fe3_sp)
%   YAl_sp  = Al_sp  / (Al_sp + Cr_sp + Fe3_sp)
%   YFe3_sp = Fe3_sp / (Al_sp + Cr_sp + Fe3_sp)
%
% 3) Fe-Mg exchange coefficient
%   KD_ol_sp = (XMg_ol * XFe2_sp) / (XFe2_ol * XMg_sp)
%
% 4) Ferric-iron correction adopted by Fujii (1978)
%   lnK5 = 4
%   lnKD_star = ln(KD_ol_sp) - 4*YFe3_sp
%
% 5) Temperature equation (Fujii, 1978, Eq. 8)
%   numerator   = 775 + 2010*YCr_sp
%   denominator = lnKD_star - 0.006 - 0.003*YCr_sp
%
%   T(K) = numerator / denominator
%   T(degreeC) = T(K) - 273.15
%
% Equivalently:
%   T(K) = (775 + 2010*YCr_sp) / ...
%          (ln(KD_ol_sp) - 4*YFe3_sp - 0.006 - 0.003*YCr_sp)
%
% Pressure does not appear in this equation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Fujii1978(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Olivine and Spinel tables
%   P_kbar         : finite non-negative scalar or vector; retained for
%                    launcher compatibility and output only
%
% Output:
%   results : table containing one row per supplied pressure value for every
%             user-selected Olivine-Spinel pair
%

%% Input validation
if nargin < 2
    error('Fujii1978 requires (rawdata_struct, P_kbar).');
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
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end
if ~isfield(rawdata_struct, 'Spinel') || ~istable(rawdata_struct.Spinel)
    error('rawdata_struct must contain table: rawdata_struct.Spinel');
end

dataset_ol = rawdata_struct.Olivine;
dataset_sp = rawdata_struct.Spinel;

validateRequiredVariables(dataset_ol, dataset_sp);

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Approximate interval represented by Fujii's natural low-temperature
% reference pairs and high-temperature experiments.
referenceT_min_degC = 550;
referenceT_max_degC = 1350;

% The pressure-independence message is displayed only once per function call.
pressureNoteIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Olivine) ===');

while true
    % ----- Olivine selection -----
    dataCodes_ol = dataset_ol{:, 1};

    [selectedIdx_ol, ok] = listdlg( ...
        'PromptString', 'Please select the Olivine data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_ol, ...
        'ListSize', [320 320]);

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
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_sp = dataset_sp(selectedIdx_sp, :);

    nanInputNames = findNaNInputs(selectedData_ol, selectedData_sp);
    validateCompositionInputs(selectedData_ol, selectedData_sp);

    row = calcTemp(selectedData_ol, selectedData_sp, P_kbar);

    % Store the selected identifiers for traceability.
    row.dataCode_ol = repmat(string(selectedCode_ol), height(row), 1);
    row.dataCode_sp = repmat(string(selectedCode_sp), height(row), 1);
    row = movevars(row, {'dataCode_ol','dataCode_sp'}, 'Before', 1);

    % Store this result as one table block.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        finiteT = row.T_deg(isfinite(row.T_deg));
        if isempty(finiteT)
            temperatureText = 'non-finite result';
        elseif max(finiteT) == min(finiteT)
            temperatureText = [num2str(finiteT(1)) ' degreeC at all pressure points'];
        else
            temperatureText = [num2str(row.T_deg(1)) ' to ' ...
                num2str(row.T_deg(end)) ' degreeC'];
        end
        disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
            ': ' temperatureText]);
    end

    % Fujii's equation has no explicit pressure term. State this once so a
    % pressure-range calculation is not misinterpreted as pressure-sensitive.
    if ~pressureNoteIssued
        fprintf(2, ...
            ['NOTE: Fujii (1978) contains no explicit pressure term. ' ...
             'P_kbar is retained for interface compatibility and output only; ' ...
             'temperature is identical at all supplied pressure points for a given pair.\n']);
        pressureNoteIssued = true;
    end

    % Warn when finite temperatures fall outside the approximate interval
    % represented by the low-T natural and high-T experimental reference data.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideReference = finiteTemperature & ...
        (row.T_deg < referenceT_min_degC | row.T_deg > referenceT_max_degC);

    if any(temperatureOutsideReference)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature lies outside the approximate ' ...
             '550-1350 degreeC interval represented by the Fujii (1978) ' ...
             'natural and experimental reference data. %d of %d finite ' ...
             'temperature point(s) are outside; finite range = %.4g-%.4g ' ...
             'degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideReference), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, but the calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
                    '(%d of %d points; NaN: %d, Inf: %d).\n' ...
                    '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    invalidDenominator = isfinite(row.denominator) & (row.denominator <= 0);
    if any(invalidDenominator)
        fprintf(2, ...
            ['WARNING: The Fujii (1978) temperature denominator is zero or ' ...
             'negative for %s & %s. Such results are non-physical and should ' ...
             'not be interpreted.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Fujii1978', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all table blocks once after selection is complete.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(data_olivine, data_spinel)
% Validate that all variables required by Fujii (1978) are available.

olivineVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
spinelVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', ...
    'Fe3_cation_apfu', 'Cr_cation_apfu', 'Al_cation_apfu'};

missingNames = strings(0, 1);

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    if ~ismember(variableName, data_olivine.Properties.VariableNames)
        missingNames(end + 1, 1) = "Olivine." + string(variableName); %#ok<AGROW>
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    if ~ismember(variableName, data_spinel.Properties.VariableNames)
        missingNames(end + 1, 1) = "Spinel." + string(variableName); %#ok<AGROW>
    end
end

if ~isempty(missingNames)
    error(['Fujii1978: required table variable(s) are missing: ' ...
        char(strjoin(missingNames, ', ')) '.']);
end

end

function nanInputNames = findNaNInputs(data_olivine, data_spinel)
% Return names of required input variables containing NaN.

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

function validateCompositionInputs(data_olivine, data_spinel)
% Validate sign constraints required by the Fujii calculation.
%
% Mg and Fe2+ must be > 0 because they enter ratios and logarithms.
% Fe3+, Cr, and Al may be zero, but finite negative values are invalid.
% NaN is intentionally allowed to propagate to the output.

positiveOlivineVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
positiveSpinelVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
nonnegativeSpinelVariables = {'Fe3_cation_apfu', ...
    'Cr_cation_apfu', 'Al_cation_apfu'};

invalidPositiveNames = strings(0, 1);
invalidNonnegativeNames = strings(0, 1);

for i = 1:numel(positiveOlivineVariables)
    variableName = positiveOlivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) <= 0)
        invalidPositiveNames(end + 1, 1) = ...
            "Olivine." + string(variableName); %#ok<AGROW>
    end
end

for i = 1:numel(positiveSpinelVariables)
    variableName = positiveSpinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) <= 0)
        invalidPositiveNames(end + 1, 1) = ...
            "Spinel." + string(variableName); %#ok<AGROW>
    end
end

for i = 1:numel(nonnegativeSpinelVariables)
    variableName = nonnegativeSpinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        invalidNonnegativeNames(end + 1, 1) = ...
            "Spinel." + string(variableName); %#ok<AGROW>
    end
end

if ~isempty(invalidPositiveNames)
    error(['Fujii1978: Mg and Fe2+ inputs must be > 0. ' ...
        'Zero or negative finite value(s) were found in: ' ...
        char(strjoin(invalidPositiveNames, ', ')) '.']);
end

if ~isempty(invalidNonnegativeNames)
    error(['Fujii1978: Fe3+, Cr, and Al inputs must be >= 0. ' ...
        'Negative finite value(s) were found in: ' ...
        char(strjoin(invalidNonnegativeNames, ', ')) '.']);
end

trivalentSum = data_spinel.Al_cation_apfu + ...
    data_spinel.Cr_cation_apfu + data_spinel.Fe3_cation_apfu;

if any(isfinite(trivalentSum(:)) & trivalentSum(:) <= 0)
    error(['Fujii1978: the Spinel trivalent-cation sum ' ...
        '(Al + Cr + Fe3+) must be > 0.']);
end

end

function row = calcTemp(data_olivine, data_spinel, P_kbar)
% Calculate Fujii (1978) temperature for one Olivine-Spinel pair.
%
% P_kbar is retained to match the common thermometer interface. It does not
% enter the published Fujii temperature equation.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% --- Extract and expand selected cation data ---
Mg_ol  = repmat(data_olivine.Mg_cation_apfu, nP, 1);
Fe2_ol = repmat(data_olivine.Fe_cation_apfu, nP, 1);

Mg_sp  = repmat(data_spinel.Mg_cation_apfu, nP, 1);
Fe2_sp = repmat(data_spinel.Fe_cation_apfu, nP, 1);
Fe3_sp = repmat(data_spinel.Fe3_cation_apfu, nP, 1);
Cr_sp  = repmat(data_spinel.Cr_cation_apfu, nP, 1);
Al_sp  = repmat(data_spinel.Al_cation_apfu, nP, 1);

% --- Divalent-cation fractions ---
XMg_ol  = Mg_ol  ./ (Mg_ol + Fe2_ol);
XFe2_ol = Fe2_ol ./ (Mg_ol + Fe2_ol);

XMg_sp  = Mg_sp  ./ (Mg_sp + Fe2_sp);
XFe2_sp = Fe2_sp ./ (Mg_sp + Fe2_sp);

% --- Trivalent-cation fractions in Spinel ---
trivalentSum_sp = Al_sp + Cr_sp + Fe3_sp;
YCr_sp  = Cr_sp  ./ trivalentSum_sp;
YAl_sp  = Al_sp  ./ trivalentSum_sp;
YFe3_sp = Fe3_sp ./ trivalentSum_sp;

% --- Fe-Mg exchange coefficient ---
KD_ol_sp = (XMg_ol .* XFe2_sp) ./ (XFe2_ol .* XMg_sp);
lnKD = log(KD_ol_sp);

% --- Fujii ferric correction ---
lnK5 = 4.0;
lnKD_star = lnKD - lnK5 .* YFe3_sp;

% --- Fujii (1978) empirical temperature equation ---
numerator = 775 + 2010 .* YCr_sp;
denominator = lnKD_star - 0.006 - 0.003 .* YCr_sp;

T_K = numerator ./ denominator;
T_deg = T_K - 273.15;

% --- Pack outputs ---
row.XMg_ol = XMg_ol;
row.XFe2_ol = XFe2_ol;
row.XMg_sp = XMg_sp;
row.XFe2_sp = XFe2_sp;

row.YCr_sp = YCr_sp;
row.YAl_sp = YAl_sp;
row.YFe3_sp = YFe3_sp;
row.trivalentSum_sp = trivalentSum_sp;

row.KD_ol_sp = KD_ol_sp;
row.lnKD = lnKD;
row.lnK5 = repmat(lnK5, nP, 1);
row.lnKD_star = lnKD_star;

row.numerator = numerator;
row.denominator = denominator;
row.T_K = T_K;
row.T_deg = T_deg;

end
