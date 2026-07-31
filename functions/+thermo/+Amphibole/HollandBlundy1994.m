function results = HollandBlundy1994(rawdata_struct, P_kbar)
% functions/+thermo/+Amphibole/HollandBlundy1994.m
% Tested with MATLAB R2024b
%
% Amphibole-Plagioclase thermometers
% Holland, T.J.B. and Blundy, J.D. (1994)
% Contributions to Mineralogy and Petrology, 116, 433–447
% DOI: https://doi.org/10.1007/BF00310910
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Amphibole analysis and one
% Plagioclase analysis and calculates temperatures using the two
% Holland and Blundy (1994) thermometers:
%
%   Thermometer A : edenite-tremolite thermometer
%                   (quartz-bearing / silica-saturated assemblages only)
%
%   Thermometer B : edenite-richterite thermometer
%                   (assemblages with or without quartz)
%
% The function accepts either a scalar pressure or a pressure vector. It is
% therefore compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for each pressure value
% for every user-selected Amphibole-Plagioclase pair.
%
% This implementation retains the practical site-allocation scheme used in
% the supplied HollandBlundy1994.m script:
%   - T site total = 8 apfu
%   - B(M4) site total = 2 apfu
%   - A site total = 1 apfu
%
% The calculated site fractions are then used in the published thermometer
% equations. For strict consistency with the original calibration, users
% should also consider the amphibole Fe3+ recalculation and site-allocation
% procedure presented in Appendix B of Holland and Blundy (1994), p. 445.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Holland and Blundy (1994) report that the two thermometers perform with a
% typical uncertainty of approximately +/-35–40 degreeC over a broad range
% of natural and experimental amphibole-plagioclase compositions. The
% abstract gives an overall range of 400–1000 degreeC and 1–15 kbar
% (p. 433), but the thermometer-specific restrictions in the main text are
% narrower and should be used for application screening:
%
% Thermometer A (edenite-tremolite; p. 438):
%   Temperature              : 400–900 degreeC
%   Pressure                 : 1–15 kbar (overall calibration range, p. 433)
%   Silica condition         : quartz-bearing / silica-saturated only
%   Amphibole Na on A site   : X_Na_A > 0.02
%   Amphibole octahedral Al  : Al(VI) < 1.8 apfu
%   Amphibole Si             : 6.0–7.7 apfu
%   Plagioclase              : X_An < 0.90
%
% Thermometer B (edenite-richterite; pp. 439–440):
%   Temperature              : 500–900 degreeC
%   Pressure                 : 1–15 kbar (overall calibration range, p. 433)
%   Silica condition         : may be used with or without quartz
%   Amphibole Na on M4 site  : X_Na_M4 > 0.03
%   Amphibole octahedral Al  : Al(VI) < 1.8 apfu
%   Amphibole Si             : 6.0–7.7 apfu
%   Plagioclase              : 0.10 < X_An < 0.90
%
% Additional cautions from the original paper:
% - Thermometer A should not be used for silica-undersaturated rocks or
%   magmas because it may return anomalously high temperatures (p. 441).
% - Thermometer B is poorly constrained for nearly pure albite and may have
%   little value in low-temperature rocks containing very albitic
%   plagioclase (pp. 439–440).
% - Application to kaersutite and Ti-rich richterite is discouraged because
%   Ti-bearing substitutions are not represented adequately by the adopted
%   amphibole model (p. 434).
% - Amphibole and plagioclase must represent an equilibrium pair. Incorrect
%   Fe3+/Fe2+ estimation, disequilibrium pairing, or compositions outside
%   the calibration dataset may produce unreliable temperatures (p. 441).
% - The two thermometer equations should not be intersected to estimate
%   pressure, and their difference should not be used quantitatively to
%   calculate silica activity (p. 441).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 1–15 kbar,
%   2) finite Thermometer A temperatures are outside 400–900 degreeC,
%   3) finite Thermometer B temperatures are outside 500–900 degreeC,
%   4) a calculation input contains NaN, or
%   5) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole   : table
%   rawdata_struct.Plagioclase : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following variable names.
%
% Required Amphibole variables (normalized cations; preferably 23 O basis):
%   Si_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu          % total Fe
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%
% Optional Amphibole variables:
%   K_cation_apfu
%   Mn_cation_apfu
%   Cr_cation_apfu
%   Fe3_cation_apfu
%
% Required Plagioclase variables (apfu cations):
%   Si_cation_apfu
%   Al_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%
% Optional Plagioclase variable:
%   K_cation_apfu
%
% Finite mineral-composition values must be non-negative. Negative values
% and Inf values stop the calculation. NaN values are retained as missing
% values, propagated through calculations that depend on them, and reported
% by non-stopping fprintf warnings. An optional variable that is completely
% absent from a table is assigned its stated default value of zero; an
% existing optional variable containing NaN is never converted to zero.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Thermometer A:
%   T_A(K) =
%     (-76.95 + 0.79*P + Y_ab + 39.4*X_Na_A + 22.4*X_K_A ...
%      + (41.5 - 2.89*P)*X_Al_M2) ...
%     / (-0.0650 - R*ln((27*X_vA*X_Si_T1*X_ab_plag) ...
%                       /(256*X_Na_A*X_Al_T1)))
%
%   Y_ab = 0                                  for X_ab > 0.5
%        = 12.0*(1 - X_ab)^2 - 3.0           otherwise
%
% Thermometer B:
%   T_B(K) =
%     (78.44 + Y_ab_an - 33.6*X_Na_M4 ...
%      - (66.8 - 2.92*P)*X_Al_M2 ...
%      + 78.5*X_Al_T1 + 9.4*X_Na_A) ...
%     / (0.0721 - R*ln((27*X_Na_M4*X_Si_T1*X_an_plag) ...
%                      /(64*X_Ca_M4*X_Al_T1*X_ab_plag)))
%
%   Y_ab_an = 3.0                             for X_ab > 0.5
%           = 12.0*(2*X_ab - 1) + 3.0        otherwise
%
% Units:
%   T in K
%   P in kbar
%   R = 0.0083144 kJ K^-1 mol^-1
%
% -------------------------------------------------------------------------
% Syntax:
%   results = HollandBlundy1994(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Amphibole and Plagioclase tables
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Amphibole-Plagioclase pair
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values. Pressure may be either a scalar or a vector.
if nargin < 2
    error('HollandBlundy1994 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || ~isreal(P_kbar) || isempty(P_kbar) || ...
        ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The tables are not
% modified here; relevant variables are read after the user selects a pair.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Amphibole') || ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end
if ~isfield(rawdata_struct, 'Plagioclase') || ~istable(rawdata_struct.Plagioclase)
    error('rawdata_struct must contain table: rawdata_struct.Plagioclase');
end

dataset_amp = rawdata_struct.Amphibole;
dataset_plag = rawdata_struct.Plagioclase;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated table concatenation inside the interactive loop is avoided because
% it repeatedly reallocates and copies the entire output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Published pressure and thermometer-specific temperature limits.
calibrationP_min_kbar = 1;
calibrationP_max_kbar = 15;
calibrationTA_min_degC = 400;
calibrationTA_max_degC = 900;
calibrationTB_min_degC = 500;
calibrationTB_max_degC = 900;

% Pressure is common to all selected mineral pairs in this function call.
% The pressure warning is therefore printed only once, after the first
% calculation has been displayed.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % ----- Amphibole selection -----
    % Assumption: the first column stores an identifier displayed to the user.
    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, ok] = listdlg( ...
        'PromptString', 'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_amp)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

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

    % ----- Calculation -----
    % Amphibole and plagioclase are selected independently; row indices are
    % not assumed to correspond between the two datasets.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);
    selectedData_plag = dataset_plag(selectedIdx_plag, :);

    % Reject negative or infinite values while intentionally permitting NaN.
    % Existing NaN values are retained and allowed to propagate.
    validateNonNegativeInputs(selectedData_amp, selectedData_plag);

    % Record input names containing NaN so that a non-stopping warning can be
    % printed immediately after the calculated temperature.
    nanInputNames = findNaNInputs(selectedData_amp, selectedData_plag);

    row = calcTemp(selectedData_amp, selectedData_plag, P_kbar);

    % Store the selected identifiers for every pressure row.
    row.dataCode_amphibole = ...
        repmat(string(selectedCode_amp), height(row), 1);
    row.dataCode_plagioclase = ...
        repmat(string(selectedCode_plag), height(row), 1);
    row = movevars(row, ...
        {'dataCode_amphibole', 'dataCode_plagioclase'}, 'Before', 1);

    % Store this calculation as one table block. The buffer is enlarged only
    % when its current capacity has been exhausted, rather than on every loop.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo calculated temperatures for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_amp)) ' & ' ...
            char(string(selectedCode_plag)) ...
            ': T_A = ' num2str(row.T_A_deg) ...
            ' degreeC, T_B = ' num2str(row.T_B_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_amp)) ' & ' ...
            char(string(selectedCode_plag)) ...
            ': T_A = ' num2str(row.T_A_deg(1)) ' to ' ...
            num2str(row.T_A_deg(end)) ...
            ' degreeC, T_B = ' num2str(row.T_B_deg(1)) ' to ' ...
            num2str(row.T_B_deg(end)) ' degreeC']);
    end

    % Warn once when any input pressure lies outside 1–15 kbar.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the published calibration ' ...
             'range of Holland and Blundy (1994): 1–15 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g–%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when finite Thermometer A temperatures lie outside 400–900 degreeC.
    finiteTemperatureA = isfinite(row.T_A_deg);
    temperatureAOutsideCalibration = finiteTemperatureA & ...
        (row.T_A_deg < calibrationTA_min_degC | ...
         row.T_A_deg > calibrationTA_max_degC);

    if any(temperatureAOutsideCalibration)
        finiteValuesA = row.T_A_deg(finiteTemperatureA);
        fprintf(2, ...
            ['WARNING: Thermometer A temperature is outside the published ' ...
             'application range of Holland and Blundy (1994): ' ...
             '400–900 degreeC (p. 438). ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureAOutsideCalibration), ...
            sum(finiteTemperatureA), ...
            min(finiteValuesA), ...
            max(finiteValuesA), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    % Warn when finite Thermometer B temperatures lie outside 500–900 degreeC.
    finiteTemperatureB = isfinite(row.T_B_deg);
    temperatureBOutsideCalibration = finiteTemperatureB & ...
        (row.T_B_deg < calibrationTB_min_degC | ...
         row.T_B_deg > calibrationTB_max_degC);

    if any(temperatureBOutsideCalibration)
        finiteValuesB = row.T_B_deg(finiteTemperatureB);
        fprintf(2, ...
            ['WARNING: Thermometer B temperature is outside the published ' ...
             'application range of Holland and Blundy (1994): ' ...
             '500–900 degreeC (p. 440). ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'calculated finite range = %.4g–%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureBOutsideCalibration), ...
            sum(finiteTemperatureB), ...
            min(finiteValuesB), ...
            max(finiteValuesB), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    % Report NaN inputs without interrupting the calculation.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and the calculation was continued; ' ...
             'dependent output values may also be NaN.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report non-finite Thermometer A results.
    invalidTemperatureA = ~isfinite(row.T_A_deg);
    if any(invalidTemperatureA)
        fprintf(2, ...
            ['WARNING: Non-finite Thermometer A values were calculated for ' ...
             '%s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            sum(invalidTemperatureA), ...
            numel(row.T_A_deg), ...
            sum(isnan(row.T_A_deg)), ...
            sum(isinf(row.T_A_deg)));
    end

    % Retain and report non-finite Thermometer B results.
    invalidTemperatureB = ~isfinite(row.T_B_deg);
    if any(invalidTemperatureB)
        fprintf(2, ...
            ['WARNING: Non-finite Thermometer B values were calculated for ' ...
             '%s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            sum(invalidTemperatureB), ...
            numel(row.T_B_deg), ...
            sum(isnan(row.T_B_deg)), ...
            sum(isinf(row.T_B_deg)));
    end

    disp('--------------------------------------------------');

    % Ask whether another mineral pair should be calculated.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'HollandBlundy1994', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after all selections have been
% completed. Return an empty table when no calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_amphibole, data_plagioclase)
% findNaNInputs
% Return names of supplied composition variables containing NaN. Missing
% optional variables are not reported because they are assigned default zero.

amphiboleRequired = {'Si_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu'};
amphiboleOptional = {'K_cation_apfu', 'Mn_cation_apfu', ...
    'Cr_cation_apfu', 'Fe3_cation_apfu'};
plagioclaseRequired = {'Si_cation_apfu', 'Al_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu'};
plagioclaseOptional = {'K_cation_apfu'};

maximumNames = numel(amphiboleRequired) + numel(amphiboleOptional) + ...
    numel(plagioclaseRequired) + numel(plagioclaseOptional);
nanInputBuffer = strings(maximumNames, 1);
nNanInputs = 0;

for i = 1:numel(amphiboleRequired)
    variableName = amphiboleRequired{i};
    variableValue = data_amphibole.(variableName);
    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = ...
            "Amphibole." + string(variableName);
    end
end

for i = 1:numel(amphiboleOptional)
    variableName = amphiboleOptional{i};
    if ismember(variableName, data_amphibole.Properties.VariableNames)
        variableValue = data_amphibole.(variableName);
        if any(isnan(variableValue(:)))
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = ...
                "Amphibole." + string(variableName);
        end
    end
end

for i = 1:numel(plagioclaseRequired)
    variableName = plagioclaseRequired{i};
    variableValue = data_plagioclase.(variableName);
    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = ...
            "Plagioclase." + string(variableName);
    end
end

for i = 1:numel(plagioclaseOptional)
    variableName = plagioclaseOptional{i};
    if ismember(variableName, data_plagioclase.Properties.VariableNames)
        variableValue = data_plagioclase.(variableName);
        if any(isnan(variableValue(:)))
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = ...
                "Plagioclase." + string(variableName);
        end
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_amphibole, data_plagioclase)
% validateNonNegativeInputs
% Require all supplied composition values to be numeric scalars that are
% either non-negative finite values or NaN. NaN is deliberately accepted so
% that it propagates through the calculation and triggers fprintf warnings.

amphiboleRequired = {'Si_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu'};
amphiboleOptional = {'K_cation_apfu', 'Mn_cation_apfu', ...
    'Cr_cation_apfu', 'Fe3_cation_apfu'};
plagioclaseRequired = {'Si_cation_apfu', 'Al_cation_apfu', ...
    'Ca_cation_apfu', 'Na_cation_apfu'};
plagioclaseOptional = {'K_cation_apfu'};

maximumNames = numel(amphiboleRequired) + numel(amphiboleOptional) + ...
    numel(plagioclaseRequired) + numel(plagioclaseOptional);
invalidInputBuffer = strings(maximumNames, 1);
nInvalidInputs = 0;

for i = 1:numel(amphiboleRequired)
    variableName = amphiboleRequired{i};
    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        error('Amphibole table must contain variable: %s', variableName);
    end

    variableValue = data_amphibole.(variableName);
    if ~isnumeric(variableValue) || ~isreal(variableValue) || ~isscalar(variableValue)
        error('Amphibole variable %s must be a numeric scalar.', variableName);
    end
    if any(isinf(variableValue(:)) | variableValue(:) < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Amphibole." + string(variableName);
    end
end

for i = 1:numel(amphiboleOptional)
    variableName = amphiboleOptional{i};
    if ismember(variableName, data_amphibole.Properties.VariableNames)
        variableValue = data_amphibole.(variableName);
        if ~isnumeric(variableValue) || ~isreal(variableValue) || ~isscalar(variableValue)
            error('Amphibole variable %s must be a numeric scalar.', variableName);
        end
        if any(isinf(variableValue(:)) | variableValue(:) < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = ...
                "Amphibole." + string(variableName);
        end
    end
end

for i = 1:numel(plagioclaseRequired)
    variableName = plagioclaseRequired{i};
    if ~ismember(variableName, data_plagioclase.Properties.VariableNames)
        error('Plagioclase table must contain variable: %s', variableName);
    end

    variableValue = data_plagioclase.(variableName);
    if ~isnumeric(variableValue) || ~isreal(variableValue) || ~isscalar(variableValue)
        error('Plagioclase variable %s must be a numeric scalar.', variableName);
    end
    if any(isinf(variableValue(:)) | variableValue(:) < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Plagioclase." + string(variableName);
    end
end

for i = 1:numel(plagioclaseOptional)
    variableName = plagioclaseOptional{i};
    if ismember(variableName, data_plagioclase.Properties.VariableNames)
        variableValue = data_plagioclase.(variableName);
        if ~isnumeric(variableValue) || ~isreal(variableValue) || ~isscalar(variableValue)
            error('Plagioclase variable %s must be a numeric scalar.', variableName);
        end
        if any(isinf(variableValue(:)) | variableValue(:) < 0)
            nInvalidInputs = nInvalidInputs + 1;
            invalidInputBuffer(nInvalidInputs) = ...
                "Plagioclase." + string(variableName);
        end
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['HollandBlundy1994: mineral-composition values must be ' ...
           'non-negative finite values or NaN. Negative or Inf value(s) ' ...
           'were found in: ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_amphibole, data_plagioclase, P_kbar)
% calcTemp
% Compute both Holland and Blundy (1994) temperatures for one selected
% Amphibole-Plagioclase pair over one or more pressure values.

P_kbar = P_kbar(:);
nP = numel(P_kbar);
R_kJ = 0.0083144;

% --- Prepare selected mineral rows ---
amphibole = prepareAmphiboleRow(data_amphibole);
plagioclase = preparePlagioclaseRow(data_plagioclase);

% --- Site allocation ---
site_amphibole = calcAmphiboleSites(amphibole);
site_plagioclase = calcPlagioclaseSites(plagioclase);

% --- Holland and Blundy variables from the site allocation ---
X_Na_A = site_amphibole.Na_A;
X_K_A = site_amphibole.K_A;
X_vA = site_amphibole.vA;

X_Ca_M4 = site_amphibole.Ca_B ./ 2;
X_Na_M4 = site_amphibole.Na_B ./ 2;

X_Al_T1 = site_amphibole.Al_4 ./ 4;
X_Si_T1 = 1 - X_Al_T1;
X_Al_M2 = site_amphibole.Al_6 ./ 2;

X_ab_plag = site_plagioclase.Xab;
X_an_plag = site_plagioclase.Xan;

% --- Plagioclase non-ideality terms ---
if isnan(X_ab_plag)
    Y_ab = NaN;
    Y_ab_an = NaN;
elseif X_ab_plag > 0.5
    Y_ab = 0.0;
    Y_ab_an = 3.0;
else
    Y_ab = 12.0 .* (1 - X_ab_plag).^2 - 3.0;
    Y_ab_an = 12.0 .* (2 .* X_ab_plag - 1) + 3.0;
end

% --- Thermometer A ---
% Numerator is pressure-dependent and is evaluated for all pressure values.
num_A = -76.95 + 0.79 .* P_kbar + Y_ab + 39.4 .* X_Na_A + ...
    22.4 .* X_K_A + (41.5 - 2.89 .* P_kbar) .* X_Al_M2;

lnArg_A = NaN(nP, 1);
denom_A = NaN(nP, 1);
T_A_K = NaN(nP, 1);

if X_vA > 0 && X_Si_T1 > 0 && X_ab_plag > 0 && ...
        X_Na_A > 0 && X_Al_T1 > 0
    lnArgAValue = (27 .* X_vA .* X_Si_T1 .* X_ab_plag) ./ ...
        (256 .* X_Na_A .* X_Al_T1);

    if isfinite(lnArgAValue) && lnArgAValue > 0
        lnArg_A(:) = lnArgAValue;
        denomAValue = -0.0650 - R_kJ .* log(lnArgAValue);
        denom_A(:) = denomAValue;

        if isfinite(denomAValue) && abs(denomAValue) > 1e-12
            T_A_K = num_A ./ denomAValue;
        end
    end
end

% --- Thermometer B ---
num_B = 78.44 + Y_ab_an - 33.6 .* X_Na_M4 - ...
    (66.8 - 2.92 .* P_kbar) .* X_Al_M2 + ...
    78.5 .* X_Al_T1 + 9.4 .* X_Na_A;

lnArg_B = NaN(nP, 1);
denom_B = NaN(nP, 1);
T_B_K = NaN(nP, 1);

if X_Na_M4 > 0 && X_Si_T1 > 0 && X_an_plag > 0 && ...
        X_Ca_M4 > 0 && X_Al_T1 > 0 && X_ab_plag > 0
    lnArgBValue = ...
        (27 .* X_Na_M4 .* X_Si_T1 .* X_an_plag) ./ ...
        (64 .* X_Ca_M4 .* X_Al_T1 .* X_ab_plag);

    if isfinite(lnArgBValue) && lnArgBValue > 0
        lnArg_B(:) = lnArgBValue;
        denomBValue = 0.0721 - R_kJ .* log(lnArgBValue);
        denom_B(:) = denomBValue;

        if isfinite(denomBValue) && abs(denomBValue) > 1e-12
            T_B_K = num_B ./ denomBValue;
        end
    end
end

T_A_deg = T_A_K - 273.15;
T_B_deg = T_B_K - 273.15;

% --- Composition-based applicability flags ---
% These flags reproduce the composition restrictions encoded in the supplied
% script. Temperature and pressure screening is reported separately by
% fprintf warnings in the main function.
isApplicable_A = ...
    (X_Na_A > 0.02) && ...
    (site_amphibole.Al_6 < 1.8) && ...
    (amphibole.Si >= 6.0) && (amphibole.Si <= 7.7) && ...
    (X_an_plag < 0.90);

isApplicable_B = ...
    (X_an_plag > 0.1) && (X_an_plag < 0.9) && ...
    (X_Na_M4 > 0.03) && ...
    (site_amphibole.Al_6 < 1.8) && ...
    (amphibole.Si >= 6.0) && (amphibole.Si <= 7.7);

% --- Pack outputs ---
% Scalar composition and site-allocation values are replicated to match the
% number of pressure rows. This permits both fixed-pressure and pressure-range
% workflows without changing the output schema.
row = table();
row.P_kbar = P_kbar;

row.Si_amp = repmat(amphibole.Si, nP, 1);
row.Ti_amp = repmat(amphibole.Ti, nP, 1);
row.Al_amp = repmat(amphibole.Al, nP, 1);
row.Fe2_amp = repmat(amphibole.Fe2, nP, 1);
row.Fe3_amp = repmat(amphibole.Fe3, nP, 1);
row.Mg_amp = repmat(amphibole.Mg, nP, 1);
row.Ca_amp = repmat(amphibole.Ca, nP, 1);
row.Na_amp = repmat(amphibole.Na, nP, 1);
row.K_amp = repmat(amphibole.K, nP, 1);
row.Mn_amp = repmat(amphibole.Mn, nP, 1);
row.Cr_amp = repmat(amphibole.Cr, nP, 1);

row.Al_4_amp = repmat(site_amphibole.Al_4, nP, 1);
row.Al_6_amp = repmat(site_amphibole.Al_6, nP, 1);
row.Ca_B_amp = repmat(site_amphibole.Ca_B, nP, 1);
row.Na_B_amp = repmat(site_amphibole.Na_B, nP, 1);
row.Na_A_amp = repmat(site_amphibole.Na_A, nP, 1);
row.K_A_amp = repmat(site_amphibole.K_A, nP, 1);
row.vA_amp = repmat(site_amphibole.vA, nP, 1);

row.Si_plag = repmat(plagioclase.Si, nP, 1);
row.Al_plag = repmat(plagioclase.Al, nP, 1);
row.Ca_plag = repmat(plagioclase.Ca, nP, 1);
row.Na_plag = repmat(plagioclase.Na, nP, 1);
row.K_plag = repmat(plagioclase.K, nP, 1);

row.Xab_plag = repmat(X_ab_plag, nP, 1);
row.Xan_plag = repmat(X_an_plag, nP, 1);
row.Xor_plag = repmat(site_plagioclase.Xor, nP, 1);

row.Y_ab_kJ = repmat(Y_ab, nP, 1);
row.Y_ab_an_kJ = repmat(Y_ab_an, nP, 1);

row.X_Na_A = repmat(X_Na_A, nP, 1);
row.X_K_A = repmat(X_K_A, nP, 1);
row.X_vA = repmat(X_vA, nP, 1);
row.X_Ca_M4 = repmat(X_Ca_M4, nP, 1);
row.X_Na_M4 = repmat(X_Na_M4, nP, 1);
row.X_Al_T1 = repmat(X_Al_T1, nP, 1);
row.X_Si_T1 = repmat(X_Si_T1, nP, 1);
row.X_Al_M2 = repmat(X_Al_M2, nP, 1);

row.lnArg_A = lnArg_A;
row.lnArg_B = lnArg_B;
row.num_A = num_A;
row.denom_A = denom_A;
row.num_B = num_B;
row.denom_B = denom_B;

row.T_A_K = T_A_K;
row.T_A_deg = T_A_deg;
row.T_B_K = T_B_K;
row.T_B_deg = T_B_deg;

row.isApplicable_A = repmat(isApplicable_A, nP, 1);
row.isApplicable_B = repmat(isApplicable_B, nP, 1);

end

function amphibole = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Extract one amphibole analysis. Missing optional variables default to zero,
% whereas NaN in an existing variable is retained unchanged.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amphibole = struct();
amphibole.Si = getRequiredVariable( ...
    data_amphibole, 'Si_cation_apfu', 'Amphibole');
amphibole.Ti = getRequiredVariable( ...
    data_amphibole, 'Ti_cation_apfu', 'Amphibole');
amphibole.Al = getRequiredVariable( ...
    data_amphibole, 'Al_cation_apfu', 'Amphibole');
amphibole.FeT = getRequiredVariable( ...
    data_amphibole, 'Fe_cation_apfu', 'Amphibole');
amphibole.Mg = getRequiredVariable( ...
    data_amphibole, 'Mg_cation_apfu', 'Amphibole');
amphibole.Ca = getRequiredVariable( ...
    data_amphibole, 'Ca_cation_apfu', 'Amphibole');
amphibole.Na = getRequiredVariable( ...
    data_amphibole, 'Na_cation_apfu', 'Amphibole');

amphibole.K = getOptionalVariable( ...
    data_amphibole, 'K_cation_apfu', 0);
amphibole.Mn = getOptionalVariable( ...
    data_amphibole, 'Mn_cation_apfu', 0);
amphibole.Cr = getOptionalVariable( ...
    data_amphibole, 'Cr_cation_apfu', 0);
amphibole.Fe3 = getOptionalVariable( ...
    data_amphibole, 'Fe3_cation_apfu', 0);

if isfinite(amphibole.Fe3) && isfinite(amphibole.FeT) && ...
        amphibole.Fe3 > amphibole.FeT
    error('Amphibole contains Fe3_cation_apfu > Fe_cation_apfu.');
end

amphibole.Fe2 = amphibole.FeT - amphibole.Fe3;

end

function plagioclase = preparePlagioclaseRow(data_plagioclase)
% preparePlagioclaseRow
% Extract one plagioclase analysis. Missing K defaults to zero, whereas NaN
% in an existing K variable is retained unchanged.

if height(data_plagioclase) ~= 1
    error('Plagioclase input must be a 1-row table.');
end

plagioclase = struct();
plagioclase.Si = getRequiredVariable( ...
    data_plagioclase, 'Si_cation_apfu', 'Plagioclase');
plagioclase.Al = getRequiredVariable( ...
    data_plagioclase, 'Al_cation_apfu', 'Plagioclase');
plagioclase.Ca = getRequiredVariable( ...
    data_plagioclase, 'Ca_cation_apfu', 'Plagioclase');
plagioclase.Na = getRequiredVariable( ...
    data_plagioclase, 'Na_cation_apfu', 'Plagioclase');
plagioclase.K = getOptionalVariable( ...
    data_plagioclase, 'K_cation_apfu', 0);

end

function site = calcAmphiboleSites(amphibole)
% calcAmphiboleSites
% Apply the practical T-, B(M4)-, and A-site allocation retained from the
% supplied script. Explicit NaN branches prevent missing values from being
% replaced by zero during min/max operations.

site = struct();

% --- T site total = 8 apfu ---
if isnan(amphibole.Si) || isnan(amphibole.Al)
    site.Al_4 = NaN;
    site.Al_6 = NaN;
else
    T_deficit = max(0, 8 - amphibole.Si);
    site.Al_4 = min(amphibole.Al, T_deficit);
    site.Al_6 = amphibole.Al - site.Al_4;

    if site.Al_6 < -1e-10
        error('Negative octahedral Al calculated. Check cation normalization.');
    end
    site.Al_6 = max(0, site.Al_6);
end

% --- B(M4) site total = 2 apfu ---
if isnan(amphibole.Ca)
    site.Ca_B = NaN;
else
    site.Ca_B = min(amphibole.Ca, 2);
end

if isnan(amphibole.Na) || isnan(site.Ca_B)
    site.Na_B = NaN;
else
    site.Na_B = min(amphibole.Na, max(0, 2 - site.Ca_B));
end

% --- A site total = 1 apfu ---
if isnan(amphibole.Na) || isnan(site.Na_B)
    site.Na_A = NaN;
else
    site.Na_A = max(0, amphibole.Na - site.Na_B);
end

site.K_A = amphibole.K;
site.A_occ = site.Na_A + site.K_A;

if isfinite(site.A_occ) && site.A_occ > 1 + 1e-8
    error('Calculated A-site occupancy exceeds 1. Check cation normalization.');
end

if isnan(site.A_occ)
    site.vA = NaN;
else
    site.vA = max(0, 1 - site.A_occ);
end

end

function site = calcPlagioclaseSites(plagioclase)
% calcPlagioclaseSites
% Calculate albite, anorthite, and orthoclase fractions. A zero denominator
% is allowed to produce NaN so that the non-finite warning system handles the
% mathematically undefined result without replacing it by zero.

site = struct();
alkaliEarthSum = plagioclase.Ca + plagioclase.Na + plagioclase.K;

site.Xab = plagioclase.Na ./ alkaliEarthSum;
site.Xan = plagioclase.Ca ./ alkaliEarthSum;
site.Xor = plagioclase.K ./ alkaliEarthSum;

end

function value = getRequiredVariable(tableRow, variableName, mineralLabel)
% getRequiredVariable
% Return one required numeric scalar from a one-row table.

if ~ismember(variableName, tableRow.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tableRow.(variableName);

if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
    error('%s variable %s must be a numeric scalar.', ...
        mineralLabel, variableName);
end

end

function value = getOptionalVariable(tableRow, variableName, defaultValue)
% getOptionalVariable
% Return an optional variable when present. Use defaultValue only when the
% variable is absent; an existing NaN value is returned unchanged.

if ismember(variableName, tableRow.Properties.VariableNames)
    value = tableRow.(variableName);

    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value)
        error('Variable %s must be a numeric scalar.', variableName);
    end
else
    value = defaultValue;
end

end
