function results = Harley1984(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/Harley1984.m
% Tested with MATLAB R2024b
%
% Garnet-orthopyroxene Fe2+-Mg exchange thermometer
% Harley, S.L. (1984)
% Contributions to Mineralogy and Petrology, 86, 359-373
% DOI: https://doi.org/10.1007/BF01187140
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Orthopyroxene analysis (selected by the user from tables) and calculates
% temperature using the Harley (1984) garnet-orthopyroxene Fe2+-Mg
% exchange thermometer.
%
% Pressure may be supplied as either a scalar or a vector. Therefore, this
% function can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. For each selected Grt-Opx pair, the output table
% contains one row per pressure value.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Grt-Opx pair, and appends results into
% a single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Harley (1984) experimentally investigated Fe-Mg partitioning between
% garnet and aluminous orthopyroxene in the following synthetic systems
% (abstract, p. 359; Introduction, p. 360):
%
%   FMAS experiments (Ca-free base calibration):
%     Temperature : 800-1200 degreeC
%     Pressure    : 5-30 kbar
%     Bulk XMg    : Mg/(Mg+Fe) = 0.15-0.80
%
%   CFMAS experiments (Ca-bearing calibration used for Eq. 11):
%     Temperature : 900-1200 degreeC
%     Pressure    : 7.5-30 kbar
%     Bulk XMg    : Mg/(Mg+Fe) = 0.50, with variable Ca contents
%
% Equation (11), which includes the garnet Ca correction implemented here,
% is presented on p. 368. Because this implementation uses Eq. (11), the
% narrower CFMAS interval of 900-1200 degreeC and 7.5-30 kbar is used for
% non-stopping calibration-range warnings.
%
% Important application notes from Harley (1984):
%   1) The thermometer is intended principally for garnet peridotites and
%      granulites (abstract, p. 359).
%   2) The symmetric regular-solution treatment of Ca-Mg-Fe garnet was
%      developed for 0 < Xgr_Gt < 0.20 (p. 361). Application to garnet with
%      substantially higher grossular fraction is an extrapolation.
%   3) Within the experimental uncertainty, orthopyroxene was treated as
%      macroscopically ideal (pp. 359-361). The Appendix shows that the
%      neglected non-ideal contribution is within the analytical and
%      experimental uncertainty over the principal experimental region,
%      approximately XMg_Opx = 0.20-0.80 and XAl_M1 = 0.05-0.15
%      (pp. 371-372).
%   4) Equation (10) reproduces the adjusted FMAS run temperatures to about
%      +/-60 degreeC (p. 367). Most CFMAS temperatures calculated with
%      Eq. (11) lie within about +/-40 degreeC of nominal run temperatures,
%      but the Ca interaction term is 1400 +/- 500 cal/mol site and carries
%      substantial uncertainty (pp. 367-368). These experimental residuals
%      are not a universal uncertainty for every natural sample.
%   5) Precision decreases at lower temperature. An analytical error of
%      only +/-0.1 in KD can produce approximately +/-40-60 degreeC error
%      in calculated temperature (p. 368). Small Mg analytical errors are
%      especially important for Fe-rich garnet (p. 370).
%   6) Mn, Cr3+, and Fe3+ were not included in the experiments. Mn enters
%      garnet preferentially; ignoring it gives an erroneously low
%      temperature. Harley suggests either combining Mn with Ca as a first
%      approximation or not applying the thermometer to Mn-rich garnet.
%      Most granulite and peridotite garnets have XMn_Gt < 0.04, whereas
%      metamorphosed iron formations may be substantially richer in Mn
%      (p. 369). The optional Ca+Mn approximation is not automatically
%      applied here; this implementation reproduces Eq. (11) using Ca.
%   7) KD is an Fe2+-Mg exchange coefficient. If all natural-rock Fe is
%      treated as Fe2+ because Fe3+ has not been estimated, the calculated
%      value should be regarded as a minimum temperature. A consistent
%      Fe3+ estimate can increase temperature by approximately 20-100
%      degreeC (p. 369). If Fe_cation_apfu and Fe3_cation_apfu are already
%      separated, Fe3+ must not be added back to Fe2+; this implementation
%      therefore uses Fe_cation_apfu alone in KD and Xgr_Gt.
%   8) Natural Grt-Opx temperatures were commonly 50-130 degreeC below
%      Grt-Cpx temperatures above approximately 1000 degreeC. Possible
%      causes include calibration errors, different Fe-Mg re-equilibration
%      rates, and failure to correct for Fe3+ (pp. 369-370).
%   9) Core-rim zoning, corona growth, recrystallisation, and continued
%      Fe-Mg exchange during cooling can cause the calculated temperature
%      to record a later event rather than peak metamorphism. Garnet and
%      orthopyroxene zones must be paired petrographically, and their
%      chemical zoning should be examined before interpreting a temperature
%      (pp. 369-370).
%  10) The calibration assumes a known pressure. The pressure term is not
%      negligible and P_kbar must represent the pressure appropriate to the
%      selected Grt-Opx equilibration stage.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside 7.5-30 kbar,
%   2) a finite calculated temperature is outside 900-1200 degreeC,
%   3) a required thermometer input contains NaN, or
%   4) a calculated temperature is NaN or Inf. In this case, the actual
%      thermometer inputs and identified calculation-domain causes are also
%      displayed.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Opx    : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns contain
% normalized cation data.
%
% Variables used directly by the thermometer:
%   Garnet table:
%     Fe_cation_apfu         % Fe2+ in garnet
%     Mg_cation_apfu
%     Ca_cation_apfu
%
%   Orthopyroxene table:
%     Fe_cation_apfu         % Fe2+ in orthopyroxene
%     Mg_cation_apfu
%
% Optional variables retained in the output when available:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%   Ca_cation_apfu in Opx
%
% IMPORTANT Fe note:
% The exchange reaction is defined using Fe2+. Therefore,
% Fe3_cation_apfu is retained for reference but is not added to
% Fe_cation_apfu. If Fe_cation_apfu contains total Fe rather than Fe2+, Fe2+
% should be estimated before this function is used.
%
% Negative finite values in variables used by the thermometer are not
% permitted. Zero values are retained; if they make the equation
% mathematically undefined, the affected result is returned as NaN and a
% non-stopping warning is printed. NaN values are never replaced by zero;
% they propagate through the calculation and remain in the output.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Fe2+-Mg distribution coefficient:
%
%   KD = (Fe2+/Mg)_Grt / (Fe2+/Mg)_Opx
%
% Mineral indices:
%
%   XMg_Gt  = Mg / (Mg + Fe2+)
%   XMg_Opx = Mg / (Mg + Fe2+)
%   Xgr_Gt  = Ca / (Ca + Mg + Fe2+)
%
% Harley (1984), Eq. (11), p. 368:
%
%          3740 + 1400*Xgr_Gt + 22.86*P_kbar
%   T(K) = ---------------------------------------
%                    R*ln(KD) + 1.96
%
% where R = 1.9872 cal/mol/K. Temperature is converted using
% T(degreeC) = T(K) - 273.15.
%
% The previous implementation added Fe3_cation_apfu to Fe_cation_apfu and
% calculated XMg_Gt on a Ca+Mg+Fe basis. Those behaviors are not used here:
% Fe3+ is excluded from KD, and XMg_Gt follows Harley's Mg/(Mg+Fe) binary
% definition while Xgr_Gt uses Ca/(Ca+Mg+Fe).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Harley1984(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Opx tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Grt-Opx pair
%

%% Input validation
if nargin < 2
    error('Harley1984 requires (rawdata_struct, P_kbar).');
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

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end

dataset_grt = rawdata_struct.Garnet;
dataset_opx = rawdata_struct.Opx;

validateRequiredVariables(dataset_grt, dataset_opx);

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Result blocks are buffered and concatenated once after the interactive
% loop, avoiding repeated reallocation of the complete output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct CFMAS calibration limits for the Ca-bearing Eq. (11).
calibrationT_min_degC = 900;
calibrationT_max_degC = 1200;
calibrationP_min_kbar = 7.5;
calibrationP_max_kbar = 30;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
    dataCodes_grt = dataset_grt{:, 1};

    [selectedIdx_grt, ok] = listdlg( ...
        'PromptString', 'Please select the Garnet data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_grt)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_grt)
        disp('Selection canceled');
        break;
    end

    selectedCode_grt = dataCodes_grt(selectedIdx_grt);
    disp(['Garnet selected: ' char(string(selectedCode_grt))]);

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

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_opx = dataset_opx(selectedIdx_opx, :);

    validateNonNegativeInputs(selectedData_grt, selectedData_opx);
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_opx);

    row = calcTemp(selectedData_grt, selectedData_opx, P_kbar);

    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_opx'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo calculated temperature values.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_opx)) ...
            ': ' num2str(row.T_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_opx)) ...
            ': ' num2str(row.T_C(1)) ' to ' num2str(row.T_C(end)) ' degreeC']);
    end

    % Pressure warning is printed once because the same pressure vector is
    % used for every mineral pair in this function call.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the direct experimental ' ...
             'calibration range of the Ca-bearing Harley (1984) Eq. (11): ' ...
             '7.5-30 kbar. %d of %d pressure point(s) are outside the ' ...
             'range; input range = %.4g-%.4g kbar (pp. 359-360).\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn for finite temperatures outside the CFMAS experimental interval.
    finiteTemperature = isfinite(row.T_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.T_C < calibrationT_min_degC | ...
         row.T_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.T_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the direct ' ...
             'experimental calibration range of the Ca-bearing Harley ' ...
             '(1984) Eq. (11): 900-1200 degreeC. %d of %d finite ' ...
             'temperature point(s) are outside the range; calculated ' ...
             'finite range = %.4g-%.4g degreeC for %s & %s ' ...
             '(pp. 359-360, 368).\n'], ...
            sum(temperatureOutsideCalibration), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)));
    end

    % Report explicitly stored NaN thermometer inputs. Calculation continues
    % and NaN values are not replaced by zero.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, and NaN was not replaced by zero.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Preserve non-finite results and print both raw inputs and identified
    % calculation-domain causes.
    invalidTemperature = ~isfinite(row.T_C);
    if any(invalidTemperature)
        nonFiniteCauses = findNonFiniteCauses(row);

        fprintf(2, ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
                    '(%d of %d points; NaN: %d, Inf: %d).\n' ...
                    '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)), ...
            sum(invalidTemperature), ...
            numel(row.T_C), ...
            sum(isnan(row.T_C)), ...
            sum(isinf(row.T_C)));

        fprintf(2, ['         Thermometer inputs used: ' ...
                    'Garnet.Fe_cation_apfu=%s, Garnet.Mg_cation_apfu=%s, ' ...
                    'Garnet.Ca_cation_apfu=%s, ' ...
                    'Opx.Fe_cation_apfu=%s, Opx.Mg_cation_apfu=%s.\n'], ...
            formatNumericValue(row.Fe2_grt(1)), ...
            formatNumericValue(row.Mg_grt(1)), ...
            formatNumericValue(row.Ca_grt(1)), ...
            formatNumericValue(row.Fe2_opx(1)), ...
            formatNumericValue(row.Mg_opx(1)));

        if isempty(nonFiniteCauses)
            fprintf(2, ['         No explicit NaN, zero, Inf, or invalid intermediate value ' ...
                        'was identified; inspect the stored intermediate variables.\n']);
        else
            fprintf(2, '         Identified cause(s): %s.\n', ...
                char(strjoin(nonFiniteCauses, '; ')));
        end
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Harley1984', ...
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

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(dataset_grt, dataset_opx)
% validateRequiredVariables
% Verify all columns required by the Harley (1984) Eq. (11).

requiredGrt = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu'};
requiredOpx = {'Fe_cation_apfu', 'Mg_cation_apfu'};

missingNames = strings(numel(requiredGrt) + numel(requiredOpx), 1);
nMissing = 0;

for i = 1:numel(requiredGrt)
    if ~ismember(requiredGrt{i}, dataset_grt.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Garnet." + string(requiredGrt{i});
    end
end

for i = 1:numel(requiredOpx)
    if ~ismember(requiredOpx{i}, dataset_opx.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Opx." + string(requiredOpx{i});
    end
end

if nMissing > 0
    missingNames = missingNames(1:nMissing);
    error(['Harley1984: required variable(s) are missing: ' ...
        char(strjoin(missingNames, ', ')) '.']);
end

end


function nanInputNames = findNaNInputs(data_grt, data_opx)
% findNaNInputs
% Return names of required thermometer inputs containing NaN.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu'};
opxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

nanInputNames = strings(numel(grtVariables) + numel(opxVariables), 1);
nNaN = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    variableValue = data_grt.(variableName);
    if any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanInputNames(nNaN) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(opxVariables)
    variableName = opxVariables{i};
    variableValue = data_opx.(variableName);
    if any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanInputNames(nNaN) = "Opx." + string(variableName);
    end
end

nanInputNames = nanInputNames(1:nNaN);

end


function validateNonNegativeInputs(data_grt, data_opx)
% validateNonNegativeInputs
% Stop when a required thermometer input is negative or infinite. Zero and
% NaN remain available to the calculation and result diagnostics.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Ca_cation_apfu'};
opxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

invalidInputNames = strings(numel(grtVariables) + numel(opxVariables), 1);
nInvalid = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    variableValue = data_grt.(variableName);
    if ~isnumeric(variableValue) || ~isscalar(variableValue)
        error('Garnet.%s must be a numeric scalar in the selected row.', ...
            variableName);
    end
    if isinf(variableValue) || variableValue < 0
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(opxVariables)
    variableName = opxVariables{i};
    variableValue = data_opx.(variableName);
    if ~isnumeric(variableValue) || ~isscalar(variableValue)
        error('Opx.%s must be a numeric scalar in the selected row.', ...
            variableName);
    end
    if isinf(variableValue) || variableValue < 0
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Opx." + string(variableName);
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNames(1:nInvalid);
    error(['Harley1984: thermometer inputs must be finite or NaN and ' ...
        'must be >= 0. Negative or infinite value(s) were found in: ' ...
        char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function row = calcTemp(data_grt, data_opx, P_kbar)
% calcTemp
% Compute Harley (1984) Grt-Opx temperatures for one mineral pair and a
% scalar or vector of pressures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;

% Gas constant used by Harley (1984), in cal/mol/K.
R_cal = 1.9872;
row.R_cal = repmat(R_cal, nP, 1);

% --- Extract and expand garnet cations ---
Fe2_grt = repmat(getRequiredValue(data_grt, 'Fe_cation_apfu', 'Garnet'), nP, 1);
Fe3_grt = repmat(getOptionalValue(data_grt, 'Fe3_cation_apfu', 0, 'Garnet'), nP, 1);
Mg_grt  = repmat(getRequiredValue(data_grt, 'Mg_cation_apfu', 'Garnet'), nP, 1);
Mn_grt  = repmat(getOptionalValue(data_grt, 'Mn_cation_apfu', 0, 'Garnet'), nP, 1);
Ca_grt  = repmat(getRequiredValue(data_grt, 'Ca_cation_apfu', 'Garnet'), nP, 1);
Ti_grt  = repmat(getOptionalValue(data_grt, 'Ti_cation_apfu', 0, 'Garnet'), nP, 1);
Al_grt  = repmat(getOptionalValue(data_grt, 'Al_cation_apfu', 0, 'Garnet'), nP, 1);
Si_grt  = repmat(getOptionalValue(data_grt, 'Si_cation_apfu', 0, 'Garnet'), nP, 1);
Na_grt  = repmat(getOptionalValue(data_grt, 'Na_cation_apfu', 0, 'Garnet'), nP, 1);
K_grt   = repmat(getOptionalValue(data_grt, 'K_cation_apfu', 0, 'Garnet'), nP, 1);

% --- Extract and expand orthopyroxene cations ---
Fe2_opx = repmat(getRequiredValue(data_opx, 'Fe_cation_apfu', 'Opx'), nP, 1);
Fe3_opx = repmat(getOptionalValue(data_opx, 'Fe3_cation_apfu', 0, 'Opx'), nP, 1);
Mg_opx  = repmat(getRequiredValue(data_opx, 'Mg_cation_apfu', 'Opx'), nP, 1);
Mn_opx  = repmat(getOptionalValue(data_opx, 'Mn_cation_apfu', 0, 'Opx'), nP, 1);
Ca_opx  = repmat(getOptionalValue(data_opx, 'Ca_cation_apfu', 0, 'Opx'), nP, 1);
Ti_opx  = repmat(getOptionalValue(data_opx, 'Ti_cation_apfu', 0, 'Opx'), nP, 1);
Al_opx  = repmat(getOptionalValue(data_opx, 'Al_cation_apfu', 0, 'Opx'), nP, 1);
Si_opx  = repmat(getOptionalValue(data_opx, 'Si_cation_apfu', 0, 'Opx'), nP, 1);
Na_opx  = repmat(getOptionalValue(data_opx, 'Na_cation_apfu', 0, 'Opx'), nP, 1);
K_opx   = repmat(getOptionalValue(data_opx, 'K_cation_apfu', 0, 'Opx'), nP, 1);

% --- Fe2+-Mg mineral indices and exchange coefficient ---
% Fe3+ is deliberately excluded from every term used by Eq. (11).
sum_grt_FeMg = Fe2_grt + Mg_grt;
sum_grt_FeMgCa = Fe2_grt + Mg_grt + Ca_grt;
sum_opx_FeMg = Fe2_opx + Mg_opx;

XMg_grt = Mg_grt ./ sum_grt_FeMg;
Xgr_grt = Ca_grt ./ sum_grt_FeMgCa;
XMg_opx = Mg_opx ./ sum_opx_FeMg;

FeMg_grt = Fe2_grt ./ Mg_grt;
FeMg_opx = Fe2_opx ./ Mg_opx;
KD = FeMg_grt ./ FeMg_opx;
lnKD = log(KD);

% --- Harley (1984), Eq. (11), p. 368 ---
temperatureNumerator = 3740 + 1400 .* Xgr_grt + 22.86 .* P_kbar;
temperatureDenominator = R_cal .* lnKD + 1.96;

% Retain NaN values and return NaN for mathematically invalid domains.
calculationDomainValid = isfinite(KD) & KD > 0 ...
    & isfinite(sum_grt_FeMg) & sum_grt_FeMg > 0 ...
    & isfinite(sum_grt_FeMgCa) & sum_grt_FeMgCa > 0 ...
    & isfinite(sum_opx_FeMg) & sum_opx_FeMg > 0 ...
    & isfinite(XMg_grt) & XMg_grt >= 0 & XMg_grt <= 1 ...
    & isfinite(Xgr_grt) & Xgr_grt >= 0 & Xgr_grt <= 1 ...
    & isfinite(XMg_opx) & XMg_opx >= 0 & XMg_opx <= 1 ...
    & isfinite(temperatureNumerator) ...
    & isfinite(temperatureDenominator) & temperatureDenominator > 0;

T_K = NaN(nP, 1);
T_K(calculationDomainValid) = ...
    temperatureNumerator(calculationDomainValid) ./ ...
    temperatureDenominator(calculationDomainValid);
T_C = T_K - 273.15;

% --- Pack outputs ---
% FeUsed is retained for compatibility with the original output variable
% set, but correctly equals Fe2+ rather than Fe2+ + Fe3+.
row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.FeUsed_grt = Fe2_grt;
row.Mg_grt = Mg_grt;
row.Mn_grt = Mn_grt;
row.Ca_grt = Ca_grt;
row.Ti_grt = Ti_grt;
row.Al_grt = Al_grt;
row.Si_grt = Si_grt;
row.Na_grt = Na_grt;
row.K_grt = K_grt;

row.Fe2_opx = Fe2_opx;
row.Fe3_opx = Fe3_opx;
row.FeUsed_opx = Fe2_opx;
row.Mg_opx = Mg_opx;
row.Mn_opx = Mn_opx;
row.Ca_opx = Ca_opx;
row.Ti_opx = Ti_opx;
row.Al_opx = Al_opx;
row.Si_opx = Si_opx;
row.Na_opx = Na_opx;
row.K_opx = K_opx;

row.sum_grt_FeMg = sum_grt_FeMg;
row.sum_grt_FeMgCa = sum_grt_FeMgCa;
row.sum_opx_FeMg = sum_opx_FeMg;
row.XMg_grt = XMg_grt;
row.Xgr_grt = Xgr_grt;
row.XMg_opx = XMg_opx;
row.FeMg_grt = FeMg_grt;
row.FeMg_opx = FeMg_opx;
row.KD = KD;
row.lnKD = lnKD;
row.temperatureNumerator = temperatureNumerator;
row.temperatureDenominator = temperatureDenominator;
row.denominator = temperatureDenominator;
row.calculationDomainValid = calculationDomainValid;
row.T_K = T_K;
row.T_C = T_C;

end


function nonFiniteCauses = findNonFiniteCauses(row)
% findNonFiniteCauses
% Identify raw-input and derived-variable conditions that can produce a NaN
% or Inf temperature. These diagnostics do not alter stored values.

maximumCauses = 24;
nonFiniteCauses = strings(maximumCauses, 1);
nCauses = 0;

inputValues = {row.Fe2_grt, row.Mg_grt, row.Ca_grt, ...
    row.Fe2_opx, row.Mg_opx};
inputNames = {"Garnet.Fe_cation_apfu", "Garnet.Mg_cation_apfu", ...
    "Garnet.Ca_cation_apfu", "Opx.Fe_cation_apfu", ...
    "Opx.Mg_cation_apfu"};
zeroIsInvalid = [true, true, false, true, true];

for i = 1:numel(inputValues)
    value = inputValues{i};
    if any(isnan(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is NaN";
    elseif any(isinf(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is Inf";
    elseif zeroIsInvalid(i) && any(value(:) == 0)
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is zero";
    end
end

derivedValues = {row.sum_grt_FeMg, row.sum_grt_FeMgCa, ...
    row.sum_opx_FeMg, row.KD, row.lnKD, row.XMg_grt, row.Xgr_grt, ...
    row.XMg_opx, row.temperatureNumerator, row.temperatureDenominator};
derivedNames = {"garnet Fe-Mg sum", "garnet Fe-Mg-Ca sum", ...
    "orthopyroxene Fe-Mg sum", "KD", "lnKD", "XMg_grt", ...
    "Xgr_grt", "XMg_opx", "temperature numerator", ...
    "temperature denominator"};

for i = 1:numel(derivedValues)
    value = derivedValues{i};
    if any(~isfinite(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is non-finite";
    elseif any(value(:) <= 0) && ...
            (derivedNames{i} == "garnet Fe-Mg sum" || ...
             derivedNames{i} == "garnet Fe-Mg-Ca sum" || ...
             derivedNames{i} == "orthopyroxene Fe-Mg sum" || ...
             derivedNames{i} == "KD" || ...
             derivedNames{i} == "temperature denominator")
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is zero or negative";
    elseif any(value(:) < 0) && ...
            (derivedNames{i} == "XMg_grt" || ...
             derivedNames{i} == "Xgr_grt" || ...
             derivedNames{i} == "XMg_opx")
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is negative";
    end
end

if any(~row.calculationDomainValid)
    nCauses = nCauses + 1;
    nonFiniteCauses(nCauses) = "Harley (1984) calculation domain is invalid";
end

nonFiniteCauses = nonFiniteCauses(1:nCauses);

end


function textValue = formatNumericValue(value)
% formatNumericValue
% Format a scalar numeric value for a compact diagnostic message.

if isnan(value)
    textValue = 'NaN';
elseif isinf(value) && value > 0
    textValue = 'Inf';
elseif isinf(value) && value < 0
    textValue = '-Inf';
else
    textValue = sprintf('%.8g', value);
end

end


function value = getRequiredValue(data_tbl, variableName, mineralLabel)
% getRequiredValue
% Extract a required scalar numeric value. NaN is returned unchanged.

if ~ismember(variableName, data_tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = data_tbl.(variableName);
if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end


function value = getOptionalValue(data_tbl, variableName, defaultValue, mineralLabel)
% getOptionalValue
% Use defaultValue only when an optional column is absent. An explicitly
% stored NaN is returned unchanged and is never converted to zero.

if ismember(variableName, data_tbl.Properties.VariableNames)
    value = data_tbl.(variableName);
    if ~isnumeric(value) || ~isscalar(value)
        error('%s.%s must be a numeric scalar in the selected row.', ...
            mineralLabel, variableName);
    end
else
    value = defaultValue;
end

end
