function results = HynesForest1988(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/HynesForest1988.m
% Tested with MATLAB R2024b
%
% Empirical garnet-muscovite Fe-Mg exchange thermometer
% Hynes, A. and Forest, R.C. (1988)
% Journal of Metamorphic Geology, 6, 297-309
% DOI: https://doi.org/10.1111/j.1525-1314.1988.tb00422.x
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% (muscovite) analysis and calculates temperature using the preferred
% empirical calibration of Hynes & Forest (1988), their equation (9).
%
% The function accepts pressure as either a scalar or a vector. Therefore,
% it can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. Equation (9) itself has no pressure term, so the
% same calculated temperature is returned for every pressure value. One
% output row is nevertheless produced for every pressure value to preserve
% the common thermometer interface and the input pressure history.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Hynes & Forest (1988) did not derive equation (9) from a dedicated
% experimental P-T calibration. It was empirically fitted to natural,
% regionally metamorphosed metapelites containing coexisting garnet,
% biotite, and muscovite. Reference temperatures were calculated with the
% Ferry & Spear (1978) garnet-biotite thermometer (abstract, p. 297;
% calibration dataset, pp. 299-303).
%
%   Rock type    : regionally metamorphosed metapelites
%   Mineral pair : coexisting garnet and phengitic muscovite
%   Pressure     : 3-7 kbar for the natural calibration dataset (p. 299)
%   Temperature  : no single formal numerical calibration range stated
%   Direct use   : 370-460 degreeC in the Ptarmigan Creek application
%                  (p. 308)
%   Main target  : low-grade metapelites, especially where garnet occurs
%                  below the biotite isograd (p. 297)
%
% The 370-460 degreeC interval is the directly demonstrated Ptarmigan Creek
% application range, not a formally stated full calibration interval. The
% natural regression dataset also included higher-grade garnet-biotite-
% muscovite metapelites, but the paper does not give one numerical minimum
% and maximum temperature for that complete dataset.
%
% Important application cautions stated in the original paper:
% - Equation (9) assumes no pressure dependence. Pressure and ferric-iron
%   effects could not be separated reliably, so the authors did not favour
%   the tentative pressure-dependent equation (10) (pp. 303-307).
% - All Fe was treated as ferrous. Treating actual Fe3+ in muscovite as
%   Fe2+ can raise the calculated temperature. The calibration rocks were
%   generally ilmenite-bearing and magnetite-free; magnetite- or hematite-
%   bearing suites may give erroneous results (pp. 305-307).
% - Fe and Mg concentrations in regional-metamorphic muscovite can be low,
%   so analytical bias and imprecision may overwhelm the temperature signal
%   (pp. 298-299).
% - Low-grade garnets commonly contain substantial Mn and Ca. The preferred
%   equation uses the Ganguly & Saxena (1984) garnet model modified for a
%   non-linear Ca effect after Hoinkes (1986), but those mixing models were
%   derived under conditions different from this application (pp. 301-303).
% - For zoned garnets, the calibration used higher-grade compositions; the
%   authors' own calculations used garnet rims only (Table 1 notes,
%   pp. 303-306). Garnet and muscovite should represent the same equilibrium
%   stage.
% - Differences from garnet-biotite temperatures were commonly about
%   15 degreeC, but this is only a minimum estimate of absolute error. The
%   authors state that absolute errors may be as high as 100 degreeC and
%   that calculated temperatures may be too low (p. 308).
%
% This implementation therefore issues non-stopping fprintf messages when:
%   1) input pressure is outside 3-7 kbar,
%   2) a finite temperature is outside the directly demonstrated
%      370-460 degreeC Ptarmigan Creek application range,
%   3) a required calculation input is NaN, or
%   4) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Mica   : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The following normalized cation
% variables are required because they enter equation (9) directly:
%
%   Garnet table variables:
%     Fe_cation_apfu
%     Mg_cation_apfu
%     Mn_cation_apfu
%     Ca_cation_apfu
%
%   Mica table variables:
%     Fe_cation_apfu
%     Mg_cation_apfu
%
% Optional calculation variable:
%   Fe3_cation_apfu
%
% In this implementation, Fe_used is calculated as:
%   Fe_used = Fe_cation_apfu + Fe3_cation_apfu
%
% when Fe3_cation_apfu is present. This reproduces the paper's simplifying
% treatment of all analysed Fe as ferrous only if Fe_cation_apfu stores Fe2+
% and Fe3_cation_apfu stores a separate, non-overlapping ferric component.
% If Fe_cation_apfu already represents total Fe, Fe3_cation_apfu must not be
% added again. Users must confirm the convention of their input tables.
%
% Optional output-only variables retained when present:
%   Ti_cation_apfu, Al_cation_apfu, Si_cation_apfu,
%   K_cation_apfu, Na_cation_apfu
%
% All finite mineral-composition values used in the calculation must be
% greater than or equal to zero. Negative finite values and Inf stop the
% calculation with an error. Zero is allowed, although it may generate a
% non-finite mathematical result that is retained and reported. NaN is
% never replaced by zero: it is propagated and reported by non-stopping
% fprintf messages. Missing Mn or Ca columns are not interpreted as zero.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Garnet mole fractions:
%   X_i_Grt = i / (Fe_used + Mg + Mn + Ca)
%
% Muscovite Fe-Mg fractions:
%   X_Fe_Mus = Fe_used / (Fe_used + Mg)
%   X_Mg_Mus = Mg      / (Fe_used + Mg)
%
% Exchange coefficient:
%   Kd = (X_Mg_Mus * X_Fe_Grt) / (X_Fe_Mus * X_Mg_Grt)
%
% Garnet Fe-Mg interaction parameter:
%   W_FeMg = 200  * X_Mg_Grt / (X_Mg_Grt + X_Fe_Grt) ...
%          + 2500 * X_Fe_Grt / (X_Mg_Grt + X_Fe_Grt)
%
% Corrected lnK:
%   lnK = lnKd ...
%       + (0.8*W_FeMg - W_FeMg*(X_Fe_Grt-X_Mg_Grt) ...
%          - 3000*X_Mn_Grt) / (R*T) ...
%       - 2.978*X_Ca_Grt*(844/T) ...
%       + 5.906*(X_Ca_Grt*(844/T))^2
%
% Preferred equation (9):
%   T(K) = 4790 / (lnK + 4.13)
%
% Because lnK contains T-dependent terms, equation (9) is solved
% numerically. R = 1.987 cal mol^-1 K^-1, consistent with the interaction
% parameters used in this formulation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = HynesForest1988(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Mica tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Mica pair
%

%% Input validation
if nargin < 2
    error('HynesForest1988 requires (rawdata_struct, P_kbar).');
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
if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_grt = rawdata_struct.Garnet;
dataset_mica = rawdata_struct.Mica;

requiredGarnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
requiredMicaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

validateRequiredVariables(dataset_grt, requiredGarnetVariables, 'Garnet');
validateRequiredVariables(dataset_mica, requiredMicaVariables, 'Mica');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each calculation as one table block and concatenate the blocks only
% once after the interactive loop. This avoids changing and copying the
% complete results table on every iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationP_min_kbar = 3;
calibrationP_max_kbar = 7;
demonstratedT_min_degC = 370;
demonstratedT_max_degC = 460;

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

    % ----- Mica selection -----
    disp('=== Step 4: Selecting a data code from the list (Mica) ===');

    dataCodes_mica = dataset_mica{:, 1};

    [selectedIdx_mica, ok] = listdlg( ...
        'PromptString', 'Please select the Mica data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_mica)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_mica)
        disp('Selection canceled');
        break;
    end

    selectedCode_mica = dataCodes_mica(selectedIdx_mica);
    disp(['Mica selected: ' char(string(selectedCode_mica))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_mica = dataset_mica(selectedIdx_mica, :);

    nanInputNames = findNaNInputs(selectedData_grt, selectedData_mica);
    validateNonNegativeInputs(selectedData_grt, selectedData_mica);

    row = calcTemp(selectedData_grt, selectedData_mica, P_kbar);

    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_mica = repmat(string(selectedCode_mica), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_mica'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': Hynes-Forest 1988 Eq. (9) = ' ...
            num2str(row.TEq9_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': Hynes-Forest 1988 Eq. (9) = ' ...
            num2str(row.TEq9_C(1)) ' to ' ...
            num2str(row.TEq9_C(end)) ' degreeC']);
    end

    % Warn once when pressure lies outside the natural calibration dataset.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the empirical pressure range ' ...
             'represented by the Hynes & Forest (1988) natural calibration ' ...
             'dataset: 3-7 kbar (p. 299). Equation (9) has no pressure term. ' ...
             '%d of %d pressure point(s) are outside the range; input range = ' ...
             '%.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), numel(P_kbar), ...
            min(P_kbar), max(P_kbar));
        pressureWarningIssued = true;
    end

    % The paper states no formal numerical T calibration range. This warning
    % therefore refers specifically to the directly demonstrated low-grade
    % Ptarmigan Creek application range of 370-460 degreeC.
    finiteTemperature = isfinite(row.TEq9_C);
    temperatureOutsideDemonstratedRange = finiteTemperature & ...
        (row.TEq9_C < demonstratedT_min_degC | ...
         row.TEq9_C > demonstratedT_max_degC);

    if any(temperatureOutsideDemonstratedRange)
        finiteValues = row.TEq9_C(finiteTemperature);
        fprintf(2, ...
            ['CAUTION: Calculated temperature is outside the directly demonstrated ' ...
             'Ptarmigan Creek application range of Hynes & Forest (1988): ' ...
             '370-460 degreeC (p. 308). This is not a formally stated full ' ...
             'calibration range. %d of %d finite temperature point(s) are outside ' ...
             'this interval; calculated finite range = %.4g-%.4g degreeC for ' ...
             '%s & %s.\n'], ...
            sum(temperatureOutsideDemonstratedRange), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)));
    end

    % Warn without stopping when a calculation input contains NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; it was not replaced by zero.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report non-finite results.
    invalidTemperature = ~isfinite(row.TEq9_C);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)), ...
            sum(invalidTemperature), numel(row.TEq9_C), ...
            sum(isnan(row.TEq9_C)), sum(isinf(row.TEq9_C)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'HynesForest1988', ...
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
function validateRequiredVariables(dataset, requiredVariables, mineralLabel)
% validateRequiredVariables
% Confirm that required calculation columns exist and contain numeric data.

missingVariables = requiredVariables( ...
    ~ismember(requiredVariables, dataset.Properties.VariableNames));

if ~isempty(missingVariables)
    error('%s table is missing required variable(s): %s.', ...
        mineralLabel, char(strjoin(string(missingVariables), ', ')));
end

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};
    if ~isnumeric(dataset.(variableName))
        error('%s.%s must be numeric.', mineralLabel, variableName);
    end
end

if ismember('Fe3_cation_apfu', dataset.Properties.VariableNames) && ...
        ~isnumeric(dataset.Fe3_cation_apfu)
    error('%s.Fe3_cation_apfu must be numeric when present.', mineralLabel);
end

end

function nanInputNames = findNaNInputs(data_garnet, data_mica)
% findNaNInputs
% Return calculation input names containing NaN. The name buffer is
% preallocated and NaN is never changed to zero.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu', 'Fe3_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Fe3_cation_apfu'};

nameBuffer = strings(numel(garnetVariables) + numel(micaVariables), 1);
nNames = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if any(isnan(variableValue(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
    if ismember(variableName, data_mica.Properties.VariableNames)
        variableValue = data_mica.(variableName);
        if any(isnan(variableValue(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Mica." + string(variableName);
        end
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_garnet, data_mica)
% validateNonNegativeInputs
% Reject negative finite values and Inf in calculation variables. Zero and
% NaN are allowed; a resulting non-finite temperature is reported later.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu', 'Fe3_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', 'Fe3_cation_apfu'};

maxNames = numel(garnetVariables) + numel(micaVariables);
negativeNameBuffer = strings(maxNames, 1);
infiniteNameBuffer = strings(maxNames, 1);
nNegative = 0;
nInfinite = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nNegative = nNegative + 1;
            negativeNameBuffer(nNegative) = "Garnet." + string(variableName);
        end
        if any(isinf(variableValue(:)))
            nInfinite = nInfinite + 1;
            infiniteNameBuffer(nInfinite) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
    if ismember(variableName, data_mica.Properties.VariableNames)
        variableValue = data_mica.(variableName);
        if any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nNegative = nNegative + 1;
            negativeNameBuffer(nNegative) = "Mica." + string(variableName);
        end
        if any(isinf(variableValue(:)))
            nInfinite = nInfinite + 1;
            infiniteNameBuffer(nInfinite) = "Mica." + string(variableName);
        end
    end
end

if nNegative > 0
    invalidNames = negativeNameBuffer(1:nNegative);
    error(['HynesForest1988: mineral-composition values used by the ' ...
           'thermometer must be >= 0. Negative finite value(s) were found ' ...
           'in: ' char(strjoin(invalidNames, ', ')) '.']);
end

if nInfinite > 0
    invalidNames = infiniteNameBuffer(1:nInfinite);
    error(['HynesForest1988: Inf is not permitted in calculation ' ...
           'inputs. Inf was found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_garnet, data_mica, P_kbar)
% calcTemp
% Calculate equation (9) for one Garnet-Mica pair and return one row for
% each input pressure. Pressure is retained but is not used by equation (9).

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;

R_cal = 1.987;
row.R_cal = repmat(R_cal, nP, 1);

grt = prepareMineralRow(data_garnet, 'Garnet');
mica = prepareMineralRow(data_mica, 'Mica');

% Mineral compositions are scalars for the selected pair. Any NaN present
% in a calculation input propagates through the following expressions.
sumDiv_grt_scalar = grt.FeUsed + grt.Mg + grt.Mn + grt.Ca;
sumFeMg_mica_scalar = mica.FeUsed + mica.Mg;

XFe_grt_scalar = grt.FeUsed ./ sumDiv_grt_scalar;
XMg_grt_scalar = grt.Mg ./ sumDiv_grt_scalar;
XMn_grt_scalar = grt.Mn ./ sumDiv_grt_scalar;
XCa_grt_scalar = grt.Ca ./ sumDiv_grt_scalar;

XFe_mica_scalar = mica.FeUsed ./ sumFeMg_mica_scalar;
XMg_mica_scalar = mica.Mg ./ sumFeMg_mica_scalar;

Kd_scalar = (XMg_mica_scalar .* XFe_grt_scalar) ./ ...
    (XFe_mica_scalar .* XMg_grt_scalar);
lnKd_scalar = log(Kd_scalar);

sumFeMg_grt_scalar = XFe_grt_scalar + XMg_grt_scalar;
W_FeMg_scalar = 200 .* (XMg_grt_scalar ./ sumFeMg_grt_scalar) ...
    + 2500 .* (XFe_grt_scalar ./ sumFeMg_grt_scalar);

[T_K_scalar, lnK_scalar, solverFlagScalar] = solveTempHF1988( ...
    lnKd_scalar, XFe_grt_scalar, XMg_grt_scalar, XMn_grt_scalar, ...
    XCa_grt_scalar, W_FeMg_scalar, R_cal);
T_C_scalar = T_K_scalar - 273.15;

% Replicate all pressure-independent outputs to the pressure-vector length.
row.Fe2_grt = repmat(grt.Fe2, nP, 1);
row.Fe3_grt = repmat(grt.Fe3, nP, 1);
row.FeUsed_grt = repmat(grt.FeUsed, nP, 1);
row.Mg_grt = repmat(grt.Mg, nP, 1);
row.Mn_grt = repmat(grt.Mn, nP, 1);
row.Ca_grt = repmat(grt.Ca, nP, 1);
row.Al_grt = repmat(grt.Al, nP, 1);
row.Si_grt = repmat(grt.Si, nP, 1);

row.Fe2_mica = repmat(mica.Fe2, nP, 1);
row.Fe3_mica = repmat(mica.Fe3, nP, 1);
row.FeUsed_mica = repmat(mica.FeUsed, nP, 1);
row.Mg_mica = repmat(mica.Mg, nP, 1);
row.Mn_mica = repmat(mica.Mn, nP, 1);
row.Ca_mica = repmat(mica.Ca, nP, 1);
row.Ti_mica = repmat(mica.Ti, nP, 1);
row.Al_mica = repmat(mica.Al, nP, 1);
row.Si_mica = repmat(mica.Si, nP, 1);
row.K_mica = repmat(mica.K, nP, 1);
row.Na_mica = repmat(mica.Na, nP, 1);

row.sumDiv_grt = repmat(sumDiv_grt_scalar, nP, 1);
row.sumFeMg_mica = repmat(sumFeMg_mica_scalar, nP, 1);
row.XFe_grt = repmat(XFe_grt_scalar, nP, 1);
row.XMg_grt = repmat(XMg_grt_scalar, nP, 1);
row.XMn_grt = repmat(XMn_grt_scalar, nP, 1);
row.XCa_grt = repmat(XCa_grt_scalar, nP, 1);
row.XFe_mica = repmat(XFe_mica_scalar, nP, 1);
row.XMg_mica = repmat(XMg_mica_scalar, nP, 1);

row.W_FeMg = repmat(W_FeMg_scalar, nP, 1);
row.Kd = repmat(Kd_scalar, nP, 1);
row.lnKd = repmat(lnKd_scalar, nP, 1);
row.lnK = repmat(lnK_scalar, nP, 1);
row.solverFlag = repmat(string(solverFlagScalar), nP, 1);

row.TEq9_K = repmat(T_K_scalar, nP, 1);
row.TEq9_C = repmat(T_C_scalar, nP, 1);

end

function [Tsol, lnKsol, solverFlag] = solveTempHF1988( ...
        lnKd, XFe, XMg, XMn, XCa, W_FeMg, R_cal)
% solveTempHF1988
% Solve implicit equation (9) for temperature in Kelvin.

Tsol = NaN;
lnKsol = NaN;
solverFlag = "failed";

if any(~isfinite([lnKd, XFe, XMg, XMn, XCa, W_FeMg, R_cal]))
    solverFlag = "invalid-input";
    return;
end

residualFun = @(T) calcResidual(T, lnKd, XFe, XMg, XMn, XCa, ...
    W_FeMg, R_cal);

% Search a broad metamorphic-temperature interval for a bracketed root.
Tgrid = linspace(273.15, 1200, 400);
residGrid = NaN(size(Tgrid));

for i = 1:numel(Tgrid)
    residGrid(i) = residualFun(Tgrid(i));
end

validMask = isfinite(residGrid);
if sum(validMask) < 2
    return;
end

Tvalid = Tgrid(validMask);
Rvalid = residGrid(validMask);
signChangeIdx = find(Rvalid(1:end-1) .* Rvalid(2:end) <= 0, 1, 'first');

if ~isempty(signChangeIdx)
    Tlo = Tvalid(signChangeIdx);
    Thi = Tvalid(signChangeIdx + 1);

    try
        Tcandidate = fzero(residualFun, [Tlo, Thi]);
        if isfinite(Tcandidate) && Tcandidate > 0
            lnKcandidate = calcLnK(Tcandidate, lnKd, XFe, XMg, XMn, ...
                XCa, W_FeMg, R_cal);
            if isfinite(lnKcandidate)
                Tsol = Tcandidate;
                lnKsol = lnKcandidate;
                solverFlag = "fzero-bracket";
                return;
            end
        end
    catch
    end
end

% Fallback 1: start fzero from the smallest absolute residual on the grid.
[~, idxMin] = min(abs(Rvalid));
Tinit = Tvalid(idxMin);

try
    Tcandidate = fzero(residualFun, Tinit);
    if isfinite(Tcandidate) && Tcandidate > 0
        lnKcandidate = calcLnK(Tcandidate, lnKd, XFe, XMg, XMn, ...
            XCa, W_FeMg, R_cal);
        if isfinite(lnKcandidate)
            Tsol = Tcandidate;
            lnKsol = lnKcandidate;
            solverFlag = "fzero-single";
            return;
        end
    end
catch
end

% Fallback 2: minimize the absolute residual within the search interval.
try
    objectiveFun = @(T) abs(residualFun(T));
    Tcandidate = fminbnd(objectiveFun, 273.15, 1200);
    residualCandidate = residualFun(Tcandidate);

    if isfinite(Tcandidate) && isfinite(residualCandidate) && ...
            abs(residualCandidate) < 1e-6
        lnKcandidate = calcLnK(Tcandidate, lnKd, XFe, XMg, XMn, ...
            XCa, W_FeMg, R_cal);
        if isfinite(lnKcandidate)
            Tsol = Tcandidate;
            lnKsol = lnKcandidate;
            solverFlag = "fminbnd";
        end
    end
catch
end

end

function residual = calcResidual(T, lnKd, XFe, XMg, XMn, XCa, ...
        W_FeMg, R_cal)
% calcResidual
% Difference between trial T and the right-hand side of equation (9).

if ~(isscalar(T) && isfinite(T) && T > 0)
    residual = NaN;
    return;
end

lnK = calcLnK(T, lnKd, XFe, XMg, XMn, XCa, W_FeMg, R_cal);
rhs = 4790 ./ (lnK + 4.13);
residual = T - rhs;

end

function lnK = calcLnK(T, lnKd, XFe, XMg, XMn, XCa, W_FeMg, R_cal)
% calcLnK
% Calculate the temperature-dependent corrected lnK in equation (9).

if ~(isscalar(T) && isfinite(T) && T > 0)
    lnK = NaN;
    return;
end

y = 844 ./ T;

lnK = lnKd ...
    + (0.8 .* W_FeMg - W_FeMg .* (XFe - XMg) - 3000 .* XMn) ...
      ./ (R_cal .* T) ...
    - 2.978 .* XCa .* y ...
    + 5.906 .* (XCa .* y) .^ 2;

end

function mineral = prepareMineralRow(dataTable, mineralLabel)
% prepareMineralRow
% Extract one-row mineral data. NaN in a present calculation variable is
% retained. Missing output-only variables are represented by NaN.

if height(dataTable) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();
mineral.Fe2 = getRequiredScalar(dataTable, 'Fe_cation_apfu', mineralLabel);
mineral.Mg = getRequiredScalar(dataTable, 'Mg_cation_apfu', mineralLabel);

% The paper treated all analysed Fe as ferrous. A missing Fe3 column means
% no separately stored ferric component; a present NaN remains NaN.
mineral.Fe3 = getOptionalScalar(dataTable, 'Fe3_cation_apfu', ...
    mineralLabel, 0);
mineral.FeUsed = mineral.Fe2 + mineral.Fe3;

if strcmp(mineralLabel, 'Garnet')
    mineral.Mn = getRequiredScalar(dataTable, 'Mn_cation_apfu', mineralLabel);
    mineral.Ca = getRequiredScalar(dataTable, 'Ca_cation_apfu', mineralLabel);
else
    mineral.Mn = getOptionalScalar(dataTable, 'Mn_cation_apfu', ...
        mineralLabel, NaN);
    mineral.Ca = getOptionalScalar(dataTable, 'Ca_cation_apfu', ...
        mineralLabel, NaN);
end

mineral.Ti = getOptionalScalar(dataTable, 'Ti_cation_apfu', ...
    mineralLabel, NaN);
mineral.Al = getOptionalScalar(dataTable, 'Al_cation_apfu', ...
    mineralLabel, NaN);
mineral.Si = getOptionalScalar(dataTable, 'Si_cation_apfu', ...
    mineralLabel, NaN);
mineral.K = getOptionalScalar(dataTable, 'K_cation_apfu', ...
    mineralLabel, NaN);
mineral.Na = getOptionalScalar(dataTable, 'Na_cation_apfu', ...
    mineralLabel, NaN);

end

function value = getRequiredScalar(dataTable, variableName, mineralLabel)
% getRequiredScalar
% Return a required scalar without changing NaN to zero.

value = dataTable.(variableName);
if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end

function value = getOptionalScalar(dataTable, variableName, ...
        mineralLabel, missingValue)
% getOptionalScalar
% Return a present optional scalar unchanged. Use the supplied missingValue
% only when the column itself is absent; a present NaN remains NaN.

if ismember(variableName, dataTable.Properties.VariableNames)
    value = dataTable.(variableName);
    if ~isnumeric(value) || ~isscalar(value)
        error('%s.%s must be a numeric scalar in the selected row.', ...
            mineralLabel, variableName);
    end
else
    value = missingValue;
end

end
