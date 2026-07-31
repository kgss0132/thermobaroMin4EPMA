function results = KleemannReinhardt1994(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/KleemannReinhardt1994.m
% Tested with MATLAB R2024b
%
% Garnet-biotite Fe-Mg exchange thermometer
% Kleemann, U. and Reinhardt, J. (1994)
% European Journal of Mineralogy, 6, 925-941
% DOI: https://doi.org/10.1127/ejm/6/6/0925
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% (biotite) analysis and calculates temperature using equation (22) of
% Kleemann & Reinhardt (1994). The formulation explicitly includes the
% effects of octahedral Al and Ti in biotite.
%
% The function accepts pressure as either a scalar or a vector. Therefore,
% it can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for every pressure
% value for each user-selected Garnet-Mica pair.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Kleemann & Reinhardt (1994) calibrated their thermometer using two
% published experimental datasets (Experimental data base, pp. 926-927):
%
%   Ferry & Spear (1978):
%     Temperature : 550-800 degreeC
%     Pressure    : 2.07 kbar
%     Data        : 16 experiments using synthetic almandine-pyrope garnet
%                   and annite-phlogopite biotite compositions
%
%   Perchuk & Lavrent'eva (1983):
%     Temperature : 575-950 degreeC
%     Pressure    : 6 kbar
%     Data        : 35 Garnet-Biotite pairs with reported biotite AlVI
%     Garnet      : Mg/(Mg+Fe) = 0.05-0.76
%     Biotite     : Mg/(Mg+Fe) = 0.22-0.83
%
% Thus, the combined experimental temperature envelope is 550-950 degreeC.
% Pressure was tested at two discrete values, 2.07 and 6 kbar; it was not
% continuously calibrated over all pressures between these values.
%
% Important application cautions stated in the original paper:
% - Natural-rock tests were most successful for greenschist- to upper-
%   amphibolite-facies metapelites (pp. 937-939).
% - Retrograde Fe-Mg exchange between garnet and biotite is identified as
%   the most serious source of error. Mineral pairs should represent the
%   same equilibrium stage and preferably be adjacent or in contact;
%   garnet core/rim and biotite textural position must be considered
%   (p. 937).
% - Many granulite-facies rocks gave temperatures too low to represent the
%   thermal peak. The paper notes limited experimental constraint above
%   800 degreeC and emphasizes partial re-equilibration during cooling
%   (p. 937).
% - The effects of Fe3+ in garnet and biotite and of F-Cl substitution for
%   OH in biotite were not included in the calibration (p. 939).
% - The Ti term was not independently fitted from the incompletely reported
%   Perchuk & Lavrent'eva biotite compositions; it incorporates parameters
%   derived from Sengupta et al. (1990) (pp. 928 and 933).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside 2.07-6 kbar,
%   2) a finite calculated temperature is outside 550-950 degreeC,
%   3) a finite calculated temperature is above 800 degreeC,
%   4) a required thermometer input is NaN, or
%   5) a calculated temperature is NaN or Inf.
%
% IMPORTANT IMPLEMENTATION LIMITATION
% Equation (22) requires the non-ideal Ca-Mg-Fe-Mn garnet activity model of
% Berman (1990). The supplied source file did not contain that model and
% instead used an ideal-garnet approximation. This modified file preserves
% that explicit approximation: dW_H_Grt = dW_S_Grt = dW_V_Grt = 0.
% Consequently, results are not an exact reproduction of the complete
% Kleemann & Reinhardt (1994) calibration for non-ideal garnet compositions.
% A non-stopping fprintf notice is printed once per function call.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Mica   : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The following normalized cation
% variables are required because they enter the calculation directly:
%
%   Garnet table variables:
%     Fe_cation_apfu       % Fe2+ used in the Fe2+-Mg exchange reaction
%     Mg_cation_apfu
%     Mn_cation_apfu
%     Ca_cation_apfu
%
%   Mica table variables (biotite normalized to 11 O equivalent):
%     Fe_cation_apfu       % Fe2+ used in the Fe2+-Mg exchange reaction
%     Mg_cation_apfu
%     Mn_cation_apfu
%     Ti_cation_apfu
%     Al_cation_apfu
%     Si_cation_apfu
%
% Optional variables retained in the output when present:
%   Fe3_cation_apfu, K_cation_apfu, Na_cation_apfu, Ca_cation_apfu
%
% Fe3_cation_apfu is not added to Fe_cation_apfu because the original
% calibration is formulated for Fe2+-Mg exchange and explicitly states
% that the effect of Fe3+ was not considered (p. 939).
%
% All finite mineral-composition values used by the thermometer must be
% greater than or equal to zero. Negative finite values and Inf values stop
% the calculation with an error. NaN values are retained as missing values,
% propagated through the calculation, and reported by non-stopping fprintf
% warnings. Missing required columns are not interpreted as zero.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Exchange coefficient:
%   K_D = (Mg/Fe2+)_Grt / (Mg/Fe2+)_Bt
%
% Garnet site fractions:
%   XAlm = Fe2+ / (Fe2+ + Mg + Mn + Ca)
%   XPrp = Mg   / (Fe2+ + Mg + Mn + Ca)
%   XSps = Mn   / (Fe2+ + Mg + Mn + Ca)
%   XGrs = Ca   / (Fe2+ + Mg + Mn + Ca)
%
% Biotite octahedral fractions (Appendix, p. 940):
%   AlIV  = max(4 - Si, 0)
%   AlVI  = max(Al_total - AlIV, 0)
%   OctSum = AlVI + Ti + Mg + Fe2+ + Mn
%   X_Al_Bt = AlVI / OctSum
%   X_Ti_Bt = Ti   / OctSum
%
% Kleemann & Reinhardt (1994), equation (22):
%
%   T(K) =
%   ( 20253 + dW_H_Grt + 77785*X_Al_Bt - 18138*X_Ti_Bt ...
%     + (dV_0 + dW_V_Grt)*P_bar ) ...
%   / ...
%   ( 10.66 - R*ln(K_D) + dW_S_Grt ...
%     + 94.1*X_Al_Bt - 11.7*X_Ti_Bt )
%
% R and all energetic constants are used consistently in J-based units:
%   R = 8.314462618 J mol^-1 K^-1
%
% The standard volume term dV_0 is calculated from equations (7), (8), and
% (10) on p. 932 rather than mixing cal- and J-based constants.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = KleemannReinhardt1994(rawdata_struct, P_kbar)
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
    error('KleemannReinhardt1994 requires (rawdata_struct, P_kbar).');
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
requiredMicaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Si_cation_apfu'};

validateRequiredVariables(dataset_grt, requiredGarnetVariables, 'Garnet');
validateRequiredVariables(dataset_mica, requiredMicaVariables, 'Mica');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each calculation as one table block and concatenate the blocks only
% once after the interactive loop. This avoids repeatedly growing and
% copying the complete results table on every loop iteration.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Combined experimental limits from the two calibration datasets.
calibrationT_min_degC = 550;
calibrationT_max_degC = 950;
highTemperatureCaution_degC = 800;
calibrationP_min_kbar = 2.07;
calibrationP_max_kbar = 6;

pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;
implementationNoticeIssued = false;

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

    % NaN is reported but deliberately retained and propagated.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_mica);

    % Reject negative finite values and Inf, but permit zero and NaN.
    validateNonNegativeInputs(selectedData_grt, selectedData_mica);

    row = calcTemp(selectedData_grt, selectedData_mica, P_kbar);

    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_mica = repmat(string(selectedCode_mica), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_mica'}, 'Before', 1);

    % Append this table block to the preallocated cell buffer. The capacity
    % is doubled only when required, rather than changing the complete
    % results table size on every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': KR94 = ' num2str(row.TKR94_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': KR94 = ' num2str(row.TKR94_C(1)) ' to ' ...
            num2str(row.TKR94_C(end)) ' degreeC']);
    end

    % Report the inherited ideal-garnet approximation once per call.
    if ~implementationNoticeIssued
        fprintf(2, ...
            ['CAUTION: KleemannReinhardt1994 currently uses an ideal-garnet ' ...
             'approximation because the Berman (1990) garnet activity model was ' ...
             'not implemented in the supplied source file. Results are not an exact ' ...
             'reproduction of the complete KR94 calibration.\n']);
        implementationNoticeIssued = true;
    end

    % Pressure warning: the experiments were performed at two discrete
    % pressure values, 2.07 and 6 kbar.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental pressure envelope ' ...
             'of Kleemann & Reinhardt (1994): 2.07-6 kbar. The calibration used ' ...
             'experiments at two discrete pressures (2.07 and 6 kbar), not a ' ...
             'continuous pressure series (pp. 926-927). %d of %d pressure point(s) ' ...
             'are outside the envelope; input range = %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), numel(P_kbar), ...
            min(P_kbar), max(P_kbar));
        pressureWarningIssued = true;
    end

    % Temperature warning for the combined experimental envelope.
    finiteTemperature = isfinite(row.TKR94_C);
    temperatureOutsideCalibration = finiteTemperature & ...
        (row.TKR94_C < calibrationT_min_degC | ...
         row.TKR94_C > calibrationT_max_degC);

    if any(temperatureOutsideCalibration)
        finiteValues = row.TKR94_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the combined experimental ' ...
             'calibration envelope of Kleemann & Reinhardt (1994): 550-950 degreeC ' ...
             '(pp. 926-927). %d of %d finite temperature point(s) are outside the ' ...
             'range; calculated finite range = %.4g-%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideCalibration), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)));
    end

    % Additional high-temperature caution from the natural-rock discussion.
    highTemperature = finiteTemperature & ...
        row.TKR94_C > highTemperatureCaution_degC;
    if any(highTemperature)
        fprintf(2, ...
            ['CAUTION: %d calculated temperature point(s) exceed 800 degreeC for ' ...
             '%s & %s. Kleemann & Reinhardt (1994, p. 937) note limited high-' ...
             'temperature constraint and that granulite-facies samples may record ' ...
             'retrograde re-equilibration rather than peak temperature.\n'], ...
            sum(highTemperature), char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)));
    end

    % Print a non-stopping warning when required inputs contain NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         NaN was retained and propagated; it was not replaced by zero.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain NaN/Inf results and report them without stopping.
    invalidTemperature = ~isfinite(row.TKR94_C);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)), ...
            sum(invalidTemperature), numel(row.TKR94_C), ...
            sum(isnan(row.TKR94_C)), sum(isinf(row.TKR94_C)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'KleemannReinhardt1994', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate the buffered blocks once. Return an empty table if the user
% canceled before completing a calculation.
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
% Confirm that every calculation variable exists and is numeric. Missing
% required variables are not silently replaced by zero.

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

end

function nanInputNames = findNaNInputs(data_garnet, data_mica)
% findNaNInputs
% Return required thermometer input names containing NaN. The output buffer
% is preallocated so its size is not increased on each loop iteration.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Si_cation_apfu'};

nameBuffer = strings(numel(garnetVariables) + numel(micaVariables), 1);
nNames = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    variableValue = data_garnet.(variableName);
    if any(isnan(variableValue(:)))
        nNames = nNames + 1;
        nameBuffer(nNames) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
    variableValue = data_mica.(variableName);
    if any(isnan(variableValue(:)))
        nNames = nNames + 1;
        nameBuffer(nNames) = "Mica." + string(variableName);
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_garnet, data_mica)
% validateNonNegativeInputs
% Reject negative finite values and Inf in required calculation variables.
% Zero and NaN are allowed; NaN propagates and is reported non-destructively.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Mn_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Si_cation_apfu'};

maxNames = numel(garnetVariables) + numel(micaVariables);
negativeNameBuffer = strings(maxNames, 1);
infiniteNameBuffer = strings(maxNames, 1);
nNegative = 0;
nInfinite = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
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

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
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

if nNegative > 0
    invalidNames = negativeNameBuffer(1:nNegative);
    error(['KleemannReinhardt1994: mineral-composition values used by ' ...
           'the thermometer must be >= 0. Negative finite value(s) were ' ...
           'found in: ' char(strjoin(invalidNames, ', ')) '.']);
end

if nInfinite > 0
    invalidNames = infiniteNameBuffer(1:nInfinite);
    error(['KleemannReinhardt1994: Inf is not permitted in required ' ...
           'thermometer inputs. Inf was found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_garnet, data_mica, P_kbar)
% calcTemp
% Compute Kleemann & Reinhardt (1994) temperatures for one garnet row, one
% mica row, and a scalar or vector of pressure values.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

row = table();

% All constants in equation (22) are treated in J-based units.
R_J = 8.314462618;

row.P_kbar = P_kbar;
row.P_bar = P_bar;
row.R_J = repmat(R_J, nP, 1);

% Prepare one-row mineral data. Optional missing values remain NaN.
grt = prepareMineralRow(data_garnet, 'Garnet');
mica = prepareMineralRow(data_mica, 'Mica');

% Fe_cation_apfu is treated as Fe2+; Fe3+ is retained for reporting only.
MgFe_grt_scalar = grt.Mg ./ grt.Fe2;
MgFe_mica_scalar = mica.Mg ./ mica.Fe2;
K_D_scalar = MgFe_grt_scalar ./ MgFe_mica_scalar;
lnK_D_scalar = log(K_D_scalar);

% Garnet mole fractions on the divalent site basis.
[XFe_grt_scalar, XMg_grt_scalar, XMn_grt_scalar, ...
    XCa_grt_scalar, Xsum_grt_scalar] = calcGarnetSiteFractions(grt);

% Biotite octahedral fractions.
[X_Al_Bt_scalar, X_Ti_Bt_scalar, AlIV_Bt_scalar, ...
    AlVI_Bt_scalar, OctSum_Bt_scalar] = calcBiotiteSiteFractions(mica);

% Standard molar volume change from equations (7), (8), and (10), p. 932.
dV_0_scalar = calcStandardVolumeChange(AlVI_Bt_scalar);

% Garnet non-ideality terms. The supplied source used the explicit ideal-
% garnet approximation retained by this modified file.
[dW_H_Grt_scalar, dW_S_Grt_scalar, dW_V_Grt_scalar] = ...
    calcKR94GarnetNonIdealTerms(XFe_grt_scalar, XMg_grt_scalar, ...
        XMn_grt_scalar, XCa_grt_scalar);

% Replicate scalar mineral terms to match the pressure-vector length.
MgFe_grt = repmat(MgFe_grt_scalar, nP, 1);
MgFe_mica = repmat(MgFe_mica_scalar, nP, 1);
K_D = repmat(K_D_scalar, nP, 1);
lnK_D = repmat(lnK_D_scalar, nP, 1);

XFe_grt = repmat(XFe_grt_scalar, nP, 1);
XMg_grt = repmat(XMg_grt_scalar, nP, 1);
XMn_grt = repmat(XMn_grt_scalar, nP, 1);
XCa_grt = repmat(XCa_grt_scalar, nP, 1);
Xsum_grt = repmat(Xsum_grt_scalar, nP, 1);

X_Al_Bt = repmat(X_Al_Bt_scalar, nP, 1);
X_Ti_Bt = repmat(X_Ti_Bt_scalar, nP, 1);
AlIV_Bt = repmat(AlIV_Bt_scalar, nP, 1);
AlVI_Bt = repmat(AlVI_Bt_scalar, nP, 1);
OctSum_Bt = repmat(OctSum_Bt_scalar, nP, 1);

dV_0 = repmat(dV_0_scalar, nP, 1);
dW_H_Grt = repmat(dW_H_Grt_scalar, nP, 1);
dW_S_Grt = repmat(dW_S_Grt_scalar, nP, 1);
dW_V_Grt = repmat(dW_V_Grt_scalar, nP, 1);

% Equation (22), Kleemann & Reinhardt (1994), p. 934.
numerator = 20253 + dW_H_Grt + 77785 .* X_Al_Bt ...
    - 18138 .* X_Ti_Bt + (dV_0 + dW_V_Grt) .* P_bar;

denominator = 10.66 - R_J .* lnK_D + dW_S_Grt ...
    + 94.1 .* X_Al_Bt - 11.7 .* X_Ti_Bt;

% Do not replace NaN/Inf with zero or stop the calculation. Non-finite
% results are retained in the output and reported by the caller.
TKR94_K = numerator ./ denominator;
TKR94_C = TKR94_K - 273.15;

% Pack pressure-independent mineral values as nP-by-1 columns.
row.Fe2_grt = repmat(grt.Fe2, nP, 1);
row.Fe3_grt = repmat(grt.Fe3, nP, 1);
row.FeUsed_grt = repmat(grt.Fe2, nP, 1);
row.Mg_grt = repmat(grt.Mg, nP, 1);
row.Mn_grt = repmat(grt.Mn, nP, 1);
row.Ca_grt = repmat(grt.Ca, nP, 1);
row.Al_grt = repmat(grt.Al, nP, 1);
row.Si_grt = repmat(grt.Si, nP, 1);

row.Fe2_mica = repmat(mica.Fe2, nP, 1);
row.Fe3_mica = repmat(mica.Fe3, nP, 1);
row.FeUsed_mica = repmat(mica.Fe2, nP, 1);
row.Mg_mica = repmat(mica.Mg, nP, 1);
row.Mn_mica = repmat(mica.Mn, nP, 1);
row.Ca_mica = repmat(mica.Ca, nP, 1);
row.Ti_mica = repmat(mica.Ti, nP, 1);
row.Al_mica = repmat(mica.Al, nP, 1);
row.Si_mica = repmat(mica.Si, nP, 1);
row.K_mica = repmat(mica.K, nP, 1);
row.Na_mica = repmat(mica.Na, nP, 1);

row.XFe_grt = XFe_grt;
row.XMg_grt = XMg_grt;
row.XMn_grt = XMn_grt;
row.XCa_grt = XCa_grt;
row.Xsum_grt = Xsum_grt;

row.AlIV_Bt = AlIV_Bt;
row.AlVI_Bt = AlVI_Bt;
row.OctSum_Bt = OctSum_Bt;
row.X_Al_Bt = X_Al_Bt;
row.X_Ti_Bt = X_Ti_Bt;

row.MgFe_grt = MgFe_grt;
row.MgFe_mica = MgFe_mica;
row.K_D = K_D;
row.lnK_D = lnK_D;

row.dW_H_Grt = dW_H_Grt;
row.dW_S_Grt = dW_S_Grt;
row.dW_V_Grt = dW_V_Grt;
row.dV_0 = dV_0;
row.numerator = numerator;
row.denominator = denominator;

row.TKR94_K = TKR94_K;
row.TKR94_C = TKR94_C;

end

function mineral = prepareMineralRow(dataTable, mineralLabel)
% prepareMineralRow
% Extract one-row cation data. Missing optional variables are represented by
% NaN, not zero. Required variables were checked before selection.

if height(dataTable) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();
mineral.Fe2 = getRequiredScalar(dataTable, 'Fe_cation_apfu', mineralLabel);
mineral.Mg = getRequiredScalar(dataTable, 'Mg_cation_apfu', mineralLabel);
mineral.Mn = getRequiredScalar(dataTable, 'Mn_cation_apfu', mineralLabel);

if strcmp(mineralLabel, 'Mica')
    mineral.Ca = getOptionalScalarOrNaN(dataTable, 'Ca_cation_apfu', mineralLabel);
    mineral.Ti = getRequiredScalar(dataTable, 'Ti_cation_apfu', mineralLabel);
    mineral.Al = getRequiredScalar(dataTable, 'Al_cation_apfu', mineralLabel);
    mineral.Si = getRequiredScalar(dataTable, 'Si_cation_apfu', mineralLabel);
else
    mineral.Ca = getRequiredScalar(dataTable, 'Ca_cation_apfu', mineralLabel);
    mineral.Ti = getOptionalScalarOrNaN(dataTable, 'Ti_cation_apfu', mineralLabel);
    mineral.Al = getOptionalScalarOrNaN(dataTable, 'Al_cation_apfu', mineralLabel);
    mineral.Si = getOptionalScalarOrNaN(dataTable, 'Si_cation_apfu', mineralLabel);
end

mineral.Fe3 = getOptionalScalarOrNaN(dataTable, 'Fe3_cation_apfu', mineralLabel);
mineral.K = getOptionalScalarOrNaN(dataTable, 'K_cation_apfu', mineralLabel);
mineral.Na = getOptionalScalarOrNaN(dataTable, 'Na_cation_apfu', mineralLabel);

end

function [XFe, XMg, XMn, XCa, Xsum] = calcGarnetSiteFractions(grt)
% calcGarnetSiteFractions
% Calculate almandine, pyrope, spessartine, and grossular fractions using
% Fe2+ only, as defined in the Appendix of Kleemann & Reinhardt (1994).

Xsum = grt.Fe2 + grt.Mg + grt.Mn + grt.Ca;
XFe = grt.Fe2 ./ Xsum;
XMg = grt.Mg ./ Xsum;
XMn = grt.Mn ./ Xsum;
XCa = grt.Ca ./ Xsum;

end

function [X_Al_Bt, X_Ti_Bt, AlIV, AlVI, OctSum] = ...
        calcBiotiteSiteFractions(mica)
% calcBiotiteSiteFractions
% Calculate octahedral Al and Ti fractions on the basis defined on p. 940.
% Explicit conditional clipping is used so that NaN remains NaN.

AlIV = 4 - mica.Si;
if isfinite(AlIV) && AlIV < 0
    AlIV = 0;
end

AlVI = mica.Al - AlIV;
if isfinite(AlVI) && AlVI < 0
    AlVI = 0;
end

OctSum = mica.Fe2 + mica.Mg + mica.Mn + mica.Ti + AlVI;
X_Al_Bt = AlVI ./ OctSum;
X_Ti_Bt = mica.Ti ./ OctSum;

end

function dV_0 = calcStandardVolumeChange(AlVI)
% calcStandardVolumeChange
% Calculate the standard molar volume change of the exchange reaction using
% equations (7), (8), and (10) of Kleemann & Reinhardt (1994), p. 932.

z = AlVI;

% Garnet end-member molar volumes (cm3/mol) used by KR94.
V_Prp = 113.20;
V_Alm = 115.11;

% Biotite end-member unit-cell volumes from equation (7), followed by the
% conversion to molar volume in equation (8).
Vcell_FeBt_A3 = 514.2 - 14.3 .* z;
Vcell_MgBt_A3 = 514.2 - 17.3 - 14.3 .* z + 7.3 .* z;
V_FeBt = 0.6022 .* Vcell_FeBt_A3 ./ 2;
V_MgBt = 0.6022 .* Vcell_MgBt_A3 ./ 2;

% 0.1 converts cm3 mol^-1 to J mol^-1 bar^-1.
dV_0 = 0.1 .* ((3 - z) .* V_Prp + 3 .* V_FeBt ...
    - (3 - z) .* V_Alm - 3 .* V_MgBt) ./ (3 .* (3 - z));

end

function [dW_H_Grt, dW_S_Grt, dW_V_Grt] = ...
        calcKR94GarnetNonIdealTerms(XFe, XMg, XMn, XCa)
% calcKR94GarnetNonIdealTerms
% Return the explicitly documented ideal-garnet approximation inherited
% from the supplied source file. Replace this function with a verified
% implementation of the Berman (1990) model for a full KR94 calculation.

if any(isnan([XFe, XMg, XMn, XCa]))
    dW_H_Grt = NaN;
    dW_S_Grt = NaN;
    dW_V_Grt = NaN;
    return;
end

dW_H_Grt = 0;
dW_S_Grt = 0;
dW_V_Grt = 0;

end

function value = getRequiredScalar(dataTable, variableName, mineralLabel)
% getRequiredScalar
% Return one required scalar without changing NaN to zero.

value = dataTable.(variableName);
if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end

function value = getOptionalScalarOrNaN(dataTable, variableName, mineralLabel)
% getOptionalScalarOrNaN
% Return one optional scalar; use NaN when the column is absent.

if ismember(variableName, dataTable.Properties.VariableNames)
    value = dataTable.(variableName);
    if ~isnumeric(value) || ~isscalar(value)
        error('%s.%s must be a numeric scalar in the selected row.', ...
            mineralLabel, variableName);
    end
    if isinf(value) || value < 0
        error('%s.%s must be NaN or a finite non-negative value.', ...
            mineralLabel, variableName);
    end
else
    value = NaN;
end

end
