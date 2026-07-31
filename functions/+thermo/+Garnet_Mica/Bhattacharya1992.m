function results = Bhattacharya1992(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/Bhattacharya1992.m
% Tested with MATLAB R2024b
%
% Garnet-biotite Fe-Mg exchange thermometer
% Bhattacharya, A., Mohanty, L., Maji, A., Sen, S.K. and Raith, M. (1992)
% Contributions to Mineralogy and Petrology, 111, 87-93
% DOI: https://doi.org/10.1007/BF00296580
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% analysis (selected by the user from tables) and calculates temperature
% using the Bhattacharya et al. (1992) garnet-biotite Fe-Mg exchange
% thermometer.
%
% Two temperature formulations are returned:
%   1) T_HW : based on the Hackler and Wood (1989) garnet mixing model
%   2) T_GS : based on the Ganguly and Saxena (1984) garnet mixing model
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Garnet-Mica pair, and appends results
% into a single output table.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Bhattacharya et al. (1992) reformulated the garnet-biotite thermometer
% using the experimental Fe-Mg partitioning data of Ferry and Spear (1978)
% and Perchuk and Lavrent'eva (1983). The experimental data used in Table 1
% cover the following conditions:
%
%   Temperature : approximately 575-950 degreeC
%   Pressure    : 2.07 and 6 kbar (two discrete experimental pressures)
%   XFe_biotite : approximately 0.17-0.78
%   XMg_garnet  : approximately 0.05-0.77
%
% The abstract states a temperature range of 550-950 degreeC (p. 87),
% whereas the tabulated experiments used for the formulation span
% approximately 575-950 degreeC (Table 1, pp. 90-91). In addition, the
% parameters at 575 degreeC (848 K) and 950 degreeC (1223 K) have large
% uncertainties and were excluded from the principal least-squares fit;
% the better-constrained fitted interval is 600-900 degreeC (p. 88).
%
% The pressure calibration is not a densely sampled continuous interval.
% It combines experiments conducted at 2.07 kbar and 6 kbar (pp. 87-88;
% Table 1, pp. 90-91). Calculations between these pressures are therefore
% interpolations between two experimental pressure conditions, and
% calculations outside 2.07-6 kbar are pressure extrapolations.
%
% The Ca-bearing garnet formulation adopts an approximation intended for
% XCa_garnet < 0.25. The adopted approximations change calculated
% temperatures by approximately 20 degreeC for garnets with
% XCa_garnet <= 0.20 (p. 91). This approximately 20 degreeC value is not a
% total uncertainty for the thermometer.
%
% Important limitations for natural samples are discussed on p. 93:
%   1) Mn in garnet affects Fe-Mg partitioning, but the required Mn-bearing
%      mixing properties were not sufficiently known for a quantitative
%      correction.
%   2) Natural biotite may contain Ti4+, Ti3+, Al3+, Fe3+, F-, Cl-, coupled
%      substitutions, and vacancies. Their effects were not fully calibrated.
%   3) Previously proposed empirical Ti-Al corrections were not considered
%      generally valid.
%   4) The thermometer should be applied to an equilibrated garnet-biotite
%      pair. The Mica table must contain biotite/phlogopite-annite data, not
%      an arbitrary mica composition.
%
% The exchange reaction is formulated for Fe2+-Mg. Consequently,
% Fe3_cation_apfu is retained in the output for inspection but is not added
% to Fe_cation_apfu in the mole fractions or K_D calculation. The effects of
% Fe3+ in natural biotite are among the uncalibrated complications discussed
% by Bhattacharya et al. (1992, p. 93).
%
% This implementation issues non-stopping warnings using fprintf when:
%   1) input pressure is outside 2.07-6 kbar,
%   2) a finite calculated T_HW or T_GS is outside 575-950 degreeC,
%   3) an input used by the thermometer contains NaN, or
%   4) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Mica   : table containing biotite analyses
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns should contain
% normalized cation data.
%
% Required variables in both Garnet and Mica tables:
%   Fe_cation_apfu          % Fe2+ used in the exchange calculation
%   Mg_cation_apfu
%
% Optional variables:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Ca_cation_apfu and Mn_cation_apfu in garnet are used by the published
% formulation. If either variable is absent, it is set to zero. If it exists
% and contains NaN, NaN is retained and propagated through the calculation.
%
% All finite mineral-composition values must be non-negative. Negative finite
% values are prohibited. NaN values are retained as missing values, propagated
% when they occur in calculation inputs, and reported by non-stopping
% warnings. A zero Fe or Mg value is allowed by the input check but may
% produce NaN or Inf, which is retained and reported after calculation.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Bhattacharya et al. (1992) incorporated non-ideal Fe-Mg mixing in both
% garnet and biotite. The two implemented formulations are:
%
% T(HW)= [20286 + 0.0193P
%          - {2080(XMgGt)^2 - 6350(XFeGt)^2
%             - 13807(XCaGt)(1-XMnGt)
%             + 8540(XFeGt)(XMgGt)(1-XCaGt-XMnGt)
%             + 4441(2XBtMg-1)}]
%         / [13.138 + 8.3143 lnK_D + 6.276XCaGt(1-XMnGt)]
%
% T(GS)= [13538 + 0.0193P
%          - {837(XMgGt)^2 - 10460(XFeGt)^2
%             - 13807XCaGt(1-XMnGt)
%             + 19246XFeGtXMgGt(1-XMnGt)
%             + 5649XCaGt(XMgGt-XFeGt)
%             + 7972(2XBtMg-1)}]
%         / [6.778 + 8.3143 lnK_D + 6.276XCaGt(1-XMnGt)]
%
% where:
%   P      : pressure in bar
%   lnK_D  : ln[(Fe2+/Mg)_Gt / (Fe2+/Mg)_Bt]
%   XFeGt  : Fe2+ / (Fe2+ + Mg + Ca + Mn) in garnet
%   XMgGt  : Mg / (Fe2+ + Mg + Ca + Mn) in garnet
%   XCaGt  : Ca / (Fe2+ + Mg + Ca + Mn) in garnet
%   XMnGt  : Mn / (Fe2+ + Mg + Ca + Mn) in garnet
%   XBtMg  : Mg / (Fe2+ + Mg) in biotite
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Bhattacharya1992(rawdata_struct, P_kbar)
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
% Basic argument checks prevent silent failures caused by missing inputs or
% invalid pressure values. Pressure is converted to a column vector so that
% fixed-pressure and pressure-range launchers use the same calculation path.
if nargin < 2
    error('Bhattacharya1992 requires (rawdata_struct, P_kbar).');
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
% modified; only the selected rows and relevant variables are read.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_grt = rawdata_struct.Garnet;
dataset_mica = rawdata_struct.Mica;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Calculation results are stored as table blocks and concatenated once after
% the interactive loop. This avoids repeatedly resizing the results table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Experimental limits used for non-stopping warnings. See the calibration
% notes above for the endpoint uncertainty and discrete pressure conditions.
calibrationT_min_degC = 575;
calibrationT_max_degC = 950;
calibrationP_min_kbar = 2.07;
calibrationP_max_kbar = 6;

% Pressure is shared by all selected mineral pairs in this function call, so
% a pressure-range warning is printed only once after the first calculation.
pressureOutsideCalibration = ...
    P_kbar < calibrationP_min_kbar | ...
    P_kbar > calibrationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
    % Assumption: the first column stores an identifier shown to the user.
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
    % Garnet and mica are selected independently; row indices are not assumed
    % to correspond between the two datasets.
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_mica = dataset_mica(selectedIdx_mica, :);

    % NaN values in calculation inputs do not stop the calculation. Their
    % variable names are collected for the warning printed after the result.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_mica);

    % Negative finite values are prohibited. NaN is deliberately allowed so
    % that it remains missing and propagates through the calculation.
    validateNonnegativeInputs(selectedData_grt, selectedData_mica);

    row = calcTemp(selectedData_grt, selectedData_mica, P_kbar);

    % Repeat the selected identifiers for every pressure row.
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_mica = repmat(string(selectedCode_mica), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_mica'}, 'Before', 1);

    % Store the result as one table block. The buffer grows geometrically only
    % when its current capacity is exhausted, rather than on every iteration.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperatures for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': HW = ' num2str(row.THW_C) ' degreeC, GS = ' ...
            num2str(row.TGS_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_mica)) ...
            ': HW = ' num2str(row.THW_C(1)) ' to ' ...
            num2str(row.THW_C(end)) ' degreeC, GS = ' ...
            num2str(row.TGS_C(1)) ' to ' num2str(row.TGS_C(end)) ' degreeC']);
    end

    % Warn once if any input pressure is outside 2.07-6 kbar. The calculation
    % is not stopped. The two endpoints represent discrete experimental
    % conditions rather than a densely sampled continuous calibration.
    if any(pressureOutsideCalibration) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the experimental pressure envelope ' ...
             'of Bhattacharya et al. (1992): 2.07-6 kbar. The calibration uses ' ...
             'two discrete pressures (2.07 and 6 kbar). %d of %d pressure point(s) ' ...
             'are outside the envelope; input range = %.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideCalibration), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Print separate range warnings for the HW and GS formulations. NaN and
    % Inf are handled by the non-finite-result warning below.
    printTemperatureCalibrationWarning( ...
        row.THW_C, 'T_HW', calibrationT_min_degC, calibrationT_max_degC, ...
        selectedCode_grt, selectedCode_mica);
    printTemperatureCalibrationWarning( ...
        row.TGS_C, 'T_GS', calibrationT_min_degC, calibrationT_max_degC, ...
        selectedCode_grt, selectedCode_mica);

    % Report NaN calculation inputs without stopping. fprintf is used for all
    % non-stopping warning messages, as in the reference Ballhaus script.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
             '         The calculation was continued, and NaN was retained as a missing value.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTHW = ~isfinite(row.THW_C);
    invalidTGS = ~isfinite(row.TGS_C);

    % Retain and report non-finite output values. This also catches cases such
    % as division by zero when a zero Fe or Mg input was supplied.
    if any(invalidTHW) || any(invalidTGS)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
             '(T_HW: %d of %d; T_GS: %d of %d; total NaN: %d; total Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_mica)), ...
            sum(invalidTHW), numel(row.THW_C), ...
            sum(invalidTGS), numel(row.TGS_C), ...
            sum(isnan(row.THW_C)) + sum(isnan(row.TGS_C)), ...
            sum(isinf(row.THW_C)) + sum(isinf(row.TGS_C)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Bhattacharya1992', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered table blocks once. Return an empty table if the
% user canceled before performing a calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_garnet, data_mica)
% findNaNInputs
% Return the names of thermometer input variables that contain NaN. Missing
% optional Ca or Mn variables are interpreted as zero and are not reported.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
micaVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

maximumNames = numel(garnetVariables) + numel(micaVariables);
nameBuffer = strings(maximumNames, 1);
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

function validateNonnegativeInputs(data_garnet, data_mica)
% validateNonnegativeInputs
% Stop when a finite mineral-composition value is negative. NaN is allowed
% and propagated, and zero is allowed as requested. Positive or zero values
% must still be numeric scalars in each selected one-row table.

garnetVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Fe3_cation_apfu', 'Mn_cation_apfu', 'Ca_cation_apfu', ...
    'Ti_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu', ...
    'K_cation_apfu', 'Na_cation_apfu'};
micaVariables = garnetVariables;

maximumNames = numel(garnetVariables) + numel(micaVariables);
nameBuffer = strings(maximumNames, 1);
nNames = 0;

for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if ~isnumeric(variableValue) || ~isscalar(variableValue) || ...
                any(isinf(variableValue(:))) || ...
                any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nNames = nNames + 1;
            nameBuffer(nNames) = "Garnet." + string(variableName);
        end
    end
end

for i = 1:numel(micaVariables)
    variableName = micaVariables{i};
    if ismember(variableName, data_mica.Properties.VariableNames)
        variableValue = data_mica.(variableName);
        if ~isnumeric(variableValue) || ~isscalar(variableValue) || ...
                any(isinf(variableValue(:))) || ...
                any(isfinite(variableValue(:)) & variableValue(:) < 0)
            nNames = nNames + 1;
            nameBuffer(nNames) = "Mica." + string(variableName);
        end
    end
end

invalidInputNames = nameBuffer(1:nNames);

if ~isempty(invalidInputNames)
    error(['Bhattacharya1992: mineral-composition values must be finite ' ...
           'non-negative scalars or NaN. Invalid value(s) were found in: ' ...
           char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(data_garnet, data_mica, P_kbar)
% calcTemp
% Compute Bhattacharya et al. (1992) HW and GS temperatures for one garnet
% row and one mica row over a scalar or vector of pressure values.
%
% Inputs:
%   data_garnet : 1-row table containing garnet cations
%   data_mica   : 1-row table containing biotite cations
%   P_kbar      : finite non-negative pressure scalar or vector in kbar
%
% Output:
%   row : table containing one row per pressure value

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

% Prepare one-row mineral data. Missing optional variables become zero;
% existing NaN values remain NaN.
garnet = prepareMineralRow(data_garnet, 'Garnet');
mica = prepareMineralRow(data_mica, 'Mica');

% Repeat composition values for every pressure point. Preconstructing all
% columns at their final height avoids changing their sizes within a loop.
Fe2_grt = repmat(garnet.Fe2, nP, 1);
Fe3_grt = repmat(garnet.Fe3, nP, 1);
FeUsed_grt = Fe2_grt;
Mg_grt = repmat(garnet.Mg, nP, 1);
Mn_grt = repmat(garnet.Mn, nP, 1);
Ca_grt = repmat(garnet.Ca, nP, 1);
Al_grt = repmat(garnet.Al, nP, 1);
Si_grt = repmat(garnet.Si, nP, 1);

Fe2_mica = repmat(mica.Fe2, nP, 1);
Fe3_mica = repmat(mica.Fe3, nP, 1);
FeUsed_mica = Fe2_mica;
Mg_mica = repmat(mica.Mg, nP, 1);
Mn_mica = repmat(mica.Mn, nP, 1);
Ca_mica = repmat(mica.Ca, nP, 1);
Ti_mica = repmat(mica.Ti, nP, 1);
Al_mica = repmat(mica.Al, nP, 1);
Si_mica = repmat(mica.Si, nP, 1);
K_mica = repmat(mica.K, nP, 1);
Na_mica = repmat(mica.Na, nP, 1);

% Garnet fractions use Fe2+ + Mg + Ca + Mn. Biotite fractions use Fe2+ + Mg.
% NaN values are intentionally allowed to propagate through these operations.
grtDivSum = FeUsed_grt + Mg_grt + Ca_grt + Mn_grt;
micaFeMgSum = FeUsed_mica + Mg_mica;

XFe_grt = FeUsed_grt ./ grtDivSum;
XMg_grt = Mg_grt ./ grtDivSum;
XCa_grt = Ca_grt ./ grtDivSum;
XMn_grt = Mn_grt ./ grtDivSum;

XFe_mica = FeUsed_mica ./ micaFeMgSum;
XMg_mica = Mg_mica ./ micaFeMgSum;

% Fe-Mg distribution coefficient.
FeMg_grt = FeUsed_grt ./ Mg_grt;
FeMg_mica = FeUsed_mica ./ Mg_mica;
K_D = FeMg_grt ./ FeMg_mica;
lnK_D = log(K_D);

% Hackler and Wood (1989) garnet mixing formulation.
numHW_inside = ...
    2080 .* (XMg_grt .^ 2) ...
    - 6350 .* (XFe_grt .^ 2) ...
    - 13807 .* XCa_grt .* (1 - XMn_grt) ...
    + 8540 .* XFe_grt .* XMg_grt .* ...
      (1 - XCa_grt - XMn_grt) ...
    + 4441 .* (2 .* XMg_mica - 1);

numHW = 20286 + 0.0193 .* P_bar - numHW_inside;
denHW = 13.138 + 8.3143 .* lnK_D + ...
    6.276 .* XCa_grt .* (1 - XMn_grt);

THW_K = numHW ./ denHW;
THW_C = THW_K - 273.15;

% Ganguly and Saxena (1984) garnet mixing formulation.
numGS_inside = ...
    837 .* (XMg_grt .^ 2) ...
    - 10460 .* (XFe_grt .^ 2) ...
    - 13807 .* XCa_grt .* (1 - XMn_grt) ...
    + 19246 .* XFe_grt .* XMg_grt .* (1 - XMn_grt) ...
    + 5649 .* XCa_grt .* (XMg_grt - XFe_grt) ...
    + 7972 .* (2 .* XMg_mica - 1);

numGS = 13538 + 0.0193 .* P_bar - numGS_inside;
denGS = 6.778 + 8.3143 .* lnK_D + ...
    6.276 .* XCa_grt .* (1 - XMn_grt);

TGS_K = numGS ./ denGS;
TGS_C = TGS_K - 273.15;

% Pack all columns at once. Each variable has nP rows, so fixed-pressure and
% pressure-range calculations return the same stable set of table variables.
row = table( ...
    P_kbar, P_bar, ...
    Fe2_grt, Fe3_grt, FeUsed_grt, Mg_grt, Mn_grt, Ca_grt, Al_grt, Si_grt, ...
    Fe2_mica, Fe3_mica, FeUsed_mica, Mg_mica, Mn_mica, Ca_mica, ...
    Ti_mica, Al_mica, Si_mica, K_mica, Na_mica, ...
    XFe_grt, XMg_grt, XCa_grt, XMn_grt, XMg_mica, XFe_mica, ...
    FeMg_grt, FeMg_mica, K_D, lnK_D, ...
    numHW, denHW, THW_K, THW_C, ...
    numGS, denGS, TGS_K, TGS_C, ...
    'VariableNames', { ...
    'P_kbar', 'P_bar', ...
    'Fe2_grt', 'Fe3_grt', 'FeUsed_grt', 'Mg_grt', 'Mn_grt', 'Ca_grt', ...
    'Al_grt', 'Si_grt', ...
    'Fe2_mica', 'Fe3_mica', 'FeUsed_mica', 'Mg_mica', 'Mn_mica', ...
    'Ca_mica', 'Ti_mica', 'Al_mica', 'Si_mica', 'K_mica', 'Na_mica', ...
    'XFe_grt', 'XMg_grt', 'XCa_grt', 'XMn_grt', 'XMg_mica', 'XFe_mica', ...
    'FeMg_grt', 'FeMg_mica', 'K_D', 'lnK_D', ...
    'numHW', 'denHW', 'THW_K', 'THW_C', ...
    'numGS', 'denGS', 'TGS_K', 'TGS_C'});

end

function mineral = prepareMineralRow(data_table, mineralLabel)
% prepareMineralRow
% Extract scalar cation values from a one-row mineral table. A missing
% optional variable is assigned zero; an existing NaN value remains NaN.

if height(data_table) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();

mineral.Fe2 = getVarOrError( ...
    data_table, 'Fe_cation_apfu', mineralLabel);
mineral.Mg = getVarOrError( ...
    data_table, 'Mg_cation_apfu', mineralLabel);

mineral.Fe3 = getVarOrZero(data_table, 'Fe3_cation_apfu', mineralLabel);
mineral.Mn = getVarOrZero(data_table, 'Mn_cation_apfu', mineralLabel);
mineral.Ca = getVarOrZero(data_table, 'Ca_cation_apfu', mineralLabel);
mineral.Ti = getVarOrZero(data_table, 'Ti_cation_apfu', mineralLabel);
mineral.Al = getVarOrZero(data_table, 'Al_cation_apfu', mineralLabel);
mineral.Si = getVarOrZero(data_table, 'Si_cation_apfu', mineralLabel);
mineral.K = getVarOrZero(data_table, 'K_cation_apfu', mineralLabel);
mineral.Na = getVarOrZero(data_table, 'Na_cation_apfu', mineralLabel);

end

function value = getVarOrError(data_table, variableName, mineralLabel)
% getVarOrError
% Read a required scalar variable. NaN and zero are allowed, negative finite
% values and Inf are rejected.

if ~ismember(variableName, data_table.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = data_table.(variableName);

if ~isnumeric(value) || ~isscalar(value) || isinf(value) || ...
        (isfinite(value) && value < 0)
    error(['%s.%s must be a finite non-negative numeric scalar or NaN.'], ...
        mineralLabel, variableName);
end

end

function value = getVarOrZero(data_table, variableName, mineralLabel)
% getVarOrZero
% Read an optional scalar variable. Missing variables become zero. Existing
% NaN values remain NaN; negative finite values and Inf are rejected.

if ismember(variableName, data_table.Properties.VariableNames)
    value = data_table.(variableName);

    if ~isnumeric(value) || ~isscalar(value) || isinf(value) || ...
            (isfinite(value) && value < 0)
        error(['%s.%s must be a finite non-negative numeric scalar or NaN.'], ...
            mineralLabel, variableName);
    end
else
    value = 0;
end

end

function printTemperatureCalibrationWarning( ...
        temperature_degC, formulationName, calibrationMin_degC, ...
        calibrationMax_degC, selectedCode_grt, selectedCode_mica)
% printTemperatureCalibrationWarning
% Print a non-stopping fprintf message when a finite calculated temperature
% lies outside the experimental temperature envelope.

finiteTemperature = isfinite(temperature_degC);
outsideCalibration = finiteTemperature & ...
    (temperature_degC < calibrationMin_degC | ...
     temperature_degC > calibrationMax_degC);

if any(outsideCalibration)
    finiteValues = temperature_degC(finiteTemperature);
    fprintf(2, ...
        ['WARNING: Calculated %s temperature is outside the experimental ' ...
         'temperature envelope of Bhattacharya et al. (1992): 575-950 degreeC. ' ...
         '%d of %d finite temperature point(s) are outside the envelope; ' ...
         'calculated finite range = %.4g-%.4g degreeC for %s & %s.\n'], ...
        formulationName, ...
        sum(outsideCalibration), ...
        sum(finiteTemperature), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(string(selectedCode_grt)), ...
        char(string(selectedCode_mica)));
end

end
