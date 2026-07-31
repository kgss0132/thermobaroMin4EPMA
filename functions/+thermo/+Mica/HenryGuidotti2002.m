function results = HenryGuidotti2002(rawdata_struct, P_kbar)
% functions/+thermo/+Mica/HenryGuidotti2002.m
% Tested with MATLAB R2024b
%
% Ti-in-biotite thermometer
% Henry, D.J. and Guidotti, C.V. (2002)
% American Mineralogist, 87, 375–382
% DOI: https://doi.org/10.2138/am-2002-0401
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Biotite analysis from
% rawdata_struct.Mica and calculates temperature using the empirical
% Ti-saturation surface of Henry and Guidotti (2002).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Biotite analysis. Each calculation is
% stored temporarily as one table block, and all blocks are concatenated only
% once after the interactive loop has finished.
%
% Both startThermoCalc_fixedP and startThermoCalc_rangeP are supported.
% P_kbar may be a finite non-negative scalar or vector. Pressure does not
% appear explicitly in the published surface-fit equation, so the calculated
% temperature is repeated for every supplied pressure value while the
% corresponding P_kbar values are retained in the output table.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION CONDITIONS AND APPLICATION NOTES
%
% Henry and Guidotti (2002) developed an isobaric Ti-saturation surface from
% natural biotites in graphitic, peraluminous metapelites from western Maine.
% The principal calibration conditions are:
%
%   Reference P : approximately 3.3 kbar
%   Principal T : approximately 500–690 degreeC
%   Extended T  : approximately 490–780 degreeC when supplementary data are
%                 included; the highest-temperature supplementary data were
%                 equilibrated at approximately 6 kbar
%   XMg         : the calibration dataset does not extend below 0.30
%   Ti          : geothermometric inference should not be made above
%                 0.50 apfu on a 22 O basis
%   Rock type   : graphitic, peraluminous metapelites
%   Assemblage  : quartz + an aluminous phase (chlorite, staurolite, or
%                 sillimanite) + ilmenite and/or rutile + graphite
%
% The approximately 3.3 kbar metamorphic framework and equilibrium criteria
% are described on pp. 375–378. The surface-fit equation and coefficients
% are given on p. 378 (Eq. 1 and Table 4). Figures 4 and 5 on pp. 378–380
% show that the calibration dataset does not include XMg < 0.30 or Ti >
% 0.50 apfu, and the text states that geothermometric inferences should not
% be made above Ti = 0.50 apfu. Chemical equilibrium, re-equilibration,
% alteration, and inclusion-thermometry limitations are discussed on
% pp. 380–381.
%
% IMPORTANT APPLICATION LIMITATIONS
% - In the strictest sense, the published surface is an approximately
%   3.3 kbar isobaric Ti-saturation surface. The paper does not define a
%   finite pressure interval around 3.3 kbar. Consequently, this
%   implementation flags every input pressure that is not numerically equal
%   to 3.3 kbar. This is a caution, not a stopping condition.
% - Pressure affects Ti solubility even though pressure is absent from the
%   fitted equation. Increasing pressure generally lowers Ti solubility in
%   biotite (discussion on pp. 375–376 and 380–381).
% - The principal western Maine data span approximately 500–690 degreeC at
%   approximately 3.3 kbar. Supplementary datasets extend the empirical
%   temperature envelope to approximately 490–780 degreeC, but the
%   highest-temperature supplementary data are from approximately 6 kbar
%   and may slightly underestimate the 3.3 kbar Ti-saturation level (p. 378).
% - The thermometer should be restricted to graphitic, peraluminous
%   metapelites with the assemblage listed above. Other bulk compositions or
%   assemblages may not satisfy the required Si-, Al-, and Ti-saturation
%   conditions.
% - Chemical equilibrium must be evaluated independently. Green,
%   chloritized, compositionally heterogeneous, locally re-equilibrated, or
%   thermally overprinted biotite may yield anomalously low or otherwise
%   misleading temperatures (pp. 380–381).
% - For inclusion thermometry, the inclusion must have satisfied the same
%   assemblage criteria at entrapment. Mg-Fe exchange with a mafic host such
%   as garnet can modify XMg during cooling (p. 381).
% - Ti and the biotite formula must be normalized on a 22 O basis. The
%   original dataset contained relatively low and approximately constant
%   Fe3+ (~12% of total Fe). Fe_cation_apfu should therefore represent the Fe
%   value used in the 22 O-normalized formula consistently with the original
%   calibration. A separately supplied Fe3_cation_apfu value is retained for
%   reporting but is not added to Fe_cation_apfu in this implementation.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure differs from the approximately 3.3 kbar reference,
%   2) finite XMg is below 0.30 or above 1.00,
%   3) finite Ti is <= 0 or above 0.50 apfu,
%   4) a finite calculated temperature lies outside the extended empirical
%      envelope of approximately 490–780 degreeC,
%   5) a finite calculated temperature lies within 490–780 degreeC but
%      outside the principal approximately 500–690 degreeC dataset,
%   6) a required thermometer input is NaN,
%   7) the cubic argument is non-positive, or
%   8) the calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Mica : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must contain
% normalized cation data.
%
% Required variables used in the thermometer:
%   Mg_cation_apfu
%   Fe_cation_apfu
%   Ti_cation_apfu       % Ti apfu, normalized to 22 O
%
% Optional variables retained in the output:
%   Al_cation_apfu
%   Si_cation_apfu
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ca_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Missing optional variables are stored as NaN, not zero. NaN values in the
% required variables are retained, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Finite negative values are not
% allowed in any mineral-composition variable read by this function. Zero is
% retained; if it makes the equation undefined, the resulting NaN/Inf is kept
% and reported rather than replaced.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Henry and Guidotti (2002) fitted the Ti-saturation surface as:
%
%   ln(z) = a + b*x^3 + c*y^3
%
% where:
%   z = Ti (apfu normalized to 22 O)
%   x = T (degreeC)
%   y = XMg = Mg / (Mg + Fe)
%
% Rearranged for temperature:
%
%   T(degreeC) = {[ln(Ti) - a - c*(XMg)^3] / b}^(1/3)
%
% with:
%   a = -2.3353
%   b = 4.3430e-9
%   c = -1.6718
%
% Only positive, finite cubic arguments are treated as valid temperature
% solutions. A non-positive cubic argument is retained in the output, but
% the calculated temperature is set to NaN. A negative real cube root is not
% interpreted as a physically meaningful temperature.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = HenryGuidotti2002(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing a Mica table (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Biotite analysis
%

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid launcher pressure values.
if nargin < 2
    error('HenryGuidotti2002 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation dataset
% Extract the required Mica table. The table itself is not modified.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_bt = rawdata_struct.Mica;

if height(dataset_bt) == 0
    error('rawdata_struct.Mica must contain at least one analysis row.');
end

requiredVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', 'Ti_cation_apfu'};
missingRequiredVariables = requiredVariables( ...
    ~ismember(requiredVariables, dataset_bt.Properties.VariableNames));

if ~isempty(missingRequiredVariables)
    error(['rawdata_struct.Mica is missing required variable(s): ' ...
        char(strjoin(string(missingRequiredVariables), ', ')) '.']);
end

% Prepare display strings once instead of rebuilding them in every loop.
dataCodes_bt = dataset_bt{:, 1};
dataCodeList_bt = cellstr(string(dataCodes_bt));

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Each calculation result is stored temporarily as one table block.
% Repeated concatenation of the full results table inside the loop is avoided.
% The buffer follows the Ballhaus1991 implementation and doubles only when
% its current capacity is exhausted.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Empirical limits and reference conditions from Henry and Guidotti (2002).
referenceP_kbar = 3.3;
referencePressureTolerance_kbar = 1e-9;
principalT_min_degC = 500;
principalT_max_degC = 690;
extendedT_min_degC = 490;
extendedT_max_degC = 780;
plotXMg_min = 0.30;
plotXMg_max = 1.00;
plotTi_min_apfu = 0.00;
plotTi_max_apfu = 0.50;

% Pressure is common to every selected analysis in this function call, so
% the reference-pressure warning is printed only once.
pressureDiffersFromReference = ...
    abs(P_kbar - referenceP_kbar) > referencePressureTolerance_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
% The loop continues until the user cancels the dialog or selects Finish.
disp('=== Step 3: Selecting a data code from the list (Biotite) ===');

while true
    % ----- Biotite selection -----
    [selectedIdx_bt, ok] = listdlg( ...
        'PromptString', 'Please select the Biotite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodeList_bt, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_bt)
        disp('Selection canceled');
        break;
    end

    selectedCode_bt = dataCodes_bt(selectedIdx_bt);
    disp(['Biotite selected: ' char(string(selectedCode_bt))]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    selectedData_bt = dataset_bt(selectedIdx_bt, :);

    % Identify NaN values in the variables actually used by the thermometer.
    % NaN does not stop the calculation.
    nanInputNames = findNaNInputs(selectedData_bt);

    % Reject finite negative values. Zero and NaN are retained so that their
    % mathematical consequences remain visible and diagnosable.
    validateNonNegativeInputs(selectedData_bt);

    row = calcTemp(selectedData_bt, P_kbar);

    % Store the selected identifier once for each pressure row.
    row.dataCode_bt = repmat(string(selectedCode_bt), height(row), 1);
    row = movevars(row, {'dataCode_bt'}, 'Before', 1);

    % Store this calculation as one buffered table block.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_bt)) ...
            ': HenryGuidotti2002 = ' ...
            num2str(row.THenryGuidotti2002_C) ' degreeC']);
    else
        disp([char(string(selectedCode_bt)) ...
            ': HenryGuidotti2002 = ' ...
            num2str(row.THenryGuidotti2002_C(1)) ' to ' ...
            num2str(row.THenryGuidotti2002_C(end)) ' degreeC']);
    end

    % Warn once if any supplied pressure differs from the approximately
    % 3.3 kbar reference pressure. The paper does not define a finite pressure
    % interval, so no broader tolerance is invented here.
    if any(pressureDiffersFromReference) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure differs from the approximately 3.3 kbar ' ...
             'reference pressure of the Henry and Guidotti (2002) isobaric ' ...
             'Ti-saturation surface. %d of %d pressure point(s) differ from ' ...
             '3.3 kbar; input range = %.4g–%.4g kbar. The calculation was ' ...
             'continued, but pressure-dependent bias may occur.\n'], ...
            sum(pressureDiffersFromReference), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when finite XMg lies outside the compositional region represented
    % by the published dataset.
    finiteXMg = isfinite(row.XMg);
    XMgOutsidePlotRange = finiteXMg & ...
        (row.XMg < plotXMg_min | row.XMg > plotXMg_max);

    if any(XMgOutsidePlotRange)
        finiteValues = row.XMg(finiteXMg);
        fprintf(2, ...
            ['WARNING: XMg is outside the compositional region represented by ' ...
             'Henry and Guidotti (2002): the calibration dataset does not ' ...
             'extend below XMg = 0.30. %d of %d finite point(s) are outside ' ...
             '0.30–1.00; finite XMg range = %.4g–%.4g for %s.\n'], ...
            sum(XMgOutsidePlotRange), ...
            sum(finiteXMg), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_bt)));
    end

    % Warn when finite Ti is non-positive or exceeds the published upper
    % limit for geothermometric inference.
    finiteTi = isfinite(row.Ti_apfu_22O);
    TiOutsidePlotRange = finiteTi & ...
        (row.Ti_apfu_22O <= plotTi_min_apfu | ...
         row.Ti_apfu_22O > plotTi_max_apfu);

    if any(TiOutsidePlotRange)
        finiteValues = row.Ti_apfu_22O(finiteTi);
        fprintf(2, ...
            ['WARNING: Ti in biotite is outside the usable range represented ' ...
             'by Henry and Guidotti (2002): Ti must be > 0 for ln(Ti), and ' ...
             'geothermometric inference should not be made above 0.50 apfu ' ...
             'on a 22 O basis. %d of %d finite point(s) are outside the ' ...
             'usable range; finite Ti range = %.4g–%.4g apfu for %s.\n'], ...
            sum(TiOutsidePlotRange), ...
            sum(finiteTi), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_bt)));
    end

    % Warn when finite calculated temperatures fall outside the full empirical
    % envelope represented by the principal and supplementary datasets.
    finiteTemperature = isfinite(row.THenryGuidotti2002_C);
    temperatureOutsideExtendedRange = finiteTemperature & ...
        (row.THenryGuidotti2002_C < extendedT_min_degC | ...
         row.THenryGuidotti2002_C > extendedT_max_degC);

    if any(temperatureOutsideExtendedRange)
        finiteValues = row.THenryGuidotti2002_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximate ' ...
             'extended empirical temperature envelope of Henry and Guidotti ' ...
             '(2002): 490–780 degreeC. %d of %d finite temperature point(s) ' ...
             'are outside the envelope; calculated finite range = %.4g–%.4g ' ...
             'degreeC for %s.\n'], ...
            sum(temperatureOutsideExtendedRange), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_bt)));
    end

    % Distinguish temperatures that are within the extended envelope but
    % outside the principal approximately 3.3 kbar western Maine dataset.
    temperatureOutsidePrincipalRange = finiteTemperature & ...
        ~temperatureOutsideExtendedRange & ...
        (row.THenryGuidotti2002_C < principalT_min_degC | ...
         row.THenryGuidotti2002_C > principalT_max_degC);

    if any(temperatureOutsidePrincipalRange)
        finiteValues = row.THenryGuidotti2002_C( ...
            temperatureOutsidePrincipalRange);
        fprintf(2, ...
            ['WARNING: Calculated temperature lies outside the principal ' ...
             'approximately 500–690 degreeC western Maine dataset but within ' ...
             'the broader 490–780 degreeC envelope that includes supplementary ' ...
             'data. %d point(s); affected finite range = %.4g–%.4g degreeC ' ...
             'for %s. High-temperature supplementary data include samples at ' ...
             'approximately 6 kbar.\n'], ...
            sum(temperatureOutsidePrincipalRange), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_bt)));
    end

    % Display the exact required input variables that contained NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'temperature may therefore be NaN.\n'], ...
            char(string(selectedCode_bt)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % A non-positive finite cubic argument has no valid positive-temperature
    % solution in this implementation. Retain the argument and keep T as NaN.
    invalidCubicArgument = isfinite(row.cubic_argument) & ...
        row.cubic_argument <= 0;

    if any(invalidCubicArgument)
        finiteArguments = row.cubic_argument(invalidCubicArgument);
        fprintf(2, ...
            ['WARNING: The Henry and Guidotti (2002) equation produced a ' ...
             'non-positive cubic argument for %s (%d of %d points; range = ' ...
             '%.4g–%.4g). No valid positive-temperature solution was assigned, ' ...
             'and the corresponding temperature remains NaN.\n'], ...
            char(string(selectedCode_bt)), ...
            sum(invalidCubicArgument), ...
            numel(row.cubic_argument), ...
            min(finiteArguments), ...
            max(finiteArguments));
    end

    % Retain and report any NaN/Inf result instead of replacing it or stopping.
    invalidTemperature = ~isfinite(row.THenryGuidotti2002_C);

    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_bt)), ...
            sum(invalidTemperature), ...
            numel(row.THenryGuidotti2002_C), ...
            sum(isnan(row.THenryGuidotti2002_C)), ...
            sum(isinf(row.THenryGuidotti2002_C)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another biotite analysis.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'HenryGuidotti2002', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after all selections finish.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_biotite)
% findNaNInputs
% Return names of required thermometer variables containing NaN. The output
% is selected from a fixed-size logical mask and does not grow in the loop.

requiredVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', 'Ti_cation_apfu'};
qualifiedNames = "Mica." + string(requiredVariables(:));
containsNaN = false(numel(requiredVariables), 1);

for i = 1:numel(requiredVariables)
    variableValue = data_biotite.(requiredVariables{i});
    containsNaN(i) = any(isnan(variableValue(:)));
end

nanInputNames = qualifiedNames(containsNaN);

end

function validateNonNegativeInputs(data_biotite)
% validateNonNegativeInputs
% Stop calculation when a supplied mineral-composition value is negative,
% nonnumeric, nonscalar, or Inf. NaN and zero are intentionally allowed.

variablesToCheck = { ...
    'Mg_cation_apfu', 'Fe_cation_apfu', 'Ti_cation_apfu', ...
    'Al_cation_apfu', 'Si_cation_apfu', 'Fe3_cation_apfu', ...
    'Mn_cation_apfu', 'Ca_cation_apfu', 'K_cation_apfu', ...
    'Na_cation_apfu'};

variableExists = ismember(variablesToCheck, ...
    data_biotite.Properties.VariableNames);
invalidShapeOrType = false(numel(variablesToCheck), 1);
containsInf = false(numel(variablesToCheck), 1);
containsNegative = false(numel(variablesToCheck), 1);

for i = 1:numel(variablesToCheck)
    if ~variableExists(i)
        continue;
    end

    variableValue = data_biotite.(variablesToCheck{i});
    invalidShapeOrType(i) = ~isnumeric(variableValue) || ~isscalar(variableValue);

    if ~invalidShapeOrType(i)
        containsInf(i) = isinf(variableValue);
        containsNegative(i) = isfinite(variableValue) && variableValue < 0;
    end
end

if any(invalidShapeOrType)
    invalidNames = string(variablesToCheck(invalidShapeOrType));
    error(['HenryGuidotti2002: mineral-composition variables must ' ...
        'contain one numeric scalar per selected row. Invalid variable(s): ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

if any(containsInf)
    invalidNames = string(variablesToCheck(containsInf));
    error(['HenryGuidotti2002: Inf is not permitted in ' ...
        'mineral-composition variables. Invalid variable(s): ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

if any(containsNegative)
    invalidNames = string(variablesToCheck(containsNegative));
    error(['HenryGuidotti2002: mineral-composition values must be >= 0 ' ...
        'or NaN. Negative finite value(s) were found in: ' ...
        char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_biotite, P_kbar)
% calcTemp
% Calculate the Henry and Guidotti (2002) Ti-in-biotite temperature for one
% selected biotite row and every supplied pressure value.
%
% Pressure is retained for launcher compatibility and calibration screening,
% but it is not used explicitly in the published surface-fit equation.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

row = table();

% --- Henry and Guidotti (2002) empirical constants ---
a_const = -2.3353;
b_const = 4.3430e-9;
c_const = -1.6718;
referenceP_kbar = 3.3;

% --- Extract one-row mineral data ---
bt = prepareMineralRow(data_biotite, 'Mica');

% Replicate pressure-independent composition values to match the supplied
% pressure vector. NaN remains NaN and is never converted to zero.
Mg_bt = repmat(bt.Mg, nP, 1);
Fe_bt = repmat(bt.Fe, nP, 1);
Ti_bt = repmat(bt.Ti, nP, 1);

Al_bt = repmat(bt.Al, nP, 1);
Si_bt = repmat(bt.Si, nP, 1);
Fe3_bt = repmat(bt.Fe3, nP, 1);
Mn_bt = repmat(bt.Mn, nP, 1);
Ca_bt = repmat(bt.Ca, nP, 1);
K_bt = repmat(bt.K, nP, 1);
Na_bt = repmat(bt.Na, nP, 1);

% --- Core thermometer parameters ---
XMg = Mg_bt ./ (Mg_bt + Fe_bt);
Ti_apfu_22O = Ti_bt;

lnTi = log(Ti_apfu_22O);
cubic_argument = ...
    (lnTi - a_const - c_const .* (XMg .^ 3)) ./ b_const;

% Only a positive finite cubic argument is treated as a valid temperature
% solution. This prevents a negative real cube root from being reported as a
% physically meaningful negative temperature.
THenryGuidotti2002_C = NaN(nP, 1);
validCubicArgument = isfinite(cubic_argument) & cubic_argument > 0;
THenryGuidotti2002_C(validCubicArgument) = ...
    cubic_argument(validCubicArgument) .^ (1 / 3);

T_K = THenryGuidotti2002_C + 273.15;

% --- Calibration and reference-condition flags ---
is_reference_pressure_match = abs(P_kbar - referenceP_kbar) <= 1e-9;
is_XMg_in_plot_range = isfinite(XMg) & XMg >= 0.30 & XMg <= 1.00;
is_Ti_in_plot_range = isfinite(Ti_apfu_22O) & ...
    Ti_apfu_22O > 0 & Ti_apfu_22O <= 0.50;
is_T_in_principal_range = isfinite(THenryGuidotti2002_C) & ...
    THenryGuidotti2002_C >= 500 & THenryGuidotti2002_C <= 690;
is_T_in_extended_range = isfinite(THenryGuidotti2002_C) & ...
    THenryGuidotti2002_C >= 490 & THenryGuidotti2002_C <= 780;

% --- Pack outputs ---
row.P_kbar = P_kbar;
row.P_reference_kbar = repmat(referenceP_kbar, nP, 1);
row.is_reference_pressure_match = is_reference_pressure_match;

row.Mg_bt = Mg_bt;
row.Fe_bt = Fe_bt;
row.Fe2_bt = Fe_bt; % Backward-compatible alias for the original script.
row.Fe3_bt = Fe3_bt;
row.Al_bt = Al_bt;
row.Si_bt = Si_bt;
row.Mn_bt = Mn_bt;
row.Ca_bt = Ca_bt;
row.Ti_bt = Ti_bt;
row.K_bt = K_bt;
row.Na_bt = Na_bt;

row.XMg = XMg;
row.Ti_apfu_22O = Ti_apfu_22O;
row.lnTi = lnTi;
row.cubic_argument = cubic_argument;

row.T_K = T_K;
row.T_deg = THenryGuidotti2002_C;
row.THenryGuidotti2002_C = THenryGuidotti2002_C;

row.is_XMg_in_plot_range = is_XMg_in_plot_range;
row.is_Ti_in_plot_range = is_Ti_in_plot_range;
row.is_T_in_principal_range = is_T_in_principal_range;
row.is_T_in_extended_range = is_T_in_extended_range;

end

function mineral = prepareMineralRow(data_tbl, mineralLabel)
% prepareMineralRow
% Extract cation values from one selected biotite row. Required and optional
% NaN values are retained. Missing optional variables are represented by NaN.

if height(data_tbl) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();

mineral.Mg = getRequiredVar(data_tbl, 'Mg_cation_apfu', mineralLabel);
mineral.Fe = getRequiredVar(data_tbl, 'Fe_cation_apfu', mineralLabel);
mineral.Ti = getRequiredVar(data_tbl, 'Ti_cation_apfu', mineralLabel);

mineral.Al = getOptionalVar(data_tbl, 'Al_cation_apfu', mineralLabel);
mineral.Si = getOptionalVar(data_tbl, 'Si_cation_apfu', mineralLabel);
mineral.Fe3 = getOptionalVar(data_tbl, 'Fe3_cation_apfu', mineralLabel);
mineral.Mn = getOptionalVar(data_tbl, 'Mn_cation_apfu', mineralLabel);
mineral.Ca = getOptionalVar(data_tbl, 'Ca_cation_apfu', mineralLabel);
mineral.K = getOptionalVar(data_tbl, 'K_cation_apfu', mineralLabel);
mineral.Na = getOptionalVar(data_tbl, 'Na_cation_apfu', mineralLabel);

end

function value = getRequiredVar(tbl, varName, mineralLabel)
% getRequiredVar
% Read one required scalar value. NaN is allowed and retained; Inf and
% negative finite values are rejected.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must contain one numeric scalar.', mineralLabel, varName);
end
if isinf(value)
    error('%s.%s must not contain Inf.', mineralLabel, varName);
end
if isfinite(value) && value < 0
    error('%s.%s must be >= 0 or NaN.', mineralLabel, varName);
end

end

function value = getOptionalVar(tbl, varName, mineralLabel)
% getOptionalVar
% Read one optional scalar value. Missing optional variables are represented
% by NaN so missing information is not silently converted to zero.

if ~ismember(varName, tbl.Properties.VariableNames)
    value = NaN;
    return;
end

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must contain one numeric scalar.', mineralLabel, varName);
end
if isinf(value)
    error('%s.%s must not contain Inf.', mineralLabel, varName);
end
if isfinite(value) && value < 0
    error('%s.%s must be >= 0 or NaN.', mineralLabel, varName);
end

end
