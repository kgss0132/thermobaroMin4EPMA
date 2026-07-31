function results = Putirka2008olbaro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Olivine/Putirka2008olbaro.m
% Tested with MATLAB R2024b
%
% Silica-activity barometer for Olivine-Orthopyroxene-saturated liquids
% Putirka, K.D. (2008), Equation (42)
% Reviews in Mineralogy and Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Liquid analysis and calculates
% pressure using the silica-activity barometer of Putirka (2008), Equation
% (42). The equation uses liquid composition and independently supplied
% temperature; no Olivine or Orthopyroxene composition is used numerically.
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Liquid row, one output row is
% returned for every input temperature value.
%
% rawdata_struct is accepted for compatibility with the common barometer
% launcher interface. Liquid composition is read through
% liquid.readLiquidExcel(), following the original implementation.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% IMPORTANT: Equation (42) is calibrated only for liquids that were in
% equilibrium with BOTH olivine and orthopyroxene, with or without additional
% phases. The relevant silica-buffering reaction and barometer development
% are discussed on pp. 96-98, especially Equations (40)-(42) on pp. 96-97
% and Figure 11 on p. 98. The presence of olivine alone is not sufficient.
% This function cannot test Ol-Opx-Liquid saturation or equilibrium because
% it uses liquid composition only.
%
% The global calibration dataset spans (p. 97):
%
%   Temperature : 825-2000 degreeC
%   Pressure    : 0.001-70 kbar
%   SiO2_liq    : 31.5-70 wt%
%   Dataset     : 510 experiments; hydrous compositions are included
%
% Equation (42) reproduces pressure with (pp. 97-98; Fig. 11):
%
%   R^2                  : 0.91
%   Standard error (SEE) : +/-2.87 kbar (approximately +/-2.9 kbar)
%
% The calibration used a global partial-melting database and iteratively
% removed data outside +/-3 standard deviations. Putirka (2008) notes that
% smaller-data-set calibrations provided little improvement in calibration
% error and performed poorly for independent data (p. 97).
%
% Silica activity alone explains approximately 80% of pressure variation for
% experiments with liquid Mg# > 0.75, but Equation (42) was calibrated using
% all 510 experiments without an Mg# restriction (p. 97). Mg# > 0.75 is
% therefore not imposed as a hard applicability limit here.
%
% Liquid cation fractions must be calculated exactly as defined by Putirka
% (2008): oxide wt% values are converted to cation proportions and normalized
% on an anhydrous basis (pp. 68-71; Table 1). Input oxide wt% values must not
% first be renormalized to 100 wt%. H2O is excluded from the cation-fraction
% denominator. In this implementation, F and Cl are also excluded from
% cationTotal_liq and from NaN-input warnings.
%
% The liquid composition must represent the liquid that equilibrated with
% olivine and orthopyroxene. A whole-rock composition is an appropriate liquid
% proxy only when petrologically justified. Crystal accumulation, fractional
% crystallization, magma mixing, alteration, or disequilibrium can yield
% misleading pressures; these general treatment errors are discussed on
% pp. 107-108.
%
% The reported +/-2.87 kbar SEE is calibration-model error. It does not
% include uncertainty in the supplied temperature, liquid analysis,
% phase-equilibrium interpretation, or choice of whole rock versus glass.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 825-2000 degreeC,
%   2) finite calculated pressure is outside 0.001-70 kbar,
%   3) finite SiO2_liq is outside 31.5-70 wt%,
%   4) a calculation input contains NaN, or
%   5) a calculated pressure is NaN or Inf.
%
% Finite input values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if zero makes an inverse power,
% logarithm, or ratio undefined, the resulting NaN or Inf is retained and
% reported.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% Liquid composition is read through liquid.readLiquidExcel(). The selected
% Liquid table should contain wt% oxide columns for the following values used
% in cationTotal_liq:
%
%   SiO2, TiO2, Al2O3, FeO or FeOt, MnO, MgO, CaO, Na2O, K2O,
%   V2O3, Cr2O3, NiO, P2O5, SO3, and Fe2O3
%
% Missing or non-numeric calculation inputs are represented by NaN and are
% never replaced by zero. Because the cation fractions are normalized to
% cationTotal_liq, NaN in any included oxide propagates through the activity
% calculation and pressure result.
%
% F and Cl may be retained in the output for traceability, but they are not
% included in cationTotal_liq, are not validated as calculation inputs, and
% are excluded from NaN-input warnings.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (Putirka, 2008, Eq. 42; p. 97)
%
%   P(kbar) =
%       231.5
%       + 0.186*T(degreeC)
%       + 0.1244*T(degreeC)*ln(aSiO2_liq)
%       - 528.5*(aSiO2_liq)^0.5
%       + 103.3*XTiO2_liq
%       + 69.9*(XNaO0.5_liq + XKO0.5_liq)
%       + 77.3*[XAlO1.5_liq/(XAlO1.5_liq + XSiO2_liq)]
%
% where the Beattie (1993) silica activity is (Eq. 41; p. 97):
%
%   aSiO2_liq =
%       (3*XSiO2_liq)^(-2)
%       *(1 - XAlO1.5_liq)^(7/2)
%       *(1 - XTiO2_liq)^7
%
% All X terms are liquid cation fractions (Table 1, pp. 69-70).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008olbaro(rawdata_struct, T_degreeC)
%   results = Putirka2008olbaro(rawdata_struct, T_degreeC, ...
%       'LiquidRow', rowNumber)
%
% Inputs:
%   rawdata_struct : struct accepted for barometer-interface compatibility
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Optional name-value input:
%   LiquidRow      : positive integer row number in the selected Liquid
%                    table. If omitted, the Liquid row is selected using a
%                    list dialog and repeated selections are permitted.
%
% Output:
%   results : table containing one row per temperature value for every
%             selected Liquid row.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('Putirka2008olbaro requires (rawdata_struct, T_degreeC).');
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

ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
    x >= 1 && x == fix(x)));
ip.parse(varargin{:});
liquidRowOption = ip.Results.LiquidRow;

%% 1) Retrieve Liquid dataset
% Read the Liquid dataset once. Selected rows are copied without modifying
% the source table.
disp('=== Step 1: Preparing Liquid dataset ===');

MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();
if isempty(liqAll)
    error('Selected Liquid dataset is empty.');
end

if ~isempty(liquidRowOption) && liquidRowOption > height(liqAll)
    error(['Requested LiquidRow (%d) exceeds the number of rows in the ' ...
           'selected Liquid dataset (%d).'], ...
        liquidRowOption, height(liqAll));
end

liquidItems = buildLiquidList(liqAll);
disp('=== Preparing Liquid dataset has been finished ===');

%% 2) Initialize output container and calibration limits
% Store each selected-Liquid result in a cell buffer and concatenate once
% after the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max(16, height(liqAll));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degreeC = 825;
calibrationT_max_degreeC = 2000;
calibrationP_min_kbar = 0.001;
calibrationP_max_kbar = 70;
calibrationSiO2_min_wtpercent = 31.5;
calibrationSiO2_max_wtpercent = 70;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);
temperatureWarningIssued = false;
phaseAssemblageCautionIssued = false;
fixedLiquidRowProcessed = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive Liquid selection + calculation
% When LiquidRow is omitted, the loop continues until the user cancels the
% selection dialog or chooses Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Liquid) ===');

while true
    % ----- Liquid selection -----
    if isempty(liquidRowOption)
        [selectedIdx_liq, ok] = listdlg( ...
            'PromptString', 'Please select the Liquid data you would like to use:', ...
            'SelectionMode', 'single', ...
            'ListString', liquidItems, ...
            'ListSize', [520 360]);

        if ~ok || isempty(selectedIdx_liq)
            disp('Selection canceled');
            break;
        end
    else
        if fixedLiquidRowProcessed
            break;
        end
        selectedIdx_liq = liquidRowOption;
        fixedLiquidRowProcessed = true;
    end

    selectedData_liq = liqAll(selectedIdx_liq, :);
    selectedLiquidLabel = getLiquidLabel(selectedData_liq, selectedIdx_liq);
    disp(['Liquid selected: ' selectedLiquidLabel]);

    % ----- Calculation -----
    disp('=== Step 4: Checking inputs and calculating pressure ===');

    liquidComposition = prepareLiquidComposition(selectedData_liq);

    % Check NaN only in variables used by Equation (42) or its cation-total
    % denominator. F and Cl are intentionally excluded.
    nanInputNames = findNaNInputs(liquidComposition, T_degreeC);

    % Reject Inf and finite negative calculation inputs. Zero and NaN remain.
    validateNonNegativeInputs(liquidComposition);

    row = calcPressure(liquidComposition, T_degreeC, MWinfo);

    % Repeat Liquid identifiers for every temperature in this calculation.
    nRows = height(row);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);

    preferredFrontVariables = ...
        {'dataRow_liq', 'liq_Index', 'liq_Experiment', 'liq_Citation'};
    frontVariables = preferredFrontVariables( ...
        ismember(preferredFrontVariables, row.Properties.VariableNames));
    if ~isempty(frontVariables)
        row = movevars(row, frontVariables, 'Before', 1);
    end

    % Store one block per selected Liquid row. The buffer expands only when
    % its current capacity is exhausted, not on every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([selectedLiquidLabel ': P = ' num2str(row.P_kbar) ' kbar']);
    else
        disp([selectedLiquidLabel ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Print the phase-assemblage requirement once per function call.
    if ~phaseAssemblageCautionIssued
        fprintf(2, ...
            ['CAUTION: Putirka (2008) Equation (42) was calibrated for liquids ' ...
             'in equilibrium with both olivine and orthopyroxene (+/- other ' ...
             'phases; pp. 96-98). This function uses liquid composition only ' ...
             'and cannot verify that phase-equilibrium requirement.\n']);
        phaseAssemblageCautionIssued = true;
    end

    % Input temperature is common to all selected Liquid rows, so print this
    % warning only once after the first completed calculation.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the 825-2000 degreeC ' ...
             'calibration range reported for Putirka (2008) Equation (42) ' ...
             '(p. 97). %d of %d finite temperature point(s) are outside the ' ...
             'range; input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures fall outside the calibration
    % pressure range. Calculated values are retained unchanged.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the 0.001-70 kbar ' ...
             'calibration range reported for Putirka (2008) Equation (42) ' ...
             '(p. 97). %d of %d finite pressure point(s) are outside the ' ...
             'range; calculated finite range = %.4g-%.4g kbar for %s.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            selectedLiquidLabel);
    end

    % Warn when finite SiO2 lies outside the reported calibration interval.
    SiO2_value = liquidComposition.SiO2;
    if isfinite(SiO2_value) && ...
            (SiO2_value < calibrationSiO2_min_wtpercent || ...
             SiO2_value > calibrationSiO2_max_wtpercent)
        fprintf(2, ...
            ['WARNING: Liquid SiO2 = %.4g wt%% is outside the 31.5-70 wt%% ' ...
             'range of the Equation (42) dataset reported by Putirka (2008; ' ...
             'p. 97) for %s.\n'], ...
            SiO2_value, selectedLiquidLabel);
    end

    % List exact calculation inputs containing NaN. F and Cl are absent from
    % this list and do not affect cationTotal_liq.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'pressure may remain NaN. F and Cl are excluded from both ' ...
             'cationTotal_liq and this NaN warning.\n'], ...
            selectedLiquidLabel, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            selectedLiquidLabel, ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressures are retained for diagnosis but are outside
    % the physical and calibration domains.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s ' ...
             '(%d of %d points). The values were retained for diagnostic ' ...
             'purposes.\n'], ...
            selectedLiquidLabel, ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    if ~isempty(liquidRowOption)
        break;
    end

    userAction = questdlg( ...
        'Continue with another Liquid selection?', ...
        'Putirka2008olbaro', ...
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

results.Properties.UserData = struct('liquid', metaLiq);
disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function liquidComposition = prepareLiquidComposition(data_liq)
% prepareLiquidComposition
% Extract one-row Liquid oxide data without replacing NaN by zero. FeO is
% used when available; FeOt is used only when FeO is NaN. F and Cl are read
% for traceability but excluded from calculation and validation.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liquidComposition = struct();
liquidComposition.SiO2 = getLiquidOxide(data_liq, 'SiO2');
liquidComposition.TiO2 = getLiquidOxide(data_liq, 'TiO2');
liquidComposition.Al2O3 = getLiquidOxide(data_liq, 'Al2O3');

FeO_direct = getLiquidOxide(data_liq, 'FeO');
FeOt = getLiquidOxide(data_liq, 'FeOt');
if isnan(FeO_direct)
    liquidComposition.FeO = FeOt;
    liquidComposition.FeO_source = "FeOt";
else
    liquidComposition.FeO = FeO_direct;
    liquidComposition.FeO_source = "FeO";
end

liquidComposition.MnO = getLiquidOxide(data_liq, 'MnO');
liquidComposition.MgO = getLiquidOxide(data_liq, 'MgO');
liquidComposition.CaO = getLiquidOxide(data_liq, 'CaO');
liquidComposition.Na2O = getLiquidOxide(data_liq, 'Na2O');
liquidComposition.K2O = getLiquidOxide(data_liq, 'K2O');
liquidComposition.V2O3 = getLiquidOxide(data_liq, 'V2O3');
liquidComposition.Cr2O3 = getLiquidOxide(data_liq, 'Cr2O3');
liquidComposition.NiO = getLiquidOxide(data_liq, 'NiO');
liquidComposition.P2O5 = getLiquidOxide(data_liq, 'P2O5');
liquidComposition.SO3 = getLiquidOxide(data_liq, 'SO3');
liquidComposition.Fe2O3 = getLiquidOxide(data_liq, 'Fe2O3');

% F and Cl are retained only for traceability.
liquidComposition.F = getLiquidOxide(data_liq, 'F');
liquidComposition.Cl = getLiquidOxide(data_liq, 'Cl');

end

function nanInputNames = findNaNInputs(liquidComposition, T_degreeC)
% findNaNInputs
% Return names of pressure-calculation inputs containing NaN. F and Cl are
% intentionally excluded. NaN values are not changed and do not stop the
% calculation.

calculationFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', ...
    'Fe2O3'};

maxNames = numel(calculationFields) + 1;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

for i = 1:numel(calculationFields)
    fieldName = calculationFields{i};
    fieldValue = liquidComposition.(fieldName);
    if isnan(fieldValue)
        nNanInputs = nNanInputs + 1;
        if strcmp(fieldName, 'FeO')
            nanInputBuffer(nNanInputs) = "Liq.FeO/FeOt";
        else
            nanInputBuffer(nNanInputs) = "Liq." + string(fieldName);
        end
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(liquidComposition)
% validateNonNegativeInputs
% Reject Inf and finite negative values in variables used by Equation (42)
% or cationTotal_liq. Zero and NaN are intentionally allowed and retained.
% F and Cl are excluded because they are not calculation inputs.

calculationFields = {'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', ...
    'CaO', 'Na2O', 'K2O', 'V2O3', 'Cr2O3', 'NiO', 'P2O5', 'SO3', ...
    'Fe2O3'};

invalidInputBuffer = strings(numel(calculationFields), 1);
nInvalidInputs = 0;

for i = 1:numel(calculationFields)
    fieldName = calculationFields{i};
    fieldValue = liquidComposition.(fieldName);
    if isinf(fieldValue) || (isfinite(fieldValue) && fieldValue < 0)
        nInvalidInputs = nInvalidInputs + 1;
        if strcmp(fieldName, 'FeO')
            invalidInputBuffer(nInvalidInputs) = "Liq.FeO/FeOt";
        else
            invalidInputBuffer(nInvalidInputs) = "Liq." + string(fieldName);
        end
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['Putirka2008olbaro: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(liquidComposition, T_degreeC, MWinfo)
% calcPressure
% Compute pressure for one Liquid row at one or more input temperatures.
% NaN values are never replaced and propagate naturally through
% cationTotal_liq, cation fractions, silica activity, and pressure.
%
% Inputs:
%   liquidComposition : scalar oxide-composition struct
%   T_degreeC         : scalar or vector temperature in degreeC
%   MWinfo            : molar-weight and cation-number information
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Convert wt% oxides to cation proportions. F, Cl, and H2O are omitted.
n = struct();
n.SiO2 = liquidComposition.SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
n.TiO2 = liquidComposition.TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
n.Al2O3 = liquidComposition.Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n.FeO = liquidComposition.FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
n.MnO = liquidComposition.MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
n.MgO = liquidComposition.MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
n.CaO = liquidComposition.CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
n.Na2O = liquidComposition.Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
n.K2O = liquidComposition.K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
n.V2O3 = liquidComposition.V2O3 .* MWinfo.Cat.V2O3 ./ MWinfo.MW.V2O3;
n.Cr2O3 = liquidComposition.Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;
n.NiO = liquidComposition.NiO .* MWinfo.Cat.NiO ./ MWinfo.MW.NiO;
n.P2O5 = liquidComposition.P2O5 .* MWinfo.Cat.P2O5 ./ MWinfo.MW.P2O5;
n.SO3 = liquidComposition.SO3 .* MWinfo.Cat.SO3 ./ MWinfo.MW.SO3;
n.Fe2O3 = liquidComposition.Fe2O3 .* MWinfo.Cat.Fe2O3 ./ MWinfo.MW.Fe2O3;

% Do not use sum(...,'omitnan'): every NaN calculation input must propagate.
cationTotal_liq_scalar = ...
    n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + n.MgO + ...
    n.CaO + n.Na2O + n.K2O + n.V2O3 + n.Cr2O3 + n.NiO + ...
    n.P2O5 + n.SO3 + n.Fe2O3;

% Liquid cation fractions. A zero or NaN total naturally produces NaN/Inf.
XSiO2_scalar = n.SiO2 ./ cationTotal_liq_scalar;
XTiO2_scalar = n.TiO2 ./ cationTotal_liq_scalar;
XAlO1_5_scalar = n.Al2O3 ./ cationTotal_liq_scalar;
XNaO0_5_scalar = n.Na2O ./ cationTotal_liq_scalar;
XKO0_5_scalar = n.K2O ./ cationTotal_liq_scalar;

% Additional diagnostic Mg# using FeO or FeOt as selected above.
MgNumber_liq_scalar = n.MgO ./ (n.MgO + n.FeO);

% Beattie (1993) silica activity as used in Putirka (2008), Equation (41).
aSiO2_liq_scalar = ...
    ((3 .* XSiO2_scalar) .^ (-2)) .* ...
    ((1 - XAlO1_5_scalar) .^ (7 ./ 2)) .* ...
    ((1 - XTiO2_scalar) .^ 7);

% Expand composition-dependent scalars to the temperature-vector length.
SiO2_liq = repmat(liquidComposition.SiO2, nT, 1);
TiO2_liq = repmat(liquidComposition.TiO2, nT, 1);
Al2O3_liq = repmat(liquidComposition.Al2O3, nT, 1);
FeO_liq = repmat(liquidComposition.FeO, nT, 1);
MnO_liq = repmat(liquidComposition.MnO, nT, 1);
MgO_liq = repmat(liquidComposition.MgO, nT, 1);
CaO_liq = repmat(liquidComposition.CaO, nT, 1);
Na2O_liq = repmat(liquidComposition.Na2O, nT, 1);
K2O_liq = repmat(liquidComposition.K2O, nT, 1);
V2O3_liq = repmat(liquidComposition.V2O3, nT, 1);
Cr2O3_liq = repmat(liquidComposition.Cr2O3, nT, 1);
NiO_liq = repmat(liquidComposition.NiO, nT, 1);
P2O5_liq = repmat(liquidComposition.P2O5, nT, 1);
SO3_liq = repmat(liquidComposition.SO3, nT, 1);
Fe2O3_liq = repmat(liquidComposition.Fe2O3, nT, 1);
F_liq = repmat(liquidComposition.F, nT, 1);
Cl_liq = repmat(liquidComposition.Cl, nT, 1);
FeO_source = repmat(string(liquidComposition.FeO_source), nT, 1);

cationTotal_liq = repmat(cationTotal_liq_scalar, nT, 1);
XSiO2_liq = repmat(XSiO2_scalar, nT, 1);
XTiO2_liq = repmat(XTiO2_scalar, nT, 1);
XAlO1_5_liq = repmat(XAlO1_5_scalar, nT, 1);
XNaO0_5_liq = repmat(XNaO0_5_scalar, nT, 1);
XKO0_5_liq = repmat(XKO0_5_scalar, nT, 1);
MgNumber_liq = repmat(MgNumber_liq_scalar, nT, 1);
aSiO2_liq = repmat(aSiO2_liq_scalar, nT, 1);

% Pressure calculation. No finite-value guard is applied: NaN and Inf are
% retained for diagnosis. Equation (42) uses temperature in degreeC.
P_kbar = ...
    231.5 ...
    + 0.186 .* T_degreeC ...
    + 0.1244 .* T_degreeC .* log(aSiO2_liq) ...
    - 528.5 .* sqrt(aSiO2_liq) ...
    + 103.3 .* XTiO2_liq ...
    + 69.9 .* (XNaO0_5_liq + XKO0_5_liq) ...
    + 77.3 .* (XAlO1_5_liq ./ (XAlO1_5_liq + XSiO2_liq));

P_GPa = P_kbar ./ 10;

% Diagnostic applicability flags. These do not verify Ol-Opx-Liquid
% equilibrium and therefore do not by themselves establish applicability.
isWithinCalibrationTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 825 & T_degreeC <= 2000;
isWithinCalibrationPRange = ...
    isfinite(P_kbar) & P_kbar >= 0.001 & P_kbar <= 70;
isWithinCalibrationSiO2Range = ...
    isfinite(SiO2_liq) & SiO2_liq >= 31.5 & SiO2_liq <= 70;
isHighMgNumberSubset = ...
    isfinite(MgNumber_liq) & MgNumber_liq > 0.75;
isSilicaActivityFinitePositive = ...
    isfinite(aSiO2_liq) & aSiO2_liq > 0;
isOlOpxLiquidEquilibriumVerified = false(nT, 1);

% Pack outputs using vectors of equal height.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.SiO2_liq = SiO2_liq;
row.TiO2_liq = TiO2_liq;
row.Al2O3_liq = Al2O3_liq;
row.FeO_liq = FeO_liq;
row.FeO_source = FeO_source;
row.MnO_liq = MnO_liq;
row.MgO_liq = MgO_liq;
row.CaO_liq = CaO_liq;
row.Na2O_liq = Na2O_liq;
row.K2O_liq = K2O_liq;
row.V2O3_liq = V2O3_liq;
row.Cr2O3_liq = Cr2O3_liq;
row.NiO_liq = NiO_liq;
row.P2O5_liq = P2O5_liq;
row.SO3_liq = SO3_liq;
row.Fe2O3_liq = Fe2O3_liq;
row.F_liq = F_liq;
row.Cl_liq = Cl_liq;

row.cationTotal_liq = cationTotal_liq;
row.XSiO2_liq = XSiO2_liq;
row.XTiO2_liq = XTiO2_liq;
row.XAlO1_5_liq = XAlO1_5_liq;
row.XNaO0_5_liq = XNaO0_5_liq;
row.XKO0_5_liq = XKO0_5_liq;
row.MgNumber_liq = MgNumber_liq;
row.aSiO2_liq = aSiO2_liq;

% P_kbar is the primary launcher-compatible output. Equation-specific aliases
% are retained for backward compatibility with the original implementation.
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.P_Eq42_kbar = P_kbar;
row.P_Eq42_GPa = P_GPa;
row.P_uncertainty_1sigma_kbar = repmat(2.87, nT, 1);

row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinCalibrationPRange = isWithinCalibrationPRange;
row.isWithinCalibrationSiO2Range = isWithinCalibrationSiO2Range;
row.isHighMgNumberSubset = isHighMgNumberSubset;
row.isSilicaActivityFinitePositive = isSilicaActivityFinitePositive;
row.isOlOpxLiquidEquilibriumVerified = isOlOpxLiquidEquilibriumVerified;

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Repeat available Liquid identifiers to match the number of temperature rows.

nRows = height(row);
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    row.liq_Index = repmat(string(getTableScalarText(data_liq, 'Index')), nRows, 1);
end
if any(strcmp(variableNames, 'Experiment'))
    row.liq_Experiment = ...
        repmat(string(getTableScalarText(data_liq, 'Experiment')), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    row.liq_Citation = ...
        repmat(string(getTableScalarText(data_liq, 'Citation')), nRows, 1);
end

end

function items = buildLiquidList(liq)
% buildLiquidList
% Construct preallocated display labels for the Liquid selection dialog.

nRows = height(liq);
items = cell(nRows, 1);
variableNames = liq.Properties.VariableNames;

hasIndex = any(strcmp(variableNames, 'Index'));
hasExperiment = any(strcmp(variableNames, 'Experiment'));
hasCitation = any(strcmp(variableNames, 'Citation'));

for i = 1:nRows
    rowData = liq(i, :);
    labelParts = cell(4, 1);
    nLabelParts = 1;
    labelParts{nLabelParts} = ['Row ' num2str(i)];

    if hasIndex
        nLabelParts = nLabelParts + 1;
        labelParts{nLabelParts} = ...
            ['Index=' getTableScalarText(rowData, 'Index')];
    end
    if hasExperiment
        nLabelParts = nLabelParts + 1;
        labelParts{nLabelParts} = getTableScalarText(rowData, 'Experiment');
    end
    if hasCitation
        nLabelParts = nLabelParts + 1;
        labelParts{nLabelParts} = getTableScalarText(rowData, 'Citation');
    end

    items{i} = strjoin(labelParts(1:nLabelParts), ' | ');
end

end

function label = getLiquidLabel(data_liq, rowNumber)
% getLiquidLabel
% Return a compact label for messages and warnings.

label = ['Liquid Row ' num2str(rowNumber)];
variableNames = data_liq.Properties.VariableNames;

if any(strcmp(variableNames, 'Index'))
    label = [label ' (Index=' getTableScalarText(data_liq, 'Index') ')'];
end

end

function textValue = getTableScalarText(tbl, variableName)
% getTableScalarText
% Convert the first value of a one-row table variable to display text.

rawValue = tbl.(variableName);
textValue = '';

if isempty(rawValue)
    return;
end

if iscell(rawValue)
    if isempty(rawValue{1})
        return;
    end
    firstValue = rawValue{1};
else
    firstValue = rawValue(1);
end

if isstring(firstValue)
    if ismissing(firstValue)
        return;
    end
    textValue = char(firstValue);
elseif ischar(firstValue)
    textValue = firstValue;
elseif iscategorical(firstValue)
    textValue = char(string(firstValue));
elseif isnumeric(firstValue) || islogical(firstValue)
    textValue = char(string(firstValue));
else
    textValue = char(string(firstValue));
end

end

function value = getLiquidOxide(data_liq, oxide)
% getLiquidOxide
% Retrieve a scalar Liquid oxide value. Missing, empty, or non-numeric data
% are represented by NaN and are never converted to zero.

columnName = findOxideColumn(data_liq.Properties.VariableNames, oxide);
if isempty(columnName)
    value = NaN;
    return;
end

rawValue = data_liq.(columnName);
value = toScalarDouble(rawValue);

end

function columnName = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match oxide columns while ignoring spaces, underscores, and hyphens.
% Both '<oxide>' and '<oxide>Value' naming styles are accepted.

canonicalNames = cell(size(variableNames));
for i = 1:numel(variableNames)
    canonicalNames{i} = canonicalizeName(variableNames{i});
end

canonicalOxide = canonicalizeName(oxide);
targetNames = {[canonicalOxide 'value'], canonicalOxide};

columnName = '';
for i = 1:numel(targetNames)
    matchedIndex = find(strcmp(canonicalNames, targetNames{i}), 1, 'first');
    if ~isempty(matchedIndex)
        columnName = variableNames{matchedIndex};
        return;
    end
end

end

function name = canonicalizeName(name)
% canonicalizeName
% Convert a variable or oxide name to a comparison-safe lowercase string.

name = lower(char(string(name)));
name = strrep(name, ' ', '');
name = strrep(name, '_', '');
name = strrep(name, '-', '');

end

function value = toScalarDouble(rawValue)
% toScalarDouble
% Convert the first value of a one-row table variable to double. Missing,
% empty, or non-convertible values remain NaN.

value = NaN;

if isempty(rawValue)
    return;
end

if isnumeric(rawValue) || islogical(rawValue)
    value = double(rawValue(1));
    return;
end

if isstring(rawValue)
    if ismissing(rawValue(1))
        return;
    end
    value = str2double(rawValue(1));
    return;
end

if ischar(rawValue)
    value = str2double(string(rawValue));
    return;
end

if iscategorical(rawValue)
    value = str2double(string(rawValue(1)));
    return;
end

if iscell(rawValue)
    if isempty(rawValue{1})
        return;
    end

    firstValue = rawValue{1};
    if isnumeric(firstValue) || islogical(firstValue)
        value = double(firstValue(1));
    elseif isstring(firstValue) || ischar(firstValue) || ...
            iscategorical(firstValue)
        value = str2double(string(firstValue));
    end
end

end
