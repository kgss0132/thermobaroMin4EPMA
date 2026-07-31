function results = WittEickschenSeck1991Opx(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/WittEickschenSeck1991Opx.m
% Tested with MATLAB R2024b
%
% Empirical Cr-Al-in-Orthopyroxene thermometer for spinel peridotite
% Witt-Eickschen, G. and Seck, H.A. (1991)
% Contributions to Mineralogy and Petrology, 106, 431-439
% DOI: https://doi.org/10.1007/BF00321986
%
% -------------------------------------------------------------------------
% OVERVIEW
%
% This function interactively selects one Orthopyroxene analysis and
% calculates temperature using Equation (7) of Witt-Eickschen and Seck
% (1991):
%
%   T(degreeC) = 636.54 ...
%              + 2088.21  * XM1Al_Opx ...
%              + 14527.32 * XM1Cr_Opx
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Opx analysis, one output row is returned for every pressure
% supplied in P_kbar. Equation (7) contains no pressure term, so all pressure
% rows for one selected Opx analysis contain the same calculated
% temperature. This interface is compatible with both
% startThermoCalc_fixedP and startThermoCalc_rangeP.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION CONTEXT AND APPLICATION NOTES
%
% Witt-Eickschen and Seck (1991) examined more than 100 natural spinel
% peridotites from the Rhenish Volcanic Province. Their empirical
% temperatures were referenced to the Brey and Kohler (1990) Ca-in-Opx
% thermometer. This calibration context and Equations (6) and (7) are
% summarized in the abstract on p. 431.
%
% IMPORTANT APPLICATION LIMITATIONS
%
%   1) Equation (7) is intended for Orthopyroxene from equilibrated natural
%      spinel peridotite. It must not be treated as a universal
%      single-Orthopyroxene thermometer for garnet peridotite, plagioclase
%      peridotite, pyroxenite, mafic rocks, or magmatic Orthopyroxene
%      outside the natural spinel-peridotite calibration domain
%      (abstract, p. 431; article, pp. 431-439).
%
%   2) The authors explicitly describe Equation (7) as convenient for
%      practical use but applicable only over a limited range of XM1Al_Opx
%      and XM1Cr_Opx. A single numerical rectangular range for these two
%      variables is not formally stated in the abstract or as a universal
%      limit for Equation (7) (abstract, p. 431).
%
%   3) The original paper does not define a single numerical pressure
%      calibration interval for Equation (7), and the equation contains no
%      pressure term. P_kbar is therefore retained only for interface
%      compatibility and traceability. A pressure-range validity test cannot
%      be made from a published Equation (7)-specific numerical interval.
%
%   4) The original paper also does not state a single formal numerical
%      temperature calibration interval for Equation (7). Consequently, a
%      published temperature-range test cannot be applied independently of
%      the Opx composition.
%
%   5) Witt-Eickschen and Seck (1991) note that recent thermal disturbance
%      or rapid cooling can produce disagreement between Ca-in-Opx and
%      Al-in-Opx temperatures because steady-state equilibrium may not have
%      been attained (abstract, p. 431). Zoned grains, exsolution,
%      porphyroclast-neoblast differences, reaction textures, and altered
%      domains therefore require careful screening.
%
%   6) Equation (7) uses only Opx composition, but its natural calibration
%      assumes the spinel-peridotite mineral assemblage and chemical
%      controls represented by the calibration samples. Textural and
%      petrological evidence for equilibrium remains essential.
%
% PRACTICAL SOFTWARE SCREENING USED HERE
%
% The supplied original source code used the following provisional
% composition-screening limits:
%
%   XM1Al_Opx : 0.00-0.20
%   XM1Cr_Opx : 0.00-0.05
%
% These values are retained only as implementation-specific, non-stopping
% screening limits. They must not be described as formal numerical
% calibration limits reported by Witt-Eickschen and Seck (1991).
%
% For the requested temperature screening, the algebraic temperature
% envelope implied by that provisional rectangular composition window is:
%
%   Temperature : 636.54-1780.548 degreeC
%
% This very broad interval is not a published calibration range. It is used
% only to flag results that necessarily lie outside the retained provisional
% composition window. The individual XM1Al_Opx and XM1Cr_Opx messages are
% more informative than this derived temperature envelope.
%
% This implementation issues non-stopping fprintf messages when:
%   1) the absence of a published numerical pressure range must be stated,
%   2) the absence of a formal numerical temperature interval must be
%      stated,
%   3) a finite result lies outside the provisional algebraic temperature
%      envelope above,
%   4) finite XM1Al_Opx or XM1Cr_Opx lies outside the provisional
%      implementation-specific screening window,
%   5) an explicitly stored calculation input is NaN,
%   6) Opx Na >= Cr + Ti + Fe3 triggers the Carswell-style allocation
%      caution,
%   7) a derived site quantity is invalid, or
%   8) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
%
% rawdata_struct must contain:
%   rawdata_struct.Opx : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns should contain
% normalized Opx cations, preferably on a 6-oxygen basis.
%
% Required Opx variables:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   Ti_cation_apfu
%   Cr_cation_apfu
%
% Optional Opx variables:
%   Mn_cation_apfu
%   Fe3_cation_apfu
%
% Existing NaN values are retained and propagated; they are never replaced
% by zero. Zero is assigned only when an optional column itself is absent.
% All finite stored mineral-composition inputs must be greater than or equal
% to zero. Negative finite values and Inf stop the calculation.
%
% Equation (7) uses Al, Cr, Ti, and Na directly. Fe3 is used only in the
% Carswell-style allocation diagnostic. Other Opx variables are retained in
% the output for traceability and interface consistency.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and liquid NaN warnings is not applicable.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% The M1-Al term follows the Carswell (1989)-type expression cited by
% Witt-Eickschen and Seck (1991):
%
%   XM1Al_Opx = 0.5 * (Al - Cr - 2*Ti + Na)
%
% Chromium is treated as occupying the Opx M1 site:
%
%   XM1Cr_Opx = Cr
%
% Witt-Eickschen and Seck (1991), Equation (7):
%
%   T(degreeC) = 636.54 ...
%              + 2088.21  * XM1Al_Opx ...
%              + 14527.32 * XM1Cr_Opx
%
% Temperature is returned in both degreeC and Kelvin:
%
%   T_K = T_degreeC + 273.15
%
% -------------------------------------------------------------------------
% Syntax:
%   results = WittEickschenSeck1991Opx(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing an Opx table
%   P_kbar         : finite non-negative numeric scalar or vector; stored in
%                    the output but not used by Equation (7)
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Opx analysis. T_K, T_degreeC, and T_deg are supplied for
%             downstream launcher and plotting compatibility.
%

%% Input validation
if nargin < 2
    error(['WittEickschenSeck1991Opx requires ' ...
           '(rawdata_struct, P_kbar).']);
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
disp('=== Step 1: Preparing Orthopyroxene dataset ===');

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end
if isempty(rawdata_struct.Opx)
    error('rawdata_struct.Opx is empty.');
end

dataset_opx = rawdata_struct.Opx;
validateRequiredColumns(dataset_opx);

disp('=== Preparing Orthopyroxene dataset has been finished ===');

%% 2) Initialize output container and screening settings
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Provisional implementation-specific composition-screening limits retained
% from the supplied original source. These are not formal published limits.
screeningXAlM1_min = 0.00;
screeningXAlM1_max = 0.20;
screeningXCrM1_min = 0.00;
screeningXCrM1_max = 0.05;

% Algebraic temperature envelope implied by the provisional rectangular
% composition-screening window above. This is not a published calibration
% range.
screeningT_min_degreeC = ...
    636.54 ...
    + 2088.21 .* screeningXAlM1_min ...
    + 14527.32 .* screeningXCrM1_min;

screeningT_max_degreeC = ...
    636.54 ...
    + 2088.21 .* screeningXAlM1_max ...
    + 14527.32 .* screeningXCrM1_max;

pressureRangeMessageIssued = false;
temperatureRangeMessageIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive selection loop and calculation
dataCodes_opx = dataset_opx{:, 1};
displayCodes_opx = cellstr(string(dataCodes_opx));

disp('=== Step 3: Selecting a data code from the list (Opx) ===');

while true
    % ----- Orthopyroxene selection -----
    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', ...
            'Please select the Opx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_opx, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = string(dataCodes_opx(selectedIdx_opx));
    selectedData_opx = dataset_opx(selectedIdx_opx, :);

    disp(['Opx selected: ' char(selectedCode_opx)]);

    % ----- Calculation -----
    disp('=== Step 4: Calculating the temperature ===');

    nanInputNames = findNaNInputs(selectedData_opx);
    validateNonNegativeInputs(selectedData_opx);

    row = calcTemp(selectedData_opx, P_kbar);

    nRows = height(row);
    row.dataCode_opx = repmat(selectedCode_opx, nRows, 1);
    row = movevars(row, {'dataCode_opx'}, 'Before', 1);

    % Store one completed table block. The complete results table is not
    % enlarged on every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = ...
            [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % ----- Immediate result display -----
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    printTemperatureSummary(selectedCode_opx, row.T_degreeC);

    % Equation (7) has no published numerical pressure interval.
    if ~pressureRangeMessageIssued
        fprintf(2, ...
            ['CAUTION: Witt-Eickschen and Seck (1991) do not state a ' ...
             'numerical pressure calibration range specifically for ' ...
             'Equation (7), and the equation contains no pressure term ' ...
             '(p. 431; article pp. 431-439). P_kbar is stored for ' ...
             'interface compatibility but is not used. Input range = ' ...
             '%.4g-%.4g kbar. A published pressure-range validity test ' ...
             'cannot be performed.\n'], ...
            min(P_kbar), ...
            max(P_kbar));
        pressureRangeMessageIssued = true;
    end

    % The paper does not define a single formal numerical T interval.
    if ~temperatureRangeMessageIssued
        fprintf(2, ...
            ['CAUTION: Witt-Eickschen and Seck (1991) do not state a ' ...
             'single formal numerical temperature calibration interval ' ...
             'for Equation (7). Temperature screening in this function ' ...
             'uses only the algebraic envelope %.4g-%.4g degreeC implied ' ...
             'by the provisional implementation-specific XM1Al_Opx and ' ...
             'XM1Cr_Opx windows; it is not a published calibration ' ...
             'range.\n'], ...
            screeningT_min_degreeC, ...
            screeningT_max_degreeC);
        temperatureRangeMessageIssued = true;
    end

    printTemperatureScreeningWarning( ...
        row.T_degreeC, ...
        screeningT_min_degreeC, ...
        screeningT_max_degreeC, ...
        selectedCode_opx);

    printCompositionScreeningWarnings( ...
        row, ...
        screeningXAlM1_min, ...
        screeningXAlM1_max, ...
        screeningXCrM1_min, ...
        screeningXCrM1_max, ...
        selectedCode_opx);

    % Existing NaN values are retained and never replaced by zero.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Equation (7) input or ' ...
             'allocation-diagnostic input(s) for %s: %s.\n' ...
             '         NaN values were retained and propagated; the ' ...
             'calculated temperature may be NaN.\n'], ...
            char(selectedCode_opx), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Carswell-style M1-Al diagnostic retained from the supplied source.
    if all(isfinite([ ...
            row.Na_opx(1), ...
            row.Cr_opx(1), ...
            row.Ti_opx(1), ...
            row.Fe3_opx(1)])) && ...
            row.Na_opx(1) >= ...
            (row.Cr_opx(1) + row.Ti_opx(1) + row.Fe3_opx(1))

        fprintf(2, ...
            ['WARNING: Opx Na >= Cr + Ti + Fe3 for %s. The implemented ' ...
             'Carswell-style XM1Al_Opx estimate may be inappropriate. ' ...
             'The calculated result has been retained.\n'], ...
            char(selectedCode_opx));
    end

    invalidTerms = findInvalidEquationTerms(row);
    if ~isempty(invalidTerms)
        fprintf(2, ...
            ['WARNING: Invalid Witt-Eickschen-Seck Equation (7) input or ' ...
             'derived term(s) were found for %s: %s.\n' ...
             '         XM1Al_Opx and XM1Cr_Opx must be finite and within ' ...
             'their physical site-fraction bounds of 0-1. Affected ' ...
             'temperature values were retained as NaN.\n'], ...
            char(selectedCode_opx), ...
            char(strjoin(invalidTerms, ', ')));
    end

    printNonFiniteTemperatureWarning( ...
        row.T_degreeC, selectedCode_opx);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'WittEickschenSeck1991Opx', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once after the selection loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'primaryTemperatureEquation', ...
        'Witt-Eickschen and Seck (1991) Equation (7)', ...
    'calibrationBasis', ...
        ['Natural spinel peridotites; reference temperatures from the ' ...
         'Brey and Kohler (1990) Ca-in-Opx thermometer'], ...
    'pressureUsedInEquation', false, ...
    'publishedPressureCalibrationRange', ...
        'Not numerically specified for Equation (7)', ...
    'publishedTemperatureCalibrationRange', ...
        'No single formal numerical interval stated for Equation (7)', ...
    'provisionalXM1AlOpxScreening', ...
        [screeningXAlM1_min, screeningXAlM1_max], ...
    'provisionalXM1CrOpxScreening', ...
        [screeningXCrM1_min, screeningXCrM1_max], ...
    'derivedProvisionalTemperatureEnvelope_degreeC', ...
        [screeningT_min_degreeC, screeningT_max_degreeC], ...
    'screeningRangeStatus', ...
        'Implementation-specific; not a formal published calibration range');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_opx)
% findNaNInputs
% Return names of explicitly stored values used directly by Equation (7) or
% by the Carswell-style allocation diagnostic that contain NaN.
%
% Equation (7) uses Al, Cr, Ti, and Na. Fe3 is included because it enters the
% diagnostic condition Na >= Cr + Ti + Fe3. Missing Fe3 is assigned zero,
% but an explicitly stored Fe3 NaN is retained and reported.

activeVariables = { ...
    'Al_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Na_cation_apfu', ...
    'Fe3_cation_apfu'};

nameBuffer = strings(numel(activeVariables), 1);
nNames = 0;

for i = 1:numel(activeVariables)
    variableName = activeVariables{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        value = data_opx.(variableName);
        if isnumeric(value) && any(isnan(value(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Opx." + string(variableName);
        end
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_opx)
% validateNonNegativeInputs
% Stop when a stored Opx composition value used or retained by this
% implementation is negative or infinite. Zero is allowed. NaN is
% deliberately allowed and propagated.

variablesToCheck = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Fe3_cation_apfu'};

nameBuffer = strings(numel(variablesToCheck), 1);
nInvalid = 0;

for i = 1:numel(variablesToCheck)
    variableName = variablesToCheck{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        value = data_opx.(variableName);
        validateScalarVariable(value, 'Opx', variableName);

        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            nameBuffer(nInvalid) = "Opx." + string(variableName);
        end
    end
end

if nInvalid > 0
    invalidNames = nameBuffer(1:nInvalid);
    error(['WittEickschenSeck1991Opx: Opx composition inputs must ' ...
           'not be negative or infinite. Invalid value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, P_kbar)
% calcTemp
% Calculate Witt-Eickschen and Seck (1991) Equation (7) for one selected
% Opx composition and repeat the pressure-independent result for every
% supplied pressure. Existing NaN values and invalid derived terms are
% retained as NaN.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

opx = prepareOpxRow(data_opx, 'Opx');
site_opx = calcOpxSiteFractions(opx);

% Equation (7) is evaluated only when both site terms are finite and within
% physical site-fraction bounds. Invalid derived terms remain NaN.
if isfinite(site_opx.XM1Al_Opx) && ...
        isfinite(site_opx.XM1Cr_Opx)

    T_scalar_raw_degreeC = ...
        636.54 ...
        + 2088.21 .* site_opx.XM1Al_Opx ...
        + 14527.32 .* site_opx.XM1Cr_Opx;
else
    T_scalar_raw_degreeC = NaN;
end

T_scalar_degreeC = T_scalar_raw_degreeC;
T_scalar_K = T_scalar_degreeC + 273.15;

% Non-positive Kelvin is physically invalid and is retained as NaN.
if ~isfinite(T_scalar_K) || T_scalar_K <= 0
    T_scalar_degreeC = NaN;
    T_scalar_K = NaN;
end

T_raw_degreeC = repmat(T_scalar_raw_degreeC, nP, 1);
T_degreeC = repmat(T_scalar_degreeC, nP, 1);
T_K = repmat(T_scalar_K, nP, 1);

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PressureUsedInEquation = false(nP, 1);
row.PrimaryEquation = ...
    repmat("WittEickschenSeck1991_Eq7", nP, 1);

% Opx inputs. Fe_cation_apfu is retained using the historical Fe2_opx output
% name for backward compatibility, although Fe does not enter Equation (7).
row.Si_opx = repmat(opx.Si, nP, 1);
row.Al_opx = repmat(opx.Al, nP, 1);
row.Fe2_opx = repmat(opx.Fe, nP, 1);
row.Fe3_opx = repmat(opx.Fe3, nP, 1);
row.Mg_opx = repmat(opx.Mg, nP, 1);
row.Ca_opx = repmat(opx.Ca, nP, 1);
row.Na_opx = repmat(opx.Na, nP, 1);
row.Mn_opx = repmat(opx.Mn, nP, 1);
row.Ti_opx = repmat(opx.Ti, nP, 1);
row.Cr_opx = repmat(opx.Cr, nP, 1);

% Site terms and raw values are both retained for diagnostics.
row.XM1Opx_Al_raw = ...
    repmat(site_opx.XM1Al_Opx_raw, nP, 1);
row.XM1Opx_Cr_raw = ...
    repmat(site_opx.XM1Cr_Opx_raw, nP, 1);

row.XM1Opx_Al = repmat(site_opx.XM1Al_Opx, nP, 1);
row.XM1Opx_Cr = repmat(site_opx.XM1Cr_Opx, nP, 1);

% Standardized aliases.
row.XM1Al_Opx = row.XM1Opx_Al;
row.XM1Cr_Opx = row.XM1Opx_Cr;

row.T_Eq7_raw_degreeC = T_raw_degreeC;
row.T_K = T_K;
row.T_degreeC = T_degreeC;
row.T_deg = T_degreeC;

end

function opx = prepareOpxRow(data_opx, mineralLabel)
% prepareOpxRow
% Extract one Opx composition. Existing NaN values remain NaN. Missing
% optional columns are assigned zero.

if height(data_opx) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

opx = struct();

opx.Si = getVarRequired( ...
    data_opx, 'Si_cation_apfu', mineralLabel);
opx.Al = getVarRequired( ...
    data_opx, 'Al_cation_apfu', mineralLabel);
opx.Fe = getVarRequired( ...
    data_opx, 'Fe_cation_apfu', mineralLabel);
opx.Mg = getVarRequired( ...
    data_opx, 'Mg_cation_apfu', mineralLabel);
opx.Ca = getVarRequired( ...
    data_opx, 'Ca_cation_apfu', mineralLabel);
opx.Na = getVarRequired( ...
    data_opx, 'Na_cation_apfu', mineralLabel);
opx.Ti = getVarRequired( ...
    data_opx, 'Ti_cation_apfu', mineralLabel);
opx.Cr = getVarRequired( ...
    data_opx, 'Cr_cation_apfu', mineralLabel);

opx.Mn = getVarOptional( ...
    data_opx, 'Mn_cation_apfu', 0, mineralLabel);
opx.Fe3 = getVarOptional( ...
    data_opx, 'Fe3_cation_apfu', 0, mineralLabel);

end

function site = calcOpxSiteFractions(opx)
% calcOpxSiteFractions
% Calculate the Opx M1-Al and M1-Cr quantities used by Equation (7).
% Invalid or non-finite derived quantities are represented by NaN rather
% than stopping the complete calculation.

site = struct();

% Carswell (1989)-type M1-Al expression cited by Witt-Eickschen and Seck
% (1991).
if all(isfinite([opx.Al, opx.Cr, opx.Ti, opx.Na]))
    XM1Al_Opx_raw = ...
        0.5 .* (opx.Al - opx.Cr - 2 .* opx.Ti + opx.Na);
else
    XM1Al_Opx_raw = NaN;
end

if isfinite(XM1Al_Opx_raw) && ...
        XM1Al_Opx_raw >= 0 && XM1Al_Opx_raw <= 1
    XM1Al_Opx = XM1Al_Opx_raw;
else
    XM1Al_Opx = NaN;
end

% Chromium is treated as an M1-site quantity in the supplied Equation (7)
% implementation.
XM1Cr_Opx_raw = opx.Cr;

if isfinite(XM1Cr_Opx_raw) && ...
        XM1Cr_Opx_raw >= 0 && XM1Cr_Opx_raw <= 1
    XM1Cr_Opx = XM1Cr_Opx_raw;
else
    XM1Cr_Opx = NaN;
end

site.XM1Al_Opx_raw = XM1Al_Opx_raw;
site.XM1Cr_Opx_raw = XM1Cr_Opx_raw;
site.XM1Al_Opx = XM1Al_Opx;
site.XM1Cr_Opx = XM1Cr_Opx;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify non-finite or physically invalid site terms and temperatures.

termBuffer = strings(8, 1);
nTerms = 0;

if ~isfinite(row.XM1Opx_Al_raw(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XM1Al_Opx raw value";
elseif row.XM1Opx_Al_raw(1) < 0 || ...
        row.XM1Opx_Al_raw(1) > 1
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XM1Al_Opx outside 0-1";
end

if ~isfinite(row.XM1Opx_Cr_raw(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XM1Cr_Opx raw value";
elseif row.XM1Opx_Cr_raw(1) < 0 || ...
        row.XM1Opx_Cr_raw(1) > 1
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XM1Cr_Opx outside 0-1";
end

if ~isfinite(row.XM1Opx_Al(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XM1Al_Opx used by Eq. (7)";
end

if ~isfinite(row.XM1Opx_Cr(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XM1Cr_Opx used by Eq. (7)";
end

if any(~isfinite(row.T_K) | row.T_K <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "T_K";
end

invalidTerms = unique(termBuffer(1:nTerms), 'stable');

end

function printTemperatureSummary(selectedCode_opx, temperatureValues)
% printTemperatureSummary
% Display one temperature or first-to-last values for a pressure vector.

if isscalar(temperatureValues)
    disp([char(selectedCode_opx) ': ' ...
        num2str(temperatureValues) ' degreeC']);
else
    disp([char(selectedCode_opx) ': ' ...
        num2str(temperatureValues(1)) ' to ' ...
        num2str(temperatureValues(end)) ' degreeC']);
end

end

function printTemperatureScreeningWarning( ...
        temperatureValues, minimumTemperature, maximumTemperature, ...
        selectedCode_opx)
% printTemperatureScreeningWarning
% Warn when finite temperatures lie outside the algebraic envelope implied
% by the provisional implementation-specific composition windows. Results
% are retained.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);

    fprintf(2, ...
        ['WARNING: Calculated Equation (7) temperature is outside the ' ...
         'provisional algebraic screening envelope of %.4g-%.4g degreeC. ' ...
         '%d of %d finite point(s) are outside; calculated finite range ' ...
         '= %.4g-%.4g degreeC for %s. This envelope is derived from the ' ...
         'implementation-specific XM1Al_Opx and XM1Cr_Opx screening ' ...
         'windows and is not a formal range published by ' ...
         'Witt-Eickschen and Seck (1991).\n'], ...
        minimumTemperature, ...
        maximumTemperature, ...
        sum(outsideMask), ...
        sum(finiteMask), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(selectedCode_opx));
end

end

function printCompositionScreeningWarnings( ...
        row, XAlMinimum, XAlMaximum, XCrMinimum, XCrMaximum, ...
        selectedCode_opx)
% printCompositionScreeningWarnings
% Apply provisional implementation-specific composition screening. These
% limits are retained from the supplied source and are not formal numerical
% calibration limits from the original paper.

XAlValue = row.XM1Opx_Al_raw(1);
XCrValue = row.XM1Opx_Cr_raw(1);

if isfinite(XAlValue) && ...
        (XAlValue < XAlMinimum || XAlValue > XAlMaximum)

    fprintf(2, ...
        ['WARNING: XM1Al_Opx = %.6g for %s is outside the provisional ' ...
         'implementation-specific screening interval %.4g-%.4g. ' ...
         'Witt-Eickschen and Seck (1991) state that Equation (7) is valid ' ...
         'only over a limited XM1Al_Opx-XM1Cr_Opx domain, but this ' ...
         'numerical interval is not presented here as a formal published ' ...
         'calibration range. The result has been retained.\n'], ...
        XAlValue, ...
        char(selectedCode_opx), ...
        XAlMinimum, ...
        XAlMaximum);
end

if isfinite(XCrValue) && ...
        (XCrValue < XCrMinimum || XCrValue > XCrMaximum)

    fprintf(2, ...
        ['WARNING: XM1Cr_Opx = %.6g for %s is outside the provisional ' ...
         'implementation-specific screening interval %.4g-%.4g. ' ...
         'Witt-Eickschen and Seck (1991) state that Equation (7) is valid ' ...
         'only over a limited XM1Al_Opx-XM1Cr_Opx domain, but this ' ...
         'numerical interval is not presented here as a formal published ' ...
         'calibration range. The result has been retained.\n'], ...
        XCrValue, ...
        char(selectedCode_opx), ...
        XCrMinimum, ...
        XCrMaximum);
end

end

function printNonFiniteTemperatureWarning( ...
        temperatureValues, selectedCode_opx)
% printNonFiniteTemperatureWarning
% Report NaN and Inf temperatures without modifying or deleting them.

invalidMask = ~isfinite(temperatureValues);

if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Witt-Eickschen-Seck Equation (7) ' ...
         'temperature values were calculated for %s (%d of %d points; ' ...
         'NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the ' ...
         'calculation has not been stopped.\n'], ...
        char(selectedCode_opx), ...
        sum(invalidMask), ...
        numel(temperatureValues), ...
        sum(isnan(temperatureValues)), ...
        sum(isinf(temperatureValues)));
end

end

function validateRequiredColumns(tbl)
% validateRequiredColumns
% Verify all normalized Opx cation columns required by the interface and
% retained output.

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};

    if ~ismember(variableName, tbl.Properties.VariableNames)
        error('Opx table must contain variable: %s', variableName);
    end
end

end

function value = getVarRequired(tbl, variableName, mineralLabel)
% getVarRequired
% Read a required numeric scalar while retaining NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', ...
        mineralLabel, variableName);
end

value = tbl.(variableName);
validateScalarVariable(value, mineralLabel, variableName);

end

function value = getVarOptional( ...
        tbl, variableName, defaultValue, mineralLabel)
% getVarOptional
% Read an optional numeric scalar. An absent column receives defaultValue;
% an explicitly stored NaN remains NaN.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);
    validateScalarVariable(value, mineralLabel, variableName);
else
    value = defaultValue;
end

end

function validateScalarVariable(value, mineralLabel, variableName)
% validateScalarVariable
% Require one numeric scalar. NaN is allowed; negative finite values and Inf
% are handled by validateNonNegativeInputs.

if ~isnumeric(value) || ~isscalar(value)
    error('%s variable %s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end
