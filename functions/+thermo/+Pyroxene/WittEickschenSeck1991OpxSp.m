function results = WittEickschenSeck1991OpxSp(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/WittEickschenSeck1991OpxSp.m
% Tested with MATLAB R2024b
%
% Empirical Orthopyroxene-Spinel-Olivine thermometer for spinel peridotite
% Witt-Eickschen, G. and Seck, H.A. (1991)
% Contributions to Mineralogy and Petrology, 106, 431-439
% DOI: https://doi.org/10.1007/BF00321986
%
% -------------------------------------------------------------------------
% OVERVIEW
%
% This function interactively selects one Orthopyroxene analysis, one
% Spinel analysis, and one Olivine analysis and calculates temperature using
% Equation (6) of Witt-Eickschen and Seck (1991).
%
% The equilibrium coefficient used in the original paper is:
%
%          (XFo_Ol)^2 * XM1Al_Opx
%   KD = --------------------------------
%        (XAl_Sp)^2 * (XMg_Sp)^2 * XM1Mg_Opx
%
% IMPORTANT:
% The Spinel Mg fraction is squared in the published KD definition. The
% supplied original script used only one power of XMg_Sp. This modified
% implementation uses (XMg_Sp)^2, following Witt-Eickschen and Seck (1991,
% abstract and Equation (6), p. 431).
%
% Equation (6):
%
%   T(degreeC) = 2248.25 ...
%              + 991.58  * lnKD ...
%              + 153.32  * (lnKD)^2 ...
%              + 539.05  * YCr_Sp ...
%              - 2005.74 * (YCr_Sp)^2
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Opx-Sp-Ol triplet, one output row is returned for every
% pressure supplied in P_kbar. Equation (6) contains no pressure term, so
% all pressure rows for one triplet contain the same calculated temperature.
% This interface is compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP.
%
% The function is designed for repeated calculations. Each completed result
% is stored as one table block in a preallocated cell buffer, and all blocks
% are concatenated only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION CONTEXT AND APPLICATION NOTES
%
% Witt-Eickschen and Seck (1991) evaluated Ca and Al solubilities in
% Orthopyroxene from more than 100 natural spinel peridotites from the
% Rhenish Volcanic Province. Equation (6) was fitted to temperatures
% calculated using the Brey and Kohler (1990) Ca-in-Opx thermometer. The
% calibration basis, KD definition, and Equations (6) and (7) are summarized
% in the abstract on p. 431; the complete article spans pp. 431-439.
%
% IMPORTANT APPLICATION LIMITATIONS
%
%   1) Equation (6) is intended for texturally and chemically equilibrated
%      Orthopyroxene-Spinel-Olivine assemblages in natural spinel
%      peridotite. It must not be treated as a universal thermometer for
%      garnet peridotite, plagioclase peridotite, pyroxenite, mafic rocks,
%      or unrelated magmatic minerals outside the calibration assemblage
%      (abstract, p. 431; article, pp. 431-439).
%
%   2) The selected Opx, Spinel, and Olivine must represent the same
%      equilibrium generation. Do not arbitrarily combine cores, rims,
%      porphyroclasts, neoblasts, exsolution-related domains, or minerals
%      affected by different reaction or metasomatic stages.
%
%   3) Witt-Eickschen and Seck (1991) note that recent thermal perturbation
%      or rapid cooling may prevent steady-state equilibrium and produce
%      discrepancies between Ca-in-Opx and Al-in-Opx temperatures
%      (abstract, p. 431). Mineral zoning, exsolution, and reaction textures
%      therefore require careful screening.
%
%   4) Equation (6) contains no explicit pressure term. The paper does not
%      state a single numerical pressure calibration interval specifically
%      for Equation (6). P_kbar is retained only for interface compatibility
%      and traceability. Absence of a pressure term must not be interpreted
%      as universal pressure independence outside the spinel-peridotite
%      calibration domain.
%
%   5) The paper does not state a single formal numerical temperature
%      calibration interval for Equation (6). A published rectangular
%      temperature-range test therefore cannot be performed independently
%      of the mineral compositions and KD.
%
%   6) Spinel Fe3 affects XAl_Sp, XMg_Sp, and YCr_Sp through the calculated
%      Fe2 and trivalent-site denominator. An explicitly stored NaN Fe3
%      value is retained rather than replaced by zero. If the Fe3 column is
%      absent, Fe3 is assumed to be zero.
%
%   7) Fe_cation_apfu is treated as total Fe for Opx, Spinel, and Olivine.
%      When Fe3_cation_apfu is available:
%
%        Fe2_cation_apfu = Fe_cation_apfu - Fe3_cation_apfu
%
%      If Fe3_cation_apfu is absent, Fe3 is assumed to be zero. A derived
%      negative Fe2 value is prohibited.
%
% PRACTICAL SOFTWARE SCREENING USED HERE
%
% The supplied original source code used the following approximate lnKD
% screening interval:
%
%   lnKD : -3.5 to -0.3
%
% This interval is retained only as an implementation-specific,
% non-stopping screening limit. It is not presented here as a formal
% numerical calibration range reported by Witt-Eickschen and Seck (1991).
%
% Because no formal numerical Equation (6)-specific temperature or pressure
% range is stated in the accessible original-paper description, this
% implementation does not invent numerical T-P calibration limits. Instead,
% it prints a one-time CAUTION describing that limitation. It additionally
% issues non-stopping fprintf messages when:
%
%   1) finite lnKD lies outside the provisional -3.5 to -0.3 interval,
%   2) an explicitly stored calculation input is NaN,
%   3) Opx Na >= Cr + Ti + Fe3 triggers the Carswell-style allocation
%      caution,
%   4) a required denominator, site term, KD, logarithm, or temperature is
%      invalid, or
%   5) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
%
% rawdata_struct must contain:
%   rawdata_struct.Opx     : table
%   rawdata_struct.Spinel  : table
%   rawdata_struct.Olivine : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialogs. Remaining columns should contain
% consistently normalized cations.
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
% Required Spinel variables:
%   Al_cation_apfu
%   Cr_cation_apfu
%   Mg_cation_apfu
%   Fe_cation_apfu
%
% Optional Spinel variable:
%   Fe3_cation_apfu
%
% Required Olivine variables:
%   Fe_cation_apfu
%   Mg_cation_apfu
%
% Optional Olivine variables:
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Ni_cation_apfu
%   Fe3_cation_apfu
%
% Existing NaN values are retained and propagated; they are never replaced
% by zero. Zero is assigned only when an optional column itself is absent.
% All finite stored mineral-composition inputs must be greater than or equal
% to zero. Negative finite values and Inf stop the calculation. A derived
% Fe2 value below zero also stops the calculation.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and liquid NaN warnings is not applicable.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% 1) Orthopyroxene M1 terms
%
%   XM1Al_Opx = 0.5 * (Al - Cr - 2*Ti + Na)
%
%   M1_fixed = XM1Al_Opx + Cr + Ti + Fe3
%   M1_remaining = 1 - M1_fixed
%
%   XM1Mg_Opx is estimated by filling the remaining M1 capacity with Mg and
%   Fe2 in proportion to Mg/(Mg + Fe2) and Fe2/(Mg + Fe2).
%
% 2) Spinel terms
%
%   XAl_Sp = Al / (Al + Cr + Fe3)
%   XMg_Sp = Mg / (Mg + Fe2)
%   YCr_Sp = Cr / (Al + Cr + Fe3)
%
% 3) Olivine forsterite component
%
%   XFo_Ol = Mg / (Mg + Fe2 + Mn + Ca + Ni)
%
% 4) Equilibrium coefficient
%
%          (XFo_Ol)^2 * XM1Al_Opx
%   KD = --------------------------------
%        (XAl_Sp)^2 * (XMg_Sp)^2 * XM1Mg_Opx
%
% 5) Temperature
%
%   T(degreeC) = 2248.25 ...
%              + 991.58  * lnKD ...
%              + 153.32  * (lnKD)^2 ...
%              + 539.05  * YCr_Sp ...
%              - 2005.74 * (YCr_Sp)^2
%
%   T_K = T_degreeC + 273.15
%
% -------------------------------------------------------------------------
% Syntax:
%   results = WittEickschenSeck1991OpxSp(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Opx, Spinel, and Olivine tables
%   P_kbar         : finite non-negative numeric scalar or vector; stored in
%                    the output but not used by Equation (6)
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Opx-Sp-Ol triplet. T_K, T_degreeC, and T_deg are supplied for
%             downstream launcher and plotting compatibility.
%

%% Input validation
if nargin < 2
    error(['WittEickschenSeck1991OpxSp requires ' ...
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

%% 1) Retrieve cation datasets
disp('=== Step 1: Preparing mineral datasets ===');

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end
if ~isfield(rawdata_struct, 'Spinel') || ~istable(rawdata_struct.Spinel)
    error('rawdata_struct must contain table: rawdata_struct.Spinel');
end
if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end
if isempty(rawdata_struct.Opx)
    error('rawdata_struct.Opx is empty.');
end
if isempty(rawdata_struct.Spinel)
    error('rawdata_struct.Spinel is empty.');
end
if isempty(rawdata_struct.Olivine)
    error('rawdata_struct.Olivine is empty.');
end

dataset_opx = rawdata_struct.Opx;
dataset_sp = rawdata_struct.Spinel;
dataset_ol = rawdata_struct.Olivine;

validateRequiredColumns(dataset_opx, 'Opx');
validateRequiredColumns(dataset_sp, 'Spinel');
validateRequiredColumns(dataset_ol, 'Olivine');

disp('=== Preparing mineral datasets has been finished ===');

%% 2) Initialize output container and screening settings
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Provisional implementation-specific screening interval retained from the
% supplied original source. It is not a formal published calibration range.
screeningLnKD_min = -3.5;
screeningLnKD_max = -0.3;

pressureRangeMessageIssued = false;
temperatureRangeMessageIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-6) Interactive selection loop and calculation
dataCodes_opx = dataset_opx{:, 1};
displayCodes_opx = cellstr(string(dataCodes_opx));

dataCodes_sp = dataset_sp{:, 1};
displayCodes_sp = cellstr(string(dataCodes_sp));

dataCodes_ol = dataset_ol{:, 1};
displayCodes_ol = cellstr(string(dataCodes_ol));

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

    % ----- Spinel selection -----
    disp('=== Step 4: Selecting a data code from the list (Spinel) ===');

    [selectedIdx_sp, ok] = listdlg( ...
        'PromptString', ...
            'Please select the Spinel data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_sp, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_sp)
        disp('Selection canceled');
        break;
    end

    selectedCode_sp = string(dataCodes_sp(selectedIdx_sp));
    selectedData_sp = dataset_sp(selectedIdx_sp, :);
    disp(['Spinel selected: ' char(selectedCode_sp)]);

    % ----- Olivine selection -----
    disp('=== Step 5: Selecting a data code from the list (Olivine) ===');

    [selectedIdx_ol, ok] = listdlg( ...
        'PromptString', ...
            'Please select the Olivine data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', displayCodes_ol, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_ol)
        disp('Selection canceled');
        break;
    end

    selectedCode_ol = string(dataCodes_ol(selectedIdx_ol));
    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    disp(['Olivine selected: ' char(selectedCode_ol)]);

    % ----- Calculation -----
    disp('=== Step 6: Calculating the temperature ===');

    nanInputNames = findNaNInputs( ...
        selectedData_opx, selectedData_sp, selectedData_ol);

    validateNonNegativeInputs( ...
        selectedData_opx, selectedData_sp, selectedData_ol);

    row = calcTemp( ...
        selectedData_opx, selectedData_sp, selectedData_ol, P_kbar);

    nRows = height(row);
    row.dataCode_opx = repmat(selectedCode_opx, nRows, 1);
    row.dataCode_sp = repmat(selectedCode_sp, nRows, 1);
    row.dataCode_ol = repmat(selectedCode_ol, nRows, 1);
    row = movevars( ...
        row, ...
        {'dataCode_opx', 'dataCode_sp', 'dataCode_ol'}, ...
        'Before', 1);

    % Store one completed result block. The complete output table is not
    % enlarged on every interactive iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = ...
            [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % ----- Immediate result display -----
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    printTemperatureSummary( ...
        selectedCode_opx, selectedCode_sp, selectedCode_ol, ...
        row.T_degreeC);

    % Equation (6) has no published numerical pressure interval.
    if ~pressureRangeMessageIssued
        fprintf(2, ...
            ['CAUTION: Witt-Eickschen and Seck (1991) do not state a ' ...
             'numerical pressure calibration range specifically for ' ...
             'Equation (6), and the equation contains no pressure term ' ...
             '(p. 431; article pp. 431-439). P_kbar is stored for ' ...
             'interface compatibility but is not used. Input range = ' ...
             '%.4g-%.4g kbar. A published pressure-range validity test ' ...
             'cannot be performed.\n'], ...
            min(P_kbar), ...
            max(P_kbar));
        pressureRangeMessageIssued = true;
    end

    % Equation (6) has no single formal numerical temperature interval.
    if ~temperatureRangeMessageIssued
        fprintf(2, ...
            ['CAUTION: Witt-Eickschen and Seck (1991) do not state a ' ...
             'single formal numerical temperature calibration interval ' ...
             'specifically for Equation (6) (p. 431; article ' ...
             'pp. 431-439). No invented numerical temperature limit has ' ...
             'been applied. Mineral assemblage, equilibrium, and ' ...
             'composition must be evaluated directly.\n']);
        temperatureRangeMessageIssued = true;
    end

    printLnKDScreeningWarning( ...
        row.lnKD, ...
        screeningLnKD_min, ...
        screeningLnKD_max, ...
        selectedCode_opx, selectedCode_sp, selectedCode_ol);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Equation (6) input(s) for ' ...
             '%s & %s & %s: %s.\n' ...
             '         Existing NaN values were retained and propagated; ' ...
             'the calculated temperature may be NaN.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_sp), ...
            char(selectedCode_ol), ...
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
            ['WARNING: Invalid Witt-Eickschen-Seck Equation (6) site, ' ...
             'fraction, KD, logarithm, or temperature term(s) were found ' ...
             'for %s & %s & %s: %s.\n' ...
             '         Affected values were retained as NaN, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_sp), ...
            char(selectedCode_ol), ...
            char(strjoin(invalidTerms, ', ')));
    end

    printNonFiniteTemperatureWarning( ...
        row.T_degreeC, ...
        selectedCode_opx, selectedCode_sp, selectedCode_ol);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'WittEickschenSeck1991OpxSp', ...
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
        'Witt-Eickschen and Seck (1991) Equation (6)', ...
    'calibrationBasis', ...
        ['More than 100 natural spinel peridotites from the Rhenish ' ...
         'Volcanic Province; reference temperatures from the Brey and ' ...
         'Kohler (1990) Ca-in-Opx thermometer'], ...
    'publishedKDDefinition', ...
        ['(XFo_Ol^2 * XM1Al_Opx) / ' ...
         '(XAl_Sp^2 * XMg_Sp^2 * XM1Mg_Opx)'], ...
    'pressureUsedInEquation', false, ...
    'publishedPressureCalibrationRange', ...
        'Not numerically specified for Equation (6)', ...
    'publishedTemperatureCalibrationRange', ...
        'No single formal numerical interval stated for Equation (6)', ...
    'provisionalLnKDScreening', ...
        [screeningLnKD_min, screeningLnKD_max], ...
    'screeningRangeStatus', ...
        'Implementation-specific; not a formal published calibration range');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_opx, data_sp, data_ol)
% findNaNInputs
% Return names of explicitly stored mineral variables used by Equation (6)
% or its site-allocation procedure that contain NaN. Missing optional
% columns are not reported because they receive the documented default zero.

opxVariables = { ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Fe3_cation_apfu'};

spinelVariables = { ...
    'Al_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Fe3_cation_apfu'};

olivineVariables = { ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Ni_cation_apfu', ...
    'Fe3_cation_apfu'};

maxNames = ...
    numel(opxVariables) + ...
    numel(spinelVariables) + ...
    numel(olivineVariables);

nameBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(opxVariables)
    variableName = opxVariables{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        value = data_opx.(variableName);
        if isnumeric(value) && any(isnan(value(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Opx." + string(variableName);
        end
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};

    if ismember(variableName, data_sp.Properties.VariableNames)
        value = data_sp.(variableName);
        if isnumeric(value) && any(isnan(value(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Spinel." + string(variableName);
        end
    end
end

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};

    if ismember(variableName, data_ol.Properties.VariableNames)
        value = data_ol.(variableName);
        if isnumeric(value) && any(isnan(value(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Olivine." + string(variableName);
        end
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(data_opx, data_sp, data_ol)
% validateNonNegativeInputs
% Stop when a stored mineral-composition value is negative or infinite.
% Zero is allowed. NaN is deliberately allowed and propagated.

opxVariables = { ...
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

spinelVariables = { ...
    'Al_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Fe3_cation_apfu'};

olivineVariables = { ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Ni_cation_apfu', ...
    'Fe3_cation_apfu'};

maxNames = ...
    numel(opxVariables) + ...
    numel(spinelVariables) + ...
    numel(olivineVariables);

nameBuffer = strings(maxNames, 1);
nInvalid = 0;

for i = 1:numel(opxVariables)
    variableName = opxVariables{i};

    if ismember(variableName, data_opx.Properties.VariableNames)
        value = data_opx.(variableName);
        validateScalarVariable(value, 'Opx', variableName);

        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            nameBuffer(nInvalid) = "Opx." + string(variableName);
        end
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};

    if ismember(variableName, data_sp.Properties.VariableNames)
        value = data_sp.(variableName);
        validateScalarVariable(value, 'Spinel', variableName);

        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            nameBuffer(nInvalid) = "Spinel." + string(variableName);
        end
    end
end

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};

    if ismember(variableName, data_ol.Properties.VariableNames)
        value = data_ol.(variableName);
        validateScalarVariable(value, 'Olivine', variableName);

        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            nameBuffer(nInvalid) = "Olivine." + string(variableName);
        end
    end
end

if nInvalid > 0
    invalidNames = nameBuffer(1:nInvalid);
    error(['WittEickschenSeck1991OpxSp: mineral-composition inputs ' ...
           'must not be negative or infinite. Invalid value(s) were ' ...
           'found in: ' char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_sp, data_ol, P_kbar)
% calcTemp
% Calculate Witt-Eickschen and Seck (1991) Equation (6) for one selected
% Opx-Sp-Ol triplet and repeat the pressure-independent result for every
% supplied pressure. Existing NaN values and invalid derived terms are
% retained as NaN.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

opx = prepareOpxRow(data_opx, 'Opx');
sp = prepareSpinelRow(data_sp, 'Spinel');
ol = prepareOlivineRow(data_ol, 'Olivine');

site_opx = calcOpxSiteFractions(opx);
spinelTerms = calcSpinelTerms(sp);
olivineTerms = calcOlivineTerms(ol);

% Published KD definition:
%   KD = (XFo_Ol^2 * XM1Al_Opx) /
%        (XAl_Sp^2 * XMg_Sp^2 * XM1Mg_Opx)
KD_numerator = NaN;
KD_denominator = NaN;
KD = NaN;
lnKD = NaN;

if all(isfinite([ ...
        olivineTerms.XFo_Ol, ...
        site_opx.XM1Al_Opx, ...
        spinelTerms.XAl_Sp, ...
        spinelTerms.XMg_Sp, ...
        site_opx.XM1Mg_Opx]))

    KD_numerator = ...
        (olivineTerms.XFo_Ol .^ 2) .* ...
        site_opx.XM1Al_Opx;

    KD_denominator = ...
        (spinelTerms.XAl_Sp .^ 2) .* ...
        (spinelTerms.XMg_Sp .^ 2) .* ...
        site_opx.XM1Mg_Opx;

    if isfinite(KD_numerator) && KD_numerator > 0 && ...
            isfinite(KD_denominator) && KD_denominator > 0

        KD = KD_numerator ./ KD_denominator;

        if isfinite(KD) && KD > 0
            lnKD = log(KD);
        end
    end
end

% Witt-Eickschen and Seck (1991), Equation (6).
if isfinite(lnKD) && isfinite(spinelTerms.YCr_Sp)
    T_scalar_raw_degreeC = ...
        2248.25 ...
        + 991.58 .* lnKD ...
        + 153.32 .* (lnKD .^ 2) ...
        + 539.05 .* spinelTerms.YCr_Sp ...
        - 2005.74 .* (spinelTerms.YCr_Sp .^ 2);
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
    repmat("WittEickschenSeck1991_Eq6", nP, 1);

% Orthopyroxene composition.
row.Si_opx = repmat(opx.Si, nP, 1);
row.Al_opx = repmat(opx.Al, nP, 1);
row.Fe_total_opx = repmat(opx.Fe_total, nP, 1);
row.Fe2_opx = repmat(opx.Fe2, nP, 1);
row.Fe3_opx = repmat(opx.Fe3, nP, 1);
row.Mg_opx = repmat(opx.Mg, nP, 1);
row.Ca_opx = repmat(opx.Ca, nP, 1);
row.Na_opx = repmat(opx.Na, nP, 1);
row.Mn_opx = repmat(opx.Mn, nP, 1);
row.Ti_opx = repmat(opx.Ti, nP, 1);
row.Cr_opx = repmat(opx.Cr, nP, 1);

% Spinel composition.
row.Al_sp = repmat(sp.Al, nP, 1);
row.Cr_sp = repmat(sp.Cr, nP, 1);
row.Fe_total_sp = repmat(sp.Fe_total, nP, 1);
row.Fe2_sp = repmat(sp.Fe2, nP, 1);
row.Fe3_sp = repmat(sp.Fe3, nP, 1);
row.Mg_sp = repmat(sp.Mg, nP, 1);

% Olivine composition.
row.Mg_ol = repmat(ol.Mg, nP, 1);
row.Fe_total_ol = repmat(ol.Fe_total, nP, 1);
row.Fe2_ol = repmat(ol.Fe2, nP, 1);
row.Fe3_ol = repmat(ol.Fe3, nP, 1);
row.Mn_ol = repmat(ol.Mn, nP, 1);
row.Ca_ol = repmat(ol.Ca, nP, 1);
row.Ni_ol = repmat(ol.Ni, nP, 1);

% Orthopyroxene site allocation.
row.XM1Opx_Al_raw = ...
    repmat(site_opx.XM1Al_Opx_raw, nP, 1);
row.XM1Opx_Al = ...
    repmat(site_opx.XM1Al_Opx, nP, 1);
row.M1_fixed_opx = ...
    repmat(site_opx.M1_fixed_Opx, nP, 1);
row.M1_remaining_opx = ...
    repmat(site_opx.M1_remaining_Opx, nP, 1);
row.MgFe_total_opx = ...
    repmat(site_opx.MgFe_total_Opx, nP, 1);
row.Mg_fraction_opx = ...
    repmat(site_opx.Mg_fraction_Opx, nP, 1);
row.Fe2_fraction_opx = ...
    repmat(site_opx.Fe2_fraction_Opx, nP, 1);
row.Mg_M1_opx = ...
    repmat(site_opx.Mg_M1_Opx, nP, 1);
row.Fe2_M1_opx = ...
    repmat(site_opx.Fe2_M1_Opx, nP, 1);
row.M1_total_opx = ...
    repmat(site_opx.M1_total_Opx, nP, 1);
row.XM1Opx_Mg = ...
    repmat(site_opx.XM1Mg_Opx, nP, 1);

% Standardized aliases.
row.XM1Al_Opx = row.XM1Opx_Al;
row.XM1Mg_Opx = row.XM1Opx_Mg;

% Spinel terms.
row.denom_trivalent_sp = ...
    repmat(spinelTerms.denom_trivalent_Sp, nP, 1);
row.denom_divalent_sp = ...
    repmat(spinelTerms.denom_divalent_Sp, nP, 1);
row.XSp_Al = repmat(spinelTerms.XAl_Sp, nP, 1);
row.XSp_Mg = repmat(spinelTerms.XMg_Sp, nP, 1);
row.YCr_sp = repmat(spinelTerms.YCr_Sp, nP, 1);

% Standardized aliases.
row.XAl_Sp = row.XSp_Al;
row.XMg_Sp = row.XSp_Mg;
row.YCr_Sp = row.YCr_sp;

% Olivine term.
row.denom_divalent_ol = ...
    repmat(olivineTerms.denom_divalent_Ol, nP, 1);
row.XOl_Fo = repmat(olivineTerms.XFo_Ol, nP, 1);
row.XFo_Ol = row.XOl_Fo;

% KD, logarithm, and temperature.
row.KD_numerator = repmat(KD_numerator, nP, 1);
row.KD_denominator = repmat(KD_denominator, nP, 1);
row.KD = repmat(KD, nP, 1);
row.lnKD = repmat(lnKD, nP, 1);

row.T_Eq6_raw_degreeC = T_raw_degreeC;
row.T_K = T_K;
row.T_degreeC = T_degreeC;
row.T_deg = T_degreeC;

end

function opx = prepareOpxRow(data_opx, mineralLabel)
% prepareOpxRow
% Extract one Opx composition. Existing NaN values remain NaN. Missing
% optional columns are assigned zero. Fe_cation_apfu is treated as total Fe.

if height(data_opx) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

opx = struct();

opx.Si = getVarRequired( ...
    data_opx, 'Si_cation_apfu', mineralLabel);
opx.Al = getVarRequired( ...
    data_opx, 'Al_cation_apfu', mineralLabel);
opx.Fe_total = getVarRequired( ...
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

opx.Fe2 = deriveFe2(opx.Fe_total, opx.Fe3, mineralLabel);

end

function sp = prepareSpinelRow(data_sp, mineralLabel)
% prepareSpinelRow
% Extract one Spinel composition. Existing NaN values remain NaN. Missing
% Fe3 is assigned zero. Fe_cation_apfu is treated as total Fe.

if height(data_sp) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

sp = struct();

sp.Al = getVarRequired( ...
    data_sp, 'Al_cation_apfu', mineralLabel);
sp.Cr = getVarRequired( ...
    data_sp, 'Cr_cation_apfu', mineralLabel);
sp.Mg = getVarRequired( ...
    data_sp, 'Mg_cation_apfu', mineralLabel);
sp.Fe_total = getVarRequired( ...
    data_sp, 'Fe_cation_apfu', mineralLabel);
sp.Fe3 = getVarOptional( ...
    data_sp, 'Fe3_cation_apfu', 0, mineralLabel);

sp.Fe2 = deriveFe2(sp.Fe_total, sp.Fe3, mineralLabel);

end

function ol = prepareOlivineRow(data_ol, mineralLabel)
% prepareOlivineRow
% Extract one Olivine composition. Existing NaN values remain NaN. Missing
% optional columns are assigned zero. Fe_cation_apfu is treated as total Fe.

if height(data_ol) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

ol = struct();

ol.Fe_total = getVarRequired( ...
    data_ol, 'Fe_cation_apfu', mineralLabel);
ol.Mg = getVarRequired( ...
    data_ol, 'Mg_cation_apfu', mineralLabel);

ol.Mn = getVarOptional( ...
    data_ol, 'Mn_cation_apfu', 0, mineralLabel);
ol.Ca = getVarOptional( ...
    data_ol, 'Ca_cation_apfu', 0, mineralLabel);
ol.Ni = getVarOptional( ...
    data_ol, 'Ni_cation_apfu', 0, mineralLabel);
ol.Fe3 = getVarOptional( ...
    data_ol, 'Fe3_cation_apfu', 0, mineralLabel);

ol.Fe2 = deriveFe2(ol.Fe_total, ol.Fe3, mineralLabel);

end

function Fe2 = deriveFe2(Fe_total, Fe3, mineralLabel)
% deriveFe2
% Derive Fe2 from total Fe and Fe3. NaN is propagated. A finite negative
% derived Fe2 value is prohibited.

Fe2 = Fe_total - Fe3;

if isfinite(Fe2) && Fe2 < 0
    error(['%s derived Fe2_cation_apfu is negative because ' ...
           'Fe3_cation_apfu exceeds Fe_cation_apfu.'], ...
        mineralLabel);
end
if isinf(Fe2)
    error('%s derived Fe2_cation_apfu is infinite.', mineralLabel);
end

end

function site = calcOpxSiteFractions(opx)
% calcOpxSiteFractions
% Calculate Opx M1 quantities used by Equation (6). Invalid derived terms
% become NaN so that the output row can be retained and diagnosed.

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

if all(isfinite([XM1Al_Opx, opx.Cr, opx.Ti, opx.Fe3]))
    M1_fixed_raw = XM1Al_Opx + opx.Cr + opx.Ti + opx.Fe3;

    if M1_fixed_raw >= 0 && M1_fixed_raw <= 1
        M1_fixed_Opx = M1_fixed_raw;
        M1_remaining_Opx = 1 - M1_fixed_Opx;
    else
        M1_fixed_Opx = NaN;
        M1_remaining_Opx = NaN;
    end
else
    M1_fixed_Opx = NaN;
    M1_remaining_Opx = NaN;
end

MgFe_total_Opx = opx.Mg + opx.Fe2;

if isfinite(MgFe_total_Opx) && MgFe_total_Opx > 0 && ...
        isfinite(opx.Mg) && isfinite(opx.Fe2) && ...
        isfinite(M1_remaining_Opx)

    Mg_fraction_Opx = opx.Mg ./ MgFe_total_Opx;
    Fe2_fraction_Opx = opx.Fe2 ./ MgFe_total_Opx;

    Mg_M1_Opx = M1_remaining_Opx .* Mg_fraction_Opx;
    Fe2_M1_Opx = M1_remaining_Opx .* Fe2_fraction_Opx;

    M1_total_Opx = ...
        M1_fixed_Opx + Mg_M1_Opx + Fe2_M1_Opx;

    if isfinite(M1_total_Opx) && M1_total_Opx > 0
        XM1Mg_Opx = Mg_M1_Opx ./ M1_total_Opx;
    else
        M1_total_Opx = NaN;
        XM1Mg_Opx = NaN;
    end
else
    Mg_fraction_Opx = NaN;
    Fe2_fraction_Opx = NaN;
    Mg_M1_Opx = NaN;
    Fe2_M1_Opx = NaN;
    M1_total_Opx = NaN;
    XM1Mg_Opx = NaN;
end

site.XM1Al_Opx_raw = XM1Al_Opx_raw;
site.XM1Al_Opx = XM1Al_Opx;
site.M1_fixed_Opx = M1_fixed_Opx;
site.M1_remaining_Opx = M1_remaining_Opx;
site.MgFe_total_Opx = MgFe_total_Opx;
site.Mg_fraction_Opx = Mg_fraction_Opx;
site.Fe2_fraction_Opx = Fe2_fraction_Opx;
site.Mg_M1_Opx = Mg_M1_Opx;
site.Fe2_M1_Opx = Fe2_M1_Opx;
site.M1_total_Opx = M1_total_Opx;
site.XM1Mg_Opx = XM1Mg_Opx;

end

function terms = calcSpinelTerms(sp)
% calcSpinelTerms
% Calculate Spinel fractions used by Equation (6). Invalid denominators
% produce NaN rather than stopping the complete calculation.

terms = struct();

denom_trivalent_Sp = sp.Al + sp.Cr + sp.Fe3;
denom_divalent_Sp = sp.Mg + sp.Fe2;

if isfinite(denom_trivalent_Sp) && denom_trivalent_Sp > 0 && ...
        all(isfinite([sp.Al, sp.Cr]))
    XAl_Sp = sp.Al ./ denom_trivalent_Sp;
    YCr_Sp = sp.Cr ./ denom_trivalent_Sp;
else
    XAl_Sp = NaN;
    YCr_Sp = NaN;
end

if isfinite(denom_divalent_Sp) && denom_divalent_Sp > 0 && ...
        isfinite(sp.Mg)
    XMg_Sp = sp.Mg ./ denom_divalent_Sp;
else
    XMg_Sp = NaN;
end

if ~(isfinite(XAl_Sp) && XAl_Sp >= 0 && XAl_Sp <= 1)
    XAl_Sp = NaN;
end
if ~(isfinite(XMg_Sp) && XMg_Sp >= 0 && XMg_Sp <= 1)
    XMg_Sp = NaN;
end
if ~(isfinite(YCr_Sp) && YCr_Sp >= 0 && YCr_Sp <= 1)
    YCr_Sp = NaN;
end

terms.denom_trivalent_Sp = denom_trivalent_Sp;
terms.denom_divalent_Sp = denom_divalent_Sp;
terms.XAl_Sp = XAl_Sp;
terms.XMg_Sp = XMg_Sp;
terms.YCr_Sp = YCr_Sp;

end

function terms = calcOlivineTerms(ol)
% calcOlivineTerms
% Calculate the Olivine forsterite component used by Equation (6). Invalid
% denominators produce NaN rather than stopping the complete calculation.

terms = struct();

denom_divalent_Ol = ...
    ol.Mg + ol.Fe2 + ol.Mn + ol.Ca + ol.Ni;

if isfinite(denom_divalent_Ol) && denom_divalent_Ol > 0 && ...
        isfinite(ol.Mg)

    XFo_Ol = ol.Mg ./ denom_divalent_Ol;
else
    XFo_Ol = NaN;
end

if ~(isfinite(XFo_Ol) && XFo_Ol >= 0 && XFo_Ol <= 1)
    XFo_Ol = NaN;
end

terms.denom_divalent_Ol = denom_divalent_Ol;
terms.XFo_Ol = XFo_Ol;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid site-allocation, fraction, KD, logarithm, and temperature
% terms. Logarithms may be negative; only finiteness is required.

termBuffer = strings(24, 1);
nTerms = 0;

positiveChecks = { ...
    'XM1Al_Opx', row.XM1Opx_Al(1); ...
    'XM1Mg_Opx', row.XM1Opx_Mg(1); ...
    'XAl_Sp', row.XSp_Al(1); ...
    'XMg_Sp', row.XSp_Mg(1); ...
    'XFo_Ol', row.XOl_Fo(1); ...
    'KD numerator', row.KD_numerator(1); ...
    'KD denominator', row.KD_denominator(1); ...
    'KD', row.KD(1)};

for i = 1:size(positiveChecks, 1)
    label = positiveChecks{i, 1};
    value = positiveChecks{i, 2};

    if ~isfinite(value) || value <= 0
        nTerms = nTerms + 1;
        termBuffer(nTerms) = string(label);
    end
end

boundedChecks = { ...
    'XM1Al_Opx raw', row.XM1Opx_Al_raw(1), 0, 1; ...
    'M1 fixed Opx', row.M1_fixed_opx(1), 0, 1; ...
    'M1 remaining Opx', row.M1_remaining_opx(1), 0, 1; ...
    'YCr_Sp', row.YCr_sp(1), 0, 1};

for i = 1:size(boundedChecks, 1)
    label = boundedChecks{i, 1};
    value = boundedChecks{i, 2};
    minimumValue = boundedChecks{i, 3};
    maximumValue = boundedChecks{i, 4};

    if ~isfinite(value) || ...
            value < minimumValue || value > maximumValue

        nTerms = nTerms + 1;
        termBuffer(nTerms) = string(label);
    end
end

if ~isfinite(row.lnKD(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "lnKD";
end

if any(~isfinite(row.T_K) | row.T_K <= 0)
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "T_K";
end

invalidTerms = unique(termBuffer(1:nTerms), 'stable');

end

function printTemperatureSummary( ...
        selectedCode_opx, selectedCode_sp, selectedCode_ol, ...
        temperatureValues)
% printTemperatureSummary
% Display one temperature or first-to-last values for a pressure vector.

label = [ ...
    char(selectedCode_opx) ' & ' ...
    char(selectedCode_sp) ' & ' ...
    char(selectedCode_ol)];

if numel(temperatureValues) == 1
    disp([label ': ' num2str(temperatureValues) ' degreeC']);
else
    disp([label ': ' ...
        num2str(temperatureValues(1)) ' to ' ...
        num2str(temperatureValues(end)) ' degreeC']);
end

end

function printLnKDScreeningWarning( ...
        lnKDValues, minimumLnKD, maximumLnKD, ...
        selectedCode_opx, selectedCode_sp, selectedCode_ol)
% printLnKDScreeningWarning
% Warn when finite lnKD lies outside the provisional implementation-specific
% interval retained from the supplied source. Results are retained.

finiteMask = isfinite(lnKDValues);
outsideMask = finiteMask & ...
    (lnKDValues < minimumLnKD | lnKDValues > maximumLnKD);

if any(outsideMask)
    finiteValues = lnKDValues(finiteMask);

    fprintf(2, ...
        ['WARNING: lnKD is outside the provisional implementation-specific ' ...
         'screening interval %.4g to %.4g. %d of %d finite point(s) are ' ...
         'outside; calculated finite lnKD range = %.6g to %.6g for ' ...
         '%s & %s & %s. This interval is not presented as a formal ' ...
         'published calibration range. The result has been retained.\n'], ...
        minimumLnKD, ...
        maximumLnKD, ...
        sum(outsideMask), ...
        sum(finiteMask), ...
        min(finiteValues), ...
        max(finiteValues), ...
        char(selectedCode_opx), ...
        char(selectedCode_sp), ...
        char(selectedCode_ol));
end

end

function printNonFiniteTemperatureWarning( ...
        temperatureValues, ...
        selectedCode_opx, selectedCode_sp, selectedCode_ol)
% printNonFiniteTemperatureWarning
% Report NaN and Inf temperatures without modifying or deleting them.

invalidMask = ~isfinite(temperatureValues);

if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Witt-Eickschen-Seck Equation (6) ' ...
         'temperature values were calculated for %s & %s & %s ' ...
         '(%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the ' ...
         'calculation has not been stopped.\n'], ...
        char(selectedCode_opx), ...
        char(selectedCode_sp), ...
        char(selectedCode_ol), ...
        sum(invalidMask), ...
        numel(temperatureValues), ...
        sum(isnan(temperatureValues)), ...
        sum(isinf(temperatureValues)));
end

end

function validateRequiredColumns(tbl, mineralLabel)
% validateRequiredColumns
% Verify normalized cation columns required for each mineral.

switch mineralLabel
    case 'Opx'
        requiredVariables = { ...
            'Si_cation_apfu', ...
            'Al_cation_apfu', ...
            'Fe_cation_apfu', ...
            'Mg_cation_apfu', ...
            'Ca_cation_apfu', ...
            'Na_cation_apfu', ...
            'Ti_cation_apfu', ...
            'Cr_cation_apfu'};

    case 'Spinel'
        requiredVariables = { ...
            'Al_cation_apfu', ...
            'Cr_cation_apfu', ...
            'Mg_cation_apfu', ...
            'Fe_cation_apfu'};

    case 'Olivine'
        requiredVariables = { ...
            'Fe_cation_apfu', ...
            'Mg_cation_apfu'};

    otherwise
        error('Unknown mineral label: %s', mineralLabel);
end

for i = 1:numel(requiredVariables)
    variableName = requiredVariables{i};

    if ~ismember(variableName, tbl.Properties.VariableNames)
        error('%s table must contain variable: %s', ...
            mineralLabel, variableName);
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
