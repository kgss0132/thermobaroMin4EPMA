function results = WoodBanno1973baro(rawdata_struct, T_degreeC)
% functions/+baro/+Pyroxene/WoodBanno1973baro.m
% Tested with MATLAB R2024b
%
% Garnet-Orthopyroxene barometer
% Wood, B.J. and Banno, S. (1973)
% Contributions to Mineralogy and Petrology, 42, 109-124
% DOI: https://doi.org/10.1007/BF00371501
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis and one
% Garnet analysis selected from input tables and calculates pressure using
% the multicomponent garnet-orthopyroxene formulation of Wood and Banno
% (1973), Eq. (17).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Opx-Garnet pair, one output row is
% returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Wood and Banno (1973) derived the barometer by extrapolating experimental
% equilibria in the simple MgSiO3-Al2O3 system to natural Ca-Fe-Mg-Al
% garnet-orthopyroxene assemblages using simplified solution models
% (pp. 109-114). The method requires an independently known temperature
% (p. 109 and p. 113).
%
% The direct comparisons with multicomponent experiments listed in Table 2
% cover approximately:
%
%   Temperature : 1200-1400 degreeC
%   Pressure    : 25.9-40.5 kbar
%   Opx Al2O3   : 2.2-5.5 wt%
%   XMg_garnet  : 0.67-0.85
%   XAl_M1_opx  : 0.046-0.115
%   XMg_M1_opx  : 0.84-0.93
%   XMg_M2_opx  : 0.76-0.95
%
% These values are taken from Table 2 on p. 116 and are used here as the
% direct experimental verification envelope. Figure 1 on p. 113 shows a
% substantially broader modeled P-T-composition space, but that figure is an
% extrapolation of the simple-system model and must not be interpreted as a
% validated calibration range.
%
% IMPORTANT SYSTEMATIC OFFSET:
% Wood and Banno (1973) state that pressure estimates based on the Boyd and
% England simple-system data may be approximately 3-5 kbar too high and that
% the calculated pressures in Table 2 should be lowered by this amount
% (p. 116). This implementation reports the uncorrected Eq. (17) result and
% issues a caution. No automatic subtraction is applied.
%
% IMPORTANT MODEL ASSUMPTIONS:
%   1) Garnet and orthopyroxene solution behavior is represented by simple
%      idealized models (pp. 110-114).
%   2) Fe-rich orthopyroxenes may deviate from the ideal two-site treatment,
%      increasing uncertainty (p. 114).
%   3) Minor Cr, Fe3+, and Ti and analytical uncertainty in Si affect the
%      assignment of Al between tetrahedral and M1 sites. The authors
%      recommend considering the average of tetrahedral and M1 Al estimates
%      when they differ (p. 115).
%   4) Thermal expansion and compressibility of the volume terms were
%      neglected; 298 K and 1 bar values were used (p. 112).
%   5) Table 1 provides the volume term only for 2-30% Al occupancy of Opx
%      M1 (p. 112). Values outside this range require extrapolation.
%   6) The original script's random Fe-Mg distribution over available M1
%      and M2 capacity is retained here. Wood and Banno discuss
%      temperature-dependent intracrystalline Fe-Mg ordering on pp. 113-114;
%      this simplification should be considered when interpreting results,
%      especially for Fe-rich or low-temperature orthopyroxene.
%
% Garnet and orthopyroxene must be a chemically equilibrated pair. Zoning,
% exsolution, alteration, inherited grains, or non-equilibrium pairing can
% produce misleading pressures.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 1200-1400 degreeC,
%   2) finite calculated pressure is outside 25.9-40.5 kbar,
%   3) Opx M1 Al is outside the 2-30% Table 1 volume-data range,
%   4) required calculation input contains NaN,
%   5) site allocation is invalid or incomplete, or
%   6) calculated pressure is NaN, Inf, or negative.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx                    : table
%   rawdata_struct.Garnet or .Grt         : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% Required normalized-cation variables in both mineral tables:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu                       % Fe2+ used directly in the equation
%   Mg_cation_apfu
%   Ca_cation_apfu
%
% Optional normalized-cation variables:
%   Na_cation_apfu
%   Mn_cation_apfu
%   Ti_cation_apfu
%   Cr_cation_apfu
%   Fe3_cation_apfu                   % independent Fe3+ value; not subtracted
%                                      from Fe_cation_apfu
%
% IMPORTANT Fe convention: Fe_cation_apfu is used directly as Fe2+.
% Fe3_cation_apfu is a separate optional value and is never subtracted
% from Fe_cation_apfu in this implementation.
%
% Optional variables that are absent are represented by NaN, not by zero.
% Because several optional Opx variables and Garnet Mn/Fe3 enter the site or
% activity calculation, their absence may propagate NaN to the result.
%
% Finite calculation inputs must be non-negative. NaN is allowed, retained,
% propagated through the calculation, and reported by fprintf warnings.
% Inf and finite negative values are prohibited. Zero is retained; if zero
% makes a site fraction, ratio, or logarithm undefined, the resulting
% NaN/Inf is retained and reported.
%
% No liquid composition is used by this barometer. Therefore, exclusion of
% Liq F and Cl from cationTotal_liq and from NaN-input warnings is not
% applicable to this function.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
% Wood and Banno (1973), Eq. (17), p. 114:
%
%   P = 1 + [
%       R*T*ln( XMg_M1_opx * XMg_M2_opx^2 * XAl_M1_opx
%               / XMg_garnet^3 )
%       + 4207 - 2.69*T
%       ] / dV
%
% where:
%   T                : temperature in Kelvin
%   R                : 1.987 cal mol^-1 K^-1
%   XMg_M1_opx       : Mg fraction on Opx M1
%   XMg_M2_opx       : Mg fraction on Opx M2
%   XAl_M1_opx       : Al fraction on Opx M1
%   XMg_garnet       : Mg / (Mg + Fe2 + Ca + Mn) in garnet
%   dV               : volume term interpolated from Table 1, p. 112
%
% Table 1 volume data:
%
%   Al in Opx M1 (%) :   2      5      10      15      20      30
%   dV (cm3/mol)     : -7.90  -8.05  -8.33   -8.57   -8.80   -9.15
%
% The numerator is in cal/mol and dV is in cm3/mol. The conversion
% 1 cal = 41.84 cm3 bar is applied, and pressure is reported in kbar.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = WoodBanno1973baro(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing Opx and Garnet/Grt tables
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Opx-Garnet pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('WoodBanno1973baro requires (rawdata_struct, T_degreeC).');
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

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The source tables are
% not modified; selected rows are read during each calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end

if isfield(rawdata_struct, 'Garnet') && istable(rawdata_struct.Garnet)
    dataset_grt = rawdata_struct.Garnet;
elseif isfield(rawdata_struct, 'Grt') && istable(rawdata_struct.Grt)
    dataset_grt = rawdata_struct.Grt;
else
    error('rawdata_struct must contain table: rawdata_struct.Garnet or rawdata_struct.Grt');
end

dataset_opx = rawdata_struct.Opx;

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu'};

validateRequiredVariables(dataset_opx, requiredVariables, 'Opx');
validateRequiredVariables(dataset_grt, requiredVariables, 'Garnet');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Direct multicomponent experimental verification envelope from Table 2.
experimentalT_min_degreeC = 1200;
experimentalT_max_degreeC = 1400;
experimentalP_min_kbar = 25.9;
experimentalP_max_kbar = 40.5;

temperatureOutsideExperimentalRange = isfinite(T_degreeC) & ...
    (T_degreeC < experimentalT_min_degreeC | ...
     T_degreeC > experimentalT_max_degreeC);

temperatureWarningIssued = false;
systematicOffsetCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Opx) ===');

while true
    % ----- Opx selection -----
    % The first table column is used only as the displayed data identifier.
    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Opx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Opx selected: ' char(string(selectedCode_opx))]);

    % ----- Garnet selection -----
    disp('=== Step 4: Selecting a data code from the list (Garnet) ===');

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

    % ----- Calculation -----
    % Opx and Garnet are selected independently; row numbers are not assumed
    % to correspond between the two tables.
    disp('=== Step 5: Calculating the pressure ===');

    selectedData_opx = dataset_opx(selectedIdx_opx, :);
    selectedData_grt = dataset_grt(selectedIdx_grt, :);

    % Check NaN only in variables directly used by the pressure calculation.
    % Missing optional variables are treated as NaN and reported.
    nanInputNames = findNaNInputs( ...
        selectedData_opx, selectedData_grt, T_degreeC);

    % Reject only Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(selectedData_opx, selectedData_grt);

    row = calcPressure(selectedData_opx, selectedData_grt, T_degreeC);

    % Repeat identifiers for all temperatures in the current calculation.
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row.dataCode_garnet = repmat(string(selectedCode_grt), height(row), 1);
    row = movevars(row, {'dataCode_opx','dataCode_garnet'}, 'Before', 1);

    % Store one block per selected mineral pair. Expand the cell buffer only
    % when its current capacity has been exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed finite pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    finitePressureValues = row.P_kbar(isfinite(row.P_kbar));
    if isempty(finitePressureValues)
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_grt)) ...
            ': P = all values are NaN or Inf']);
    elseif numel(row.P_kbar) == 1
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_grt)) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_opx)) ' & ' ...
            char(string(selectedCode_grt)) ': finite P range = ' ...
            num2str(min(finitePressureValues)) ' to ' ...
            num2str(max(finitePressureValues)) ' kbar']);
    end

    % Print the systematic-offset caution once per function call.
    if ~systematicOffsetCautionIssued
        fprintf(2, ...
            ['CAUTION: Wood and Banno (1973) state that pressures based on the ' ...
             'Boyd and England simple-system data may be approximately 3-5 kbar ' ...
             'too high and that the Table 2 calculated pressures should be lowered ' ...
             'by this amount (p. 116). The uncorrected Eq. (17) result is retained.\n']);
        systematicOffsetCautionIssued = true;
    end

    % Input temperature is common to all selected pairs, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideExperimentalRange) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the direct multicomponent ' ...
             'experimental verification range of Wood and Banno (1973): ' ...
             '1200-1400 degreeC (Table 2, p. 116). %d of %d finite temperature ' ...
             'point(s) are outside the range; input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideExperimentalRange), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures lie outside the direct
    % multicomponent experimental verification range.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideExperimentalRange = finitePressure & ...
        (row.P_kbar < experimentalP_min_kbar | ...
         row.P_kbar > experimentalP_max_kbar);

    if any(pressureOutsideExperimentalRange)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the direct multicomponent ' ...
             'experimental verification range of Wood and Banno (1973): ' ...
             '25.9-40.5 kbar (Table 2, p. 116). %d of %d finite pressure point(s) ' ...
             'are outside the range; calculated finite range = %.4g-%.4g kbar ' ...
             'for %s & %s.\n'], ...
            sum(pressureOutsideExperimentalRange), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)));
    end

    % Warn when Table 1 volume data must be extrapolated.
    finiteAlM1 = isfinite(row.AlM1_percent_opx(1));
    if finiteAlM1 && ...
            (row.AlM1_percent_opx(1) < 2 || row.AlM1_percent_opx(1) > 30)
        fprintf(2, ...
            ['WARNING: Opx M1 Al = %.4g%% is outside the 2-30%% range of the ' ...
             'volume data in Wood and Banno (1973), Table 1, p. 112. The dV term ' ...
             'was linearly extrapolated for %s & %s.\n'], ...
            row.AlM1_percent_opx(1), ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)));
    end

    % Report incomplete or physically inconsistent site allocation without
    % stopping the calculation. Invalid derived values remain NaN.
    if ~row.isOpxSiteAllocationValid(1)
        fprintf(2, ...
            ['WARNING: Opx site allocation is incomplete or invalid for %s & %s. ' ...
             'Possible causes include NaN inputs, AlVI < 0, fixed M1 or M2 ' ...
             'occupancy > 1, or Mg + Fe2 <= 0. Derived NaN values were retained.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)));
    end

    % List the exact calculation inputs containing NaN. For vector
    % temperature input, the NaN element indices are included.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s & %s: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'pressure may remain NaN.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressure is retained for diagnostics but is outside
    % the experimental range and physical/useful domain.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & %s ' ...
             '(%d of %d points). The values were retained for diagnostic purposes.\n'], ...
            char(string(selectedCode_opx)), ...
            char(string(selectedCode_grt)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'WoodBanno1973baro', ...
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

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(tbl, requiredVariables, mineralLabel)
% validateRequiredVariables
% Confirm that required table variables exist before interactive selection.

missingMask = ~ismember(requiredVariables, tbl.Properties.VariableNames);
if any(missingMask)
    missingVariables = string(requiredVariables(missingMask));
    error('%s table is missing required variable(s): %s', ...
        mineralLabel, char(strjoin(missingVariables, ', ')));
end

end

function nanInputNames = findNaNInputs(data_opx, data_grt, T_degreeC)
% findNaNInputs
% Return names of pressure-equation inputs containing NaN. Missing optional
% variables are represented by NaN. F and Cl are not relevant because no
% liquid composition is used.

maxNames = 16;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

opxVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Fe3_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

for i = 1:numel(opxVariables)
    variableName = opxVariables{i};
    variableValue = getVarOrNaN(data_opx, variableName);
    if isnan(variableValue)
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Opx." + string(variableName);
    end
end

grtVariables = { ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Mn_cation_apfu'};

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    variableValue = getVarOrNaN(data_grt, variableName);
    if isnan(variableValue)
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = "Garnet." + string(variableName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_opx, data_grt)
% validateNonNegativeInputs
% Reject Inf and finite negative values in variables directly used by the
% pressure calculation. Zero and NaN are intentionally allowed and retained.

opxVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Fe3_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu'};

grtVariables = { ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Mn_cation_apfu'};

maxNames = numel(opxVariables) + numel(grtVariables);
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

for i = 1:numel(opxVariables)
    variableName = opxVariables{i};
    variableValue = getVarOrNaN(data_opx, variableName);

    if isinf(variableValue) || ...
            (isfinite(variableValue) && variableValue < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Opx." + string(variableName);
    end
end

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    variableValue = getVarOrNaN(data_grt, variableName);

    if isinf(variableValue) || ...
            (isfinite(variableValue) && variableValue < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Garnet." + string(variableName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['WoodBanno1973baro: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_opx, data_grt, T_degreeC)
% calcPressure
% Compute pressure for one Opx row and one Garnet row at one or more input
% temperatures. NaN values propagate through the calculation.
%
% Inputs:
%   data_opx   : 1-row Opx table
%   data_grt   : 1-row Garnet table
%   T_degreeC : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Physical and conversion constants.
R_cal = 1.987;
calToCm3Bar = 41.84;

% Extract one-row cation data. Missing optional variables are NaN and are
% never replaced by zero.
opx = prepareMineralRow(data_opx, 'Opx');
grt = prepareMineralRow(data_grt, 'Garnet');

% Opx site allocation. The retained implementation distributes Fe2 and Mg
% randomly over the available M1 and M2 capacities.
site_opx = calcOpxSiteFractions(opx);

XMg_M1_scalar = site_opx.XMg_M1;
XMg_M2_scalar = site_opx.XMg_M2;
XAl_M1_scalar = site_opx.XAl_M1;

% Garnet Mg fraction on the divalent site. A NaN or zero denominator
% naturally produces NaN/Inf and is not replaced.
garnetDivalentSum_scalar = grt.Mg + grt.Fe2 + grt.Ca + grt.Mn;
XMg_garnet_scalar = grt.Mg ./ garnetDivalentSum_scalar;

% Volume term from Wood and Banno (1973), Table 1, p. 112. Linear
% extrapolation is retained outside the tabulated 2-30% range and is
% reported by a warning in the main function.
AlM1_percent_scalar = 100 .* XAl_M1_scalar;
AlM1_percent_grid = [2 5 10 15 20 30];
dV_grid = [-7.90 -8.05 -8.33 -8.57 -8.80 -9.15];
dV_cm3mol_scalar = interp1( ...
    AlM1_percent_grid, dV_grid, AlM1_percent_scalar, ...
    'linear', 'extrap');

% Eq. (17) activity term. Evaluate the logarithm only for a finite,
% strictly positive value to prevent complex output. Invalid values are
% represented by NaN and propagated.
Kterm_scalar = ...
    (XMg_M1_scalar .* (XMg_M2_scalar.^2) .* XAl_M1_scalar) ./ ...
    (XMg_garnet_scalar.^3);

lnKterm_scalar = NaN;
if isfinite(Kterm_scalar) && Kterm_scalar > 0
    lnKterm_scalar = log(Kterm_scalar);
end

% Expand composition-dependent scalars to the temperature-vector length.
R_output = repmat(R_cal, nT, 1);
calToCm3Bar_output = repmat(calToCm3Bar, nT, 1);

Si_opx = repmat(opx.Si, nT, 1);
Al_opx = repmat(opx.Al, nT, 1);
FeT_opx = repmat(opx.FeT, nT, 1);
Fe2_opx = repmat(opx.Fe2, nT, 1);
Fe3_opx = repmat(opx.Fe3, nT, 1);
Mg_opx = repmat(opx.Mg, nT, 1);
Ca_opx = repmat(opx.Ca, nT, 1);
Na_opx = repmat(opx.Na, nT, 1);
Mn_opx = repmat(opx.Mn, nT, 1);
Ti_opx = repmat(opx.Ti, nT, 1);
Cr_opx = repmat(opx.Cr, nT, 1);

Si_grt = repmat(grt.Si, nT, 1);
Al_grt = repmat(grt.Al, nT, 1);
FeT_grt = repmat(grt.FeT, nT, 1);
Fe2_grt = repmat(grt.Fe2, nT, 1);
Fe3_grt = repmat(grt.Fe3, nT, 1);
Mg_grt = repmat(grt.Mg, nT, 1);
Ca_grt = repmat(grt.Ca, nT, 1);
Na_grt = repmat(grt.Na, nT, 1);
Mn_grt = repmat(grt.Mn, nT, 1);
Ti_grt = repmat(grt.Ti, nT, 1);
Cr_grt = repmat(grt.Cr, nT, 1);

AlIV_opx = repmat(site_opx.AlIV, nT, 1);
AlVI_opx = repmat(site_opx.AlVI, nT, 1);
XAl_T_opx = repmat(site_opx.XAl_T, nT, 1);
XMg_M1_opx = repmat(XMg_M1_scalar, nT, 1);
XMg_M2_opx = repmat(XMg_M2_scalar, nT, 1);
XAl_M1_opx = repmat(XAl_M1_scalar, nT, 1);
M1_total_opx = repmat(site_opx.M1_total, nT, 1);
M2_total_opx = repmat(site_opx.M2_total, nT, 1);
isOpxSiteAllocationValid = ...
    repmat(site_opx.isValid, nT, 1);

garnetDivalentSum = ...
    repmat(garnetDivalentSum_scalar, nT, 1);
XMg_garnet = repmat(XMg_garnet_scalar, nT, 1);
AlM1_percent_opx = repmat(AlM1_percent_scalar, nT, 1);
dV_cm3mol = repmat(dV_cm3mol_scalar, nT, 1);
Kterm = repmat(Kterm_scalar, nT, 1);
lnKterm = repmat(lnKterm_scalar, nT, 1);

% Pressure calculation. NaN inputs remain NaN and propagate to P.
G_cal = R_cal .* T_K .* lnKterm + 4207 - 2.69 .* T_K;
P_bar = 1 + (G_cal .* calToCm3Bar ./ dV_cm3mol);
P_kbar = P_bar ./ 1000;

% Applicability flags are diagnostic only and do not stop calculation.
isWithinExperimentalTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 1200 & T_degreeC <= 1400;

isWithinExperimentalPRange = ...
    isfinite(P_kbar) & P_kbar >= 25.9 & P_kbar <= 40.5;

isWithinTable1AlRange = ...
    isfinite(AlM1_percent_opx) & ...
    AlM1_percent_opx >= 2 & AlM1_percent_opx <= 30;

isWithinTable2CompositionRange = ...
    isfinite(XMg_garnet) & XMg_garnet >= 0.67 & XMg_garnet <= 0.85 & ...
    isfinite(XAl_M1_opx) & XAl_M1_opx >= 0.046 & XAl_M1_opx <= 0.115 & ...
    isfinite(XMg_M1_opx) & XMg_M1_opx >= 0.84 & XMg_M1_opx <= 0.93 & ...
    isfinite(XMg_M2_opx) & XMg_M2_opx >= 0.76 & XMg_M2_opx <= 0.95;

% Pack outputs using pre-sized vectors of equal height.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.R_cal_molK = R_output;
row.calToCm3Bar = calToCm3Bar_output;

row.Si_opx = Si_opx;
row.Al_opx = Al_opx;
row.FeT_opx = FeT_opx;
row.Fe2_opx = Fe2_opx;
row.Fe3_opx = Fe3_opx;
row.Mg_opx = Mg_opx;
row.Ca_opx = Ca_opx;
row.Na_opx = Na_opx;
row.Mn_opx = Mn_opx;
row.Ti_opx = Ti_opx;
row.Cr_opx = Cr_opx;

row.Si_grt = Si_grt;
row.Al_grt = Al_grt;
row.FeT_grt = FeT_grt;
row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.Mg_grt = Mg_grt;
row.Ca_grt = Ca_grt;
row.Na_grt = Na_grt;
row.Mn_grt = Mn_grt;
row.Ti_grt = Ti_grt;
row.Cr_grt = Cr_grt;

row.AlIV_opx = AlIV_opx;
row.AlVI_opx = AlVI_opx;
row.XAl_T_opx = XAl_T_opx;
row.XMg_M1_opx = XMg_M1_opx;
row.XMg_M2_opx = XMg_M2_opx;
row.XAl_M1_opx = XAl_M1_opx;
row.M1_total_opx = M1_total_opx;
row.M2_total_opx = M2_total_opx;
row.isOpxSiteAllocationValid = isOpxSiteAllocationValid;

row.garnetDivalentSum = garnetDivalentSum;
row.XMg1_gar = XMg_garnet;
row.XMg_garnet = XMg_garnet;
row.AlM1_percent_opx = AlM1_percent_opx;
row.dV_cm3mol = dV_cm3mol;
row.Kterm = Kterm;
row.lnKterm = lnKterm;
row.G_cal = G_cal;
row.P_bar = P_bar;
row.P_kbar = P_kbar;

row.isWithinExperimentalTRange = isWithinExperimentalTRange;
row.isWithinExperimentalPRange = isWithinExperimentalPRange;
row.isWithinTable1AlRange = isWithinTable1AlRange;
row.isWithinTable2CompositionRange = ...
    isWithinTable2CompositionRange;

% Systematic offset discussed on p. 116; values are metadata only and are
% not subtracted automatically from P_kbar.
row.reportedPossibleOverestimateMin_kbar = repmat(3, nT, 1);
row.reportedPossibleOverestimateMax_kbar = repmat(5, nT, 1);

end

function mineral = prepareMineralRow(data_mineral, mineralLabel)
% prepareMineralRow
% Extract one-row cation data. Required variables must exist; missing
% optional variables are retained as NaN, never replaced by zero.

if height(data_mineral) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();
mineral.Si = getVarOrError( ...
    data_mineral, 'Si_cation_apfu', mineralLabel);
mineral.Al = getVarOrError( ...
    data_mineral, 'Al_cation_apfu', mineralLabel);
% Fe_cation_apfu is the Fe2+ value used directly by the original
% WoodBanno1973baro implementation and by the thermobaroMin cation tables.
% Fe3_cation_apfu is retained independently for Opx M1-site allocation and
% is not subtracted from Fe_cation_apfu.
mineral.Fe2 = getVarOrError( ...
    data_mineral, 'Fe_cation_apfu', mineralLabel);
mineral.Mg = getVarOrError( ...
    data_mineral, 'Mg_cation_apfu', mineralLabel);
mineral.Ca = getVarOrError( ...
    data_mineral, 'Ca_cation_apfu', mineralLabel);

mineral.Na = getVarOrNaN(data_mineral, 'Na_cation_apfu');
mineral.Mn = getVarOrNaN(data_mineral, 'Mn_cation_apfu');
mineral.Ti = getVarOrNaN(data_mineral, 'Ti_cation_apfu');
mineral.Cr = getVarOrNaN(data_mineral, 'Cr_cation_apfu');
mineral.Fe3 = getVarOrNaN(data_mineral, 'Fe3_cation_apfu');

% Total Fe is a diagnostic output only. If Fe3 is NaN, FeT remains NaN;
% this does not change the Fe2+ value used in the barometer calculation.
mineral.FeT = mineral.Fe2 + mineral.Fe3;

fieldNames = fieldnames(mineral);
for i = 1:numel(fieldNames)
    value = mineral.(fieldNames{i});
    if ~isnumeric(value) || ~isscalar(value)
        error('%s variable %s must be a numeric scalar in a 1-row table.', ...
            mineralLabel, fieldNames{i});
    end
    if isinf(value) || (isfinite(value) && value < 0)
        error('%s contains an invalid negative or Inf value for %s.', ...
            mineralLabel, fieldNames{i});
    end
end

end

function site = calcOpxSiteFractions(opx)
% calcOpxSiteFractions
% Calculate Opx site fractions without replacing NaN by zero. Invalid
% derived site capacities are represented by NaN and flagged.

site = struct();
site.AlIV = NaN;
site.AlVI = NaN;
site.XAl_T = NaN;
site.M1_fixed = NaN;
site.M2_fixed = NaN;
site.M1_remaining = NaN;
site.M2_remaining = NaN;
site.Mg_M1 = NaN;
site.Fe2_M1 = NaN;
site.Mg_M2 = NaN;
site.Fe2_M2 = NaN;
site.M1_total = NaN;
site.M2_total = NaN;
site.XMg_M1 = NaN;
site.XMg_M2 = NaN;
site.XAl_M1 = NaN;
site.isValid = true;

% T-site Al from Si deficiency relative to two tetrahedral cations.
if isfinite(opx.Si)
    site.AlIV = max(0, 2 - opx.Si);
end

% M1 Al. Do not clamp a negative derived value to zero.
if isfinite(opx.Al) && isfinite(site.AlIV)
    calculatedAlVI = opx.Al - site.AlIV;
    if calculatedAlVI >= 0
        site.AlVI = calculatedAlVI;
    else
        site.isValid = false;
    end
end

% For diagnostics, tetrahedral Al is expressed on the coupled one-site
% scale used by the MgAl2SiO6 substitution.
site.XAl_T = site.AlIV;

if all(isfinite([site.AlVI, opx.Cr, opx.Ti, opx.Fe3]))
    site.M1_fixed = site.AlVI + opx.Cr + opx.Ti + opx.Fe3;
end

if all(isfinite([opx.Ca, opx.Na, opx.Mn]))
    site.M2_fixed = opx.Ca + opx.Na + opx.Mn;
end

if isfinite(site.M1_fixed)
    if site.M1_fixed >= 0 && site.M1_fixed <= 1 + 1e-8
        site.M1_remaining = max(0, 1 - site.M1_fixed);
    else
        site.isValid = false;
    end
end

if isfinite(site.M2_fixed)
    if site.M2_fixed >= 0 && site.M2_fixed <= 1 + 1e-8
        site.M2_remaining = max(0, 1 - site.M2_fixed);
    else
        site.isValid = false;
    end
end

MgFe_total = opx.Mg + opx.Fe2;
if isfinite(MgFe_total) && MgFe_total > 0 && ...
        isfinite(site.M1_remaining) && isfinite(site.M2_remaining)

    Mg_fraction = opx.Mg ./ MgFe_total;
    Fe_fraction = opx.Fe2 ./ MgFe_total;

    site.Mg_M1 = site.M1_remaining .* Mg_fraction;
    site.Fe2_M1 = site.M1_remaining .* Fe_fraction;
    site.Mg_M2 = site.M2_remaining .* Mg_fraction;
    site.Fe2_M2 = site.M2_remaining .* Fe_fraction;

    site.M1_total = ...
        site.M1_fixed + site.Mg_M1 + site.Fe2_M1;
    site.M2_total = ...
        site.M2_fixed + site.Mg_M2 + site.Fe2_M2;

    if isfinite(site.M1_total) && site.M1_total > 0 && ...
            isfinite(site.M2_total) && site.M2_total > 0
        site.XMg_M1 = site.Mg_M1 ./ site.M1_total;
        site.XMg_M2 = site.Mg_M2 ./ site.M2_total;
        site.XAl_M1 = site.AlVI ./ site.M1_total;
    else
        site.isValid = false;
    end
else
    if isfinite(MgFe_total) && MgFe_total <= 0
        site.isValid = false;
    end
end

requiredDerivedValues = [ ...
    site.XMg_M1, site.XMg_M2, site.XAl_M1, ...
    site.M1_total, site.M2_total];

if any(~isfinite(requiredDerivedValues))
    site.isValid = false;
end

end

function value = getVarOrError(tbl, variableName, mineralLabel)
% getVarOrError
% Retrieve a required scalar table variable without altering NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);
if iscell(value)
    value = value{1};
end

if ~isnumeric(value) || ~isscalar(value)
    error('%s variable %s must be a numeric scalar in a 1-row table.', ...
        mineralLabel, variableName);
end

end

function value = getVarOrNaN(tbl, variableName)
% getVarOrNaN
% Retrieve an optional scalar variable. Missing optional variables are
% represented by NaN, never by zero.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);
    if iscell(value)
        value = value{1};
    end

    if ~isnumeric(value) || ~isscalar(value)
        error('Variable %s must be a numeric scalar in a 1-row table.', ...
            variableName);
    end
else
    value = NaN;
end

end
