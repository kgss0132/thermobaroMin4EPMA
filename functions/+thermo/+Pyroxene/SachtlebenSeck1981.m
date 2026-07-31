function results = SachtlebenSeck1981(rawdata_struct, P_kbar)
% functions/+thermo/+Pyroxene/SachtlebenSeck1981.m
% Tested with MATLAB R2024b
%
% Empirical Orthopyroxene-Spinel-Olivine thermometer for spinel peridotite
% Sachtleben, Th. and Seck, H.A. (1981)
% Contributions to Mineralogy and Petrology, 78, 157-165
% DOI: https://doi.org/10.1007/BF00373777
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Orthopyroxene analysis, one Spinel
% analysis, and one Olivine analysis and calculates temperature using the
% empirical thermometer of Sachtleben and Seck (1981).
%
% The function accepts either a scalar pressure or a pressure vector. For
% each selected Opx-Sp-Ol triplet, one output row is returned for every
% pressure supplied in P_kbar. The empirical equation has no pressure term,
% so all pressure rows for one triplet have the same calculated temperature.
% This interface is compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP.
%
% The function is designed for repeated calculations. Each completed result
% is stored in a preallocated cell buffer, and all table blocks are combined
% only once after the interactive loop has finished.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Sachtleben and Seck (1981) derived the thermometer from natural
% Orthopyroxene-Spinel pairs in two suites of Westeifel spinel peridotite.
% The empirical equation is presented on p. 163.
%
% The thermometer was not calibrated directly against a single experimental
% P-T grid. Instead, the regression was fitted to natural Opx-Sp
% compositions and temperatures estimated from Ca solubility in Opx using
% the 15-kbar experimental data of Lindsley and Dixon (1976). The reference
% temperatures are discussed on pp. 160-163.
%
% The two Westeifel suites yielded approximately:
%   Ia suite : 945-980 degreeC
%   Ib suite : 1150-1165 degreeC
%
% In the Discussion, the authors state that, provided equilibrium was
% attained, the phase relations and thermometer should be generally
% applicable to peridotitic rocks equilibrated between approximately
% 950 and 1150 degreeC (p. 163). This published 950-1150 degreeC interval
% is used here for non-stopping temperature-range warnings.
%
% No numerical pressure calibration interval is stated for the empirical
% thermometer. The equation contains no pressure term. The use of 15-kbar
% Ca-in-Opx reference data during construction of the temperature scale must
% not be interpreted as a formal pressure calibration range of exactly
% 15 kbar. Consequently, this implementation stores P_kbar and accepts a
% pressure vector, but it cannot determine whether an input pressure is
% inside or outside a published numerical pressure range. A non-stopping
% fprintf caution is printed once per function call.
%
% Important application cautions from the original paper:
%
%   1) The thermometer is intended for equilibrated spinel peridotite
%      assemblages containing coexisting Opx, Spinel, and Olivine.
%
%   2) The selected analyses must represent mutually equilibrated mineral
%      domains. The Westeifel application was supported by homogeneous
%      mineral compositions, nearly parallel non-crossing Opx-Sp tie lines,
%      and the absence of obvious disequilibrium textures (pp. 163-164).
%
%   3) Porphyroclastic rocks require careful textural and chemical screening.
%      Opx porphyroclast cores, rims, neoblasts, and different Spinel
%      generations may record different stages of cooling. Arbitrary
%      combinations can produce meaningless temperatures (p. 164).
%
%   4) Al-Cr exchange between Opx and Spinel may close at a different
%      temperature from Ca-Mg exchange between the pyroxenes. The empirical
%      regression assumes that the Opx-Sp equilibrium was not blocked at a
%      substantially higher temperature than the reference two-pyroxene
%      exchange equilibrium (pp. 163-164).
%
%   5) Opx Al solubility is chemically controlled by the composition of the
%      coexisting Spinel. The Cr correction must therefore be retained.
%      Ignoring Spinel Cr produced temperatures more than 200 degreeC too
%      high for Cr-rich pairs in the studied suites (pp. 162-163).
%
%   6) The derivation assumes that Al in Opx is assigned to the M1 site.
%      The authors note that some Al may occupy M2, which introduces a
%      site-allocation uncertainty (p. 163).
%
%   7) Fe-Mg ordering in Opx also affects the equilibrium constant. For the
%      Westeifel data, the authors estimated an approximately common
%      temperature effect not exceeding about 25 degreeC, but this estimate
%      need not apply to compositions far from those studied (pp. 162-163).
%
%   8) Spinel Al and Cr analyses require care. The authors reported that
%      different Cr-spinel standards could produce approximately 1 wt%
%      differences in Al2O3 and Cr2O3, complicating interlaboratory
%      comparison (p. 163).
%
% This implementation issues non-stopping fprintf messages when:
%   1) no numerical pressure-range test can be made from the original paper,
%   2) a finite calculated temperature is outside 950-1150 degreeC,
%   3) an explicitly stored calculation input is NaN,
%   4) the Opx M1-site allocation or an equation term is invalid,
%   5) the Carswell-style M1-Al screening relation is triggered, or
%   6) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Opx     : table
%   rawdata_struct.Spinel  : table
%   rawdata_struct.Olivine : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialogs. The mineral data must be normalized
% cations on an appropriate formula basis.
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
% Optional Opx variables; absent columns are assigned zero:
%   Mn_cation_apfu
%   Fe3_cation_apfu
%
% Required Spinel variables:
%   Al_cation_apfu
%   Cr_cation_apfu
%   Mg_cation_apfu
%   Fe_cation_apfu
%
% Optional Spinel variable; an absent column is assigned zero:
%   Fe3_cation_apfu
%
% Required Olivine variables:
%   Fe_cation_apfu
%   Mg_cation_apfu
%
% Optional Olivine variables; absent columns are assigned zero:
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Ni_cation_apfu
%   Fe3_cation_apfu
%
% Existing NaN values are retained and propagated; they are never replaced
% by zero. Zero is assigned only when an optional column itself is absent.
% All finite mineral-composition inputs must be greater than or equal to
% zero. Negative finite values and Inf stop the calculation.
%
% This thermometer does not use a liquid composition. Therefore exclusion
% of F and Cl from cationTotal_liq and liquid NaN warnings is not applicable.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (AS IMPLEMENTED)
%
% The implementation retains the KD definition in the supplied source code:
%
%   KD = ((XMg_Ol)^2 * XM1Al_Opx) ...
%        / ((XM1Mg_Opx * XMg_Sp * XAl_Sp)^2)
%
%   YCr_Sp = Cr_Sp / (Cr_Sp + Al_Sp + Fe3_Sp)
%
% Sachtleben and Seck (1981), empirical equation, p. 163:
%
%   T(degreeC) = ...
%       (4.59 + ln(KD) - 1.552 * YCr_Sp) / 0.0025
%
% Implemented compositional terms:
%
%   XM1Al_Opx = 0.5 * (Al_Opx - Cr_Opx - 2*Ti_Opx + Na_Opx)
%
%   M1_fixed = XM1Al_Opx + Cr_Opx + Ti_Opx + Fe3_Opx
%
%   The remaining M1 capacity is divided between Mg and Fe2+ according to
%   their relative abundances to estimate XM1Mg_Opx.
%
%   XAl_Sp = Al_Sp / (Al_Sp + Cr_Sp + Fe3_Sp)
%   XMg_Sp = Mg_Sp / (Mg_Sp + Fe2_Sp)
%
%   XMg_Ol = Mg_Ol / (Mg_Ol + Fe2_Ol + Mn_Ol + Ca_Ol + Ni_Ol)
%
% Temperature is calculated directly in degreeC. Kelvin is returned as:
%   T_K = T_degreeC + 273.15
%
% -------------------------------------------------------------------------
% Syntax:
%   results = SachtlebenSeck1981(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Opx, Spinel, and Olivine tables
%   P_kbar         : finite non-negative numeric scalar or vector; stored in
%                    the output but not used by the empirical equation
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Opx-Sp-Ol triplet. T_degreeC and T_deg are standardized aliases
%             for downstream launcher and plotting compatibility.
%

%% Input validation
if nargin < 2
    error('SachtlebenSeck1981 requires (rawdata_struct, P_kbar).');
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
disp('=== Step 1: Preparing cation datasets ===');

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

disp('=== Preparing cation datasets has been finished ===');

%% 2) Initialize output container and application limits
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

applicationT_min_degreeC = 950;
applicationT_max_degreeC = 1150;

% The paper provides no numerical pressure calibration interval.
pressureCautionIssued = false;

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
        'PromptString', 'Please select the Opx data you would like to use:', ...
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
        'PromptString', 'Please select the Spinel data you would like to use:', ...
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
        'PromptString', 'Please select the Olivine data you would like to use:', ...
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
        row, {'dataCode_opx', 'dataCode_sp', 'dataCode_ol'}, 'Before', 1);

    % Store one table block. The complete output table is not enlarged on
    % every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    printTemperatureSummary( ...
        selectedCode_opx, selectedCode_sp, selectedCode_ol, row.T_deg);

    % The original paper gives no numerical pressure calibration interval.
    if ~pressureCautionIssued
        fprintf(2, ...
            ['CAUTION: Sachtleben and Seck (1981) do not state a numerical ' ...
             'pressure calibration range for this empirical thermometer, ' ...
             'and the equation contains no pressure term. P_kbar is stored ' ...
             'for interface compatibility but is not used. Input range = ' ...
             '%.4g-%.4g kbar. A pressure-range validity test cannot be ' ...
             'performed from the original paper.\n'], ...
            min(P_kbar), ...
            max(P_kbar));
        pressureCautionIssued = true;
    end

    printTemperatureRangeWarning( ...
        row.T_deg, ...
        applicationT_min_degreeC, ...
        applicationT_max_degreeC, ...
        selectedCode_opx, selectedCode_sp, selectedCode_ol);

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & ' ...
             '%s & %s: %s.\n' ...
             '         NaN values were retained and propagated; the ' ...
             'calculated temperature may be NaN.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_sp), ...
            char(selectedCode_ol), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Carswell-style Opx M1-Al screening used in the supplied implementation.
    if isfinite(row.Na_opx(1)) && isfinite(row.Cr_opx(1)) && ...
            isfinite(row.Ti_opx(1)) && isfinite(row.Fe3_opx(1)) && ...
            row.Na_opx(1) >= ...
            (row.Cr_opx(1) + row.Ti_opx(1) + row.Fe3_opx(1))
        fprintf(2, ...
            ['WARNING: Opx Na >= Cr + Ti + Fe3 for %s. The implemented ' ...
             'Carswell-style XM1Al_Opx estimate may be inappropriate. ' ...
             'The calculated result is retained.\n'], ...
            char(selectedCode_opx));
    end

    invalidTerms = findInvalidEquationTerms(row);
    if ~isempty(invalidTerms)
        fprintf(2, ...
            ['WARNING: Invalid Sachtleben-Seck equation or site-allocation ' ...
             'term(s) were found for %s & %s & %s: %s.\n' ...
             '         Required site fractions, KD, and logarithm arguments ' ...
             'must be finite and physically valid. Affected temperatures ' ...
             'were retained as NaN.\n'], ...
            char(selectedCode_opx), ...
            char(selectedCode_sp), ...
            char(selectedCode_ol), ...
            char(strjoin(invalidTerms, ', ')));
    end

    printNonFiniteTemperatureWarning( ...
        row.T_deg, selectedCode_opx, selectedCode_sp, selectedCode_ol);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'SachtlebenSeck1981', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered table blocks only once.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'primaryTemperatureEquation', 'Sachtleben and Seck (1981)', ...
    'applicationTemperature_degreeC', [950, 1150], ...
    'pressureUsedInEquation', false, ...
    'pressureCalibrationRange', 'Not numerically specified in original paper');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_opx, data_sp, data_ol)
% findNaNInputs
% Return names of explicitly stored variables used by the implemented site
% allocation or thermometer that contain NaN. Missing optional columns are
% assigned zero and therefore are not reported as stored NaN values.

opxVariables = { ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Na_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Fe3_cation_apfu'};

spVariables = { ...
    'Al_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Fe3_cation_apfu'};

olVariables = { ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Ni_cation_apfu'};

maxNames = numel(opxVariables) + numel(spVariables) + numel(olVariables);
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

for i = 1:numel(spVariables)
    variableName = spVariables{i};
    if ismember(variableName, data_sp.Properties.VariableNames)
        value = data_sp.(variableName);
        if isnumeric(value) && any(isnan(value(:)))
            nNames = nNames + 1;
            nameBuffer(nNames) = "Spinel." + string(variableName);
        end
    end
end

for i = 1:numel(olVariables)
    variableName = olVariables{i};
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
% Stop when a stored mineral-composition value used or retained by the
% implementation is negative or infinite. Zero is allowed. NaN is allowed
% and propagated.

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

spVariables = { ...
    'Al_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Fe3_cation_apfu'};

olVariables = { ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Ni_cation_apfu', ...
    'Fe3_cation_apfu'};

maxNames = numel(opxVariables) + numel(spVariables) + numel(olVariables);
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

for i = 1:numel(spVariables)
    variableName = spVariables{i};
    if ismember(variableName, data_sp.Properties.VariableNames)
        value = data_sp.(variableName);
        validateScalarVariable(value, 'Spinel', variableName);
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            nameBuffer(nInvalid) = "Spinel." + string(variableName);
        end
    end
end

for i = 1:numel(olVariables)
    variableName = olVariables{i};
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
    error(['SachtlebenSeck1981: mineral-composition inputs must not be ' ...
           'negative or infinite. Invalid value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_opx, data_sp, data_ol, P_kbar)
% calcTemp
% Calculate one Sachtleben-Seck temperature for one selected Opx-Sp-Ol
% triplet and repeat it for each supplied pressure. NaN inputs and invalid
% derived terms are retained as NaN.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

opx = prepareOpxRow(data_opx, 'Opx');
sp = prepareSpinelRow(data_sp, 'Spinel');
ol = prepareOlivineRow(data_ol, 'Olivine');

site_opx = calcOpxSiteFractions(opx);

% Spinel trivalent sublattice fractions.
denom_sp_tri = sp.Al + sp.Cr + sp.Fe3;
if isfinite(denom_sp_tri) && denom_sp_tri > 0
    XAl_Sp = sp.Al ./ denom_sp_tri;
    YCr_sp = sp.Cr ./ denom_sp_tri;
else
    XAl_Sp = NaN;
    YCr_sp = NaN;
end

% Spinel divalent sublattice fraction.
denom_sp_div = sp.Mg + sp.Fe2;
if isfinite(denom_sp_div) && denom_sp_div > 0
    XMg_Sp = sp.Mg ./ denom_sp_div;
else
    XMg_Sp = NaN;
end

% Olivine Mg fraction.
denom_ol_div = ol.Mg + ol.Fe2 + ol.Mn + ol.Ca + ol.Ni;
if isfinite(denom_ol_div) && denom_ol_div > 0
    XMg_Ol = ol.Mg ./ denom_ol_div;
else
    XMg_Ol = NaN;
end

% Equilibrium constant retained from the supplied implementation.
if isPositiveFinite(XMg_Ol) && ...
        isPositiveFinite(site_opx.XM1Al_Opx) && ...
        isPositiveFinite(site_opx.XM1Mg_Opx) && ...
        isPositiveFinite(XMg_Sp) && ...
        isPositiveFinite(XAl_Sp)

    KD = ((XMg_Ol .^ 2) .* site_opx.XM1Al_Opx) ./ ...
        ((site_opx.XM1Mg_Opx .* XMg_Sp .* XAl_Sp) .^ 2);
else
    KD = NaN;
end

lnKD = safeLogPositive(KD);

if isfinite(lnKD) && isfinite(YCr_sp)
    T_scalar_degreeC = ...
        (4.59 + lnKD - 1.552 .* YCr_sp) ./ 0.0025;
else
    T_scalar_degreeC = NaN;
end

if ~isfinite(T_scalar_degreeC)
    T_scalar_degreeC = NaN;
end

T_degreeC = repmat(T_scalar_degreeC, nP, 1);
T_K = T_degreeC + 273.15;

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PressureUsedInEquation = false(nP, 1);
row.PrimaryEquation = repmat("SachtlebenSeck1981", nP, 1);

% Orthopyroxene inputs.
row.Si_opx = repmat(opx.Si, nP, 1);
row.Al_opx = repmat(opx.Al, nP, 1);
row.Fe2_opx = repmat(opx.Fe2, nP, 1);
row.Fe3_opx = repmat(opx.Fe3, nP, 1);
row.Mg_opx = repmat(opx.Mg, nP, 1);
row.Ca_opx = repmat(opx.Ca, nP, 1);
row.Na_opx = repmat(opx.Na, nP, 1);
row.Mn_opx = repmat(opx.Mn, nP, 1);
row.Ti_opx = repmat(opx.Ti, nP, 1);
row.Cr_opx = repmat(opx.Cr, nP, 1);
row.cationSum_opx = repmat(opx.cationSum, nP, 1);

% Spinel inputs.
row.Al_sp = repmat(sp.Al, nP, 1);
row.Cr_sp = repmat(sp.Cr, nP, 1);
row.Fe2_sp = repmat(sp.Fe2, nP, 1);
row.Fe3_sp = repmat(sp.Fe3, nP, 1);
row.Mg_sp = repmat(sp.Mg, nP, 1);
row.cationSum_sp = repmat(sp.cationSum, nP, 1);

% Olivine inputs.
row.Mg_ol = repmat(ol.Mg, nP, 1);
row.Fe2_ol = repmat(ol.Fe2, nP, 1);
row.Fe3_ol = repmat(ol.Fe3, nP, 1);
row.Mn_ol = repmat(ol.Mn, nP, 1);
row.Ca_ol = repmat(ol.Ca, nP, 1);
row.Ni_ol = repmat(ol.Ni, nP, 1);
row.cationSum_ol = repmat(ol.cationSum, nP, 1);

% Site allocations and equation terms.
row.XM1Al_Opx = repmat(site_opx.XM1Al_Opx, nP, 1);
row.XM1Mg_Opx = repmat(site_opx.XM1Mg_Opx, nP, 1);
row.M1_fixed_opx = repmat(site_opx.M1_fixed, nP, 1);
row.M1_remaining_opx = repmat(site_opx.M1_remaining, nP, 1);
row.M1_total_opx = repmat(site_opx.M1_total, nP, 1);

row.denom_sp_tri = repmat(denom_sp_tri, nP, 1);
row.denom_sp_div = repmat(denom_sp_div, nP, 1);
row.denom_ol_div = repmat(denom_ol_div, nP, 1);

row.XAl_Sp = repmat(XAl_Sp, nP, 1);
row.XMg_Sp = repmat(XMg_Sp, nP, 1);
row.YCr_sp = repmat(YCr_sp, nP, 1);
row.XMg_Ol = repmat(XMg_Ol, nP, 1);

row.KD = repmat(KD, nP, 1);
row.lnKD = repmat(lnKD, nP, 1);

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

opx.Si = getVarRequired(data_opx, 'Si_cation_apfu', mineralLabel);
opx.Al = getVarRequired(data_opx, 'Al_cation_apfu', mineralLabel);
opx.Fe2 = getVarRequired(data_opx, 'Fe_cation_apfu', mineralLabel);
opx.Mg = getVarRequired(data_opx, 'Mg_cation_apfu', mineralLabel);
opx.Ca = getVarRequired(data_opx, 'Ca_cation_apfu', mineralLabel);
opx.Na = getVarRequired(data_opx, 'Na_cation_apfu', mineralLabel);
opx.Ti = getVarRequired(data_opx, 'Ti_cation_apfu', mineralLabel);
opx.Cr = getVarRequired(data_opx, 'Cr_cation_apfu', mineralLabel);

opx.Mn = getVarOptional( ...
    data_opx, 'Mn_cation_apfu', 0, mineralLabel);
opx.Fe3 = getVarOptional( ...
    data_opx, 'Fe3_cation_apfu', 0, mineralLabel);

opx.cationSum = ...
    opx.Si + opx.Al + opx.Fe2 + opx.Fe3 + opx.Mg ...
    + opx.Ca + opx.Na + opx.Mn + opx.Ti + opx.Cr;

end

function sp = prepareSpinelRow(data_sp, mineralLabel)
% prepareSpinelRow
% Extract one Spinel composition. Existing NaN values remain NaN. A missing
% Fe3 column is assigned zero.

if height(data_sp) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

sp = struct();

sp.Al = getVarRequired(data_sp, 'Al_cation_apfu', mineralLabel);
sp.Cr = getVarRequired(data_sp, 'Cr_cation_apfu', mineralLabel);
sp.Mg = getVarRequired(data_sp, 'Mg_cation_apfu', mineralLabel);
sp.Fe2 = getVarRequired(data_sp, 'Fe_cation_apfu', mineralLabel);
sp.Fe3 = getVarOptional( ...
    data_sp, 'Fe3_cation_apfu', 0, mineralLabel);

sp.cationSum = sp.Al + sp.Cr + sp.Mg + sp.Fe2 + sp.Fe3;

end

function ol = prepareOlivineRow(data_ol, mineralLabel)
% prepareOlivineRow
% Extract one Olivine composition. Existing NaN values remain NaN. Missing
% optional columns are assigned zero.

if height(data_ol) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

ol = struct();

ol.Fe2 = getVarRequired(data_ol, 'Fe_cation_apfu', mineralLabel);
ol.Mg = getVarRequired(data_ol, 'Mg_cation_apfu', mineralLabel);

ol.Mn = getVarOptional( ...
    data_ol, 'Mn_cation_apfu', 0, mineralLabel);
ol.Ca = getVarOptional( ...
    data_ol, 'Ca_cation_apfu', 0, mineralLabel);
ol.Ni = getVarOptional( ...
    data_ol, 'Ni_cation_apfu', 0, mineralLabel);
ol.Fe3 = getVarOptional( ...
    data_ol, 'Fe3_cation_apfu', 0, mineralLabel);

ol.cationSum = ol.Fe2 + ol.Fe3 + ol.Mg + ol.Mn + ol.Ca + ol.Ni;

end

function site = calcOpxSiteFractions(opx)
% calcOpxSiteFractions
% Calculate the M1-Al and M1-Mg terms retained from the supplied
% implementation. Invalid or non-finite site allocations are represented by
% NaN rather than stopping the complete calculation.

site = struct();

if all(isfinite([opx.Al, opx.Cr, opx.Ti, opx.Na]))
    XM1Al_Opx = ...
        0.5 .* (opx.Al - opx.Cr - 2 .* opx.Ti + opx.Na);
else
    XM1Al_Opx = NaN;
end

if isfinite(XM1Al_Opx) && XM1Al_Opx >= 0 && XM1Al_Opx <= 1 && ...
        all(isfinite([opx.Cr, opx.Ti, opx.Fe3]))

    M1_fixed = XM1Al_Opx + opx.Cr + opx.Ti + opx.Fe3;
else
    M1_fixed = NaN;
end

if isfinite(M1_fixed) && M1_fixed >= 0 && M1_fixed <= 1
    M1_remaining = 1 - M1_fixed;
else
    M1_remaining = NaN;
end

MgFe_total = opx.Mg + opx.Fe2;
if isfinite(MgFe_total) && MgFe_total > 0 && ...
        isfinite(opx.Mg) && isfinite(opx.Fe2) && ...
        isfinite(M1_remaining)

    Mg_fraction = opx.Mg ./ MgFe_total;
    Fe_fraction = opx.Fe2 ./ MgFe_total;

    Mg_M1 = M1_remaining .* Mg_fraction;
    Fe2_M1 = M1_remaining .* Fe_fraction;

    M1_total = M1_fixed + Mg_M1 + Fe2_M1;

    if isfinite(M1_total) && M1_total > 0
        XM1Mg_Opx = Mg_M1 ./ M1_total;
    else
        XM1Mg_Opx = NaN;
        M1_total = NaN;
    end
else
    XM1Mg_Opx = NaN;
    M1_total = NaN;
end

site.XM1Al_Opx = XM1Al_Opx;
site.XM1Mg_Opx = XM1Mg_Opx;
site.M1_fixed = M1_fixed;
site.M1_remaining = M1_remaining;
site.M1_total = M1_total;

end

function invalidTerms = findInvalidEquationTerms(row)
% findInvalidEquationTerms
% Identify invalid site-allocation, fraction, KD, logarithm, and temperature
% terms. Only the first row is needed because the equation is pressure
% independent.

termBuffer = strings(12, 1);
nTerms = 0;

if ~isfinite(row.XM1Al_Opx(1)) || row.XM1Al_Opx(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XM1Al_Opx";
end
if ~isfinite(row.XM1Mg_Opx(1)) || row.XM1Mg_Opx(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XM1Mg_Opx";
end
if ~isfinite(row.M1_fixed_opx(1)) || ...
        row.M1_fixed_opx(1) < 0 || row.M1_fixed_opx(1) > 1
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "M1_fixed_opx";
end
if ~isfinite(row.XAl_Sp(1)) || row.XAl_Sp(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XAl_Sp";
end
if ~isfinite(row.XMg_Sp(1)) || row.XMg_Sp(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XMg_Sp";
end
if ~isfinite(row.YCr_sp(1)) || row.YCr_sp(1) < 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "YCr_sp";
end
if ~isfinite(row.XMg_Ol(1)) || row.XMg_Ol(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "XMg_Ol";
end
if ~isfinite(row.KD(1)) || row.KD(1) <= 0
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "KD";
end
if ~isfinite(row.lnKD(1))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "ln(KD)";
end
if any(~isfinite(row.T_degreeC))
    nTerms = nTerms + 1;
    termBuffer(nTerms) = "T_degreeC";
end

invalidTerms = termBuffer(1:nTerms);

end

function validateRequiredColumns(tbl, mineralLabel)
% validateRequiredColumns
% Verify the required columns for each mineral table.

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
    error('%s table must contain variable: %s', mineralLabel, variableName);
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
% are handled separately.

if ~isnumeric(value) || ~isscalar(value)
    error('%s variable %s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end

function result = isPositiveFinite(value)
% isPositiveFinite
% Return true only for a finite value greater than zero.

result = isfinite(value) && value > 0;

end

function value = safeLogPositive(value)
% safeLogPositive
% Return ln(value) only when value is finite and strictly positive.

if isPositiveFinite(value)
    value = log(value);
else
    value = NaN;
end

end

function printTemperatureSummary( ...
        selectedCode_opx, selectedCode_sp, selectedCode_ol, temperatureValues)
% printTemperatureSummary
% Display one temperature or the first-to-last values for a pressure vector.

if isscalar(temperatureValues)
    disp([char(selectedCode_opx) ' & ' char(selectedCode_sp) ' & ' ...
        char(selectedCode_ol) ': ' num2str(temperatureValues) ' degreeC']);
else
    disp([char(selectedCode_opx) ' & ' char(selectedCode_sp) ' & ' ...
        char(selectedCode_ol) ': ' num2str(temperatureValues(1)) ' to ' ...
        num2str(temperatureValues(end)) ' degreeC']);
end

end

function printTemperatureRangeWarning( ...
        temperatureValues, minimumTemperature, maximumTemperature, ...
        selectedCode_opx, selectedCode_sp, selectedCode_ol)
% printTemperatureRangeWarning
% Warn when a finite temperature lies outside the published general
% application interval of approximately 950-1150 degreeC.

finiteMask = isfinite(temperatureValues);
outsideMask = finiteMask & ...
    (temperatureValues < minimumTemperature | ...
     temperatureValues > maximumTemperature);

if any(outsideMask)
    finiteValues = temperatureValues(finiteMask);
    fprintf(2, ...
        ['WARNING: Calculated Sachtleben-Seck temperature is outside the ' ...
         'approximately %.4g-%.4g degreeC general application interval ' ...
         'stated for equilibrated peridotitic rocks by Sachtleben and Seck ' ...
         '(1981, p. 163). %d of %d finite point(s) are outside the range; ' ...
         'calculated finite range = %.4g-%.4g degreeC for %s & %s & %s.\n'], ...
        minimumTemperature, ...
        maximumTemperature, ...
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
        temperatureValues, selectedCode_opx, selectedCode_sp, selectedCode_ol)
% printNonFiniteTemperatureWarning
% Report NaN and Inf temperatures without modifying or deleting them.

invalidMask = ~isfinite(temperatureValues);
if any(invalidMask)
    fprintf(2, ...
        ['WARNING: Non-finite Sachtleben-Seck temperature values were ' ...
         'calculated for %s & %s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
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
