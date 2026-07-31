function results = Zhang2020baro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Quartz/Zhang2020baro.m
% Tested with MATLAB R2024b
%
% Quartz-Liquid Ti barometer based on Zhang et al. (2020), Eq. (8)
% Zhang, C., Li, X., Almeev, R.R., Horn, I., Behrens, H. and Holtz, F. (2020)
% Earth and Planetary Science Letters, 538, 116213
% DOI: https://doi.org/10.1016/j.epsl.2020.116213
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Quartz analysis and combines it
% with one Liquid analysis to calculate pressure from the Ti partitioning
% relation of Zhang et al. (2020), Eq. (8).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Quartz-Liquid pair, one output row
% is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Zhang et al. (2020) calibrated Ti partitioning between quartz and a
% coexisting high-silica aluminosilicate melt at:
%
%   Temperature : 700-900 degreeC
%   Pressure    : 0.5-4 kbar
%   System      : hydrous, high-silica, quartz-bearing rhyolitic melt
%
% These are the direct experimental conditions stated in the Abstract and
% listed in Table 1 (pp. 1 and 3). Calculated temperatures or pressures
% outside these limits are extrapolations from the directly calibrated
% Quartz-Liquid dataset. This implementation therefore issues non-stopping
% fprintf warnings outside 700-900 degreeC or 0.5-4 kbar.
%
% The starting compositions were strongly peraluminous synthetic Qz-Ab-Or
% glasses, approximately Qz70-Ab20-Or10, with an aluminum saturation index
% of approximately 1.6-1.8 (p. 2). The experiments were H2O-saturated and
% were designed for silicic magmas stored at shallow crustal depths
% (pp. 1-4). Application to mafic, intermediate, peralkaline, markedly
% metaluminous, hydrothermal, or metamorphic quartz-fluid systems is not a
% direct use of this calibration.
%
% The TiO2-solubility model incorporated into Eq. (8) was proposed for
% silicic melts with SiO2 > 70 wt%, after normalization to 100 wt% on an
% anhydrous basis, over 700-1000 degreeC and 0.5-10 kbar (p. 9). However,
% the wider 700-1000 degreeC and 0.5-10 kbar limits include literature data
% used for the melt-solubility model and are not the direct 0.5-4 kbar,
% 700-900 degreeC calibration range of the complete Quartz-Liquid Eq. (8).
%
% Equation (8) requires quartz and liquid that coexisted and approached
% equilibrium. In natural applications, Zhang et al. (2020) used Ti in
% quartz together with Ti in melt inclusions hosted by quartz (pp. 9-11).
% Pairing quartz with an unrelated matrix glass, a different crystallization
% stage, or a composition modified after entrapment can produce misleading
% pressures.
%
% Equation (8) assumes ideal behavior of TiO2 in the liquid, i.e.
% gamma_TiO2_liq = 1, and requires either pressure or temperature to be
% independently known (p. 9). Here temperature is supplied externally and
% pressure is calculated.
%
% Pressure is highly sensitive to the measured Ti ratio. The natural
% examples have C_Ti_Qtz/C_Ti_Liq values of approximately 0.09-0.13; at
% approximately 800 degreeC this small ratio interval corresponds roughly
% to 5-1 kbar (pp. 1 and 10). High-precision Ti analyses are therefore
% required. Zhang et al. report approximate quartz-Ti detection limits of
% 38 ppm by EPMA and 10 ppm by fs-LA-ICP-MS for their analytical conditions
% (p. 4). Rutile or glass inclusions in the analyzed quartz volume must be
% avoided; experiments containing rutile inclusions in quartz were excluded
% from the calibration (p. 2).
%
% With an assumed temperature uncertainty of +/-25 degreeC and sufficiently
% precise Ti analyses, the Abstract reports pressure precision of about
% +/-0.2 kbar (p. 1). This value must not be treated as a universal error for
% disequilibrium pairs, low-Ti measurements, compositionally inappropriate
% melts, or extrapolated P-T conditions.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 700-900 degreeC,
%   2) finite calculated pressure is outside 0.5-4 kbar,
%   3) finite anhydrous-normalized liquid SiO2 is <= 70 wt%,
%   4) a required or present calculation input contains NaN,
%   5) calculated pressure is NaN or Inf, or
%   6) the rearranged pressure base is finite and negative.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Quartz : table
%
% The FIRST column of the Quartz table is treated as an identifier
% ("data code") displayed in the selection dialog.
%
% Required Quartz variable:
%   Ti_cation_apfu
%
% Liquid data are loaded through liquid.readLiquidExcel(). The following
% liquid oxide variables are required:
%   SiO2, TiO2, Al2O3, Na2O, K2O
%
% The following liquid oxide variables are used when present:
%   FeOt, FeO, Fe2O3, MnO, MgO, CaO, V2O3, Cr2O3, NiO, P2O5, SO3
%
% If an optional oxide column is absent, its contribution is treated as
% absent from the analytical dataset. If the column exists but its selected
% value is NaN, the NaN is retained and propagated through cationTotal_liq,
% FM, and pressure. NaN is never converted to zero.
%
% F and Cl are intentionally excluded from cationTotal_liq and are excluded
% from the NaN-input warning check, as requested for the liquid treatment.
%
% Finite calculation inputs must be non-negative. NaN is allowed, retained
% as missing, propagated through the calculation, and reported by
% non-stopping fprintf warnings. Inf and finite negative values are rejected.
% Zero is retained; if it makes a logarithm, ratio, or pressure inversion
% undefined, the resulting NaN or Inf is retained and reported.
%
% The optional name-value argument 'LiquidRow' selects a row from the liquid
% table. If omitted, row 1 is used. This behavior is retained from the
% original implementation.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
% Zhang et al. (2020), Eq. (8), p. 9:
%
%   log10(C_Ti_Qtz / C_Ti_Liq) =
%     -1.1963 + (1058.1 - 520.4*P^0.2)/T - 0.1155*FM
%
% Rearranged for pressure:
%
%   P(kbar) =
%     ((1058.1 - T(K)*[log10(C_Ti_Qtz/C_Ti_Liq)
%       + 1.1963 + 0.1155*FM]) / 520.4)^5
%
% where:
%
%   FM = (XNa + XK + 2*XCa + 2*XMg + 2*XFe) / (XSi*XAl)
%
% and X denotes the cation molar fraction in the liquid (Eqs. 5 and 8,
% p. 9). C_Ti_Qtz and C_Ti_Liq are elemental Ti concentrations in ppm.
%
% Quartz Ti is converted from Ti_cation_apfu to elemental Ti ppm using the
% low-concentration SiO2 approximation:
%
%   C_Ti_Qtz(ppm) = Ti_apfu * MW_Ti / MW_SiO2 * 1e6
%
% Liquid TiO2 wt% is converted to elemental Ti ppm as:
%
%   C_Ti_Liq(ppm) = TiO2(wt%) * MW_Ti / MW_TiO2 * 1e4
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Zhang2020baro(rawdata_struct, T_degreeC)
%   results = Zhang2020baro(rawdata_struct, T_degreeC, ...
%       'LiquidRow', liquidRowNumber)
%
% Inputs:
%   rawdata_struct : struct containing a Quartz table
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Quartz-Liquid pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('Zhang2020baro requires (rawdata_struct, T_degreeC).');
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

%% 1) Retrieve Quartz dataset
% Extract the required table from the input struct. The source table is not
% modified; selected rows are read during each calculation.
disp('=== Step 1: Preparing Quartz dataset ===');

if ~isfield(rawdata_struct, 'Quartz') || ~istable(rawdata_struct.Quartz)
    error('rawdata_struct must contain table: rawdata_struct.Quartz');
end

dataset_qtz = rawdata_struct.Quartz;

if ~ismember('Ti_cation_apfu', dataset_qtz.Properties.VariableNames)
    error('Quartz table must contain variable: Ti_cation_apfu');
end

disp('=== Preparing Quartz dataset has been finished ===');

%% 2) Parse options and load Liquid dataset
% Keep the original LiquidRow option while validating it as a positive
% integer row number.
disp('=== Step 2: Preparing Liquid dataset ===');

ip = inputParser;
ip.addParameter('LiquidRow', [], @(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && fix(x) == x));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();

if isempty(liqAll)
    error('Selected liquid dataset is empty.');
end

if isempty(liquidRowOpt)
    idxLiq = 1;
    if height(liqAll) > 1
        fprintf(2, ...
            ['WARNING: Liquid dataset has multiple rows (%d). Row 1 is used because ' ...
             'the LiquidRow option was not specified. Confirm that this liquid is ' ...
             'the equilibrium partner of the selected quartz.\n'], ...
            height(liqAll));
    end
else
    idxLiq = liquidRowOpt;
    if idxLiq > height(liqAll)
        error(['Requested LiquidRow (%d) exceeds the number of rows in the ' ...
               'selected liquid dataset (%d).'], idxLiq, height(liqAll));
    end
end

selectedData_liq = liqAll(idxLiq, :);
preparedLiq = prepareLiquidRow(selectedData_liq, MWinfo);

disp(['Liquid selected: Row ' num2str(idxLiq)]);
disp('=== Preparing Liquid dataset has been finished ===');

%% 3) Initialize output container and applicability limits
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 3: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degreeC = 700;
calibrationT_max_degreeC = 900;
calibrationP_min_kbar = 0.5;
calibrationP_max_kbar = 4;

isTemperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);
temperatureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 4-5) Interactive Quartz selection and pressure calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 4: Selecting a data code from the list (Quartz) ===');

while true
    % ----- Quartz selection -----
    % The first table column is used only as the displayed data identifier.
    dataCodes_qtz = dataset_qtz{:, 1};

    [selectedIdx_qtz, ok] = listdlg( ...
        'PromptString', 'Please select the Quartz data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_qtz)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_qtz)
        disp('Selection canceled');
        break;
    end

    selectedCode_qtz = dataCodes_qtz(selectedIdx_qtz);
    disp(['Quartz selected: ' char(string(selectedCode_qtz))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the pressure ===');

    selectedData_qtz = dataset_qtz(selectedIdx_qtz, :);
    preparedQtz = prepareQuartzRow(selectedData_qtz);

    % List NaN only for values actually used by the pressure calculation.
    % F and Cl are neither used nor checked.
    nanInputNames = findNaNInputs(preparedQtz, preparedLiq, T_degreeC);

    row = calcPressure(preparedQtz, preparedLiq, T_degreeC, MWinfo);

    % Repeat identifiers for all temperatures in the current calculation.
    row.dataCode_qtz = repmat(string(selectedCode_qtz), height(row), 1);
    row.dataRow_liq = repmat(idxLiq, height(row), 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_qtz', 'dataRow_liq'}, 'Before', 1);

    % Store one result block per selected Quartz-Liquid pair. Expand the
    % preallocated cell buffer only when its capacity is exhausted, rather
    % than changing table size on every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_qtz)) ' & Liquid row ' ...
            num2str(idxLiq) ': P = ' num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_qtz)) ' & Liquid row ' ...
            num2str(idxLiq) ': P = ' num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Input temperature is common to all selected Quartz-Liquid pairs, so
    % this warning is printed only once after the first completed result.
    if any(isTemperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the direct Zhang et al. (2020) ' ...
             'Quartz-Liquid calibration range of 700-900 degreeC (Abstract, p. 1; ' ...
             'Table 1, p. 3). %d of %d finite temperature point(s) are outside ' ...
             'the range; input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(isTemperatureOutsideCalibration), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures are outside the direct
    % experimental calibration range.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the direct Zhang et al. (2020) ' ...
             'Quartz-Liquid calibration range of 0.5-4 kbar (Abstract, p. 1; ' ...
             'Table 1, p. 3). %d of %d finite pressure point(s) are outside ' ...
             'the range; calculated finite range = %.4g-%.4g kbar for %s and ' ...
             'Liquid row %d.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_qtz)), ...
            idxLiq);
    end

    % The melt-solubility model was proposed for >70 wt% SiO2 after
    % anhydrous normalization. Equality at 70 wt% is outside the stated
    % strictly-greater-than domain.
    finiteSiO2anh = isfinite(row.SiO2_anhydrous_wt(1));
    if finiteSiO2anh && row.SiO2_anhydrous_wt(1) <= 70
        fprintf(2, ...
            ['WARNING: Anhydrous-normalized liquid SiO2 = %.4g wt%% is outside the ' ...
             'stated silicic-melt domain SiO2 > 70 wt%% for the incorporated TiO2 ' ...
             'solubility model of Zhang et al. (2020; p. 9), for %s and Liquid ' ...
             'row %d.\n'], ...
            row.SiO2_anhydrous_wt(1), ...
            char(string(selectedCode_qtz)), ...
            idxLiq);
    end

    % List the exact calculation inputs containing NaN. NaN values are not
    % replaced and calculation continues.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Zhang2020baro input(s) for %s and ' ...
             'Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'pressure may remain NaN. F and Cl are excluded from this check.\n'], ...
            char(string(selectedCode_qtz)), ...
            idxLiq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % A negative rearranged base has no valid non-negative real pressure for
    % the P^0.2 term. The output pressure is retained as NaN.
    negativePbase = isfinite(row.Pbase_Eq8) & row.Pbase_Eq8 < 0;
    if any(negativePbase)
        fprintf(2, ...
            ['WARNING: The rearranged Eq. (8) pressure base is negative for %s and ' ...
             'Liquid row %d (%d of %d points). The corresponding pressure values ' ...
             'were retained as NaN because a non-negative real P^0.2 solution ' ...
             'does not exist.\n'], ...
            char(string(selectedCode_qtz)), ...
            idxLiq, ...
            sum(negativePbase), ...
            numel(row.Pbase_Eq8));
    end

    % Retain and report all non-finite calculated pressures.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s and ' ...
             'Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_qtz)), ...
            idxLiq, ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Quartz selection using the same Liquid dataset?', ...
        'Zhang2020baro', ...
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
function qtz = prepareQuartzRow(data_qtz)
% prepareQuartzRow
% Extract one-row Quartz data without replacing NaN by zero. Inf and finite
% negative values are prohibited; zero and NaN are retained.

if height(data_qtz) ~= 1
    error('Quartz input must be a 1-row table.');
end

qtz = struct();
qtz.Ti_cation_apfu = getTableScalar(data_qtz, 'Ti_cation_apfu', 'Quartz');

value = qtz.Ti_cation_apfu;
if isinf(value) || (isfinite(value) && value < 0)
    error(['Quartz.Ti_cation_apfu must be non-negative. NaN is allowed, ' ...
           'but Inf and finite negative values are prohibited.']);
end

end

function liq = prepareLiquidRow(data_liq, MWinfo)
% prepareLiquidRow
% Extract one-row Liquid oxide data. Required columns must exist. Missing
% optional columns contribute zero because they are absent from the supplied
% analytical dataset; an existing NaN value is retained and propagated.
% F and Cl are deliberately not read or used.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liq = struct();

% Required oxides. Existing NaN values are allowed and retained.
[liq.SiO2, liq.name_SiO2] = readLiquidOxide(data_liq, 'SiO2', true);
[liq.TiO2, liq.name_TiO2] = readLiquidOxide(data_liq, 'TiO2', true);
[liq.Al2O3, liq.name_Al2O3] = readLiquidOxide(data_liq, 'Al2O3', true);
[liq.Na2O, liq.name_Na2O] = readLiquidOxide(data_liq, 'Na2O', true);
[liq.K2O, liq.name_K2O] = readLiquidOxide(data_liq, 'K2O', true);

% Optional oxides. Missing columns are represented by zero contribution;
% present NaN values are retained as NaN.
[liq.FeOt, liq.name_FeOt, liq.has_FeOt] = ...
    readLiquidOxide(data_liq, 'FeOt', false);
[liq.FeO, liq.name_FeO, liq.has_FeO] = ...
    readLiquidOxide(data_liq, 'FeO', false);
[liq.Fe2O3, liq.name_Fe2O3, liq.has_Fe2O3] = ...
    readLiquidOxide(data_liq, 'Fe2O3', false);
[liq.MnO, liq.name_MnO, liq.has_MnO] = ...
    readLiquidOxide(data_liq, 'MnO', false);
[liq.MgO, liq.name_MgO, liq.has_MgO] = ...
    readLiquidOxide(data_liq, 'MgO', false);
[liq.CaO, liq.name_CaO, liq.has_CaO] = ...
    readLiquidOxide(data_liq, 'CaO', false);
[liq.V2O3, liq.name_V2O3, liq.has_V2O3] = ...
    readLiquidOxide(data_liq, 'V2O3', false);
[liq.Cr2O3, liq.name_Cr2O3, liq.has_Cr2O3] = ...
    readLiquidOxide(data_liq, 'Cr2O3', false);
[liq.NiO, liq.name_NiO, liq.has_NiO] = ...
    readLiquidOxide(data_liq, 'NiO', false);
[liq.P2O5, liq.name_P2O5, liq.has_P2O5] = ...
    readLiquidOxide(data_liq, 'P2O5', false);
[liq.SO3, liq.name_SO3, liq.has_SO3] = ...
    readLiquidOxide(data_liq, 'SO3', false);

% Choose the Fe representation without double counting. If FeOt is present,
% it is treated as total Fe expressed as FeO and separate FeO/Fe2O3 columns
% are not used. Otherwise, available FeO and Fe2O3 components are combined.
MW = MWinfo.MW;
if liq.has_FeOt
    liq.nFe = liq.FeOt ./ MW.FeO;
    liq.Fe_input_names = string(liq.name_FeOt);
    liq.Fe_oxide_total_wt = liq.FeOt;
elseif liq.has_FeO || liq.has_Fe2O3
    liq.nFe = liq.FeO ./ MW.FeO + ...
        2 .* liq.Fe2O3 ./ MW.Fe2O3;
    feNames = strings(2, 1);
    nFeNames = 0;
    if liq.has_FeO
        nFeNames = nFeNames + 1;
        feNames(nFeNames) = string(liq.name_FeO);
    end
    if liq.has_Fe2O3
        nFeNames = nFeNames + 1;
        feNames(nFeNames) = string(liq.name_Fe2O3);
    end
    liq.Fe_input_names = feNames(1:nFeNames);
    liq.Fe_oxide_total_wt = liq.FeO + liq.Fe2O3;
else
    liq.nFe = 0;
    liq.Fe_input_names = strings(0, 1);
    liq.Fe_oxide_total_wt = 0;
end

% Validate all values that are used in the calculation. Missing optional
% columns have already been assigned zero and therefore do not trigger an
% error. Existing NaN remains allowed.
calculationValues = [ ...
    liq.SiO2; liq.TiO2; liq.Al2O3; liq.Na2O; liq.K2O; ...
    liq.Fe_oxide_total_wt; liq.MnO; liq.MgO; liq.CaO; ...
    liq.V2O3; liq.Cr2O3; liq.NiO; liq.P2O5; liq.SO3];

if any(isinf(calculationValues) | ...
        (isfinite(calculationValues) & calculationValues < 0))
    error(['Liquid calculation inputs must be non-negative. NaN is allowed, ' ...
           'but Inf and finite negative values are prohibited.']);
end

% Record names and values actually used for the NaN diagnostic. F and Cl are
% intentionally absent. Optional columns are listed only when present.
maxNames = 16;
inputNameBuffer = strings(maxNames, 1);
inputValueBuffer = NaN(maxNames, 1);
nInputs = 0;

requiredNames = {liq.name_SiO2, liq.name_TiO2, liq.name_Al2O3, ...
    liq.name_Na2O, liq.name_K2O};
requiredValues = [liq.SiO2, liq.TiO2, liq.Al2O3, liq.Na2O, liq.K2O];
for i = 1:numel(requiredNames)
    nInputs = nInputs + 1;
    inputNameBuffer(nInputs) = "Liquid." + string(requiredNames{i});
    inputValueBuffer(nInputs) = requiredValues(i);
end

for i = 1:numel(liq.Fe_input_names)
    nInputs = nInputs + 1;
    inputNameBuffer(nInputs) = "Liquid." + liq.Fe_input_names(i);
    if liq.has_FeOt
        inputValueBuffer(nInputs) = liq.FeOt;
    elseif liq.Fe_input_names(i) == string(liq.name_FeO)
        inputValueBuffer(nInputs) = liq.FeO;
    else
        inputValueBuffer(nInputs) = liq.Fe2O3;
    end
end

optionalFlags = [liq.has_MnO, liq.has_MgO, liq.has_CaO, liq.has_V2O3, ...
    liq.has_Cr2O3, liq.has_NiO, liq.has_P2O5, liq.has_SO3];
optionalNames = {liq.name_MnO, liq.name_MgO, liq.name_CaO, liq.name_V2O3, ...
    liq.name_Cr2O3, liq.name_NiO, liq.name_P2O5, liq.name_SO3};
optionalValues = [liq.MnO, liq.MgO, liq.CaO, liq.V2O3, ...
    liq.Cr2O3, liq.NiO, liq.P2O5, liq.SO3];

for i = 1:numel(optionalFlags)
    if optionalFlags(i)
        nInputs = nInputs + 1;
        inputNameBuffer(nInputs) = "Liquid." + string(optionalNames{i});
        inputValueBuffer(nInputs) = optionalValues(i);
    end
end

liq.calculationInputNames = inputNameBuffer(1:nInputs);
liq.calculationInputValues = inputValueBuffer(1:nInputs);

end

function nanInputNames = findNaNInputs(qtz, liq, T_degreeC)
% findNaNInputs
% Return names of pressure-equation inputs containing NaN. NaN values do not
% cause an error and are not replaced by zero. F and Cl are not checked.

maxNames = numel(liq.calculationInputNames) + 2;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = "T_degreeC(indices=" + indexText + ")";
end

if isnan(qtz.Ti_cation_apfu)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Quartz.Ti_cation_apfu";
end

for i = 1:numel(liq.calculationInputNames)
    if isnan(liq.calculationInputValues(i))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = liq.calculationInputNames(i);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function row = calcPressure(qtz, liq, T_degreeC, MWinfo)
% calcPressure
% Compute pressure for one Quartz row and one Liquid row at one or more input
% temperatures. NaN values propagate through the calculation.
%
% Inputs:
%   qtz        : prepared one-row Quartz structure
%   liq        : prepared one-row Liquid structure
%   T_degreeC  : scalar or vector temperature in degreeC
%   MWinfo     : molecular-weight structure from liquid.getMolarWeights
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

MW = MWinfo.MW;
MW_Ti = 47.867;

% Quartz Ti conversion. NaN and zero propagate naturally.
Ti_cation_apfu_qtz_scalar = qtz.Ti_cation_apfu;
Ti_Qtz_ppm_scalar = Ti_cation_apfu_qtz_scalar .* ...
    (MW_Ti ./ MW.SiO2) .* 1e6;

% Liquid Ti conversion from TiO2 wt% to elemental Ti ppm.
Ti_Liq_ppm_scalar = liq.TiO2 .* ...
    (MW_Ti ./ MW.TiO2) .* 1e4;

% Cation molar amounts. F and Cl are intentionally excluded.
nSi = liq.SiO2 ./ MW.SiO2;
nTi = liq.TiO2 ./ MW.TiO2;
nAl = 2 .* liq.Al2O3 ./ MW.Al2O3;
nFe = liq.nFe;
nMn = liq.MnO ./ MW.MnO;
nMg = liq.MgO ./ MW.MgO;
nCa = liq.CaO ./ MW.CaO;
nNa = 2 .* liq.Na2O ./ MW.Na2O;
nK = 2 .* liq.K2O ./ MW.K2O;
nV = 2 .* liq.V2O3 ./ MW.V2O3;
nCr = 2 .* liq.Cr2O3 ./ MW.Cr2O3;
nNi = liq.NiO ./ MW.NiO;
nP = 2 .* liq.P2O5 ./ MW.P2O5;
nS = liq.SO3 ./ MW.SO3;

cationTotal_liq_scalar = nSi + nTi + nAl + nFe + nMn + nMg + ...
    nCa + nNa + nK + nV + nCr + nNi + nP + nS;

XSi_liq_scalar = nSi ./ cationTotal_liq_scalar;
XAl_liq_scalar = nAl ./ cationTotal_liq_scalar;
XNa_liq_scalar = nNa ./ cationTotal_liq_scalar;
XK_liq_scalar = nK ./ cationTotal_liq_scalar;
XCa_liq_scalar = nCa ./ cationTotal_liq_scalar;
XMg_liq_scalar = nMg ./ cationTotal_liq_scalar;
XFe_liq_scalar = nFe ./ cationTotal_liq_scalar;

FM_scalar = (XNa_liq_scalar + XK_liq_scalar + ...
    2 .* XCa_liq_scalar + 2 .* XMg_liq_scalar + ...
    2 .* XFe_liq_scalar) ./ ...
    (XSi_liq_scalar .* XAl_liq_scalar);

% Anhydrous oxide normalization for the stated >70 wt% SiO2 domain.
% F and Cl are excluded. Fe is included once according to the selected Fe
% representation used above.
anhydrousOxideTotal_wt_scalar = liq.SiO2 + liq.TiO2 + liq.Al2O3 + ...
    liq.Na2O + liq.K2O + liq.Fe_oxide_total_wt + liq.MnO + liq.MgO + ...
    liq.CaO + liq.V2O3 + liq.Cr2O3 + liq.NiO + liq.P2O5 + liq.SO3;
SiO2_anhydrous_wt_scalar = ...
    100 .* liq.SiO2 ./ anhydrousOxideTotal_wt_scalar;

% Expand composition-dependent scalars to the temperature-vector length.
Ti_cation_apfu_qtz = repmat(Ti_cation_apfu_qtz_scalar, nT, 1);
Ti_Qtz_ppm = repmat(Ti_Qtz_ppm_scalar, nT, 1);
Ti_Liq_ppm = repmat(Ti_Liq_ppm_scalar, nT, 1);

SiO2_liq = repmat(liq.SiO2, nT, 1);
TiO2_liq = repmat(liq.TiO2, nT, 1);
Al2O3_liq = repmat(liq.Al2O3, nT, 1);
FeOt_liq = repmat(liq.FeOt, nT, 1);
FeO_liq = repmat(liq.FeO, nT, 1);
Fe2O3_liq = repmat(liq.Fe2O3, nT, 1);
MnO_liq = repmat(liq.MnO, nT, 1);
MgO_liq = repmat(liq.MgO, nT, 1);
CaO_liq = repmat(liq.CaO, nT, 1);
Na2O_liq = repmat(liq.Na2O, nT, 1);
K2O_liq = repmat(liq.K2O, nT, 1);
V2O3_liq = repmat(liq.V2O3, nT, 1);
Cr2O3_liq = repmat(liq.Cr2O3, nT, 1);
NiO_liq = repmat(liq.NiO, nT, 1);
P2O5_liq = repmat(liq.P2O5, nT, 1);
SO3_liq = repmat(liq.SO3, nT, 1);

cationTotal_liq = repmat(cationTotal_liq_scalar, nT, 1);
XSi_liq = repmat(XSi_liq_scalar, nT, 1);
XAl_liq = repmat(XAl_liq_scalar, nT, 1);
XNa_liq = repmat(XNa_liq_scalar, nT, 1);
XK_liq = repmat(XK_liq_scalar, nT, 1);
XCa_liq = repmat(XCa_liq_scalar, nT, 1);
XMg_liq = repmat(XMg_liq_scalar, nT, 1);
XFe_liq = repmat(XFe_liq_scalar, nT, 1);
FM = repmat(FM_scalar, nT, 1);
anhydrousOxideTotal_wt = repmat(anhydrousOxideTotal_wt_scalar, nT, 1);
SiO2_anhydrous_wt = repmat(SiO2_anhydrous_wt_scalar, nT, 1);

% Zhang et al. (2020), Eq. (8), rearranged for pressure. No finite-value
% guard is used before the arithmetic: NaN inputs remain NaN. Zero Ti values
% are retained and may produce NaN or Inf through the ratio/logarithm.
log_TiQtz_over_TiLiq = log10(Ti_Qtz_ppm ./ Ti_Liq_ppm);
A_Eq8 = log_TiQtz_over_TiLiq + 1.1963 + 0.1155 .* FM;
Pbase_Eq8 = (1058.1 - T_K .* A_Eq8) ./ 520.4;
P_kbar = Pbase_Eq8 .^ 5;

% A finite negative Pbase is outside the real non-negative pressure domain
% implied by P^0.2. Preserve this diagnostic condition as NaN pressure.
negativePbase = isfinite(Pbase_Eq8) & Pbase_Eq8 < 0;
P_kbar(negativePbase) = NaN;
P_GPa = P_kbar ./ 10;

% Applicability flags are diagnostic only and do not alter results.
isWithinDirectCalibrationTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 700 & T_degreeC <= 900;
isWithinDirectCalibrationPRange = ...
    isfinite(P_kbar) & P_kbar >= 0.5 & P_kbar <= 4;
isHighSilicaMelt = ...
    isfinite(SiO2_anhydrous_wt) & SiO2_anhydrous_wt > 70;
isWithinDirectCalibrationRange = ...
    isWithinDirectCalibrationTRange & ...
    isWithinDirectCalibrationPRange & ...
    isHighSilicaMelt;

% Pack outputs using pre-sized vectors of equal height.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.Ti_cation_apfu_qtz = Ti_cation_apfu_qtz;
row.Ti_Qtz_ppm = Ti_Qtz_ppm;
row.Ti_Liq_ppm = Ti_Liq_ppm;

row.SiO2_liq = SiO2_liq;
row.TiO2_liq = TiO2_liq;
row.Al2O3_liq = Al2O3_liq;
row.FeOt_liq = FeOt_liq;
row.FeO_liq = FeO_liq;
row.Fe2O3_liq = Fe2O3_liq;
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

row.cationTotal_liq = cationTotal_liq;
row.XSi_liq = XSi_liq;
row.XAl_liq = XAl_liq;
row.XNa_liq = XNa_liq;
row.XK_liq = XK_liq;
row.XCa_liq = XCa_liq;
row.XMg_liq = XMg_liq;
row.XFe_liq = XFe_liq;
row.FM = FM;

row.anhydrousOxideTotal_wt = anhydrousOxideTotal_wt;
row.SiO2_anhydrous_wt = SiO2_anhydrous_wt;

row.log_TiQtz_over_TiLiq = log_TiQtz_over_TiLiq;
row.A_Eq8 = A_Eq8;
row.Pbase_Eq8 = Pbase_Eq8;

% Standard launcher-compatible output names.
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% Backward-compatible aliases retained from the original implementation.
row.P_Zhang2020_Eq8_kbar = P_kbar;
row.P_Zhang2020_Eq8_GPa = P_GPa;

row.isWithinDirectCalibrationTRange = ...
    isWithinDirectCalibrationTRange;
row.isWithinDirectCalibrationPRange = ...
    isWithinDirectCalibrationPRange;
row.isHighSilicaMelt = isHighSilicaMelt;
row.isWithinDirectCalibrationRange = ...
    isWithinDirectCalibrationRange;

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Repeat available liquid identifiers to match the number of temperature
% rows in the calculated output block.

variableNames = data_liq.Properties.VariableNames;
nRows = height(row);

if any(strcmp(variableNames, 'Index'))
    row.liq_Index = repmat(data_liq.('Index'), nRows, 1);
end
if any(strcmp(variableNames, 'Experiment'))
    row.liq_Experiment = repmat(string(data_liq.('Experiment')), nRows, 1);
end
if any(strcmp(variableNames, 'Citation'))
    row.liq_Citation = repmat(string(data_liq.('Citation')), nRows, 1);
end
if any(strcmp(variableNames, 'Sample'))
    row.liq_Sample = repmat(string(data_liq.('Sample')), nRows, 1);
end

end

function [value, actualName, isPresent] = readLiquidOxide(data_liq, oxide, isRequired)
% readLiquidOxide
% Retrieve one liquid oxide value. Existing NaN is retained. Missing optional
% oxides contribute zero; missing required oxides cause an error.

actualName = findOxideColumn(data_liq.Properties.VariableNames, oxide);
isPresent = ~isempty(actualName);

if ~isPresent
    if isRequired
        error('Selected liquid row must contain variable: %s', oxide);
    end
    value = 0;
    actualName = oxide;
    return
end

value = toScalarDouble(data_liq.(actualName));

if isinf(value) || (isfinite(value) && value < 0)
    error(['Selected liquid row contains an invalid Inf or finite negative ' ...
           'value for %s.'], actualName);
end

end

function name = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match common oxide-column naming variants, including an optional "value"
% suffix, while preserving the original table variable name.

canonicalNames = cell(size(variableNames));

for i = 1:numel(variableNames)
    textValue = lower(variableNames{i});
    textValue = strrep(textValue, ' ', '');
    textValue = strrep(textValue, '_', '');
    textValue = strrep(textValue, '-', '');
    canonicalNames{i} = textValue;
end

canonicalOxide = lower(oxide);
canonicalOxide = strrep(canonicalOxide, ' ', '');
canonicalOxide = strrep(canonicalOxide, '_', '');
canonicalOxide = strrep(canonicalOxide, '-', '');

targets = {[canonicalOxide 'value'], canonicalOxide};
name = '';

for i = 1:numel(targets)
    idx = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(idx)
        name = variableNames{idx};
        return
    end
end

end

function value = toScalarDouble(raw)
% toScalarDouble
% Convert a one-row table value to scalar double. Missing or non-convertible
% content is represented by NaN and is never converted to zero.

value = NaN;

if isempty(raw)
    return
end

if isnumeric(raw)
    value = double(raw(1));
    return
end

if islogical(raw)
    value = double(raw(1));
    return
end

if isstring(raw)
    if ismissing(raw(1))
        return
    end
    value = str2double(raw(1));
    return
end

if ischar(raw)
    value = str2double(string(raw));
    return
end

if iscell(raw)
    if isempty(raw{1})
        return
    end

    cellValue = raw{1};
    if isnumeric(cellValue) || islogical(cellValue)
        value = double(cellValue(1));
        return
    end
    if isstring(cellValue) || ischar(cellValue)
        value = str2double(string(cellValue));
        return
    end
end

end

function value = getTableScalar(tbl, variableName, tableLabel)
% getTableScalar
% Retrieve a required scalar table variable without altering NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', tableLabel, variableName);
end

value = tbl.(variableName);
if ~isscalar(value)
    error('Variable %s must be scalar in a 1-row table.', variableName);
end

value = toScalarDouble(value);

end
