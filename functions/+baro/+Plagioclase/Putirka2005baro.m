function results = Putirka2005baro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Plagioclase_Liquid/Putirka2005baro.m
% Tested with MATLAB R2024b
%
% Plagioclase-Liquid barometer, Model C
% Putirka, K.D. (2005)
% American Mineralogist, 90, 336-346
% DOI: https://doi.org/10.2138/am.2005.1449
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Plagioclase analysis with one row
% from a user-selected Liquid dataset and calculates pressure using Model C
% of Putirka (2005), Table 2 (p. 342).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Plagioclase-Liquid pair, one
% output row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Putirka (2005) calibrated Model C using plagioclase-saturated partial-
% melting experiments. From the experiments marked as calibration data for
% "This study" in Table 1 (p. 338), the approximate calibration-data
% envelope is:
%
%   Temperature : approximately 850-1350 degreeC (1123-1623 K)
%   Pressure    : approximately 0.001-20 kbar
%
% The complete experimental database examined in the paper spans a broader
% range of approximately 0.001-27 kbar, 850-1430 degreeC, and 42-73 wt% SiO2
% in the liquid (p. 337). These broader limits are NOT the strict Model C
% calibration range. This implementation therefore uses 850-1350 degreeC
% and 0.001-20 kbar only as non-stopping calibration-envelope warnings.
%
% Model C has a standard error of estimate (SEE) of 1.8 kbar for the
% calibration data and 2.2 kbar for independent test data (pp. 343-344;
% Fig. 4c on p. 344). Mean pressures for nominally isobaric experimental
% groups were recovered to approximately +/-1.0 kbar (p. 343; Fig. 4d on
% p. 344).
%
% Individual low-pressure estimates may scatter by several kbar and may be
% negative even for samples expected to have crystallized near 1 atm. For
% Kilauea lava-tube samples, individual calculated pressures ranged from
% -5.4 to 3.1 kbar or from -2.5 to 2.9 kbar, depending on the thermometer
% used to supply T. Averaging multiple nominally isobaric equilibrium
% Plagioclase-Liquid pairs produced a mean much closer to the expected
% surface pressure (pp. 343-345). A single calculated pressure should
% therefore not be interpreted more precisely than the stated model error.
%
% Plagioclase and Liquid must represent an equilibrium pair. Crystal cores,
% rims, inherited crystals, and matrix glass should not be paired unless
% textural and compositional evidence indicates that they coexisted at the
% same stage of crystallization. The barometer is especially sensitive to
% Na because Ab_pl and Na_liq occur in logarithmic and linear terms. Putirka
% (2005) excludes experiments with suspected Na2O loss and emphasizes Na
% loss during experiments and microbeam analysis of hydrous glasses
% (pp. 337-338 and p. 341).
%
% Liquid components must be calculated as anhydrous cation fractions.
% Following Putirka (2005; pp. 341-342), the cation total used here contains
% only:
%
%   SiO2, TiO2, AlO1.5, FeO, MnO, MgO, CaO, NaO0.5, KO0.5, CrO1.5
%
% Oxide wt% values are not renormalized before conversion to cation
% fractions. H2O is excluded. F and Cl are also excluded from
% cationTotal_liq and are excluded from the NaN-input warning. Other
% components not listed above are retained only as optional diagnostic
% outputs when present and are not included in cationTotal_liq.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside approximately 850-1350 degreeC,
%   2) finite calculated pressure is outside approximately 0.001-20 kbar,
%   3) a calculation input contains NaN,
%   4) a calculated pressure is NaN or Inf, or
%   5) a finite negative pressure is calculated.
%
% Finite inputs used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% input values are rejected. Zero is retained; if it makes a logarithm,
% ratio, oxygen normalization, or cation normalization undefined, the
% resulting NaN or Inf is retained and reported.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Plagioclase : table
%
% The FIRST column of the Plagioclase table is treated as an identifier
% ("data code") displayed in the selection dialog. Oxide columns may be
% named either, for example, "SiO2" or "SiO2_value". The following oxide
% inputs are read for the plagioclase calculation and diagnostic 8-oxygen
% normalization:
%
%   SiO2, TiO2, Al2O3, FeO or FeOt, MnO, MgO,
%   CaO, Na2O, K2O, Cr2O3
%
% The Liquid dataset is selected by liquid.readLiquidExcel(). The following
% liquid oxides are used in cationTotal_liq and the pressure calculation:
%
%   SiO2, TiO2, Al2O3, FeO or FeOt, MnO, MgO,
%   CaO, Na2O, K2O, Cr2O3
%
% Missing columns and missing/non-numeric cell contents are represented as
% NaN rather than zero. NaN values are propagated and listed in an fprintf
% warning. F and Cl are not part of this list and are not checked for NaN.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
% Putirka (2005), Table 2, Model C (p. 342):
%
%   P(kbar) = -42.2
%             + 4.94e-2*T(K)
%             + 1.16e-2*T(K)*ln[(Ab_pl*Al_liq*Ca_liq) /
%                                (An_pl*Na_liq*Si_liq)]
%             - 382.3*(Si_liq)^2
%             + 514.2*(Si_liq)^3
%             - 19.6*ln(Ab_pl)
%             - 139.8*Ca_liq
%             + 287.2*Na_liq
%             + 163.9*K_liq
%
% where:
%
%   T is in Kelvin and P is in kbar.
%   An_pl = Ca / (Ca + Na + K) in plagioclase.
%   Ab_pl = Na / (Ca + Na + K) in plagioclase.
%   Or_pl = K  / (Ca + Na + K) in plagioclase.
%   Si_liq, Al_liq, Ca_liq, Na_liq, and K_liq are anhydrous
%   liquid cation fractions.
%   Natural logarithms are used.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2005baro(rawdata_struct, T_degreeC)
%   results = Putirka2005baro(rawdata_struct, T_degreeC, ...
%                             'LiquidRow', liquidRowNumber)
%
% Inputs:
%   rawdata_struct : struct containing a Plagioclase table
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Name-value option:
%   'LiquidRow' : positive integer scalar or [] (default []). If empty, row
%                 1 of the selected Liquid dataset is used. A non-stopping
%                 fprintf warning is printed when the dataset has more than
%                 one row.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Plagioclase-Liquid pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('Putirka2005baro requires (rawdata_struct, T_degreeC).');
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

%% Options
ip = inputParser;
ip.addParameter('LiquidRow', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && ...
                       x >= 1 && fix(x) == x));
ip.parse(varargin{:});
liquidRowOpt = ip.Results.LiquidRow;

%% 1) Retrieve Plagioclase and Liquid datasets
% Extract the required Plagioclase table and load the Liquid dataset. Source
% tables are not modified; selected rows are prepared during calculation.
disp('=== Step 1: Preparing Plagioclase and Liquid datasets ===');

if ~isfield(rawdata_struct, 'Plagioclase') || ...
        ~istable(rawdata_struct.Plagioclase)
    error('rawdata_struct must contain table: rawdata_struct.Plagioclase');
end

dataset_pl = rawdata_struct.Plagioclase;
if isempty(dataset_pl)
    error('rawdata_struct.Plagioclase is empty.');
end

MWinfo = liquid.getMolarWeights();
[liqAll, metaLiq] = liquid.readLiquidExcel();
if isempty(liqAll)
    error('Selected Liquid dataset is empty.');
end

if isempty(liquidRowOpt)
    selectedIdx_liq = 1;
    if height(liqAll) > 1
        fprintf(2, ...
            ['WARNING: The selected Liquid dataset contains %d rows. ' ...
             'Putirka2005baro is using row 1. Use the ''LiquidRow'' option ' ...
             'to select a different row explicitly.\n'], ...
            height(liqAll));
    end
else
    selectedIdx_liq = liquidRowOpt;
    if selectedIdx_liq > height(liqAll)
        error(['Requested LiquidRow (%d) exceeds the number of rows in the ' ...
               'selected Liquid dataset (%d).'], ...
            selectedIdx_liq, height(liqAll));
    end
end

selectedData_liq = liqAll(selectedIdx_liq, :);
disp(['Liquid selected: Row ' num2str(selectedIdx_liq)]);
disp('=== Preparing Plagioclase and Liquid datasets has been finished ===');

%% 2) Initialize output container and warning limits
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Approximate Model C calibration-data envelope derived from Table 1.
calibrationT_min_degreeC = 850;
calibrationT_max_degreeC = 1350;
calibrationP_min_kbar = 0.001;
calibrationP_max_kbar = 20;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);
temperatureWarningIssued = false;
modelCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive Plagioclase selection + calculation
% The loop continues until the user cancels the selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Plagioclase) ===');

while true
    % ----- Plagioclase selection -----
    % The first table column is used only as the displayed data identifier.
    dataCodes_pl = string(dataset_pl{:, 1});

    [selectedIdx_pl, ok] = listdlg( ...
        'PromptString', 'Please select the Plagioclase data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(dataCodes_pl), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_pl)
        disp('Selection canceled');
        break;
    end

    selectedCode_pl = dataCodes_pl(selectedIdx_pl);
    disp(['Plagioclase selected: ' char(string(selectedCode_pl))]);

    % ----- Prepare inputs and calculate -----
    disp('=== Step 4: Calculating the pressure ===');

    selectedData_pl = dataset_pl(selectedIdx_pl, :);

    % Read all calculation inputs without replacing NaN by zero.
    pl = preparePlagioclaseRow(selectedData_pl, MWinfo);
    liq = prepareLiquidRow(selectedData_liq, MWinfo);

    % List NaN inputs used by the calculation. F and Cl are deliberately
    % absent from this check because they are not used in cationTotal_liq.
    nanInputNames = findNaNInputs(pl, liq, T_degreeC);

    % Reject only Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(pl, liq);

    row = calcBaro(pl, liq, T_degreeC);

    % Repeat identifiers for every input temperature in the current block.
    nRows = height(row);
    row.dataCode_pl = repmat(string(selectedCode_pl), nRows, 1);
    row.dataRow_liq = repmat(selectedIdx_liq, nRows, 1);
    row = attachLiquidIDs(row, selectedData_liq);
    row = movevars(row, {'dataCode_pl', 'dataRow_liq'}, 'Before', 1);

    % Store one block per selected pair. The buffer expands geometrically,
    % rather than changing size on every loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_pl)) ' & Liquid row ' ...
            num2str(selectedIdx_liq) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        finitePressureForDisplay = row.P_kbar(isfinite(row.P_kbar));
        if isempty(finitePressureForDisplay)
            disp([char(string(selectedCode_pl)) ' & Liquid row ' ...
                num2str(selectedIdx_liq) ': P = non-finite for all ' ...
                num2str(height(row)) ' temperature point(s)']);
        else
            disp([char(string(selectedCode_pl)) ' & Liquid row ' ...
                num2str(selectedIdx_liq) ': finite P range = ' ...
                num2str(min(finitePressureForDisplay)) ' to ' ...
                num2str(max(finitePressureForDisplay)) ' kbar']);
        end
    end

    % Print the principal model limitations once per function call.
    if ~modelCautionIssued
        fprintf(2, ...
            ['CAUTION: Putirka (2005) Model C has an SEE of 1.8 kbar for ' ...
             'calibration data and 2.2 kbar for independent test data ' ...
             '(pp. 343-344). Individual low-pressure results can scatter ' ...
             'by several kbar and can be negative. Use equilibrium ' ...
             'Plagioclase-Liquid pairs and, where appropriate, average ' ...
             'multiple nominally isobaric pairs (pp. 343-345).\n']);
        modelCautionIssued = true;
    end

    % Input temperature is common to all selected pairs, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the approximate ' ...
             '850-1350 degreeC Model C calibration-data envelope derived ' ...
             'from Table 1 of Putirka (2005; p. 338). %d of %d finite ' ...
             'temperature point(s) are outside the envelope; input finite ' ...
             'range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures fall outside the approximate
    % calibration-data pressure envelope.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the approximate ' ...
             '0.001-20 kbar Model C calibration-data envelope derived ' ...
             'from Table 1 of Putirka (2005; p. 338). %d of %d finite ' ...
             'pressure point(s) are outside the envelope; calculated ' ...
             'finite range = %.4g-%.4g kbar for %s & Liquid row %d.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_pl)), ...
            selectedIdx_liq);
    end

    % List the exact calculation-input names containing NaN. For vector
    % temperature input, the NaN element indices are included.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Putirka (2005) Model C input(s) ' ...
             'for %s & Liquid row %d: %s.\n' ...
             '         NaN values were retained and propagated; they were ' ...
             'not replaced by zero. F and Cl are excluded from this check.\n'], ...
            char(string(selectedCode_pl)), ...
            selectedIdx_liq, ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s ' ...
             '& Liquid row %d (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_pl)), ...
            selectedIdx_liq, ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressure is retained for diagnosis but is outside the
    % physical and calibration domains.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & ' ...
             'Liquid row %d (%d of %d points). The value(s) were retained ' ...
             'for diagnostic purposes. Putirka (2005) documents comparable ' ...
             'scatter at low pressure (pp. 343-345).\n'], ...
            char(string(selectedCode_pl)), ...
            selectedIdx_liq, ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Plagioclase selection (same Liquid row)?', ...
        'Putirka2005baro', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered blocks once. If the user canceled before any
% calculation, return an empty table.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

results.Properties.UserData = struct( ...
    'liquid', metaLiq, ...
    'liquidRow', selectedIdx_liq, ...
    'reference', 'Putirka (2005), Model C', ...
    'doi', '10.2138/am.2005.1449');

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(pl, liq, T_degreeC)
% findNaNInputs
% Return the names of calculation inputs containing NaN. NaN values do not
% cause an error and are not replaced by zero. Liquid F and Cl are not
% included because they are excluded from cationTotal_liq and Model C.

maxNames = 1 + numel(pl.usedFields) + numel(liq.usedFields);
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

for i = 1:numel(pl.usedFields)
    fieldName = pl.usedFields{i};
    if isnan(pl.oxide.(fieldName))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = pl.label.(fieldName);
    end
end

for i = 1:numel(liq.usedFields)
    fieldName = liq.usedFields{i};
    if isnan(liq.oxide.(fieldName))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = liq.label.(fieldName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(pl, liq)
% validateNonNegativeInputs
% Reject Inf and finite negative values in calculation inputs. Zero and NaN
% are intentionally allowed and retained. Liquid F and Cl are not checked
% because they are not used in cationTotal_liq or Model C.

maxNames = numel(pl.usedFields) + numel(liq.usedFields);
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

for i = 1:numel(pl.usedFields)
    fieldName = pl.usedFields{i};
    value = pl.oxide.(fieldName);
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = pl.label.(fieldName);
    end
end

for i = 1:numel(liq.usedFields)
    fieldName = liq.usedFields{i};
    value = liq.oxide.(fieldName);
    if isinf(value) || (isfinite(value) && value < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = liq.label.(fieldName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['Putirka2005baro: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcBaro(pl, liq, T_degreeC)
% calcBaro
% Compute Putirka (2005) Model C pressure for one Plagioclase row and one
% Liquid row at one or more input temperatures. NaN values propagate through
% all normalizations and the pressure equation.
%
% Inputs:
%   pl          : prepared one-row Plagioclase structure
%   liq         : prepared one-row Liquid structure
%   T_degreeC   : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Expand Plagioclase composition-dependent scalars to the temperature-vector
% length. No NaN substitutions are performed.
XSi_pl = repmat(pl.XSi, nT, 1);
XTi_pl = repmat(pl.XTi, nT, 1);
XAl_pl = repmat(pl.XAl, nT, 1);
XFe_pl = repmat(pl.XFe, nT, 1);
XMn_pl = repmat(pl.XMn, nT, 1);
XMg_pl = repmat(pl.XMg, nT, 1);
XCa_pl = repmat(pl.XCa, nT, 1);
XNa_pl = repmat(pl.XNa, nT, 1);
XK_pl = repmat(pl.XK, nT, 1);
XCr_pl = repmat(pl.XCr, nT, 1);
cationSum_pl = repmat(pl.cationSum, nT, 1);
oxygenSum_pl = repmat(pl.oxygenSum, nT, 1);

An_pl = repmat(pl.An, nT, 1);
Ab_pl = repmat(pl.Ab, nT, 1);
Or_pl = repmat(pl.Or, nT, 1);

% Expand raw Liquid oxides and cation fractions.
SiO2_liq = repmat(liq.oxide.SiO2, nT, 1);
TiO2_liq = repmat(liq.oxide.TiO2, nT, 1);
Al2O3_liq = repmat(liq.oxide.Al2O3, nT, 1);
FeO_liq = repmat(liq.oxide.FeO, nT, 1);
MnO_liq = repmat(liq.oxide.MnO, nT, 1);
MgO_liq = repmat(liq.oxide.MgO, nT, 1);
CaO_liq = repmat(liq.oxide.CaO, nT, 1);
Na2O_liq = repmat(liq.oxide.Na2O, nT, 1);
K2O_liq = repmat(liq.oxide.K2O, nT, 1);
Cr2O3_liq = repmat(liq.oxide.Cr2O3, nT, 1);

F_liq = repmat(liq.auxiliary.F, nT, 1);
Cl_liq = repmat(liq.auxiliary.Cl, nT, 1);
H2O_liq = repmat(liq.auxiliary.H2O, nT, 1);
V2O3_liq = repmat(liq.auxiliary.V2O3, nT, 1);
NiO_liq = repmat(liq.auxiliary.NiO, nT, 1);
P2O5_liq = repmat(liq.auxiliary.P2O5, nT, 1);
SO3_liq = repmat(liq.auxiliary.SO3, nT, 1);
Fe2O3_liq = repmat(liq.auxiliary.Fe2O3, nT, 1);

cationTotal_liq = repmat(liq.cationTotal, nT, 1);
XSiO2_liq = repmat(liq.XSiO2, nT, 1);
XTiO2_liq = repmat(liq.XTiO2, nT, 1);
XAlO1_5_liq = repmat(liq.XAlO1_5, nT, 1);
XFeO_liq = repmat(liq.XFeO, nT, 1);
XMnO_liq = repmat(liq.XMnO, nT, 1);
XMgO_liq = repmat(liq.XMgO, nT, 1);
XCaO_liq = repmat(liq.XCaO, nT, 1);
XNaO0_5_liq = repmat(liq.XNaO0_5, nT, 1);
XKO0_5_liq = repmat(liq.XKO0_5, nT, 1);
XCrO1_5_liq = repmat(liq.XCrO1_5, nT, 1);

% Putirka (2005) Model C terms. Direct evaluation is intentional: NaN
% propagates; zero denominators and log(0) may produce NaN or Inf, which are
% retained and reported by the calling function.
exchangeRatio_ModelC = ...
    (Ab_pl .* XAlO1_5_liq .* XCaO_liq) ./ ...
    (An_pl .* XNaO0_5_liq .* XSiO2_liq);

lnK_ModelC = log(exchangeRatio_ModelC);
lnAb_pl = log(Ab_pl);

P_kbar = ...
    -42.2 ...
    + 4.94e-2 .* T_K ...
    + 1.16e-2 .* T_K .* lnK_ModelC ...
    - 382.3 .* (XSiO2_liq .^ 2) ...
    + 514.2 .* (XSiO2_liq .^ 3) ...
    - 19.6 .* lnAb_pl ...
    - 139.8 .* XCaO_liq ...
    + 287.2 .* XNaO0_5_liq ...
    + 163.9 .* XKO0_5_liq;

P_GPa = P_kbar ./ 10;

% Diagnostic applicability flags. These flags do not stop the calculation.
isWithinEquationDomain = ...
    isfinite(T_K) & T_K > 0 & ...
    isfinite(An_pl) & An_pl > 0 & ...
    isfinite(Ab_pl) & Ab_pl > 0 & ...
    isfinite(XSiO2_liq) & XSiO2_liq > 0 & ...
    isfinite(XAlO1_5_liq) & XAlO1_5_liq > 0 & ...
    isfinite(XCaO_liq) & XCaO_liq > 0 & ...
    isfinite(XNaO0_5_liq) & XNaO0_5_liq > 0 & ...
    isfinite(XKO0_5_liq) & XKO0_5_liq >= 0 & ...
    isfinite(cationTotal_liq) & cationTotal_liq > 0;

isWithinCalibrationTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 850 & T_degreeC <= 1350;

isWithinCalibrationPRange = ...
    isfinite(P_kbar) & P_kbar >= 0.001 & P_kbar <= 20;

isWithinBroadLiquidSiO2Range = ...
    isfinite(SiO2_liq) & SiO2_liq >= 42 & SiO2_liq <= 73;

% Pack outputs using pre-sized vectors of equal height.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;

row.XSi_pl = XSi_pl;
row.XTi_pl = XTi_pl;
row.XAl_pl = XAl_pl;
row.XFe_pl = XFe_pl;
row.XMn_pl = XMn_pl;
row.XMg_pl = XMg_pl;
row.XCa_pl = XCa_pl;
row.XNa_pl = XNa_pl;
row.XK_pl = XK_pl;
row.XCr_pl = XCr_pl;
row.cationSum_pl = cationSum_pl;
row.oxygenSum_pl = oxygenSum_pl;

row.An_pl = An_pl;
row.Ab_pl = Ab_pl;
row.Or_pl = Or_pl;

row.SiO2_liq = SiO2_liq;
row.TiO2_liq = TiO2_liq;
row.Al2O3_liq = Al2O3_liq;
row.FeO_liq = FeO_liq;
row.MnO_liq = MnO_liq;
row.MgO_liq = MgO_liq;
row.CaO_liq = CaO_liq;
row.Na2O_liq = Na2O_liq;
row.K2O_liq = K2O_liq;
row.Cr2O3_liq = Cr2O3_liq;

% Optional diagnostic oxides. None are included in cationTotal_liq.
row.F_liq = F_liq;
row.Cl_liq = Cl_liq;
row.H2O_liq_wt = H2O_liq;
row.V2O3_liq = V2O3_liq;
row.NiO_liq = NiO_liq;
row.P2O5_liq = P2O5_liq;
row.SO3_liq = SO3_liq;
row.Fe2O3_liq = Fe2O3_liq;

row.cationTotal_liq = cationTotal_liq;
row.XSiO2_liq = XSiO2_liq;
row.XTiO2_liq = XTiO2_liq;
row.XAlO1_5_liq = XAlO1_5_liq;
row.XFeO_liq = XFeO_liq;
row.XMnO_liq = XMnO_liq;
row.XMgO_liq = XMgO_liq;
row.XCaO_liq = XCaO_liq;
row.XNaO0_5_liq = XNaO0_5_liq;
row.XKO0_5_liq = XKO0_5_liq;
row.XCrO1_5_liq = XCrO1_5_liq;

row.term_An_pl = An_pl;
row.term_Ab_pl = Ab_pl;
row.term_Si_liq = XSiO2_liq;
row.term_Al_liq = XAlO1_5_liq;
row.term_Ca_liq = XCaO_liq;
row.term_Na_liq = XNaO0_5_liq;
row.term_K_liq = XKO0_5_liq;
row.exchangeRatio_ModelC = exchangeRatio_ModelC;
row.lnK_ModelC = lnK_ModelC;
row.lnAb_pl = lnAb_pl;

% Primary output name for launcher compatibility, followed by backward-
% compatible aliases retained from the original implementation.
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;
row.PModelC_kbar = P_kbar;
row.PModelC_GPa = P_GPa;

row.P_uncertainty_calibrationSEE_kbar = repmat(1.8, nT, 1);
row.P_uncertainty_testSEE_kbar = repmat(2.2, nT, 1);

row.isWithinEquationDomain = isWithinEquationDomain;
row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinCalibrationPRange = isWithinCalibrationPRange;
row.isWithinBroadLiquidSiO2Range = isWithinBroadLiquidSiO2Range;

% Backward-compatible / concise diagnostic aliases.
row.isRecommended_T_range = isWithinCalibrationTRange;
row.isApplicable_composition = isWithinEquationDomain;

end

function pl = preparePlagioclaseRow(data_pl, MWinfo)
% preparePlagioclaseRow
% Read one Plagioclase row, preserve NaN, and calculate cations normalized
% to 8 oxygens. Missing columns are represented by NaN rather than zero.

if height(data_pl) ~= 1
    error('Plagioclase input must be a 1-row table.');
end

pl = struct();
pl.usedFields = { ...
    'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};

aliasSets = { ...
    {'SiO2'}, ...
    {'TiO2'}, ...
    {'Al2O3'}, ...
    {'FeO', 'FeOt'}, ...
    {'MnO'}, ...
    {'MgO'}, ...
    {'CaO'}, ...
    {'Na2O'}, ...
    {'K2O'}, ...
    {'Cr2O3'}};

for i = 1:numel(pl.usedFields)
    fieldName = pl.usedFields{i};
    [value, label] = getOxideOrNaN( ...
        data_pl, aliasSets{i}, 'Plagioclase');
    pl.oxide.(fieldName) = value;
    pl.label.(fieldName) = label;
end

% Molar oxide proportions.
molProp.SiO2 = pl.oxide.SiO2 ./ MWinfo.MW.SiO2;
molProp.TiO2 = pl.oxide.TiO2 ./ MWinfo.MW.TiO2;
molProp.Al2O3 = pl.oxide.Al2O3 ./ MWinfo.MW.Al2O3;
molProp.FeO = pl.oxide.FeO ./ MWinfo.MW.FeO;
molProp.MnO = pl.oxide.MnO ./ MWinfo.MW.MnO;
molProp.MgO = pl.oxide.MgO ./ MWinfo.MW.MgO;
molProp.CaO = pl.oxide.CaO ./ MWinfo.MW.CaO;
molProp.Na2O = pl.oxide.Na2O ./ MWinfo.MW.Na2O;
molProp.K2O = pl.oxide.K2O ./ MWinfo.MW.K2O;
molProp.Cr2O3 = pl.oxide.Cr2O3 ./ MWinfo.MW.Cr2O3;

% Eight-oxygen normalization. NaN and zero totals are not intercepted;
% resulting NaN or Inf values propagate naturally.
pl.oxygenSum = ...
    2 .* molProp.SiO2 + ...
    2 .* molProp.TiO2 + ...
    3 .* molProp.Al2O3 + ...
    molProp.FeO + ...
    molProp.MnO + ...
    molProp.MgO + ...
    molProp.CaO + ...
    molProp.Na2O + ...
    molProp.K2O + ...
    3 .* molProp.Cr2O3;

oxygenRenormalizationFactor = 8 ./ pl.oxygenSum;

pl.XSi = molProp.SiO2 .* oxygenRenormalizationFactor;
pl.XTi = molProp.TiO2 .* oxygenRenormalizationFactor;
pl.XAl = 2 .* molProp.Al2O3 .* oxygenRenormalizationFactor;
pl.XFe = molProp.FeO .* oxygenRenormalizationFactor;
pl.XMn = molProp.MnO .* oxygenRenormalizationFactor;
pl.XMg = molProp.MgO .* oxygenRenormalizationFactor;
pl.XCa = molProp.CaO .* oxygenRenormalizationFactor;
pl.XNa = 2 .* molProp.Na2O .* oxygenRenormalizationFactor;
pl.XK = 2 .* molProp.K2O .* oxygenRenormalizationFactor;
pl.XCr = 2 .* molProp.Cr2O3 .* oxygenRenormalizationFactor;

pl.cationSum = ...
    pl.XSi + pl.XTi + pl.XAl + pl.XFe + pl.XMn + ...
    pl.XMg + pl.XCa + pl.XNa + pl.XK + pl.XCr;

feldsparCationSum = pl.XCa + pl.XNa + pl.XK;
pl.An = pl.XCa ./ feldsparCationSum;
pl.Ab = pl.XNa ./ feldsparCationSum;
pl.Or = pl.XK ./ feldsparCationSum;

end

function liq = prepareLiquidRow(data_liq, MWinfo)
% prepareLiquidRow
% Read one Liquid row and calculate the anhydrous cation fractions required
% by Putirka (2005). NaN is retained. F and Cl are read only as optional
% diagnostics and are excluded from cationTotal_liq and NaN warnings.

if height(data_liq) ~= 1
    error('Liquid input must be a 1-row table.');
end

liq = struct();
liq.usedFields = { ...
    'SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', ...
    'MgO', 'CaO', 'Na2O', 'K2O', 'Cr2O3'};

aliasSets = { ...
    {'SiO2'}, ...
    {'TiO2'}, ...
    {'Al2O3'}, ...
    {'FeO', 'FeOt'}, ...
    {'MnO'}, ...
    {'MgO'}, ...
    {'CaO'}, ...
    {'Na2O'}, ...
    {'K2O'}, ...
    {'Cr2O3'}};

for i = 1:numel(liq.usedFields)
    fieldName = liq.usedFields{i};
    [value, label] = getOxideOrNaN( ...
        data_liq, aliasSets{i}, 'Liquid');
    liq.oxide.(fieldName) = value;
    liq.label.(fieldName) = label;
end

% Optional diagnostic components. They are intentionally not part of the
% Model C cation total or NaN-input warning.
liq.auxiliary.F = getAuxiliaryOxideOrNaN(data_liq, {'F'});
liq.auxiliary.Cl = getAuxiliaryOxideOrNaN(data_liq, {'Cl'});
liq.auxiliary.H2O = getAuxiliaryOxideOrNaN(data_liq, {'H2O'});
liq.auxiliary.V2O3 = getAuxiliaryOxideOrNaN(data_liq, {'V2O3'});
liq.auxiliary.NiO = getAuxiliaryOxideOrNaN(data_liq, {'NiO'});
liq.auxiliary.P2O5 = getAuxiliaryOxideOrNaN(data_liq, {'P2O5'});
liq.auxiliary.SO3 = getAuxiliaryOxideOrNaN(data_liq, {'SO3'});
liq.auxiliary.Fe2O3 = getAuxiliaryOxideOrNaN(data_liq, {'Fe2O3'});

% Cation proportions corresponding to SiO2, TiO2, AlO1.5, FeO, MnO,
% MgO, CaO, NaO0.5, KO0.5, and CrO1.5. F, Cl, H2O, and all other
% auxiliary components are excluded.
n.SiO2 = liq.oxide.SiO2 .* MWinfo.Cat.SiO2 ./ MWinfo.MW.SiO2;
n.TiO2 = liq.oxide.TiO2 .* MWinfo.Cat.TiO2 ./ MWinfo.MW.TiO2;
n.Al2O3 = liq.oxide.Al2O3 .* MWinfo.Cat.Al2O3 ./ MWinfo.MW.Al2O3;
n.FeO = liq.oxide.FeO .* MWinfo.Cat.FeO ./ MWinfo.MW.FeO;
n.MnO = liq.oxide.MnO .* MWinfo.Cat.MnO ./ MWinfo.MW.MnO;
n.MgO = liq.oxide.MgO .* MWinfo.Cat.MgO ./ MWinfo.MW.MgO;
n.CaO = liq.oxide.CaO .* MWinfo.Cat.CaO ./ MWinfo.MW.CaO;
n.Na2O = liq.oxide.Na2O .* MWinfo.Cat.Na2O ./ MWinfo.MW.Na2O;
n.K2O = liq.oxide.K2O .* MWinfo.Cat.K2O ./ MWinfo.MW.K2O;
n.Cr2O3 = liq.oxide.Cr2O3 .* MWinfo.Cat.Cr2O3 ./ MWinfo.MW.Cr2O3;

liq.cationTotal = ...
    n.SiO2 + n.TiO2 + n.Al2O3 + n.FeO + n.MnO + ...
    n.MgO + n.CaO + n.Na2O + n.K2O + n.Cr2O3;

% No finite-value guard is used. A NaN input produces a NaN total and NaN
% fractions; a zero total produces NaN/Inf through direct division.
liq.XSiO2 = n.SiO2 ./ liq.cationTotal;
liq.XTiO2 = n.TiO2 ./ liq.cationTotal;
liq.XAlO1_5 = n.Al2O3 ./ liq.cationTotal;
liq.XFeO = n.FeO ./ liq.cationTotal;
liq.XMnO = n.MnO ./ liq.cationTotal;
liq.XMgO = n.MgO ./ liq.cationTotal;
liq.XCaO = n.CaO ./ liq.cationTotal;
liq.XNaO0_5 = n.Na2O ./ liq.cationTotal;
liq.XKO0_5 = n.K2O ./ liq.cationTotal;
liq.XCrO1_5 = n.Cr2O3 ./ liq.cationTotal;

end

function row = attachLiquidIDs(row, data_liq)
% attachLiquidIDs
% Repeat common Liquid identifiers so every output row has matching height.

variableNames = data_liq.Properties.VariableNames;
nRows = height(row);

if ismember('Index', variableNames)
    rawIndex = data_liq.('Index');
    if isnumeric(rawIndex) || islogical(rawIndex)
        row.liq_Index = repmat(rawIndex(1), nRows, 1);
    else
        row.liq_Index = repmat(string(rawIndex(1)), nRows, 1);
    end
end

if ismember('Experiment', variableNames)
    rawExperiment = data_liq.('Experiment');
    row.liq_Experiment = repmat(string(rawExperiment(1)), nRows, 1);
end

if ismember('Citation', variableNames)
    rawCitation = data_liq.('Citation');
    row.liq_Citation = repmat(string(rawCitation(1)), nRows, 1);
end

end

function [value, label] = getOxideOrNaN(data_tbl, aliases, phaseLabel)
% getOxideOrNaN
% Retrieve one scalar oxide value using one or more accepted aliases.
% Missing columns or non-numeric/missing contents return NaN, never zero.

if ischar(aliases) || isstring(aliases)
    aliases = cellstr(string(aliases));
end

variableName = '';
for i = 1:numel(aliases)
    variableName = findOxideColumn( ...
        data_tbl.Properties.VariableNames, aliases{i});
    if ~isempty(variableName)
        break;
    end
end

if isempty(variableName)
    value = NaN;
    aliasText = strjoin(string(aliases), '/');
    label = string(phaseLabel) + "." + aliasText + "(missing)";
    return
end

value = toScalarDouble(data_tbl.(variableName));
label = string(phaseLabel) + "." + string(variableName);

end

function value = getAuxiliaryOxideOrNaN(data_tbl, aliases)
% getAuxiliaryOxideOrNaN
% Retrieve an optional diagnostic oxide. Missing or non-numeric contents are
% represented by NaN. These values are not used in Model C.

[ value, ~ ] = getOxideOrNaN(data_tbl, aliases, 'Liquid');

end

function variableName = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match common oxide-column variants such as SiO2 and SiO2_value.

canonicalNames = cell(size(variableNames));
for i = 1:numel(variableNames)
    canonicalNames{i} = canonicalizeName(variableNames{i});
end

canonicalOxide = canonicalizeName(oxide);
targets = {[canonicalOxide 'value'], canonicalOxide};

variableName = '';
for i = 1:numel(targets)
    index = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(index)
        variableName = variableNames{index};
        return
    end
end

end

function canonicalName = canonicalizeName(name)
% canonicalizeName
% Lower-case a name and remove common separators.

canonicalName = lower(char(string(name)));
canonicalName = strrep(canonicalName, ' ', '');
canonicalName = strrep(canonicalName, '_', '');
canonicalName = strrep(canonicalName, '-', '');

end

function value = toScalarDouble(raw)
% toScalarDouble
% Convert a one-row table entry to a scalar double. Missing or non-numeric
% values return NaN. Explicit NaN and Inf are preserved.

value = NaN;

if isempty(raw)
    return
end

if iscell(raw)
    if isempty(raw{1})
        return
    end
    raw = raw{1};
end

if isnumeric(raw) || islogical(raw)
    if ~isscalar(raw)
        error('Oxide input must be scalar in a 1-row table.');
    end
    value = double(raw);
    return
end

if isstring(raw)
    if numel(raw) ~= 1 || ismissing(raw)
        return
    end
    convertedValue = str2double(raw);
    if ~isnan(convertedValue) || strcmpi(strtrim(raw), 'NaN')
        value = convertedValue;
    end
    return
end

if ischar(raw)
    convertedValue = str2double(string(raw));
    if ~isnan(convertedValue) || strcmpi(strtrim(string(raw)), 'NaN')
        value = convertedValue;
    end
    return
end

if iscategorical(raw)
    if numel(raw) ~= 1 || isundefined(raw)
        return
    end
    convertedValue = str2double(string(raw));
    if ~isnan(convertedValue) || strcmpi(strtrim(string(raw)), 'NaN')
        value = convertedValue;
    end
end

end
