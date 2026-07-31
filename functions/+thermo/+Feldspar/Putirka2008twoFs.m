function results = Putirka2008twoFs(rawdata_struct, P_kbar, varargin)
% functions/+thermo/+Feldsper/Putirka2008twoFs.m
% Tested with MATLAB R2024b
%
% Two-feldspar thermometers
% Putirka, K.D. (2008)
% Reviews in Mineralogy & Geochemistry, 69, 61-120
% DOI: https://doi.org/10.2138/rmg.2008.69.3
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Plagioclase analysis and one alkali
% Feldspar analysis and calculates temperature using Equations (27a) and
% (27b) of Putirka (2008).
%
% The function accepts pressure as either a scalar or vector. Therefore, it
% can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for each pressure value
% for every user-selected Plagioclase-Feldspar pair.
%
% The function is designed for repeated calculations. Each result block is
% stored in a preallocated cell buffer and all blocks are concatenated only
% once after the interactive loop finishes.
%
% Equation (27a) is assigned to the common output variables T_degreeC, T_deg,
% and T_K so that the standard fixed-pressure and pressure-range launchers
% can use the result directly. Equation-specific outputs for both equations
% are retained separately.
%
% The varargin input is retained for interface compatibility with the
% original function. No optional parameter is currently used.
%
% -------------------------------------------------------------------------
% CALIBRATION CONTEXT, PRACTICAL RANGE, AND APPLICATION NOTES
%
% Putirka (2008) developed the two-feldspar thermometers from only 41
% experimental observations compiled from seven studies. The limited number
% of experimental data is explicitly identified as an important limitation
% of two-feldspar thermometer development (pp. 83-84).
%
% Equation (27a):
%   Calibration data : 30 of the 41 experimental observations
%   Calibration error: +/-23 degreeC
%   Test-data error  : +/-44 degreeC
%
% Equation (27b):
%   Calibration data : global calibration using all 41 observations
%   Model error      : +/-30 degreeC
%
% The equations, calibration-data selection, and error statistics are given
% on p. 84. Figure 5h on p. 81 illustrates calibration and test results over
% an approximate experimental temperature coverage of about 550-1050
% degreeC. Putirka (2008) does not state 550-1050 degreeC as a strict formal
% validity interval; it is used here only as a practical Figure-5h screening
% range for non-stopping fprintf warnings.
%
% Putirka (2008) does not report a single numerical pressure calibration
% interval specifically for Equations (27a) and (27b). Pressure therefore
% cannot be screened against a formally published equation-specific range.
% This implementation does not invent pressure limits; it prints one
% non-stopping fprintf message explaining this limitation.
%
% Important application notes:
%   1) The selected Plagioclase and alkali Feldspar must represent an
%      equilibrium pair from the same crystallization or re-equilibration
%      stage. Core-rim mismatches, inherited crystals, exsolution, alteration,
%      and diffusional re-equilibration can invalidate the result.
%   2) Temperature estimates are highly sensitive to small changes in XAb,
%      XAn, and XOr. Putirka (2008) discusses this compositional sensitivity
%      and the problems caused by adjusting experimental compositions on
%      pp. 83-84.
%   3) Mineral components must be calculated in the same manner as in the
%      calibration. General cation-fraction and feldspar-component methods
%      are described on pp. 68-71, and the XSi(afs) definition is stated on
%      p. 84.
%   4) Putirka (2008) recommends comparison among independent feldspar-based
%      thermometers and saturation temperatures. In the Toba application,
%      disequilibrium may cause two-feldspar temperatures to be approximately
%      40-70 degreeC too low (discussion on pp. 101 and 104).
%   5) The approximately 550-1050 degreeC interval used below is only the
%      plotted experimental coverage in Figure 5h, not a formally stated
%      strict calibration range.
%
% This implementation issues non-stopping fprintf messages when:
%   1) pressure applicability cannot be screened because no numerical range
%      is reported specifically for Equations (27a) and (27b),
%   2) a finite temperature from either equation lies outside the practical
%      Figure-5h screening interval of 550-1050 degreeC,
%   3) an existing oxide input used by the calculation contains NaN, or
%   4) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Plagioclase : table
%
% The alkali-feldspar table may use either of these field names:
%   rawdata_struct.Feldspar    : preferred spelling
%   rawdata_struct.Feldsper    : accepted legacy spelling
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. Oxide columns are located by names
% equivalent to the following, optionally with a "_value" suffix:
%
%   SiO2, Al2O3, CaO, Na2O, K2O
%
% SiO2 and Al2O3 columns must exist. CaO, Na2O, and K2O are optional only in
% the sense that a completely absent column is assigned zero. If an existing
% input cell contains NaN, that NaN is retained rather than converted to
% zero. It propagates through component calculations and is reported by a
% non-stopping fprintf message.
%
% All finite oxide values used by the calculation must be greater than or
% equal to zero. A finite value below zero or an Inf value stops the
% calculation with an error. Zero is allowed as an analytical value, but a
% zero that makes a logarithmic term undefined produces NaN and a
% non-stopping result warning.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (Putirka, 2008, Eq. 27a-b, p. 84)
%
% Feldspar components:
%   XAn = XCa / (XCa + XNa + XK)
%   XAb = XNa / (XCa + XNa + XK)
%   XOr = XK  / (XCa + XNa + XK)
%
% Equation (27a):
%   10^4 / T(degreeC) =
%       9.8
%     - 0.098 P(kbar)
%     - 2.46 ln(XAb_afs / XAb_plag)
%     - 14.2 XSi_afs
%     + 423 XCa_afs
%     - 2.42 ln(XAn_afs)
%     - 11.4 (XAn_plag * XAb_plag)
%
% Equation (27b):
%   T(degreeC) =
%     [-442 - 3.72 P(kbar)] /
%     [-0.11
%      + 0.1 ln(XAb_afs / XAb_plag)
%      - 3.27 XAn_afs
%      + 0.098 ln(XAn_afs)
%      + 0.52 (XAn_plag * XAb_plag)]
%
% XSi_afs and XCa_afs are cation fractions calculated from the complete
% Si-Al-Ca-Na-K cation sum. XAn, XAb, and XOr are normalized on the
% Ca-Na-K sum.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Putirka2008twoFs(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Plagioclase and alkali-Feldspar tables
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every selected
%             Plagioclase-Feldspar pair
%

%% Input validation
if nargin < 2
    error('Putirka2008twoFs requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || ~isreal(P_kbar) || isempty(P_kbar) || ...
        ~isvector(P_kbar) || any(~isfinite(P_kbar(:))) || ...
        any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve feldspar datasets
% Extract the two tables and retain their original data. The preferred
% alkali-feldspar field is Feldspar; Feldsper is accepted for compatibility
% with legacy workbook or package naming.
disp('=== Step 1: Preparing feldspar datasets ===');

if ~isfield(rawdata_struct, 'Plagioclase') || ...
        ~istable(rawdata_struct.Plagioclase)
    error('rawdata_struct must contain table: rawdata_struct.Plagioclase');
end

dataset_pl = rawdata_struct.Plagioclase;

if isfield(rawdata_struct, 'Feldspar') && ...
        istable(rawdata_struct.Feldspar)
    dataset_afs = rawdata_struct.Feldspar;
    alkaliFeldsparField = 'Feldspar';
elseif isfield(rawdata_struct, 'Feldsper') && ...
        istable(rawdata_struct.Feldsper)
    dataset_afs = rawdata_struct.Feldsper;
    alkaliFeldsparField = 'Feldsper';
else
    error(['rawdata_struct must contain an alkali-feldspar table named ' ...
           'rawdata_struct.Feldspar or rawdata_struct.Feldsper.']);
end

validateRequiredColumns(dataset_pl, 'Plagioclase');
validateRequiredColumns(dataset_afs, alkaliFeldsparField);

MWinfo = liquid.getMolarWeights();

disp(['Alkali-feldspar field used: rawdata_struct.' alkaliFeldsparField]);
disp('=== Preparing feldspar datasets has been finished ===');

%% 2) Initialize output container and application checks
% Repeated table concatenation inside the interactive loop is avoided. Each
% result is stored in a preallocated cell buffer and concatenated once.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Figure 5h practical screening interval. This is not described by Putirka
% (2008) as a strict formal calibration range.
practicalT_min_degC = 550;
practicalT_max_degC = 1050;
pressureRangeMessageIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop and calculation
disp('=== Step 3: Selecting a data code from the list (Plagioclase) ===');

while true
    % ----- Plagioclase selection -----
    dataCodes_pl = dataset_pl{:, 1};

    [selectedIdx_pl, ok] = listdlg( ...
        'PromptString', ...
        'Please select the Plagioclase data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_pl)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_pl)
        disp('Selection canceled');
        break;
    end

    selectedCode_pl = dataCodes_pl(selectedIdx_pl);
    selectedData_pl = dataset_pl(selectedIdx_pl, :);
    disp(['Plagioclase selected: ' char(string(selectedCode_pl))]);

    % ----- Alkali-feldspar selection -----
    disp('=== Step 4: Selecting a data code from the list (Feldspar) ===');

    dataCodes_afs = dataset_afs{:, 1};

    [selectedIdx_afs, ok] = listdlg( ...
        'PromptString', ...
        'Please select the alkali Feldspar data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_afs)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_afs)
        disp('Selection canceled');
        break;
    end

    selectedCode_afs = dataCodes_afs(selectedIdx_afs);
    selectedData_afs = dataset_afs(selectedIdx_afs, :);
    disp(['Feldspar selected: ' char(string(selectedCode_afs))]);

    % ----- Input checks -----
    % Negative finite values and Inf are rejected. Existing NaN values are
    % permitted and retained so that they propagate into dependent results.
    validateInputValues(selectedData_pl, selectedData_afs);
    nanInputNames = findNaNInputs(selectedData_pl, selectedData_afs);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    row = calcTemp(selectedData_pl, selectedData_afs, P_kbar, MWinfo);

    % Add identifiers to every pressure row for traceability.
    nRows = height(row);
    row.dataCode_pl = repmat(string(selectedCode_pl), nRows, 1);
    row.dataCode_afs = repmat(string(selectedCode_afs), nRows, 1);
    row = movevars(row, {'dataCode_pl', 'dataCode_afs'}, 'Before', 1);

    % Store the complete table block. Capacity is doubled only when needed.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % ----- Immediate result display -----
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_pl)) ' & ' ...
            char(string(selectedCode_afs)) ...
            ': Eq. 27a = ' num2str(row.T_Eq27a_C) ...
            ' degreeC, Eq. 27b = ' num2str(row.T_Eq27b_C) ' degreeC']);
    else
        disp([char(string(selectedCode_pl)) ' & ' ...
            char(string(selectedCode_afs)) ...
            ': Eq. 27a = ' num2str(row.T_Eq27a_C(1)) ' to ' ...
            num2str(row.T_Eq27a_C(end)) ...
            ' degreeC; Eq. 27b = ' num2str(row.T_Eq27b_C(1)) ' to ' ...
            num2str(row.T_Eq27b_C(end)) ' degreeC']);
    end

    % ----- Pressure applicability message -----
    % No equation-specific numerical pressure range is published.
    if ~pressureRangeMessageIssued
        fprintf(2, ...
            ['WARNING: Putirka (2008) does not report a numerical pressure ' ...
             'calibration range specifically for Equations (27a) and (27b) ' ...
             '(pp. 83-84). The input pressure therefore cannot be screened ' ...
             'against a published equation-specific range; input range = ' ...
             '%.4g-%.4g kbar.\n'], ...
            min(P_kbar), max(P_kbar));
        pressureRangeMessageIssued = true;
    end

    % ----- Practical temperature-range messages -----
    finiteTemperature27a = isfinite(row.T_Eq27a_C);
    outsideTemperature27a = finiteTemperature27a & ...
        (row.T_Eq27a_C < practicalT_min_degC | ...
         row.T_Eq27a_C > practicalT_max_degC);

    if any(outsideTemperature27a)
        finiteValues27a = row.T_Eq27a_C(finiteTemperature27a);
        fprintf(2, ...
            ['WARNING: Equation (27a) calculated temperature is outside the ' ...
             'approximate 550-1050 degreeC experimental coverage illustrated ' ...
             'in Putirka (2008), Figure 5h (p. 81). This interval is a ' ...
             'practical screening range, not a strict formal calibration ' ...
             'limit. %d of %d finite point(s) are outside; calculated finite ' ...
             'range = %.4g-%.4g degreeC for %s & %s.\n'], ...
            sum(outsideTemperature27a), ...
            sum(finiteTemperature27a), ...
            min(finiteValues27a), ...
            max(finiteValues27a), ...
            char(string(selectedCode_pl)), ...
            char(string(selectedCode_afs)));
    end

    finiteTemperature27b = isfinite(row.T_Eq27b_C);
    outsideTemperature27b = finiteTemperature27b & ...
        (row.T_Eq27b_C < practicalT_min_degC | ...
         row.T_Eq27b_C > practicalT_max_degC);

    if any(outsideTemperature27b)
        finiteValues27b = row.T_Eq27b_C(finiteTemperature27b);
        fprintf(2, ...
            ['WARNING: Equation (27b) calculated temperature is outside the ' ...
             'approximate 550-1050 degreeC experimental coverage illustrated ' ...
             'in Putirka (2008), Figure 5h (p. 81). This interval is a ' ...
             'practical screening range, not a strict formal calibration ' ...
             'limit. %d of %d finite point(s) are outside; calculated finite ' ...
             'range = %.4g-%.4g degreeC for %s & %s.\n'], ...
            sum(outsideTemperature27b), ...
            sum(finiteTemperature27b), ...
            min(finiteValues27b), ...
            max(finiteValues27b), ...
            char(string(selectedCode_pl)), ...
            char(string(selectedCode_afs)));
    end

    % ----- NaN input message -----
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s & ' ...
             '%s: %s.\n' ...
             '         NaN was retained and the calculation was continued; ' ...
             'dependent output values may also be NaN.\n'], ...
            char(string(selectedCode_pl)), ...
            char(string(selectedCode_afs)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % ----- Non-finite result messages -----
    invalidTemperature27a = ~isfinite(row.T_Eq27a_C);
    if any(invalidTemperature27a)
        fprintf(2, ...
            ['WARNING: Non-finite Equation (27a) temperature values were ' ...
             'calculated for %s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped. Zero or NaN in a logarithmic ' ...
             'component may be responsible.\n'], ...
            char(string(selectedCode_pl)), ...
            char(string(selectedCode_afs)), ...
            sum(invalidTemperature27a), ...
            numel(row.T_Eq27a_C), ...
            sum(isnan(row.T_Eq27a_C)), ...
            sum(isinf(row.T_Eq27a_C)));
    end

    invalidTemperature27b = ~isfinite(row.T_Eq27b_C);
    if any(invalidTemperature27b)
        fprintf(2, ...
            ['WARNING: Non-finite Equation (27b) temperature values were ' ...
             'calculated for %s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped. Zero or NaN in a logarithmic ' ...
             'component may be responsible.\n'], ...
            char(string(selectedCode_pl)), ...
            char(string(selectedCode_afs)), ...
            sum(invalidTemperature27b), ...
            numel(row.T_Eq27b_C), ...
            sum(isnan(row.T_Eq27b_C)), ...
            sum(isinf(row.T_Eq27b_C)));
    end

    disp('--------------------------------------------------');

    % Ask whether another pair should be calculated.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Putirka2008twoFs', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate buffered result blocks only once. Return an empty table when no
% calculation was performed.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredColumns(dataTable, phaseName)
% validateRequiredColumns
% Confirm that the two required oxide columns exist before selection begins.

requiredOxides = {'SiO2', 'Al2O3'};
missingBuffer = strings(numel(requiredOxides), 1);
nMissing = 0;

for i = 1:numel(requiredOxides)
    oxide = requiredOxides{i};
    name = findOxideColumn(dataTable.Properties.VariableNames, oxide);
    if isempty(name)
        nMissing = nMissing + 1;
        missingBuffer(nMissing) = string(oxide);
    end
end

if nMissing > 0
    missingNames = missingBuffer(1:nMissing);
    error('%s table is missing required oxide column(s): %s.', ...
        phaseName, char(strjoin(missingNames, ', ')));
end

end

function nanInputNames = findNaNInputs(data_pl, data_afs)
% findNaNInputs
% Return names of existing oxide inputs containing NaN. Missing optional
% columns are not reported because they are intentionally assigned zero.

oxideNames = {'SiO2', 'Al2O3', 'CaO', 'Na2O', 'K2O'};
maxNames = 2 .* numel(oxideNames);
nanNameBuffer = strings(maxNames, 1);
nNames = 0;

for i = 1:numel(oxideNames)
    oxide = oxideNames{i};
    name = findOxideColumn(data_pl.Properties.VariableNames, oxide);
    if ~isempty(name)
        value = toScalarDouble(data_pl.(name));
        if isnan(value)
            nNames = nNames + 1;
            nanNameBuffer(nNames) = "Plagioclase." + string(name);
        end
    end
end

for i = 1:numel(oxideNames)
    oxide = oxideNames{i};
    name = findOxideColumn(data_afs.Properties.VariableNames, oxide);
    if ~isempty(name)
        value = toScalarDouble(data_afs.(name));
        if isnan(value)
            nNames = nNames + 1;
            nanNameBuffer(nNames) = "Feldspar." + string(name);
        end
    end
end

nanInputNames = nanNameBuffer(1:nNames);

end

function validateInputValues(data_pl, data_afs)
% validateInputValues
% Stop when an existing oxide input is Inf or a finite value below zero.
% Zero is allowed. NaN is deliberately allowed and retained.

oxideNames = {'SiO2', 'Al2O3', 'CaO', 'Na2O', 'K2O'};
maxNames = 2 .* numel(oxideNames);
invalidNameBuffer = strings(maxNames, 1);
nInvalid = 0;

for i = 1:numel(oxideNames)
    oxide = oxideNames{i};
    name = findOxideColumn(data_pl.Properties.VariableNames, oxide);
    if ~isempty(name)
        value = toScalarDouble(data_pl.(name));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            invalidNameBuffer(nInvalid) = "Plagioclase." + string(name);
        end
    end
end

for i = 1:numel(oxideNames)
    oxide = oxideNames{i};
    name = findOxideColumn(data_afs.Properties.VariableNames, oxide);
    if ~isempty(name)
        value = toScalarDouble(data_afs.(name));
        if isinf(value) || (isfinite(value) && value < 0)
            nInvalid = nInvalid + 1;
            invalidNameBuffer(nInvalid) = "Feldspar." + string(name);
        end
    end
end

if nInvalid > 0
    invalidNames = invalidNameBuffer(1:nInvalid);
    error(['Putirka2008twoFs: oxide inputs must be non-negative. ' ...
           'Negative finite or infinite value(s) were found in: ' ...
           char(strjoin(invalidNames, ', ')) '.']);
end

end

function row = calcTemp(data_pl, data_afs, P_kbar, MWinfo)
% calcTemp
% Calculate Equations (27a) and (27b) for one selected Plagioclase-Feldspar
% pair over a scalar or vector of pressure values.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

pl = prepareFeldsparRow(data_pl, MWinfo);
afs = prepareFeldsparRow(data_afs, MWinfo);

% Logarithmic terms are evaluated only for finite, strictly positive
% arguments. Zero and NaN are retained as invalid and produce NaN results.
logAbRatio = NaN;
if isfinite(afs.XAb) && afs.XAb > 0 && ...
        isfinite(pl.XAb) && pl.XAb > 0
    logAbRatio = log(afs.XAb ./ pl.XAb);
end

logAnAfs = NaN;
if isfinite(afs.XAn) && afs.XAn > 0
    logAnAfs = log(afs.XAn);
end

% Evaluate the pressure-dependent denominators as nP-by-1 vectors.
den27a = ...
    9.8 ...
    - 0.098 .* P_kbar ...
    - 2.46 .* logAbRatio ...
    - 14.2 .* afs.XSi ...
    + 423 .* afs.XCa ...
    - 2.42 .* logAnAfs ...
    - 11.4 .* (pl.XAn .* pl.XAb);

den27b = repmat( ...
    -0.11 ...
    + 0.1 .* logAbRatio ...
    - 3.27 .* afs.XAn ...
    + 0.098 .* logAnAfs ...
    + 0.52 .* (pl.XAn .* pl.XAb), ...
    nP, 1);

T_Eq27a_C = nan(nP, 1);
valid27a = isfinite(den27a) & den27a > 0;
T_Eq27a_C(valid27a) = 1e4 ./ den27a(valid27a);

num27b = -442 - 3.72 .* P_kbar;
T_Eq27b_C = nan(nP, 1);
valid27b = isfinite(den27b) & den27b ~= 0 & isfinite(num27b);
T_Eq27b_C(valid27b) = num27b(valid27b) ./ den27b(valid27b);

T_Eq27a_K = T_Eq27a_C + 273.15;
T_Eq27b_K = T_Eq27b_C + 273.15;

% Common launcher variables use Equation (27a) as the primary temperature.
T_degreeC = T_Eq27a_C;
T_deg = T_degreeC;
T_K = T_Eq27a_K;

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% Expand pressure-independent composition terms to nP rows.
row.XSi_pl = repmat(pl.XSi, nP, 1);
row.XCa_pl = repmat(pl.XCa, nP, 1);
row.XNa_pl = repmat(pl.XNa, nP, 1);
row.XK_pl = repmat(pl.XK, nP, 1);
row.XAn_pl = repmat(pl.XAn, nP, 1);
row.XAb_pl = repmat(pl.XAb, nP, 1);
row.XOr_pl = repmat(pl.XOr, nP, 1);

row.XSi_afs = repmat(afs.XSi, nP, 1);
row.XCa_afs = repmat(afs.XCa, nP, 1);
row.XNa_afs = repmat(afs.XNa, nP, 1);
row.XK_afs = repmat(afs.XK, nP, 1);
row.XAn_afs = repmat(afs.XAn, nP, 1);
row.XAb_afs = repmat(afs.XAb, nP, 1);
row.XOr_afs = repmat(afs.XOr, nP, 1);

row.log_XAb_afs_over_XAb_pl = repmat(logAbRatio, nP, 1);
row.log_XAn_afs = repmat(logAnAfs, nP, 1);
row.den_Eq27a = den27a;
row.den_Eq27b = den27b;

row.T_Eq27a_C = T_Eq27a_C;
row.T_Eq27a_K = T_Eq27a_K;
row.T_Eq27b_C = T_Eq27b_C;
row.T_Eq27b_K = T_Eq27b_K;

row.T_degreeC = T_degreeC;
row.T_deg = T_deg;
row.T_K = T_K;

end

function fs = prepareFeldsparRow(data_fs, MWinfo)
% prepareFeldsparRow
% Convert one feldspar analysis to cation fractions and An-Ab-Or fractions.
% Existing NaN values are retained. Missing optional CaO, Na2O, and K2O
% columns are assigned zero.

SiO2 = getMineralOxideRequired(data_fs, 'SiO2');
Al2O3 = getMineralOxideRequired(data_fs, 'Al2O3');
CaO = getMineralOxideOptional(data_fs, 'CaO', 0);
Na2O = getMineralOxideOptional(data_fs, 'Na2O', 0);
K2O = getMineralOxideOptional(data_fs, 'K2O', 0);

nSi = SiO2 ./ MWinfo.MW.SiO2 .* MWinfo.Cat.SiO2;
nAl = Al2O3 ./ MWinfo.MW.Al2O3 .* MWinfo.Cat.Al2O3;
nCa = CaO ./ MWinfo.MW.CaO .* MWinfo.Cat.CaO;
nNa = Na2O ./ MWinfo.MW.Na2O .* MWinfo.Cat.Na2O;
nK = K2O ./ MWinfo.MW.K2O .* MWinfo.Cat.K2O;

catSum = nSi + nAl + nCa + nNa + nK;
alkSum = nCa + nNa + nK;

% A zero or NaN normalization sum is allowed to propagate as NaN. It is not
% replaced by zero and does not stop the calculation.
fs.XSi = nSi ./ catSum;
fs.XAl = nAl ./ catSum;
fs.XCa = nCa ./ catSum;
fs.XNa = nNa ./ catSum;
fs.XK = nK ./ catSum;

fs.XAn = nCa ./ alkSum;
fs.XAb = nNa ./ alkSum;
fs.XOr = nK ./ alkSum;

end

function value = getMineralOxideRequired(dataTable, oxide)
% getMineralOxideRequired
% Read a required oxide column. The column must exist, but an existing NaN
% value is returned unchanged.

name = findOxideColumn(dataTable.Properties.VariableNames, oxide);
if isempty(name)
    error('Selected feldspar row must contain variable: %s', oxide);
end

value = toScalarDouble(dataTable.(name));

end

function value = getMineralOxideOptional(dataTable, oxide, defaultValue)
% getMineralOxideOptional
% Return defaultValue only when the optional column is absent. NaN in an
% existing column remains NaN.

name = findOxideColumn(dataTable.Properties.VariableNames, oxide);
if isempty(name)
    value = defaultValue;
else
    value = toScalarDouble(dataTable.(name));
end

end

function name = findOxideColumn(variableNames, oxide)
% findOxideColumn
% Match oxide names after removing spaces, underscores, and hyphens. A
% column with the suffix "value" is preferred over the bare oxide name.

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
    index = find(strcmp(canonicalNames, targets{i}), 1, 'first');
    if ~isempty(index)
        name = variableNames{index};
        return;
    end
end

end

function value = toScalarDouble(rawValue)
% toScalarDouble
% Convert the first table value to double. Empty, missing, non-numeric, or
% explicitly NaN content is represented as NaN and is never replaced by zero.

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

if iscell(rawValue)
    if isempty(rawValue{1})
        return;
    end

    firstValue = rawValue{1};
    if isnumeric(firstValue) || islogical(firstValue)
        value = double(firstValue(1));
    elseif isstring(firstValue)
        if ~ismissing(firstValue(1))
            value = str2double(firstValue(1));
        end
    elseif ischar(firstValue)
        value = str2double(string(firstValue));
    end
end

end
