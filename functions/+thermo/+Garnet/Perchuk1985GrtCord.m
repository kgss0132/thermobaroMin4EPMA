function results = Perchuk1985GrtCord(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet/Perchuk1985GrtCord.m
% Designed for MATLAB R2024b
%
% Garnet-cordierite Fe-Mg exchange geothermometer
% Perchuk, L.L. et al. (1985)
% Journal of Metamorphic Geology, 3, 265-310
% DOI: https://doi.org/10.1111/j.1525-1314.1985.tb00321.x
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Cordierite
% analysis (selected independently by the user from tables) and calculates
% temperature using the garnet-cordierite Fe-Mg exchange geothermometer of
% Perchuk et al. (1985).
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Grt-Crd pair and appends results into
% a single output table.
%
% Both fixed-pressure and pressure-range launchers are supported:
%   startThermoCalc_fixedP
%   startThermoCalc_rangeP
%
% -------------------------------------------------------------------------
% CALIBRATION, DEMONSTRATED APPLICATION RANGE, AND APPLICATION NOTES
%
% Perchuk et al. (1985) do not state a simple rectangular experimental
% calibration range for this thermometer. The garnet-cordierite thermometer
% is presented as Eq. (A18) in the Appendix on p. 305, where the authors
% report an estimated precision of +/-12 degreeC. Experimental calibration
% is attributed to Perchuk and Lavrent'eva (1983).
%
% The following range represents Grt-Crd applications demonstrated in the
% 1985 paper, rather than a formal calibration range:
%
%   Temperature : approximately 590-840 degreeC
%   Pressure    : approximately 3-7.5 kbar
%
% Examples include the retrograde paths discussed for Ald-49, Sut-72, and
% Sut-2/2 on p. 289, with the broader P-T interpretation summarized in
% Table 8 and Fig. 20 on pp. 295-296.
%
% Important application notes:
%   1) Apply the thermometer to texturally and chemically equilibrated
%      coexisting garnet and cordierite. Cordierite may form during either
%      prograde or retrograde metamorphism, so unrelated core-rim or
%      inclusion-matrix combinations should not be paired (p. 283).
%   2) Garnet and other Fe-Mg minerals commonly preserve diffusion zoning.
%      Rim compositions may record retrograde Fe-Mg re-equilibration rather
%      than peak metamorphic temperature (pp. 274-279 and 289-295).
%   3) The thermometer was developed for Ca-poor almandine-pyrope garnet
%      associated with cordierite. The Appendix states that the broader
%      garnet activity model is valid for 0 < XCa < 0.3, but Eq. (A18)
%      itself contains no Ca correction because cordierite-associated
%      garnets were generally very Ca-poor (pp. 303 and 306).
%   4) Altered or pinitized/sericitized cordierite should not be used because
%      alteration can modify the Fe-Mg ratio (cordierite alteration is
%      described, for example, on pp. 274-277).
%   5) Equation (A18) describes Fe2+-Mg exchange. This implementation uses
%      Fe_cation_apfu as Fe2+ and does not add Fe3_cation_apfu to the Fe-Mg
%      exchange term. Fe3_cation_apfu, when present, is retained only as an
%      output for inspection.
%   6) Pressure is required in the equation. P_kbar is converted internally
%      to bar because Eq. (A18) uses P in bar (p. 305).
%
% This implementation therefore issues non-stopping fprintf warnings when:
%   1) input pressure is outside the approximately 3-7.5 kbar range
%      demonstrated in the paper,
%   2) a finite calculated temperature is outside approximately
%      590-840 degreeC,
%   3) XCa in garnet is greater than or equal to 0.3,
%   4) a required calculation input contains NaN or zero, or
%   5) the calculated temperature is NaN or Inf.
%
% The numerical T-P warnings identify extrapolation beyond the examples in
% Perchuk et al. (1985); they must not be described as strict experimental
% calibration limits.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain one Garnet table:
%   rawdata_struct.Garnet
% or
%   rawdata_struct.Grt
%
% and one Cordierite table:
%   rawdata_struct.Cordierite
% or
%   rawdata_struct.Cord
% or
%   rawdata_struct.Crd
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. Remaining columns contain normalized
% cation data as apfu.
%
% Required Garnet variables:
%   Fe_cation_apfu       % Fe2+ used in the exchange thermometer
%   Mg_cation_apfu
%   Ca_cation_apfu       % used for compositional screening
%
% Required Cordierite variables:
%   Fe_cation_apfu       % Fe2+ used in the exchange thermometer
%   Mg_cation_apfu
%
% Optional variables retained in the output:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Si_cation_apfu
%   Al_cation_apfu
%
% Missing optional variables are represented by NaN, not zero. NaN values
% are retained as missing values and propagated through calculations that
% use them. All finite cation values must be non-negative. Negative values
% and Inf stop the calculation with an error; zero and NaN do not stop the
% calculation but are reported by fprintf warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Perchuk et al. (1985), Eq. (A18), p. 305:
%
%   T(K) = (3087 + 0.018 * P_bar) / (ln(KD) + 1.342)
%
% where
%
%   KD = (Fe2+/Mg)_garnet / (Fe2+/Mg)_cordierite
%
% and
%
%   P_bar = P_kbar * 1000
%   T(degreeC) = T(K) - 273.15
%
% Notes:
% - Fe_cation_apfu is treated as Fe2+.
% - Fe3_cation_apfu is not included in KD.
% - If KD is non-positive/non-finite, or if ln(KD) + 1.342 is non-positive
%   or non-finite, the calculated temperature is returned as NaN.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Perchuk1985GrtCord(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Cordierite tables
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Grt-Crd pair.
%

%% Input validation
% Accept a scalar for fixed-pressure calculations or a vector for
% pressure-range calculations.
if nargin < 2
    error('Perchuk1985GrtCord requires (rawdata_struct, P_kbar).');
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
% Accept the alternative field names used elsewhere in the project.
disp('=== Step 1: Preparing cation dataset ===');

if isfield(rawdata_struct, 'Garnet') && istable(rawdata_struct.Garnet)
    dataset_gt = rawdata_struct.Garnet;
elseif isfield(rawdata_struct, 'Grt') && istable(rawdata_struct.Grt)
    dataset_gt = rawdata_struct.Grt;
else
    error(['rawdata_struct must contain a garnet table as either ' ...
           'rawdata_struct.Garnet or rawdata_struct.Grt.']);
end

if isfield(rawdata_struct, 'Cordierite') && istable(rawdata_struct.Cordierite)
    dataset_crd = rawdata_struct.Cordierite;
elseif isfield(rawdata_struct, 'Cord') && istable(rawdata_struct.Cord)
    dataset_crd = rawdata_struct.Cord;
elseif isfield(rawdata_struct, 'Crd') && istable(rawdata_struct.Crd)
    dataset_crd = rawdata_struct.Crd;
else
    error(['rawdata_struct must contain a cordierite table as either ' ...
           'rawdata_struct.Cordierite, rawdata_struct.Cord, or ' ...
           'rawdata_struct.Crd.']);
end

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each result in a preallocated cell buffer and concatenate only once
% after the interactive loop. This avoids repeatedly resizing the full table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Demonstrated application range in Perchuk et al. (1985), not a formal
% rectangular experimental calibration range.
applicationT_min_degC = 590;
applicationT_max_degC = 840;
applicationP_min_kbar = 3;
applicationP_max_kbar = 7.5;

pressureOutsideApplication = ...
    P_kbar < applicationP_min_kbar | ...
    P_kbar > applicationP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
    dataCodes_gt = dataset_gt{:, 1};

    [selectedIdx_gt, ok] = listdlg( ...
        'PromptString', 'Please select the Garnet data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_gt)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_gt)
        disp('Selection canceled');
        break;
    end

    selectedCode_gt = dataCodes_gt(selectedIdx_gt);
    selectedCodeGtText = char(string(selectedCode_gt));
    disp(['Garnet selected: ' selectedCodeGtText]);

    % ----- Cordierite selection -----
    disp('=== Step 4: Selecting a data code from the list (Cordierite) ===');

    dataCodes_crd = dataset_crd{:, 1};

    [selectedIdx_crd, ok] = listdlg( ...
        'PromptString', 'Please select the Cordierite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_crd)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_crd)
        disp('Selection canceled');
        break;
    end

    selectedCode_crd = dataCodes_crd(selectedIdx_crd);
    selectedCodeCrdText = char(string(selectedCode_crd));
    disp(['Cordierite selected: ' selectedCodeCrdText]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_gt = dataset_gt(selectedIdx_gt, :);
    selectedData_crd = dataset_crd(selectedIdx_crd, :);

    % Resolve required and optional input variables. Optional missing values
    % are represented by NaN and are never silently changed to zero.
    gt = prepareGarnetRow(selectedData_gt);
    crd = prepareCordieriteRow(selectedData_crd);

    % Reject negative finite values and Inf, but retain NaN and zero.
    validateNonNegativeInputs(gt, crd);

    % Identify required inputs that will propagate NaN or produce invalid
    % Fe/Mg ratios. These warnings are printed after the result.
    nanInputNames = findNaNInputs(gt, crd);
    zeroInputNames = findZeroInputs(gt, crd);

    row = calcTemp(gt, crd, P_kbar);

    % Store selected identifiers for every pressure point.
    row.dataCode_gt = repmat(string(selectedCode_gt), height(row), 1);
    row.dataCode_crd = repmat(string(selectedCode_crd), height(row), 1);
    row = movevars(row, {'dataCode_gt', 'dataCode_crd'}, 'Before', 1);

    % Append the table block to the preallocated buffer.
    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the computed temperature range.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    disp([selectedCodeGtText ' & ' selectedCodeCrdText ': ' ...
        formatTemperatureRange(row.T_deg) ' degreeC']);

    % Pressure warning is issued only once because the same pressure vector
    % is used for every selected mineral pair in this function call.
    if any(pressureOutsideApplication) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximately 3-7.5 kbar ' ...
             'range demonstrated for Grt-Crd applications in Perchuk et al. ' ...
             '(1985). This is not a formal experimental calibration limit. ' ...
             '%d of %d pressure point(s) are outside; input range = ' ...
             '%.4g-%.4g kbar.\n'], ...
            sum(pressureOutsideApplication), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn when finite calculated temperatures lie outside the range
    % demonstrated by the Grt-Crd examples in the paper.
    finiteTemperature = isfinite(row.T_deg);
    temperatureOutsideApplication = finiteTemperature & ...
        (row.T_deg < applicationT_min_degC | ...
         row.T_deg > applicationT_max_degC);

    if any(temperatureOutsideApplication)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the approximately ' ...
             '590-840 degreeC range demonstrated for Grt-Crd applications in ' ...
             'Perchuk et al. (1985). This is not a formal experimental ' ...
             'calibration limit. %d of %d finite temperature point(s) are ' ...
             'outside; calculated finite range = %.4g-%.4g degreeC for ' ...
             '%s & %s.\n'], ...
            sum(temperatureOutsideApplication), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            selectedCodeGtText, ...
            selectedCodeCrdText);
    end

    % Warn if the garnet composition lies beyond the XCa range stated for
    % the broader activity model in the Appendix.
    finiteXCa = row.XCa_g(isfinite(row.XCa_g));
    if ~isempty(finiteXCa) && any(finiteXCa >= 0.3)
        fprintf(2, ...
            ['WARNING: Garnet XCa is greater than or equal to 0.3 for %s ' ...
             '(XCa = %.4g). The broader garnet activity model in Perchuk ' ...
             'et al. (1985) is stated for 0 < XCa < 0.3, while the Grt-Crd ' ...
             'thermometer itself was applied to generally Ca-poor garnet.\n'], ...
            selectedCodeGtText, ...
            finiteXCa(1));
    end

    % Report retained NaN and zero values without stopping calculation.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in required thermometer input(s) for ' ...
             '%s & %s: %s.\n' ...
             '         NaN values were retained, and affected temperature ' ...
             'results may remain NaN.\n'], ...
            selectedCodeGtText, ...
            selectedCodeCrdText, ...
            char(strjoin(nanInputNames, ', ')));
    end

    if ~isempty(zeroInputNames)
        fprintf(2, ...
            ['WARNING: Zero was found in required thermometer input(s) for ' ...
             '%s & %s: %s.\n' ...
             '         Zero was retained, but Fe/Mg, KD, or ln(KD) may be ' ...
             'undefined and temperature may be NaN.\n'], ...
            selectedCodeGtText, ...
            selectedCodeCrdText, ...
            char(strjoin(zeroInputNames, ', ')));
    end

    % Report retained non-finite outputs.
    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for ' ...
             '%s & %s (%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            selectedCodeGtText, ...
            selectedCodeCrdText, ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another mineral pair.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Perchuk1985GrtCord', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered table blocks once.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(gt, crd)
% findNaNInputs
% Return names of required calculation/screening inputs containing NaN.

variableLabels = { ...
    'Garnet.Fe_cation_apfu', gt.Fe2; ...
    'Garnet.Mg_cation_apfu', gt.Mg; ...
    'Garnet.Ca_cation_apfu', gt.Ca; ...
    'Cordierite.Fe_cation_apfu', crd.Fe2; ...
    'Cordierite.Mg_cation_apfu', crd.Mg};

nameBuffer = strings(size(variableLabels, 1), 1);
nNames = 0;

for i = 1:size(variableLabels, 1)
    if isnan(variableLabels{i, 2})
        nNames = nNames + 1;
        nameBuffer(nNames) = string(variableLabels{i, 1});
    end
end

nanInputNames = nameBuffer(1:nNames);

end

function zeroInputNames = findZeroInputs(gt, crd)
% findZeroInputs
% Return names of required calculation/screening inputs containing zero.

variableLabels = { ...
    'Garnet.Fe_cation_apfu', gt.Fe2; ...
    'Garnet.Mg_cation_apfu', gt.Mg; ...
    'Garnet.Ca_cation_apfu', gt.Ca; ...
    'Cordierite.Fe_cation_apfu', crd.Fe2; ...
    'Cordierite.Mg_cation_apfu', crd.Mg};

nameBuffer = strings(size(variableLabels, 1), 1);
nNames = 0;

for i = 1:size(variableLabels, 1)
    value = variableLabels{i, 2};
    if isfinite(value) && value == 0
        nNames = nNames + 1;
        nameBuffer(nNames) = string(variableLabels{i, 1});
    end
end

zeroInputNames = nameBuffer(1:nNames);

end

function validateNonNegativeInputs(gt, crd)
% validateNonNegativeInputs
% Stop calculation for negative finite values, Inf, nonnumeric values, or
% nonscalar values. NaN and finite zero are intentionally allowed.

gtFields = fieldnames(gt);
crdFields = fieldnames(crd);
nameBuffer = strings(numel(gtFields) + numel(crdFields), 1);
nNames = 0;

for i = 1:numel(gtFields)
    fieldName = gtFields{i};
    value = gt.(fieldName);

    if ~isnumeric(value) || ~isscalar(value) || isinf(value) || ...
            (isfinite(value) && value < 0)
        nNames = nNames + 1;
        nameBuffer(nNames) = "Garnet." + string(fieldName);
    end
end

for i = 1:numel(crdFields)
    fieldName = crdFields{i};
    value = crd.(fieldName);

    if ~isnumeric(value) || ~isscalar(value) || isinf(value) || ...
            (isfinite(value) && value < 0)
        nNames = nNames + 1;
        nameBuffer(nNames) = "Cordierite." + string(fieldName);
    end
end

invalidInputNames = nameBuffer(1:nNames);

if ~isempty(invalidInputNames)
    error(['Perchuk1985GrtCord: cation values must be numeric scalars ' ...
           'that are NaN or finite and non-negative. Invalid value(s) were ' ...
           'found in: ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcTemp(gt, crd, P_kbar)
% calcTemp
% Compute temperatures for one Grt-Crd pair at one or more pressures.
%
% Inputs:
%   gt       : scalar struct containing Garnet cation values
%   crd      : scalar struct containing Cordierite cation values
%   P_kbar   : pressure in kbar; scalar or column vector
%
% Output:
%   row      : table containing one row per pressure value.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

row = table();

% Fe_cation_apfu is used as Fe2+ in the Fe-Mg exchange calculation.
Fe2_g = gt.Fe2;
Mg_g = gt.Mg;
Fe2_crd = crd.Fe2;
Mg_crd = crd.Mg;

% Calculate composition indices. NaN and zero propagate naturally.
FeMg_garnet = Fe2_g ./ Mg_g;
FeMg_cordierite = Fe2_crd ./ Mg_crd;
KD = FeMg_garnet ./ FeMg_cordierite;
lnKD = log(KD);
denominator = lnKD + 1.342;

% Garnet Ca and Mn fractions are retained for compositional screening.
sumDivalent_g = gt.Fe2 + gt.Mg + gt.Mn + gt.Ca;
XCa_g = gt.Ca ./ sumDivalent_g;
XMn_g = gt.Mn ./ sumDivalent_g;

Mg_number_garnet = gt.Mg ./ (gt.Mg + gt.Fe2);
Mg_number_cordierite = crd.Mg ./ (crd.Mg + crd.Fe2);

% Equation (A18). A non-positive/non-finite denominator is considered
% physically invalid and returned as NaN for every pressure point.
if ~isfinite(denominator) || denominator <= 0
    T_K = nan(nP, 1);
else
    T_K = (3087 + 0.018 .* P_bar) ./ denominator;
end

T_deg = T_K - 273.15;

% Repeat scalar mineral and exchange values for every pressure point.
row.P_kbar = P_kbar;
row.P_bar = P_bar;

row.Fe2_g = repmat(gt.Fe2, nP, 1);
row.Fe_raw_g = repmat(gt.Fe2, nP, 1);
row.Fe3_g = repmat(gt.Fe3, nP, 1);
row.Mg_g = repmat(gt.Mg, nP, 1);
row.Mn_g = repmat(gt.Mn, nP, 1);
row.Ca_g = repmat(gt.Ca, nP, 1);
row.Si_g = repmat(gt.Si, nP, 1);
row.Al_g = repmat(gt.Al, nP, 1);

row.Fe2_crd = repmat(crd.Fe2, nP, 1);
row.Fe_raw_crd = repmat(crd.Fe2, nP, 1);
row.Fe3_crd = repmat(crd.Fe3, nP, 1);
row.Mg_crd = repmat(crd.Mg, nP, 1);
row.Mn_crd = repmat(crd.Mn, nP, 1);
row.Si_crd = repmat(crd.Si, nP, 1);
row.Al_crd = repmat(crd.Al, nP, 1);

row.FeMg_garnet = repmat(FeMg_garnet, nP, 1);
row.FeMg_cordierite = repmat(FeMg_cordierite, nP, 1);
row.KD = repmat(KD, nP, 1);
row.lnKD = repmat(lnKD, nP, 1);
row.denominator = repmat(denominator, nP, 1);
row.XCa_g = repmat(XCa_g, nP, 1);
row.XMn_g = repmat(XMn_g, nP, 1);

row.Mg_number_garnet = repmat(Mg_number_garnet, nP, 1);
row.Mg_number_cordierite = repmat(Mg_number_cordierite, nP, 1);

row.T_K = T_K;
row.T_deg = T_deg;

row.is_positive_FeMg_garnet = repmat( ...
    isfinite(FeMg_garnet) && FeMg_garnet > 0, nP, 1);
row.is_positive_FeMg_cordierite = repmat( ...
    isfinite(FeMg_cordierite) && FeMg_cordierite > 0, nP, 1);
row.is_positive_KD = repmat(isfinite(KD) && KD > 0, nP, 1);
row.is_positive_denominator = repmat( ...
    isfinite(denominator) && denominator > 0, nP, 1);
row.is_XCa_g_reasonable = repmat( ...
    isfinite(XCa_g) && XCa_g >= 0 && XCa_g < 0.3, nP, 1);
row.is_XMn_g_low = repmat( ...
    isfinite(XMn_g) && XMn_g < 0.1, nP, 1);

end

function gt = prepareGarnetRow(data_gt)
% prepareGarnetRow
% Extract one selected Garnet row. Required variables must exist, but their
% values may be NaN. Optional missing variables are represented by NaN.

if height(data_gt) ~= 1
    error('Garnet input must be a 1-row table.');
end

gt = struct();

gt.Fe2 = getVarOrError(data_gt, 'Fe_cation_apfu', 'Garnet');
gt.Mg = getVarOrError(data_gt, 'Mg_cation_apfu', 'Garnet');
gt.Ca = getVarOrError(data_gt, 'Ca_cation_apfu', 'Garnet');

gt.Fe3 = getVarOrNaN(data_gt, 'Fe3_cation_apfu');
gt.Mn = getVarOrNaN(data_gt, 'Mn_cation_apfu');
gt.Si = getVarOrNaN(data_gt, 'Si_cation_apfu');
gt.Al = getVarOrNaN(data_gt, 'Al_cation_apfu');

end

function crd = prepareCordieriteRow(data_crd)
% prepareCordieriteRow
% Extract one selected Cordierite row. Required variables must exist, but
% their values may be NaN. Optional missing variables become NaN.

if height(data_crd) ~= 1
    error('Cordierite input must be a 1-row table.');
end

crd = struct();

crd.Fe2 = getVarOrError(data_crd, 'Fe_cation_apfu', 'Cordierite');
crd.Mg = getVarOrError(data_crd, 'Mg_cation_apfu', 'Cordierite');

crd.Fe3 = getVarOrNaN(data_crd, 'Fe3_cation_apfu');
crd.Mn = getVarOrNaN(data_crd, 'Mn_cation_apfu');
crd.Si = getVarOrNaN(data_crd, 'Si_cation_apfu');
crd.Al = getVarOrNaN(data_crd, 'Al_cation_apfu');

end

function value = getVarOrError(tbl, varName, mineralLabel)
% getVarOrError
% Return a required numeric scalar. NaN is allowed and retained.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s in %s table must be a numeric scalar.', ...
        varName, mineralLabel);
end

end

function value = getVarOrNaN(tbl, varName)
% getVarOrNaN
% Return an optional numeric scalar, or NaN when the variable is absent.

if ~ismember(varName, tbl.Properties.VariableNames)
    value = NaN;
    return;
end

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s must be a numeric scalar.', varName);
end

end

function txt = formatTemperatureRange(values)
% formatTemperatureRange
% Format scalar or vector temperature output for the command window.

if isempty(values)
    txt = 'NaN';
elseif isscalar(values)
    txt = numberOrNaN(values(1));
else
    txt = [numberOrNaN(values(1)) ' to ' numberOrNaN(values(end))];
end

end

function txt = numberOrNaN(value)
% numberOrNaN
% Convert one numeric value to text while retaining NaN/Inf labels.

if isnan(value)
    txt = 'NaN';
elseif isinf(value)
    if value > 0
        txt = 'Inf';
    else
        txt = '-Inf';
    end
else
    txt = num2str(value);
end

end
