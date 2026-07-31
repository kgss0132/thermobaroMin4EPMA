function results = Xie1997(rawdata_struct, P_kbar)
% functions/+thermo/+Chlorite/Xie1997.m
% Designed for MATLAB R2024b; static checks performed for this revision
%
% Fe/(Fe+Mg)-corrected Al-in-Chlorite thermometer
% Xie, X., Byerly, G.R., Ferrell, R.E. Jr. (1997)
% Contributions to Mineralogy and Petrology, 126, 275-291
% DOI: https://doi.org/10.1007/s004100050250
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Chlorite analysis at a time and
% calculates temperature using the Fe/(Fe+Mg)-corrected Al(IV) thermometer
% proposed by Xie et al. (1997).
%
% Xie et al. (1997) did not establish a wholly independent experimental
% thermometer. They modified the Cathelineau (1988) Al(IV) thermometer to
% reduce the influence of Fe-Mg and host-rock compositional variation.
%
% A finite non-negative scalar or vector pressure is accepted. Therefore,
% the function can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for every supplied
% pressure value and every user-selected Chlorite analysis.
%
% IMPORTANT:
% The Xie et al. (1997) equation contains no pressure term and the paper
% defines no formal pressure calibration range. Pressure is retained only
% for traceability. Consequently, a pressure vector produces identical
% temperature estimates repeated for each pressure value. This behavior is
% reported by fprintf.
%
% -------------------------------------------------------------------------
% FIELD-CALIBRATION CONTEXT AND APPLICATION NOTES
%
% Xie et al. (1997) emphasize that Al(IV) in Chlorite is controlled not only
% by temperature, but also by host-rock composition, coexisting Al-bearing
% minerals, and crystal-chemical constraints. In the Barberton greenstone
% belt, Chlorites formed under broadly similar greenschist-facies conditions
% produced uncorrected Cathelineau (1988) temperatures of approximately
% 109-449 degreeC, which the authors considered largely unrealistic because
% the calculated variation tracks rock type and host-rock MgO
% (pp. 287-289; Figs. 10-11).
%
% The modified thermometer is therefore most appropriate for:
%
%   - IIb trioctahedral Chlorite,
%   - basaltic to dacitic host rocks,
%   - low- to very-low-grade metamorphic or hydrothermal conditions,
%   - Chlorite that is texturally and chemically related to the event being
%     evaluated and is not a mixed-layer or mixed-phase analysis.
%
% Temperature context:
%
%   - The underlying Cathelineau (1988) field calibration is based mainly on
%     Los Azufres andesitic samples with measured temperatures of
%     approximately 130-310 degreeC (discussion on p. 287).
%
%   - Xie et al. (1997) obtained approximately 320 degreeC for Barberton
%     basaltic to dacitic samples, consistent with greenschist-facies
%     mineral assemblages and independent thermometers (pp. 288-290).
%
%   - Xie et al. (1997) do not define a new formal temperature calibration
%     range. This implementation therefore distinguishes:
%
%       130-310 degreeC : direct field-calibration interval inherited from
%                         Cathelineau (1988)
%       130-350 degreeC : practical reference interval that includes the
%                         approximately 320 degreeC Barberton application
%
%     The 130-350 degreeC interval is a screening range used here for
%     warnings; it is not a formally published Xie et al. calibration range.
%
% Pressure context:
%
%   - No formal pressure calibration range or pressure correction is given
%     by Xie et al. (1997). Pressure is not used in the equation
%     (pp. 287-290).
%
% Compositional and petrological limitations:
%
%   - The Los Azufres Chlorites used for the underlying calibration have
%     Fe/(Fe+Mg) approximately 0.24-0.37, with an average near 0.31.
%     Barberton basaltic Chlorites considered most comparable have
%     Fe/(Fe+Mg) approximately 0.32-0.48 (pp. 288-289; Figs. 10-11).
%     The combined 0.24-0.48 interval is used here as a direct-comparison
%     reference range, not as a universal calibration range.
%
%   - The full Barberton dataset spans Fe/(Fe+Mg) approximately 0.12-0.80,
%     but that broad observed interval reflects host-rock control and must
%     not be interpreted as the thermometer calibration range
%     (abstract and pp. 278-289).
%
%   - The slope of the Al(IV)-Fe/(Fe+Mg) correction may itself depend on
%     temperature; coefficient 1.33 should not be assumed universal for all
%     thermal conditions and rock types (p. 289).
%
%   - For Barberton greenschist-facies trioctahedral Chlorite, approximate
%     full-formula compositional limits are Mg = 1.5-9.2 apfu and
%     Al(IV) = 1.0-3.2 apfu on O20(OH)16 (p. 287). On the 14-oxygen
%     half-formula basis used here, these correspond approximately to:
%
%       Mg     = 0.75-4.60 apfu
%       Al(IV) = 0.50-1.60 apfu
%
%     These are compositional reference limits, not temperature-calibration
%     limits.
%
%   - The empirical thermometer should not be applied indiscriminately to
%     ultramafic, strongly silicified, B-Si metasomatized, carbonatized, or
%     locally disequilibrated Chlorites. Chlorite replacing igneous biotite
%     or hornblende, Chlorite in quartz veins, and trace Chlorite may record
%     local rather than bulk-rock equilibrium (pp. 280-282 and 288-290).
%
%   - At low temperatures, corrensite, Chlorite-smectite mixed layers,
%     saponite, and fine inclusions can contaminate EPMA analyses. EPMA alone
%     may not identify these problems, especially in diagenetic materials
%     below approximately 150-200 degreeC (pp. 287-288).
%
%   - Structural formulae in Xie et al. (1997) were calculated on an
%     O20(OH)16 basis with total Fe treated as Fe2+. The authors note that
%     this approximation is generally reasonable for metamorphic IIb
%     trioctahedral Chlorite, but Fe-rich or oxidized Chlorite may be more
%     sensitive to Fe3+ (pp. 278-279).
%
% This implementation issues non-stopping fprintf messages when:
%   1) pressure is supplied, because no formal pressure calibration range or
%      pressure term exists,
%   2) a finite calculated temperature is outside 130-310 degreeC,
%   3) a finite calculated temperature is outside the broader 130-350
%      degreeC practical reference interval,
%   4) Fe/(Fe+Mg) is outside the direct-comparison range 0.24-0.48,
%   5) Mg or Al(IV) is outside the Barberton compositional reference range,
%   6) an input contains NaN, or
%   7) a non-finite temperature is returned.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain either:
%   rawdata_struct.Chlorite : table
% or
%   rawdata_struct.Chl      : table
%
% The FIRST column is treated as an identifier ("data code") displayed in
% the selection dialog.
%
% Required variables, normalized on a 14-oxygen Chlorite basis:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu       % total Fe; treated as Fe2+ following the paper
%   Mg_cation_apfu
%
% Optional variable:
%   Mn_cation_apfu
%
% Optional trace variables retained in the output but not used directly by
% the Xie et al. temperature equation:
%   Ti_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%   Cr_cation_apfu
%
% All finite cation values must be greater than or equal to zero. Negative
% finite values are prohibited. Existing NaN values are retained and never
% replaced by zero. Optional variables are assigned zero only when their
% columns are absent.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% On the 14-oxygen basis:
%
%   Al(IV)_raw = 4 - Si
%
% For compatibility with the supplied original Xie1997.m implementation:
%
%   Al(IV) = min[Al_total, max(0, Al(IV)_raw)]
%   Al(VI) = Al_total - Al(IV)
%
% Fe ratio:
%
%   XFe = Fe/(Fe + Mg)
%
% Uncorrected Cathelineau (1988)-type temperature:
%
%   T_uncorrected(degreeC) = 321.98*Al(IV) - 61.92
%
% Xie et al. (1997) correction:
%
%   Al(IV)_corr = Al(IV) - 1.33*[XFe - 0.31]
%
%   T_Xie1997(degreeC) = 321.98*Al(IV)_corr - 61.92
%
% This unified expression is algebraically equivalent to the two piecewise
% equations printed by Xie et al. (1997, p. 289).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Xie1997(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing a Chlorite or Chl table
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per supplied pressure value for every
%             user-selected Chlorite analysis
%

%% Input validation
if nargin < 2
    error('Xie1997 requires (rawdata_struct, P_kbar).');
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
disp('=== Step 1: Preparing cation dataset ===');

if isfield(rawdata_struct, 'Chlorite') && istable(rawdata_struct.Chlorite)
    dataset_chl = rawdata_struct.Chlorite;
elseif isfield(rawdata_struct, 'Chl') && istable(rawdata_struct.Chl)
    dataset_chl = rawdata_struct.Chl;
else
    error(['rawdata_struct must contain a Chlorite table as either ' ...
           'rawdata_struct.Chlorite or rawdata_struct.Chl']);
end

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max(16, height(dataset_chl));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

baseCalibrationT_min_degC = 130;
baseCalibrationT_max_degC = 310;
practicalReferenceT_min_degC = 130;
practicalReferenceT_max_degC = 350;
XFe_reference_min = 0.24;
XFe_reference_max = 0.48;
Mg_reference_min_14O = 0.75;
Mg_reference_max_14O = 4.60;
AlIV_reference_min_14O = 0.50;
AlIV_reference_max_14O = 1.60;

pressureCautionIssued = false;
applicationCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Chlorite) ===');

while true
    dataCodes_chl = dataset_chl{:, 1};

    [selectedIdx_chl, ok] = listdlg( ...
        'PromptString', 'Please select the Chlorite data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_chl)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_chl)
        disp('Selection canceled');
        break;
    end

    selectedCode_chl = dataCodes_chl(selectedIdx_chl);
    disp(['Chlorite selected: ' char(string(selectedCode_chl))]);

    disp('=== Step 4: Calculating the temperature ===');

    selectedData_chl = dataset_chl(selectedIdx_chl, :);

    nanInputNames = findNaNInputs(selectedData_chl);
    validateNonNegativeInputs(selectedData_chl);

    row = calcTemp(selectedData_chl, P_kbar);

    row.dataCode_chl = repmat(string(selectedCode_chl), height(row), 1);
    row = movevars(row, {'dataCode_chl'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    disp([char(string(selectedCode_chl)) ': T_Xie1997 = ' ...
        formatFiniteRange(row.T_deg) ...
        ' degreeC, T_Cathelineau1988_uncorrected = ' ...
        formatFiniteRange(row.T_uncorrected_deg) ' degreeC']);

    if ~pressureCautionIssued
        fprintf(2, ...
            ['WARNING: Xie et al. (1997) define no formal pressure ' ...
             'calibration range and their thermometer equation contains no ' ...
             'pressure term. P_kbar is stored only for traceability, so ' ...
             'temperature is identical at every supplied pressure point ' ...
             '(pp. 287-290).\n']);
        pressureCautionIssued = true;
    end

    if ~applicationCautionIssued
        fprintf(2, ...
            ['WARNING: Xie et al. (1997) recommend cautious application to ' ...
             'IIb trioctahedral Chlorite from basaltic to dacitic rocks ' ...
             'under low- to very-low-grade conditions. This script cannot ' ...
             'verify host-rock type, Chlorite polytype, mixed layering, ' ...
             'local equilibrium, replacement origin, or contamination by ' ...
             'fine inclusions. These conditions must be assessed ' ...
             'independently (pp. 278-290).\n']);
        applicationCautionIssued = true;
    end

    finiteTemperature = isfinite(row.T_deg);

    outsidePracticalReference = finiteTemperature & ...
        (row.T_deg < practicalReferenceT_min_degC | ...
         row.T_deg > practicalReferenceT_max_degC);

    outsideBaseCalibration = finiteTemperature & ...
        (row.T_deg < baseCalibrationT_min_degC | ...
         row.T_deg > baseCalibrationT_max_degC);

    if any(outsidePracticalReference)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated Xie et al. (1997) temperature for %s is ' ...
             'outside the 130-350 degreeC practical reference interval ' ...
             'used by this implementation. %d of %d finite pressure-row ' ...
             'value(s) are outside; calculated finite range = %.4g-%.4g ' ...
             'degreeC. This is a strong extrapolation. Xie et al. (1997) ' ...
             'did not publish a formal independent calibration range ' ...
             '(pp. 287-290).\n'], ...
            char(string(selectedCode_chl)), ...
            sum(outsidePracticalReference), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues));
    elseif any(outsideBaseCalibration)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated Xie et al. (1997) temperature for %s is ' ...
             'outside the approximately 130-310 degreeC direct field-' ...
             'calibration interval inherited from Cathelineau (1988), but ' ...
             'remains within the broader 130-350 degreeC practical ' ...
             'reference interval that includes the approximately 320 ' ...
             'degreeC Barberton application. Calculated finite range = ' ...
             '%.4g-%.4g degreeC (pp. 287-290).\n'], ...
            char(string(selectedCode_chl)), ...
            min(finiteValues), ...
            max(finiteValues));
    end

    finiteXFe = isfinite(row.XFe_Fe_over_FeMg);
    outsideXFeReference = finiteXFe & ...
        (row.XFe_Fe_over_FeMg < XFe_reference_min | ...
         row.XFe_Fe_over_FeMg > XFe_reference_max);

    if any(outsideXFeReference)
        XFeValues = row.XFe_Fe_over_FeMg(finiteXFe);
        fprintf(2, ...
            ['WARNING: Fe/(Fe+Mg) for %s is outside the approximately ' ...
             '0.24-0.48 direct-comparison range represented by Los ' ...
             'Azufres and Barberton basaltic Chlorites; calculated finite ' ...
             'range = %.4g-%.4g. The broader Barberton observed range ' ...
             '0.12-0.80 is not a thermometer calibration range ' ...
             '(pp. 278-289).\n'], ...
            char(string(selectedCode_chl)), ...
            min(XFeValues), ...
            max(XFeValues));
    end

    finiteMg = isfinite(row.Mg);
    outsideMgReference = finiteMg & ...
        (row.Mg < Mg_reference_min_14O | ...
         row.Mg > Mg_reference_max_14O);

    finiteAlIV = isfinite(row.Al_IV);
    outsideAlIVReference = finiteAlIV & ...
        (row.Al_IV < AlIV_reference_min_14O | ...
         row.Al_IV > AlIV_reference_max_14O);

    if any(outsideMgReference) || any(outsideAlIVReference)
        fprintf(2, ...
            ['WARNING: Chlorite composition for %s is outside the ' ...
             'approximate Barberton trioctahedral-Chlorite reference ' ...
             'limits on a 14-oxygen basis (Mg = 0.75-4.60 apfu; ' ...
             'Al(IV) = 0.50-1.60 apfu). These are compositional reference ' ...
             'limits, not formal thermometer calibration limits ' ...
             '(p. 287).\n'], ...
            char(string(selectedCode_chl)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Xie et al. (1997) input(s) ' ...
             'for %s: %s.\n' ...
             '         Existing NaN values were retained and were not ' ...
             'replaced by zero. Dependent calculated values may remain ' ...
             'NaN, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_chl)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite Xie et al. (1997) temperature values ' ...
             'were returned for %s (%d of %d pressure rows; NaN: %d, ' ...
             'Inf: %d). These values remain in the output table, and the ' ...
             'calculation has not been stopped.\n'], ...
            char(string(selectedCode_chl)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Xie1997', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_chl)
% findNaNInputs
% Return names of existing Chlorite cation variables containing NaN.

variableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu', ...
    'Cr_cation_apfu'};

displayNames = "Chlorite." + string(variableNames(:));
nanMask = false(numel(variableNames), 1);

for i = 1:numel(variableNames)
    variableName = variableNames{i};

    if ismember(variableName, data_chl.Properties.VariableNames)
        variableValue = data_chl.(variableName);

        if isnumeric(variableValue)
            nanMask(i) = any(isnan(variableValue(:)));
        end
    end
end

nanInputNames = displayNames(nanMask);

end

function validateNonNegativeInputs(data_chl)
% validateNonNegativeInputs
% Reject finite negative values, infinite values, non-numeric cation
% columns, and non-scalar values. Zero and NaN are allowed.

variableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu', ...
    'Cr_cation_apfu'};

displayNames = "Chlorite." + string(variableNames(:));
negativeMask = false(numel(variableNames), 1);
infiniteMask = false(numel(variableNames), 1);
invalidShapeMask = false(numel(variableNames), 1);

for i = 1:numel(variableNames)
    variableName = variableNames{i};

    if ismember(variableName, data_chl.Properties.VariableNames)
        variableValue = data_chl.(variableName);

        if ~isnumeric(variableValue) || ~isscalar(variableValue)
            invalidShapeMask(i) = true;
            continue;
        end

        negativeMask(i) = ...
            any(isfinite(variableValue(:)) & variableValue(:) < 0);
        infiniteMask(i) = any(isinf(variableValue(:)));
    end
end

if any(invalidShapeMask)
    error(['Xie1997: cation variables must be numeric scalars. ' ...
           'Invalid variable(s): ' ...
           char(strjoin(displayNames(invalidShapeMask), ', ')) '.']);
end

if any(negativeMask)
    error(['Xie1997: finite cation values must be greater than or ' ...
           'equal to zero. Negative value(s) were found in: ' ...
           char(strjoin(displayNames(negativeMask), ', ')) '.']);
end

if any(infiniteMask)
    error(['Xie1997: infinite cation values are not allowed in: ' ...
           char(strjoin(displayNames(infiniteMask), ', ')) '.']);
end

end

function row = calcTemp(data_chl, P_kbar)
% calcTemp
% Calculate one scalar Xie et al. temperature and repeat its outputs for
% every supplied pressure value. Pressure is not used in the equation.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

chl = prepareChloriteRow(data_chl);
site = calcChloriteSites(chl);

if isfinite(chl.Fe2) && isfinite(chl.Mg) && ...
        (chl.Fe2 + chl.Mg) > 0
    XFe = chl.Fe2 ./ (chl.Fe2 + chl.Mg);
    MgNumber = chl.Mg ./ (chl.Mg + chl.Fe2);
else
    XFe = NaN;
    MgNumber = NaN;
end

if isfinite(chl.Fe2) && isfinite(chl.Mg) && isfinite(chl.Mn) && ...
        (chl.Fe2 + chl.Mg + chl.Mn) > 0
    FeRatioFeMgMn = chl.Fe2 ./ (chl.Fe2 + chl.Mg + chl.Mn);
else
    FeRatioFeMgMn = NaN;
end

T_uncorrected_deg_scalar = 321.98 .* site.Al_IV - 61.92;
T_uncorrected_K_scalar = T_uncorrected_deg_scalar + 273.15;

Al_IV_corr_scalar = site.Al_IV - 1.33 .* (XFe - 0.31);

T_deg_scalar = 321.98 .* Al_IV_corr_scalar - 61.92;
T_K_scalar = T_deg_scalar + 273.15;

deltaTScalar = T_deg_scalar - T_uncorrected_deg_scalar;

is14OTetraReasonableScalar = ...
    isfinite(chl.Si) && chl.Si >= 2.0 && chl.Si <= 4.0;
isAlIVInBGBReferenceScalar = ...
    isfinite(site.Al_IV) && site.Al_IV >= 0.50 && site.Al_IV <= 1.60;
isMgInBGBReferenceScalar = ...
    isfinite(chl.Mg) && chl.Mg >= 0.75 && chl.Mg <= 4.60;
isAlIVCorrFiniteScalar = isfinite(Al_IV_corr_scalar);
isXFeDefinedScalar = isfinite(XFe) && XFe >= 0.0 && XFe <= 1.0;
isXFeInDirectReferenceScalar = ...
    isfinite(XFe) && XFe >= 0.24 && XFe <= 0.48;
isInBaseCalibrationTScalar = ...
    isfinite(T_deg_scalar) && ...
    T_deg_scalar >= 130 && T_deg_scalar <= 310;
isInPracticalReferenceTScalar = ...
    isfinite(T_deg_scalar) && ...
    T_deg_scalar >= 130 && T_deg_scalar <= 350;
isOctahedralTotalFiniteScalar = isfinite(site.Sum_VI);

recommendedNumericalScreenScalar = ...
    is14OTetraReasonableScalar && ...
    isAlIVInBGBReferenceScalar && ...
    isMgInBGBReferenceScalar && ...
    isAlIVCorrFiniteScalar && ...
    isXFeDefinedScalar && ...
    isXFeInDirectReferenceScalar && ...
    isInPracticalReferenceTScalar && ...
    isOctahedralTotalFiniteScalar;

row = table();

row.P_kbar = P_kbar;
row.pressure_used_in_equation = false(nP, 1);
row.formal_pressure_calibration_range_defined = false(nP, 1);
row.formal_independent_temperature_calibration_range_defined = ...
    false(nP, 1);

row.Si = repmat(chl.Si, nP, 1);
row.Al = repmat(chl.Al, nP, 1);
row.FeT = repmat(chl.FeT, nP, 1);
row.Fe2 = repmat(chl.Fe2, nP, 1);
row.Mg = repmat(chl.Mg, nP, 1);
row.Mn = repmat(chl.Mn, nP, 1);
row.Ti = repmat(chl.Ti, nP, 1);
row.Ca = repmat(chl.Ca, nP, 1);
row.Na = repmat(chl.Na, nP, 1);
row.K = repmat(chl.K, nP, 1);
row.Cr = repmat(chl.Cr, nP, 1);

row.Al_IV_raw = repmat(site.Al_IV_raw, nP, 1);
row.Al_IV = repmat(site.Al_IV, nP, 1);
row.Al_VI = repmat(site.Al_VI, nP, 1);
row.Sum_VI = repmat(site.Sum_VI, nP, 1);
row.VAC = repmat(site.VAC, nP, 1);

row.XFe_Fe_over_FeMg = repmat(XFe, nP, 1);
row.Mg_number = repmat(MgNumber, nP, 1);
row.Fe_ratio_Fe_over_FeMgMn = repmat(FeRatioFeMgMn, nP, 1);

row.Al_IV_corr = repmat(Al_IV_corr_scalar, nP, 1);

row.T_uncorrected_deg = repmat(T_uncorrected_deg_scalar, nP, 1);
row.T_uncorrected_K = repmat(T_uncorrected_K_scalar, nP, 1);
row.T_deg = repmat(T_deg_scalar, nP, 1);
row.T_K = repmat(T_K_scalar, nP, 1);
row.deltaT_corrected_minus_uncorrected = repmat(deltaTScalar, nP, 1);

row.is_14O_tetra_reasonable = ...
    repmat(is14OTetraReasonableScalar, nP, 1);
row.is_AlIV_in_Barberton_reference = ...
    repmat(isAlIVInBGBReferenceScalar, nP, 1);
row.is_Mg_in_Barberton_reference = ...
    repmat(isMgInBGBReferenceScalar, nP, 1);
row.is_AlIVcorr_finite = repmat(isAlIVCorrFiniteScalar, nP, 1);
row.is_XFe_defined = repmat(isXFeDefinedScalar, nP, 1);
row.is_XFe_in_direct_reference_range = ...
    repmat(isXFeInDirectReferenceScalar, nP, 1);
row.is_in_Cathelineau1988_base_T_range = ...
    repmat(isInBaseCalibrationTScalar, nP, 1);
row.is_in_Xie1997_practical_T_reference = ...
    repmat(isInPracticalReferenceTScalar, nP, 1);
row.is_octahedral_total_finite = ...
    repmat(isOctahedralTotalFiniteScalar, nP, 1);

row.requires_basaltic_to_dacitic_host_confirmation = true(nP, 1);
row.requires_IIb_trioctahedral_confirmation = true(nP, 1);
row.requires_mixed_layer_exclusion = true(nP, 1);
row.requires_local_equilibrium_confirmation = true(nP, 1);
row.total_Fe_assumed_Fe2 = true(nP, 1);

row.recommended_by_Xie1997_numerical_screen = ...
    repmat(recommendedNumericalScreenScalar, nP, 1);

end

function chl = prepareChloriteRow(data_chl)
% prepareChloriteRow
% Extract one Chlorite analysis. Existing NaN values are retained. Optional
% variables are set to zero only when their columns are absent.

if height(data_chl) ~= 1
    error('Chlorite input must be a 1-row table.');
end

chl = struct();

chl.Si = getVarOrError(data_chl, 'Si_cation_apfu', 'Chlorite');
chl.Al = getVarOrError(data_chl, 'Al_cation_apfu', 'Chlorite');
chl.FeT = getVarOrError(data_chl, 'Fe_cation_apfu', 'Chlorite');
chl.Mg = getVarOrError(data_chl, 'Mg_cation_apfu', 'Chlorite');

chl.Mn = getVarOrZeroIfMissing(data_chl, 'Mn_cation_apfu');
chl.Ti = getVarOrZeroIfMissing(data_chl, 'Ti_cation_apfu');
chl.Ca = getVarOrZeroIfMissing(data_chl, 'Ca_cation_apfu');
chl.Na = getVarOrZeroIfMissing(data_chl, 'Na_cation_apfu');
chl.K = getVarOrZeroIfMissing(data_chl, 'K_cation_apfu');
chl.Cr = getVarOrZeroIfMissing(data_chl, 'Cr_cation_apfu');

% Xie et al. (1997) treat total Fe as Fe2+ for the structural formula.
chl.Fe2 = chl.FeT;

end

function site = calcChloriteSites(chl)
% calcChloriteSites
% Approximate site allocation on a 14-oxygen basis. NaN inputs remain NaN.

site = struct();

site.Al_IV_raw = 4 - chl.Si;

if isfinite(site.Al_IV_raw) && isfinite(chl.Al)
    site.Al_IV = min(chl.Al, max(0, site.Al_IV_raw));
    site.Al_VI = chl.Al - site.Al_IV;
else
    site.Al_IV = NaN;
    site.Al_VI = NaN;
end

if isfinite(site.Al_VI) && site.Al_VI < -1.0e-10
    site.Al_IV = NaN;
    site.Al_VI = NaN;
end

if all(isfinite([site.Al_VI, chl.Fe2, chl.Mg, chl.Mn]))
    site.Sum_VI = site.Al_VI + chl.Fe2 + chl.Mg + chl.Mn;
    site.VAC = 6 - site.Sum_VI;
else
    site.Sum_VI = NaN;
    site.VAC = NaN;
end

end

function value = getVarOrError(tbl, varName, mineralLabel)
% getVarOrError
% Retrieve a required numeric scalar. NaN is retained. Infinite and finite
% negative values are prohibited.

if ~ismember(varName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, varName);
end

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s must be a numeric scalar in a 1-row table.', varName);
end

if isinf(value)
    error('Variable %s contains an infinite value.', varName);
end

if isfinite(value) && value < 0
    error('Variable %s contains a negative value.', varName);
end

end

function value = getVarOrZeroIfMissing(tbl, varName)
% getVarOrZeroIfMissing
% Retrieve an optional numeric scalar. Assign zero only if the column is
% absent. Existing NaN is retained.

if ismember(varName, tbl.Properties.VariableNames)
    value = tbl.(varName);

    if ~isnumeric(value) || ~isscalar(value)
        error('Variable %s must be a numeric scalar in a 1-row table.', varName);
    end

    if isinf(value)
        error('Variable %s contains an infinite value.', varName);
    end

    if isfinite(value) && value < 0
        error('Variable %s contains a negative value.', varName);
    end
else
    value = 0;
end

end

function textValue = formatFiniteRange(values)
% formatFiniteRange
% Format a scalar or a finite numerical range for command-window output.

values = values(:);
finiteValues = values(isfinite(values));

if isempty(finiteValues)
    textValue = 'NaN';
elseif isscalar(finiteValues) || ...
        abs(max(finiteValues) - min(finiteValues)) <= ...
        eps(max(1, max(abs(finiteValues))))
    textValue = num2str(finiteValues(1));
else
    textValue = [num2str(min(finiteValues)) ' to ' ...
                 num2str(max(finiteValues))];
end

end
