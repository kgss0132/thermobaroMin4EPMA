function results = Ridolfi2010(rawdata_struct, P_kbar)
% functions/+thermo/+Amphibole/Ridolfi2010.m
% Tested with MATLAB R2024b
%
% Single-amphibole Si*-sensitive thermometer
% Ridolfi, F., Renzulli, A., Puerini, M. (2010)
% Contributions to Mineralogy and Petrology, 160, 45–66
% DOI: https://doi.org/10.1007/s00410-009-0465-7
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Amphibole analysis and calculates
% temperature using Equation (1) of Ridolfi et al. (2010), which is based on
% the amphibole Si* compositional index.
%
% The function accepts either a scalar pressure or a pressure vector and is
% therefore compatible with both startThermoCalc_fixedP and
% startThermoCalc_rangeP. Because Equation (1) is pressure independent, the
% same calculated temperature is repeated for every supplied pressure value.
%
% The function is designed for repeated calculations. Each selected
% Amphibole analysis produces one output row per pressure value, and all
% result blocks are concatenated only once after the interactive loop.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Ridolfi et al. (2010) developed the thermometer for calcic amphiboles in
% calc-alkaline magmas, particularly subduction-related volcanic systems.
% The complete experimental compilation considered in the paper spans:
%
%   Temperature : 550–1120 degreeC
%   Pressure    : <1200 MPa (<12 kbar)
%   Redox range : approximately DeltaNNO = -1 to +5
%
% These broad limits are reported in the abstract on p. 45. However, not all
% compiled experiments were used in the regression. The experimental data
% were filtered for phase equilibrium, texture, formula balance, melt
% composition, and amphibole composition (pp. 47–52).
%
% Equation (1) was calibrated using the selected "consistent" experimental
% amphiboles. Figure 4a on p. 53 indicates an approximate calibration-data
% temperature coverage of about 800–1050 degreeC. Ridolfi et al. (2010) do
% not state these values as strict validity limits; they are used here only
% as approximate calibration-coverage limits for non-stopping warnings.
%
% The principal application developed in the paper concerns calc-alkaline
% crustal magmatic systems below approximately 1000 MPa (10 kbar; final
% remarks, p. 61). Equation (1) itself contains no pressure term. The authors
% also tested it on high-Al amphiboles synthesized at mantle pressures of
% 1–3 GPa and obtained a larger standard error of approximately 56 degreeC
% (p. 53). This implementation therefore warns outside 0–10 kbar, but does
% not stop the calculation and does not alter the calculated temperature.
%
% Important compositional and textural considerations include:
%
%   1) The intended materials are calcic amphiboles from calc-alkaline
%      magmas, especially magnesiohornblende, tschermakitic pargasite, and
%      magnesiohastingsite compositions (pp. 45, 50–51, and 57).
%
%   2) The selected "consistent" amphiboles have:
%        Al# = [6]Al / AlT <= 0.21
%      whereas many amphiboles with Al# > 0.21 were excluded from the
%      regression because of anomalous amphibole/melt compositions and
%      probable excess-water effects (pp. 50–53).
%
%   3) The volcanic calcic amphiboles considered by the authors are Mg-rich,
%      with Mg/(Mg + Fe2+) > 0.5. This implementation reports a warning when
%      a finite calculated Mg number is <= 0.5 (pp. 50–51 and 56).
%
%   4) Application to suspected xenocrysts, patchily zoned amphiboles,
%      microlites, rapidly quenched crystals, skeletal/swallow-tail/needle
%      crystals, fluid-rich late-stage amphiboles, and vein amphiboles is
%      not recommended, even when Al# is below 0.21 (pp. 47 and 56).
%      Homogeneous crystal interiors may retain pre-breakdown compositions,
%      whereas reaction rims should be treated cautiously (p. 47).
%
%   5) Ridolfi et al. (2010) calculated amphibole formulae following the IMA
%      calcic-amphibole procedure and estimated Fe3+/Fe2+ by charge balance
%      after normalization to 13 T+C-site cations (p. 50). This script keeps
%      the site-allocation approximation used in the original Ridolfi2010.m
%      implementation for reproducibility. Results may therefore differ from
%      the authors' AMP-TB spreadsheet when cation normalization or Fe3+
%      estimation differs.
%
%   6) Reported performance of Equation (1) is approximately:
%        standard error of estimate : +/-22 degreeC
%        maximum reported error     : +/-57 degreeC
%        coefficient of determination: R^2 = 0.84
%      These statistics are given with Equation (1) on p. 53.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) input pressure is outside the approximate principal application range
%      of 0–10 kbar,
%   2) a finite calculated temperature is outside the approximate Equation
%      (1) calibration-data coverage of 800–1050 degreeC,
%   3) the simple calcic-amphibole screening fails,
%   4) finite Al# is > 0.21,
%   5) finite Mg/(Mg + Fe2+) is <= 0.5,
%   6) a required thermometer input contains NaN, or
%   7) a calculated temperature is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain either:
%   rawdata_struct.Amphibole : table
% or
%   rawdata_struct.Amp       : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% Required variables:
%   Si_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu         % total Fe cations
%   Mg_cation_apfu
%   Ca_cation_apfu
%   Na_cation_apfu
%
% Optional variables retained from the original implementation:
%   K_cation_apfu
%   Mn_cation_apfu
%   Cr_cation_apfu
%   Fe3_cation_apfu
%
% If an optional COLUMN is absent, the original implementation's assumption
% of zero is retained. If a column exists but its selected VALUE is NaN, the
% NaN is preserved, propagated through all dependent calculations, and never
% replaced by zero.
%
% All finite required and optional mineral-composition values must be >= 0.
% Finite negative values and Inf are prohibited. Zero and NaN are allowed.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
%   T(degreeC) = -151.487 * Si_star + 2041
%
%   Si_star = Si + Al_4/15 - 2*Ti_4 - Al_6/2 - Ti_6/1.8
%             + Fe_3/9 + Fe_2/3.3 + Mg/26 + Ca_B/5 + Na_B/1.3
%             - Na_A/15 + vA/2.3
%
% where:
%   Al_4 = tetrahedral Al
%   Ti_4 = tetrahedral Ti
%   Al_6 = octahedral Al
%   Ti_6 = octahedral Ti
%   Ca_B, Na_B = B-site Ca and Na
%   Na_A = A-site Na
%   vA = A-site vacancy
%
% Site-allocation approximation retained from the original script:
%   T site total = 8 apfu
%     Al_4 = min(Al_total, max(0, 8 - Si))
%     Ti_4 = min(Ti_total, max(0, 8 - Si - Al_4))
%
%   C site:
%     Al_6 = Al_total - Al_4
%     Ti_6 = Ti_total - Ti_4
%
%   B site total = 2 apfu
%     Ca_B = min(Ca, 2)
%     Na_B = min(Na, max(0, 2 - Ca_B))
%
%   A site total = 1 apfu
%     Na_A = max(0, Na - Na_B)
%     K_A  = K
%     vA   = max(0, 1 - Na_A - K_A)
%
% Pressure is accepted and stored for interface compatibility, but it is not
% used in Equation (1). T_degreeC is supplied for common P-T plotting, while
% T_deg is retained as a legacy alias with identical values.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Ridolfi2010(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing an Amphibole or Amp table
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Amphibole analysis
%

%% Input validation
if nargin < 2
    error('Ridolfi2010 requires (rawdata_struct, P_kbar).');
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

if isfield(rawdata_struct, 'Amphibole') && istable(rawdata_struct.Amphibole)
    dataset_amp = rawdata_struct.Amphibole;
elseif isfield(rawdata_struct, 'Amp') && istable(rawdata_struct.Amp)
    dataset_amp = rawdata_struct.Amp;
else
    error(['rawdata_struct must contain amphibole table as either ' ...
           'rawdata_struct.Amphibole or rawdata_struct.Amp']);
end

requiredVariables = { ...
    'Si_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu'};

missingVariables = requiredVariables( ...
    ~ismember(requiredVariables, dataset_amp.Properties.VariableNames));

if ~isempty(missingVariables)
    error('Missing required Amphibole column(s): %s', ...
        strjoin(missingVariables, ', '));
end

if width(dataset_amp) < 1
    error('The Amphibole table must contain at least one column.');
end
if height(dataset_amp) < 1
    error('The Amphibole table contains no analyses.');
end

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each calculation as a table block and concatenate only once after
% the loop. This avoids repeated growth and copying of the results table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = max(16, height(dataset_amp));
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Approximate Equation (1) calibration-data coverage from Fig. 4a (p. 53).
approxT_min_degC = 800;
approxT_max_degC = 1050;

% Principal crustal application range discussed in the final remarks (p. 61).
% Equation (1) is pressure independent; these are warning limits only.
approxP_min_kbar = 0;
approxP_max_kbar = 10;

pressureOutsideApproxRange = ...
    P_kbar < approxP_min_kbar | P_kbar > approxP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3–5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % The first column is used as the displayed analysis identifier.
    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, ok] = listdlg( ...
        'PromptString', ...
        'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_amp)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

    disp('=== Step 4: Checking the selected Amphibole data ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);

    % NaN values are retained. Only variables that can influence the
    % thermometer calculation are included in this warning list.
    nanInputNames = findNaNInputs(selectedData_amp);

    % Finite negative values and Inf stop the calculation. NaN and zero are
    % intentionally allowed.
    validateCompositionInputs(selectedData_amp);

    disp('=== Step 5: Calculating the temperature ===');

    row = calcTemp(selectedData_amp, P_kbar, ...
        approxT_min_degC, approxT_max_degC, ...
        approxP_min_kbar, approxP_max_kbar);

    row.dataCode_amp = repmat(string(selectedCode_amp), height(row), 1);
    row = movevars(row, 'dataCode_amp', 'Before', 1);

    nResultBlocks = nResultBlocks + 1;

    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end

    resultBlocks{nResultBlocks} = row;

    % Echo the calculated temperature for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    disp([char(string(selectedCode_amp)) ': ' ...
        formatTemperatureRange(row.T_degreeC) ' degreeC']);

    % Pressure is not used by Equation (1), but warn once when the supplied
    % pressure lies outside the principal crustal application range.
    if any(pressureOutsideApproxRange) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the approximate principal ' ...
             'application range discussed by Ridolfi et al. (2010): about ' ...
             '0–10 kbar (<1000 MPa). %d of %d pressure point(s) are outside ' ...
             'this range; input range = %.4g–%.4g kbar. Equation (1) is ' ...
             'pressure independent, so pressure does not change the calculated ' ...
             'temperature. This is an application-range warning, not a strict ' ...
             'pressure validity boundary for the equation.\n'], ...
            sum(pressureOutsideApproxRange), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn for finite temperatures outside the approximate calibration-data
    % coverage shown by Figure 4a. NaN and Inf are handled separately below.
    printTemperatureRangeWarning( ...
        row.T_degreeC, selectedCode_amp, ...
        approxT_min_degC, approxT_max_degC);

    % Simple compositional screening warnings. These flags do not modify or
    % remove the calculated result.
    finiteCa = isfinite(row.Ca_B);
    if any(finiteCa & ~row.is_calcic_candidate)
        fprintf(2, ...
            ['WARNING: The selected analysis %s fails the simple calcic-amphibole ' ...
             'screen used in this implementation (Ca_B < 1.5 apfu). Ridolfi ' ...
             'et al. (2010) developed the calibration primarily for calcic ' ...
             'amphiboles in calc-alkaline magmas. Formal classification requires ' ...
             'a complete amphibole formula and site allocation.\n'], ...
            char(string(selectedCode_amp)));
    end

    finiteAlNumber = isfinite(row.Al_number);
    if any(finiteAlNumber & row.Al_number > 0.21)
        fprintf(2, ...
            ['WARNING: The selected analysis %s has Al# = [6]Al/AlT > 0.21. ' ...
             'Ridolfi et al. (2010) excluded many such "inconsistent" ' ...
             'amphiboles from the Equation (1) regression (pp. 50–53). ' ...
             'The calculated temperature is retained.\n'], ...
            char(string(selectedCode_amp)));
    end

    finiteMgNumber = isfinite(row.Mg_number);
    if any(finiteMgNumber & row.Mg_number <= 0.5)
        fprintf(2, ...
            ['WARNING: The selected analysis %s has Mg/(Mg + Fe2+) <= 0.5. ' ...
             'The volcanic calcic amphiboles emphasized by Ridolfi et al. ' ...
             '(2010) are Mg-rich. The calculated temperature is retained.\n'], ...
            char(string(selectedCode_amp)));
    end

    % Print a non-stopping warning when any required thermometer input is NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the thermometer input(s) for %s: %s.\n' ...
             '         NaN values were retained and propagated; they were not ' ...
             'replaced by zero. The calculated temperature may therefore be NaN.\n'], ...
            char(string(selectedCode_amp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report NaN/Inf output values without stopping calculation.
    printNonfiniteTemperatureWarning( ...
        row.T_degreeC, selectedCode_amp);

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Ridolfi2010', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered table blocks once after the loop.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_amphibole)
% findNaNInputs
% Return names of thermometer input variables that contain NaN. Missing
% optional columns are not listed because the original zero assumption is
% retained only when the entire optional column is absent.

thermometerVariables = { ...
    'Si_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu', ...
    'Fe3_cation_apfu'};

nVariables = numel(thermometerVariables);
nanFlags = false(nVariables, 1);
allNames = strings(nVariables, 1);

for i = 1:nVariables
    variableName = thermometerVariables{i};
    allNames(i) = "Amphibole." + string(variableName);

    if ismember(variableName, data_amphibole.Properties.VariableNames)
        variableValue = data_amphibole.(variableName);
        if isnumeric(variableValue) && any(isnan(variableValue(:)))
            nanFlags(i) = true;
        end
    end
end

nanInputNames = allNames(nanFlags);

end

function validateCompositionInputs(data_amphibole)
% validateCompositionInputs
% Required and present optional variables must be numeric. Inf and finite
% negative values are prohibited. NaN and zero are allowed.

variablesToCheck = { ...
    'Si_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Fe3_cation_apfu'};

nVariables = numel(variablesToCheck);
nonNumericFlags = false(nVariables, 1);
infFlags = false(nVariables, 1);
negativeFlags = false(nVariables, 1);
allNames = strings(nVariables, 1);

for i = 1:nVariables
    variableName = variablesToCheck{i};
    allNames(i) = "Amphibole." + string(variableName);

    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        continue;
    end

    variableValue = data_amphibole.(variableName);

    if ~isnumeric(variableValue)
        nonNumericFlags(i) = true;
        continue;
    end

    if any(isinf(variableValue(:)))
        infFlags(i) = true;
    end

    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        negativeFlags(i) = true;
    end
end

if any(nonNumericFlags)
    error(['Ridolfi2010: required/present mineral-composition inputs must ' ...
           'be numeric. Non-numeric variable(s): ' ...
           char(strjoin(allNames(nonNumericFlags), ', ')) '.']);
end

if any(infFlags)
    error(['Ridolfi2010: Inf is not allowed in mineral-composition inputs. ' ...
           'Inf was found in: ' ...
           char(strjoin(allNames(infFlags), ', ')) '.']);
end

if any(negativeFlags)
    error(['Ridolfi2010: finite negative mineral-composition values are ' ...
           'not allowed. Negative value(s) were found in: ' ...
           char(strjoin(allNames(negativeFlags), ', ')) '. ' ...
           'Zero and NaN values are allowed.']);
end

end

function row = calcTemp(data_amphibole, P_kbar, ...
        approxT_min_degC, approxT_max_degC, ...
        approxP_min_kbar, approxP_max_kbar)
% calcTemp
% Calculate Ridolfi et al. (2010) Equation (1) for one selected Amphibole
% analysis over a scalar or vector of pressures. Pressure is stored but does
% not enter the thermometer equation.

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

% Extract one amphibole analysis and allocate its cations to sites using the
% same approximation as the original Ridolfi2010.m script.
amp = prepareAmphiboleRow(data_amphibole);
site = calcAmphiboleSites(amp);

% Ridolfi et al. (2010), Equation (1) Si* index.
Si_star_scalar = amp.Si ...
    + site.Al_4 ./ 15 ...
    - 2 .* site.Ti_4 ...
    - site.Al_6 ./ 2 ...
    - site.Ti_6 ./ 1.8 ...
    + amp.Fe3 ./ 9 ...
    + amp.Fe2 ./ 3.3 ...
    + amp.Mg ./ 26 ...
    + site.Ca_B ./ 5 ...
    + site.Na_B ./ 1.3 ...
    - site.Na_A ./ 15 ...
    + site.vA ./ 2.3;

% Temperature equation. The coefficients are unchanged from the original
% Ridolfi2010.m implementation.
T_degreeC_scalar = -151.487 .* Si_star_scalar + 2041;
T_K_scalar = T_degreeC_scalar + 273.15;

% Diagnostic Mg and Al numbers. NaN is retained when the denominator is NaN
% or zero.
mgDenominator = amp.Mg + amp.Fe2;
if isfinite(mgDenominator) && mgDenominator > 0
    Mg_number_scalar = amp.Mg ./ mgDenominator;
else
    Mg_number_scalar = NaN;
end

alDenominator = site.Al_4 + site.Al_6;
if isfinite(alDenominator) && alDenominator > 0
    Al_number_scalar = site.Al_6 ./ alDenominator;
else
    Al_number_scalar = NaN;
end

isCalcicScalar = isfinite(site.Ca_B) && site.Ca_B >= 1.5;
isRecommendedScalar = isCalcicScalar && ...
    isfinite(Mg_number_scalar) && Mg_number_scalar > 0.5 && ...
    isfinite(Al_number_scalar) && Al_number_scalar <= 0.21;

% Replicate the pressure-independent composition and temperature results so
% that the output has one row for every supplied pressure value.
Si = repmat(amp.Si, nP, 1);
Ti = repmat(amp.Ti, nP, 1);
Al = repmat(amp.Al, nP, 1);
FeT = repmat(amp.FeT, nP, 1);
Fe2 = repmat(amp.Fe2, nP, 1);
Fe3 = repmat(amp.Fe3, nP, 1);
Mg = repmat(amp.Mg, nP, 1);
Ca = repmat(amp.Ca, nP, 1);
Na = repmat(amp.Na, nP, 1);
K = repmat(amp.K, nP, 1);
Mn = repmat(amp.Mn, nP, 1);
Cr = repmat(amp.Cr, nP, 1);

Al_4 = repmat(site.Al_4, nP, 1);
Ti_4 = repmat(site.Ti_4, nP, 1);
Al_6 = repmat(site.Al_6, nP, 1);
Ti_6 = repmat(site.Ti_6, nP, 1);
Ca_B = repmat(site.Ca_B, nP, 1);
Na_B = repmat(site.Na_B, nP, 1);
Na_A = repmat(site.Na_A, nP, 1);
K_A = repmat(site.K_A, nP, 1);
A_occ = repmat(site.A_occ, nP, 1);
vA = repmat(site.vA, nP, 1);

Mg_number = repmat(Mg_number_scalar, nP, 1);
Al_number = repmat(Al_number_scalar, nP, 1);
Si_star = repmat(Si_star_scalar, nP, 1);

T_degreeC = repmat(T_degreeC_scalar, nP, 1);
T_deg = T_degreeC;
T_K = repmat(T_K_scalar, nP, 1);

is_calcic_candidate = repmat(isCalcicScalar, nP, 1);
recommended_by_Ridolfi2010 = repmat(isRecommendedScalar, nP, 1);

isPWithinApproxRange = isfinite(P_kbar) & ...
    P_kbar >= approxP_min_kbar & P_kbar <= approxP_max_kbar;
isTWithinApproxRange = isfinite(T_degreeC) & ...
    T_degreeC >= approxT_min_degC & T_degreeC <= approxT_max_degC;

% Pack outputs. T_degreeC is the common plotting variable expected by the
% range-pressure launcher; T_deg is retained for backward compatibility.
row = table( ...
    P_kbar, P_GPa, ...
    Si, Ti, Al, FeT, Fe2, Fe3, Mg, Ca, Na, K, Mn, Cr, ...
    Al_4, Ti_4, Al_6, Ti_6, ...
    Ca_B, Na_B, Na_A, K_A, A_occ, vA, ...
    Mg_number, Al_number, Si_star, ...
    T_degreeC, T_deg, T_K, ...
    is_calcic_candidate, recommended_by_Ridolfi2010, ...
    isPWithinApproxRange, isTWithinApproxRange, ...
    'VariableNames', { ...
    'P_kbar', 'P_GPa', ...
    'Si', 'Ti', 'Al', 'FeT', 'Fe2', 'Fe3', 'Mg', 'Ca', 'Na', 'K', 'Mn', 'Cr', ...
    'Al_4', 'Ti_4', 'Al_6', 'Ti_6', ...
    'Ca_B', 'Na_B', 'Na_A', 'K_A', 'A_occ', 'vA', ...
    'Mg_number', 'Al_number', 'Si_star', ...
    'T_degreeC', 'T_deg', 'T_K', ...
    'is_calcic_candidate', 'recommended_by_Ridolfi2010', ...
    'isPWithinApproxRange', 'isTWithinApproxRange'});

end

function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Extract one-row amphibole cation data. Missing optional columns retain the
% original assumption of zero; NaN values in existing columns are preserved.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();

amp.Si = getRequiredVariable(data_amphibole, ...
    'Si_cation_apfu', 'Amphibole');
amp.Ti = getRequiredVariable(data_amphibole, ...
    'Ti_cation_apfu', 'Amphibole');
amp.Al = getRequiredVariable(data_amphibole, ...
    'Al_cation_apfu', 'Amphibole');
amp.FeT = getRequiredVariable(data_amphibole, ...
    'Fe_cation_apfu', 'Amphibole');
amp.Mg = getRequiredVariable(data_amphibole, ...
    'Mg_cation_apfu', 'Amphibole');
amp.Ca = getRequiredVariable(data_amphibole, ...
    'Ca_cation_apfu', 'Amphibole');
amp.Na = getRequiredVariable(data_amphibole, ...
    'Na_cation_apfu', 'Amphibole');

amp.K = getOptionalVariable(data_amphibole, 'K_cation_apfu', 0);
amp.Mn = getOptionalVariable(data_amphibole, 'Mn_cation_apfu', 0);
amp.Cr = getOptionalVariable(data_amphibole, 'Cr_cation_apfu', 0);
amp.Fe3 = getOptionalVariable(data_amphibole, 'Fe3_cation_apfu', 0);

% Fe2+ is calculated from total Fe and Fe3+. NaN propagates naturally.
if isfinite(amp.Fe3) && isfinite(amp.FeT) && amp.Fe3 > amp.FeT
    error(['Ridolfi2010: Fe3_cation_apfu exceeds Fe_cation_apfu. ' ...
           'This would produce a negative Fe2+ value.']);
end

amp.Fe2 = amp.FeT - amp.Fe3;

if isfinite(amp.Fe2) && amp.Fe2 < 0
    error('Ridolfi2010: calculated Fe2+ is negative.');
end

end

function site = calcAmphiboleSites(amp)
% calcAmphiboleSites
% Approximate amphibole site allocation retained from the original script.
% Helper functions explicitly preserve NaN during min/max operations.

site = struct();

% --- T site total = 8 apfu ---
T_deficit = maxPreserveNaN(0, 8 - amp.Si);
site.Al_4 = minPreserveNaN(amp.Al, T_deficit);

T_remaining = maxPreserveNaN(0, 8 - amp.Si - site.Al_4);
site.Ti_4 = minPreserveNaN(amp.Ti, T_remaining);

site.Al_6 = amp.Al - site.Al_4;
site.Ti_6 = amp.Ti - site.Ti_4;

if isfinite(site.Al_6) && site.Al_6 < -1e-10
    error('Negative octahedral Al calculated. Check cation normalization.');
end
if isfinite(site.Ti_6) && site.Ti_6 < -1e-10
    error('Negative octahedral Ti calculated. Check cation normalization.');
end

site.Al_6 = maxPreserveNaN(0, site.Al_6);
site.Ti_6 = maxPreserveNaN(0, site.Ti_6);

% --- B site total = 2 apfu ---
site.Ca_B = minPreserveNaN(amp.Ca, 2);
B_remaining = maxPreserveNaN(0, 2 - site.Ca_B);
site.Na_B = minPreserveNaN(amp.Na, B_remaining);

% --- A site total = 1 apfu ---
site.Na_A = maxPreserveNaN(0, amp.Na - site.Na_B);
site.K_A = amp.K;
site.A_occ = site.Na_A + site.K_A;

if isfinite(site.A_occ) && site.A_occ > 1 + 1e-8
    error(['Calculated A-site occupancy exceeds 1. ' ...
           'Check cation normalization.']);
end

site.vA = maxPreserveNaN(0, 1 - site.A_occ);

end

function value = getRequiredVariable(tbl, variableName, mineralLabel)
% getRequiredVariable
% Retrieve a required scalar numeric value. NaN is retained.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);

if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s must be a numeric scalar in a 1-row table.', ...
        variableName);
end

end

function value = getOptionalVariable(tbl, variableName, defaultValue)
% getOptionalVariable
% Retrieve an optional scalar numeric value. If the column is absent, use the
% stated default. If the column exists and contains NaN, retain the NaN.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);

    if ~isnumeric(value) || ~isscalar(value)
        error('Variable %s must be a numeric scalar in a 1-row table.', ...
            variableName);
    end
else
    value = defaultValue;
end

end

function value = minPreserveNaN(a, b)
% minPreserveNaN
% Scalar minimum that returns NaN when either input is NaN.

if isnan(a) || isnan(b)
    value = NaN;
else
    value = min(a, b);
end

end

function value = maxPreserveNaN(a, b)
% maxPreserveNaN
% Scalar maximum that returns NaN when either input is NaN.

if isnan(a) || isnan(b)
    value = NaN;
else
    value = max(a, b);
end

end

function printTemperatureRangeWarning(temperatureValues, selectedCode, ...
        approxT_min_degC, approxT_max_degC)
% printTemperatureRangeWarning
% Print a non-stopping warning when a finite calculated temperature lies
% outside the approximate Equation (1) calibration-data coverage.

finiteTemperature = isfinite(temperatureValues);
outsideApproxRange = finiteTemperature & ...
    (temperatureValues < approxT_min_degC | ...
     temperatureValues > approxT_max_degC);

if any(outsideApproxRange)
    finiteValues = temperatureValues(finiteTemperature);
    fprintf(2, ...
        ['WARNING: Calculated temperature for %s is outside the approximate ' ...
         'Equation (1) calibration-data coverage shown by Ridolfi et al. ' ...
         '(2010), Figure 4a: about %.0f–%.0f degreeC. %d of %d finite ' ...
         'temperature point(s) are outside this range; calculated finite ' ...
         'range = %.4g–%.4g degreeC. These are approximate data-coverage ' ...
         'limits, not strict validity limits stated by the authors.\n'], ...
        char(string(selectedCode)), ...
        approxT_min_degC, ...
        approxT_max_degC, ...
        sum(outsideApproxRange), ...
        sum(finiteTemperature), ...
        min(finiteValues), ...
        max(finiteValues));
end

end

function printNonfiniteTemperatureWarning(temperatureValues, selectedCode)
% printNonfiniteTemperatureWarning
% Print a non-stopping warning for NaN or Inf calculated temperatures.

invalidTemperature = ~isfinite(temperatureValues);

if any(invalidTemperature)
    fprintf(2, ...
        ['WARNING: Non-finite temperature value(s) were calculated for %s ' ...
         '(%d of %d points; NaN: %d, Inf: %d).\n' ...
         '         These values remain in the output table, and the calculation ' ...
         'has not been stopped.\n'], ...
        char(string(selectedCode)), ...
        sum(invalidTemperature), ...
        numel(temperatureValues), ...
        sum(isnan(temperatureValues)), ...
        sum(isinf(temperatureValues)));
end

end

function rangeText = formatTemperatureRange(temperatureValues)
% formatTemperatureRange
% Format one or more temperature values for display without altering data.

if numel(temperatureValues) == 1
    rangeText = num2str(temperatureValues, '%.1f');
    return;
end

finiteValues = temperatureValues(isfinite(temperatureValues));

if isempty(finiteValues)
    rangeText = 'NaN to NaN';
    return;
end

minimumValue = min(finiteValues);
maximumValue = max(finiteValues);

if abs(maximumValue - minimumValue) < 1e-12
    rangeText = num2str(minimumValue, '%.1f');
else
    rangeText = [num2str(minimumValue, '%.1f') ' to ' ...
        num2str(maximumValue, '%.1f')];
end

if any(~isfinite(temperatureValues))
    rangeText = [rangeText ' (including non-finite value(s))'];
end

end
