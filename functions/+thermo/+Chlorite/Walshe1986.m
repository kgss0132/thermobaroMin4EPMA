function results = Walshe1986(rawdata_struct, P_kbar)
% functions/+thermo/+Chlorite/Walshe1986.m
% Designed for MATLAB R2024b; static checks performed for this revision
%
% Six-component single-Chlorite thermodynamic thermometer/model
% Walshe, J.L. (1986)
% Economic Geology, 81, 681-703
% DOI: https://doi.org/10.2113/gsecongeo.81.3.681
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Chlorite analysis at a time and
% estimates temperature using the six-component Chlorite solid-solution
% model of Walshe (1986).
%
% The implementation follows the supplied Walshe1986.m formulation:
%   - Walshe (1986), Table 10A: structural-site calculation
%   - Walshe (1986), Table 10B: component mole fractions X1-X6
%   - Walshe (1986), Table 5: activity-composition relations
%   - reactions (23) and (32): two geothermometer constraints
%   - Gibbs-Duhem derivative consistency: third nonlinear constraint
%
% The model solves simultaneously for:
%   T(K) : temperature
%   x    : amount of initial Fe converted to Fe3+ in component 4
%   X6   : mole fraction of component 6
%
% A finite non-negative scalar or vector pressure is accepted. Therefore,
% the function can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for every supplied
% pressure value and every user-selected Chlorite analysis.
%
% IMPORTANT:
% The supplied formulation does not apply an explicit pressure correction
% to either logK relation. Pressure is retained for traceability only.
% Consequently, a pressure vector produces identical temperature estimates
% repeated for each pressure value. This behavior is reported by fprintf.
%
% -------------------------------------------------------------------------
% MODEL-CONSTRAINT RANGE AND APPLICATION NOTES
%
% Walshe (1986) did not publish a single formal temperature-pressure
% calibration range comparable to an experimental empirical thermometer.
% The model was constructed and evaluated mainly from low-pressure
% hydrothermal and geothermal Chlorites.
%
% Main temperature constraints used in model development include:
%
%   Salton Sea : 190-322 degreeC
%   Broadlands : approximately 290 degreeC
%   Creede     : approximately 268 degreeC
%   Quaama     : approximately 325-338 degreeC
%
% The principal direct parameter-constraint interval is therefore treated
% here as approximately 190-325 degreeC. These calibration examples and
% thermodynamic constraints are described on pp. 685-696, especially
% Tables 4, 6, 9, and 12 and Figs. 2-8.
%
% Published applications in Walshe (1986) produced temperatures of
% approximately 150-374 degreeC:
%
%   Jabiluka : 150-248 degreeC (Table 13, pp. 698-699)
%   Juno     : 160-374 degreeC (Table 14, pp. 699-700)
%
% The 150-374 degreeC interval is an application-example range, not a formal
% calibration range. This implementation therefore distinguishes:
%
%   190-325 degreeC : principal direct model-constraint interval
%   150-374 degreeC : broader published application-example interval
%
% Walshe (1986) also reports the following limitations:
%
%   - Above approximately 330 degreeC, the model may produce artificially
%     high oxygen-fugacity estimates (concluding comments, p. 702). This
%     caution directly concerns redox outputs, but also identifies a region
%     where model extrapolation becomes more important.
%
%   - No formal pressure calibration range is defined. Most calculations
%     use the pure-water liquid-vapor pressure or assumed pressures of
%     approximately 1 kbar, with limited comparisons at 2 kbar
%     (Tables 9, 13, and 14, pp. 690, 698, and 700).
%
%   - No attempt was made to evaluate the pressure effect on the alpha and
%     beta parameters controlling non-ideal behavior of component 6
%     (p. 695). Reaction (23) may require a pressure correction, but the
%     available volume data produced unresolved inconsistencies
%     (concluding comments, p. 702).
%
%   - The inverse calculation assumes equilibrium among Chlorite, quartz,
%     and an aqueous phase at a known or assumed pressure (pp. 683 and
%     693-695). Quartz-absent, silica-undersaturated, water-poor, or
%     disequilibrium assemblages require additional caution.
%
%   - Fine-grained Chlorite can be intergrown with other phases, and poor
%     quality probe analyses may cause significant errors (p. 681 and
%     concluding comments, p. 702). Serpentine, mixed-layer material, and
%     smectite interlayering affected some Salton Sea data (pp. 685-687).
%
%   - Exact solutions may not be obtainable for Fe-rich Chlorites
%     (concluding comments, p. 702). Solver non-convergence or a
%     non-physical solution must therefore be retained as NaN and not
%     interpreted as a valid temperature.
%
%   - The model system is SiO2-Al2O3-MgO-FeO-Fe2O3-H2O. Mn is treated as
%     Mg-equivalent in the paper (Table 9 note, p. 690). Large Ti, Ca, Na,
%     K, Cr, Ni, or other excluded components may place a Chlorite outside
%     the six-component model.
%
%   - Calculated Fe3+ and H2O are model-dependent inverse estimates rather
%     than direct measurements. Agreement with measured Fe3+ is variable,
%     particularly in some Jabiluka samples (pp. 698-699).
%
%   - Temperature precision estimated from three Broadlands probe analyses
%     is approximately +/-10-15 degreeC, but Walshe states that definitive
%     model accuracy cannot yet be established (concluding comments,
%     pp. 701-702). This precision is not an absolute accuracy estimate.
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure exceeds the low-pressure conditions directly examined
%      in the paper (>2 kbar practical reference screen),
%   2) a finite temperature is outside 190-325 degreeC,
%   3) a finite temperature is outside the broader 150-374 degreeC
%      published application-example interval,
%   4) a finite temperature exceeds approximately 330 degreeC,
%   5) a calculation input contains NaN,
%   6) fsolve does not converge to a physical low-residual solution, or
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
% Required variables, normalized as Chlorite cations on a 14-oxygen basis:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu       % total Fe
%   Mg_cation_apfu
%
% Optional calculation variable:
%   Mn_cation_apfu       % treated as Mg-equivalent when present
%
% Optional trace variables retained in the output but not used directly by
% the supplied six-component calculation:
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
% If NaN occurs in Si, Al, total Fe, Mg, or an existing Mn column, the
% nonlinear solver is not called. A NaN-filled result is returned and a
% non-stopping fprintf message is displayed.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Initial composition:
%   Si0 = Si
%   Al0 = Al
%   Fe0 = total Fe
%   Mg0 = Mg + Mn
%   CAT = Si0 + Al0 + Fe0 + Mg0
%
% Unknowns:
%   T(K), x, X6
%
% Walshe (1986), Table 10A:
%   F = 28/(28 + x)
%
%   Si_T  = F*Si0
%   Al_T  = 4 - F*(Si0 + x/2)
%   Fe3_T = F*x/2
%
%   Al_O  = F*(Si0 + Al0 + x/2) - 4
%   Fe3_O = F*x/2 + X6
%   Fe2_O = F*(Fe0 - x) - X6
%   Mg_O  = F*Mg0
%
% Walshe (1986), Table 10B is used to calculate X1-X6.
%
% Activity-composition relations:
%   Components 1-3 use ideal random site mixing (Table 5).
%   a4 = X4
%   a5 = X5
%   a6 = gamma6*X6
%
%   log10(gamma6) =
%       alpha*(1-X6) + 0.5*beta*(1-X6)^2
%
%   alpha = -1.353e4*exp(-5026/T)
%   beta  =  3.560e4*exp(-4845/T)
%
% Reaction (23):
%   log10(K23) = 1626/T - 6.542
%
% Reaction (32):
%   The supplied original Walshe1986.m approximates the equilibrium
%   relation by linear regression in 1/T through Walshe's Table 2 values.
%   This approximation is retained for reproducibility and is explicitly
%   reported in the output diagnostics.
%
% The two equilibrium residuals and the Gibbs-Duhem derivative consistency
% residual are solved with fsolve. Optimization Toolbox is required.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Walshe1986(rawdata_struct, P_kbar)
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
    error('Walshe1986 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || isempty(P_kbar) || ~isvector(P_kbar) || ...
        any(~isfinite(P_kbar(:))) || any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

if exist('fsolve', 'file') ~= 2
    error(['Walshe1986 requires fsolve (Optimization Toolbox). ' ...
           'fsolve was not found on the MATLAB path.']);
end

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

primaryT_min_degC = 190;
primaryT_max_degC = 325;
applicationT_min_degC = 150;
applicationT_max_degC = 374;
redoxCautionT_degC = 330;
lowPressureReferenceMax_kbar = 2;

pressureOutsideReference = P_kbar > lowPressureReferenceMax_kbar;
pressureWarningIssued = false;
pressureImplementationCautionIssued = false;
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
    disp([char(string(selectedCode_chl)) ': T = ' ...
        formatFiniteRange(row.T_deg) ...
        ' degreeC, success = ' char(string(row.success(1)))]);

    if ~pressureImplementationCautionIssued
        fprintf(2, ...
            ['WARNING: This Walshe1986 implementation stores P_kbar ' ...
             'for traceability but does not apply an explicit pressure ' ...
             'correction to logK. Therefore, temperature is identical at ' ...
             'all supplied pressure points. Walshe (1986) requires known ' ...
             'or assumed pressure and discusses unresolved pressure effects ' ...
             '(pp. 683, 695, and 702).\n']);
        pressureImplementationCautionIssued = true;
    end

    if any(pressureOutsideReference) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: %d of %d pressure point(s) exceed 2 kbar; input ' ...
             'range = %.4g-%.4g kbar. Walshe (1986) defines no formal ' ...
             'pressure calibration range, and the model was evaluated ' ...
             'mainly at pure-water liquid-vapor pressures and assumed ' ...
             'pressures near 1 kbar, with limited comparisons at 2 kbar ' ...
             '(Tables 9, 13, and 14; pp. 690, 698, and 700). Values above ' ...
             '2 kbar are retained as pressure extrapolations.\n'], ...
            sum(pressureOutsideReference), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    if ~applicationCautionIssued
        fprintf(2, ...
            ['WARNING: Walshe (1986) assumes equilibrium among Chlorite, ' ...
             'quartz, and an aqueous phase. This script cannot verify ' ...
             'quartz saturation, water activity, equilibrium, true 14-A ' ...
             'Chlorite identity, mixed layering, fine intergrowths, or ' ...
             'whether excluded non-six-component cations are minor. These ' ...
             'conditions must be assessed independently (pp. 681, 683, ' ...
             '685-687, 693-695, and 702).\n']);
        applicationCautionIssued = true;
    end

    finiteTemperature = isfinite(row.T_deg);

    outsidePublishedApplication = finiteTemperature & ...
        (row.T_deg < applicationT_min_degC | ...
         row.T_deg > applicationT_max_degC);

    outsidePrimaryConstraint = finiteTemperature & ...
        (row.T_deg < primaryT_min_degC | ...
         row.T_deg > primaryT_max_degC);

    if any(outsidePublishedApplication)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature for %s is outside the ' ...
             'approximately 150-374 degreeC range represented by the ' ...
             'published Walshe (1986) applications. %d of %d finite ' ...
             'pressure-row value(s) are outside; calculated finite range ' ...
             '= %.4g-%.4g degreeC. This is a strong temperature ' ...
             'extrapolation, not a formal calibration-range test.\n'], ...
            char(string(selectedCode_chl)), ...
            sum(outsidePublishedApplication), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues));
    elseif any(outsidePrimaryConstraint)
        finiteValues = row.T_deg(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature for %s is outside the ' ...
             'approximately 190-325 degreeC principal model-constraint ' ...
             'interval, although it remains within the broader published ' ...
             'application-example range of about 150-374 degreeC. ' ...
             'Calculated finite range = %.4g-%.4g degreeC ' ...
             '(pp. 685-700).\n'], ...
            char(string(selectedCode_chl)), ...
            min(finiteValues), ...
            max(finiteValues));
    end

    aboveRedoxCaution = finiteTemperature & ...
        row.T_deg > redoxCautionT_degC;

    if any(aboveRedoxCaution)
        fprintf(2, ...
            ['WARNING: Calculated temperature exceeds approximately ' ...
             '330 degreeC for %s. Walshe (1986, p. 702) reports that the ' ...
             'model may produce artificially high oxygen-fugacity ' ...
             'estimates above about 330 degreeC. Temperature is retained, ' ...
             'but Fe3+, H2O, and redox interpretations require additional ' ...
             'caution.\n'], ...
            char(string(selectedCode_chl)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Walshe (1986) input(s) for ' ...
             '%s: %s.\n' ...
             '         Existing NaN values were retained and were not ' ...
             'replaced by zero. If a required model input is NaN, the ' ...
             'solver is skipped and dependent outputs remain NaN.\n'], ...
            char(string(selectedCode_chl)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    if any(~row.success)
        fprintf(2, ...
            ['WARNING: Walshe (1986) did not obtain a converged physical ' ...
             'low-residual solution for %s. exitflag = %g; resnorm = %.4g; ' ...
             'message = %s. The diagnostic values and any NaN results are ' ...
             'retained, and the interactive calculation continues. ' ...
             'Walshe (1986, p. 702) notes that exact solutions may not be ' ...
             'obtainable for Fe-rich Chlorites.\n'], ...
            char(string(selectedCode_chl)), ...
            row.exitflag(1), ...
            row.resnorm(1), ...
            char(row.message(1)));
    end

    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite Walshe (1986) temperature values were ' ...
             'returned for %s (%d of %d pressure rows; NaN: %d, Inf: %d). ' ...
             'These values remain in the output table, and calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_chl)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Walshe1986', ...
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
    error(['Walshe1986: cation variables must be numeric scalars. ' ...
           'Invalid variable(s): ' ...
           char(strjoin(displayNames(invalidShapeMask), ', ')) '.']);
end

if any(negativeMask)
    error(['Walshe1986: finite cation values must be greater than or ' ...
           'equal to zero. Negative value(s) were found in: ' ...
           char(strjoin(displayNames(negativeMask), ', ')) '.']);
end

if any(infiniteMask)
    error(['Walshe1986: infinite cation values are not allowed in: ' ...
           char(strjoin(displayNames(infiniteMask), ', ')) '.']);
end

end

function row = calcTemp(data_chl, P_kbar)
% calcTemp
% Calculate one scalar Walshe-model solution and repeat its outputs for
% every supplied pressure value. Pressure is not used by the supplied
% equilibrium calculation.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

chl = prepareChloriteRow(data_chl);

comp0 = struct();
comp0.Si0 = chl.Si;
comp0.Al0 = chl.Al;
comp0.Fe0 = chl.FeT;
comp0.Mg0 = chl.Mg + chl.Mn;
comp0.Mg_only = chl.Mg;
comp0.Mn_only = chl.Mn;
comp0.CAT = comp0.Si0 + comp0.Al0 + comp0.Fe0 + comp0.Mg0;

sol = solveWalsheModel(comp0);
sol.site = ensureSiteFields(sol.site);

T_deg_scalar = sol.T_K - 273.15;
sumX_scalar = sum(sol.X);
isSumXCloseScalar = ...
    isfinite(sumX_scalar) && abs(sumX_scalar - 1) < 1.0e-6;

isInPrimaryTScalar = ...
    isfinite(T_deg_scalar) && ...
    T_deg_scalar >= 190 && T_deg_scalar <= 325;

isInPublishedApplicationTScalar = ...
    isfinite(T_deg_scalar) && ...
    T_deg_scalar >= 150 && T_deg_scalar <= 374;

isInLowPressureReference = P_kbar <= 2;
aboveRedoxCautionScalar = ...
    isfinite(T_deg_scalar) && T_deg_scalar > 330;

if isfinite(chl.Mg) && isfinite(chl.FeT) && ...
        (chl.Mg + chl.FeT) > 0
    MgNumberScalar = chl.Mg ./ (chl.Mg + chl.FeT);
else
    MgNumberScalar = NaN;
end

recommendedNumericalScreen = ...
    repmat(sol.success, nP, 1) & ...
    repmat(sol.isPhysical, nP, 1) & ...
    repmat(isSumXCloseScalar, nP, 1) & ...
    repmat(isInPrimaryTScalar, nP, 1) & ...
    isInLowPressureReference;

row = table();

row.P_kbar = P_kbar;
row.pressure_used_in_solver = false(nP, 1);
row.pressure_correction_applied = false(nP, 1);

row.Si = repmat(chl.Si, nP, 1);
row.Al = repmat(chl.Al, nP, 1);
row.FeT = repmat(chl.FeT, nP, 1);
row.Fe_total = repmat(chl.FeT, nP, 1);
row.Mg = repmat(chl.Mg, nP, 1);
row.Mn = repmat(chl.Mn, nP, 1);
row.Ti = repmat(chl.Ti, nP, 1);
row.Ca = repmat(chl.Ca, nP, 1);
row.Na = repmat(chl.Na, nP, 1);
row.K = repmat(chl.K, nP, 1);
row.Cr = repmat(chl.Cr, nP, 1);

row.success = repmat(sol.success, nP, 1);
row.exitflag = repmat(sol.exitflag, nP, 1);
row.message = repmat(string(sol.message), nP, 1);
row.resnorm = repmat(sol.resnorm, nP, 1);

row.T_deg = repmat(T_deg_scalar, nP, 1);
row.T_K = repmat(sol.T_K, nP, 1);

row.x = repmat(sol.x, nP, 1);
row.X6 = repmat(sol.X6, nP, 1);
row.alpha = repmat(sol.alpha, nP, 1);
row.beta = repmat(sol.beta, nP, 1);

row.Si_T = repmat(sol.site.Si_T, nP, 1);
row.Al_T = repmat(sol.site.Al_T, nP, 1);
row.Fe3_T = repmat(sol.site.Fe3_T, nP, 1);
row.Al_O = repmat(sol.site.Al_O, nP, 1);
row.Fe3_O = repmat(sol.site.Fe3_O, nP, 1);
row.Fe2_O = repmat(sol.site.Fe2_O, nP, 1);
row.Mg_O = repmat(sol.site.Mg_O, nP, 1);
row.Vac_O = repmat(sol.site.Vac_O, nP, 1);
row.OH_apfu = repmat(sol.site.OH, nP, 1);

row.Tetrahedral_total = repmat(sol.site.Tet_total, nP, 1);
row.Octahedral_total = repmat(sol.site.Oct_total, nP, 1);

row.Fe3_total_apfu = ...
    repmat(sol.site.Fe3_T + sol.site.Fe3_O, nP, 1);
row.Fe2_total_apfu = repmat(sol.site.Fe2_O, nP, 1);

row.X1 = repmat(sol.X(1), nP, 1);
row.X2 = repmat(sol.X(2), nP, 1);
row.X3 = repmat(sol.X(3), nP, 1);
row.X4 = repmat(sol.X(4), nP, 1);
row.X5 = repmat(sol.X(5), nP, 1);
row.X6_component = repmat(sol.X(6), nP, 1);

row.a1 = repmat(sol.a(1), nP, 1);
row.a2 = repmat(sol.a(2), nP, 1);
row.a3 = repmat(sol.a(3), nP, 1);
row.a4 = repmat(sol.a(4), nP, 1);
row.a5 = repmat(sol.a(5), nP, 1);
row.a6 = repmat(sol.a(6), nP, 1);

row.log10K23_calc = repmat(sol.log10K23_calc, nP, 1);
row.log10K23_eq = repmat(sol.log10K23_eq, nP, 1);
row.log10K32_calc = repmat(sol.log10K32_calc, nP, 1);
row.log10K32_eq = repmat(sol.log10K32_eq, nP, 1);

row.residual1 = repmat(sol.residuals(1), nP, 1);
row.residual2 = repmat(sol.residuals(2), nP, 1);
row.residual3 = repmat(sol.residuals(3), nP, 1);

row.is_physical = repmat(sol.isPhysical, nP, 1);
row.sum_X = repmat(sumX_scalar, nP, 1);
row.is_sumX_close_to_1 = repmat(isSumXCloseScalar, nP, 1);

row.Mg_number_totalFeBasis = repmat(MgNumberScalar, nP, 1);

row.is_in_Walshe1986_primary_T_constraint_range = ...
    repmat(isInPrimaryTScalar, nP, 1);
row.is_in_Walshe1986_published_application_T_range = ...
    repmat(isInPublishedApplicationTScalar, nP, 1);
row.is_in_Walshe1986_lowP_reference = isInLowPressureReference;
row.is_above_Walshe1986_redox_caution_T = ...
    repmat(aboveRedoxCautionScalar, nP, 1);

row.formal_temperature_calibration_range_defined = false(nP, 1);
row.formal_pressure_calibration_range_defined = false(nP, 1);
row.requires_quartz_aqueous_equilibrium_confirmation = true(nP, 1);
row.requires_true_chlorite_phase_confirmation = true(nP, 1);
row.requires_nonSixComponent_assessment = true(nP, 1);
row.reaction32_uses_table2_regression_approximation = true(nP, 1);

row.recommended_by_Walshe1986_numerical_screen = ...
    recommendedNumericalScreen;

end

function sol = solveWalsheModel(comp0)
% solveWalsheModel
% Solve simultaneously for temperature, ferric-iron conversion variable x,
% and component-6 mole fraction X6. If any required calculation input is
% NaN, return a NaN-filled solution rather than passing NaN into fsolve.

requiredValues = [comp0.Si0, comp0.Al0, comp0.Fe0, comp0.Mg0, comp0.CAT];
if any(~isfinite(requiredValues))
    sol = initSolution();
    sol.message = ['Required Walshe-model input contains NaN or another ' ...
                   'non-finite value. The result was retained as NaN.'];
    return;
end

% reaction (23): explicit calibration
A23_log10 = 1626;
B23_log10 = -6.542;

% reaction (32): regression from Walshe Table 2 values
Tfit_C = [25; 100; 200; 300; 350];
Tfit_K = Tfit_C + 273.15;
log10K32_fit = [-2.53; -1.54; -0.63; 0.02; 0.29];
pfit = polyfit(1 ./ Tfit_K, log10K32_fit, 1);
A32_log10 = pfit(1);
B32_log10 = pfit(2);

A23_ln = log(10) * A23_log10;
A32_ln = log(10) * A32_log10;

% Initial guesses. The 50-600 degreeC transform used below is only a
% numerical search bound; it is not a published Walshe (1986) calibration
% range.
T0_deg = (max(0, 4 - comp0.Si0) + 8.26e-2) ./ 4.71e-3;
T0_deg = min(max(T0_deg, 120), 380);
T0_K = T0_deg + 273.15;

x0 = min(max(0.05 * comp0.Fe0, 1e-4), max(0.8 * comp0.Fe0, 1e-3));
X60 = 0.10;

p0 = [ ...
    invSigmoid((T0_K - 323.15) ./ (873.15 - 323.15)); ...
    invSigmoid(x0 ./ max(comp0.Fe0, 1e-6)); ...
    invSigmoid(X60 ./ 0.95)];

opts = optimoptions('fsolve', ...
    'Display', 'off', ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12, ...
    'MaxIterations', 800, ...
    'MaxFunctionEvaluations', 5000);

try
    [pbest, ~, exitflag] = fsolve( ...
        @(p) residualWalshe(p, comp0, A23_log10, B23_log10, A32_log10, B32_log10, A23_ln, A32_ln), ...
        p0, opts);
catch ME
    sol = initSolution();
    sol.message = ['fsolve failed: ' ME.message];
    return;
end

[residuals, state] = residualWalshe(pbest, comp0, A23_log10, B23_log10, A32_log10, B32_log10, A23_ln, A32_ln);

sol = initSolution();
sol.exitflag = exitflag;
sol.message = state.message;
sol.T_K = state.T_K;
sol.x = state.x;
sol.X6 = state.X6;
sol.alpha = state.alpha;
sol.beta = state.beta;
sol.site = state.site;
sol.X = state.X;
sol.a = state.a;
sol.log10K23_calc = state.log10K23_calc;
sol.log10K23_eq = state.log10K23_eq;
sol.log10K32_calc = state.log10K32_calc;
sol.log10K32_eq = state.log10K32_eq;
sol.residuals = residuals(:).';
sol.resnorm = norm(residuals);
sol.isPhysical = state.isPhysical;

sol.success = (exitflag > 0) && state.isPhysical && all(isfinite(residuals)) && norm(residuals) < 1e-6;

end

function [r, state] = residualWalshe(p, comp0, A23_log10, B23_log10, A32_log10, B32_log10, A23_ln, A32_ln)

state = initState();

[T_K, x, X6] = unpackParams(p, comp0.Fe0);

state.T_K = T_K;
state.x = x;
state.X6 = X6;

try
    alpha = -1.353e4 * exp(-5026 ./ T_K);
    beta = 3.560e4 * exp(-4845 ./ T_K);

    state.alpha = alpha;
    state.beta = beta;

    s = buildState(comp0, T_K, x, X6, alpha, beta);

    state.site = s.site;
    state.X = s.X;
    state.a = s.a;
    state.isPhysical = s.isPhysical;

    if ~s.isPhysical
        r = [100; 100; 100];
        state.message = 'Non-physical state.';
        return;
    end

    lnK23_calc = log(s.a(5)) + (10 / 6) * log(s.a(1)) - 2 * log(s.a(2));
    lnK32_calc = (10 / 7) * log(s.a(6)) + (5 / 6) * log(s.a(1)) ...
               - log(s.a(2)) - (3 / 7) * log(s.a(3)) - (5 / 7) * log(s.a(4));

    lnK23_eq = log(10) * (A23_log10 ./ T_K + B23_log10);
    lnK32_eq = log(10) * (A32_log10 ./ T_K + B32_log10);

    r1 = lnK23_calc - lnK23_eq;
    r2 = lnK32_calc - lnK32_eq;

    d = calcDerivativeTerms(comp0, T_K, x, X6, alpha, beta);

    term23 = d.dlnadx(5) + (10 / 6) * d.dlnadx(1) - 2 * d.dlnadx(2);
    dTdx = -(T_K ^ 2 ./ A23_ln) .* term23;

    term32 = (10 / 7) * d.dlnadx(6) + (5 / 6) * d.dlnadx(1) ...
           - d.dlnadx(2) - (3 / 7) * d.dlnadx(3) - (5 / 7) * d.dlnadx(4);

    r3 = term32 + (A32_ln ./ T_K .^ 2) .* dTdx;

    r = [r1; r2; r3];

    state.log10K23_calc = lnK23_calc ./ log(10);
    state.log10K23_eq = lnK23_eq ./ log(10);
    state.log10K32_calc = lnK32_calc ./ log(10);
    state.log10K32_eq = lnK32_eq ./ log(10);
    state.message = 'Residuals evaluated successfully.';

catch ME
    r = [100; 100; 100];
    state.message = ['Residual evaluation failed: ' ME.message];
end

end

function d = calcDerivativeTerms(comp0, T_K, x, X6, alpha, beta)

base = buildState(comp0, T_K, x, X6, alpha, beta);
if ~base.isPhysical
    error('Base state is non-physical.');
end

hx = max(1e-6, 1e-4 * max(abs(x), 1));
h6 = max(1e-6, 1e-4 * max(abs(X6), 1));

xp = min(x + hx, 0.999999 * comp0.Fe0);
xm = max(x - hx, 1e-12);
if xp <= xm
    xp = x + 1e-6;
    xm = x - 1e-6;
end

x6p = min(X6 + h6, 0.949999);
x6m = max(X6 - h6, 1e-12);
if x6p <= x6m
    x6p = X6 + 1e-6;
    x6m = X6 - 1e-6;
end

sp = buildState(comp0, T_K, xp, X6, alpha, beta);
sm = buildState(comp0, T_K, xm, X6, alpha, beta);
s6p = buildState(comp0, T_K, x, x6p, alpha, beta);
s6m = buildState(comp0, T_K, x, x6m, alpha, beta);

if ~sp.isPhysical || ~sm.isPhysical || ~s6p.isPhysical || ~s6m.isPhysical
    error('Derivative stencil is non-physical.');
end

dln_dx_partial = (log(sp.a) - log(sm.a)) ./ (xp - xm);
dln_dX6_partial = (log(s6p.a) - log(s6m.a)) ./ (x6p - x6m);

numer = sum(base.X .* dln_dx_partial);
denom = sum(base.X .* dln_dX6_partial);

if ~isfinite(numer) || ~isfinite(denom) || abs(denom) < 1e-12
    error('Degenerate Gibbs-Duhem derivative.');
end

dX6dx = -numer ./ denom;

d = struct();
d.dX6dx = dX6dx;
d.dlnadx_partial = dln_dx_partial;
d.dlnadX6_partial = dln_dX6_partial;
d.dlnadx = dln_dx_partial + dln_dX6_partial .* dX6dx;

end

function s = buildState(comp0, T_K, x, X6, alpha, beta)

site = buildSiteState(comp0, x, X6);
X = calcComponentFractions(comp0, x, X6);
xf = calcSiteFractions(site);

k1 = 1.0;
k2 = 59.720;
k3 = 59.720;
k6 = 729.0;

a1 = k1 .* (xf.XMg_O .^ 6) .* (xf.XSi_T .^ 4);
a2 = k2 .* (xf.XMg_O .^ 5) .* (xf.XAl_O .^ 2) .* (xf.XSi_T .^ 3) .* xf.XAl_T;
a3 = k3 .* (xf.XFe2_O .^ 5) .* (xf.XAl_O .^ 2) .* (xf.XSi_T .^ 3) .* xf.XAl_T;

a4 = X(4);
a5 = X(5);

log10gamma6 = alpha .* (1 - X6) + 0.5 .* beta .* (1 - X6) .^ 2;
gamma6 = 10 .^ log10gamma6;
a6 = gamma6 .* X(6);

a = [a1, a2, a3, a4, a5, a6];

isPhysical = true;

if any(~isfinite(struct2array(site)))
    isPhysical = false;
end
if any(~isfinite(X))
    isPhysical = false;
end
if any(~isfinite(a))
    isPhysical = false;
end

% Site occupancies must be non-negative
fieldsSite = {'Si_T', 'Al_T', 'Fe3_T', 'Al_O', 'Fe3_O', 'Fe2_O', 'Mg_O', 'Vac_O', 'OH'};
for i = 1:numel(fieldsSite)
    if site.(fieldsSite{i}) < -1e-10
        isPhysical = false;
    end
end

% Exact site totals required by Table 10A
if abs(site.Tet_total - 4) > 1e-8
    isPhysical = false;
end
if abs((site.Oct_total + site.Vac_O) - 6) > 1e-8
    isPhysical = false;
end

% Activities used in logs must be positive
if any(a <= 0)
    isPhysical = false;
end

% Walshe components 4, 5, 6 must be positive because activity = mole fraction
if X(4) <= 0 || X(5) <= 0 || X(6) <= 0
    isPhysical = false;
end

s = struct();
s.site = site;
s.X = X;
s.a = a;
s.isPhysical = isPhysical;

end

function site = buildSiteState(comp0, x, X6)

% Walshe (1986) Table 10A
F = 28 ./ (28 + x);

site = struct();
site.Si_T = F .* comp0.Si0;
site.Al_T = 4 - F .* (comp0.Si0 + x ./ 2);
site.Fe3_T = F .* x ./ 2;

site.Al_O = F .* (comp0.Si0 + comp0.Al0 + x ./ 2) - 4;
site.Fe3_O = F .* x ./ 2 + X6;
site.Fe2_O = F .* (comp0.Fe0 - x) - X6;
site.Mg_O = F .* comp0.Mg0;

site.Tet_total = site.Si_T + site.Al_T + site.Fe3_T;
site.Oct_total = site.Al_O + site.Fe3_O + site.Fe2_O + site.Mg_O;
site.Vac_O = 6 - site.Oct_total;

% Walshe component 6 has one OH less than the normal 8 OH basis
site.OH = 8 - X6;

end

function X = calcComponentFractions(comp0, x, X6)

% Walshe (1986) Table 10B
F = 28 ./ (28 + x);
CAT = comp0.CAT;

X1 = F .* (5 ./ 6) .* (comp0.Mg0 ./ 5 - comp0.Al0 ./ 2 - CAT + comp0.Fe0 ./ 5 - 7 .* x ./ 10) + 50 ./ 6;
X2 = F .* (comp0.Al0 ./ 2 + CAT - comp0.Fe0 ./ 5 + 7 .* x ./ 10) - 10;
X3 = F .* (comp0.Fe0 ./ 5 - 7 .* x ./ 10) - X6;
X4 = F .* x ./ 2;
X5 = 5 - F .* CAT ./ 2;
X6c = X6;

X = [X1, X2, X3, X4, X5, X6c];

end

function xf = calcSiteFractions(site)

% Walshe Table 5 assumes mixing over 6 octahedral and 4 tetrahedral sites
xf = struct();

xf.XSi_T = site.Si_T ./ 4;
xf.XAl_T = site.Al_T ./ 4;
xf.XFe3_T = site.Fe3_T ./ 4;

xf.XAl_O = site.Al_O ./ 6;
xf.XMg_O = site.Mg_O ./ 6;
xf.XFe2_O = site.Fe2_O ./ 6;
xf.XFe3_O = site.Fe3_O ./ 6;
xf.XVac_O = site.Vac_O ./ 6;

end

function [T_K, x, X6] = unpackParams(p, Fe0)

s1 = sigmoid(p(1));
s2 = sigmoid(p(2));
s3 = sigmoid(p(3));

T_K = 323.15 + (873.15 - 323.15) .* s1;
x = max(Fe0, 1e-6) .* s2;
x = min(x, 0.999999 .* max(Fe0, 1e-6));
X6 = 0.95 .* s3;

end

function y = sigmoid(x)
y = 1 ./ (1 + exp(-x));
end

function z = invSigmoid(y)
y = min(max(y, 1e-8), 1 - 1e-8);
z = log(y ./ (1 - y));
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

function site = ensureSiteFields(site)
% ensureSiteFields
% Ensure diagnostic site fields exist even when the solver returns early.

requiredFields = { ...
    'Si_T', 'Al_T', 'Fe3_T', ...
    'Al_O', 'Fe3_O', 'Fe2_O', 'Mg_O', ...
    'Tet_total', 'Oct_total', 'Vac_O', 'OH'};

if ~isstruct(site)
    site = struct();
end

for i = 1:numel(requiredFields)
    if ~isfield(site, requiredFields{i})
        site.(requiredFields{i}) = NaN;
    end
end

end

function sol = initSolution()
% initSolution
% Return a scalar solution structure filled with diagnostic NaN values.

sol = struct();
sol.success = false;
sol.exitflag = NaN;
sol.message = '';
sol.resnorm = NaN;

sol.T_K = NaN;
sol.x = NaN;
sol.X6 = NaN;
sol.alpha = NaN;
sol.beta = NaN;

sol.site = ensureSiteFields(struct());
sol.X = nan(1, 6);
sol.a = nan(1, 6);

sol.log10K23_calc = NaN;
sol.log10K23_eq = NaN;
sol.log10K32_calc = NaN;
sol.log10K32_eq = NaN;

sol.residuals = [NaN, NaN, NaN];
sol.isPhysical = false;

end

function state = initState()
% initState
% Return a scalar residual-evaluation state filled with diagnostic NaNs.

state = struct();
state.message = '';
state.T_K = NaN;
state.x = NaN;
state.X6 = NaN;
state.alpha = NaN;
state.beta = NaN;
state.site = ensureSiteFields(struct());
state.X = nan(1, 6);
state.a = nan(1, 6);
state.log10K23_calc = NaN;
state.log10K23_eq = NaN;
state.log10K32_calc = NaN;
state.log10K32_eq = NaN;
state.isPhysical = false;

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
