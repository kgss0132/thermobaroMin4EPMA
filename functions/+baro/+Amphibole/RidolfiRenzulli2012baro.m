function results = RidolfiRenzulli2012baro(rawdata_struct, T_degreeC, varargin)
% functions/+baro/+Amphibole/RidolfiRenzulli2012baro.m
% Tested with MATLAB R2024b
%
% Single-amphibole empirical barometer
% Ridolfi, F. and Renzulli, A. (2012)
% Contributions to Mineralogy and Petrology, 163, 877-895
% DOI: https://doi.org/10.1007/s00410-011-0704-6
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Amphibole analysis and calculates
% pressure using Equations 1a-1e and the empirical pressure-selection
% procedure of Ridolfi and Renzulli (2012).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. Equations 1a-1e do not contain temperature; the same
% composition-dependent pressure is repeated for every input temperature,
% while temperature-dependent applicability flags are evaluated separately
% for each temperature value.
%
% The launcher-compatible pressure output is P_kbar. Equation-specific
% pressures, the intermediate P2 value, APE, and the authors' final pressure
% are also retained in MPa, kbar, and GPa where appropriate.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Calibration dataset:
% Ridolfi and Renzulli (2012) selected 61 equilibrium experimental calcic
% Amphiboles crystallized from calc-alkaline and alkaline melts over:
%
%   Temperature : 800-1130 degreeC
%   Pressure    : 130-2200 MPa (1.3-22 kbar)
%
% These limits and the dataset-selection criteria are reported on pp.
% 878-880 and in Table 1. The P-T calibration data occupy a horn-shaped,
% non-rectangular field in Fig. 3a (p. 884); the simple rectangular limits
% alone do not define a uniformly validated domain.
%
% Individual pressure equations have different calibration intervals
% (Table 3, p. 887):
%
%   Eq. 1a : 130-2200 MPa
%   Eq. 1b : 130-500 MPa
%   Eq. 1c : 130-500 MPa
%   Eq. 1d : 400-1500 MPa
%   Eq. 1e : 930-2200 MPa
%
% The equations must not be treated as interchangeable independent
% barometers. The authors require all five pressures to be calculated and
% combined using the empirical decision procedure on pp. 891-892.
% Equations 1b and 1c develop large errors above approximately 400 MPa,
% Equation 1e is unreliable below approximately 1.3 GPa, and Equation 1d
% may give negative pressures for shallow-crustal Amphiboles (pp. 889-891).
%
% Uncertainty:
% The complete pressure-selection method has an overall standard error of
% approximately 11.5% (Fig. 5f; pp. 888-891). If APE is 50% or greater,
% the final pressure is the mean of P1a and P2; the authors state that the
% uncertainty of such averaged pressures cannot be determined exactly but
% is expected to be lower than approximately 20% (pp. 891-892).
%
% Amphibole type and texture:
% The model is intended for Mg-rich calcic Amphiboles from igneous rocks.
% The experimental crystals have B-site Ca >= 1.5 apfu and
% Mg/(Mg + Fe2+) >= 0.5 and include magnesiohornblende, tschermakitic
% pargasite, pargasite, magnesiohastingsite, and kaersutite (pp. 881-883).
% This implementation provides a numerical B-site-Ca screen and a
% conservative Mg/(Mg + FeT) proxy. Because total Fe does not uniquely give
% Fe2+, the original Mg/(Mg + Fe2+) criterion cannot be verified exactly.
%
% The method is limited to magmatic Amphibole phenocrysts, including
% crystals with dehydration or breakdown rims, and to euhedral homogeneous
% crystals from sub-volcanic or plutonic rocks. It should not be applied to
% hydrothermal Amphibole veins, microlites, fast-grown or quenched zones,
% strongly heterogeneous disequilibrium domains, or subsolidus/metamorphic
% Amphiboles (pp. 891-892).
%
% Amphibole normalization and composition:
% Amphibole formulae must be calculated by the 13-cation method. SiO2,
% TiO2, Al2O3, total Fe as FeO, MgO, CaO, Na2O, and K2O are required;
% Cr2O3 and MnO are optional. Natural compositions should approximately
% match the experimental fields in Fig. 1 and the oxide ranges in Table 2
% (pp. 881-883, 891):
%
%   SiO2  : 38.8-47.9 wt%
%   TiO2  : 1.1-6.4 wt%
%   Al2O3 : 7.0-15.9 wt%
%   FeOt  : 5.9-16.9 wt%
%   MnO   : 0.1-0.6 wt%   (optional)
%   MgO   : 9.7-18.0 wt%
%   CaO   : 9.9-12.4 wt%
%   Na2O  : 1.3-3.1 wt%
%   K2O   : 0.1-2.03 wt%
%
% These oxide ranges are experimental-data envelopes, not absolute
% rejection limits. Large extrapolations should nevertheless be treated
% cautiously because equation uncertainties may increase unpredictably.
%
% P-T consistency:
% The authors recommend plotting calculated P-T results in Fig. 3 and
% accepting values within or near the experimental fields after considering
% expected uncertainty (pp. 891-892). This implementation evaluates an
% approximate polygon from the whole-data coordinates listed in the Fig. 3a
% caption and issues a non-stopping warning outside that polygon.
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 800-1130 degreeC,
%   2) finite final pressure is outside 130-2200 MPa,
%   3) a finite P-T pair is outside the approximate Fig. 3a field,
%   4) equation-specific pressures are outside their calibration intervals,
%   5) Amphibole composition is outside the Table 2 experimental ranges,
%   6) the numerical calcic or conservative Mg-rich screen is not passed,
%   7) APE is 50% or greater or is non-finite,
%   8) an equation input contains NaN,
%   9) a calculated pressure is NaN or Inf, or
%  10) a negative finite pressure is calculated.
%
% All warnings are diagnostic and do not stop the calculation. Passing all
% numerical checks does not prove magmatic origin, equilibrium, suitable
% texture, or complete applicability.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole : table
%
% The FIRST Amphibole-table column is treated as the data identifier shown
% in the selection dialog. The selected row must contain oxide wt% columns
% compatible with the following names:
%
%   Required:
%     SiO2, TiO2, Al2O3, MgO, CaO, Na2O, K2O
%     FeOt, or FeO with optional Fe2O3
%
%   Optional:
%     MnO, Cr2O3, F, Cl
%
% Direct FeOt is preferred. If FeOt is absent, FeO and Fe2O3 are converted
% to total Fe as FeO. Missing optional MnO and Cr2O3 columns are treated as
% omitted components (zero), consistent with the authors' statement that
% they are optional. An explicit NaN in any existing column is retained as
% NaN and is never replaced by zero.
%
% Finite calculation inputs must be non-negative. NaN is allowed, retained,
% propagated, and reported. Inf and finite negative values are rejected.
% Zero is retained; if it makes normalization or a later ratio undefined,
% the affected result remains NaN or Inf and a warning is printed.
%
% No Liquid composition is used by this barometer. Therefore, exclusion of
% Liq F and Cl from cationTotal_liq and Liquid NaN diagnostics is not
% applicable to this function.
%
% -------------------------------------------------------------------------
% PRESSURE FORMULATION
%
% Equations 1a-1e use total Si, Ti, Al, Fe, Mg, Ca, Na, and K apfu on the
% 13-cation basis (Table 3, p. 887). Pressure is in MPa.
%
% Authors' pressure-selection procedure (pp. 891-892):
%
%   DPdb = P1d - P1b
%   XPae = (P1a - P1e) / P1a
%
%   if P1b < 335 MPa      : P2 = P1b
%   elseif P1c < 415 MPa  : P2 = P1c
%   elseif P1d < 470 MPa  : P2 = P1c
%   elseif DPdb > 500 MPa : P2 = P1e
%   elseif DPdb > 250 MPa : P2 = P1d
%   elseif DPdb < 100 MPa : P2 = P1c
%   elseif XPae < -0.45   : P2 = P1c
%   otherwise             : P2 = mean(P1b, P1c, P1e)
%
%   APE = abs((P1a - P2) * 200 / (P1a + P2))
%
%   if APE < 50% : Pfinal = P2
%   otherwise    : Pfinal = mean(P2, P1a)
%
% NaN values are not omitted from either mean. If any required input is
% NaN, the affected pressures and final result remain NaN.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = RidolfiRenzulli2012baro(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing an Amphibole table
%   T_degreeC      : non-negative numeric scalar or vector; NaN allowed
%
% Output:
%   results : table containing one row per temperature value for every
%             selected Amphibole analysis
%

%% Input validation
if nargin < 2
    error('RidolfiRenzulli2012baro requires (rawdata_struct, T_degreeC).');
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
if ~isfield(rawdata_struct, 'Amphibole') || ...
        ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end

% varargin is retained for backward-compatible function signatures. No
% optional parameters are currently defined.
if ~isempty(varargin)
    fprintf(2, ['CAUTION: RidolfiRenzulli2012baro received optional input ' ...
        'arguments, but this implementation defines no optional parameters. ' ...
        'The extra arguments are ignored.\n']);
end

T_degreeC = T_degreeC(:);

%% 1) Retrieve Amphibole dataset
disp('=== Step 1: Preparing cation dataset ===');
dataset_amp = rawdata_struct.Amphibole;
disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store result blocks in a cell buffer and concatenate only once after the
% interactive loop, avoiding repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

calibrationT_min_degreeC = 800;
calibrationT_max_degreeC = 1130;
calibrationP_min_MPa = 130;
calibrationP_max_MPa = 2200;

% Approximate whole-data field coordinates listed in the Fig. 3a caption
% (p. 884), ordered around the experimental P-T envelope.
ptField_T_degreeC = [800; 825; 850; 930; 1040; 1130; ...
    1050; 1000; 950; 900; 850];
ptField_P_MPa = [200; 130; 130; 226; 490; 1500; ...
    2000; 2200; 1000; 500; 300];

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);

temperatureWarningIssued = false;
applicationCautionIssued = false;

% Table 2 experimental Amphibole oxide ranges, wt%.
oxideRangeNames = {'SiO2','TiO2','Al2O3','FeOt','MnO', ...
    'MgO','CaO','Na2O','K2O'};
oxideRangeMin = [38.8, 1.1, 7.0, 5.9, 0.1, 9.7, 9.9, 1.3, 0.1];
oxideRangeMax = [47.9, 6.4, 15.9, 16.9, 0.6, 18.0, 12.4, 3.1, 2.03];

disp('=== Preparing output container has been finished ===');

%% 3-4) Interactive selection loop and calculation
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, ok] = listdlg( ...
        'PromptString', ...
        'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_amp)), ...
        'ListSize', [420 360]);

    if ~ok || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

    disp('=== Step 4: Calculating the pressure ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);
    amp = prepareAmphibole13Cat(selectedData_amp);
    row = calcPressure(amp, T_degreeC, ...
        ptField_T_degreeC, ptField_P_MPa);

    row.dataCode_amphibole = ...
        repmat(string(selectedCode_amp), height(row), 1);
    row = movevars(row, 'dataCode_amphibole', 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = ...
            [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    fprintf('%s: P1a = %.6g, P1b = %.6g, P1c = %.6g, P1d = %.6g, P1e = %.6g MPa; Pfinal = %.6g MPa', ...
        char(string(selectedCode_amp)), ...
        row.P1a_MPa(1), row.P1b_MPa(1), row.P1c_MPa(1), ...
        row.P1d_MPa(1), row.P1e_MPa(1), row.P_RR2012_MPa(1));

    finiteTemperatureValues = row.T_degreeC(isfinite(row.T_degreeC));
    if isempty(finiteTemperatureValues)
        fprintf('; all %d input temperature value(s) are NaN\n', height(row));
    elseif height(row) == 1
        fprintf(' at T = %.6g degreeC\n', row.T_degreeC);
    else
        fprintf(' for finite T = %.6g to %.6g degreeC (%d rows)\n', ...
            min(finiteTemperatureValues), ...
            max(finiteTemperatureValues), ...
            height(row));
    end

    if ~applicationCautionIssued
        fprintf(2, ['CAUTION: Ridolfi and Renzulli (2012) require all five ' ...
            'pressure equations and the empirical selection procedure on ' ...
            'pp. 891-892. Individual equations must not be interpreted as ' ...
            'independent wide-range barometers. The method is intended for ' ...
            'magmatic Mg-rich calcic Amphiboles and should not be applied ' ...
            'to hydrothermal veins, microlites, quenched zones, strongly ' ...
            'disequilibrium domains, or subsolidus/metamorphic Amphiboles.\n']);
        applicationCautionIssued = true;
    end

    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteT = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ['WARNING: Input temperature is outside the 800-1130 ' ...
            'degreeC experimental calibration range of Ridolfi and ' ...
            'Renzulli (2012; pp. 878-880; Table 1). %d of %d finite ' ...
            'temperature point(s) are outside; finite input range = ' ...
            '%.6g-%.6g degreeC. Equations 1a-1e do not contain a ' ...
            'temperature term, but P-T applicability must still be checked.\n'], ...
            sum(temperatureOutsideCalibration), numel(finiteT), ...
            min(finiteT), max(finiteT));
        temperatureWarningIssued = true;
    end

    % Report exact NaN input names. Missing optional MnO and Cr2O3 are
    % omitted components and are not reported; explicit NaN values are.
    if ~isempty(amp.nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the Ridolfi and Renzulli ' ...
            '(2012) input(s) for %s: %s.\n' ...
            '         NaN values were retained and were not replaced by ' ...
            'zero; affected normalized cations and pressures remain NaN.\n'], ...
            char(string(selectedCode_amp)), ...
            char(strjoin(amp.nanInputNames, ', ')));
    end

    nanTemperatureIndices = find(isnan(T_degreeC));
    if ~isempty(nanTemperatureIndices)
        fprintf(2, ['WARNING: T_degreeC contains NaN at index/indices %s ' ...
            'for %s. Temperature is not used in Equations 1a-1e, so the ' ...
            'composition-dependent pressure is retained; temperature-based ' ...
            'applicability flags are false at those rows.\n'], ...
            char(strjoin(string(nanTemperatureIndices.'), ',')), ...
            char(string(selectedCode_amp)));
    end

    % Table 2 oxide-envelope diagnostics.
    oxideValues = [amp.raw.SiO2, amp.raw.TiO2, amp.raw.Al2O3, ...
        amp.raw.FeOt, amp.raw.MnO, amp.raw.MgO, amp.raw.CaO, ...
        amp.raw.Na2O, amp.raw.K2O];
    for iRange = 1:numel(oxideRangeNames)
        value = oxideValues(iRange);
        isOptionalOmitted = ...
            strcmp(oxideRangeNames{iRange}, 'MnO') && ...
            ~amp.raw.MnO_isPresent;
        if ~isOptionalOmitted && isfinite(value) && ...
                (value < oxideRangeMin(iRange) || ...
                 value > oxideRangeMax(iRange))
            fprintf(2, ['WARNING: Amphibole %s = %.6g wt%% is outside the ' ...
                'experimental range %.6g-%.6g wt%% in Ridolfi and ' ...
                'Renzulli (2012; Table 2, p. 883) for %s. This is a ' ...
                'non-stopping extrapolation warning.\n'], ...
                oxideRangeNames{iRange}, value, ...
                oxideRangeMin(iRange), oxideRangeMax(iRange), ...
                char(string(selectedCode_amp)));
        end
    end

    if isfinite(row.Ca_B_amp_13cat(1)) && ...
            row.Ca_B_amp_13cat(1) < 1.5
        fprintf(2, ['WARNING: Estimated B-site Ca = %.6g apfu is below ' ...
            'the calcic-Amphibole criterion BCa >= 1.5 used by Ridolfi ' ...
            'and Renzulli (2012; pp. 881-883) for %s.\n'], ...
            row.Ca_B_amp_13cat(1), char(string(selectedCode_amp)));
    end

    if isfinite(row.MgNumber_totalFe_proxy(1)) && ...
            row.MgNumber_totalFe_proxy(1) < 0.5
        fprintf(2, ['WARNING: Mg/(Mg+FeT) = %.6g is below 0.5 for %s. ' ...
            'The original Mg-rich criterion is Mg/(Mg+Fe2+) >= 0.5 ' ...
            '(Ridolfi and Renzulli, 2012; pp. 881-883). Because only ' ...
            'total Fe is available, this is a conservative proxy and ' ...
            'does not exactly determine the original criterion.\n'], ...
            row.MgNumber_totalFe_proxy(1), ...
            char(string(selectedCode_amp)));
    end

    % Equation-specific range diagnostics.
    printEquationRangeWarning(row.P1a_MPa(1), 'Equation 1a', ...
        130, 2200, selectedCode_amp);
    printEquationRangeWarning(row.P1b_MPa(1), 'Equation 1b', ...
        130, 500, selectedCode_amp);
    printEquationRangeWarning(row.P1c_MPa(1), 'Equation 1c', ...
        130, 500, selectedCode_amp);
    printEquationRangeWarning(row.P1d_MPa(1), 'Equation 1d', ...
        400, 1500, selectedCode_amp);
    printEquationRangeWarning(row.P1e_MPa(1), 'Equation 1e', ...
        930, 2200, selectedCode_amp);

    finiteFinalPressure = isfinite(row.P_RR2012_MPa);
    outsideFinalPressureRange = finiteFinalPressure & ...
        (row.P_RR2012_MPa < calibrationP_min_MPa | ...
         row.P_RR2012_MPa > calibrationP_max_MPa);

    if any(outsideFinalPressureRange)
        finiteP = row.P_RR2012_MPa(finiteFinalPressure);
        fprintf(2, ['WARNING: Final pressure is outside the 130-2200 MPa ' ...
            '(1.3-22 kbar) experimental calibration range of Ridolfi and ' ...
            'Renzulli (2012; pp. 878-880; Table 1) for %s. %d of %d ' ...
            'finite row(s) are outside; finite range = %.6g-%.6g MPa. ' ...
            'Values are retained.\n'], ...
            char(string(selectedCode_amp)), ...
            sum(outsideFinalPressureRange), sum(finiteFinalPressure), ...
            min(finiteP), max(finiteP));
    end

    finitePT = isfinite(row.T_degreeC) & isfinite(row.P_RR2012_MPa);
    outsidePTField = finitePT & ~row.isWithinApproxExperimentalPTField;
    if any(outsidePTField)
        fprintf(2, ['WARNING: %d of %d finite P-T row(s) for %s lie ' ...
            'outside the approximate whole-data field reconstructed from ' ...
            'the Fig. 3a caption of Ridolfi and Renzulli (2012; p. 884). ' ...
            'The authors recommend accepting calculated data only within ' ...
            'or near the experimental fields after considering uncertainty.\n'], ...
            sum(outsidePTField), sum(finitePT), ...
            char(string(selectedCode_amp)));
    end

    if isfinite(row.APE_percent(1)) && row.APE_percent(1) >= 50
        fprintf(2, ['WARNING: APE = %.6g%% is >= 50%% for %s. The final ' ...
            'pressure is therefore the mean of P1a and P2. Ridolfi and ' ...
            'Renzulli (2012; pp. 891-892) state that uncertainty for this ' ...
            'averaged result cannot be determined exactly and is expected ' ...
            'to be lower than approximately 20%%.\n'], ...
            row.APE_percent(1), char(string(selectedCode_amp)));
    elseif ~isfinite(row.APE_percent(1))
        fprintf(2, ['WARNING: APE is NaN or Inf for %s. The pressure-' ...
            'selection diagnostic is non-finite and the final result must ' ...
            'be treated cautiously.\n'], char(string(selectedCode_amp)));
    end

    % Retain and report non-finite and negative pressures.
    pressureVariableNames = {'P1a_MPa','P1b_MPa','P1c_MPa', ...
        'P1d_MPa','P1e_MPa','P2_MPa','P_RR2012_MPa'};
    for iPressure = 1:numel(pressureVariableNames)
        variableName = pressureVariableNames{iPressure};
        pressureValues = row.(variableName);

        invalidPressure = ~isfinite(pressureValues);
        if any(invalidPressure)
            fprintf(2, ['WARNING: %s contains non-finite pressure for %s ' ...
                '(%d of %d row(s); NaN: %d, Inf: %d). Values remain in ' ...
                'the output table.\n'], ...
                variableName, char(string(selectedCode_amp)), ...
                sum(invalidPressure), numel(pressureValues), ...
                sum(isnan(pressureValues)), sum(isinf(pressureValues)));
        end

        negativePressure = isfinite(pressureValues) & pressureValues < 0;
        if any(negativePressure)
            fprintf(2, ['WARNING: %s contains negative finite pressure for ' ...
                '%s (%d of %d row(s)). Negative values are retained for ' ...
                'diagnostic purposes and are not physical negative pressure.\n'], ...
                variableName, char(string(selectedCode_amp)), ...
                sum(negativePressure), numel(pressureValues));
        end
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another Amphibole selection?', ...
        'RidolfiRenzulli2012baro', ...
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

disp('=== RidolfiRenzulli2012baro finished ===');

end

%% ---- local functions ----
function amp = prepareAmphibole13Cat(data_amp)
% prepareAmphibole13Cat
% Extract one Amphibole oxide analysis and calculate 13-cation-normalized
% Si, Ti, Al, total Fe, Mg, Ca, Na, and K. Explicit NaN values are retained.
% Missing optional MnO and Cr2O3 columns are treated as omitted components.

if height(data_amp) ~= 1
    error('Amphibole input must be a 1-row table.');
end

raw = struct();
raw.SiO2 = getMineralOxideRequired(data_amp, {'SiO2'}, 'Amphibole.SiO2');
raw.TiO2 = getMineralOxideRequired(data_amp, {'TiO2'}, 'Amphibole.TiO2');
raw.Al2O3 = getMineralOxideRequired(data_amp, {'Al2O3'}, 'Amphibole.Al2O3');
raw.MgO = getMineralOxideRequired(data_amp, {'MgO'}, 'Amphibole.MgO');
raw.CaO = getMineralOxideRequired(data_amp, {'CaO'}, 'Amphibole.CaO');
raw.Na2O = getMineralOxideRequired(data_amp, {'Na2O'}, 'Amphibole.Na2O');
raw.K2O = getMineralOxideRequired(data_amp, {'K2O'}, 'Amphibole.K2O');

[raw.MnO, mnPresent] = getMineralOxideOptional( ...
    data_amp, {'MnO'}, 'Amphibole.MnO', 0);
[raw.Cr2O3, crPresent] = getMineralOxideOptional( ...
    data_amp, {'Cr2O3'}, 'Amphibole.Cr2O3', 0);
raw.MnO_isPresent = mnPresent;
raw.Cr2O3_isPresent = crPresent;

[raw.FeOt, feSource, feInputNames] = extractTotalFeAsFeO(data_amp);

MW = struct();
MW.SiO2 = 60.0843;
MW.TiO2 = 79.866;
MW.Al2O3 = 101.9613;
MW.FeO = 71.844;
MW.Fe2O3 = 159.6882;
MW.MnO = 70.9374;
MW.MgO = 40.3044;
MW.CaO = 56.0774;
MW.Na2O = 61.9789;
MW.K2O = 94.196;
MW.Cr2O3 = 151.9904;

cat = struct();
cat.Si = raw.SiO2 ./ MW.SiO2;
cat.Ti = raw.TiO2 ./ MW.TiO2;
cat.Al = 2 .* raw.Al2O3 ./ MW.Al2O3;
cat.Fe = raw.FeOt ./ MW.FeO;
cat.Mn = raw.MnO ./ MW.MnO;
cat.Mg = raw.MgO ./ MW.MgO;
cat.Ca = raw.CaO ./ MW.CaO;
cat.Na = 2 .* raw.Na2O ./ MW.Na2O;
cat.K = 2 .* raw.K2O ./ MW.K2O;
cat.Cr = 2 .* raw.Cr2O3 ./ MW.Cr2O3;

% 13-CNK normalization denominator excludes Ca, Na, and K.
denominator13 = cat.Si + cat.Ti + cat.Al + cat.Fe + ...
    cat.Mn + cat.Mg + cat.Cr;

if isfinite(denominator13) && denominator13 > 0
    normFactor = 13 ./ denominator13;
else
    normFactor = NaN;
end

amp = struct();
amp.Si = cat.Si .* normFactor;
amp.Ti = cat.Ti .* normFactor;
amp.Al = cat.Al .* normFactor;
amp.Fe = cat.Fe .* normFactor;
amp.Mn = cat.Mn .* normFactor;
amp.Mg = cat.Mg .* normFactor;
amp.Ca = cat.Ca .* normFactor;
amp.Na = cat.Na .* normFactor;
amp.K = cat.K .* normFactor;
amp.Cr = cat.Cr .* normFactor;
amp.normFactor = normFactor;
amp.denominator13 = denominator13;
amp.raw = raw;
amp.feSource = feSource;

% Practical site-allocation diagnostics. These do not replace a complete
% amphibole-classification routine.
if isfinite(amp.Ca)
    amp.Ca_B = min(amp.Ca, 2);
else
    amp.Ca_B = NaN;
end

if isfinite(amp.Na) && isfinite(amp.Ca_B)
    amp.Na_B = min(amp.Na, max(0, 2 - amp.Ca_B));
else
    amp.Na_B = NaN;
end

amp.CaNa_B = amp.Ca_B + amp.Na_B;
amp.isCalcic_numeric = ...
    isfinite(amp.Ca_B) && amp.Ca_B >= 1.5;

if isfinite(amp.Mg) && isfinite(amp.Fe) && (amp.Mg + amp.Fe) > 0
    amp.MgNumber_totalFe = amp.Mg ./ (amp.Mg + amp.Fe);
else
    amp.MgNumber_totalFe = NaN;
end

% Build a fixed-capacity NaN input list. Missing optional MnO and Cr2O3 do
% not enter the list; explicit NaN values in existing optional columns do.
maxNames = 12;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

requiredNames = {'SiO2','TiO2','Al2O3','MgO','CaO','Na2O','K2O'};
for i = 1:numel(requiredNames)
    fieldName = requiredNames{i};
    if isnan(raw.(fieldName))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = ...
            "Amphibole." + string(fieldName);
    end
end

if isnan(raw.FeOt)
    for i = 1:numel(feInputNames)
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = string(feInputNames{i});
    end
end

if mnPresent && isnan(raw.MnO)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Amphibole.MnO";
end
if crPresent && isnan(raw.Cr2O3)
    nNanInputs = nNanInputs + 1;
    nanInputBuffer(nNanInputs) = "Amphibole.Cr2O3";
end

amp.nanInputNames = nanInputBuffer(1:nNanInputs);

end

function [FeOt, sourceLabel, inputNames] = extractTotalFeAsFeO(data_amp)
% extractTotalFeAsFeO
% Prefer a direct total-Fe-as-FeO column. If absent, combine FeO and Fe2O3.
% Missing Fe2O3 is treated as an omitted zero component; explicit NaN is
% retained. At least FeOt, FeO, or Fe2O3 must be present.

variableNames = data_amp.Properties.VariableNames;

feOtColumn = findOxideColumn(variableNames, ...
    {'FeOt','FeOT','FeOtotal','FeO_Total','FeOtot'});
if ~isempty(feOtColumn)
    FeOt = convertToScalarDouble(data_amp.(feOtColumn));
    validateNonNegativeScalar(FeOt, 'Amphibole.FeOt');
    sourceLabel = "FeOt";
    inputNames = {'Amphibole.FeOt'};
    return;
end

feOColumn = findOxideColumn(variableNames, {'FeO'});
fe2O3Column = findOxideColumn(variableNames, {'Fe2O3'});

if isempty(feOColumn) && isempty(fe2O3Column)
    error(['Selected Amphibole row must contain FeOt, FeO, or Fe2O3 ' ...
        'for total-iron calculation.']);
end

if isempty(feOColumn)
    FeO = 0;
    feOPresent = false;
else
    FeO = convertToScalarDouble(data_amp.(feOColumn));
    validateNonNegativeScalar(FeO, 'Amphibole.FeO');
    feOPresent = true;
end

if isempty(fe2O3Column)
    Fe2O3 = 0;
    fe2O3Present = false;
else
    Fe2O3 = convertToScalarDouble(data_amp.(fe2O3Column));
    validateNonNegativeScalar(Fe2O3, 'Amphibole.Fe2O3');
    fe2O3Present = true;
end

MW_FeO = 71.844;
MW_Fe2O3 = 159.6882;
FeOt = FeO + Fe2O3 .* (2 .* MW_FeO ./ MW_Fe2O3);
sourceLabel = "FeO+Fe2O3";

inputNameBuffer = cell(2, 1);
nInputNames = 0;
if feOPresent && isnan(FeO)
    nInputNames = nInputNames + 1;
    inputNameBuffer{nInputNames} = 'Amphibole.FeO';
end
if fe2O3Present && isnan(Fe2O3)
    nInputNames = nInputNames + 1;
    inputNameBuffer{nInputNames} = 'Amphibole.Fe2O3';
end
if nInputNames == 0 && isnan(FeOt)
    inputNames = {'Amphibole.totalFe'};
else
    inputNames = inputNameBuffer(1:nInputNames);
end

end

function row = calcPressure(amp, T_degreeC, ptField_T, ptField_P)
% calcPressure
% Calculate Equations 1a-1e and apply the authors' pressure-selection
% procedure. Composition-dependent scalar results are expanded to the input
% temperature-vector length. NaN values propagate naturally.

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

Si = amp.Si;
Ti = amp.Ti;
Al = amp.Al;
Fe = amp.Fe;
Mg = amp.Mg;
Ca = amp.Ca;
Na = amp.Na;
K = amp.K;

lnP1a_scalar = ...
    125.93 ...
    - 9.5876 .* Si ...
    -10.116  .* Ti ...
    - 8.1735 .* Al ...
    - 9.2261 .* Fe ...
    - 8.7934 .* Mg ...
    - 1.6659 .* Ca ...
    + 2.4835 .* Na ...
    + 2.5192 .* K;
P1a_scalar = exp(lnP1a_scalar);

lnP1b_scalar = ...
     38.723 ...
    - 2.6957 .* Si ...
    - 2.3565 .* Ti ...
    - 1.3006 .* Al ...
    - 2.7780 .* Fe ...
    - 2.4838 .* Mg ...
    - 0.6614 .* Ca ...
    - 0.2705 .* Na ...
    + 0.1117 .* K;
P1b_scalar = exp(lnP1b_scalar);

P1c_scalar = ...
    24023 ...
    -1925.3 .* Si ...
    -1720.6 .* Ti ...
    -1478.5 .* Al ...
    -1843.2 .* Fe ...
    -1746.9 .* Mg ...
    - 158.28 .* Ca ...
    -  40.444 .* Na ...
    + 253.52 .* K;

P1d_scalar = ...
    26106 ...
    -1991.9 .* Si ...
    -3035.0 .* Ti ...
    -1472.2 .* Al ...
    -2454.8 .* Fe ...
    -2125.8 .* Mg ...
    - 830.64 .* Ca ...
    +2708.8 .* Na ...
    +2204.1 .* K;

lnP1e_scalar = ...
     26.543 ...
    - 1.2085 .* Si ...
    - 3.8593 .* Ti ...
    - 1.1054 .* Al ...
    - 2.9068 .* Fe ...
    - 2.6483 .* Mg ...
    + 0.5134 .* Ca ...
    + 2.9752 .* Na ...
    + 1.8147 .* K;
P1e_scalar = exp(lnP1e_scalar);

DPdb_scalar = P1d_scalar - P1b_scalar;
XPae_scalar = (P1a_scalar - P1e_scalar) ./ P1a_scalar;

% Authors' sequential pressure-selection procedure. Comparisons involving
% NaN are false; the final otherwise branch then returns a NaN mean without
% omitting missing values.
if P1b_scalar < 335
    P2_scalar = P1b_scalar;
    P2_source_scalar = "P1b";
elseif P1c_scalar < 415
    P2_scalar = P1c_scalar;
    P2_source_scalar = "P1c";
elseif P1d_scalar < 470
    P2_scalar = P1c_scalar;
    P2_source_scalar = "P1c_from_P1d_lt470";
elseif DPdb_scalar > 500
    P2_scalar = P1e_scalar;
    P2_source_scalar = "P1e";
elseif DPdb_scalar > 250
    P2_scalar = P1d_scalar;
    P2_source_scalar = "P1d";
elseif DPdb_scalar < 100
    P2_scalar = P1c_scalar;
    P2_source_scalar = "P1c_from_DPdb_lt100";
elseif XPae_scalar < -0.45
    P2_scalar = P1c_scalar;
    P2_source_scalar = "P1c_from_XPae";
else
    P2_scalar = mean([P1b_scalar, P1c_scalar, P1e_scalar]);
    P2_source_scalar = "mean_P1b_P1c_P1e";
end

APE_scalar = abs((P1a_scalar - P2_scalar) .* 200 ./ ...
    (P1a_scalar + P2_scalar));

if APE_scalar < 50
    Pfinal_scalar = P2_scalar;
    finalPressureSource_scalar = "P2";
else
    Pfinal_scalar = mean([P2_scalar, P1a_scalar]);
    finalPressureSource_scalar = "mean_P2_P1a";
end

% Expand scalar composition and pressure values to all input temperatures.
lnP1a = repmat(lnP1a_scalar, nT, 1);
lnP1b = repmat(lnP1b_scalar, nT, 1);
lnP1e = repmat(lnP1e_scalar, nT, 1);
P1a = repmat(P1a_scalar, nT, 1);
P1b = repmat(P1b_scalar, nT, 1);
P1c = repmat(P1c_scalar, nT, 1);
P1d = repmat(P1d_scalar, nT, 1);
P1e = repmat(P1e_scalar, nT, 1);
DPdb = repmat(DPdb_scalar, nT, 1);
XPae = repmat(XPae_scalar, nT, 1);
P2 = repmat(P2_scalar, nT, 1);
APE = repmat(APE_scalar, nT, 1);
Pfinal = repmat(Pfinal_scalar, nT, 1);
P2_source = repmat(P2_source_scalar, nT, 1);
finalPressureSource = repmat(finalPressureSource_scalar, nT, 1);

isWithinCalibrationTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 800 & T_degreeC <= 1130;
isWithinCalibrationPRange = ...
    isfinite(Pfinal) & Pfinal >= 130 & Pfinal <= 2200;

isWithinApproxExperimentalPTField = false(nT, 1);
finitePT = isfinite(T_degreeC) & isfinite(Pfinal);
if any(finitePT)
    isWithinApproxExperimentalPTField(finitePT) = inpolygon( ...
        T_degreeC(finitePT), Pfinal(finitePT), ptField_T, ptField_P);
end

allEquationInputsFinite_scalar = all(isfinite( ...
    [Si, Ti, Al, Fe, Mg, Ca, Na, K]));

isApplicable_numeric = ...
    repmat(allEquationInputsFinite_scalar, nT, 1) & ...
    repmat(amp.isCalcic_numeric, nT, 1) & ...
    isWithinCalibrationTRange & ...
    isWithinCalibrationPRange & ...
    isWithinApproxExperimentalPTField;

row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.temperatureUsedInPressureEquation = false(nT, 1);

row.Si_amp_13cat = repmat(amp.Si, nT, 1);
row.Ti_amp_13cat = repmat(amp.Ti, nT, 1);
row.Al_amp_13cat = repmat(amp.Al, nT, 1);
row.FeT_amp_13cat = repmat(amp.Fe, nT, 1);
row.Mn_amp_13cat = repmat(amp.Mn, nT, 1);
row.Mg_amp_13cat = repmat(amp.Mg, nT, 1);
row.Ca_amp_13cat = repmat(amp.Ca, nT, 1);
row.Na_amp_13cat = repmat(amp.Na, nT, 1);
row.K_amp_13cat = repmat(amp.K, nT, 1);
row.Cr_amp_13cat = repmat(amp.Cr, nT, 1);
row.normFactor_13cat = repmat(amp.normFactor, nT, 1);
row.denominator_13cat = repmat(amp.denominator13, nT, 1);
row.Ca_B_amp_13cat = repmat(amp.Ca_B, nT, 1);
row.Na_B_amp_13cat = repmat(amp.Na_B, nT, 1);
row.CaNa_B_amp_13cat = repmat(amp.CaNa_B, nT, 1);
row.MgNumber_totalFe_proxy = repmat(amp.MgNumber_totalFe, nT, 1);

row.SiO2_amp = repmat(amp.raw.SiO2, nT, 1);
row.TiO2_amp = repmat(amp.raw.TiO2, nT, 1);
row.Al2O3_amp = repmat(amp.raw.Al2O3, nT, 1);
row.FeOt_amp = repmat(amp.raw.FeOt, nT, 1);
row.MnO_amp = repmat(amp.raw.MnO, nT, 1);
row.MgO_amp = repmat(amp.raw.MgO, nT, 1);
row.CaO_amp = repmat(amp.raw.CaO, nT, 1);
row.Na2O_amp = repmat(amp.raw.Na2O, nT, 1);
row.K2O_amp = repmat(amp.raw.K2O, nT, 1);
row.Cr2O3_amp = repmat(amp.raw.Cr2O3, nT, 1);
row.Fe_input_source = repmat(string(amp.feSource), nT, 1);

row.lnP1a = lnP1a;
row.lnP1b = lnP1b;
row.lnP1e = lnP1e;

row.P1a_MPa = P1a;
row.P1b_MPa = P1b;
row.P1c_MPa = P1c;
row.P1d_MPa = P1d;
row.P1e_MPa = P1e;

row.P1a_kbar = P1a ./ 100;
row.P1b_kbar = P1b ./ 100;
row.P1c_kbar = P1c ./ 100;
row.P1d_kbar = P1d ./ 100;
row.P1e_kbar = P1e ./ 100;

row.P1a_GPa = P1a ./ 1000;
row.P1b_GPa = P1b ./ 1000;
row.P1c_GPa = P1c ./ 1000;
row.P1d_GPa = P1d ./ 1000;
row.P1e_GPa = P1e ./ 1000;

row.DPdb_MPa = DPdb;
row.XPae = XPae;
row.APE_percent = APE;
row.P2_MPa = P2;
row.P2_source = P2_source;

row.P_kbar = Pfinal ./ 100;
row.P_RR2012_MPa = Pfinal;
row.P_RR2012_kbar = Pfinal ./ 100;
row.P_RR2012_GPa = Pfinal ./ 1000;
row.finalPressureSource = finalPressureSource;

row.P_RR2012_uncertainty_percent = repmat(11.5, nT, 1);
row.P_RR2012_1sigma_MPa = abs(Pfinal) .* 0.115;
row.P_RR2012_APE50_expectedMaxUncertainty_percent = ...
    repmat(20.0, nT, 1);

row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinCalibrationPRange = isWithinCalibrationPRange;
row.isWithinApproxExperimentalPTField = ...
    isWithinApproxExperimentalPTField;
row.isCalcicAmphibole_numeric = ...
    repmat(amp.isCalcic_numeric, nT, 1);
row.isMgRich_totalFe_conservative_proxy = ...
    repmat(isfinite(amp.MgNumber_totalFe) && ...
    amp.MgNumber_totalFe >= 0.5, nT, 1);
row.isAPEbelow50 = isfinite(APE) & APE < 50;
row.isApplicable_numeric = isApplicable_numeric;

% Backward-compatible general flag.
row.isApplicable = isApplicable_numeric;

end

function printEquationRangeWarning( ...
        pressure_MPa, equationLabel, rangeMin_MPa, rangeMax_MPa, ...
        selectedCode_amp)
% printEquationRangeWarning
% Print a non-stopping warning when an equation-specific finite result lies
% outside the pressure interval used to calibrate that equation.

if isfinite(pressure_MPa) && ...
        (pressure_MPa < rangeMin_MPa || pressure_MPa > rangeMax_MPa)
    fprintf(2, ['WARNING: %s gives %.6g MPa for %s, outside its ' ...
        'calibration interval %.6g-%.6g MPa in Ridolfi and Renzulli ' ...
        '(2012; Table 3, p. 887). The authors still evaluate all five ' ...
        'equations as part of the final selection procedure; do not ' ...
        'interpret this equation independently.\n'], ...
        equationLabel, pressure_MPa, char(string(selectedCode_amp)), ...
        rangeMin_MPa, rangeMax_MPa);
end

end

function value = getMineralOxideRequired( ...
        data_tbl, candidates, displayName)
% getMineralOxideRequired
% Retrieve a required one-row oxide value. Missing columns cause an error;
% explicit NaN is retained. Inf and finite negative values are rejected.

columnName = findOxideColumn( ...
    data_tbl.Properties.VariableNames, candidates);
if isempty(columnName)
    error('Selected Amphibole row must contain variable: %s', displayName);
end

value = convertToScalarDouble(data_tbl.(columnName));
validateNonNegativeScalar(value, displayName);

end

function [value, isPresent] = getMineralOxideOptional( ...
        data_tbl, candidates, displayName, missingDefault)
% getMineralOxideOptional
% Retrieve an optional oxide. A missing column uses the supplied omitted-
% component default. Explicit NaN is retained and never replaced by zero.

columnName = findOxideColumn( ...
    data_tbl.Properties.VariableNames, candidates);
if isempty(columnName)
    value = missingDefault;
    isPresent = false;
    return;
end

value = convertToScalarDouble(data_tbl.(columnName));
validateNonNegativeScalar(value, displayName);
isPresent = true;

end

function columnName = findOxideColumn(variableNames, candidates)
% findOxideColumn
% Find a table column using case-insensitive canonical name matching after
% removing spaces, underscores, hyphens, and periods.

canonicalVariables = strings(numel(variableNames), 1);
for i = 1:numel(variableNames)
    canonicalVariables(i) = canonicalizeName(variableNames{i});
end

columnName = '';
for i = 1:numel(candidates)
    canonicalCandidate = canonicalizeName(candidates{i});
    possibleNames = [canonicalCandidate; canonicalCandidate + "value"];

    for j = 1:numel(possibleNames)
        matchIndex = find( ...
            canonicalVariables == possibleNames(j), 1, 'first');
        if ~isempty(matchIndex)
            columnName = variableNames{matchIndex};
            return;
        end
    end
end

end

function canonicalName = canonicalizeName(inputName)
% canonicalizeName
% Convert a variable name to a lower-case comparison key.

canonicalName = lower(string(inputName));
canonicalName = replace(canonicalName, " ", "");
canonicalName = replace(canonicalName, "_", "");
canonicalName = replace(canonicalName, "-", "");
canonicalName = replace(canonicalName, ".", "");

end

function value = convertToScalarDouble(rawValue)
% convertToScalarDouble
% Convert one scalar numeric, logical, string, char, or cell value to double.
% Missing or non-convertible scalar content is returned as NaN.

if isempty(rawValue)
    value = NaN;
    return;
end
if numel(rawValue) ~= 1
    error('Selected table variable must contain one scalar value.');
end

if isnumeric(rawValue)
    value = double(rawValue);
    return;
end
if islogical(rawValue)
    value = double(rawValue);
    return;
end
if isstring(rawValue)
    if ismissing(rawValue)
        value = NaN;
    else
        value = str2double(rawValue);
    end
    return;
end
if ischar(rawValue)
    value = str2double(string(rawValue));
    return;
end
if iscell(rawValue)
    cellValue = rawValue{1};
    if isempty(cellValue)
        value = NaN;
    elseif isnumeric(cellValue) || islogical(cellValue)
        if numel(cellValue) ~= 1
            error('Selected cell value must contain one scalar value.');
        end
        value = double(cellValue);
    else
        value = str2double(string(cellValue));
    end
    return;
end

value = NaN;

end

function validateNonNegativeScalar(value, displayName)
% validateNonNegativeScalar
% Allow zero and NaN; reject Inf and finite negative values.

if isinf(value) || (isfinite(value) && value < 0)
    error('%s must be non-negative or NaN; Inf is prohibited.', displayName);
end

end
