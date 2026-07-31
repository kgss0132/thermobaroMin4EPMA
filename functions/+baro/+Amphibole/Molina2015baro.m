function results = Molina2015baro(rawdata_struct, T_degreeC)
% functions/+baro/+Amphibole/Molina2015baro.m
% Tested with MATLAB R2024b
%
% Plagioclase-Amphibole Al-Si partitioning barometer
% Molina, J.F., Moreno, J.A., Castro, A., Rodriguez, C., and Fershtater, G.B. (2015)
% Lithos, 232, 286-305
% DOI: https://doi.org/10.1016/j.lithos.2015.06.027
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Amphibole analysis and one
% Plagioclase analysis and calculates pressure using the Al-Si partitioning
% barometer of Molina et al. (2015).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Amphibole-Plagioclase pair, one
% output row is returned for every input temperature value.
%
% Temperature enters the published pressure equation directly. The selected
% Amphibole and Plagioclase analyses must therefore represent an equilibrium
% pair, and the supplied temperature must represent the same equilibration
% event.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Molina et al. (2015) calibrated the barometer using 92 selected
% Amphibole-Plagioclase experimental pairs. The equation is presented in
% Section 6.1 on p. 295 and is summarized in the abstract on p. 286.
%
% Experimental-data ranges for the selected barometer dataset (Table 2,
% p. 289) are approximately:
%
%   Temperature : 650-1050 degreeC
%   Pressure    : 1-15 kbar
%   XAb_plag    : 0.108-0.840
%   Si_amp      : 5.791-7.062 apfu
%   XAl_T1_amp  : 0.234-0.552
%   RTlnD       : -2066 to 10975 J mol^-1
%
% The 1-15 kbar interval is the pressure envelope represented by the selected
% experimental dataset, rather than a separately stated universal validity
% range. Values outside this interval are extrapolations and trigger a
% non-stopping warning in this implementation.
%
% Required Amphibole conditions (13-CNK normalization, 23 oxygens):
%
%   Temperature : 650-1050 degreeC
%   Ti_amp      : > 0.02 apfu
%   AlVI_amp    : > 0.05 apfu
%
% These restrictions are stated in Section 6.1 on pp. 294-295, in the
% applications discussion on p. 298, and in the Conclusions on p. 304.
%
% IMPORTANT NORMALIZATION REQUIREMENT:
% Amphibole formulae and Fe3+/Fe2+ ratios in the calibration were calculated
% on a 23-oxygen basis using the 13-CNK method: 13 cations exclusive of Ca,
% Na, and K (Section 5.1, p. 292). This function reads previously normalized
% cation data and cannot verify that the upstream calculation used the
% 13-CNK method. AlVI_amp_est is only a simple estimate from total Al and Si;
% it is not a replacement for a complete 13-CNK site allocation.
%
% Mineral-pair equilibrium:
% The barometer is based on Al-Si partitioning between coexisting phases.
% Amphibole and Plagioclase must therefore be texturally and compositionally
% equilibrated and must record the same P-T event. The experimental pairs
% used for regression were screened using the Holland and Blundy (1994)
% Amphibole-Plagioclase thermometers (data-selection discussion, pp. 287-289).
% Non-adjacent analyses, mismatched core-rim zones, inherited crystals,
% alteration, or later subsolidus re-equilibration may give misleading
% pressures.
%
% Compositional compatibility:
% The authors recommend that natural phase compositions lie within or close
% to the experimental compositional fields shown in Fig. 14 (pp. 301-302).
% Individual one-dimensional range checks cannot fully reproduce those
% multivariate fields. This function therefore reports several numerical
% range flags, but the user must evaluate the full compositional and
% petrographic compatibility independently.
%
% Rock types and assemblages:
% The barometer is applicable to a wide range of calcic-amphibole-
% plagioclase-bearing metamorphic and igneous rocks, including amphibolites,
% mafic granulites, metaluminous granitoids, and gabbros. Quartz and garnet
% are not required (abstract, p. 286; Conclusions, p. 304).
%
% Precision:
%   Calibration dataset : approximately +/-1.5 kbar at 1 sigma
%   Independent test set: approximately +/-2.3 kbar at 1 sigma
%
% These values are reported in Section 6.1 and Table 5 on pp. 295-297.
% Additional uncertainty may arise from temperature uncertainty, analytical
% error, zoning, and incorrect mineral pairing. The pressure sensitivity to
% temperature is approximately 0.032 kbar K^-1 (Clapeyron slope, p. 295).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 650-1050 degreeC,
%   2) finite calculated pressure is outside the 1-15 kbar experimental-data
%      envelope,
%   3) Ti_amp is not greater than 0.02 apfu,
%   4) estimated AlVI_amp is not greater than 0.05 apfu,
%   5) XAb, Si_amp, XAl_T1_amp, or RTlnD lies outside the corresponding
%      one-dimensional experimental-data envelope,
%   6) the equation domain is invalid because a required ratio or logarithm
%      cannot be calculated,
%   7) a calculation or screening input contains NaN,
%   8) calculated pressure is NaN or Inf, or
%   9) negative finite pressure is calculated.
%
% All warnings are diagnostic and do not stop the calculation. Full
% applicability cannot be established automatically because this function
% cannot verify 13-CNK normalization, mineral-pair equilibrium, zoning,
% alteration, or compatibility with the multivariate fields in Fig. 14.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Amphibole   : table
%   rawdata_struct.Plagioclase : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include:
%
%   Required Amphibole variables:
%     Si_cation_apfu           % used in XAl_T1 and D_AlSi
%     Ti_cation_apfu           % application screening
%     Al_cation_apfu           % used in D_AlSi and estimated AlVI
%     Fe_cation_apfu           % total Fe, retained in output
%     Mg_cation_apfu
%     Ca_cation_apfu
%     Na_cation_apfu
%
%   Optional Amphibole variables retained in the output when present:
%     Fe3_cation_apfu
%     Mn_cation_apfu
%     K_cation_apfu
%     Cr_cation_apfu
%
%   Required Plagioclase variables:
%     Si_cation_apfu           % used in D_AlSi
%     Al_cation_apfu           % used in D_AlSi
%     Ca_cation_apfu           % used in feldspar normalization
%     Na_cation_apfu           % used in XAb
%     K_cation_apfu            % used in feldspar normalization
%
% Finite input cation values must be non-negative. NaN is allowed, retained
% as missing, propagated through the calculation, and reported by fprintf.
% Inf and finite negative values are rejected. Zero is retained; if zero
% makes a denominator or logarithm invalid, the derived value and pressure
% are returned as NaN and a warning is printed.
%
% This barometer does not use a liquid (Liq) composition. Therefore, rules
% concerning exclusion of F and Cl from cationTotal_liq and their exclusion
% from Liq NaN warnings are not applicable to this function.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (as implemented here)
%
%   P(kbar) =
%     [R*T(K)*ln(D_AlSi_plg_amp) - 8.7*T(K)
%      + 23377*XAl_T1_amp + 7579*XAb_plag - 11302] / -274
%
% where:
%
%   R = 8.3144 J K^-1 mol^-1
%
%   XAl_T1_amp = (8 - Si_amp) / 4
%
%   XAb_plag = Na_plag / (Ca_plag + Na_plag + K_plag)
%   XAn_plag = Ca_plag / (Ca_plag + Na_plag + K_plag)
%   XOr_plag = K_plag  / (Ca_plag + Na_plag + K_plag)
%
%   D_AlSi_plg_amp =
%     (Al_plag / Si_plag) / (Al_amp / Si_amp)
%
% The regression coefficients and equation are given on p. 295. NaN inputs
% are not replaced by zero. Negative calculated pressure is retained for
% diagnostic purposes and reported by a non-stopping warning.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Molina2015baro(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing Amphibole and Plagioclase tables
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Amphibole-Plagioclase pair
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('Molina2015baro requires (rawdata_struct, T_degreeC).');
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
% Extract the required tables. The source tables are not modified.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Amphibole') || ~istable(rawdata_struct.Amphibole)
    error('rawdata_struct must contain table: rawdata_struct.Amphibole');
end
if ~isfield(rawdata_struct, 'Plagioclase') || ~istable(rawdata_struct.Plagioclase)
    error('rawdata_struct must contain table: rawdata_struct.Plagioclase');
end

dataset_amp = rawdata_struct.Amphibole;
dataset_plag = rawdata_struct.Plagioclase;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Experimental-data limits and formal compositional restrictions.
calibrationT_min_degreeC = 650;
calibrationT_max_degreeC = 1050;
experimentalP_min_kbar = 1;
experimentalP_max_kbar = 15;
minimumTi_amp_apfu = 0.02;
minimumAlVI_amp_apfu = 0.05;
experimentalXAb_min = 0.1083147;
experimentalXAb_max = 0.84;
experimentalSi_amp_min = 5.791408;
experimentalSi_amp_max = 7.062472;
experimentalXAlT1_min = 0.234382;
experimentalXAlT1_max = 0.5521481;
experimentalRTlnD_min_J = -2065.716;
experimentalRTlnD_max_J = 10974.73;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);

temperatureWarningIssued = false;
applicationCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Amphibole) ===');

while true
    % ----- Amphibole selection -----
    % The first table column is used only as the displayed identifier.
    dataCodes_amp = dataset_amp{:, 1};

    [selectedIdx_amp, ok] = listdlg( ...
        'PromptString', 'Please select the Amphibole data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_amp)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_amp)
        disp('Selection canceled');
        break;
    end

    selectedCode_amp = dataCodes_amp(selectedIdx_amp);
    disp(['Amphibole selected: ' char(string(selectedCode_amp))]);

    % ----- Plagioclase selection -----
    disp('=== Step 4: Selecting a data code from the list (Plagioclase) ===');

    dataCodes_plag = dataset_plag{:, 1};

    [selectedIdx_plag, ok] = listdlg( ...
        'PromptString', 'Please select the Plagioclase data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_plag)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_plag)
        disp('Selection canceled');
        break;
    end

    selectedCode_plag = dataCodes_plag(selectedIdx_plag);
    disp(['Plagioclase selected: ' char(string(selectedCode_plag))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the pressure ===');

    selectedData_amp = dataset_amp(selectedIdx_amp, :);
    selectedData_plag = dataset_plag(selectedIdx_plag, :);

    % Required columns must exist. Inf and finite negative values are
    % rejected; zero and NaN are retained.
    validateNonNegativeInputs(selectedData_amp, selectedData_plag);

    % Report only NaN variables used by the equation or its numerical
    % applicability screening. NaN values are never changed to zero.
    nanInputNames = findNaNInputs( ...
        selectedData_amp, selectedData_plag, T_degreeC);

    row = calcPressure(selectedData_amp, selectedData_plag, T_degreeC);

    % Repeat identifiers for all temperature rows in this calculation.
    row.dataCode_amphibole = ...
        repmat(string(selectedCode_amp), height(row), 1);
    row.dataCode_plagioclase = ...
        repmat(string(selectedCode_plag), height(row), 1);
    row = movevars(row, ...
        {'dataCode_amphibole','dataCode_plagioclase'}, 'Before', 1);

    % Store one block per selected pair. Expand the cell buffer only when
    % its current capacity has been exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    finitePressure = isfinite(row.P_kbar);
    finiteTemperature = isfinite(row.T_degreeC);

    if height(row) == 1
        fprintf('%s & %s: P = %.6g kbar at T = %.6g degreeC\n', ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            row.P_kbar, row.T_degreeC);
    elseif any(finitePressure)
        finitePressureValues = row.P_kbar(finitePressure);
        if any(finiteTemperature)
            finiteTemperatureValues = row.T_degreeC(finiteTemperature);
            fprintf('%s & %s: P = %.6g to %.6g kbar for finite T = %.6g to %.6g degreeC\n', ...
                char(string(selectedCode_amp)), ...
                char(string(selectedCode_plag)), ...
                min(finitePressureValues), max(finitePressureValues), ...
                min(finiteTemperatureValues), max(finiteTemperatureValues));
        else
            fprintf('%s & %s: P = %.6g to %.6g kbar; all %d input temperatures are NaN\n', ...
                char(string(selectedCode_amp)), ...
                char(string(selectedCode_plag)), ...
                min(finitePressureValues), max(finitePressureValues), ...
                height(row));
        end
    elseif any(finiteTemperature)
        finiteTemperatureValues = row.T_degreeC(finiteTemperature);
        fprintf('%s & %s: P = NaN/Inf for all %d temperature point(s), finite T = %.6g to %.6g degreeC\n', ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            height(row), ...
            min(finiteTemperatureValues), max(finiteTemperatureValues));
    else
        fprintf('%s & %s: P = NaN/Inf for all %d point(s); all input temperatures are NaN\n', ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            height(row));
    end

    % Print major limitations that cannot be assessed numerically once per
    % function call.
    if ~applicationCautionIssued
        fprintf(2, ['CAUTION: Molina et al. (2015) require Amphibole formulae ' ...
            'normalized by the 13-CNK method and an equilibrium Amphibole-' ...
            'Plagioclase pair recording the same P-T event (pp. 287-289, 292, ' ...
            '298, 304). This function cannot verify normalization method, ' ...
            'textural equilibrium, zoning, alteration, or compatibility with ' ...
            'the multivariate experimental fields in Fig. 14 (pp. 301-302).\n']);
        applicationCautionIssued = true;
    end

    % Temperature is common to all selected pairs during this function call,
    % so the range warning is printed only once.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperatureValues = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ['WARNING: Input temperature is outside the Molina et al. ' ...
            '(2015) barometer range of 650-1050 degreeC (pp. 286-289, ' ...
            '294-295, 304). %d of %d finite temperature point(s) are ' ...
            'outside the range; input finite range = %.6g-%.6g degreeC.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperatureValues), ...
            min(finiteTemperatureValues), ...
            max(finiteTemperatureValues));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures lie outside the pressure envelope
    % represented by the selected experimental dataset.
    pressureOutsideExperimentalRange = finitePressure & ...
        (row.P_kbar < experimentalP_min_kbar | ...
         row.P_kbar > experimentalP_max_kbar);

    if any(pressureOutsideExperimentalRange)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ['WARNING: Calculated pressure is outside the approximately ' ...
            '1-15 kbar experimental-data envelope of the 92 selected ' ...
            'Amphibole-Plagioclase pairs in Molina et al. (2015; Table 2, ' ...
            'p. 289). %d of %d finite pressure point(s) are outside the ' ...
            'envelope; calculated finite range = %.6g-%.6g kbar for %s & %s. ' ...
            'The values are retained in the output table.\n'], ...
            sum(pressureOutsideExperimentalRange), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    % Formal minimum Amphibole composition restrictions.
    finiteTi = isfinite(row.Ti_amp(1));
    if finiteTi && row.Ti_amp(1) <= minimumTi_amp_apfu
        fprintf(2, ['WARNING: Amphibole Ti = %.6g apfu is not greater than ' ...
            'the required minimum of 0.02 apfu for Molina et al. (2015; ' ...
            'pp. 294-295, 298, 304) for %s & %s.\n'], ...
            row.Ti_amp(1), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    finiteAlVI = isfinite(row.AlVI_amp_est(1));
    if finiteAlVI && row.AlVI_amp_est(1) <= minimumAlVI_amp_apfu
        fprintf(2, ['WARNING: Estimated Amphibole AlVI = %.6g apfu is not ' ...
            'greater than the required minimum of 0.05 apfu for Molina et al. ' ...
            '(2015; pp. 294-295, 298, 304) for %s & %s. The value used here ' ...
            'is a simple estimate and does not replace 13-CNK site allocation.\n'], ...
            row.AlVI_amp_est(1), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    % One-dimensional experimental composition envelopes. These checks are
    % descriptive and do not reproduce the multivariate fields in Fig. 14.
    finiteXAb = isfinite(row.XAb_plag(1));
    if finiteXAb && ...
            (row.XAb_plag(1) < experimentalXAb_min || ...
             row.XAb_plag(1) > experimentalXAb_max)
        fprintf(2, ['WARNING: Plagioclase XAb = %.6g is outside the ' ...
            'experimental range 0.108-0.840 in Table 2 of Molina et al. ' ...
            '(2015; p. 289) for %s & %s. This is a descriptive data ' ...
            'envelope, not an absolute compositional cutoff.\n'], ...
            row.XAb_plag(1), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    finiteSiAmp = isfinite(row.Si_amp(1));
    if finiteSiAmp && ...
            (row.Si_amp(1) < experimentalSi_amp_min || ...
             row.Si_amp(1) > experimentalSi_amp_max)
        fprintf(2, ['WARNING: Amphibole Si = %.6g apfu is outside the ' ...
            'experimental range 5.791-7.062 apfu in Table 2 of Molina et al. ' ...
            '(2015; p. 289) for %s & %s.\n'], ...
            row.Si_amp(1), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    finiteXAlT1 = isfinite(row.XAl_T1_amp(1));
    if finiteXAlT1 && ...
            (row.XAl_T1_amp(1) < experimentalXAlT1_min || ...
             row.XAl_T1_amp(1) > experimentalXAlT1_max)
        fprintf(2, ['WARNING: Amphibole XAl_T1 = %.6g is outside the ' ...
            'experimental range 0.234-0.552 in Table 2 of Molina et al. ' ...
            '(2015; p. 289) for %s & %s.\n'], ...
            row.XAl_T1_amp(1), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    finiteRTlnD = isfinite(row.RTlnD_AlSi_plg_amp_J);
    rtlnDOutsideExperimentalRange = finiteRTlnD & ...
        (row.RTlnD_AlSi_plg_amp_J < experimentalRTlnD_min_J | ...
         row.RTlnD_AlSi_plg_amp_J > experimentalRTlnD_max_J);
    if any(rtlnDOutsideExperimentalRange)
        finiteRTlnDValues = row.RTlnD_AlSi_plg_amp_J(finiteRTlnD);
        fprintf(2, ['WARNING: RTlnD_AlSi is outside the approximate ' ...
            'experimental range -2066 to 10975 J mol^-1 in Table 2 of ' ...
            'Molina et al. (2015; p. 289). %d of %d finite point(s) are ' ...
            'outside the range; finite RTlnD range = %.6g-%.6g J mol^-1 ' ...
            'for %s & %s.\n'], ...
            sum(rtlnDOutsideExperimentalRange), ...
            sum(finiteRTlnD), ...
            min(finiteRTlnDValues), ...
            max(finiteRTlnDValues), ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)));
    end

    % Report an invalid mathematical domain without stopping. Zero and NaN
    % inputs are retained, and the corresponding result remains NaN.
    if ~row.isWithinEquationDomain(1)
        fprintf(2, ['WARNING: The Molina et al. (2015) equation domain is ' ...
            'invalid for %s & %s. The calculation requires positive finite ' ...
            'Al/Si ratios in both minerals, Ca+Na+K > 0 in Plagioclase, and ' ...
            'a positive finite D_AlSi_plg_amp. Si_amp = %.6g, Al_amp = %.6g, ' ...
            'Si_plag = %.6g, Al_plag = %.6g, Ca+Na+K_plag = %.6g, ' ...
            'D_AlSi_plg_amp = %.6g. Pressure is retained as NaN.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            row.Si_amp(1), row.Al_amp(1), ...
            row.Si_plag(1), row.Al_plag(1), ...
            row.sumCNK_plag(1), row.D_AlSi_plg_amp(1));
    end

    % List exact calculation/screening inputs that contain NaN.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the barometer input(s) for ' ...
            '%s & %s: %s.\n         NaN values were retained and were not ' ...
            'replaced by zero; derived values and pressure may remain NaN.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite pressure results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ['WARNING: Non-finite pressure values were calculated for ' ...
            '%s & %s (%d of %d points; NaN: %d, Inf: %d).\n         These ' ...
            'values remain in the output table, and the calculation has not ' ...
            'been stopped.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative finite pressure is retained for diagnostic purposes.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ['WARNING: Negative finite pressure was calculated for ' ...
            '%s & %s (%d of %d points). The values were retained for ' ...
            'diagnostic purposes and are outside the experimental pressure ' ...
            'envelope.\n'], ...
            char(string(selectedCode_amp)), ...
            char(string(selectedCode_plag)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    % Ask whether to repeat with another mineral pair.
    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Molina2015baro', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered result blocks once. Return an empty table when
% the user canceled before completing any calculation.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs( ...
        data_amphibole, data_plagioclase, T_degreeC)
% findNaNInputs
% Return names of NaN inputs used by the pressure equation or numerical
% applicability screening. NaN values do not cause an error and are never
% replaced by zero.

maxNames = 9;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = ...
        "T_degreeC(indices=" + indexText + ")";
end

amphiboleVariables = { ...
    'Si_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Al_cation_apfu'};

for i = 1:numel(amphiboleVariables)
    variableName = amphiboleVariables{i};
    variableValue = data_amphibole.(variableName);
    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = ...
            "Amphibole." + string(variableName);
    end
end

plagioclaseVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu'};

for i = 1:numel(plagioclaseVariables)
    variableName = plagioclaseVariables{i};
    variableValue = data_plagioclase.(variableName);
    if any(isnan(variableValue(:)))
        nNanInputs = nNanInputs + 1;
        nanInputBuffer(nNanInputs) = ...
            "Plagioclase." + string(variableName);
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_amphibole, data_plagioclase)
% validateNonNegativeInputs
% Require all calculation/output cation variables to be scalar and reject
% Inf and finite negative values. Zero and NaN are allowed and retained.

requiredAmphiboleVariables = { ...
    'Si_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu'};

optionalAmphiboleVariables = { ...
    'Fe3_cation_apfu', ...
    'Mn_cation_apfu', ...
    'K_cation_apfu', ...
    'Cr_cation_apfu'};

requiredPlagioclaseVariables = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu'};

for i = 1:numel(requiredAmphiboleVariables)
    variableName = requiredAmphiboleVariables{i};
    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        error('Amphibole table must contain variable: %s', variableName);
    end
end

for i = 1:numel(requiredPlagioclaseVariables)
    variableName = requiredPlagioclaseVariables{i};
    if ~ismember(variableName, data_plagioclase.Properties.VariableNames)
        error('Plagioclase table must contain variable: %s', variableName);
    end
end

maxInvalidNames = numel(requiredAmphiboleVariables) + ...
    numel(optionalAmphiboleVariables) + ...
    numel(requiredPlagioclaseVariables);
invalidInputBuffer = strings(maxInvalidNames, 1);
nInvalidInputs = 0;

allAmphiboleVariables = ...
    [requiredAmphiboleVariables, optionalAmphiboleVariables];

for i = 1:numel(allAmphiboleVariables)
    variableName = allAmphiboleVariables{i};
    if ~ismember(variableName, data_amphibole.Properties.VariableNames)
        continue;
    end

    variableValue = data_amphibole.(variableName);
    if ~isscalar(variableValue)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end

    if isinf(variableValue) || ...
            (isfinite(variableValue) && variableValue < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Amphibole." + string(variableName);
    end
end

for i = 1:numel(requiredPlagioclaseVariables)
    variableName = requiredPlagioclaseVariables{i};
    variableValue = data_plagioclase.(variableName);

    if ~isscalar(variableValue)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end

    if isinf(variableValue) || ...
            (isfinite(variableValue) && variableValue < 0)
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = ...
            "Plagioclase." + string(variableName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['Molina2015baro: cation values must be non-negative. NaN is ' ...
           'allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure( ...
        data_amphibole, data_plagioclase, T_degreeC)
% calcPressure
% Compute pressure for one Amphibole row and one Plagioclase row at one or
% more input temperatures. NaN values propagate through the calculation.
%
% Inputs:
%   data_amphibole   : 1-row Amphibole table
%   data_plagioclase : 1-row Plagioclase table
%   T_degreeC        : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end
if height(data_plagioclase) ~= 1
    error('Plagioclase input must be a 1-row table.');
end

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;
R_J = 8.3144;

amp = prepareAmphiboleRow(data_amphibole);
plag = preparePlagioclaseRow(data_plagioclase);

% Simple Al-site estimates retained for screening and output. NaN values
% propagate. These values do not constitute a complete 13-CNK site allocation.
if isnan(amp.Si) || isnan(amp.Al)
    AlIV_amp_est_scalar = NaN;
    AlVI_amp_est_scalar = NaN;
else
    AlIV_amp_est_scalar = min(max(8 - amp.Si, 0), amp.Al);
    AlVI_amp_est_scalar = amp.Al - AlIV_amp_est_scalar;
end

% Amphibole T1-site Al fraction used directly by the published equation.
XAl_T1_amp_scalar = (8 - amp.Si) ./ 4;

% Plagioclase component fractions. Zero or NaN total produces NaN fractions
% rather than a stopping error.
sumCNK_plag_scalar = plag.Ca + plag.Na + plag.K;
if isfinite(sumCNK_plag_scalar) && sumCNK_plag_scalar > 0
    XAb_plag_scalar = plag.Na ./ sumCNK_plag_scalar;
    XAn_plag_scalar = plag.Ca ./ sumCNK_plag_scalar;
    XOr_plag_scalar = plag.K ./ sumCNK_plag_scalar;
else
    XAb_plag_scalar = NaN;
    XAn_plag_scalar = NaN;
    XOr_plag_scalar = NaN;
end

% Al/Si ratios and partition coefficient. Raw ratios are retained for
% diagnostics. lnD is set to NaN unless D is finite and strictly positive,
% preventing complex or infinite logarithmic results.
AlSi_plag_scalar = plag.Al ./ plag.Si;
AlSi_amp_scalar = amp.Al ./ amp.Si;
D_AlSi_plg_amp_scalar = AlSi_plag_scalar ./ AlSi_amp_scalar;

if isfinite(D_AlSi_plg_amp_scalar) && D_AlSi_plg_amp_scalar > 0
    lnD_AlSi_plg_amp_scalar = log(D_AlSi_plg_amp_scalar);
else
    lnD_AlSi_plg_amp_scalar = NaN;
end

% Temperature-dependent terms and pressure calculation. NaN temperature or
% composition inputs propagate naturally to pressure.
RTlnD_AlSi_plg_amp_J = ...
    R_J .* T_K .* lnD_AlSi_plg_amp_scalar;

P_kbar = (RTlnD_AlSi_plg_amp_J - 8.7 .* T_K + ...
    23377 .* XAl_T1_amp_scalar + ...
    7579 .* XAb_plag_scalar - 11302) ./ -274;

% Mathematical domain and application flags.
isWithinEquationDomain_scalar = ...
    isfinite(amp.Si) && amp.Si > 0 && ...
    isfinite(amp.Al) && amp.Al > 0 && ...
    isfinite(plag.Si) && plag.Si > 0 && ...
    isfinite(plag.Al) && plag.Al > 0 && ...
    isfinite(sumCNK_plag_scalar) && sumCNK_plag_scalar > 0 && ...
    isfinite(D_AlSi_plg_amp_scalar) && D_AlSi_plg_amp_scalar > 0;

isWithinCalibrationTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 650 & T_degreeC <= 1050;

isWithinExperimentalPRange = ...
    isfinite(P_kbar) & P_kbar >= 1 & P_kbar <= 15;

isTiAboveMinimum_scalar = ...
    isfinite(amp.Ti) && amp.Ti > 0.02;

isAlVIAboveMinimum_scalar = ...
    isfinite(AlVI_amp_est_scalar) && AlVI_amp_est_scalar > 0.05;

isWithinExperimentalXAbRange_scalar = ...
    isfinite(XAb_plag_scalar) && ...
    XAb_plag_scalar >= 0.1083147 && XAb_plag_scalar <= 0.84;

isWithinExperimentalSiRange_scalar = ...
    isfinite(amp.Si) && amp.Si >= 5.791408 && amp.Si <= 7.062472;

isWithinExperimentalXAlT1Range_scalar = ...
    isfinite(XAl_T1_amp_scalar) && ...
    XAl_T1_amp_scalar >= 0.234382 && XAl_T1_amp_scalar <= 0.5521481;

isWithinExperimentalRTlnDRange = ...
    isfinite(RTlnD_AlSi_plg_amp_J) & ...
    RTlnD_AlSi_plg_amp_J >= -2065.716 & ...
    RTlnD_AlSi_plg_amp_J <= 10974.73;

% This flag evaluates only numerical restrictions available to the function.
% It does not verify 13-CNK normalization or mineral-pair equilibrium.
isApplicable_numeric = ...
    isWithinCalibrationTRange & ...
    isWithinExperimentalPRange & ...
    repmat(isWithinEquationDomain_scalar, nT, 1) & ...
    repmat(isTiAboveMinimum_scalar, nT, 1) & ...
    repmat(isAlVIAboveMinimum_scalar, nT, 1);

% Expand composition-dependent scalars to the temperature-vector length.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.R_J = repmat(R_J, nT, 1);

row.Si_amp = repmat(amp.Si, nT, 1);
row.Ti_amp = repmat(amp.Ti, nT, 1);
row.Al_amp = repmat(amp.Al, nT, 1);
row.AlIV_amp_est = repmat(AlIV_amp_est_scalar, nT, 1);
row.AlVI_amp_est = repmat(AlVI_amp_est_scalar, nT, 1);
row.FeT_amp = repmat(amp.FeT, nT, 1);
row.Fe2_amp = repmat(amp.Fe2, nT, 1);
row.Fe3_amp = repmat(amp.Fe3, nT, 1);
row.Mg_amp = repmat(amp.Mg, nT, 1);
row.Ca_amp = repmat(amp.Ca, nT, 1);
row.Na_amp = repmat(amp.Na, nT, 1);
row.K_amp = repmat(amp.K, nT, 1);
row.Mn_amp = repmat(amp.Mn, nT, 1);
row.Cr_amp = repmat(amp.Cr, nT, 1);

row.Si_plag = repmat(plag.Si, nT, 1);
row.Al_plag = repmat(plag.Al, nT, 1);
row.Ca_plag = repmat(plag.Ca, nT, 1);
row.Na_plag = repmat(plag.Na, nT, 1);
row.K_plag = repmat(plag.K, nT, 1);
row.sumCNK_plag = repmat(sumCNK_plag_scalar, nT, 1);

row.XAb_plag = repmat(XAb_plag_scalar, nT, 1);
row.XAn_plag = repmat(XAn_plag_scalar, nT, 1);
row.XOr_plag = repmat(XOr_plag_scalar, nT, 1);

% Lower-case aliases follow the naming style used by other barometers.
row.Xab_plag = row.XAb_plag;
row.Xan_plag = row.XAn_plag;
row.Xor_plag = row.XOr_plag;

row.XAl_T1_amp = repmat(XAl_T1_amp_scalar, nT, 1);
row.AlSi_plag = repmat(AlSi_plag_scalar, nT, 1);
row.AlSi_amp = repmat(AlSi_amp_scalar, nT, 1);
row.D_AlSi_plg_amp = repmat(D_AlSi_plg_amp_scalar, nT, 1);
row.lnD_AlSi_plg_amp = repmat(lnD_AlSi_plg_amp_scalar, nT, 1);
row.RTlnD_AlSi_plg_amp_J = RTlnD_AlSi_plg_amp_J;

% Primary launcher-compatible pressure name plus the original alias.
row.P_kbar = P_kbar;
row.P_Molina2015_kbar = P_kbar;
row.P_uncertainty_calibration_1sigma_kbar = repmat(1.5, nT, 1);
row.P_uncertainty_test_1sigma_kbar = repmat(2.3, nT, 1);
row.P_temperature_sensitivity_kbar_per_K = repmat(0.032, nT, 1);

row.isWithinEquationDomain = ...
    repmat(isWithinEquationDomain_scalar, nT, 1);
row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinExperimentalPRange = isWithinExperimentalPRange;
row.isTiAboveMinimum = repmat(isTiAboveMinimum_scalar, nT, 1);
row.isAlVIAboveMinimum = repmat(isAlVIAboveMinimum_scalar, nT, 1);
row.isWithinExperimentalXAbRange = ...
    repmat(isWithinExperimentalXAbRange_scalar, nT, 1);
row.isWithinExperimentalSiRange = ...
    repmat(isWithinExperimentalSiRange_scalar, nT, 1);
row.isWithinExperimentalXAlT1Range = ...
    repmat(isWithinExperimentalXAlT1Range_scalar, nT, 1);
row.isWithinExperimentalRTlnDRange = ...
    isWithinExperimentalRTlnDRange;
row.requires13CNKNormalization = true(nT, 1);
row.is13CNKNormalizationVerified = false(nT, 1);
row.isApplicable_numeric = isApplicable_numeric;

% Backward-compatible name from the original implementation. This remains a
% numerical diagnostic and does not verify complete applicability.
row.isApplicable_Molina2015 = isApplicable_numeric;
row.isApplicable = isApplicable_numeric;

end

function amp = prepareAmphiboleRow(data_amphibole)
% prepareAmphiboleRow
% Extract one-row Amphibole cation data. Required variables must exist;
% missing optional variables are represented as NaN, never as zero.

if height(data_amphibole) ~= 1
    error('Amphibole input must be a 1-row table.');
end

amp = struct();
amp.Si = getVarOrError(data_amphibole, ...
    'Si_cation_apfu', 'Amphibole');
amp.Ti = getVarOrError(data_amphibole, ...
    'Ti_cation_apfu', 'Amphibole');
amp.Al = getVarOrError(data_amphibole, ...
    'Al_cation_apfu', 'Amphibole');
amp.FeT = getVarOrError(data_amphibole, ...
    'Fe_cation_apfu', 'Amphibole');
amp.Mg = getVarOrError(data_amphibole, ...
    'Mg_cation_apfu', 'Amphibole');
amp.Ca = getVarOrError(data_amphibole, ...
    'Ca_cation_apfu', 'Amphibole');
amp.Na = getVarOrError(data_amphibole, ...
    'Na_cation_apfu', 'Amphibole');

amp.Fe3 = getVarOrNaN(data_amphibole, 'Fe3_cation_apfu');
amp.Mn = getVarOrNaN(data_amphibole, 'Mn_cation_apfu');
amp.K = getVarOrNaN(data_amphibole, 'K_cation_apfu');
amp.Cr = getVarOrNaN(data_amphibole, 'Cr_cation_apfu');

% Missing Fe3 is not interpreted as zero; derived Fe2 remains NaN.
if isfinite(amp.Fe3) && isfinite(amp.FeT) && amp.Fe3 > amp.FeT
    error('Amphibole contains Fe3_cation_apfu > Fe_cation_apfu.');
end

amp.Fe2 = amp.FeT - amp.Fe3;
if isfinite(amp.Fe2) && amp.Fe2 < 0
    error('Amphibole contains a negative derived Fe2 value.');
end

end

function plag = preparePlagioclaseRow(data_plagioclase)
% preparePlagioclaseRow
% Extract one-row Plagioclase cation data without replacing NaN by zero.

if height(data_plagioclase) ~= 1
    error('Plagioclase input must be a 1-row table.');
end

plag = struct();
plag.Si = getVarOrError(data_plagioclase, ...
    'Si_cation_apfu', 'Plagioclase');
plag.Al = getVarOrError(data_plagioclase, ...
    'Al_cation_apfu', 'Plagioclase');
plag.Ca = getVarOrError(data_plagioclase, ...
    'Ca_cation_apfu', 'Plagioclase');
plag.Na = getVarOrError(data_plagioclase, ...
    'Na_cation_apfu', 'Plagioclase');
plag.K = getVarOrError(data_plagioclase, ...
    'K_cation_apfu', 'Plagioclase');

end

function value = getVarOrError(tbl, variableName, mineralLabel)
% getVarOrError
% Retrieve one required scalar variable without altering NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);
if ~isscalar(value)
    error('Variable %s must be scalar in a 1-row table.', variableName);
end

if isinf(value) || (isfinite(value) && value < 0)
    error('Variable %s must be non-negative or NaN.', variableName);
end

end

function value = getVarOrNaN(tbl, variableName)
% getVarOrNaN
% Retrieve one optional scalar variable. Missing variables and explicit NaN
% values remain NaN and are never interpreted as zero.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);

    if ~isscalar(value)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end

    if isinf(value) || (isfinite(value) && value < 0)
        error('Variable %s must be non-negative or NaN.', variableName);
    end
else
    value = NaN;
end

end
