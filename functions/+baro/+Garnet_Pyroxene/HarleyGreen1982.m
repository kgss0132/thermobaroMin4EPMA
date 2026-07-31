function results = HarleyGreen1982(rawdata_struct, T_degreeC)
% functions/+baro/+Garnet_Pyroxene/HarleyGreen1982.m
% Tested with MATLAB R2024b
%
% Garnet-Orthopyroxene barometer
% Harley, S.L. and Green, D.H. (1982)
% Nature, 300, 697-701
% DOI: https://doi.org/10.1038/300697a0
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Orthopyroxene analysis selected from the input tables and calculates
% pressure using Eq. (5) of Harley and Green (1982, p. 700).
%
% The function accepts either a scalar temperature or a temperature vector.
% It is therefore compatible with both startBaroCalc_fixedT and
% startBaroCalc_rangeT. For each selected Garnet-Orthopyroxene pair, one
% output row is returned for every input temperature value.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL RANGE AND APPLICATION NOTES
%
% Harley and Green (1982) calibrated the Garnet-Orthopyroxene barometer
% using experiments in the FMAS and CFMAS systems. The experimental range
% stated in the abstract is:
%
%   Temperature : 800-1200 degreeC
%   Pressure    : 5-20 kbar
%   Systems     : FeO-MgO-Al2O3-SiO2 and
%                 CaO-FeO-MgO-Al2O3-SiO2
%
% The discussion of experiments aimed directly at crustal granulites gives
% approximately 800-1150 degreeC and 5-20 kbar (p. 697). The broader
% 800-1200 degreeC and 5-20 kbar limits stated in the abstract are used in
% this implementation as non-stopping calibration-envelope warnings.
%
% The equation is intended for equilibrated Garnet-Orthopyroxene
% assemblages in granulite-facies rocks and garnet peridotites (pp. 697,
% 700-701). Important limitations stated by Harley and Green (1982) are:
%
%   1) High Fe3+ in Garnet or Orthopyroxene may require additional
%      non-ideal terms and may involve Fe3+-Al exchange in Orthopyroxene
%      (p. 700).
%   2) High Cr3+ in Garnet limits application of the barometer (p. 700).
%   3) High Mn2+ in the phases limits application (p. 700).
%   4) Very low Al2O3 in Orthopyroxene produces considerable pressure
%      uncertainty, particularly in high-pressure garnet peridotites and
%      low-temperature assemblages containing Fe-Ca-rich Garnet (p. 700).
%   5) Pressure is strongly temperature dependent. The Al2O3 isopleths have
%      an approximate slope of 4 kbar per 100 degreeC, so a reliable and
%      independent temperature estimate is essential (pp. 700-701).
%   6) Garnet and Orthopyroxene must represent local equilibrium. Zoned
%      corona textures can yield unrealistic pressures if adjacent,
%      compositionally corresponding domains are not analysed (p. 700).
%   7) Thermometry and barometry are assumed to record comparable closure
%      temperatures (p. 701).
%
% Equation (5) is shown with an approximate calibration uncertainty of
% +/-1 kbar on p. 700. This does not include potentially much larger errors
% caused by uncertain temperature, low Orthopyroxene Al, disequilibrium, or
% extrapolation in pressure, temperature, or composition.
%
% The original paper does not define numerical upper limits for Fe3+, Cr,
% or Mn, and does not define a numerical lower limit for Orthopyroxene Al.
% For diagnostic warnings only, this implementation flags:
%
%   Orthopyroxene X_Al_M1 < 0.02
%   Garnet or Orthopyroxene Cr > 0.05 apfu
%   Garnet or Orthopyroxene Mn > 0.05 apfu
%   Any finite positive Fe3_cation_apfu value
%
% These numerical thresholds are implementation-specific warning criteria,
% not calibration boundaries specified by Harley and Green (1982).
%
% This implementation issues non-stopping fprintf warnings when:
%   1) finite input temperature is outside 800-1200 degreeC,
%   2) finite calculated pressure is outside 5-20 kbar,
%   3) a composition-related limitation is detected,
%   4) a required calculation input contains NaN, or
%   5) a calculated pressure is NaN or Inf.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Opx    : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must include the
% following normalized-cation variables:
%
%   Garnet variables used directly in the equation:
%     Fe_cation_apfu           % treated as Fe2+ by this implementation
%     Mg_cation_apfu
%     Ca_cation_apfu
%
%   Orthopyroxene variables used directly in the equation:
%     Al_cation_apfu           % cations normalized to 6 oxygens
%     Fe_cation_apfu           % treated as Fe2+ by this implementation
%     Mg_cation_apfu
%
%   Optional variables retained in the output when present:
%     Si_cation_apfu
%     Ti_cation_apfu
%     Fe3_cation_apfu
%     Mn_cation_apfu
%     Cr_cation_apfu
%     Ca_cation_apfu           % optional for Orthopyroxene output only
%     Al_cation_apfu           % optional for Garnet output only
%
% Finite values used in the calculation must be non-negative. NaN is
% allowed, retained as missing, propagated through the calculation, and
% reported by non-stopping fprintf warnings. Inf values and finite negative
% values are rejected. Zero is retained; if it makes a logarithm, ratio, or
% denominator undefined, the resulting NaN or Inf is retained and reported.
%
% IMPORTANT: Fe_cation_apfu is treated as Fe2+ in the calculation, matching
% the original implementation. If Fe3_cation_apfu is present and positive,
% a warning is printed because Harley and Green (1982) identify elevated
% Fe3+ as a limitation (p. 700). Fe3+ is not subtracted from total Fe here.
%
% No liquid composition is used by this barometer. Therefore, exclusion of
% Liq F and Cl from cationTotal_liq and from Liq NaN warnings is not
% applicable to this function.
%
% -------------------------------------------------------------------------
% BAROMETER FORMULATION (Eq. 5, p. 700)
%
%   X_Al_M1 = Al_opx / 2
%
%   X_Fe_opx = Fe2+_opx / (Mg_opx + Fe2+_opx)
%
%   X_Ca_grt = Ca_grt / (Ca_grt + Mg_grt + Fe2+_grt)
%   X_Mg_grt = Mg_grt / (Ca_grt + Mg_grt + Fe2+_grt)
%   X_Fe_grt = Fe2+_grt / (Ca_grt + Mg_grt + Fe2+_grt)
%
%   K = X_Al_M1*(1-X_Al_M1) / (1-X_Ca_grt)^3
%
%   DeltaVr = -[183.3 + 178.98*X_Al_M1/(1-X_Al_M1)]
%              cal kbar^-1 mol^-1
%
%   P(kbar) = [(R*ln(K)-2.93)*T(K) + 5650
%              + 5157*(1-X_Al_M1)*(1-2*X_Al_M1)*X_Fe_opx
%              - 6300*(X_Ca_grt*X_Fe_grt + X_Ca_grt^2)] / DeltaVr
%
%   R = 1.987 cal mol^-1 K^-1
%
% The pressure equation is implemented without an additional leading minus
% sign, consistent with Eq. (5). The Garnet Ca correction uses X_Fe_grt,
% not X_Mg_grt, in the cross term shown in Eq. (5).
%
% -------------------------------------------------------------------------
% Syntax:
%   results = HarleyGreen1982(rawdata_struct, T_degreeC)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Opx tables
%   T_degreeC      : temperature in degreeC; non-negative numeric scalar or
%                    vector. NaN is allowed and retained; Inf is prohibited.
%
% Output:
%   results : table containing one row per temperature value for every
%             user-selected Garnet-Orthopyroxene pair.
%

%% Input validation
% Accept scalar or vector temperature input so that both fixed-temperature
% and temperature-range launchers use the same barometer implementation.
if nargin < 2
    error('HarleyGreen1982 requires (rawdata_struct, T_degreeC).');
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
% Extract the required tables from the input struct. The source tables are
% not modified; selected rows are read during each calculation.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Opx') || ~istable(rawdata_struct.Opx)
    error('rawdata_struct must contain table: rawdata_struct.Opx');
end

dataset_grt = rawdata_struct.Garnet;
dataset_opx = rawdata_struct.Opx;

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Store each selected-pair result in a cell buffer and concatenate once after
% the interactive loop. This avoids repeated growth of the output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Experimental calibration envelope stated in the paper abstract (p. 697).
calibrationT_min_degreeC = 800;
calibrationT_max_degreeC = 1200;
calibrationP_min_kbar = 5;
calibrationP_max_kbar = 20;

temperatureOutsideCalibration = isfinite(T_degreeC) & ...
    (T_degreeC < calibrationT_min_degreeC | ...
     T_degreeC > calibrationT_max_degreeC);
temperatureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
    % The first table column is used only as the displayed data identifier.
    dataCodes_grt = dataset_grt{:, 1};

    [selectedIdx_grt, ok] = listdlg( ...
        'PromptString', 'Please select the Garnet data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_grt)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_grt)
        disp('Selection canceled');
        break;
    end

    selectedCode_grt = dataCodes_grt(selectedIdx_grt);
    disp(['Garnet selected: ' char(string(selectedCode_grt))]);

    % ----- Orthopyroxene selection -----
    disp('=== Step 4: Selecting a data code from the list (Opx) ===');

    dataCodes_opx = dataset_opx{:, 1};

    [selectedIdx_opx, ok] = listdlg( ...
        'PromptString', 'Please select the Opx data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_opx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_opx)
        disp('Selection canceled');
        break;
    end

    selectedCode_opx = dataCodes_opx(selectedIdx_opx);
    disp(['Opx selected: ' char(string(selectedCode_opx))]);

    % ----- Calculation -----
    % Garnet and Orthopyroxene are selected independently; row numbers are
    % not assumed to correspond between the two tables.
    disp('=== Step 5: Calculating the pressure ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_opx = dataset_opx(selectedIdx_opx, :);

    % List NaN only in variables used directly by Eq. (5). NaN values are
    % not changed and do not stop the calculation.
    nanInputNames = findNaNInputs(selectedData_grt, selectedData_opx, T_degreeC);

    % Reject Inf and finite negative values. Zero and NaN are retained.
    validateNonNegativeInputs(selectedData_grt, selectedData_opx);

    row = calcPressure(selectedData_grt, selectedData_opx, T_degreeC);

    % Repeat identifiers for all temperatures in the current calculation.
    row.dataCode_garnet = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_opx = repmat(string(selectedCode_opx), height(row), 1);
    row = movevars(row, {'dataCode_garnet','dataCode_opx'}, 'Before', 1);

    % Store one block per selected mineral pair. Expand the cell buffer only
    % when its current capacity has been exhausted.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the computed pressure range for immediate feedback.
    disp('--------------------------------------------------');
    disp('=== Pressure was calculated: ===');

    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' ...
            char(string(selectedCode_opx)) ': P = ' ...
            num2str(row.P_kbar) ' kbar']);
    else
        disp([char(string(selectedCode_grt)) ' & ' ...
            char(string(selectedCode_opx)) ': P = ' ...
            num2str(row.P_kbar(1)) ' to ' ...
            num2str(row.P_kbar(end)) ' kbar']);
    end

    % Input temperature is common to all selected pairs, so this warning is
    % printed only once after the first completed calculation.
    if any(temperatureOutsideCalibration) && ~temperatureWarningIssued
        finiteTemperature = T_degreeC(isfinite(T_degreeC));
        fprintf(2, ...
            ['WARNING: Input temperature is outside the Harley and Green (1982) ' ...
             'experimental calibration envelope of 800-1200 degreeC (p. 697). ' ...
             '%d of %d finite temperature point(s) are outside the range; ' ...
             'input finite range = %.4g-%.4g degreeC.\n'], ...
            sum(temperatureOutsideCalibration), ...
            numel(finiteTemperature), ...
            min(finiteTemperature), ...
            max(finiteTemperature));
        temperatureWarningIssued = true;
    end

    % Warn when finite calculated pressures lie outside the experimental
    % pressure calibration envelope.
    finitePressure = isfinite(row.P_kbar);
    pressureOutsideCalibration = finitePressure & ...
        (row.P_kbar < calibrationP_min_kbar | ...
         row.P_kbar > calibrationP_max_kbar);

    if any(pressureOutsideCalibration)
        finitePressureValues = row.P_kbar(finitePressure);
        fprintf(2, ...
            ['WARNING: Calculated pressure is outside the Harley and Green (1982) ' ...
             'experimental calibration envelope of 5-20 kbar (p. 697). ' ...
             '%d of %d finite pressure point(s) are outside the range; ' ...
             'calculated finite range = %.4g-%.4g kbar for %s & %s.\n'], ...
            sum(pressureOutsideCalibration), ...
            sum(finitePressure), ...
            min(finitePressureValues), ...
            max(finitePressureValues), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)));
    end

    % Warn when the mathematical domain of Eq. (5) is not satisfied.
    equationDomainInvalid = ~row.isWithinEquationDomain;
    if any(equationDomainInvalid)
        fprintf(2, ...
            ['WARNING: One or more Eq. (5) domain checks failed for %s & %s. ' ...
             'Required conditions include finite component fractions, ' ...
             '0 < X_Al_M1 < 1, finite K > 0, and a non-zero finite DeltaVr. %d of %d result point(s) are outside ' ...
             'the equation domain. NaN or Inf results are retained.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)), ...
            sum(equationDomainInvalid), ...
            numel(equationDomainInvalid));
    end

    % Warn for low Orthopyroxene Al using the implementation-specific
    % diagnostic threshold documented above.
    X_Al_M1_scalar = row.X_Al_M1(1);
    if isfinite(X_Al_M1_scalar) && X_Al_M1_scalar < 0.02
        fprintf(2, ...
            ['WARNING: Orthopyroxene X_Al_M1 = %.4g is below the diagnostic ' ...
             'threshold 0.02 for %s & %s. Harley and Green (1982) state that ' ...
             'very low Orthopyroxene Al causes considerable pressure ' ...
             'uncertainty (p. 700); 0.02 is not a paper-defined limit.\n'], ...
            X_Al_M1_scalar, ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)));
    end

    % Warn when Cr or Mn exceeds the diagnostic thresholds. The paper does
    % not provide explicit numerical limits.
    compositionValues = struct( ...
        'Cr_grt', row.Cr_grt(1), ...
        'Cr_opx', row.Cr_opx(1), ...
        'Mn_grt', row.Mn_grt(1), ...
        'Mn_opx', row.Mn_opx(1));

    if (isfinite(compositionValues.Cr_grt) && compositionValues.Cr_grt > 0.05) || ...
            (isfinite(compositionValues.Cr_opx) && compositionValues.Cr_opx > 0.05)
        fprintf(2, ...
            ['WARNING: Cr exceeds the implementation-specific diagnostic ' ...
             'threshold 0.05 apfu for %s & %s (Cr_grt = %.4g, Cr_opx = %.4g). ' ...
             'High Cr3+ in Garnet is identified as a limitation by Harley and ' ...
             'Green (1982, p. 700); 0.05 apfu is not a paper-defined limit.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)), ...
            compositionValues.Cr_grt, ...
            compositionValues.Cr_opx);
    end

    if (isfinite(compositionValues.Mn_grt) && compositionValues.Mn_grt > 0.05) || ...
            (isfinite(compositionValues.Mn_opx) && compositionValues.Mn_opx > 0.05)
        fprintf(2, ...
            ['WARNING: Mn exceeds the implementation-specific diagnostic ' ...
             'threshold 0.05 apfu for %s & %s (Mn_grt = %.4g, Mn_opx = %.4g). ' ...
             'High Mn2+ is identified as a limitation by Harley and Green ' ...
             '(1982, p. 700); 0.05 apfu is not a paper-defined limit.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)), ...
            compositionValues.Mn_grt, ...
            compositionValues.Mn_opx);
    end

    % Warn whenever a positive Fe3+ value is supplied because the equation
    % still treats Fe_cation_apfu as Fe2+.
    Fe3_grt_scalar = row.Fe3_grt(1);
    Fe3_opx_scalar = row.Fe3_opx(1);
    if (isfinite(Fe3_grt_scalar) && Fe3_grt_scalar > 0) || ...
            (isfinite(Fe3_opx_scalar) && Fe3_opx_scalar > 0)
        fprintf(2, ...
            ['WARNING: Positive Fe3_cation_apfu was supplied for %s & %s ' ...
             '(Fe3_grt = %.4g, Fe3_opx = %.4g). Fe_cation_apfu is still ' ...
             'treated as Fe2+ in this implementation. Elevated Fe3+ is a ' ...
             'stated limitation of the barometer (Harley and Green, 1982, p. 700).\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)), ...
            Fe3_grt_scalar, ...
            Fe3_opx_scalar);
    end

    % List exact calculation inputs containing NaN. For vector temperature
    % input, NaN element indices are included.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the barometer input(s) for %s & %s: %s.\n' ...
             '         NaN values were retained and propagated; the calculated ' ...
             'pressure may remain NaN.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Retain and report all non-finite calculated results.
    invalidPressure = ~isfinite(row.P_kbar);
    if any(invalidPressure)
        fprintf(2, ...
            ['WARNING: Non-finite pressure values were calculated for %s & %s ' ...
             '(%d of %d points; NaN: %d, Inf: %d).\n' ...
             '         These values remain in the output table, and the calculation ' ...
             'has not been stopped.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)), ...
            sum(invalidPressure), ...
            numel(row.P_kbar), ...
            sum(isnan(row.P_kbar)), ...
            sum(isinf(row.P_kbar)));
    end

    % Negative calculated pressure is retained for diagnosis but is outside
    % the physical and calibration domain.
    negativePressure = isfinite(row.P_kbar) & row.P_kbar < 0;
    if any(negativePressure)
        fprintf(2, ...
            ['WARNING: Negative finite pressure was calculated for %s & %s ' ...
             '(%d of %d points). The values were retained for diagnostic purposes.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_opx)), ...
            sum(negativePressure), ...
            numel(row.P_kbar));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'HarleyGreen1982', ...
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

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function nanInputNames = findNaNInputs(data_garnet, data_opx, T_degreeC)
% findNaNInputs
% Return names of Eq. (5) inputs containing NaN. NaN values do not cause an
% error and are not replaced by zero.

maxNames = 7;
nanInputBuffer = strings(maxNames, 1);
nNanInputs = 0;

nanTemperatureIndices = find(isnan(T_degreeC));
if ~isempty(nanTemperatureIndices)
    nNanInputs = nNanInputs + 1;
    indexText = strjoin(string(nanTemperatureIndices.'), ',');
    nanInputBuffer(nNanInputs) = "T_degreeC(indices=" + indexText + ")";
end

garnetVariables = {'Ca_cation_apfu', 'Mg_cation_apfu', 'Fe_cation_apfu'};
for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ismember(variableName, data_garnet.Properties.VariableNames)
        variableValue = data_garnet.(variableName);
        if any(isnan(variableValue(:)))
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = "Garnet." + string(variableName);
        end
    end
end

opxVariables = {'Al_cation_apfu', 'Mg_cation_apfu', 'Fe_cation_apfu'};
for i = 1:numel(opxVariables)
    variableName = opxVariables{i};
    if ismember(variableName, data_opx.Properties.VariableNames)
        variableValue = data_opx.(variableName);
        if any(isnan(variableValue(:)))
            nNanInputs = nNanInputs + 1;
            nanInputBuffer(nNanInputs) = "Opx." + string(variableName);
        end
    end
end

nanInputNames = nanInputBuffer(1:nNanInputs);

end

function validateNonNegativeInputs(data_garnet, data_opx)
% validateNonNegativeInputs
% Reject Inf and finite negative values in cation variables used directly by
% Eq. (5). Zero and NaN are intentionally allowed and retained.

maxNames = 6;
invalidInputBuffer = strings(maxNames, 1);
nInvalidInputs = 0;

garnetVariables = {'Ca_cation_apfu', 'Mg_cation_apfu', 'Fe_cation_apfu'};
for i = 1:numel(garnetVariables)
    variableName = garnetVariables{i};
    if ~ismember(variableName, data_garnet.Properties.VariableNames)
        error('Garnet table must contain variable: %s', variableName);
    end
    variableValue = data_garnet.(variableName);
    if any(isinf(variableValue(:)) | ...
            (isfinite(variableValue(:)) & variableValue(:) < 0))
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = "Garnet." + string(variableName);
    end
end

opxVariables = {'Al_cation_apfu', 'Mg_cation_apfu', 'Fe_cation_apfu'};
for i = 1:numel(opxVariables)
    variableName = opxVariables{i};
    if ~ismember(variableName, data_opx.Properties.VariableNames)
        error('Opx table must contain variable: %s', variableName);
    end
    variableValue = data_opx.(variableName);
    if any(isinf(variableValue(:)) | ...
            (isfinite(variableValue(:)) & variableValue(:) < 0))
        nInvalidInputs = nInvalidInputs + 1;
        invalidInputBuffer(nInvalidInputs) = "Opx." + string(variableName);
    end
end

if nInvalidInputs > 0
    invalidInputNames = invalidInputBuffer(1:nInvalidInputs);
    error(['HarleyGreen1982: calculation inputs must be non-negative. ' ...
           'NaN is allowed, but Inf and finite negative value(s) are prohibited. ' ...
           'Invalid input(s): ' char(strjoin(invalidInputNames, ', ')) '.']);
end

end

function row = calcPressure(data_garnet, data_opx, T_degreeC)
% calcPressure
% Compute Harley and Green (1982) pressure for one Garnet row and one
% Orthopyroxene row at one or more input temperatures. NaN values propagate
% through the calculation and are never converted to zero.
%
% Inputs:
%   data_garnet : 1-row Garnet table
%   data_opx    : 1-row Orthopyroxene table
%   T_degreeC   : scalar or vector temperature in degreeC
%
% Output:
%   row : table with one row per temperature value

T_degreeC = T_degreeC(:);
nT = numel(T_degreeC);
T_K = T_degreeC + 273.15;

% Physical constant and calibration constants from Eq. (5), p. 700.
R_cal = 1.987;
termH_scalar = 5650;
P_uncertainty_scalar = 1.0;

% Extract one-row cation data. Required columns are checked explicitly;
% optional columns are returned as NaN when absent and never set to zero.
grt = prepareGarnetRow(data_garnet);
opx = prepareOpxRow(data_opx);

% Composition variables used by Eq. (5). NaN values propagate naturally.
X_Al_M1_scalar = opx.Al ./ 2;
X_Fe_opx_scalar = opx.Fe2 ./ (opx.Mg + opx.Fe2);

garnetDenominator = grt.Ca + grt.Mg + grt.Fe2;
X_Ca_grt_scalar = grt.Ca ./ garnetDenominator;
X_Mg_grt_scalar = grt.Mg ./ garnetDenominator;
X_Fe_grt_scalar = grt.Fe2 ./ garnetDenominator;

K_HG1982_scalar = ...
    (X_Al_M1_scalar .* (1 - X_Al_M1_scalar)) ./ ...
    ((1 - X_Ca_grt_scalar).^3);

% Preserve real-valued diagnostics without allowing a negative K to create
% complex output. K = 0 gives lnK = -Inf; NaN and negative K give NaN.
if isnan(K_HG1982_scalar) || K_HG1982_scalar < 0
    lnK_scalar = NaN;
else
    lnK_scalar = log(K_HG1982_scalar);
end

DeltaVr_scalar = -(183.3 + 178.98 .* ...
    (X_Al_M1_scalar ./ (1 - X_Al_M1_scalar)));

termFeOpx_scalar = 5157 .* ...
    (1 - X_Al_M1_scalar) .* ...
    (1 - 2 .* X_Al_M1_scalar) .* ...
    X_Fe_opx_scalar;

% Eq. (5) uses X_Ca_grt*X_Fe_grt + X_Ca_grt^2.
termCaGrt_scalar = -6300 .* ...
    (X_Ca_grt_scalar .* X_Fe_grt_scalar + X_Ca_grt_scalar.^2);

% Expand composition-dependent scalars to the temperature-vector length.
Si_grt = repmat(grt.Si, nT, 1);
Ti_grt = repmat(grt.Ti, nT, 1);
Al_grt = repmat(grt.Al, nT, 1);
FeT_grt = repmat(grt.FeT, nT, 1);
Fe2_grt = repmat(grt.Fe2, nT, 1);
Fe3_grt = repmat(grt.Fe3, nT, 1);
Mg_grt = repmat(grt.Mg, nT, 1);
Ca_grt = repmat(grt.Ca, nT, 1);
Mn_grt = repmat(grt.Mn, nT, 1);
Cr_grt = repmat(grt.Cr, nT, 1);

Si_opx = repmat(opx.Si, nT, 1);
Ti_opx = repmat(opx.Ti, nT, 1);
Al_opx = repmat(opx.Al, nT, 1);
FeT_opx = repmat(opx.FeT, nT, 1);
Fe2_opx = repmat(opx.Fe2, nT, 1);
Fe3_opx = repmat(opx.Fe3, nT, 1);
Mg_opx = repmat(opx.Mg, nT, 1);
Ca_opx = repmat(opx.Ca, nT, 1);
Mn_opx = repmat(opx.Mn, nT, 1);
Cr_opx = repmat(opx.Cr, nT, 1);

X_Al_M1 = repmat(X_Al_M1_scalar, nT, 1);
X_Fe_opx = repmat(X_Fe_opx_scalar, nT, 1);
X_Ca_grt = repmat(X_Ca_grt_scalar, nT, 1);
X_Mg_grt = repmat(X_Mg_grt_scalar, nT, 1);
X_Fe_grt = repmat(X_Fe_grt_scalar, nT, 1);
K_HG1982 = repmat(K_HG1982_scalar, nT, 1);
lnK = repmat(lnK_scalar, nT, 1);
DeltaVr = repmat(DeltaVr_scalar, nT, 1);
termH = repmat(termH_scalar, nT, 1);
termFeOpx = repmat(termFeOpx_scalar, nT, 1);
termCaGrt = repmat(termCaGrt_scalar, nT, 1);
R_output = repmat(R_cal, nT, 1);
P_uncertainty = repmat(P_uncertainty_scalar, nT, 1);

% Pressure calculation from Eq. (5). No finite-value replacement is used:
% NaN and Inf values propagate and remain in the output table.
termRT = (R_cal .* lnK - 2.93) .* T_K;
numerator = termRT + termH + termFeOpx + termCaGrt;
P_kbar = numerator ./ DeltaVr;

% Diagnostic flags. These flags do not alter calculated values.
isWithinEquationDomain = ...
    isfinite(X_Al_M1) & X_Al_M1 > 0 & X_Al_M1 < 1 & ...
    isfinite(X_Fe_opx) & X_Fe_opx >= 0 & X_Fe_opx <= 1 & ...
    isfinite(X_Ca_grt) & X_Ca_grt >= 0 & X_Ca_grt <= 1 & ...
    isfinite(X_Fe_grt) & X_Fe_grt >= 0 & X_Fe_grt <= 1 & ...
    isfinite(K_HG1982) & K_HG1982 > 0 & ...
    isfinite(DeltaVr) & abs(DeltaVr) > 1e-12;

isWithinCalibrationTRange = ...
    isfinite(T_degreeC) & T_degreeC >= 800 & T_degreeC <= 1200;

isWithinCalibrationPRange = ...
    isfinite(P_kbar) & P_kbar >= 5 & P_kbar <= 20;

isLowAlDiagnostic = isfinite(X_Al_M1) & X_Al_M1 < 0.02;
isHighCrDiagnostic = ...
    (isfinite(Cr_grt) & Cr_grt > 0.05) | ...
    (isfinite(Cr_opx) & Cr_opx > 0.05);
isHighMnDiagnostic = ...
    (isfinite(Mn_grt) & Mn_grt > 0.05) | ...
    (isfinite(Mn_opx) & Mn_opx > 0.05);
isFe3Present = ...
    (isfinite(Fe3_grt) & Fe3_grt > 0) | ...
    (isfinite(Fe3_opx) & Fe3_opx > 0);

isApplicable = ...
    isWithinEquationDomain & ...
    isWithinCalibrationTRange & ...
    isWithinCalibrationPRange & ...
    ~isLowAlDiagnostic & ...
    ~isHighCrDiagnostic & ...
    ~isHighMnDiagnostic & ...
    ~isFe3Present;

% Build a diagnostic warning string for backward-compatible table output.
warningText = strings(nT, 1);
for i = 1:nT
    currentWarning = "";
    if ~isWithinEquationDomain(i)
        currentWarning = currentWarning + "Outside equation domain; ";
    end
    if isfinite(T_degreeC(i)) && ~isWithinCalibrationTRange(i)
        currentWarning = currentWarning + "Temperature outside 800-1200 degreeC; ";
    end
    if isfinite(P_kbar(i)) && ~isWithinCalibrationPRange(i)
        currentWarning = currentWarning + "Pressure outside 5-20 kbar; ";
    end
    if isLowAlDiagnostic(i)
        currentWarning = currentWarning + "Very low Opx Al; ";
    end
    if isHighCrDiagnostic(i)
        currentWarning = currentWarning + "Cr diagnostic threshold exceeded; ";
    end
    if isHighMnDiagnostic(i)
        currentWarning = currentWarning + "Mn diagnostic threshold exceeded; ";
    end
    if isFe3Present(i)
        currentWarning = currentWarning + "Fe3+ present while total Fe is treated as Fe2+; ";
    end
    if isnan(P_kbar(i))
        currentWarning = currentWarning + "Calculated pressure is NaN; ";
    elseif isinf(P_kbar(i))
        currentWarning = currentWarning + "Calculated pressure is Inf; ";
    end
    warningText(i) = currentWarning;
end

% Pack outputs using pre-sized vectors of equal height.
row = table();
row.T_degreeC = T_degreeC;
row.T_K = T_K;
row.R_cal_mol_K = R_output;

row.Si_grt = Si_grt;
row.Ti_grt = Ti_grt;
row.Al_grt = Al_grt;
row.FeT_grt = FeT_grt;
row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.Mg_grt = Mg_grt;
row.Ca_grt = Ca_grt;
row.Mn_grt = Mn_grt;
row.Cr_grt = Cr_grt;

row.Si_opx = Si_opx;
row.Ti_opx = Ti_opx;
row.Al_opx = Al_opx;
row.FeT_opx = FeT_opx;
row.Fe2_opx = Fe2_opx;
row.Fe3_opx = Fe3_opx;
row.Mg_opx = Mg_opx;
row.Ca_opx = Ca_opx;
row.Mn_opx = Mn_opx;
row.Cr_opx = Cr_opx;

row.X_Al_M1 = X_Al_M1;
row.X_Fe_opx = X_Fe_opx;
row.X_Ca_grt = X_Ca_grt;
row.X_Mg_grt = X_Mg_grt;
row.X_Fe_grt = X_Fe_grt;

% Backward-compatible Garnet fraction aliases retained from the original.
row.X_Ca_ga = X_Ca_grt;
row.X_Mg_ga = X_Mg_grt;
row.X_Fe_ga = X_Fe_grt;

row.K_HG1982 = K_HG1982;
row.lnK = lnK;
row.DeltaVr_cal_kbar_mol = DeltaVr;

row.term_RT = termRT;
row.term_H = termH;
row.term_Fe_opx = termFeOpx;
row.term_Ca_grt = termCaGrt;
row.numerator = numerator;

row.P_kbar = P_kbar;
row.P_uncertainty_calibration_kbar = P_uncertainty;
row.isWithinEquationDomain = isWithinEquationDomain;
row.isWithinCalibrationTRange = isWithinCalibrationTRange;
row.isWithinCalibrationPRange = isWithinCalibrationPRange;
row.isLowAlDiagnostic = isLowAlDiagnostic;
row.isHighCrDiagnostic = isHighCrDiagnostic;
row.isHighMnDiagnostic = isHighMnDiagnostic;
row.isFe3Present = isFe3Present;
row.isApplicable = isApplicable;
row.warning = warningText;

end

function grt = prepareGarnetRow(data_garnet)
% prepareGarnetRow
% Extract one-row Garnet cation data. Required variables must exist;
% optional variables are retained as NaN when absent. Fe_cation_apfu is
% treated as Fe2+ and Fe3+ is not subtracted.

if height(data_garnet) ~= 1
    error('Garnet input must be a 1-row table.');
end

grt = struct();
grt.Si = getVarOrNaN(data_garnet, 'Si_cation_apfu');
grt.Ti = getVarOrNaN(data_garnet, 'Ti_cation_apfu');
grt.Al = getVarOrNaN(data_garnet, 'Al_cation_apfu');
grt.FeT = getVarOrError(data_garnet, 'Fe_cation_apfu', 'Garnet');
grt.Fe2 = grt.FeT;
grt.Fe3 = getVarOrNaN(data_garnet, 'Fe3_cation_apfu');
grt.Mg = getVarOrError(data_garnet, 'Mg_cation_apfu', 'Garnet');
grt.Ca = getVarOrError(data_garnet, 'Ca_cation_apfu', 'Garnet');
grt.Mn = getVarOrNaN(data_garnet, 'Mn_cation_apfu');
grt.Cr = getVarOrNaN(data_garnet, 'Cr_cation_apfu');

validateExtractedCations(grt, 'Garnet');

end

function opx = prepareOpxRow(data_opx)
% prepareOpxRow
% Extract one-row Orthopyroxene cation data. Required variables must exist;
% optional variables are retained as NaN when absent. Fe_cation_apfu is
% treated as Fe2+ and Fe3+ is not subtracted.

if height(data_opx) ~= 1
    error('Opx input must be a 1-row table.');
end

opx = struct();
opx.Si = getVarOrNaN(data_opx, 'Si_cation_apfu');
opx.Ti = getVarOrNaN(data_opx, 'Ti_cation_apfu');
opx.Al = getVarOrError(data_opx, 'Al_cation_apfu', 'Opx');
opx.FeT = getVarOrError(data_opx, 'Fe_cation_apfu', 'Opx');
opx.Fe2 = opx.FeT;
opx.Fe3 = getVarOrNaN(data_opx, 'Fe3_cation_apfu');
opx.Mg = getVarOrError(data_opx, 'Mg_cation_apfu', 'Opx');
opx.Ca = getVarOrNaN(data_opx, 'Ca_cation_apfu');
opx.Mn = getVarOrNaN(data_opx, 'Mn_cation_apfu');
opx.Cr = getVarOrNaN(data_opx, 'Cr_cation_apfu');

validateExtractedCations(opx, 'Opx');

end

function validateExtractedCations(cationStruct, mineralLabel)
% validateExtractedCations
% Reject Inf and finite negative values in extracted cation variables.
% Missing optional variables and explicit NaN values are retained.

fieldNames = fieldnames(cationStruct);
for i = 1:numel(fieldNames)
    value = cationStruct.(fieldNames{i});
    if ~isscalar(value)
        error('%s variable %s must be scalar in a 1-row table.', ...
            mineralLabel, fieldNames{i});
    end
    if isinf(value) || (isfinite(value) && value < 0)
        error('%s contains an invalid negative or Inf value for %s.', ...
            mineralLabel, fieldNames{i});
    end
end

end

function value = getVarOrError(tbl, variableName, mineralLabel)
% getVarOrError
% Retrieve a required scalar table variable without altering NaN. The name
% indicates that a missing column is an error; NaN values are permitted.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = tbl.(variableName);
if ~isscalar(value)
    error('Variable %s must be scalar in a 1-row table.', variableName);
end

end

function value = getVarOrNaN(tbl, variableName)
% getVarOrNaN
% Retrieve an optional scalar variable. Missing optional variables are
% represented by NaN, never by zero.

if ismember(variableName, tbl.Properties.VariableNames)
    value = tbl.(variableName);
    if ~isscalar(value)
        error('Variable %s must be scalar in a 1-row table.', variableName);
    end
else
    value = NaN;
end

end
