function results = Lanari2014(rawdata_struct, P_kbar)
% functions/+thermo/+Chlorite/Lanari2014.m
% Tested with MATLAB R2024b
%
% Semi-empirical single-Chlorite geothermometers Chl(1) and Chl(2)
% Lanari, P., Wagner, T. and Vidal, O. (2014)
% Contributions to Mineralogy and Petrology, 167, 968
% DOI: https://doi.org/10.1007/s00410-014-0968-8
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively selects one Chlorite analysis at a time and
% calculates temperature using the two semi-empirical geothermometers
% proposed by Lanari et al. (2014):
%
%   Chl(1): explicitly uses Fe3+ and requires reliable grain-scale
%           Fe3+/sumFe information.
%
%   Chl(2): assumes that all Fe is Fe2+ as a calibration convention and
%           includes an explicit pressure term.
%
% The function is designed for repeated calculations. After each run it asks
% whether another Chlorite analysis should be calculated and stores all
% results in a single output table.
%
% Both a scalar pressure and a pressure vector are accepted. Therefore, this
% function can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. One output row is returned for every supplied
% pressure value and every user-selected Chlorite analysis.
%
% -------------------------------------------------------------------------
% CALIBRATION RANGE AND APPLICATION NOTES
%
% Lanari et al. (2014) developed a thermodynamic activity-composition model
% for di-trioctahedral Chlorite in the FMASH system and derived two
% simplified semi-empirical geothermometers based on the equilibrium:
%
%   2 clinochlore + 3 sudoite =
%       4 amesite + 4 H2O + 7 quartz
%
% The thermodynamic model was calibrated using 271 published natural
% Chlorite analyses with independent P-T estimates. The natural calibration
% data cover approximately 100-600 degreeC and 1-20 kbar (model calibration,
% Table 3 and Fig. 4, pp. 6-7).
%
% Chl(2) is explicitly reported to predict correct temperatures over:
%
%   Temperature : 100-500 degreeC
%   Pressure    : 1-20 kbar
%
% These limits are stated in the Chl(2) calibration discussion on
% pp. 11-12 and reiterated in the Conclusions on pp. 15-16.
%
% Chl(1) does not have a separately stated formal P-T range. It was
% calibrated using a smaller natural data subset for which Fe3+/sumFe was
% measured or estimated. This implementation uses 100-500 degreeC and
% 1-20 kbar only as conservative reference-screening limits for Chl(1), not
% as separately published formal Chl(1) calibration limits.
%
% IMPORTANT APPLICATION LIMITATIONS:
%
%   - Chl(1) requires accurate Fe3+/sumFe information for individual
%     Chlorite grains. Lanari et al. (2014) describe this as a serious
%     practical limitation, especially in polymetamorphic rocks containing
%     multiple Chlorite generations (p. 11).
%
%   - Bulk-sample Mossbauer values may average grains or generations with
%     different Fe3+ contents and formation temperatures, which can increase
%     scatter in Chl(1) results (pp. 10-11).
%
%   - Chl(2) assumes all Fe is Fe2+ only as a semi-empirical calibration
%     convention. This assumption does not demonstrate that natural
%     Chlorite contains no Fe3+. The temperature-dependent Fe3+ effect is
%     interpreted as being implicitly incorporated into the Chl(2)
%     calibration (pp. 7-9 and 11).
%
%   - Both thermometers are based on Chlorite-quartz-water equilibrium and
%     assume unit water activity when deriving ln(K) (p. 9). Application to
%     quartz-absent, silica-undersaturated, water-poor, high-salinity, or
%     mixed-volatile systems requires additional caution.
%
%   - The thermodynamic framework is primarily the
%     FeO-MgO-Al2O3-SiO2-H2O system. Chlorites with large Mn, Cr, Ni, Ca,
%     Na, K, or other non-FMASH components may lie outside the calibrated
%     compositional model (abstract and model description, pp. 1-4).
%
%   - Development of the Chl(LWV) composition model was restricted to
%     Chlorites with Si < 3 apfu on a 14-oxygen basis (p. 3). This
%     implementation therefore prints a caution for Si >= 3 apfu.
%
%   - Analyses with very small sudoite component, expressed as the vacancy
%     variable z < 0.045, were removed from the model-calibration data set
%     (p. 7). Temperatures are retained but reported as compositional
%     extrapolations when finite z is below 0.045.
%
%   - Positive activities of amesite, clinochlore, and sudoite are required.
%     Invalid site fractions, non-positive activities, or physically
%     inconsistent x-y-z values make ln(K) undefined. Such outputs are
%     retained as NaN and reported by fprintf rather than stopping the full
%     interactive workflow.
%
%   - The semi-empirical thermometers use only ideal on-site activities.
%     This simplification performs well for temperature estimation but
%     cannot replace the complete non-ideal Chl(LWV) model in internally
%     consistent phase-equilibrium calculations (pp. 11-12).
%
%   - Different Chlorite generations, detrital cores, recrystallized rims,
%     and overgrowths must be distinguished petrographically. The
%     application example on pp. 13-15 resolves distinct 200-275 degreeC
%     and 300-400 degreeC domains within partially recrystallized grains.
%
%   - Figures 8 and 9 use an average temperature uncertainty of
%     approximately +/-50 degreeC because complete uncertainty information
%     was unavailable for many source studies (pp. 10-11). Individual
%     calculated values should not be interpreted at unrealistically high
%     precision.
%
% This implementation issues non-stopping fprintf messages when:
%   1) pressure is outside 1-20 kbar,
%   2) a finite Chl(2) temperature is outside 100-500 degreeC,
%   3) a finite Chl(1) temperature is outside the conservative
%      100-500 degreeC reference interval,
%   4) Si >= 3 apfu,
%   5) finite z < 0.045,
%   6) Chl(1) Fe3+ information is unavailable or non-finite,
%   7) a calculation input contains NaN,
%   8) the calculated composition or activities are invalid, or
%   9) a non-finite temperature is calculated.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain either:
%   rawdata_struct.Chlorite : table
% or
%   rawdata_struct.Chl      : table
%
% The FIRST column of the table is treated as an identifier ("data code")
% displayed in the selection dialog.
%
% Required variables, normalized on a 14-oxygen basis:
%   Si_cation_apfu
%   Al_cation_apfu
%   Fe_cation_apfu      % total Fe; used as Fe2+ by Chl(2)
%   Mg_cation_apfu
%
% Optional calculation variables:
%   Ti_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%
% Optional Fe3+ descriptors for Chl(1):
%   Fe3_cation_apfu
%   Fe2_cation_apfu
%   Fe3_ratio_chl       % Fe3+/(Fe2+ + Fe3+), from 0 to 1
%
% Fe3_cation_apfu is preferred when its column exists. If Fe2_cation_apfu
% is absent and finite Fe3_cation_apfu and Fe_cation_apfu are available,
% Fe2 is calculated as total Fe minus Fe3. If an existing Fe2, Fe3, or
% Fe3-ratio value is NaN, it remains NaN and Chl(1) is not silently
% calculated by replacing the missing value with zero.
%
% Optional trace variables retained in the output:
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Cr_cation_apfu
%
% All finite cation values must be greater than or equal to zero. Negative
% values are prohibited. NaN values are retained as missing values and are
% never replaced by zero when the corresponding column exists. Optional
% variables are assigned zero only when their columns are absent.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Structural formula and composition variables on a 14-oxygen basis:
%
%   x = Fe2/(Fe2 + Mg)
%
%   AlIV = 4 - (Si + Ti)
%   AlVI = Al_total - AlIV
%   R1   = Na + K
%   h    = 0.5*(AlVI - AlIV + Fe3 + R1)
%   z    = h
%
%   Al_M4  = 1 - Fe3
%   Al_M23 = 2*h
%   Al_M1  = AlVI - Al_M23 - Al_M4
%   y      = Al_M1
%
% Site fractions:
%
%   xFe_M23 = x*(1 - 0.5*z)
%   xMg_M23 = (1 - x)*(1 - 0.5*z)
%   xAl_M23 = 0.5*z
%
%   xFe_M1 = x*(1 - y - z)
%   xMg_M1 = (1 - x)*(1 - y - z)
%   xAl_M1 = y
%   xVa_M1 = z
%
%   xAl_T2 = 0.5*(1 + y)
%   xSi_T2 = 0.5*(1 - y)
%
% NOTE:
% The T2 fractions above follow AlIV = 1 + y in Lanari et al. (2014,
% Eq. 33, p. 8). The signs are intentionally opposite to those in the
% original unmodified Lanari2014.m supplied for this task.
%
% Ideal activities:
%
%   a_ames = xMg_M23^4*xAl_M1*xAl_T2^2
%
%   a_clin = 4*xMg_M23^4*xMg_M1*xAl_T2*xSi_T2
%
%   a_sud  = 64*xAl_M23^2*xMg_M23^2*xVa_M1*xAl_T2*xSi_T2
%
%   lnK = ln[a_ames^4/(a_clin^2*a_sud^3)]
%
% Corrected temperature equations:
%
%   Chl(1):
%     T(degreeC) = 172341/(315.149 - R*lnK) - 273.15
%
%   Chl(2):
%     T(degreeC) =
%       (203093 + 4996.99*P_kbar)/(455.782 - R*lnK) - 273.15
%
% These denominator signs follow Eqs. 37-40 on pp. 10-11. The original
% unmodified Lanari2014.m used addition instead of subtraction.
%
% R = 8.31446261815324 J mol^-1 K^-1.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Lanari2014(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing a Chlorite or Chl table
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Chlorite analysis
%

%% Input validation
if nargin < 2
    error('Lanari2014 requires (rawdata_struct, P_kbar).');
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

calibrationT_min_degC = 100;
calibrationT_max_degC = 500;
calibrationP_min_kbar = 1;
calibrationP_max_kbar = 20;
primaryModelSi_max = 3;
minimumCalibrationZ = 0.045;

pressureWarningIssued = false;
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

    disp([char(string(selectedCode_chl)) ': T_Chl2 = ' ...
        formatFiniteRange(row.T_Chl2_deg) ' degreeC']);

    if row.Chl1_available(1)
        disp([char(string(selectedCode_chl)) ': T_Chl1 = ' ...
            formatFiniteRange(row.T_Chl1_deg) ' degreeC']);
    else
        disp([char(string(selectedCode_chl)) ...
            ': T_Chl1 = NaN (usable grain-scale Fe3+ information unavailable)']);
    end

    % Pressure range warning is common to all selected analyses and is
    % therefore displayed only once per function call.
    if ~pressureWarningIssued
        outsidePressureRange = ...
            P_kbar < calibrationP_min_kbar | ...
            P_kbar > calibrationP_max_kbar;

        if any(outsidePressureRange)
            fprintf(2, ...
                ['WARNING: %d of %d pressure point(s) are outside the ' ...
                 'Lanari et al. (2014) Chl(2) calibration range of ' ...
                 '1-20 kbar (pp. 11-12 and 15-16); input range = ' ...
                 '%.4g-%.4g kbar. Chl(1) has no separately stated formal ' ...
                 'pressure range, but 1-20 kbar is retained as a ' ...
                 'conservative natural-data reference interval.\n'], ...
                sum(outsidePressureRange), ...
                numel(P_kbar), ...
                min(P_kbar), ...
                max(P_kbar));
        end

        pressureWarningIssued = true;
    end

    % Print general assumptions once. These conditions cannot be determined
    % from a single EPMA row.
    if ~applicationCautionIssued
        fprintf(2, ...
            ['WARNING: Lanari et al. (2014) Chl(1) and Chl(2) assume ' ...
             'Chlorite-quartz-water equilibrium and unit water activity. ' ...
             'The script cannot verify quartz saturation, water activity, ' ...
             'equilibrium, Chlorite generation, recrystallization state, ' ...
             'or whether non-FMASH components are minor. These conditions ' ...
             'must be assessed independently (pp. 9 and 12-15).\n']);
        applicationCautionIssued = true;
    end

    % Chl(2) temperature calibration warning.
    finiteChl2 = isfinite(row.T_Chl2_deg);
    outsideChl2Temperature = finiteChl2 & ...
        (row.T_Chl2_deg < calibrationT_min_degC | ...
         row.T_Chl2_deg > calibrationT_max_degC);

    if any(outsideChl2Temperature)
        finiteValues = row.T_Chl2_deg(finiteChl2);
        fprintf(2, ...
            ['WARNING: Chl(2) temperature is outside the Lanari et al. ' ...
             '(2014) calibration range of 100-500 degreeC ' ...
             '(pp. 11-12 and 15-16). %d of %d finite point(s) are ' ...
             'outside; calculated finite range = %.4g-%.4g degreeC ' ...
             'for %s.\n'], ...
            sum(outsideChl2Temperature), ...
            sum(finiteChl2), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    end

    % Chl(1) has no separately stated formal range. Use the same interval
    % only as a conservative reference screen.
    finiteChl1 = isfinite(row.T_Chl1_deg);
    outsideChl1Reference = finiteChl1 & ...
        (row.T_Chl1_deg < calibrationT_min_degC | ...
         row.T_Chl1_deg > calibrationT_max_degC);

    if any(outsideChl1Reference)
        finiteValues = row.T_Chl1_deg(finiteChl1);
        fprintf(2, ...
            ['WARNING: Chl(1) temperature is outside the conservative ' ...
             '100-500 degreeC natural-data reference interval used by ' ...
             'this implementation. Lanari et al. (2014) do not state a ' ...
             'separate formal Chl(1) temperature range. Calculated finite ' ...
             'range = %.4g-%.4g degreeC for %s.\n'], ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_chl)));
    end

    % Primary model composition restriction.
    finiteSi = isfinite(row.Si);
    outsidePrimarySiRange = finiteSi & row.Si >= primaryModelSi_max;

    if any(outsidePrimarySiRange)
        fprintf(2, ...
            ['WARNING: Si >=3 apfu was found for %s. Development of the ' ...
             'Lanari et al. (2014) Chl(LWV) composition model was ' ...
             'restricted to Chlorites with Si <3 apfu on a 14-oxygen ' ...
             'basis (p. 3). The result is a compositional extrapolation.\n'], ...
            char(string(selectedCode_chl)));
    end

    % Calibration-data screen for the sudoite/vacancy variable z.
    finiteZChl2 = isfinite(row.z_Chl2);
    lowZChl2 = finiteZChl2 & row.z_Chl2 < minimumCalibrationZ;

    if any(lowZChl2)
        fprintf(2, ...
            ['WARNING: Chl(2) vacancy variable z is below 0.045 for %s. ' ...
             'Analyses with z <0.045 were excluded from the Lanari et al. ' ...
             '(2014) model-calibration data set (p. 7). The temperature ' ...
             'is retained as a compositional extrapolation.\n'], ...
            char(string(selectedCode_chl)));
    end

    finiteZChl1 = isfinite(row.z_Chl1);
    lowZChl1 = finiteZChl1 & row.z_Chl1 < minimumCalibrationZ;

    if any(lowZChl1)
        fprintf(2, ...
            ['WARNING: Chl(1) vacancy variable z is below 0.045 for %s. ' ...
             'Analyses with z <0.045 were excluded from the Lanari et al. ' ...
             '(2014) model-calibration data set (p. 7). The temperature ' ...
             'is retained as a compositional extrapolation.\n'], ...
            char(string(selectedCode_chl)));
    end

    % Chl(1) Fe3+ availability.
    if ~row.Chl1_available(1)
        if ~row.has_Fe3_descriptor(1)
            fprintf(2, ...
                ['WARNING: Chl(1) was not calculated for %s because no ' ...
                 'Fe3_cation_apfu or Fe3_ratio_chl descriptor is present. ' ...
                 'Lanari et al. (2014, p. 11) require accurate ' ...
                 'grain-scale Fe3+/sumFe information for Chl(1). Chl(2) ' ...
                 'was calculated independently.\n'], ...
                char(string(selectedCode_chl)));
        else
            fprintf(2, ...
                ['WARNING: Chl(1) was not calculated for %s because the ' ...
                 'existing Fe2+, Fe3+, Fe3-ratio, or total-Fe information ' ...
                 'needed by Chl(1) contains NaN or is otherwise not usable. ' ...
                 'The missing value was not replaced by zero. Chl(2) was ' ...
                 'calculated independently.\n'], ...
                char(string(selectedCode_chl)));
        end
    end

    % Model-composition and activity checks.
    if any(~row.isValidComposition_Chl2)
        fprintf(2, ...
            ['WARNING: The Lanari et al. (2014) Chl(2) composition model ' ...
             'or ideal activities are invalid for %s. Site fractions, ' ...
             'x-y-z variables, and positive amesite-clinochlore-sudoite ' ...
             'activities are required. Dependent Chl(2) outputs remain ' ...
             'NaN and calculation continued.\n'], ...
            char(string(selectedCode_chl)));
    end

    if row.Chl1_available(1) && any(~row.isValidComposition_Chl1)
        fprintf(2, ...
            ['WARNING: The Lanari et al. (2014) Chl(1) composition model ' ...
             'or ideal activities are invalid for %s. Site fractions, ' ...
             'x-y-z variables, and positive amesite-clinochlore-sudoite ' ...
             'activities are required. Dependent Chl(1) outputs remain ' ...
             'NaN and calculation continued.\n'], ...
            char(string(selectedCode_chl)));
    end

    % Report NaN inputs after the displayed result. NaN values remain
    % unchanged and propagate only into calculations that depend on them.
    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in the Lanari et al. (2014) ' ...
             'calculation input(s) for %s: %s.\n' ...
             '         NaN values were retained and were not replaced by ' ...
             'zero. Dependent outputs may therefore be NaN.\n'], ...
            char(string(selectedCode_chl)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    invalidChl2Temperature = ~isfinite(row.T_Chl2_deg);
    if any(invalidChl2Temperature)
        fprintf(2, ...
            ['WARNING: Non-finite Chl(2) temperatures were calculated for ' ...
             '%s (%d of %d points; NaN: %d, Inf: %d). These values remain ' ...
             'in the output table, and calculation has not been stopped.\n'], ...
            char(string(selectedCode_chl)), ...
            sum(invalidChl2Temperature), ...
            numel(row.T_Chl2_deg), ...
            sum(isnan(row.T_Chl2_deg)), ...
            sum(isinf(row.T_Chl2_deg)));
    end

    if row.Chl1_available(1)
        invalidChl1Temperature = ~isfinite(row.T_Chl1_deg);
        if any(invalidChl1Temperature)
            fprintf(2, ...
                ['WARNING: Non-finite Chl(1) temperatures were calculated ' ...
                 'for %s (%d of %d points; NaN: %d, Inf: %d). These values ' ...
                 'remain in the output table, and calculation has not been ' ...
                 'stopped.\n'], ...
                char(string(selectedCode_chl)), ...
                sum(invalidChl1Temperature), ...
                numel(row.T_Chl1_deg), ...
                sum(isnan(row.T_Chl1_deg)), ...
                sum(isinf(row.T_Chl1_deg)));
        end
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Lanari2014', ...
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
% Return names of existing variables that are used directly in Chl(1) or
% Chl(2) and contain NaN.

variableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu', ...
    'Fe2_cation_apfu', ...
    'Fe3_cation_apfu', ...
    'Fe3_ratio_chl'};

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
% Reject finite negative cation values, infinite values, and non-numeric
% cation columns. Zero and NaN are allowed. Fe3_ratio_chl must lie between
% zero and one when finite.

cationVariableNames = { ...
    'Si_cation_apfu', ...
    'Al_cation_apfu', ...
    'Fe_cation_apfu', ...
    'Mg_cation_apfu', ...
    'Mn_cation_apfu', ...
    'Ti_cation_apfu', ...
    'Ca_cation_apfu', ...
    'Na_cation_apfu', ...
    'K_cation_apfu', ...
    'Cr_cation_apfu', ...
    'Fe2_cation_apfu', ...
    'Fe3_cation_apfu'};

displayNames = "Chlorite." + string(cationVariableNames(:));
negativeMask = false(numel(cationVariableNames), 1);
infiniteMask = false(numel(cationVariableNames), 1);
nonNumericMask = false(numel(cationVariableNames), 1);

for i = 1:numel(cationVariableNames)
    variableName = cationVariableNames{i};

    if ismember(variableName, data_chl.Properties.VariableNames)
        variableValue = data_chl.(variableName);

        if ~isnumeric(variableValue)
            nonNumericMask(i) = true;
            continue;
        end

        negativeMask(i) = any(isfinite(variableValue(:)) & variableValue(:) < 0);
        infiniteMask(i) = any(isinf(variableValue(:)));
    end
end

if any(nonNumericMask)
    error(['Lanari2014: cation variables must be numeric. ' ...
           'Non-numeric variable(s): ' ...
           char(strjoin(displayNames(nonNumericMask), ', ')) '.']);
end

if any(negativeMask)
    error(['Lanari2014: cation values must be greater than or equal to ' ...
           'zero. Negative value(s) were found in: ' ...
           char(strjoin(displayNames(negativeMask), ', ')) '.']);
end

if any(infiniteMask)
    error(['Lanari2014: infinite cation value(s) are not allowed in: ' ...
           char(strjoin(displayNames(infiniteMask), ', ')) '.']);
end

if ismember('Fe3_ratio_chl', data_chl.Properties.VariableNames)
    ratioValue = data_chl.Fe3_ratio_chl;

    if ~isnumeric(ratioValue) || ~isscalar(ratioValue)
        error('Lanari2014: Fe3_ratio_chl must be a numeric scalar.');
    end

    if isinf(ratioValue)
        error('Lanari2014: Fe3_ratio_chl cannot be infinite.');
    end

    if isfinite(ratioValue) && (ratioValue < 0 || ratioValue > 1)
        error('Lanari2014: finite Fe3_ratio_chl must be between 0 and 1.');
    end
end

end

function row = calcTemp(data_chl, P_kbar)
% calcTemp
% Compute Chl(1) and Chl(2) for one Chlorite composition and repeat scalar
% composition outputs for every supplied pressure value.

P_kbar = P_kbar(:);
nP = numel(P_kbar);

R = 8.31446261815324;

chl = prepareChloriteRow(data_chl);

% Chl(2): all total Fe is treated as Fe2+ by calibration convention.
modelChl2 = calcLanariModel(chl, false);

% Chl(1): use grain-scale Fe2+ and Fe3+ information only when usable.
if chl.hasFe3Info
    modelChl1 = calcLanariModel(chl, true);
else
    modelChl1 = emptyLanariModel();
end

% Correct signs from Lanari et al. (2014), Eqs. 37-40.
if modelChl1.isValidComposition
    denominatorChl1 = 315.149 - R .* modelChl1.lnK;
    T_Chl1_deg_scalar = 172341 ./ denominatorChl1 - 273.15;
else
    denominatorChl1 = NaN;
    T_Chl1_deg_scalar = NaN;
end
T_Chl1_K_scalar = T_Chl1_deg_scalar + 273.15;

if modelChl2.isValidComposition
    denominatorChl2_scalar = 455.782 - R .* modelChl2.lnK;
    T_Chl2_deg = ...
        (203093 + 4996.99 .* P_kbar) ./ denominatorChl2_scalar - 273.15;
else
    denominatorChl2_scalar = NaN;
    T_Chl2_deg = NaN(nP, 1);
end
T_Chl2_K = T_Chl2_deg + 273.15;

T_Chl1_deg = repmat(T_Chl1_deg_scalar, nP, 1);
T_Chl1_K = repmat(T_Chl1_K_scalar, nP, 1);

if isfinite(chl.Mg) && isfinite(chl.FeT) && (chl.Mg + chl.FeT) > 0
    Mg_number_totalFe_scalar = chl.Mg ./ (chl.Mg + chl.FeT);
else
    Mg_number_totalFe_scalar = NaN;
end

if chl.hasFe3Info && isfinite(chl.Fe2_input) && ...
        isfinite(chl.Fe3_input) && ...
        (chl.Fe2_input + chl.Fe3_input) > 0
    Fe3_ratio_input_scalar = ...
        chl.Fe3_input ./ (chl.Fe2_input + chl.Fe3_input);
else
    Fe3_ratio_input_scalar = NaN;
end

isInPressureRange = P_kbar >= 1 & P_kbar <= 20;
isInPrimarySiRangeScalar = isfinite(chl.Si) && chl.Si < 3;

isInZRangeChl2Scalar = ...
    isfinite(modelChl2.z) && modelChl2.z >= 0.045;

isInZRangeChl1Scalar = ...
    isfinite(modelChl1.z) && modelChl1.z >= 0.045;

isInTemperatureRangeChl2 = ...
    isfinite(T_Chl2_deg) & T_Chl2_deg >= 100 & T_Chl2_deg <= 500;

isInTemperatureReferenceChl1 = ...
    isfinite(T_Chl1_deg) & T_Chl1_deg >= 100 & T_Chl1_deg <= 500;

% These are numerical-screening flags. Quartz-water equilibrium, low
% non-FMASH content, grain generation, and recrystallization must be checked
% independently and cannot be established from one EPMA row.
recommended_by_Lanari2014_Chl2 = ...
    repmat(modelChl2.isValidComposition, nP, 1) & ...
    repmat(isInPrimarySiRangeScalar, nP, 1) & ...
    repmat(isInZRangeChl2Scalar, nP, 1) & ...
    isInPressureRange & ...
    isInTemperatureRangeChl2;

recommended_by_Lanari2014_Chl1 = ...
    repmat(chl.hasFe3Info, nP, 1) & ...
    repmat(modelChl1.isValidComposition, nP, 1) & ...
    repmat(isInPrimarySiRangeScalar, nP, 1) & ...
    repmat(isInZRangeChl1Scalar, nP, 1) & ...
    isInPressureRange & ...
    isInTemperatureReferenceChl1;

row = table();

row.P_kbar = P_kbar;
row.R = repmat(R, nP, 1);

% Raw inputs.
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

row.has_Fe3_descriptor = repmat(chl.hasFe3Descriptor, nP, 1);
row.has_Fe3_input = repmat(chl.hasFe3Info, nP, 1);
row.Fe2_input = repmat(chl.Fe2_input, nP, 1);
row.Fe3_input = repmat(chl.Fe3_input, nP, 1);
row.Fe3_ratio_input = repmat(Fe3_ratio_input_scalar, nP, 1);
row.Mg_number_totalFe = repmat(Mg_number_totalFe_scalar, nP, 1);

% Chl(1) model outputs.
row.Chl1_available = repmat(chl.hasFe3Info, nP, 1);
row.Fe2_Chl1 = repmat(modelChl1.Fe2, nP, 1);
row.Fe3_Chl1 = repmat(modelChl1.Fe3, nP, 1);

row.AlIV_Chl1 = repmat(modelChl1.AlIV, nP, 1);
row.AlVI_Chl1 = repmat(modelChl1.AlVI, nP, 1);
row.R1_Chl1 = repmat(modelChl1.R1, nP, 1);
row.h_Chl1 = repmat(modelChl1.h, nP, 1);
row.x_Chl1 = repmat(modelChl1.x, nP, 1);
row.y_Chl1 = repmat(modelChl1.y, nP, 1);
row.z_Chl1 = repmat(modelChl1.z, nP, 1);

row.Al_M4_Chl1 = repmat(modelChl1.Al_M4, nP, 1);
row.Al_M23_total_Chl1 = repmat(modelChl1.Al_M23_total, nP, 1);
row.Al_M1_Chl1 = repmat(modelChl1.Al_M1, nP, 1);

row.xFe_M23_Chl1 = repmat(modelChl1.xFe_M23, nP, 1);
row.xMg_M23_Chl1 = repmat(modelChl1.xMg_M23, nP, 1);
row.xAl_M23_Chl1 = repmat(modelChl1.xAl_M23, nP, 1);
row.xFe_M1_Chl1 = repmat(modelChl1.xFe_M1, nP, 1);
row.xMg_M1_Chl1 = repmat(modelChl1.xMg_M1, nP, 1);
row.xAl_M1_Chl1 = repmat(modelChl1.xAl_M1, nP, 1);
row.xVa_M1_Chl1 = repmat(modelChl1.xVa_M1, nP, 1);
row.xAl_T2_Chl1 = repmat(modelChl1.xAl_T2, nP, 1);
row.xSi_T2_Chl1 = repmat(modelChl1.xSi_T2, nP, 1);

row.a_ames_Chl1 = repmat(modelChl1.a_ames, nP, 1);
row.a_clin_Chl1 = repmat(modelChl1.a_clin, nP, 1);
row.a_sud_Chl1 = repmat(modelChl1.a_sud, nP, 1);
row.lnK_Chl1 = repmat(modelChl1.lnK, nP, 1);
row.denominator_Chl1 = repmat(denominatorChl1, nP, 1);
row.isValidComposition_Chl1 = ...
    repmat(modelChl1.isValidComposition, nP, 1);

row.T_Chl1_deg = T_Chl1_deg;
row.T_Chl1_K = T_Chl1_K;

% Chl(2) model outputs.
row.Fe2_Chl2_assumed = repmat(modelChl2.Fe2, nP, 1);
row.Fe3_Chl2_assumed = repmat(modelChl2.Fe3, nP, 1);

row.AlIV_Chl2 = repmat(modelChl2.AlIV, nP, 1);
row.AlVI_Chl2 = repmat(modelChl2.AlVI, nP, 1);
row.R1_Chl2 = repmat(modelChl2.R1, nP, 1);
row.h_Chl2 = repmat(modelChl2.h, nP, 1);
row.x_Chl2 = repmat(modelChl2.x, nP, 1);
row.y_Chl2 = repmat(modelChl2.y, nP, 1);
row.z_Chl2 = repmat(modelChl2.z, nP, 1);

row.Al_M4_Chl2 = repmat(modelChl2.Al_M4, nP, 1);
row.Al_M23_total_Chl2 = repmat(modelChl2.Al_M23_total, nP, 1);
row.Al_M1_Chl2 = repmat(modelChl2.Al_M1, nP, 1);

row.xFe_M23_Chl2 = repmat(modelChl2.xFe_M23, nP, 1);
row.xMg_M23_Chl2 = repmat(modelChl2.xMg_M23, nP, 1);
row.xAl_M23_Chl2 = repmat(modelChl2.xAl_M23, nP, 1);
row.xFe_M1_Chl2 = repmat(modelChl2.xFe_M1, nP, 1);
row.xMg_M1_Chl2 = repmat(modelChl2.xMg_M1, nP, 1);
row.xAl_M1_Chl2 = repmat(modelChl2.xAl_M1, nP, 1);
row.xVa_M1_Chl2 = repmat(modelChl2.xVa_M1, nP, 1);
row.xAl_T2_Chl2 = repmat(modelChl2.xAl_T2, nP, 1);
row.xSi_T2_Chl2 = repmat(modelChl2.xSi_T2, nP, 1);

row.a_ames_Chl2 = repmat(modelChl2.a_ames, nP, 1);
row.a_clin_Chl2 = repmat(modelChl2.a_clin, nP, 1);
row.a_sud_Chl2 = repmat(modelChl2.a_sud, nP, 1);
row.lnK_Chl2 = repmat(modelChl2.lnK, nP, 1);
row.denominator_Chl2 = repmat(denominatorChl2_scalar, nP, 1);
row.isValidComposition_Chl2 = ...
    repmat(modelChl2.isValidComposition, nP, 1);

row.T_Chl2_deg = T_Chl2_deg;
row.T_Chl2_K = T_Chl2_K;

% Generic temperature aliases use Chl(2), the pressure-dependent practical
% calibration that remains available when Fe3+ information is unknown.
row.T_deg = T_Chl2_deg;
row.T_K = T_Chl2_K;

% Comparison and applicability flags.
row.deltaT_Chl1_minus_Chl2_deg = T_Chl1_deg - T_Chl2_deg;

row.is_in_Lanari2014_P_range = isInPressureRange;
row.is_in_Lanari2014_Chl2_T_range = isInTemperatureRangeChl2;
row.is_in_Lanari2014_Chl1_T_reference_range = ...
    isInTemperatureReferenceChl1;
row.is_in_primary_Siel_range = ...
    repmat(isInPrimarySiRangeScalar, nP, 1);
row.is_in_z_calibration_range_Chl1 = ...
    repmat(isInZRangeChl1Scalar, nP, 1);
row.is_in_z_calibration_range_Chl2 = ...
    repmat(isInZRangeChl2Scalar, nP, 1);

row.Chl1_pressure_calibration_defined = false(nP, 1);
row.Chl2_pressure_calibration_defined = true(nP, 1);

row.requires_quartz_water_equilibrium_confirmation = true(nP, 1);
row.requires_water_activity_assessment = true(nP, 1);
row.requires_chlorite_generation_confirmation = true(nP, 1);
row.requires_nonFMASH_component_assessment = true(nP, 1);

row.recommended_by_Lanari2014_Chl1 = ...
    recommended_by_Lanari2014_Chl1;
row.recommended_by_Lanari2014_Chl2 = ...
    recommended_by_Lanari2014_Chl2;

end

function chl = prepareChloriteRow(data_chl)
% prepareChloriteRow
% Extract one Chlorite analysis from a 1-row table. Existing NaN values are
% retained. Optional values are assigned zero only when their columns are
% absent.

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

chl.hasFe3Descriptor = false;
chl.hasFe3Info = false;
chl.Fe2_input = NaN;
chl.Fe3_input = NaN;

hasFe3 = ismember('Fe3_cation_apfu', data_chl.Properties.VariableNames);
hasFe2 = ismember('Fe2_cation_apfu', data_chl.Properties.VariableNames);
hasFe3ratio = ismember('Fe3_ratio_chl', data_chl.Properties.VariableNames);

% Direct Fe3 apfu information takes priority when its column exists.
if hasFe3
    chl.hasFe3Descriptor = true;
    Fe3_input = getExistingScalarAllowNaN(data_chl, 'Fe3_cation_apfu');

    if hasFe2
        Fe2_input = getExistingScalarAllowNaN(data_chl, 'Fe2_cation_apfu');
    elseif isfinite(Fe3_input) && isfinite(chl.FeT)
        Fe2_input = chl.FeT - Fe3_input;

        if Fe2_input < -1.0e-10
            error(['Lanari2014: Fe3_cation_apfu exceeds finite ' ...
                   'Fe_cation_apfu.']);
        end

        if Fe2_input < 0
            Fe2_input = 0;
        end
    else
        Fe2_input = NaN;
    end

    if isfinite(Fe2_input) && isfinite(Fe3_input) && isfinite(chl.FeT)
        if abs((Fe2_input + Fe3_input) - chl.FeT) > 1.0e-6
            error(['Lanari2014: finite Fe2_cation_apfu + ' ...
                   'Fe3_cation_apfu must equal Fe_cation_apfu for Chl(1).']);
        end
    end

    chl.Fe2_input = Fe2_input;
    chl.Fe3_input = Fe3_input;
    chl.hasFe3Info = isfinite(Fe2_input) && isfinite(Fe3_input);

% Use Fe3 ratio only when a direct Fe3 column is absent.
elseif hasFe3ratio
    chl.hasFe3Descriptor = true;
    Fe3ratio = getExistingScalarAllowNaN(data_chl, 'Fe3_ratio_chl');

    if isfinite(Fe3ratio) && isfinite(chl.FeT)
        chl.Fe3_input = chl.FeT .* Fe3ratio;
        chl.Fe2_input = chl.FeT - chl.Fe3_input;
        chl.hasFe3Info = true;
    end
end

end

function model = calcLanariModel(chl, useFe3Info)
% calcLanariModel
% Calculate Lanari et al. (2014) composition variables, site fractions,
% ideal activities, and ln(K). Invalid derived compositions return NaN
% activities rather than stopping the complete workflow.

model = emptyLanariModel();

if useFe3Info
    if ~chl.hasFe3Info
        return;
    end

    Fe2 = chl.Fe2_input;
    Fe3 = chl.Fe3_input;
else
    Fe2 = chl.FeT;
    Fe3 = 0;
end

model.Fe2 = Fe2;
model.Fe3 = Fe3;

% Structural formula quantities. Do not clamp invalid finite values because
% doing so would conceal compositions outside the published model.
AlIV = 4 - (chl.Si + chl.Ti);
AlVI = chl.Al - AlIV;
R1 = chl.Na + chl.K;
h = 0.5 .* (AlVI - AlIV + Fe3 + R1);

Al_M4 = 1 - Fe3;
Al_M23_total = 2 .* h;
Al_M1 = AlVI - Al_M23_total - Al_M4;

x = Fe2 ./ (Fe2 + chl.Mg);
y = Al_M1;
z = h;

% Correct T2-site fractions from Eq. 33: AlIV = 1 + y.
xFe_M23 = x .* (1 - 0.5 .* z);
xMg_M23 = (1 - x) .* (1 - 0.5 .* z);
xAl_M23 = 0.5 .* z;

xFe_M1 = x .* (1 - y - z);
xMg_M1 = (1 - x) .* (1 - y - z);
xAl_M1 = y;
xVa_M1 = z;

xAl_T2 = 0.5 .* (1 + y);
xSi_T2 = 0.5 .* (1 - y);

model.AlIV = AlIV;
model.AlVI = AlVI;
model.R1 = R1;
model.h = h;

model.Al_M4 = Al_M4;
model.Al_M23_total = Al_M23_total;
model.Al_M1 = Al_M1;

model.x = x;
model.y = y;
model.z = z;

model.xFe_M23 = xFe_M23;
model.xMg_M23 = xMg_M23;
model.xAl_M23 = xAl_M23;

model.xFe_M1 = xFe_M1;
model.xMg_M1 = xMg_M1;
model.xAl_M1 = xAl_M1;
model.xVa_M1 = xVa_M1;

model.xAl_T2 = xAl_T2;
model.xSi_T2 = xSi_T2;

siteFractions = [ ...
    xFe_M23, xMg_M23, xAl_M23, ...
    xFe_M1, xMg_M1, xAl_M1, xVa_M1, ...
    xAl_T2, xSi_T2];

derivedQuantities = [ ...
    AlIV, AlVI, R1, h, Al_M4, Al_M23_total, Al_M1, x, y, z];

tolerance = 1.0e-10;

isFiniteAll = all(isfinite(siteFractions)) && ...
    all(isfinite(derivedQuantities));

isNonNegativeAll = all(siteFractions >= -tolerance);
isNotTooLarge = all(siteFractions <= 1 + tolerance);

isReasonableX = isfinite(x) && x >= -tolerance && x <= 1 + tolerance;
isReasonableY = isfinite(y) && y >= -tolerance && y <= 1 + tolerance;
isReasonableZ = isfinite(z) && z >= -tolerance && z <= 1 + tolerance;
isReasonableAlM4 = ...
    isfinite(Al_M4) && Al_M4 >= -tolerance && Al_M4 <= 1 + tolerance;

isReasonableAlIV = isfinite(AlIV) && AlIV >= -tolerance && AlIV <= 2 + tolerance;
isReasonableAlVI = isfinite(AlVI) && AlVI >= -tolerance;

isValidComposition = ...
    isFiniteAll && ...
    isNonNegativeAll && ...
    isNotTooLarge && ...
    isReasonableX && ...
    isReasonableY && ...
    isReasonableZ && ...
    isReasonableAlM4 && ...
    isReasonableAlIV && ...
    isReasonableAlVI;

if ~isValidComposition
    model.isValidComposition = false;
    return;
end

% Remove only tiny numerical negatives within tolerance after validation.
xMg_M23_activity = max(0, xMg_M23);
xAl_M23_activity = max(0, xAl_M23);
xMg_M1_activity = max(0, xMg_M1);
xAl_M1_activity = max(0, xAl_M1);
xVa_M1_activity = max(0, xVa_M1);
xAl_T2_activity = max(0, xAl_T2);
xSi_T2_activity = max(0, xSi_T2);

a_ames = ...
    (xMg_M23_activity.^4) .* ...
    xAl_M1_activity .* ...
    (xAl_T2_activity.^2);

a_clin = ...
    4 .* (xMg_M23_activity.^4) .* ...
    xMg_M1_activity .* ...
    xAl_T2_activity .* ...
    xSi_T2_activity;

a_sud = ...
    64 .* (xAl_M23_activity.^2) .* ...
    (xMg_M23_activity.^2) .* ...
    xVa_M1_activity .* ...
    xAl_T2_activity .* ...
    xSi_T2_activity;

model.a_ames = a_ames;
model.a_clin = a_clin;
model.a_sud = a_sud;

if ~isfinite(a_ames) || ~isfinite(a_clin) || ~isfinite(a_sud) || ...
        a_ames <= 0 || a_clin <= 0 || a_sud <= 0
    model.isValidComposition = false;
    return;
end

lnK = log((a_ames.^4) ./ ((a_clin.^2) .* (a_sud.^3)));

if ~isfinite(lnK)
    model.isValidComposition = false;
    return;
end

model.lnK = lnK;
model.isValidComposition = true;

end

function model = emptyLanariModel()
% emptyLanariModel
% Return a scalar model structure filled with NaN diagnostic values.

model = struct();

model.Fe2 = NaN;
model.Fe3 = NaN;

model.AlIV = NaN;
model.AlVI = NaN;
model.R1 = NaN;
model.h = NaN;

model.Al_M4 = NaN;
model.Al_M23_total = NaN;
model.Al_M1 = NaN;

model.x = NaN;
model.y = NaN;
model.z = NaN;

model.xFe_M23 = NaN;
model.xMg_M23 = NaN;
model.xAl_M23 = NaN;

model.xFe_M1 = NaN;
model.xMg_M1 = NaN;
model.xAl_M1 = NaN;
model.xVa_M1 = NaN;

model.xAl_T2 = NaN;
model.xSi_T2 = NaN;

model.a_ames = NaN;
model.a_clin = NaN;
model.a_sud = NaN;
model.lnK = NaN;

model.isValidComposition = false;

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
% Retrieve an optional numeric scalar. Assign zero only when the column is
% absent. Existing NaN values remain NaN.

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

function value = getExistingScalarAllowNaN(tbl, varName)
% getExistingScalarAllowNaN
% Retrieve an existing numeric scalar. NaN is retained. Infinite and finite
% negative values are prohibited. Fe3_ratio_chl is additionally restricted
% to the interval zero to one.

value = tbl.(varName);

if ~isnumeric(value) || ~isscalar(value)
    error('Variable %s must be a numeric scalar in a 1-row table.', varName);
end

if isinf(value)
    error('Variable %s contains an infinite value.', varName);
end

if strcmp(varName, 'Fe3_ratio_chl')
    if isfinite(value) && (value < 0 || value > 1)
        error('Variable Fe3_ratio_chl must be between 0 and 1 when finite.');
    end
elseif isfinite(value) && value < 0
    error('Variable %s contains a negative value.', varName);
end

end

function textValue = formatFiniteRange(values)
% formatFiniteRange
% Format one scalar or a finite numerical range for command-window output.

values = values(:);
finiteValues = values(isfinite(values));

if isempty(finiteValues)
    textValue = 'NaN';
elseif isscalar(finiteValues) || ...
        abs(max(finiteValues) - min(finiteValues)) <= eps(max(abs(finiteValues)))
    textValue = num2str(finiteValues(1));
else
    textValue = [num2str(min(finiteValues)) ' to ' num2str(max(finiteValues))];
end

end
