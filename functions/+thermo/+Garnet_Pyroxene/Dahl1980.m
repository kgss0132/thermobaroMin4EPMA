function results = Dahl1980(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Pyroxene/Dahl1980.m
% Tested with MATLAB R2024b
%
% Garnet-clinopyroxene Fe2+-Mg exchange thermometer
% Dahl, P.S. (1980)
% American Mineralogist, 65, 852-866
% DOI: None
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one
% Clinopyroxene analysis (selected by the user from tables) and calculates
% temperature using the Dahl (1980) garnet-clinopyroxene Fe2+-Mg exchange
% thermometer.
%
% Pressure may be supplied as either a scalar or a vector. Therefore, this
% function can be called from both startThermoCalc_fixedP and
% startThermoCalc_rangeP. For each selected Grt-Cpx pair, the output table
% contains one row per pressure value.
%
% The function is designed for repeated calculations: after each run it asks
% whether you want to compute another Grt-Cpx pair, and appends results into
% a single output table.
%
% -------------------------------------------------------------------------
% EMPIRICAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Dahl (1980) is not a broad experimental calibration. Equation (19) was
% derived by multiple linear regression of 13 natural Grt-Cpx pairs from
% mafic gneisses and iron formations in two areas of the Ruby Range,
% Montana. The regression used independently estimated peak metamorphic
% conditions (pp. 855-859):
%
%   Kelly area        : 745 +/- 50 degreeC; 7.2 +/- 1.2 kbar
%   Carter Creek area : 675 +/- 45 degreeC; 6.2 +/- 1.2 kbar
%
% Dahl (1980) did not define formal hard P-T calibration limits. This
% implementation uses the union of the quoted uncertainty intervals as a
% conservative empirical-support range for non-stopping warnings:
%
%   Temperature : 630-795 degreeC
%   Pressure    : 5.0-8.4 kbar
%
% These values are practical warning limits, not an experimental calibration
% envelope. The natural samples compared in Figs. 6 and 7 span approximately
% 4.5-8.5 kbar, but this broader interval is not a direct calibration range
% (Limitations, p. 864).
%
% The 13-sample Grt-Cpx regression data span approximately (Table 4, p. 859):
%
%   (XFe_Grt - XMg_Grt) : 0.303-0.592
%   XCa_Grt             : 0.159-0.541
%   XMn_Grt             : 0.012-0.322
%
% Important application notes from Dahl (1980):
%   1) Equation (19) was proposed mainly for upper-amphibolite to
%      lower-granulite facies rocks, especially rocks with Cpx low in Na and
%      Al and rocks rich in Mn (p. 862; conclusion on p. 865).
%   2) The calibration assumes that Grt and Cpx attained equilibrium, that
%      pairs from each area equilibrated at similar temperatures, and that
%      retrograde alteration did not erase the peak compositions. Dahl used
%      rim compositions for the calibration (pp. 853 and 855).
%   3) Composition corrections are essential. Applying an Fe-Mg exchange
%      thermometer without adequate composition correction can produce
%      false temperature variations within a single isotherm (p. 862).
%   4) The Ruby Range Cpx compositions contain very low Na, Al, Fe3+, and Ti.
%      The model does not explicitly account for variable Fe3+ in garnet;
%      variable Na, Al, Ca, or Mn in pyroxene; or non-ideal Fe-Mg mixing
%      between pyroxene M1-M2 sites (pp. 857 and 864-865).
%   5) Most analyses used total Fe because Fe3+ was small, but a sample with
%      appreciable Fe3+ was corrected. The exchange reaction itself is
%      defined for Fe2+, so independently estimated Fe2+ should be used when
%      Fe3+ is significant (p. 853; Fig. 6 caption on p. 861).
%   6) The regression explains approximately 89 percent of the observed
%      variance. Dahl notes possible over-correction in a high-Ca sample and
%      a high-Fe sample, and emphasizes the limitations of the small natural
%      dataset (p. 857).
%   7) Application outside the Ruby Range assumes that the standard Gibbs
%      free energy term remains nearly constant over a broad temperature
%      interval. Dahl argues that this is reasonable, but it is still an
%      extrapolation assumption rather than broad experimental validation
%      (pp. 862-863).
%   8) Dahl's Grt-Opx formulation is only a relative thermometer. This file
%      implements the Grt-Cpx absolute formulation in Eq. (19), not the
%      Grt-Opx formulation (pp. 863-864).
%
% This implementation issues non-stopping fprintf messages when:
%   1) input pressure is outside 5.0-8.4 kbar,
%   2) a finite calculated temperature is outside 630-795 degreeC,
%   3) a required thermometer input contains NaN, or
%   4) a calculated temperature is NaN or Inf. In this case, the actual
%      thermometer inputs and identified calculation-domain causes are also
%      displayed.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Cpx    : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns must contain
% normalized cation data.
%
% Variables used directly by the thermometer:
%   Garnet table:
%     Fe_cation_apfu         % Fe2+ in garnet
%     Mg_cation_apfu
%     Ca_cation_apfu
%     Mn_cation_apfu
%
%   Clinopyroxene table:
%     Fe_cation_apfu         % Fe2+ in clinopyroxene
%     Mg_cation_apfu
%
% Optional variables retained in the output when available:
%   Fe3_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   Na_cation_apfu
%   K_cation_apfu
%
% IMPORTANT Fe note:
% The exchange reaction is defined using Fe2+. Therefore,
% Fe3_cation_apfu is retained for reference but is not added to
% Fe_cation_apfu. If Fe_cation_apfu contains total Fe rather than Fe2+, Fe2+
% should be estimated before this function is used.
%
% Negative finite values in variables used by the thermometer are not
% permitted. Zero values are retained; if they make the equation
% mathematically undefined, the affected result is returned as NaN and a
% non-stopping warning is printed. NaN values are never replaced by zero;
% they propagate through the calculation and remain in the output.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% Fe2+-Mg distribution coefficient:
%
%   KD = (Fe2+/Mg)_Grt / (Fe2+/Mg)_Cpx
%
% Garnet X-site mole fractions:
%
%   XFe_Grt = Fe2+ / (Fe2+ + Mg + Ca + Mn)
%   XMg_Grt = Mg   / (Fe2+ + Mg + Ca + Mn)
%   XCa_Grt = Ca   / (Fe2+ + Mg + Ca + Mn)
%   XMn_Grt = Mn   / (Fe2+ + Mg + Ca + Mn)
%
% Original equation (Dahl, 1980, Eq. 19 on p. 862):
%
%          2324 + 0.022*P_bar + 1509*(XFe_Grt - XMg_Grt)
%   T(K) = ----------------------------------------------------
%                         1.987*ln(KD)
%          + 2810*XCa_Grt + 2855*XMn_Grt
%
% Written on one line:
%
%   T(K) = [2324 + 0.022*P_bar ...
%           + 1509*(XFe_Grt - XMg_Grt) ...
%           + 2810*XCa_Grt + 2855*XMn_Grt] ...
%          / [1.987*ln(KD)]
%
% Dividing the numerator coefficients by 1.987 gives an equivalent rounded
% form whose constant is approximately 1169.6 (commonly rounded to 1170),
% not 170. This implementation uses Eq. (19) directly to avoid coefficient-
% rounding and transcription errors.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Dahl1980(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Cpx tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Grt-Cpx pair
%

%% Input validation
if nargin < 2
    error('Dahl1980 requires (rawdata_struct, P_kbar).');
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
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Cpx') || ~istable(rawdata_struct.Cpx)
    error('rawdata_struct must contain table: rawdata_struct.Cpx');
end

dataset_grt = rawdata_struct.Garnet;
dataset_cpx = rawdata_struct.Cpx;

validateRequiredVariables(dataset_grt, dataset_cpx);

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Result blocks are buffered and concatenated once after the interactive
% loop, avoiding repeated reallocation of the complete output table.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Conservative empirical-support limits derived from the quoted Ruby Range
% P-T estimates and their uncertainties (Dahl, 1980, p. 855).
supportT_min_degC = 630;
supportT_max_degC = 795;
supportP_min_kbar = 5.0;
supportP_max_kbar = 8.4;

pressureOutsideSupport = ...
    P_kbar < supportP_min_kbar | P_kbar > supportP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Garnet) ===');

while true
    % ----- Garnet selection -----
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

    % ----- Clinopyroxene selection -----
    disp('=== Step 4: Selecting a data code from the list (Clinopyroxene) ===');

    dataCodes_cpx = dataset_cpx{:, 1};

    [selectedIdx_cpx, ok] = listdlg( ...
        'PromptString', 'Please select the Clinopyroxene data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_cpx)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_cpx)
        disp('Selection canceled');
        break;
    end

    selectedCode_cpx = dataCodes_cpx(selectedIdx_cpx);
    disp(['Clinopyroxene selected: ' char(string(selectedCode_cpx))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_cpx = dataset_cpx(selectedIdx_cpx, :);

    nanInputNames = findNaNInputs(selectedData_grt, selectedData_cpx);
    validateNonNegativeInputs(selectedData_grt, selectedData_cpx);

    row = calcTemp(selectedData_grt, selectedData_cpx, P_kbar);

    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_cpx = repmat(string(selectedCode_cpx), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_cpx'}, 'Before', 1);

    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo calculated temperature values.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' char(string(selectedCode_cpx)) ...
            ': ' num2str(row.T_C(1)) ' to ' num2str(row.T_C(end)) ' degreeC']);
    end

    % Pressure warning is printed only once because the same pressure vector
    % is used for every mineral pair in this function call.
    if any(pressureOutsideSupport) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the empirical P-T support range ' ...
             'derived from the Ruby Range calibration of Dahl (1980): ' ...
             '5.0-8.4 kbar. %d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar. This is a practical warning range, ' ...
             'not a formal experimental calibration limit.\n'], ...
            sum(pressureOutsideSupport), ...
            numel(P_kbar), ...
            min(P_kbar), ...
            max(P_kbar));
        pressureWarningIssued = true;
    end

    % Warn for finite temperatures outside the empirical support interval.
    finiteTemperature = isfinite(row.T_C);
    temperatureOutsideSupport = finiteTemperature & ...
        (row.T_C < supportT_min_degC | row.T_C > supportT_max_degC);

    if any(temperatureOutsideSupport)
        finiteValues = row.T_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the empirical P-T support ' ...
             'range derived from the Ruby Range calibration of Dahl (1980): ' ...
             '630-795 degreeC. %d of %d finite temperature point(s) are outside ' ...
             'the range; calculated finite range = %.4g-%.4g degreeC for %s & %s. ' ...
             'This is a practical warning range, not a formal experimental calibration limit.\n'], ...
            sum(temperatureOutsideSupport), ...
            sum(finiteTemperature), ...
            min(finiteValues), ...
            max(finiteValues), ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)));
    end

    % Report explicitly stored NaN thermometer inputs. Calculation continues
    % and NaN values are not replaced by zero.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, and NaN was not replaced by zero.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Preserve non-finite results and print both raw inputs and identified
    % calculation-domain causes.
    invalidTemperature = ~isfinite(row.T_C);
    if any(invalidTemperature)
        nonFiniteCauses = findNonFiniteCauses(row);

        fprintf(2, ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
                    '(%d of %d points; NaN: %d, Inf: %d).\n' ...
                    '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_grt)), ...
            char(string(selectedCode_cpx)), ...
            sum(invalidTemperature), ...
            numel(row.T_C), ...
            sum(isnan(row.T_C)), ...
            sum(isinf(row.T_C)));

        fprintf(2, ['         Thermometer inputs used: ' ...
                    'Garnet.Fe_cation_apfu=%s, Garnet.Mg_cation_apfu=%s, ' ...
                    'Garnet.Ca_cation_apfu=%s, Garnet.Mn_cation_apfu=%s, ' ...
                    'Cpx.Fe_cation_apfu=%s, Cpx.Mg_cation_apfu=%s.\n'], ...
            formatNumericValue(row.Fe2_grt(1)), ...
            formatNumericValue(row.Mg_grt(1)), ...
            formatNumericValue(row.Ca_grt(1)), ...
            formatNumericValue(row.Mn_grt(1)), ...
            formatNumericValue(row.Fe2_cpx(1)), ...
            formatNumericValue(row.Mg_cpx(1)));

        if isempty(nonFiniteCauses)
            fprintf(2, ['         No explicit NaN, zero, Inf, or invalid intermediate value ' ...
                        'was identified; inspect the stored intermediate variables.\n']);
        else
            fprintf(2, '         Identified cause(s): %s.\n', ...
                char(strjoin(nonFiniteCauses, '; ')));
        end
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Dahl1980', ...
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
function validateRequiredVariables(dataset_grt, dataset_cpx)
% validateRequiredVariables
% Verify all columns required by the Dahl (1980) Grt-Cpx equation.

requiredGrt = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
requiredCpx = {'Fe_cation_apfu', 'Mg_cation_apfu'};

missingNames = strings(numel(requiredGrt) + numel(requiredCpx), 1);
nMissing = 0;

for i = 1:numel(requiredGrt)
    if ~ismember(requiredGrt{i}, dataset_grt.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Garnet." + string(requiredGrt{i});
    end
end

for i = 1:numel(requiredCpx)
    if ~ismember(requiredCpx{i}, dataset_cpx.Properties.VariableNames)
        nMissing = nMissing + 1;
        missingNames(nMissing) = "Cpx." + string(requiredCpx{i});
    end
end

if nMissing > 0
    missingNames = missingNames(1:nMissing);
    error(['Dahl1980: required variable(s) are missing: ' ...
        char(strjoin(missingNames, ', ')) '.']);
end

end


function nanInputNames = findNaNInputs(data_grt, data_cpx)
% findNaNInputs
% Return names of required thermometer inputs containing NaN.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

nanInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nNaN = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    variableValue = data_grt.(variableName);
    if any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanInputNames(nNaN) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(cpxVariables)
    variableName = cpxVariables{i};
    variableValue = data_cpx.(variableName);
    if any(isnan(variableValue(:)))
        nNaN = nNaN + 1;
        nanInputNames(nNaN) = "Cpx." + string(variableName);
    end
end

nanInputNames = nanInputNames(1:nNaN);

end


function validateNonNegativeInputs(data_grt, data_cpx)
% validateNonNegativeInputs
% Stop when a finite required thermometer input is negative. Zero and NaN
% remain available to the calculation and result diagnostics.

grtVariables = {'Fe_cation_apfu', 'Mg_cation_apfu', ...
    'Ca_cation_apfu', 'Mn_cation_apfu'};
cpxVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};

invalidInputNames = strings(numel(grtVariables) + numel(cpxVariables), 1);
nInvalid = 0;

for i = 1:numel(grtVariables)
    variableName = grtVariables{i};
    variableValue = data_grt.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Garnet." + string(variableName);
    end
end

for i = 1:numel(cpxVariables)
    variableName = cpxVariables{i};
    variableValue = data_cpx.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        nInvalid = nInvalid + 1;
        invalidInputNames(nInvalid) = "Cpx." + string(variableName);
    end
end

if nInvalid > 0
    invalidInputNames = invalidInputNames(1:nInvalid);
    error(['Dahl1980: thermometer inputs must be >= 0. ' ...
        'Negative finite value(s) were found in: ' ...
        char(strjoin(invalidInputNames, ', ')) '.']);
end

end


function row = calcTemp(data_grt, data_cpx, P_kbar)
% calcTemp
% Compute Dahl (1980) Grt-Cpx temperatures for one mineral pair and a
% scalar or vector of pressures.

P_kbar = P_kbar(:);
nP = numel(P_kbar);
P_bar = P_kbar .* 1000;

row = table();
row.P_kbar = P_kbar;
row.P_bar = P_bar;

% --- Extract and expand garnet cations ---
Fe2_grt = repmat(getRequiredValue(data_grt, 'Fe_cation_apfu', 'Garnet'), nP, 1);
Fe3_grt = repmat(getOptionalValue(data_grt, 'Fe3_cation_apfu', 0, 'Garnet'), nP, 1);
Mg_grt  = repmat(getRequiredValue(data_grt, 'Mg_cation_apfu', 'Garnet'), nP, 1);
Mn_grt  = repmat(getRequiredValue(data_grt, 'Mn_cation_apfu', 'Garnet'), nP, 1);
Ca_grt  = repmat(getRequiredValue(data_grt, 'Ca_cation_apfu', 'Garnet'), nP, 1);
Al_grt  = repmat(getOptionalValue(data_grt, 'Al_cation_apfu', 0, 'Garnet'), nP, 1);
Si_grt  = repmat(getOptionalValue(data_grt, 'Si_cation_apfu', 0, 'Garnet'), nP, 1);

% --- Extract and expand clinopyroxene cations ---
Fe2_cpx = repmat(getRequiredValue(data_cpx, 'Fe_cation_apfu', 'Cpx'), nP, 1);
Fe3_cpx = repmat(getOptionalValue(data_cpx, 'Fe3_cation_apfu', 0, 'Cpx'), nP, 1);
Mg_cpx  = repmat(getRequiredValue(data_cpx, 'Mg_cation_apfu', 'Cpx'), nP, 1);
Mn_cpx  = repmat(getOptionalValue(data_cpx, 'Mn_cation_apfu', 0, 'Cpx'), nP, 1);
Ca_cpx  = repmat(getOptionalValue(data_cpx, 'Ca_cation_apfu', 0, 'Cpx'), nP, 1);
Ti_cpx  = repmat(getOptionalValue(data_cpx, 'Ti_cation_apfu', 0, 'Cpx'), nP, 1);
Al_cpx  = repmat(getOptionalValue(data_cpx, 'Al_cation_apfu', 0, 'Cpx'), nP, 1);
Si_cpx  = repmat(getOptionalValue(data_cpx, 'Si_cation_apfu', 0, 'Cpx'), nP, 1);
Na_cpx  = repmat(getOptionalValue(data_cpx, 'Na_cation_apfu', 0, 'Cpx'), nP, 1);
K_cpx   = repmat(getOptionalValue(data_cpx, 'K_cation_apfu', 0, 'Cpx'), nP, 1);

% --- Fe2+-Mg distribution coefficient ---
FeMg_grt = Fe2_grt ./ Mg_grt;
FeMg_cpx = Fe2_cpx ./ Mg_cpx;
KD = FeMg_grt ./ FeMg_cpx;
lnKD = log(KD);

% --- Garnet X-site mole fractions ---
xSiteSum_grt = Fe2_grt + Mg_grt + Mn_grt + Ca_grt;
XFe_grt = Fe2_grt ./ xSiteSum_grt;
XMg_grt = Mg_grt ./ xSiteSum_grt;
XCa_grt = Ca_grt ./ xSiteSum_grt;
XMn_grt = Mn_grt ./ xSiteSum_grt;
XFeMinusXMg_grt = XFe_grt - XMg_grt;

% --- Dahl (1980), Eq. (19), p. 862 ---
numerator = 2324 ...
    + 0.022 .* P_bar ...
    + 1509 .* XFeMinusXMg_grt ...
    + 2810 .* XCa_grt ...
    + 2855 .* XMn_grt;
denominator = 1.987 .* lnKD;
T_K = numerator ./ denominator;

% Return NaN for domain-invalid cases while preserving all intermediate
% values for diagnosis.
validDomain = isfinite(KD) & KD > 0 ...
    & isfinite(XFe_grt) & XFe_grt >= 0 ...
    & isfinite(XMg_grt) & XMg_grt >= 0 ...
    & isfinite(XCa_grt) & XCa_grt >= 0 ...
    & isfinite(XMn_grt) & XMn_grt >= 0 ...
    & isfinite(numerator) & isfinite(denominator) ...
    & denominator ~= 0;
T_K(~validDomain) = NaN;
T_C = T_K - 273.15;

% --- Pack outputs ---
% FeUsed is retained for compatibility with the original output variable
% set, but now correctly equals Fe2+ rather than Fe2+ + Fe3+.
row.Fe2_grt = Fe2_grt;
row.Fe3_grt = Fe3_grt;
row.FeUsed_grt = Fe2_grt;
row.Mg_grt = Mg_grt;
row.Mn_grt = Mn_grt;
row.Ca_grt = Ca_grt;
row.Al_grt = Al_grt;
row.Si_grt = Si_grt;

row.Fe2_cpx = Fe2_cpx;
row.Fe3_cpx = Fe3_cpx;
row.FeUsed_cpx = Fe2_cpx;
row.Mg_cpx = Mg_cpx;
row.Mn_cpx = Mn_cpx;
row.Ca_cpx = Ca_cpx;
row.Ti_cpx = Ti_cpx;
row.Al_cpx = Al_cpx;
row.Si_cpx = Si_cpx;
row.Na_cpx = Na_cpx;
row.K_cpx = K_cpx;

row.FeMg_grt = FeMg_grt;
row.FeMg_cpx = FeMg_cpx;
row.XFe_grt = XFe_grt;
row.XMg_grt = XMg_grt;
row.XCa_grt = XCa_grt;
row.XMn_grt = XMn_grt;
row.XFeMinusXMg_grt = XFeMinusXMg_grt;
row.KD = KD;
row.lnKD = lnKD;
row.numerator = numerator;
row.denominator = denominator;
row.T_K = T_K;
row.T_C = T_C;

end


function nonFiniteCauses = findNonFiniteCauses(row)
% findNonFiniteCauses
% Identify raw-input and derived-variable conditions that can produce a NaN
% or Inf temperature. These diagnostics do not alter stored values.

maximumCauses = 24;
nonFiniteCauses = strings(maximumCauses, 1);
nCauses = 0;

inputValues = {row.Fe2_grt, row.Mg_grt, row.Ca_grt, row.Mn_grt, ...
    row.Fe2_cpx, row.Mg_cpx};
inputNames = {"Garnet.Fe_cation_apfu", "Garnet.Mg_cation_apfu", ...
    "Garnet.Ca_cation_apfu", "Garnet.Mn_cation_apfu", ...
    "Cpx.Fe_cation_apfu", "Cpx.Mg_cation_apfu"};

for i = 1:numel(inputValues)
    value = inputValues{i};
    if any(isnan(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is NaN";
    elseif any(isinf(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is Inf";
    elseif any(value(:) == 0)
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = inputNames{i} + " is zero";
    end
end

derivedValues = {row.KD, row.XFe_grt, row.XMg_grt, row.XCa_grt, ...
    row.XMn_grt, row.numerator, row.denominator};
derivedNames = {"KD", "XFe_grt", "XMg_grt", "XCa_grt", ...
    "XMn_grt", "temperature numerator", "temperature denominator"};

for i = 1:numel(derivedValues)
    value = derivedValues{i};
    if any(~isfinite(value(:)))
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is non-finite";
    elseif any(value(:) == 0) && ...
            (derivedNames{i} == "KD" || derivedNames{i} == "temperature denominator")
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is zero";
    elseif any(value(:) < 0) && derivedNames{i} ~= "temperature numerator" ...
            && derivedNames{i} ~= "temperature denominator"
        nCauses = nCauses + 1;
        nonFiniteCauses(nCauses) = derivedNames{i} + " is negative";
    end
end

nonFiniteCauses = nonFiniteCauses(1:nCauses);

end


function textValue = formatNumericValue(value)
% formatNumericValue
% Format a scalar numeric value for a compact diagnostic message.

if isnan(value)
    textValue = 'NaN';
elseif isinf(value) && value > 0
    textValue = 'Inf';
elseif isinf(value) && value < 0
    textValue = '-Inf';
else
    textValue = sprintf('%.8g', value);
end

end


function value = getRequiredValue(data_tbl, variableName, mineralLabel)
% getRequiredValue
% Extract a required scalar numeric value. NaN is returned unchanged.

if ~ismember(variableName, data_tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = data_tbl.(variableName);
if ~isnumeric(value) || ~isscalar(value)
    error('%s.%s must be a numeric scalar in the selected row.', ...
        mineralLabel, variableName);
end

end


function value = getOptionalValue(data_tbl, variableName, defaultValue, mineralLabel)
% getOptionalValue
% Use defaultValue only when an optional column is absent. An explicitly
% stored NaN is returned unchanged and is never converted to zero.

if ismember(variableName, data_tbl.Properties.VariableNames)
    value = data_tbl.(variableName);
    if ~isnumeric(value) || ~isscalar(value)
        error('%s.%s must be a numeric scalar in the selected row.', ...
            mineralLabel, variableName);
    end
else
    value = defaultValue;
end

end
