function results = Jackson1969(rawdata_struct, P_kbar)
% functions/+thermo/+Olivine_Spinel/Jackson1969.m
% Target environment: MATLAB R2024b
%
% Fe-Mg exchange thermometer between Olivine and Spinel
% Jackson, E.D. (1969)
% Chemical variation in coexisting chromite and olivine in chromitite
% zones of the Stillwater Complex.
% Economic Geology Monograph 4, 41-71.
% DOI: https://doi.org/10.5382/Mono.04.03
%
% Equation and application notes used in this implementation are reproduced
% in Hartmann, L.A. and Chemale-Junior, F. (2003), Appendix 2, pp. 124-125.
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Olivine analysis and one Spinel
% analysis and calculates temperature from Fe2+-Mg exchange between the two
% minerals using the Jackson (1969) relation.
%
% The user-selection, repeated-calculation, and table-output workflow follows
% Ballhaus1991.m and Fujii1978.m so that this function can be called by the
% same fixed-pressure and pressure-range launchers.
%
% Jackson's equation contains no explicit pressure term. P_kbar is retained
% for launcher compatibility and output traceability, but changing pressure
% does not change the calculated temperature for a given Olivine-Spinel pair.
%
% -------------------------------------------------------------------------
% CALIBRATION AND APPLICATION NOTES
%
% Jackson (1969) formulated the Irvine Olivine-Spinel Fe-Mg exchange model
% using thermochemical data and coexisting chromite-olivine compositions from
% the Stillwater Complex. This is an early thermochemical formulation and does
% not have a modern experimentally defined calibration range or uncertainty.
%
% Hartmann and Chemale-Junior (2003, Appendix 2, p. 125) describe the
% following compositional conditions for their application:
%   1) only small variation in Mg/(Mg+Fe2+) of Olivine around Fo90; and
%   2) chrome-spinel compositions near the
%      (Fe2+,Mg)(Cr,Fe3+)2O4 face of the Spinel compositional prism.
%
% Important cautions:
%   1) Fe3+ in Spinel must be distinguished from Fe2+. Ferric-iron
%      estimation affects both the trivalent-cation fractions and KD.
%   2) The equation does not include Ti. Ti-rich Spinel compositions may not
%      be represented adequately by the formulation.
%   3) Jackson temperatures may be too high for Olivine-Chrome Spinel pairs
%      formed during regional metamorphism of serpentinites.
%   4) Ferric-iron-rich Spinel may yield a large and unreasonable range of
%      calculated temperatures.
%   5) Altered Spinel, Ferritchromite, Cr-magnetite, magnetite-rich rims,
%      and mismatched core-rim mineral pairs should not be interpreted
%      without independent evidence of equilibrium.
%   6) Fe-Mg exchange may continue during subsolidus cooling. A calculated
%      value can therefore record re-equilibration or closure rather than
%      magmatic crystallization or peak metamorphic temperature.
%   7) The equation contains no explicit pressure correction.
%
% This implementation does not impose numerical Fe3+, Ti, pressure, or
% temperature cutoffs because Jackson (1969) and the secondary source used
% here do not define universally applicable numerical thresholds.
%
% SOURCE-CONSISTENCY NOTE:
% Hartmann and Chemale-Junior (2003) label the printed quotient as a
% temperature in degreeC. Applying the rounded values reported in their
% Table VI to the printed equation does not exactly reproduce every listed
% Table VI temperature. This implementation follows the printed Appendix 2
% equation verbatim and does not fit or modify its coefficients to reproduce
% the rounded table values.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT
% rawdata_struct must contain:
%   rawdata_struct.Olivine : table
%   rawdata_struct.Spinel  : table
%
% The FIRST column of each table is treated as an identifier displayed in
% the selection dialog. The following normalized-cation variables are used:
%
%   Olivine:
%     Mg_cation_apfu
%     Fe_cation_apfu         % Fe2+ in Olivine
%
%   Spinel:
%     Mg_cation_apfu
%     Fe_cation_apfu         % Fe2+ in Spinel
%     Fe3_cation_apfu        % Fe3+ in Spinel
%     Cr_cation_apfu
%     Al_cation_apfu
%
% Finite Mg and Fe2+ values must be strictly positive because they occur in
% ratios and logarithms. Finite Fe3+, Cr, and Al values may be zero but must
% not be negative. The trivalent-cation sum Cr + Al + Fe3+ must be positive.
% NaN values are allowed to propagate and are reported by non-stopping
% warnings.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION
%
% 1) Divalent-cation fractions
%   XMg_ol  = Mg_ol  / (Mg_ol  + Fe2_ol)
%   XFe2_ol = Fe2_ol / (Mg_ol  + Fe2_ol)
%   XMg_sp  = Mg_sp  / (Mg_sp  + Fe2_sp)
%   XFe2_sp = Fe2_sp / (Mg_sp  + Fe2_sp)
%
% 2) Trivalent-cation fractions in Spinel
%   a = Cr_sp  / (Cr_sp + Al_sp + Fe3_sp)
%   b = Al_sp  / (Cr_sp + Al_sp + Fe3_sp)
%   c = Fe3_sp / (Cr_sp + Al_sp + Fe3_sp)
%   a + b + c = 1
%
% 3) Fe-Mg exchange coefficient
%   KD_FeMg_ol_sp = ...
%       (XMg_ol * XFe2_sp) / (XFe2_ol * XMg_sp)
%
% 4) Temperature equation reproduced by Hartmann and Chemale-Junior (2003)
%    from Jackson (1969, p. 63)
%
%   numerator = 5580*a + 1018*b - 1720*c + 2400
%
%   denominator = 0.9*a + 2.56*b - 3.08*c - 1.47 ...
%                 + 1.9871*ln(KD_FeMg_ol_sp)
%
%   T(degreeC) = numerator / denominator
%
% The coefficient 1.9871 is retained exactly as printed. Pressure does not
% appear in the equation.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Jackson1969(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Olivine and Spinel tables
%   P_kbar         : finite non-negative scalar or vector; retained for
%                    launcher compatibility and output only
%
% Output:
%   results : table containing one row per supplied pressure value for every
%             user-selected Olivine-Spinel pair
%

%% Input validation
if nargin < 2
    error('Jackson1969 requires (rawdata_struct, P_kbar).');
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

if ~isfield(rawdata_struct, 'Olivine') || ~istable(rawdata_struct.Olivine)
    error('rawdata_struct must contain table: rawdata_struct.Olivine');
end
if ~isfield(rawdata_struct, 'Spinel') || ~istable(rawdata_struct.Spinel)
    error('rawdata_struct must contain table: rawdata_struct.Spinel');
end

dataset_ol = rawdata_struct.Olivine;
dataset_sp = rawdata_struct.Spinel;

validateRequiredVariables(dataset_ol, dataset_sp);

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Display these notes only once per function call.
pressureNoteIssued = false;
applicationCautionIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop + calculation
disp('=== Step 3: Selecting a data code from the list (Olivine) ===');

while true
    % ----- Olivine selection -----
    dataCodes_ol = dataset_ol{:, 1};

    [selectedIdx_ol, ok] = listdlg( ...
        'PromptString', 'Please select the Olivine data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_ol, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_ol)
        disp('Selection canceled');
        break;
    end

    selectedCode_ol = dataCodes_ol(selectedIdx_ol);
    disp(['Olivine selected: ' char(string(selectedCode_ol))]);

    % ----- Spinel selection -----
    disp('=== Step 4: Selecting a data code from the list (Spinel) ===');

    dataCodes_sp = dataset_sp{:, 1};

    [selectedIdx_sp, ok] = listdlg( ...
        'PromptString', 'Please select the Spinel data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', dataCodes_sp, ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_sp)
        disp('Selection canceled');
        break;
    end

    selectedCode_sp = dataCodes_sp(selectedIdx_sp);
    disp(['Spinel selected: ' char(string(selectedCode_sp))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_ol = dataset_ol(selectedIdx_ol, :);
    selectedData_sp = dataset_sp(selectedIdx_sp, :);

    nanInputNames = findNaNInputs(selectedData_ol, selectedData_sp);
    validateCompositionInputs(selectedData_ol, selectedData_sp);

    row = calcTemp(selectedData_ol, selectedData_sp, P_kbar);

    % Store selected identifiers for traceability.
    row.dataCode_ol = repmat(string(selectedCode_ol), height(row), 1);
    row.dataCode_sp = repmat(string(selectedCode_sp), height(row), 1);
    row = movevars(row, {'dataCode_ol','dataCode_sp'}, 'Before', 1);

    % Store this result as one table block.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo calculated temperature.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
            ': ' num2str(row.T_deg) ' degreeC']);
    else
        finiteT = row.T_deg(isfinite(row.T_deg));
        if isempty(finiteT)
            temperatureText = 'non-finite result';
        elseif max(finiteT) == min(finiteT)
            temperatureText = [num2str(finiteT(1)) ...
                ' degreeC at all pressure points'];
        else
            temperatureText = [num2str(row.T_deg(1)) ' to ' ...
                num2str(row.T_deg(end)) ' degreeC'];
        end
        disp([char(string(selectedCode_ol)) ' & ' char(string(selectedCode_sp)) ...
            ': ' temperatureText]);
    end

    % Jackson's equation has no pressure term. State this once so that a
    % pressure-range run is not misinterpreted as pressure sensitive.
    if ~pressureNoteIssued
        fprintf(2, ...
            ['NOTE: Jackson (1969) contains no explicit pressure term. ' ...
             'P_kbar is retained for interface compatibility and output only; ' ...
             'temperature is identical at all supplied pressure points for a given pair.\n']);
        pressureNoteIssued = true;
    end

    % Display the principal literature cautions once per function call.
    if ~applicationCautionIssued
        fprintf(2, ...
            ['CAUTION: The Jackson (1969) thermometer has no modern, experimentally ' ...
             'defined calibration range. It may yield temperatures that are too high ' ...
             'for regionally metamorphosed serpentinites and a large, unreasonable ' ...
             'temperature range for ferric-iron-rich Spinel. Use only texturally ' ...
             'equilibrated Olivine-Spinel pairs and independently evaluated Fe2+/Fe3+.\n']);
        applicationCautionIssued = true;
    end

    % NaN input warning.
    if ~isempty(nanInputNames)
        fprintf(2, ['WARNING: NaN was found in the thermometer input(s) for %s & %s: %s.\n' ...
                    '         The calculation was continued, but the calculated temperature may be NaN.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    % Non-finite result warning.
    invalidTemperature = ~isfinite(row.T_deg);
    if any(invalidTemperature)
        fprintf(2, ['WARNING: Non-finite temperature values were calculated for %s & %s ' ...
                    '(%d of %d points; NaN: %d, Inf: %d).\n' ...
                    '         These values remain in the output table, and the calculation has not been stopped.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)), ...
            sum(invalidTemperature), ...
            numel(row.T_deg), ...
            sum(isnan(row.T_deg)), ...
            sum(isinf(row.T_deg)));
    end

    % The printed equation is a quotient. A zero or negative denominator
    % produces an undefined or physically questionable result.
    invalidDenominator = isfinite(row.denominator) & (row.denominator <= 0);
    if any(invalidDenominator)
        fprintf(2, ...
            ['WARNING: The Jackson (1969) temperature denominator is zero or ' ...
             'negative for %s & %s. Such results should not be interpreted.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)));
    end

    % The reproduced relation is explicitly described as degreeC. Warn for
    % values below absolute zero without otherwise imposing a calibration range.
    belowAbsoluteZero = isfinite(row.T_deg) & row.T_deg < -273.15;
    if any(belowAbsoluteZero)
        fprintf(2, ...
            ['WARNING: Calculated Jackson (1969) temperature is below absolute ' ...
             'zero for %s & %s and is non-physical.\n'], ...
            char(string(selectedCode_ol)), ...
            char(string(selectedCode_sp)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Jackson1969', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all table blocks once after selection is complete.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(data_olivine, data_spinel)
% Validate that all variables required by Jackson (1969) are available.

olivineVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
spinelVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', ...
    'Fe3_cation_apfu', 'Cr_cation_apfu', 'Al_cation_apfu'};

missingNames = strings(0, 1);

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    if ~ismember(variableName, data_olivine.Properties.VariableNames)
        missingNames(end + 1, 1) = ...
            "Olivine." + string(variableName); %#ok<AGROW>
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    if ~ismember(variableName, data_spinel.Properties.VariableNames)
        missingNames(end + 1, 1) = ...
            "Spinel." + string(variableName); %#ok<AGROW>
    end
end

if ~isempty(missingNames)
    error(['Jackson1969: required table variable(s) are missing: ' ...
        char(strjoin(missingNames, ', ')) '.']);
end

end

function nanInputNames = findNaNInputs(data_olivine, data_spinel)
% Return names of required input variables containing NaN.

olivineVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
spinelVariables = {'Mg_cation_apfu', 'Fe_cation_apfu', ...
    'Fe3_cation_apfu', 'Cr_cation_apfu', 'Al_cation_apfu'};

nanInputNames = strings(0, 1);

for i = 1:numel(olivineVariables)
    variableName = olivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isnan(variableValue(:)))
        nanInputNames(end + 1, 1) = ...
            "Olivine." + string(variableName); %#ok<AGROW>
    end
end

for i = 1:numel(spinelVariables)
    variableName = spinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isnan(variableValue(:)))
        nanInputNames(end + 1, 1) = ...
            "Spinel." + string(variableName); %#ok<AGROW>
    end
end

end

function validateCompositionInputs(data_olivine, data_spinel)
% Validate sign constraints required by the Jackson calculation.
%
% Mg and Fe2+ must be > 0 because they enter ratios and logarithms.
% Fe3+, Cr, and Al may be zero, but finite negative values are invalid.
% NaN is intentionally allowed to propagate to the output.

positiveOlivineVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
positiveSpinelVariables = {'Mg_cation_apfu', 'Fe_cation_apfu'};
nonnegativeSpinelVariables = {'Fe3_cation_apfu', ...
    'Cr_cation_apfu', 'Al_cation_apfu'};

invalidPositiveNames = strings(0, 1);
invalidNonnegativeNames = strings(0, 1);

for i = 1:numel(positiveOlivineVariables)
    variableName = positiveOlivineVariables{i};
    variableValue = data_olivine.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) <= 0)
        invalidPositiveNames(end + 1, 1) = ...
            "Olivine." + string(variableName); %#ok<AGROW>
    end
end

for i = 1:numel(positiveSpinelVariables)
    variableName = positiveSpinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) <= 0)
        invalidPositiveNames(end + 1, 1) = ...
            "Spinel." + string(variableName); %#ok<AGROW>
    end
end

for i = 1:numel(nonnegativeSpinelVariables)
    variableName = nonnegativeSpinelVariables{i};
    variableValue = data_spinel.(variableName);
    if any(isfinite(variableValue(:)) & variableValue(:) < 0)
        invalidNonnegativeNames(end + 1, 1) = ...
            "Spinel." + string(variableName); %#ok<AGROW>
    end
end

if ~isempty(invalidPositiveNames)
    error(['Jackson1969: Mg and Fe2+ inputs must be > 0. ' ...
        'Zero or negative finite value(s) were found in: ' ...
        char(strjoin(invalidPositiveNames, ', ')) '.']);
end

if ~isempty(invalidNonnegativeNames)
    error(['Jackson1969: Fe3+, Cr, and Al inputs must be >= 0. ' ...
        'Negative finite value(s) were found in: ' ...
        char(strjoin(invalidNonnegativeNames, ', ')) '.']);
end

trivalentSum = data_spinel.Cr_cation_apfu + ...
    data_spinel.Al_cation_apfu + data_spinel.Fe3_cation_apfu;

if any(isfinite(trivalentSum(:)) & trivalentSum(:) <= 0)
    error(['Jackson1969: the Spinel trivalent-cation sum ' ...
        '(Cr + Al + Fe3+) must be > 0.']);
end

end

function row = calcTemp(data_olivine, data_spinel, P_kbar)
% Calculate Jackson (1969) temperature for one Olivine-Spinel pair.
%
% P_kbar is retained to match the common thermometer interface. It does not
% enter the Jackson temperature equation reproduced in Appendix 2 of
% Hartmann and Chemale-Junior (2003).

P_kbar = P_kbar(:);
P_GPa = P_kbar ./ 10;
nP = numel(P_kbar);

row = table();
row.P_kbar = P_kbar;
row.P_GPa = P_GPa;

% --- Extract and expand selected cation data ---
Mg_ol  = repmat(data_olivine.Mg_cation_apfu, nP, 1);
Fe2_ol = repmat(data_olivine.Fe_cation_apfu, nP, 1);

Mg_sp  = repmat(data_spinel.Mg_cation_apfu, nP, 1);
Fe2_sp = repmat(data_spinel.Fe_cation_apfu, nP, 1);
Fe3_sp = repmat(data_spinel.Fe3_cation_apfu, nP, 1);
Cr_sp  = repmat(data_spinel.Cr_cation_apfu, nP, 1);
Al_sp  = repmat(data_spinel.Al_cation_apfu, nP, 1);

% --- Divalent-cation fractions ---
XMg_ol  = Mg_ol  ./ (Mg_ol + Fe2_ol);
XFe2_ol = Fe2_ol ./ (Mg_ol + Fe2_ol);

XMg_sp  = Mg_sp  ./ (Mg_sp + Fe2_sp);
XFe2_sp = Fe2_sp ./ (Mg_sp + Fe2_sp);

% --- Trivalent-cation fractions in Spinel ---
trivalentSum_sp = Cr_sp + Al_sp + Fe3_sp;
a_Cr_sp  = Cr_sp  ./ trivalentSum_sp;
b_Al_sp  = Al_sp  ./ trivalentSum_sp;
c_Fe3_sp = Fe3_sp ./ trivalentSum_sp;
sum_abc = a_Cr_sp + b_Al_sp + c_Fe3_sp;

% --- Fe-Mg exchange coefficient ---
KD_FeMg_ol_sp = ...
    (XMg_ol .* XFe2_sp) ./ (XFe2_ol .* XMg_sp);
lnKD_FeMg_ol_sp = log(KD_FeMg_ol_sp);

% --- Jackson (1969) equation as reproduced by Hartmann and
%     Chemale-Junior (2003, Appendix 2) ---
R_cal = 1.9871;

numerator = 5580 .* a_Cr_sp ...
          + 1018 .* b_Al_sp ...
          - 1720 .* c_Fe3_sp ...
          + 2400;

denominator = 0.9 .* a_Cr_sp ...
            + 2.56 .* b_Al_sp ...
            - 3.08 .* c_Fe3_sp ...
            - 1.47 ...
            + R_cal .* lnKD_FeMg_ol_sp;

% Hartmann and Chemale-Junior (2003) explicitly describe the printed
% quotient as temperature in degreeC. This implementation follows that
% description verbatim.
T_deg = numerator ./ denominator;
T_K = T_deg + 273.15;

% --- Pack outputs ---
row.XMg_ol = XMg_ol;
row.XFe2_ol = XFe2_ol;
row.Fo_ol_molpct = 100 .* XMg_ol;

row.XMg_sp = XMg_sp;
row.XFe2_sp = XFe2_sp;

row.a_Cr_sp = a_Cr_sp;
row.b_Al_sp = b_Al_sp;
row.c_Fe3_sp = c_Fe3_sp;
row.sum_abc = sum_abc;
row.trivalentSum_sp = trivalentSum_sp;

row.KD_FeMg_ol_sp = KD_FeMg_ol_sp;
row.lnKD_FeMg_ol_sp = lnKD_FeMg_ol_sp;
row.R_cal = repmat(R_cal, nP, 1);

row.numerator = numerator;
row.denominator = denominator;
row.T_K = T_K;
row.T_deg = T_deg;

end
