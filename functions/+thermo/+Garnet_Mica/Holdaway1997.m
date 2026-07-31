function results = Holdaway1997(rawdata_struct, P_kbar)
% functions/+thermo/+Garnet_Mica/Holdaway1997.m
% Tested with MATLAB R2024b
%
% Garnet-biotite Fe-Mg exchange thermometer (simplified implementation)
% Holdaway, M.J., Mukhopadhyay, B., Dyar, M.D., Guidotti, C.V. and
% Dutrow, B.L. (1997)
% American Mineralogist, 82, 582-595
% DOI: https://doi.org/10.2138/am-1997-5-618
%
% -------------------------------------------------------------------------
% OVERVIEW
% This function interactively pairs one Garnet analysis and one Mica
% analysis and calculates temperature using the Fe-Mg non-ideal core of the
% Holdaway et al. (1997) garnet-biotite thermometer.
%
% IMPORTANT: this legacy implementation is not the complete multicomponent
% Holdaway et al. (1997) calibration. It does not explicitly implement the
% garnet Ca-Mn or biotite Al-Ti terms, and stored Fe3+ values are not used in
% the calculation. The approximately +/-25 degreeC precision discussed by
% Holdaway et al. (1997) for optimum use of their complete calibration must
% therefore not be assigned automatically to results from this simplified
% code (p. 594).
%
% The function is designed for repeated calculations. For every selected
% Garnet-Mica pair, it returns one output row per supplied pressure value.
% P_kbar can therefore be either a scalar (fixed-pressure calculation) or a
% vector (pressure-range calculation), following the Ballhaus1991 interface.
%
% -------------------------------------------------------------------------
% EXPERIMENTAL CALIBRATION RANGE AND APPLICATION NOTES
%
% Holdaway et al. (1997) combined two principal experimental datasets:
%
%   Ferry-Spear experiments      : 550-800 degreeC at 2.07 kbar
%                                  (pp. 584-585)
%   Perchuk-Lavrent'eva data     : 575-950 degreeC at approximately 6 kbar
%                                  (pp. 586-587)
%
% These are discrete experimental pressure levels, not a continuously
% calibrated 2-6 kbar pressure series. The Maine natural-specimen dataset
% chiefly tests approximately 550-650 degreeC and 3-4.5 kbar (p. 582), and
% the Ontario granulite comparison tests approximately 650-700 degreeC and
% 4-6 kbar (p. 594). Thus, the best-tested natural-rock interval is much
% narrower than the full temperature envelope of the experiments.
%
% This implementation uses the conservative envelope 550-950 degreeC and
% 2-6 kbar only as a screening range. It prints non-stopping warnings when:
%   1) input pressure is outside 2-6 kbar, or
%   2) a finite calculated temperature is outside 550-950 degreeC.
% A result inside this envelope is not, by itself, proof that the mineral
% compositions, oxidation state, or equilibrium assumptions are valid.
%
% Additional cautions from Holdaway et al. (1997):
% - Fe3+ must be measured or estimated. Ignoring it can cause significant
%   error; biotite Fe3+ is especially important. Very reduced, graphite-
%   bearing rocks containing hematite-free ilmenite are preferred
%   (pp. 583-584). This simplified code does not apply an Fe3+ correction.
% - The Perchuk-Lavrent'eva products have limitations in Al, Ca, Mn, and Ti
%   characterization, and their equilibrium status requires care
%   (pp. 586-587).
% - Use the highest Mg/Fe portion of a garnet traverse as the best estimate
%   of peak conditions, avoid retrograde rims, and compare multiple nearby
%   biotites. Retrograde exchange can lower high-grade estimates by as much
%   as about 25 degreeC (pp. 591-592).
% - Independent partial-melting experiments at 825-975 degreeC and 7-13
%   kbar were not reproduced: the thermometer yielded only 666-725 degreeC.
%   The authors discussed Fe loss and disequilibrium as possible causes
%   (p. 592). High-temperature/high-pressure use is therefore poorly tested.
% - The quoted approximately +/-25 degreeC performance applies only under
%   optimum analytical, ferric-iron, compositional, and equilibrium
%   conditions; simple error propagation is closer to 50 degreeC (p. 594).
%
% A later revision (Holdaway, 2000, American Mineralogist 85, 881-892;
% DOI: https://doi.org/10.2138/am-2000-0701) recommends the internally
% consistent 5AV model and notes an error in the garnet activity equation
% of the 1997 computer program (pp. 881 and 885).
% Do not mix parameters among different Holdaway model versions.
%
% -------------------------------------------------------------------------
% REQUIRED INPUT FORMAT (IMPORTANT)
% rawdata_struct must contain:
%   rawdata_struct.Garnet : table
%   rawdata_struct.Mica   : table
%
% The FIRST column of each table is treated as an identifier ("data code")
% displayed in the selection dialog. The remaining columns should contain
% normalized cation data.
%
% Minimum required variables used by this simplified calculation:
%   Fe_cation_apfu
%   Mg_cation_apfu
%
% Optional variables retained in the output when available:
%   Fe3_cation_apfu
%   Mn_cation_apfu
%   Ca_cation_apfu
%   Ti_cation_apfu
%   Al_cation_apfu
%   Si_cation_apfu
%   K_cation_apfu
%   Na_cation_apfu
%
% Finite cation values must be non-negative. Negative values and Inf are
% rejected. NaN is retained as missing data and propagated; it is never
% replaced by zero. Missing optional variables are also represented by NaN.
% A zero or NaN in any Fe-Mg value required by the exchange coefficient
% makes the temperature undefined, so the output temperature remains NaN
% and a non-stopping warning is printed. In particular, NaN cannot become
% 0 K or -273.15 degreeC in this implementation.
%
% -------------------------------------------------------------------------
% THERMOMETER FORMULATION (as implemented here)
%
% Exchange coefficient:
%   K_D = (Mg/Fe2)_Grt / (Mg/Fe2)_Bt
%
% Temperature equation:
%   T(K) = (41952 + 0.311*P(bar) + G + B) / (10.35 - 3*R*ln(K_D))
%
% Garnet Fe-Mg Margules parameters:
%   W_FeMg_Grt = -24166 + 22.09*T - 0.034*P(bar)
%   W_MgFe_Grt =  22265 - 12.40*T + 0.050*P(bar)
%
% Garnet correction:
%   G = 2*XMg*XFe*(W_FeMg-W_MgFe)
%       + XFe^2*W_MgFe - XMg^2*W_FeMg
%
% Biotite Fe-Mg parameter and correction:
%   W_MgFe_Bt = 40719 - 30*T
%   B = W_MgFe_Bt*(XMg_Bt-XFe_Bt)
%
% The equation is solved iteratively because G and B depend on temperature.
%
% -------------------------------------------------------------------------
% Syntax:
%   results = Holdaway1997(rawdata_struct, P_kbar)
%
% Inputs:
%   rawdata_struct : struct containing Garnet and Mica tables (see above)
%   P_kbar         : pressure in kbar (finite, non-negative scalar or vector)
%
% Output:
%   results : table containing one row per pressure value for every
%             user-selected Garnet-Mica pair

%% Input validation
% Basic argument checks prevent silent failures due to missing inputs or
% invalid pressure values.
if nargin < 2
    error('Holdaway1997 requires (rawdata_struct, P_kbar).');
end
if ~isstruct(rawdata_struct)
    error('rawdata_struct must be a struct.');
end
if ~isnumeric(P_kbar) || ~isreal(P_kbar) || isempty(P_kbar) || ...
        ~isvector(P_kbar) || any(~isfinite(P_kbar(:))) || ...
        any(P_kbar(:) < 0)
    error('P_kbar must be a finite non-negative numeric scalar or vector.');
end

P_kbar = P_kbar(:);

%% 1) Retrieve cation datasets
% Extract the required tables from the input struct. The source tables are
% read only; selected rows are passed to the calculation local function.
disp('=== Step 1: Preparing cation dataset ===');

if ~isfield(rawdata_struct, 'Garnet') || ~istable(rawdata_struct.Garnet)
    error('rawdata_struct must contain table: rawdata_struct.Garnet');
end
if ~isfield(rawdata_struct, 'Mica') || ~istable(rawdata_struct.Mica)
    error('rawdata_struct must contain table: rawdata_struct.Mica');
end

dataset_grt = rawdata_struct.Garnet;
dataset_mica = rawdata_struct.Mica;

if width(dataset_grt) < 1 || height(dataset_grt) < 1
    error('rawdata_struct.Garnet must be a non-empty table.');
end
if width(dataset_mica) < 1 || height(dataset_mica) < 1
    error('rawdata_struct.Mica must be a non-empty table.');
end

validateRequiredVariables(dataset_grt, 'Garnet');
validateRequiredVariables(dataset_mica, 'Mica');

disp('=== Preparing cation dataset has been finished ===');

%% 2) Initialize output container
% Result blocks are buffered so the results table is not resized on every
% loop iteration. The buffer doubles only when its capacity is exhausted,
% and all blocks are concatenated once after the interactive loop.
disp('=== Step 2: Preparing output container ===');

initialBufferCapacity = 16;
resultBlocks = cell(initialBufferCapacity, 1);
nResultBlocks = 0;

% Conservative screening limits spanning the principal experiments.
screeningT_min_degC = 550;
screeningT_max_degC = 950;
screeningP_min_kbar = 2;
screeningP_max_kbar = 6;

pressureOutsideScreening = ...
    P_kbar < screeningP_min_kbar | P_kbar > screeningP_max_kbar;
pressureWarningIssued = false;

disp('=== Preparing output container has been finished ===');

%% 3-5) Interactive selection loop and calculation
% The loop continues until the user cancels a selection dialog or chooses
% Finish after a calculation.
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

    % ----- Mica selection -----
    disp('=== Step 4: Selecting a data code from the list (Mica) ===');
    dataCodes_mica = dataset_mica{:, 1};

    [selectedIdx_mica, ok] = listdlg( ...
        'PromptString', 'Please select the Mica data you would like to use:', ...
        'SelectionMode', 'single', ...
        'ListString', cellstr(string(dataCodes_mica)), ...
        'ListSize', [320 320]);

    if ~ok || isempty(selectedIdx_mica)
        disp('Selection canceled');
        break;
    end

    selectedCode_mica = dataCodes_mica(selectedIdx_mica);
    disp(['Mica selected: ' char(string(selectedCode_mica))]);

    % ----- Calculation -----
    disp('=== Step 5: Calculating the temperature ===');

    selectedData_grt = dataset_grt(selectedIdx_grt, :);
    selectedData_mica = dataset_mica(selectedIdx_mica, :);

    % Negative values and Inf stop the calculation. NaN and zero are kept
    % so that an undefined exchange calculation is returned as NaN and can
    % be reported with non-stopping fprintf messages below.
    validateSelectedInputs(selectedData_grt, selectedData_mica);
    [nanInputNames, zeroInputNames] = ...
        findSpecialExchangeInputs(selectedData_grt, selectedData_mica);

    row = calcTemp(selectedData_grt, selectedData_mica, P_kbar);

    % Store identifiers once per output pressure point.
    row.dataCode_grt = repmat(string(selectedCode_grt), height(row), 1);
    row.dataCode_mica = repmat(string(selectedCode_mica), height(row), 1);
    row = movevars(row, {'dataCode_grt', 'dataCode_mica'}, 'Before', 1);

    % Store this result as one table block. Buffer growth occurs only when
    % capacity is exhausted, not on every interactive-loop iteration.
    nResultBlocks = nResultBlocks + 1;
    if nResultBlocks > numel(resultBlocks)
        resultBlocks = [resultBlocks; cell(numel(resultBlocks), 1)]; %#ok<AGROW>
    end
    resultBlocks{nResultBlocks} = row;

    % Echo the calculated scalar value or pressure-series endpoint values.
    disp('--------------------------------------------------');
    disp('=== Temperature was calculated: ===');
    if height(row) == 1
        disp([char(string(selectedCode_grt)) ' & ' ...
            char(string(selectedCode_mica)) ': Holdaway1997 = ' ...
            num2str(row.T_H97_C) ' degreeC']);
    else
        disp([char(string(selectedCode_grt)) ' & ' ...
            char(string(selectedCode_mica)) ': Holdaway1997 = ' ...
            num2str(row.T_H97_C(1)) ' to ' ...
            num2str(row.T_H97_C(end)) ' degreeC']);
    end

    % Pressure warnings are issued only once because pressure is common to
    % all selected mineral pairs in this function call.
    if any(pressureOutsideScreening) && ~pressureWarningIssued
        fprintf(2, ...
            ['WARNING: Input pressure is outside the conservative screening ' ...
             'range used for Holdaway et al. (1997): 2-6 kbar. ' ...
             '%d of %d pressure point(s) are outside the range; ' ...
             'input range = %.4g-%.4g kbar. The principal experiments were ' ...
             'conducted at discrete pressures near 2.07 and 6 kbar.\n'], ...
            sum(pressureOutsideScreening), numel(P_kbar), ...
            min(P_kbar), max(P_kbar));
        pressureWarningIssued = true;
    end

    finiteTemperature = isfinite(row.T_H97_C);
    temperatureOutsideScreening = finiteTemperature & ...
        (row.T_H97_C < screeningT_min_degC | ...
         row.T_H97_C > screeningT_max_degC);

    if any(temperatureOutsideScreening)
        finiteValues = row.T_H97_C(finiteTemperature);
        fprintf(2, ...
            ['WARNING: Calculated temperature is outside the conservative ' ...
             'screening range used for Holdaway et al. (1997): ' ...
             '550-950 degreeC. %d of %d finite point(s) are outside; ' ...
             'finite calculated range = %.4g-%.4g degreeC for %s & %s.\n'], ...
            sum(temperatureOutsideScreening), sum(finiteTemperature), ...
            min(finiteValues), max(finiteValues), ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)));
    end

    if ~isempty(nanInputNames)
        fprintf(2, ...
            ['WARNING: NaN was found in an Fe-Mg thermometer input for ' ...
             '%s & %s: %s. NaN was retained; it was not replaced by zero.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)), ...
            char(strjoin(nanInputNames, ', ')));
    end

    if ~isempty(zeroInputNames)
        fprintf(2, ...
            ['WARNING: Zero was found in an Fe-Mg thermometer input for ' ...
             '%s & %s: %s. A positive Fe2 and Mg value is required to ' ...
             'define K_D, so the corresponding temperature remains NaN.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)), ...
            char(strjoin(zeroInputNames, ', ')));
    end

    invalidTemperature = ~isfinite(row.T_H97_C);
    if any(invalidTemperature)
        fprintf(2, ...
            ['WARNING: Non-finite temperature values were calculated for ' ...
             '%s & %s (%d of %d points; NaN: %d, Inf: %d). ' ...
             'They remain in the output table and calculation continues.\n'], ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)), ...
            sum(invalidTemperature), numel(row.T_H97_C), ...
            sum(isnan(row.T_H97_C)), sum(isinf(row.T_H97_C)));
    end

    finiteNotConverged = finiteTemperature & ~row.isConverged;
    if any(finiteNotConverged)
        fprintf(2, ...
            ['WARNING: The iterative Holdaway1997 solution did not meet ' ...
             'the convergence tolerance for %d of %d finite point(s) for ' ...
             '%s & %s. The final finite estimates were retained.\n'], ...
            sum(finiteNotConverged), sum(finiteTemperature), ...
            char(string(selectedCode_grt)), char(string(selectedCode_mica)));
    end

    disp('--------------------------------------------------');

    userAction = questdlg( ...
        'Continue with another data selection ?', ...
        'Holdaway1997', ...
        'Continue', 'Finish', 'Continue');

    if isempty(userAction) || strcmp(userAction, 'Finish')
        disp('Proceed to save data');
        break;
    else
        disp('Repeat again');
    end
end

% Concatenate all buffered blocks once. Return an empty table if no pair was
% calculated.
if nResultBlocks == 0
    results = table();
else
    results = vertcat(resultBlocks{1:nResultBlocks});
end

disp('=== Calculation has been finished! ===');

end

%% ---- local functions ----
function validateRequiredVariables(dataset, mineralLabel)
% validateRequiredVariables
% Confirm that the columns used by the simplified Fe-Mg equation exist.

requiredVariables = {'Fe_cation_apfu', 'Mg_cation_apfu'};
missingVariables = requiredVariables(~ismember( ...
    requiredVariables, dataset.Properties.VariableNames));

if ~isempty(missingVariables)
    error('%s table must contain variable(s): %s', mineralLabel, ...
        char(strjoin(string(missingVariables), ', ')));
end

end

function validateSelectedInputs(data_grt, data_mica)
% validateSelectedInputs
% Reject negative, infinite, complex, nonnumeric, or nonscalar cation data.
% NaN and zero are intentionally accepted here and handled as undefined
% exchange inputs by calcTemp and the caller's fprintf warnings.

if height(data_grt) ~= 1 || height(data_mica) ~= 1
    error('Selected Garnet and Mica inputs must each be a 1-row table.');
end

variableNames = {'Fe_cation_apfu', 'Fe3_cation_apfu', ...
    'Mg_cation_apfu', 'Mn_cation_apfu', 'Ca_cation_apfu', ...
    'Ti_cation_apfu', 'Al_cation_apfu', 'Si_cation_apfu', ...
    'K_cation_apfu', 'Na_cation_apfu'};

mineralTables = {data_grt, data_mica};
mineralLabels = {'Garnet', 'Mica'};

for mineralIndex = 1:numel(mineralTables)
    selectedTable = mineralTables{mineralIndex};
    for variableIndex = 1:numel(variableNames)
        variableName = variableNames{variableIndex};
        if ismember(variableName, selectedTable.Properties.VariableNames)
            value = selectedTable.(variableName);
            if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
                    isinf(value) || value < 0
                error(['%s.%s must be a real scalar that is non-negative ' ...
                    'or NaN; negative values and Inf are not allowed.'], ...
                    mineralLabels{mineralIndex}, variableName);
            end
        end
    end
end

end

function [nanInputNames, zeroInputNames] = ...
        findSpecialExchangeInputs(data_grt, data_mica)
% findSpecialExchangeInputs
% Identify NaN and zero only in the four values actually used by the
% simplified Fe-Mg exchange calculation. Fixed-size arrays avoid growth.

inputNames = ["Garnet.Fe_cation_apfu"; "Garnet.Mg_cation_apfu"; ...
    "Mica.Fe_cation_apfu"; "Mica.Mg_cation_apfu"];
inputValues = [data_grt.Fe_cation_apfu; data_grt.Mg_cation_apfu; ...
    data_mica.Fe_cation_apfu; data_mica.Mg_cation_apfu];

nanInputNames = inputNames(isnan(inputValues));
zeroInputNames = inputNames(isfinite(inputValues) & inputValues == 0);

end

function row = calcTemp(data_grt, data_mica, P_kbar)
% calcTemp
% Compute Holdaway et al. (1997) simplified Fe-Mg temperatures for one
% Garnet-Mica pair at every supplied pressure.

P_kbar = P_kbar(:);
P_bar = P_kbar .* 1000;
nP = numel(P_kbar);

R_J = 8.31441;
row = table();
row.P_kbar = P_kbar;
row.P_bar = P_bar;

% Extract one-row cation data. Missing optional variables become NaN and
% supplied NaN values remain NaN.
grt = prepareMineralRow(data_grt, 'Garnet');
mica = prepareMineralRow(data_mica, 'Mica');

% Store mineral data at every pressure point.
row.Fe2_grt = repmat(grt.Fe2, nP, 1);
row.Fe3_grt = repmat(grt.Fe3, nP, 1);
row.Mg_grt = repmat(grt.Mg, nP, 1);
row.Mn_grt = repmat(grt.Mn, nP, 1);
row.Ca_grt = repmat(grt.Ca, nP, 1);
row.Al_grt = repmat(grt.Al, nP, 1);
row.Si_grt = repmat(grt.Si, nP, 1);

row.Fe2_mica = repmat(mica.Fe2, nP, 1);
row.Fe3_mica = repmat(mica.Fe3, nP, 1);
row.Mg_mica = repmat(mica.Mg, nP, 1);
row.Mn_mica = repmat(mica.Mn, nP, 1);
row.Ca_mica = repmat(mica.Ca, nP, 1);
row.Ti_mica = repmat(mica.Ti, nP, 1);
row.Al_mica = repmat(mica.Al, nP, 1);
row.Si_mica = repmat(mica.Si, nP, 1);
row.K_mica = repmat(mica.K, nP, 1);
row.Na_mica = repmat(mica.Na, nP, 1);

% Preallocate every pressure-dependent result so no variable changes size
% during the iterative pressure loop.
XFe_grt = NaN(nP, 1);
XMg_grt = NaN(nP, 1);
XFe_mica = NaN(nP, 1);
XMg_mica = NaN(nP, 1);
MgFe_grt = NaN(nP, 1);
MgFe_mica = NaN(nP, 1);
K_D = NaN(nP, 1);
lnK_D = NaN(nP, 1);

W_FeMg_grt_Jmol = NaN(nP, 1);
W_MgFe_grt_Jmol = NaN(nP, 1);
W_MgFe_mica_Jmol = NaN(nP, 1);
G_Jmol = NaN(nP, 1);
B_Jmol = NaN(nP, 1);
T_H97_K = NaN(nP, 1);
T_H97_C = NaN(nP, 1);
nIterations = NaN(nP, 1);
isConverged = false(nP, 1);

exchangeValues = [grt.Fe2, grt.Mg, mica.Fe2, mica.Mg];
exchangeInputsUsable = all(isfinite(exchangeValues) & exchangeValues > 0);

% Fractions can be retained for finite non-negative inputs even if a zero
% later makes the exchange coefficient undefined.
if all(isfinite(exchangeValues))
    sumFeMg_grt = grt.Fe2 + grt.Mg;
    sumFeMg_mica = mica.Fe2 + mica.Mg;
    if sumFeMg_grt > 0
        XFe_grt(:) = grt.Fe2 ./ sumFeMg_grt;
        XMg_grt(:) = grt.Mg ./ sumFeMg_grt;
    end
    if sumFeMg_mica > 0
        XFe_mica(:) = mica.Fe2 ./ sumFeMg_mica;
        XMg_mica(:) = mica.Mg ./ sumFeMg_mica;
    end
end

% A positive value is required for every ratio term. If NaN or zero occurs,
% all temperature values stay preallocated as NaN; no zero-K placeholder is
% ever introduced.
if exchangeInputsUsable
    MgFe_grt(:) = grt.Mg ./ grt.Fe2;
    MgFe_mica(:) = mica.Mg ./ mica.Fe2;
    K_D(:) = MgFe_grt(1) ./ MgFe_mica(1);

    if isfinite(K_D(1)) && K_D(1) > 0
        lnK_D(:) = log(K_D(1));

        initialT_K = 600 + 273.15;
        tolerance_K = 0.02;
        maxIterations = 200;

        for pressureIndex = 1:nP
            currentT_K = initialT_K;
            lastCandidate_K = NaN;
            iterationsUsed = 0;

            for iterationIndex = 1:maxIterations
                currentParams = calcHoldaway1997Params( ...
                    XFe_grt(pressureIndex), XMg_grt(pressureIndex), ...
                    XFe_mica(pressureIndex), XMg_mica(pressureIndex), ...
                    currentT_K, P_bar(pressureIndex));

                denominator = 10.35 - 3 .* R_J .* lnK_D(pressureIndex);
                numerator = 41952 + 0.311 .* P_bar(pressureIndex) + ...
                    currentParams.G_Jmol + currentParams.B_Jmol;

                if ~isfinite(denominator) || denominator <= 0 || ...
                        ~isfinite(numerator)
                    lastCandidate_K = NaN;
                    break;
                end

                candidateT_K = numerator ./ denominator;
                if ~isfinite(candidateT_K) || candidateT_K <= 0
                    lastCandidate_K = NaN;
                    break;
                end

                iterationsUsed = iterationIndex;
                lastCandidate_K = candidateT_K;

                if abs(candidateT_K - currentT_K) < tolerance_K
                    isConverged(pressureIndex) = true;
                    break;
                end

                currentT_K = 0.5 .* (currentT_K + candidateT_K);
            end

            if isfinite(lastCandidate_K)
                T_H97_K(pressureIndex) = lastCandidate_K;
                T_H97_C(pressureIndex) = lastCandidate_K - 273.15;
                nIterations(pressureIndex) = iterationsUsed;

                finalParams = calcHoldaway1997Params( ...
                    XFe_grt(pressureIndex), XMg_grt(pressureIndex), ...
                    XFe_mica(pressureIndex), XMg_mica(pressureIndex), ...
                    lastCandidate_K, P_bar(pressureIndex));

                W_FeMg_grt_Jmol(pressureIndex) = ...
                    finalParams.W_FeMg_grt_Jmol;
                W_MgFe_grt_Jmol(pressureIndex) = ...
                    finalParams.W_MgFe_grt_Jmol;
                W_MgFe_mica_Jmol(pressureIndex) = ...
                    finalParams.W_MgFe_mica_Jmol;
                G_Jmol(pressureIndex) = finalParams.G_Jmol;
                B_Jmol(pressureIndex) = finalParams.B_Jmol;
            end
        end
    end
end

% Pack composition indices, thermodynamic terms, and temperatures.
row.XFe_grt = XFe_grt;
row.XMg_grt = XMg_grt;
row.XFe_mica = XFe_mica;
row.XMg_mica = XMg_mica;
row.MgFe_grt = MgFe_grt;
row.MgFe_mica = MgFe_mica;
row.K_D = K_D;
row.lnK_D = lnK_D;
row.W_FeMg_grt_Jmol = W_FeMg_grt_Jmol;
row.W_MgFe_grt_Jmol = W_MgFe_grt_Jmol;
row.W_MgFe_mica_Jmol = W_MgFe_mica_Jmol;
row.G_Jmol = G_Jmol;
row.B_Jmol = B_Jmol;
row.T_H97_K = T_H97_K;
row.T_H97_C = T_H97_C;
row.nIterations = nIterations;
row.isConverged = isConverged;

end

function params = calcHoldaway1997Params( ...
        XFe_grt, XMg_grt, XFe_mica, XMg_mica, T_K, P_bar)
% calcHoldaway1997Params
% Evaluate the simplified Holdaway et al. (1997) Fe-Mg non-ideal terms.

params = struct();

W_FeMg_grt_Jmol = -24166 + 22.09 .* T_K - 0.034 .* P_bar;
W_MgFe_grt_Jmol = 22265 - 12.40 .* T_K + 0.050 .* P_bar;

G_Jmol = 2 .* XMg_grt .* XFe_grt .* ...
    (W_FeMg_grt_Jmol - W_MgFe_grt_Jmol) + ...
    (XFe_grt .^ 2) .* W_MgFe_grt_Jmol - ...
    (XMg_grt .^ 2) .* W_FeMg_grt_Jmol;

W_MgFe_mica_Jmol = 40719 - 30 .* T_K;
B_Jmol = W_MgFe_mica_Jmol .* (XMg_mica - XFe_mica);

params.W_FeMg_grt_Jmol = W_FeMg_grt_Jmol;
params.W_MgFe_grt_Jmol = W_MgFe_grt_Jmol;
params.W_MgFe_mica_Jmol = W_MgFe_mica_Jmol;
params.G_Jmol = G_Jmol;
params.B_Jmol = B_Jmol;

end

function mineral = prepareMineralRow(data_tbl, mineralLabel)
% prepareMineralRow
% Extract one mineral analysis without replacing NaN or missing optional
% variables by zero.

if height(data_tbl) ~= 1
    error('%s input must be a 1-row table.', mineralLabel);
end

mineral = struct();
mineral.Fe2 = getRequiredVariable(data_tbl, 'Fe_cation_apfu', mineralLabel);
mineral.Mg = getRequiredVariable(data_tbl, 'Mg_cation_apfu', mineralLabel);

mineral.Fe3 = getOptionalVariable(data_tbl, 'Fe3_cation_apfu', mineralLabel);
mineral.Mn = getOptionalVariable(data_tbl, 'Mn_cation_apfu', mineralLabel);
mineral.Ca = getOptionalVariable(data_tbl, 'Ca_cation_apfu', mineralLabel);
mineral.Ti = getOptionalVariable(data_tbl, 'Ti_cation_apfu', mineralLabel);
mineral.Al = getOptionalVariable(data_tbl, 'Al_cation_apfu', mineralLabel);
mineral.Si = getOptionalVariable(data_tbl, 'Si_cation_apfu', mineralLabel);
mineral.K = getOptionalVariable(data_tbl, 'K_cation_apfu', mineralLabel);
mineral.Na = getOptionalVariable(data_tbl, 'Na_cation_apfu', mineralLabel);

end

function value = getRequiredVariable(tbl, variableName, mineralLabel)
% getRequiredVariable
% Read a required scalar cation value while preserving NaN.

if ~ismember(variableName, tbl.Properties.VariableNames)
    error('%s table must contain variable: %s', mineralLabel, variableName);
end

value = validateCationValue(tbl.(variableName), mineralLabel, variableName);

end

function value = getOptionalVariable(tbl, variableName, mineralLabel)
% getOptionalVariable
% Read an optional scalar cation value. Missing variables are recorded as
% NaN because absence of a measurement must not be interpreted as zero.

if ismember(variableName, tbl.Properties.VariableNames)
    value = validateCationValue(tbl.(variableName), mineralLabel, variableName);
else
    value = NaN;
end

end

function value = validateCationValue(value, mineralLabel, variableName)
% validateCationValue
% Accept non-negative finite scalars or NaN; reject negative values and Inf.

if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        isinf(value) || value < 0
    error(['%s.%s must be a real scalar that is non-negative or NaN; ' ...
        'negative values and Inf are not allowed.'], mineralLabel, variableName);
end

end
