function MWinfo = getMolarWeights()
% liquid.getMolarWeights
% Load molar weights and cation numbers from:
%   functions/+liquid/Cation_moduli.xlsx, sheet 'Molar weight'
%
% Uses ONLY oxides consistent with cationCalc4EPMA:
%   SiO2 TiO2 Al2O3 FeO MnO MgO CaO Na2O K2O V2O3 Cr2O3 NiO P2O5 SO3 F Cl Fe2O3

thisFileDir = fileparts(mfilename('fullpath')); % .../functions/+liquid
xlsxPath = fullfile(thisFileDir, 'Cation_moduli.xlsx');
if ~exist(xlsxPath, 'file')
    error('Cation_moduli.xlsx not found at: %s', xlsxPath);
end

opts = detectImportOptions(xlsxPath, 'Sheet', 'Molar weight', 'NumHeaderLines', 0);
T = readtable(xlsxPath, opts);
T.Properties.VariableNames = opts.VariableNames;

rowLabelCol = T.Properties.VariableNames{1};
labelsS = string(T.(rowLabelCol));

iMW  = find(strcmpi(labelsS, 'Molar weights'), 1, 'first');
iCat = find(strcmpi(labelsS, 'Cation#'), 1, 'first');

if isempty(iMW) || isempty(iCat)
    error('Sheet "Molar weight" must contain row labels: "Molar weights" and "Cation#".');
end

useOx = {'SiO2','TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O', ...
         'V2O3','Cr2O3','NiO','P2O5','SO3','F','Cl','Fe2O3'};

MW = struct();
Cat = struct();

for j = 1:numel(useOx)
    ox = useOx{j};

    if ~any(strcmp(T.Properties.VariableNames, ox))
        error('Cation_moduli.xlsx: column "%s" not found in sheet "Molar weight".', ox);
    end

    mwVal  = T{iMW,  ox};
    catVal = T{iCat, ox};

    if isempty(mwVal) || ~isfinite(mwVal)
        error('Cation_moduli.xlsx: MW for "%s" is missing/invalid.', ox);
    end
    if isempty(catVal) || ~isfinite(catVal)
        error('Cation_moduli.xlsx: Cation# for "%s" is missing/invalid.', ox);
    end

    MW.(ox)  = mwVal;
    Cat.(ox) = catVal;
end

MWinfo = struct();
MWinfo.file = xlsxPath;
MWinfo.MW   = MW;
MWinfo.Cat  = Cat;
MWinfo.useOxides = useOx;

end
