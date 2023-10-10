%% Load initial state from data for the large network
% Input: Inifile: Filename, string
%        N_E/I:  Number of neurons used here. 
%        EEDly, IEDly, EIDly: logical, if these delay mats will be used
% Output: InE/I: Adjusted initial conditions for E/I cells
%         DlyMat: Delay Mats used 

function [InE, InI] = LargeNW_LoadIniState(InifileName, Inifile,N_E, N_I)
if exist(InifileName) == 2
InSs = load(InifileName);
else
    InSs = Inifile;
end
if ~isfield(InSs,'VE') || ~isfield(InSs,'VI')
    disp('Not proper ini state file')
    return
end
N_Eori = length(InSs.VE); N_Iori = length(InSs.VI); 
%% get all field names
FNames = fieldnames(InSs);
DlyNames = FNames(contains(FNames,'Dly'));
for FInd = 1:length(DlyNames)
    DlyMat.(DlyNames{FInd}) = InSs.(DlyNames{FInd});
end
InSu = rmfield(InSs,DlyNames); % remove dly fields 

% getE and I fields
ENames = FNames(contains(FNames,'E'));
for FInd = 1:length(ENames)
     CurrentField = InSu.(ENames{FInd});
     InE.(ENames{FInd}) = zeros(N_E,1);
     if N_Eori <= N_E
         InE.(ENames{FInd})(1:N_Eori) = CurrentField;
     else
         InE.(ENames{FInd}) = CurrentField(1:N_E);
     end
end

INames = FNames(contains(FNames,'I'));
for FInd = 1:length(INames)
     CurrentField = InSu.(INames{FInd});
     InI.(INames{FInd}) = zeros(N_I,1);
     if N_Iori <= N_I
         InI.(INames{FInd})(1:N_Iori) = CurrentField;
     else
         InI.(INames{FInd}) = CurrentField(1:N_I);
     end
end

end
