%% HPC for large amount of IC tests
% ICId Initial Condition Id : 1-4: 0, 22.5, 45, 90 deg
% ContrastID: A. bg fr
%             B. "1/2 contrast" = bg + 1/2 * (X_E-bg, X_I-bg)
%             C. (normal) full contrast = (X_E, X_I)
%             D. super-charged = bg + 3/2 * (X_E-bg, X_I-bg)
% SaveID Just for saving tag

function [] = HPCGlobConvTest(ICId, ContrastID, ITestSize, PertSize)
CurrentFolder = pwd
%FigurePath = [CurrentFolder '/Figures/Demo022722/'];
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/GlobConv/']; % V1D2
%SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/']; % V1
addpath(SaveFolder)
% DataFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/16Function_Scheme/']; % V1D2
% addpath(DataFolder)
addpath([CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/'])

load('GlobalConvTestWS1.mat')
clear LDEFrfuncLarge LDEFrfunc
load('ICTesfuncLarge.mat')
load('ICTesfuncSmall.mat')

CurrentFolder = pwd; % Since loaded data will overwrite pwd...
SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/GlobConv/']; % V1D2
%SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/']; % V1
addpath(SaveFolder)
%% 4.1 Test a large sample
%% first prepare a driven IC
NSample = PertSize*ITestSize;
load('LDETraces_ang0.0.mat')
LDEequv = LDEEpoOutAll{1}{end};
switch ICId
    case 1
        load('LDETraces_ang0.0.mat')
        LDEICPart = LDEEpoOutAll{1}{end};
    case 2
        load('LDETraces_ang22.5.mat')
        LDEICPart = LDEEpoOutAll{1}{end};
    case 3 % rotate 90 deg: 45 angle OD
        load('LDETraces_ang0.0.mat')
        LDEICPart = LDEEpoOutAll{1}{end};
        LDEICPart = ICRot(LDEICPart,3,N_HCOut,NPixX,NPixY);
    case 4 % rotate 90 deg: 45 angle OD
        load('LDETraces_ang0.0.mat')
        LDEICPart = LDEEpoOutAll{1}{end};
        LDEICPart = ICRot(LDEICPart,2,N_HCOut,NPixX,NPixY);
end

% Then a BG IC
LDEbgIC = LDEICPart;
LDEbgIC.S(:) = 2;
LDEbgIC.C(:) = 6;
LDEbgIC.I(:) = 12;

% Compose equv and BG to new testIC
ICWeights = [1,0; 
             1, 1/2; 
             0, 1; 
             1, 3/2];
fields = fieldnames(LDEICPart);
for FInd = 1:length(fields)
    LDEIC_temp.(fields{FInd}) = ...
        ICWeights(ContrastID,1) * LDEbgIC.(fields{FInd}) + ...
        ICWeights(ContrastID,2) * LDEICPart.(fields{FInd});          
end

% Now perturbe and collect a group
IPertList = linspace(0.80,1.20,ITestSize);
PixPertScl = 0.15;
ICTestAll = cell(ITestSize,PertSize);
for ITestInd = 1:ITestSize
    for PertInd = 1:PertSize
        LDEIC_pert = LDEIC_temp;
        LDEIC_pert.I = LDEIC_pert.I * IPertList(ITestInd);
        
        % take a pixelwise peerturbation vec
        PertVec = 1 + rand(size(LDEIC_pert.I)) * PixPertScl*2 - PixPertScl;
        for FInd = 1:length(fields)
            LDEIC_pert.(fields{FInd}) = LDEIC_pert.(fields{FInd}).*PertVec;
        end
        ICTestAll{ITestInd,PertInd} = LDEIC_pert;
    end
end
ICTestAll = reshape(ICTestAll',NSample,1);

%% Simulate for each IC
LDEPertOutAll = cell(NSample,1);
EpocTest = 200;
L2Diff_OneStepAll = zeros(NSample,EpocTest+1);
DiffVecAll = zeros(NSample,EpocTest+1,NPixY*NPixX*3);

tic
for TestInd = 1:NSample
    ICUse = ICTestAll{TestInd};
    IniTest.S = symmHCs(ICUse.S,N_HCOut,NPixX,NPixY);
    IniTest.C = symmHCs(ICUse.C,N_HCOut,NPixX,NPixY);
    IniTest.I = symmHCs(ICUse.I,N_HCOut,NPixX,NPixY);
    
    % for <10 iterations, use large domain, otherwise use small

    [LargeIter,~,~,~,~] = ...
        LDEIteration_16FuncMain(...
        PixInptCtgrUse,IniTest,p,20,...
        C_SS_mean,C_CS_mean,C_IS_mean,...
        C_SC_mean,C_CC_mean,C_IC_mean,...
        C_SI_mean,C_CI_mean,C_II_mean,...
        L4EmeshXLarge,L4ImeshYLarge,LDEFrfuncLarge,...
        N_HCOut,NPixX,NPixY)  ;
    
    [SmallIter,~,~,~,~] = ...
        LDEIteration_16FuncMain(...
        PixInptCtgrUse,LargeIter{end},p,EpocTest,...
        C_SS_mean,C_CS_mean,C_IS_mean,...
        C_SC_mean,C_CC_mean,C_IC_mean,...
        C_SI_mean,C_CI_mean,C_II_mean,...
        L4EmeshXSmall,L4ImeshYSmall,LDEFrfuncSmall,...
        N_HCOut,NPixX,NPixY)  ;
    LDEPertOutAll{TestInd} = [LargeIter;SmallIter];
    
    for EpcInd = 1:EpocTest+1
        LDEPertRslt = LDEPertOutAll{TestInd}{EpcInd};
        SDiff = LDEPertRslt.S - LDEequv.S;
        SDiffHC = reshape(SDiff,N_HCOut*NPixY,N_HCOut*NPixX);
        SDiff = reshape(SDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);
        CDiff = LDEPertRslt.C - LDEequv.C;
        CDiffHC = reshape(CDiff,N_HCOut*NPixY,N_HCOut*NPixX);
        CDiff = reshape(CDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);
        IDiff = LDEPertRslt.I - LDEequv.I;
        IDiffHC = reshape(IDiff,N_HCOut*NPixY,N_HCOut*NPixX);
        IDiff = reshape(IDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);
        
        L2Diff_OneStepAll(TestInd,EpcInd) = ...
            sqrt(sum([SDiff;CDiff;IDiff].^2));
        DiffVecAll(TestInd,EpcInd,:) = [SDiff;CDiff;IDiff];
    end
    fprintf("sample %d is done.\n",TestInd)
end
toc
save([SaveFolder sprintf('GlobConv_IC%d_Ctrst%d.mat',...
    ICId,ContrastID)],'L2Diff_OneStepAll','DiffVecAll')
end

% rotate ICs
% LDEequv contains SCI
function [LDEequvUse] = ICRot(LDEequv,rotID,N_HCOut,NPixX,NPixY)
fields = fieldnames(LDEequv);
LDEequvUse = LDEequv;

[HCX, HCY] = meshgrid(1:N_HCOut,1:N_HCOut);
HCX = mod(HCX,2); HCY = mod(HCY,2); 
for FInd = 1:length(fields)
    CurrFieldVec = LDEequv.(fields{FInd});
    CurrFieldMap = reshape(CurrFieldVec,N_HCOut*NPixY,N_HCOut*NPixX);
    OneHC = rot90(CurrFieldMap(1:NPixY,1:NPixX),rotID); % get rotated one HC
    
    MapOut = zeros(size(CurrFieldMap));
    for xid = 1:N_HCOut
        for yid = 1:N_HCOut
            xmod = HCX(yid, xid); ymod = HCY(yid, xid);
            HCHold = OneHC;
            if xmod == 0
                HCHold = HCHold(:,end:-1:1);
            end
            
            if ymod == 0
                HCHold = HCHold(end:-1:1,:);
            end
            MapOut((yid-1)*NPixY+1:yid*NPixY, (xid-1)*NPixX+1:xid*NPixX) = HCHold;
        end
    end
    LDEequvUse.(fields{FInd})= reshape(MapOut,...
        length(CurrFieldVec),1);
end

end