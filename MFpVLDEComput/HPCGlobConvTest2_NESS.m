%% HPC for large amount of IC tests. Start from NESS
% AngId: 1-13 for 0:7.5:90
% ContrastID: B. "1/2 contrast" = bg + 1/2 * (X_E-bg, X_I-bg)
%             C. (normal) full contrast = (X_E, X_I)
%             D. super-charged = bg + 3/2 * (X_E-bg, X_I-bg)
% SaveID Just for saving tag

function [] = HPCGlobConvTest2_NESS(AngId, ContrastID, NSample)
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
addpath([CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper2PlotingData/'])

load('GlobalConvTestWS1.mat')
clear LDEFrfuncLarge LDEFrfunc

%% NEW unified functions
%load('ICTesfuncLarge.mat')
%load('ICTesfuncSmall.mat')
DomList =  {'Small','Larger'}; % we have small, large, larger domains
DataPoint = 'V4D2'; %'V4D1'; 'V5D1'

Large = load(sprintf('Func16%sAng%s%s.mat',DataPoint,'0.0',DomList{2}),...
    'LDEFrfunc','L4EmeshX','L4ImeshY');

Small = load(sprintf('Func16%sAng%s%s.mat',DataPoint,'0.0',DomList{1}),...
    'LDEFrfunc','L4EmeshX','L4ImeshY');

L4EmeshXAll = {Small.L4EmeshX, Large.L4EmeshX};
L4ImeshYAll = {Small.L4ImeshY, Large.L4ImeshY};
LDEFrfuncAll = {Small.LDEFrfunc, Large.LDEFrfunc};
%% rest of folders
CurrentFolder = pwd; % Since loaded data will overwrite pwd...
SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/GlobConv/']; % V1D2
%SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/']; % V1
addpath(SaveFolder)
%% 4.1 Test a large sample
%% first prepare a driven IC
AngleAll = 0:7.5:90;
AngleInpt = AngleAll(AngId);% has to be multiples of 7.5
RotInd = floor(AngleInpt/45); % 0-3
MirInd = mod(floor(AngleInpt/22.5),2); % 0: no mirror; 1: mirror
switch MirInd % find the corresponding source
    case 0
        AngSource = mod(AngleInpt,45);
    case 1
        AngSource = mod(-AngleInpt,45);
end

% Use the source, decide which response function for angle should I use
AngFuncCtgrCdid = 0:7.5:22.5;
[~,AngFuncCtgr] = min(abs(AngFuncCtgrCdid - AngSource));

AngleList = {'0.0','7.5','15.0','22.5'};
AngPrint = AngleList{AngFuncCtgr};
load(sprintf('LDETraces_ang%s.mat',AngPrint))
LDEICPart = LDEEpoOutAll{1}{end};
LDEICPart = HCRot(LDEICPart,RotInd,N_HCOut,NPixX,NPixY,MirInd); % rotate to the write direction

% Compose contrasts
ICWeights = 1/2:1/2:3/2;
fields = fieldnames(LDEICPart);
for FInd = 1:length(fields)
    LDEIC_temp.(fields{FInd}) = ...
        ICWeights(ContrastID) * LDEICPart.(fields{FInd});
end

% Now perturbe and collect a group
PertBckSz = 2; % size of the block
BlockN = floor(NPixX/PertBckSz);
% Use Kronecker Tensor Product
ICTestAll = cell(NSample,1);
for SampInd = 1:NSample
    LDEIC_pert = LDEIC_temp;
    Block = ones(PertBckSz);
    PertMtx = ones(BlockN);
    
    PertPlus  = rand(BlockN)    +1; % random in [1,2]
    Pertminus = rand(BlockN)/2+1/2; % random in [1/2,1]
    
    [PlusInd, MinusInd] = RandPermIndMat(BlockN,BlockN,randi([5 20])); % randomly use the random numbers above
    PertMtx(PlusInd)  = PertPlus(PlusInd);
    PertMtx(MinusInd) = Pertminus(MinusInd);
    PertMtxUse = kron(PertMtx,Block);
    % take a blockwise peerturbation
    
    LDEIC_pert.S = symmHCs(LDEIC_pert.S,N_HCOut,NPixX,NPixY,PertMtxUse,'prod');
    LDEIC_pert.C = symmHCs(LDEIC_pert.C,N_HCOut,NPixX,NPixY,PertMtxUse,'prod');
    ICTestAll{SampInd} = LDEIC_pert;
    
end

%% Redefine L6 and LGN inpt category
FigOn = false;
L6smear = n_S_HC*0.34; LGNsmear = n_S_HC*0.2;
TruncL6 = 1.25; TruncLGN = 1;
L6SFilt_Grating = SpatialGaussianFilt_my(OD_SMap,...
    N_HC,n_S_HC,L6smear,TruncL6,FigOn);
LGNFilt_Grating = SpatialGaussianFilt_my(OD_SMap,...
    N_HC,n_S_HC,LGNsmear,TruncLGN,FigOn);

PixL6Ctgr = LGNIndSpat(L6SFilt_Grating,1:4,NnSPixel,N_HC,N_HCOut,NPixX,NPixY);
PixLGNCtgr = LGNIndSpat(LGNFilt_Grating,1:4,NnSPixel,N_HC,N_HCOut,NPixX,NPixY);
% However, the Ctgr for pixels needs to be rotated/flipped to represent the actual gratings
MirInd = logical(MirInd);
for Ctgr = 1:4
    PixL6Ctgr(:,Ctgr) = HCRot(PixL6Ctgr(:,Ctgr),RotInd,N_HCOut,NPixX,NPixY,MirInd);
    PixLGNCtgr(:,Ctgr) = HCRot(PixLGNCtgr(:,Ctgr),RotInd,N_HCOut,NPixX,NPixY,MirInd);
end

PixInptCtgrUse  = zeros(PixNumOut,LGNlist,L6list);
for PixInd = 1:PixNumOut
    PixInptCtgrUse(PixInd,:,:) = PixLGNCtgr(PixInd,:)' *PixL6Ctgr(PixInd,:);% by multiplying both indexes
end


%% Simulate for each IC
LDEPertOutAll = cell(NSample,1);
EpocTest = 200;
L2Diff_EIWgtsAll = zeros(NSample,EpocTest+1);
DiffVecAll = zeros(NSample,EpocTest+1,NPixY*NPixX*3);
p = 0.33;

tic
for TestInd = 1:NSample
    ICUse = ICTestAll{TestInd};
    IniTest.S = symmHCs(ICUse.S,N_HCOut,NPixX,NPixY);
    IniTest.C = symmHCs(ICUse.C,N_HCOut,NPixX,NPixY);
    IniTest.I = symmHCs(ICUse.I,N_HCOut,NPixX,NPixY);
    
    % Function: small large combined
    [LDEPertOutAll{TestInd},~,~,~,~] = ...
        LDEIteration_16FuncMain_CombDom(...
        PixInptCtgrUse,IniTest,p,EpocTest,...
        C_SS_mean,C_CS_mean,C_IS_mean,...
        C_SC_mean,C_CC_mean,C_IC_mean,...
        C_SI_mean,C_CI_mean,C_II_mean,...
        L4EmeshXAll,L4ImeshYAll,LDEFrfuncAll,...
        N_HCOut,NPixX,NPixY,false,'xn')  ;
    LDEequv = LDEPertOutAll{TestInd}{end};
    
    if ~isempty(LDEequv)
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
            
            EDiff = SDiff*(1-CplxR) + CDiff*CplxR;
            L2Diff_EIWgtsAll(TestInd,EpcInd) = ...
                sqrt(sum(EDiff.^2 * 0.8^2 + IDiff.^2 * 0.2^2)/(NPixY*NPixX));
            DiffVecAll(TestInd,EpcInd,:) = [SDiff;CDiff;IDiff];
        end
    else
        sprintf('***Test %d fails.***', TestInd)
    end
    fprintf("sample %d is done.\n",TestInd)
end
toc
save([SaveFolder sprintf('Paper2GlobConv_IC%d_Ctrst%d_p%.2f_NewSmear.mat',...
    AngId,ContrastID,p)], 'LDEPertOutAll','L2Diff_EIWgtsAll','DiffVecAll')
end

% generate random binary matrix
% p row q col
% plus and minus randomly selected from n entries
function [PlusInd, MinusInd] = RandPermIndMat(p,q,n)
a = zeros(p,q);
perturbIndBoth = randperm(numel(a), n);
PMInd = logical(randi(2,size(perturbIndBoth,1),size(perturbIndBoth,2))-1);

PlusInd = a; MinusInd = a;
PlusInd(perturbIndBoth(PMInd)) = 1;
MinusInd(perturbIndBoth(~PMInd)) = 1;

PlusInd = logical(PlusInd);
MinusInd = logical(MinusInd);
end