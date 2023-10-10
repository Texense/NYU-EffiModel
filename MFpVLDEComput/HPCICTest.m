%% HPC for large amount of IC tests
% RadiusInd: 10^-RadiusInd when using
% SampleSize: specify an integer. Fix
% SaveID Just for saving tag
function [] = HPCICTest(RadiusInd,SampleSize,SaveID)
CurrentFolder = pwd
FigurePath = [CurrentFolder '/Figures/Demo022722/'];
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/ICTest/']; % V1D2
%SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/']; % V1
addpath(SaveFolder)
% DataFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/16Function_Scheme/']; % V1D2
% addpath(DataFolder)
addpath([CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/'])

load('ICTestWS.mat')
CurrentFolder = pwd
SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/ICTest/']; % V1D2
%SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/']; % V1
addpath(SaveFolder)
%% 4.1 Test a large sample
NSample = SampleSize;
Radius = 10^(-RadiusInd);
NoiseAll = normrnd(0,Radius,NPixY*NPixX*3,NSample);
%NoiseUse = [];
NoiNorm  = sqrt(sum(NoiseAll.^2,1));
NoiseAll = NoiseAll./repmat(NoiNorm,NPixY*NPixX*3,1)*Radius;
% OutBound = (NoiNorm>1e-7) | (NoiNorm<1e-10);
% NoiseAll(:,OutBound) = [];

LDEequv = LDEEpoOutAll{1}{end};

LDEPertOutAll = cell(NSample,1);
EpocTest = 10;
L2Diff_OneStepAll = zeros(NSample,EpocTest+1);
DiffVecAll = zeros(NSample,EpocTest+1,NPixY*NPixX*3);

tic
for TestInd = 1:NSample
    IniTest.S = symmHCs(LDEequv.S,N_HCOut,NPixX,NPixY,...
        NoiseAll(1:NPixY*NPixX,TestInd));
    IniTest.C = symmHCs(LDEequv.C,N_HCOut,NPixX,NPixY,...
        NoiseAll(NPixY*NPixX+1:2*NPixY*NPixX,TestInd));
    IniTest.I = symmHCs(LDEequv.I,N_HCOut,NPixX,NPixY,...
        NoiseAll(2*NPixY*NPixX+1:3*NPixY*NPixX,TestInd));
    
    [LDEPertOutAll{TestInd},~,~,~,~] = ...
        LDEIteration_16FuncMain(...
        PixInptCtgrUse,IniTest,p,EpocTest,...
        C_SS_mean,C_CS_mean,C_IS_mean,...
        C_SC_mean,C_CC_mean,C_IC_mean,...
        C_SI_mean,C_CI_mean,C_II_mean,...
        L4EmeshX,L4ImeshY,LDEFrfunc,...
        N_HCOut,NPixX,NPixY)  ;
    
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
save([SaveFolder sprintf('FewSteps_Rad%d_size%d_ID%d.mat',...
    RadiusInd,NSample,SaveID)],'L2Diff_OneStepAll','DiffVecAll')