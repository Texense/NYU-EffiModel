%% Paper3: 1st test
% Testing different combinations of x and y for E and I cells.
% Output: save the turning curve data
% Input xE, xI, yE, yI;
% Version 1: xE, xI only. yE and yI are kept as 1 from the beginning to end
% Zhuo-Cheng Xiao 06/18/2023

function [] = HPC_TestEIxy(xE,xI)
CurrentFolder = pwd
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
DataFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper3PlotingData/']; % V1D2
DataFolder1 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/GlobConv/']; % V1D2
DataFolder2 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper2PlotingData/Typical_trajs/']; % V1D2
DataFolder3 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper2PlotingData/']; % V1D2
SaveToFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper3TestingKernelData/']; % V1D2

addpath(DataFolder3)
addpath(DataFolder2)
addpath(DataFolder1)
addpath(DataFolder)
addpath([CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/'])

% FigurePaperPath = [CurrentFolder '/Figures/Demo061123/']; % V1D2
%% B1. Load LDE computations and do Fr distributions
% Should rerun Sect 1 after in case loaded file contain any paths.
DataPoint = 'V4D2'; %'V4D1'; 'V5D1'
load(sprintf("AllMFPixPara_Paper2TuneFig1%s.mat",DataPoint))

% CurrentFolder = pwd;
% FigurePaperPath = [CurrentFolder '/Figures/Demo061123/']; % V1D2
close all

%% 2. Test for different x for E and I
%first get basic function information
% Too time comsuming. --->> using HPC
fprintf('xE=%.1f, xI=%.1f\n',xE,xI)
tic
%AngleList = {'0.0'}; % only vertical for now
AngleList = {'0.0','7.5','15.0','22.5'};
DomList =  {'Small','Large','LARGER'};
%% 3. Set up initial conditions
N_HCOutX = 4; N_HCOutY = 8; PixNumOut = N_HCOutX*N_HCOutY*NPixX*NPixY;
N_HCinX = 4; N_HCinY = 4;
NPixX = 10; NPixY = 10;
%Parameters of L4 kernels
MonocuFlag = true;
AngleTestAll = 0:7.5:180; % 0:7.5:22.5; % has to be multiples of 7.5
TestNum = length(AngleTestAll);
TestExSeq = xE*ones(size(0:7.5:180)); TestEySeq = 1*ones(size(TestExSeq));
TestIxSeq = xI*ones(size(0:7.5:180)); TestIySeq = 1*ones(size(TestExSeq));
%Background
% LDEbgIC.S = 2.5*ones(N_HCOutX*NPixX*N_HCOutY*NPixY,1);
% LDEbgIC.C = 8*ones(N_HCOutX*NPixX*N_HCOutY*NPixY,1);
% LDEbgIC.I = 18*ones(N_HCOutX*NPixX*N_HCOutY*NPixY,1); %16
% LDEICOut = LDEbgIC;
%%%60 deg: by rotating 15 deg pattern
%Start from 90, then go to 22.5
ICData = load([DataFolder2 sprintf('LDETracesV4D2_NewSmear_ang0.0.mat')]);
load('LDETracesV4D2_NewSmear_ang0.0.mat')
LDEICPart = ICData.LDEEquv;

LDEICOut = HCRotXY(LDEICPart,0,N_HCinX,N_HCinY,N_HCOutX,N_HCOutY,NPixX,NPixY,0);
for TestId = 1:TestNum
    ICTestAll{TestId} = LDEICOut;
end
LDEOutAll = cell(size(ICTestAll));
NANFlag = false(size(ICTestAll));
%BGFlagAllTest = cell(size(ICTestAll));
%load equiv
%PlotUse = load([DataFolder2 sprintf('LDETraces_bg-ang%s.mat',AngPrint)],'LDEEquv');
%LDEEquv = PlotUse.LDEEquv;
%LDEEquv = HCRot(LDEEquv,RotInd,N_HCOut,NPixX,NPixY,MirInd);
%Hyperparameters of iterations
p = 0.33;% 1-p for the original input
EpocTest = 50;
ExportFlag = 'xn';

L2Diff_OneStepAll = zeros(length(ICTestAll),EpocTest+1);
DiffVecAll = zeros(length(ICTestAll),EpocTest+1,3*NPixY*NPixX);
%Since we are doing different network architectures now - need to separate inputs to the S C I functions!!!
% Determine L4 Input from a range
L4SE = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4SI = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4CE = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4CI = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4IE = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4II = zeros(N_HC*NPixX*N_HC*NPixY,1);
for PInd = 1:N_HC*NPixX*N_HC*NPixY
    L4SE(PInd) = (C_SS_Pixel_Us(PInd,:)*FrSPixVec          + C_SC_Pixel_Us(PInd,:)*FrCPixVec);
    L4SI(PInd) =  C_SI_Pixel_Us(PInd,:)*FrIPixVec ;
    L4CE(PInd) = (C_CS_Pixel_Us(PInd,:)*FrSPixVec          + C_CC_Pixel_Us(PInd,:)*FrCPixVec);
    L4CI(PInd) =  C_CI_Pixel_Us(PInd,:)*FrIPixVec ;
    L4IE(PInd) = (C_IS_Pixel_Us(PInd,:)*FrSPixVec          + C_IC_Pixel_Us(PInd,:)*FrCPixVec);
    L4II(PInd) =  C_II_Pixel_Us(PInd,:)*FrIPixVec ;
end
L4Eall = L4SE+L4CE+L4IE; L4Iall = L4SI+L4CI+L4II;
L4SEp = mean(L4SE./L4Eall);L4CEp = mean(L4CE./L4Eall);L4IEp = mean(L4IE./L4Eall);
L4SIp = mean(L4SI./L4Iall);L4CIp = mean(L4CI./L4Iall);L4IIp = mean(L4II./L4Iall);
%Compute NESSs for each angle
for TestId = 1:TestNum
    AngleInpt = AngleTestAll(TestId);

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

    AngPrint = AngleList{AngFuncCtgr};
    Angle = str2num(AngPrint); % can select from 0 7.5 15 22.5
    %load all precomputed domains
    L4EmeshXAll = cell(size(DomList));
    L4ImeshYAll = cell(size(DomList));
    LDEFrfuncAll = cell(size(DomList));
    for DomInd = 1:length(DomList)
        FuncTemp = load(sprintf('Func25%sAng%s%s.mat',DataPoint,AngPrint,DomList{DomInd}),...
            'LDEFrfunc','L4EmeshX','L4ImeshY');
        L4EmeshXAll{DomInd} = FuncTemp.L4EmeshX;
        L4ImeshYAll{DomInd} = FuncTemp.L4ImeshY;
        LDEFrfuncAll{DomInd} = FuncTemp.LDEFrfunc;
    end
    % the finest BG function as well
    FuncTemp_BG = load('FuncBG_Smaller_SM15_5.mat','LDEFrfunc','L4EmeshX','L4ImeshY');
    L4EmeshXBG = FuncTemp_BG.L4EmeshX;
    L4ImeshYBG = FuncTemp_BG.L4ImeshY;
    LDEFrfuncBG = FuncTemp_BG.LDEFrfunc;
    %% 4. Iteration
    %New: Need to extend all these to arbitrary HC numbers!
    %    N_HCOut = 4; PixNumOut = N_HCOut^2*NPixX*NPixY;
    LGNlist = 5; L6list = 5;
    OD_SMapModi = OD_SMap;
    %Needs to replace the midel row as "5" -- Background
    %OD_SMapModi = OD_SMap(1:2*n_S_HC,1:2*n_S_HC);
    if MonocuFlag
        OD_SMapModi(n_S_HC+1:2*n_S_HC,1:end) = 5;
        SaveStr = 'Monocu';
    else
        SaveStr = 'Binocu';
    end
    %LGN and L6 are blurred
    N_HCs = 3; % shrink to 2*2
    FigOn = false;
    L6smear = n_S_HC*0.25; LGNsmear = n_S_HC*0.2;
    TruncL6 = 0.33/(L6smear/n_S_HC); TruncLGN = 1;
    L6SFilt_Grating = SpatialGaussianFilt_my(OD_SMapModi,...
        N_HCs,n_S_HC,L6smear,TruncL6,FigOn);
    LGNFilt_Grating = SpatialGaussianFilt_my(OD_SMapModi,...
        N_HCs,n_S_HC,LGNsmear,TruncLGN,FigOn);

    if ~MonocuFlag
        L6SFilt_Grating = [L6SFilt_Grating,zeros(size(L6SFilt_Grating,1),1)];
        LGNFilt_Grating = [LGNFilt_Grating,zeros(size(LGNFilt_Grating,1),1)];
    end

    PixL6Ctgr = LGNIndSpat_Rec(...
        L6SFilt_Grating,1:L6list,NnSPixel,N_HC,N_HCOutX,N_HCOutY,NPixX,NPixY,false);
    PixLGNCtgr = LGNIndSpat_Rec(...
        LGNFilt_Grating,1:LGNlist,NnSPixel,N_HC,N_HCOutX,N_HCOutY,NPixX,NPixY,true);

    % symmetrize !! This is actually tricky, so I do it ad hoc
    PixL6Ctgr = Ocu_LGNL6symm(PixL6Ctgr,N_HCOutX,N_HCOutY,NPixX,NPixY);
    PixLGNCtgr = Ocu_LGNL6symm(PixLGNCtgr,N_HCOutX,N_HCOutY,NPixX,NPixY);
    %completely BG or not?
    BGFlag = false;
    if BGFlag
        PixL6Ctgr(:,end) = 1; PixL6Ctgr(:,1:end-1) = 0;
        PixLGNCtgr(:,end) = 1; PixLGNCtgr(:,1:end-1) = 0;
    end
    %However, the Ctgr for pixels needs to be rotated/flipped to represent the actual gratings
    %Instead of rotating existing ctgr fields, I should simply permute 1234
    MirInd = logical(MirInd);
    CtgrOrder = [2,1;3,4]; % for lower left HC
    CtgrOrderUse = rot90(CtgrOrder,RotInd);
    if MirInd
        %flipDim = mod(RotInd,2)+1; %1 for col 2 ofr row
        CtgrOrderUse = flip(CtgrOrderUse,1);
    end
    CtgrOrderReadout = ...
        [CtgrOrderUse(1,2);
        CtgrOrderUse(1,1);
        CtgrOrderUse(2,1);
        CtgrOrderUse(2,2)];
    %    for Ctgr = 1:4
    %         PixL6Ctgr(:,Ctgr) = ...
    %             HCRot_Rec(PixL6Ctgr(:,Ctgr),RotInd,N_HCOutX,N_HCOutY,NPixX,NPixY,MirInd);
    %         PixLGNCtgr(:,Ctgr) = ...
    %             HCRot_Rec(PixLGNCtgr(:,Ctgr),RotInd,N_HCOutX,N_HCOutY,NPixX,NPixY,MirInd);
    %    end

    % let's implement something cheap now for 45 deg
    PixL6Ctgr = PixL6Ctgr(:,[CtgrOrderReadout',5]);
    PixLGNCtgr = PixLGNCtgr(:,[CtgrOrderReadout',5]);

    PixInptCtgrUse  = zeros(PixNumOut,LGNlist,L6list);
    for PixInd = 1:PixNumOut
        PixInptCtgrUse(PixInd,:,:) = PixLGNCtgr(PixInd,:)' *PixL6Ctgr(PixInd,:);% by multiplying both indexes
    end
    %Get pixelwise connectivities
    %NOTE: NEED to modify AveSpatKer_Rec to implement the ocular modulation of connectivity kernels
    ParaEx = TestExSeq(TestId); ParaEy = TestEySeq(TestId);
    ParaIx = TestIxSeq(TestId); ParaIy = TestIySeq(TestId);

    C_SS_mean = sparse(AveSpatKer_Rec(C_SS_Pixel_Us,N_HC,...
        N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
    C_CS_mean = sparse(AveSpatKer_Rec(C_CS_Pixel_Us,N_HC,...
        N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
    C_IS_mean = sparse(AveSpatKer_Rec(C_IS_Pixel_Us,N_HC,...
        N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
    C_SC_mean = sparse(AveSpatKer_Rec(C_SC_Pixel_Us,N_HC,...
        N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
    C_CC_mean = sparse(AveSpatKer_Rec(C_CC_Pixel_Us,N_HC,...
        N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
    C_IC_mean = sparse(AveSpatKer_Rec(C_IC_Pixel_Us,N_HC,...
        N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
    C_SI_mean = sparse(AveSpatKer_Rec(C_SI_Pixel_Us,N_HC,...
        N_HCOutX,N_HCOutY,NPixX,NPixY,ParaIx,ParaIy));
    C_CI_mean = sparse(AveSpatKer_Rec(C_CI_Pixel_Us,N_HC,...
        N_HCOutX,N_HCOutY,NPixX,NPixY,ParaIx,ParaIy));
    C_II_mean = sparse(AveSpatKer_Rec(C_II_Pixel_Us,N_HC,...
        N_HCOutX,N_HCOutY,NPixX,NPixY,ParaIx,ParaIy));

    Conn_Reduc_Test = 5;
    Reduc_S = Conn_Reduc_Test*(1-CplxR); Reduc_C = Conn_Reduc_Test*CplxR;
    P_SS_rate = (N_SS-Reduc_S)/N_SS; C_SS_meanU = C_SS_mean * P_SS_rate;
    P_SC_rate = (N_SC-Reduc_C)/N_SC; C_SC_meanU = C_SC_mean * P_SC_rate;
    P_CS_rate = (N_CS-Reduc_S)/N_CS; C_CS_meanU = C_CS_mean * P_CS_rate;
    P_CC_rate = (N_CC-Reduc_C)/N_CC; C_CC_meanU = C_CC_mean * P_CC_rate;
    %Firing rate suppresion
    %     EKp.Thrsld = 50; EKp.Highist = [200,150]; EKp.HardBound = 200; EKp.Slope = 0;
    %
    %     IKp.Thrsld = 70; IKp.Highist = [100,95];  IKp.HardBound = 120; IKp.Slope = 1; % 1; 0.9(95); 0.8(9)
    %set up IC
    tic
    ICUse = ICTestAll{TestId}; %
    IniTest = ICUse;
    %IniTest.S = symmHCs(ICUse.S,N_HCOut,NPixX,NPixY);
    %IniTest.C = symmHCs(ICUse.C,N_HCOut,NPixX,NPixY);
    %IniTest.I = symmHCs(ICUse.I,N_HCOut,NPixX,NPixY);

    %start iteration
    [LDEOutAll{TestId},~,~,~,~,FuncUse,...
        NANFlag(TestId)] = ...
        LDEIteration_25FuncMain_CombDom_BG(...
        PixInptCtgrUse,IniTest,p,EpocTest,...
        C_SS_meanU,C_CS_meanU,C_IS_mean,...
        C_SC_meanU,C_CC_meanU,C_IC_mean,...
        C_SI_mean, C_CI_mean, C_II_mean,...
        L4SEp, L4SIp, ...
        L4CEp, L4CIp, ...
        L4IEp, L4IIp, ...
        L4EmeshXAll,L4ImeshYAll,LDEFrfuncAll,...
        L4EmeshXBG, L4ImeshYBG, LDEFrfuncBG, ...
        N_HCOutY,NPixX,NPixY,false,'xn')  ;
    LDEEquv = LDEOutAll{TestId}{end};
    toc
    %             RunTime = toc;
    %             sprintf('RunTime = %.2f',RunTime)

    if ~NANFlag(TestId)
        for EpcInd = 1:EpocTest+1
            LDEPertRslt = LDEOutAll{TestId}{EpcInd};
            SDiff = LDEPertRslt.S - LDEEquv.S;
            SDiffHC = reshape(SDiff,N_HCOutY*NPixY,N_HCOutX*NPixX);
            SDiff = reshape(SDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);
            CDiff = LDEPertRslt.C - LDEEquv.C;
            CDiffHC = reshape(CDiff,N_HCOutY*NPixY,N_HCOutX*NPixX);
            CDiff = reshape(CDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);
            IDiff = LDEPertRslt.I - LDEEquv.I;
            IDiffHC = reshape(IDiff,N_HCOutY*NPixY,N_HCOutX*NPixX);
            IDiff = reshape(IDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);

            L2Diff_OneStepAll(TestId,EpcInd) = ...
                sqrt(sum((SDiff*(1-CplxR)+CDiff*CplxR).^2,'all'));%sqrt(sum([SDiff;CDiff;IDiff].^2));
            DiffVecAll(TestId,EpcInd,:) = [SDiff;CDiff;IDiff];
        end
    end
end
%Comment below if don't need to save data (usually for setting up new NESSs)
SaveStrUse = sprintf('%s_L6%.2f_%.2f_Ereduc%d_xE%.1fxI%.1fyE%.1fyI%.1f',...
    SaveStr,L6smear/n_S_HC,TruncL6,Conn_Reduc_Test,ParaEx,ParaIx,ParaEy,ParaIy);

save([SaveToFolder sprintf('NESSs_%s.mat',SaveStrUse)],...
    'LDEOutAll','L2Diff_OneStepAll','DiffVecAll','AngleTestAll')
RunTime = toc;
fprintf('RunTime = %.2f',RunTime)

end
