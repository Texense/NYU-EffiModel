%% Precomputing on Planes: Give input, compute firing rates as functions of inputs
% Output: save([SaveFolder 'MFpV_LDE_' num2str(InputCtgr) '.mat'],...
%              'f_EnIOut','meanVs','SteadyIndicate','FailureIndicate','L4ERcrd','L4IRcrd')   
% Input:  Angle:ranging from 0-90 deg, although we only use 0 7.5 15 22.5
%         LGNctgr, L6ctgr: ranging from 1-4
%         LGNL6Mapctgr: 1 for linear, 2 for cosine
%         FlagLargeDom: 1 for small domain, 2 for large domain
% Version 1: L6 range switched to [10 60] Hz per L6 neuron
% Zhuo-Cheng Xiao 01/24/2022

function [] = LIFoLDEComput_BG(varargin)
CurrentFolder = pwd
%FigurePath = [CurrentFolder '/Figures'];
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/'];
addpath(SaveFolder)
DataFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/'];
addpath(DataFolder)

% FlagLargeDom determines the domain
if length(varargin)<1
    FlagLargeDom = 1;
else
    FlagLargeDom = varargin{1};
end
DomStr = {'Smaller','Small','Large'};
disp(sprintf('Computing %s domain.', DomStr{FlagLargeDom}))

DataPt = 'V4D2';
S = load(sprintf('AllMFPixPara_Paper2TuneFig1%s.mat',DataPt));
% For part one: Preparing...
N_Slgn = S.N_Slgn; N_Clgn = S.N_Clgn; N_Ilgn = S.N_Ilgn;
NS_L6 = S.NS_L6; NC_L6 = S.NC_L6; NI_L6 = S.NI_L6;
N_HC = S.N_HC; NPixX = S.NPixX; NPixY = S.NPixY;
FrSPixVec = S.FrSPixVec;
FrCPixVec = S.FrCPixVec;
FrIPixVec = S.FrIPixVec;
% For part two: parfor...
S_EE=S.S_EE; S_EI=S.S_EI; S_IE=S.S_IE; S_II=S.S_II; p_EEFail=S.p_EEFail; 
S_EL6=S.S_EL6; S_IL6=S.S_IL6; 
S_amb=S.S_amb; rS_amb=S.rS_amb; rC_amb=S.rC_amb; rI_amb=S.rI_amb;                                 
S_Elgn=S.S_Elgn; S_Ilgn=S.S_Ilgn; 
gL_E=S.gL_E; gL_I=S.gL_I; Ve=S.Ve; Vi=S.Vi;  tau_ref=S.tau_ref; 

tau_ampa_R=S.tau_ampa_R; tau_ampa_D=S.tau_ampa_D; 
tau_nmda_R=S.tau_nmda_R; tau_nmda_D=S.tau_nmda_D; 
tau_gaba_R=S.tau_gaba_R; tau_gaba_D=S.tau_gaba_D; 
rhoE_ampa=S.rhoE_ampa; rhoE_nmda=S.rhoE_nmda; 
rhoI_ampa=S.rhoI_ampa; rhoI_nmda=S.rhoI_nmda; 
LGNFreq = S.LGNFreq;
N_SS = S.N_SS;  N_SC = S.N_SC;  N_SI = S.N_SI; 
N_CS = S.N_CS;  N_CC = S.N_CC;  N_CI = S.N_CI;            
N_IS = S.N_IS;  N_IC = S.N_IC;  N_II = S.N_II;
clear S
%% Prepare for a canonical pixel
NL6S = NS_L6; NL6C = NC_L6; NL6I = NI_L6; 
rL6_One = 0.007;
rL6SU = NL6S*rL6_One; %rL6EU = rL6E(InputCtgr); 
rL6CU = NL6C*rL6_One;
rL6IU = NL6I*rL6_One;

rLGN_One = 0.02;
lgn_S = rLGN_One;
lgn_C = rLGN_One;
lgn_I = rLGN_One; 
% Determine L4 Input proportions from a range
fS = 2.5; fC = 8; fI = 16;

L4SE = (N_SS*fS + N_SC*fC);
L4SI =  N_SI*fI ;
L4CE = (N_CS*fS + N_CC*fC);
L4CI =  N_CI*fI ;
L4IE = (N_IS*fS + N_IC*fC);
L4II =  N_II*fI ;

L4Eall = L4SE+L4CE+L4IE; L4Iall = L4SI+L4CI+L4II;
L4SEp = mean(L4SE./L4Eall);L4CEp = mean(L4CE./L4Eall);L4IEp = mean(L4IE./L4Eall);
L4SIp = mean(L4SI./L4Iall);L4CIp = mean(L4CI./L4Iall);L4IIp = mean(L4II./L4Iall);

LineFit = [L4Iall/L4Eall,0];

% domain time scale: 160000 cost 3hrs:
switch FlagLargeDom
    case 1
        L4ERange =               0:L4Eall/100:L4Eall*3; % small and dense
        L4IDiffRange = -0.8*L4Iall:L4Iall/100:0.8*L4Iall;
    case 2
        L4ERange =               0:L4Eall/10 :L4Eall*8;
        L4IDiffRange = -3*L4Iall  :L4Iall/20 :3*L4Iall;        
    case 3
        L4ERange =               0:L4Eall/3  :L4Eall*30; % large and coarse
        L4IDiffRange = -15*L4Iall :L4Iall/3  :15*L4Iall;
end

%% Start parallel computation
cluster = gcp('nocreate');
if isempty(cluster)
cluster = parpool("local",[4,128]);
end
addAttachedFiles(cluster, {'AllMFPixPara_Paper2TuneFig1V4D2.mat'}); 

a0 = length(L4ERange)*length(L4IDiffRange);
f_EnIOut = cell(a0,1);
% meanVs = cell(a0,1);
% SteadyIndicate = zeros(a0,1);
% FailureIndicate = zeros(a0,1);
L4ERcrd = zeros(a0,1);
L4IRcrd = zeros(a0,1);
HyperPara = {100*1e3};
parfor  LDEInd = 1:a0
        L4EInd = ceil(LDEInd/length(L4IDiffRange));
        L4IDiffInd = mod(LDEInd,length(L4IDiffRange));
        if L4IDiffInd == 0
            L4IDiffInd = length(L4IDiffRange);
        end
        L4EU = L4ERange(L4EInd);
        L4IU = LineFit(1)*L4EU+LineFit(2) + L4IDiffRange(L4IDiffInd);
        L4ERcrd(LDEInd) = L4EU;
        L4IRcrd(LDEInd) = L4IU;
        
        if min(L4EU,L4IU)<0
            continue
        end
        %Distribute L4Input
        L4SEU = L4SEp*L4EU; L4CEU = L4CEp*L4EU; L4IEU = L4IEp*L4EU;
        L4SIU = L4SIp*L4IU; L4CIU = L4CIp*L4IU; L4IIU = L4IIp*L4IU;
         
        tic
       [f_EnIOut{LDEInd}]...
           = LIFo_SinglePixel(...
...% MF Parameters                     
                     LGNFreq, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,... %3 
                     S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                     S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                     lgn_S, lgn_C,lgn_I,N_Slgn,N_Clgn,N_Ilgn, S_Elgn,S_Ilgn,... %7
                     gL_E,gL_I,Ve,Vi, tau_ref,... %5
...% Below are LIF details
                     tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,... %7
                     rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                     HyperPara) % if varagin non empty, we only export the last state
        toc   
end

% Save MFpV Data
%CurrentFolder = pwd;
%SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/'];

save([SaveFolder ...
      sprintf('Paper2_LIFoLDE_Fig4Bg%s_%s.mat',DataPt, DomStr{FlagLargeDom})],...
    'f_EnIOut','L4ERcrd','L4IRcrd')    
end