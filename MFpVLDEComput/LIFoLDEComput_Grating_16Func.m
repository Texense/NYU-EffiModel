%% Precomputing on Planes: Give input, compute firing rates as functions of inputs
% Output: save([SaveFolder 'MFpV_LDE_' num2str(InputCtgr) '.mat'],...
%              'f_EnIOut','meanVs','SteadyIndicate','FailureIndicate','L4ERcrd','L4IRcrd')   
% Input:  Angle:ranging from 0-90 deg, although we only use 0 7.5 15 22.5
%         LGNctgr, L6ctgr: ranging from 1-4
%% NOTE: Now we are computing for binocular, and has to mix BG and drive inputs
%         LGNctgr, L6ctgr: ranging from 1-5!! 5 means background!
%         LGNL6Mapctgr: 1 for linear, 2 for cosine
%         FlagLargeDom: 1 for small domain, 2 for large domain
% Version 1: L6 range switched to [10 60] Hz per L6 neuron
% Zhuo-Cheng Xiao 01/24/2022

function [] = LIFoLDEComput_Grating_16Func(Angle, LGNctgr, L6ctgr, varargin)
CurrentFolder = pwd
%FigurePath = [CurrentFolder '/Figures'];
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/25Function_Binocular/'];
addpath(SaveFolder)
DataFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/'];
addpath(DataFolder)
DataFolder1 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/'];
addpath(DataFolder1)
% LGNL6Mapctgr determines the mapping of input to LGN/L6
if length(varargin)<1
    LGNL6Mapctgr = 1; 
else
    LGNL6Mapctgr = varargin{1};
end
% FlagLargeDom determines the domain
if length(varargin)<2
    FlagLargeDom = 1;
else
    FlagLargeDom = varargin{2};
end
DomStr = {'','Large','LARGER'};
% if ~ismember(InputCtgr,[1,2,3])
%     error('No such input category. Exit.')
% end

DataPt = 'V4D2';
S = load(sprintf('AllMFPixPara_Paper2TuneFig1%s.mat',DataPt));
% For part one: Preparing...
C_SS_Pixel_Us = S.C_SS_Pixel_Us;
C_CS_Pixel_Us = S.C_CS_Pixel_Us;
C_IS_Pixel_Us = S.C_IS_Pixel_Us;
C_SC_Pixel_Us = S.C_SC_Pixel_Us;
C_CC_Pixel_Us = S.C_CC_Pixel_Us;
C_IC_Pixel_Us = S.C_IC_Pixel_Us;
C_SI_Pixel_Us = S.C_SI_Pixel_Us;
C_CI_Pixel_Us = S.C_CI_Pixel_Us;
C_II_Pixel_Us = S.C_II_Pixel_Us;
N_Slgn = S.N_Slgn; N_Clgn = S.N_Clgn; N_Ilgn = S.N_Ilgn;
NS_L6 = S.NS_L6; NC_L6 = S.NC_L6; NI_L6 = S.NI_L6;
L6Ord_F  = S.L6Ord_F;
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
clear S
%% Prepare for a canonical pixel
Angles_4Input = Angle + [0 45 90 135];
lgn_SOnOff = [45+abs(mod(Angles_4Input,180)-90)/2;
              45-abs(mod(Angles_4Input,180)-90)/2]/1e3;              
% redefine low-up bounds for L6: [10 54]
NL6S = NS_L6; NL6C = NC_L6; NL6I = NI_L6;
L6up = 60; L6low = 6; %[10 54; 12 50]
if LGNL6Mapctgr == 1
   lgn_SOnOff = [45+abs(mod(Angles_4Input,180)-90)/2;
                 45-abs(mod(Angles_4Input,180)-90)/2]/1e3;  
   FL6_Angle1 = ((abs(mod(Angles_4Input,180)-90)/90)      *(L6up-L6low)+L6low) /1e3; % get L6 frs
   ExpTex = 'Linear';
elseif LGNL6Mapctgr == 2
   lgn_SOnOff = [45+abs(mod(Angles_4Input,180)-90)/2;
                 45-abs(mod(Angles_4Input,180)-90)/2]/1e3;  
   FL6_Angle1 = ((cosd(abs(mod(Angles_4Input,180))*2)+1)/2*(L6up-L6low)+L6low) /1e3; 
   ExpTex = 'Cosine_L6';
elseif LGNL6Mapctgr == 3
   lgn_SOnOff = [45+(cosd(abs(mod(Angles_4Input,180))*2)+1)/2*45;
                 45-(cosd(abs(mod(Angles_4Input,180))*2)+1)/2*45]/1e3;  
   FL6_Angle1 = ((cosd(abs(mod(Angles_4Input,180))*2)+1)/2*(L6up-L6low)+L6low) /1e3; 
   
   ExpTex = 'Cosine_All';   
else
   disp("illigal LGN->L6 mapping.")
%rL6E = [1.0; 1.75; 2.5]; rL6I = 3*rL6E;
end    
%% NOTE: Now we are computing for binocular, and has to mix BG and drive 
FL6_Angle1 = [FL6_Angle1, 0.007]; % background L6 frs
lgn_SOnOff = [lgn_SOnOff, [0.02;0.02]];
lgn_COnOff = [0.045*ones(2,4),[0.02;0.02]];
lgn_IOnOff = [0.045*ones(1,4),0.02]; 

rL6SU = NL6S*FL6_Angle1(L6ctgr); %rL6EU = rL6E(InputCtgr); 
rL6CU = NL6C*FL6_Angle1(L6ctgr);
rL6IU = NL6I*FL6_Angle1(L6ctgr);
lgn_SU = lgn_SOnOff(:,LGNctgr);
lgn_C = lgn_COnOff(:,LGNctgr);
lgn_I = lgn_IOnOff(LGNctgr);
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
LineFit = polyfit(L4Eall,L4Iall,1);
YRange = L4Iall - (LineFit(1)*L4Eall+LineFit(2)); % We plan to use 2 times of the range
Bdry = floor(max(abs(YRange))/100)*100;
% domain time scale: 160000 cost 3hrs:
switch FlagLargeDom
    case 1
        L4ERange = 0:100:ceil(max(L4Eall)*1.5/100)*100;
        L4IDiffRange = -20*Bdry:100:20*Bdry;
    case 2
        L4ERange = 0:600:ceil(max(L4Eall)*12/100)*100;
        L4IDiffRange = -120*Bdry:600:60*Bdry;
    case 3
        L4ERange = 0:2000:ceil(max(L4Eall)*50/100)*100;
        L4IDiffRange = -500*Bdry:2000:500*Bdry;
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
HyperPara = {50*1e3};
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
                     lgn_SU, lgn_C,lgn_I,N_Slgn,N_Clgn,N_Ilgn, S_Elgn,S_Ilgn,... %7
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
      sprintf('Paper2_LIFoLDE_Fig1%s_Ang%.1f_LGNc%d_L6c%d_L6_%d_%d_%s%s.mat',...
      DataPt, Angle, LGNctgr, L6ctgr, L6up,L6low, ExpTex,DomStr{FlagLargeDom})],...
    'f_EnIOut','L4ERcrd','L4IRcrd')    
end