%% Precomputing on Planes: Give input, compute firing rates as functions of inputs
% Output: save([SaveFolder 'MFpV_LDE_' num2str(InputCtgr) '.mat'],...
%              'f_EnIOut','meanVs','SteadyIndicate','FailureIndicate','L4ERcrd','L4IRcrd')   
% Input:  InputCtgr: Choices from 1-3. Orthogonal, suboptimal, optimal

% Version 0: Get all parameters from a para file, and I only specify two
% input.
% Zhuo-Cheng Xiao 08/20/2021

InputCtgr = 3
CurrentFolder = pwd
%FigurePath = [CurrentFolder '/Figures'];
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
SaveFolder = [CurrentFolder '/Figures/Demo082721/'];
addpath(SaveFolder)
DataFolder = [CurrentFolder '/Figures/Demo082121/'];
addpath(DataFolder)

if InputCtgr>3 || InputCtgr<1
    error('No such input category. Exit.')
end

load('AllMFPixPara.mat')
%% Prepare for a canonical pixel
% Connection within a pixel
N_PreSynPix = [mean(diag(C_SS_Pixel_Us)),mean(diag(C_CS_Pixel_Us)),mean(diag(C_IS_Pixel_Us));
               mean(diag(C_SC_Pixel_Us)),mean(diag(C_CC_Pixel_Us)),mean(diag(C_IC_Pixel_Us));
               mean(diag(C_SI_Pixel_Us)),mean(diag(C_CI_Pixel_Us)),mean(diag(C_II_Pixel_Us))];
% external drives: lgn and L6
NlgnS = N_Slgn; NlgnC = N_Clgn; NlgnI = N_Ilgn;
rL6E = [1.0; 1.75; 2.5]; rL6I = 3*rL6E;
lgn_SOnOff = [0.045, 0.0675, 0.09;
              0.045, 0.0225, 0.00];
lgn_COnOff = [0.045;
              0.045];
lgn_I =       0.045;
rL6EU = rL6E(InputCtgr); rL6IU = rL6I(InputCtgr);
lgn_SU = lgn_SOnOff(:,InputCtgr);
% Determine L4 Input from a range
L4SE = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4SI = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4CE = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4CI = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4IE = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4II = zeros(N_HC*NPixX*N_HC*NPixY,1);
for PInd = 1:N_HC*NPixX*N_HC*NPixY
L4SE(PInd) = (C_SS_Pixel_Us(PInd,:)*FrSPixVec          + C_SC_Pixel_Us(PInd,:)*FrCPixVec)...
            - C_SS_Pixel_Us(PInd,PInd)*FrSPixVec(PInd) - C_SC_Pixel_Us(PInd,PInd)*FrCPixVec(PInd);
L4SI(PInd) =  C_SI_Pixel_Us(PInd,:)*FrIPixVec - C_SI_Pixel_Us(PInd,PInd)*FrIPixVec(PInd);
L4CE(PInd) = (C_CS_Pixel_Us(PInd,:)*FrSPixVec          + C_CC_Pixel_Us(PInd,:)*FrCPixVec)...
            - C_CS_Pixel_Us(PInd,PInd)*FrSPixVec(PInd) - C_CC_Pixel_Us(PInd,PInd)*FrCPixVec(PInd);
L4CI(PInd) =  C_CI_Pixel_Us(PInd,:)*FrIPixVec - C_CI_Pixel_Us(PInd,PInd)*FrIPixVec(PInd);
L4IE(PInd) = (C_IS_Pixel_Us(PInd,:)*FrSPixVec          + C_IC_Pixel_Us(PInd,:)*FrCPixVec)...
            - C_IS_Pixel_Us(PInd,PInd)*FrSPixVec(PInd) - C_IC_Pixel_Us(PInd,PInd)*FrCPixVec(PInd);
L4II(PInd) =  C_II_Pixel_Us(PInd,:)*FrIPixVec - C_II_Pixel_Us(PInd,PInd)*FrIPixVec(PInd);
end
L4Eall = L4SE+L4CE+L4IE; L4Iall = L4SI+L4CI+L4II;
L4SEp = mean(L4SE./L4Eall);L4CEp = mean(L4CE./L4Eall);L4IEp = mean(L4IE./L4Eall);
L4SIp = mean(L4SI./L4Iall);L4CIp = mean(L4CI./L4Iall);L4IIp = mean(L4II./L4Iall);
LineFit = polyfit(L4Eall,L4Iall,1);
YRange = L4Iall - (LineFit(1)*L4Eall+LineFit(2)); % We plan to use 2 times of the range
Bdry = floor(max(abs(YRange))/100)*100;
L4ERange = floor(min(L4Eall)/100)*100:100:ceil(max(L4Eall)/100)*100;
L4IDiffRange = -2*Bdry:50:2*Bdry;
% Get a much larger domain

L4ERange = 0:200:ceil(max(L4Eall)*2/100)*100;
L4IDiffRange = -12*Bdry:200:12*Bdry;
%% Start parallel computation
cluster = parpool([4 128]);

a0 = length(L4ERange)*length(L4IDiffRange);
f_EnIOut = cell(a0,1);
meanVs = cell(a0,1);
SteadyIndicate = zeros(a0,1);
FailureIndicate = zeros(a0,1);
L4ERcrd = zeros(a0,1);
L4IRcrd = zeros(a0,1);

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
        
        if L4IU<0
            continue
        end
        %Distribute L4Input
        L4SEU = L4SEp*L4EU; L4CEU = L4CEp*L4EU; L4IEU = L4IEp*L4EU;
        L4SIU = L4SIp*L4IU; L4CIU = L4CIp*L4IU; L4IIU = L4IIp*L4IU;
        
       HyperPara = {'Traj',50,50,1,5000,'thre',0.2};
        tic
       [f_EnIOut{LDEInd},meanVs{LDEInd},~,...
        SteadyIndicate(LDEInd),FailureIndicate(LDEInd)]...
           = MFpV_SinglePixel(...% MF Parameters                     
                     N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,... %3 
                     S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                     S_EL6,S_IL6,rL6EU,rL6IU,S_amb,rE_amb,rI_amb,...%7 L6 Amb                                   
                     lgn_SU, lgn_COnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... %7
                     gL_E,gL_I,Ve,Vi, tau_ref,... %5
                     ...% Below are LIF details
                     tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,... %7
                     rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                     HyperPara);
        toc   
end

% Save MFpV Data
CurrentFolder = pwd;
SaveFolder = [CurrentFolder '/Data/LDE_Precomputing/'];
save([SaveFolder 'MFpV_LDE_thre' num2str(InputCtgr) '.mat'],...
    'f_EnIOut','meanVs','SteadyIndicate','FailureIndicate','L4ERcrd','L4IRcrd')   
