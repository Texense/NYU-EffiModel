CurrentFolder = pwd
FigurePath = [CurrentFolder '/Figures'];
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
SaveFolder = [CurrentFolder '/Figures/Demo082721/'];
addpath(SaveFolder)
DataFolder = [CurrentFolder '/Figures/Demo082121/'];
addpath(DataFolder)

load('AllMFPixPara.mat')

%% Start parallel computation
cluster = parpool([4 128]);
vRcrdThre = 0.05;

mVLIF = zeros(5,PixNum); FrLIF = zeros(5,PixNum);
f_EnIOut = cell(PixNum,1);
meanVs = cell(PixNum,1);
SteadyIndicate = zeros(PixNum,1);
FailureIndicate = zeros(PixNum,1);
parfor  PInd = 1:900
        VSonP = mVSOnPixVec(PInd); VSoffP = mVSOffPixVec(PInd);
        VConP = mVCOnPixVec(PInd); VCoffP = mVCOffPixVec(PInd);
        VIP   = mVIPixVec(PInd);
        %meanVs = [VSonP;VConP;VSoffP;VCoffP;VIP]
        f_pre = [FrSOnPixVec(PInd) ;
            FrCOnPixVec(PInd);
            FrSOffPixVec(PInd);
            FrCOffPixVec(PInd);
            FrIPixVec(PInd)];
        % Input
        NlgnS = N_Slgn; NlgnC = N_Clgn; NlgnI = N_Ilgn;
        rL6E = (L6S_Pixel(PInd) + L6C_Pixel(PInd))/2; rL6I = L6I_Pixel(PInd);
        lgn_SOnOff = [lambda_SOn_Pixel(PInd)/NlgnS;
                     lambda_SOff_Pixel(PInd)/NlgnS];
        lgn_COnOff = [lambda_COn_Pixel(PInd)/NlgnC;
                     lambda_COff_Pixel(PInd)/NlgnC];
        lgn_I =      lambda_I_Pixel/NlgnI;
        % L4 Input and parameters
        N_PreSynPix = [C_SS_Pixel_Us(PInd, PInd),C_CS_Pixel_Us(PInd, PInd),C_IS_Pixel_Us(PInd, PInd);
                       C_SC_Pixel_Us(PInd, PInd),C_CC_Pixel_Us(PInd, PInd),C_IC_Pixel_Us(PInd, PInd);
                       C_SI_Pixel_Us(PInd, PInd),C_CI_Pixel_Us(PInd, PInd),C_II_Pixel_Us(PInd, PInd)];
%         N_PreSynPix(1,1) = C_SS_Pixel_Us(PInd, PInd);%mean(diag(C_SS_Pixel_Us)); %;
%         N_PreSynPix(1,2) = C_CS_Pixel_Us(PInd, PInd);%mean(diag(C_CS_Pixel_Us)); %;
%         N_PreSynPix(1,3) = C_IS_Pixel_Us(PInd, PInd);%mean(diag(C_IS_Pixel_Us)); %;
%         N_PreSynPix(2,1) = C_SC_Pixel_Us(PInd, PInd);%mean(diag(C_SC_Pixel_Us)); %;
%         N_PreSynPix(2,2) = C_CC_Pixel_Us(PInd, PInd);%mean(diag(C_CC_Pixel_Us)); %;
%         N_PreSynPix(2,3) = C_IC_Pixel_Us(PInd, PInd);%mean(diag(C_IC_Pixel_Us)); %;
%         N_PreSynPix(3,1) = C_SI_Pixel_Us(PInd, PInd);%mean(diag(C_SI_Pixel_Us)); %;
%         N_PreSynPix(3,2) = C_CI_Pixel_Us(PInd, PInd);%mean(diag(C_CI_Pixel_Us)); %;
%         N_PreSynPix(3,3) = C_II_Pixel_Us(PInd, PInd);%mean(diag(C_II_Pixel_Us)); %;
        
        L4SE = (C_SS_Pixel_Us(PInd,:)*FrSPixVec          + C_SC_Pixel_Us(PInd,:)*FrCPixVec)...
            - C_SS_Pixel_Us(PInd,PInd)*FrSPixVec(PInd) - C_SC_Pixel_Us(PInd,PInd)*FrCPixVec(PInd);
        L4SI = C_SI_Pixel_Us(PInd,:)*FrIPixVec - C_SI_Pixel_Us(PInd,PInd)*FrIPixVec(PInd);
        L4CE = (C_CS_Pixel_Us(PInd,:)*FrSPixVec          + C_CC_Pixel_Us(PInd,:)*FrCPixVec)...
            - C_CS_Pixel_Us(PInd,PInd)*FrSPixVec(PInd) - C_CC_Pixel_Us(PInd,PInd)*FrCPixVec(PInd);
        L4CI = C_CI_Pixel_Us(PInd,:)*FrIPixVec - C_CI_Pixel_Us(PInd,PInd)*FrIPixVec(PInd);
        L4IE = (C_IS_Pixel_Us(PInd,:)*FrSPixVec          + C_IC_Pixel_Us(PInd,:)*FrCPixVec)...
            - C_IS_Pixel_Us(PInd,PInd)*FrSPixVec(PInd) - C_IC_Pixel_Us(PInd,PInd)*FrCPixVec(PInd);
        L4II = C_II_Pixel_Us(PInd,:)*FrIPixVec - C_II_Pixel_Us(PInd,PInd)*FrIPixVec(PInd);
        
        LIFSimuT = 10000;
        Fr_MFinv = f_pre;
        [mVLIF(:,PInd),FrLIF(:,PInd)] = LIF1Pixel(Fr_MFinv, N_PreSynPix, L4SE,L4SI, L4CE,L4CI, L4IE,L4II,...
                                 S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                 S_EL6,S_IL6,rL6E,rL6I,S_amb,rE_amb,rI_amb,...%7 L6 Amb                                   
                                 lgn_SOnOff,lgn_COnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,...
                                 tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                 rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                 gL_E,gL_I,Ve,Vi,LIFSimuT, dt, vRcrdThre);

        HyperPara = {'Traj',50,50,1,5000};
        tic
       [f_EnIOut{PInd},meanVs{PInd},loop,SteadyIndicate(PInd),FailureIndicate(PInd)]...
           = MFpV_SinglePixel(...
...% MF Parameters                     
                     N_PreSynPix, L4SE,L4SI, L4CE,L4CI, L4IE,L4II,... %3 
                     S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                     S_EL6,S_IL6,rL6E,rL6I,S_amb,rE_amb,rI_amb,...%7 L6 Amb                                   
                     lgn_SOnOff, lgn_COnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... %7
                     gL_E,gL_I,Ve,Vi, tau_ref,... %5
...% Below are LIF details
                     tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,... %7
                     rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                     HyperPara);
        toc
    
end
CurrentFolder = pwd;
SaveFolder = [CurrentFolder '/Figures/Demo082721/'];
save([SaveFolder 'MFpVPixVth' num2str(vRcrdThre) '.mat'],...
    'mVLIF','FrLIF','f_EnIOut','meanVs','SteadyIndicate','FailureIndicate')   