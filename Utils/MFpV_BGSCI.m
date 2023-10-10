%% mean-field stuff. We use parameters andto estimate firing rates
%% NOTE! mean Vs are now autonomus from single cell simulations
% Input: C_EE,C_EI,C_IE,C_II connectivity matrices
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        S_EL6,S_IL6,rE_L6,rI_L6 L6 Input
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rE_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulationsave('fig9MLMC.mat',fig9MLMC
%        varargin            >34: 'Mean' or 'Traj', indicate the form of the output
%                            >35: Sample Number after stopping criteria
%                            >36: Max Iteration before converged
%                            >37: Stepsize h
%                            >38: LIF simulation time
% Output:f_EnI               Estimation of firing rates, E and I; A sequences
%        meanVs              mean V of E and I
%        loop                Number of loops
%        SteadyIndicate      logical value for convergence

% Ver 2.0: New convergence conditions!
% Ver 3.0: We can now define 1. form of the output
%                            2. Sample number and max iteration number
%                            3. test different stepsizes and LIF simulation
%                            time
% Ver Test3.1: We try LSY's idea about using previous mV to stablize
% Ver Test3.3: To compensate the overstimate: add ref to MF estimate
% Ver 4:       L6 added
% iteration
function [f_EnIOut,meanVs,loop,SteadyIndicate,FailureIndicate] = ...
          MFpV_BGSCI(...
          N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,... %3 
                     S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                     S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                     lgn_S, lgn_C,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... %7
                     gL_E,gL_I,Ve,Vi, tau_ref,... %5
...% Below are LIF details
                     tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,... %7
                     rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                     HyperPara) % if varagin non empty, we only export the last state
%% Hyperparameters
N_Hyp = length(HyperPara);
if N_Hyp > 1  % Specify: The number of loops I want, after stopping criteria met    
    AveLoop = HyperPara{2};
else 
    AveLoop = 100;
end

if N_Hyp > 2 % Specify: Maximum nubmer of loops before stopping
        StopLoop = HyperPara{3};
else 
        StopLoop = AveLoop;
end

if N_Hyp > 3 % Specify: stepsize h
    h_Step = HyperPara{4};
else
    h_Step = 1;        
end

if N_Hyp > 4 % Specify: LIF simulation timef_pre
    LIFSimuT = HyperPara{5};
else
    LIFSimuT = 10*1e3;    % unit in ms    
end

if N_Hyp > 5 % Specify: Starting point of Recording
    if strcmpi(HyperPara{6},'thre')
        RecdThre = HyperPara{7};
        RecdDely = 0;
    elseif strcmpi(HyperPara{6},'delay')
        RecdThre = 0;
        RecdDely = HyperPara{7};
    else
        disp('***Wrong test mode')
        return
    end
else
    RecdThre = 0;    % unit in ms  
    RecdDely = 0;
end
dt = 0.1;
% initialize with a reasonable guess     
mVS = 0.59; 
mVC = 0.67; 
mVI = 0.67; 
f_pre = ones(3,1);
meanVs = [mVS;mVC;mVI];
f_SCIIni = MF_SCI_1Pix(N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                       S_EE,S_EI,S_IE,S_II,p_EEFail,...
                       S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                       lgn_S, lgn_C,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... %7
                       gL_E,gL_I,Ve,Vi, tau_ref, meanVs(:,end), f_pre,f_pre);
f_EnIOut = f_SCIIni;

loop = 0;
SteadyCounter = 0; % Indicate the number of loops after steady condition
SteadyIndicate = false;
TestPoints = floor(7); % How many consecutive points we test
%while( norm([mVEpre;mVIpre] - [mVE;mVI]) > 0.01 || norm(f_EnIpre - f_EnI0)>0.1) %%% relative difference for firing rates!!

FailureIndicate = 0;
Suspicious = false;
while SteadyCounter<AveLoop %the formal ending condition
% simulate one neuron with input 
if max(f_EnIOut(:,end))>501
    disp('Large Fr >501Hz. Break')
    FailureIndicate = 1;
    break 
end
mVLIF = zeros(3,1); FrLIF = zeros(3,1);
[mVLIF(1),FrLIF(1)] = MEanFieldEst_SingleCell_L6(...
                      's', f_EnIOut(:,end), ...
                      N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                      S_EE,S_EI,S_IE,S_II,p_EEFail,...
                      S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                      lgn_S, lgn_C,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,...
                      tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                      rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                      gL_E,gL_I,Ve,Vi,LIFSimuT,dt, RecdThre, RecdDely);
[mVLIF(2),FrLIF(2)] = MEanFieldEst_SingleCell_L6(...
                      'c', f_EnIOut(:,end), ...
                      N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                      S_EE,S_EI,S_IE,S_II,p_EEFail,...
                      S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                      lgn_S, lgn_C,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,...
                      tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                      rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                      gL_E,gL_I,Ve,Vi,LIFSimuT,dt, RecdThre, RecdDely);                  
[mVLIF(3),FrLIF(3)] = MEanFieldEst_SingleCell_L6(...
                      'i', f_EnIOut(:,end), ...
                      N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                      S_EE,S_EI,S_IE,S_II,p_EEFail,...
                      S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                      lgn_S, lgn_C,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,...
                      tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                      rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                      gL_E,gL_I,Ve,Vi,LIFSimuT,dt, RecdThre, RecdDely);  
                 
Fail = sum(isnan(mVLIF))>0; 
if Fail
    disp('NaN appears in mV. Break')
    FailureIndicate = 1;
    break    
end                             
%% NOW! consider the previous mVs if it already satisfies steady condition
if loop>StopLoop
    mVIn = mean(meanVs(:,end-10+1:end),2) * 0.9 + mVLIF*0.1;
else
    mVIn = mVLIF;
end

% estimate with ref now
f_EnI0 = MF_SCI_1Pix(N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                     S_EE,S_EI,S_IE,S_II,p_EEFail,...
                     S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                     lgn_S, lgn_C,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... % lgn
                     gL_E,gL_I,Ve,Vi, tau_ref, mVIn, f_SCIIni, FrLIF); 
                               
                               
Suspicious = (min(f_EnI0)<-1);                              
%f_EnI0 = max([f_EnI0,[0;0]],[],2);
%f_EnI0 = abs(f_EnI0);                                     
%% The new input!
f_SCIIni = f_EnI0*h_Step + f_SCIIni*(1-h_Step);

loop = loop+1;

f_EnIOut = [f_EnIOut,f_EnI0];
meanVs = [meanVs, mVIn];

if ~SteadyIndicate
  if loop >= TestPoints && mean(std(f_EnIOut(:,end-TestPoints+1:end),0,2) ...
                            ./mean(f_EnIOut(:,end-TestPoints+1:end),2))<0.07 % std/mean for the last 10 samples
     SteadyIndicate = true;
  end  
else 
    SteadyCounter = SteadyCounter+1;
end

if ((loop>StopLoop || SteadyIndicate) && Suspicious) 
    FailureIndicate = 1;
end

% break out if not reaching convergence after 100 iterations. Tbis number
% should be larger than end condition of SteadyCounter
if (~SteadyIndicate && loop >= StopLoop) 
    disp('Firing rates unconverged')
    break
end
end

if N_Hyp > 0 && strcmpi(HyperPara(1),'Mean') % Specify: I don't need the trajectory but the final ones    
    f_EnIOut = mean(f_EnIOut(:,end-50:end),2);
    meanVs   = mean(meanVs(:,end-50:end),2);
end
end



%% mean-field est With ref. We use parameters and mean Vs to estimate firing rates
% Input: N_PreSynPix         connectivity matrix for SCI within pixel. 3*3
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        S_EL6,S_IL6,rE_L6,rI_L6 L6 Input
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rS_amb,rC_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
%        tau_ref             ref period, in ms
%        f_pre               from previous
% Output:Fr_MFinv            Estimation of firing rates, Son, Soff, Con, Coff, I

function Fr_MFinv = MF_SCI_1Pix(N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                             S_EE,S_EI,S_IE,S_II,p_EEFail,...
                             S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                             lgn_S, lgn_C,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... % lgn
                             gL_E,gL_I,Ve,Vi, tau_ref, meanVs, f_pre, FrLIF)                             
%% PreSynaptic neuron numbers. External and internal (pixel)
N_SSi = N_PreSynPix(1,1); N_CSi = N_PreSynPix(1,2); N_ISi = N_PreSynPix(1,3); %% REally?? Row/Col
N_SCi = N_PreSynPix(2,1); N_CCi = N_PreSynPix(2,2); N_ICi = N_PreSynPix(2,3); 
N_SIi = N_PreSynPix(3,1); N_CIi = N_PreSynPix(3,2); N_IIi = N_PreSynPix(3,3); 

% Downplay current by a ref factor
% Use LIF firing rates if negative value emerges
f_preU = f_pre; 
f_preU(f_pre<0) = FrLIF(f_pre<0);
%f_preU(f_preU>200) = 200; % give a maximum
RefM = diag(1-f_preU*tau_ref/1000);
% Mats
MatSS = (S_EE*(1-p_EEFail))*N_SSi; MatCS = (S_EE*(1-p_EEFail))*N_CSi; MatIS = S_IE * N_ISi;
MatSC = (S_EE*(1-p_EEFail))*N_SCi; MatCC = (S_EE*(1-p_EEFail))*N_CCi; MatIC = S_IE * N_ICi;
MatSI = S_EI *              N_SIi; MatCI = S_EI *              N_CIi; MatII = S_II * N_IIi;

eVS  = (Ve-meanVs(1));  iVS  = (Vi-meanVs(1));
eVC  = (Ve-meanVs(2));  iVC  = (Vi-meanVs(2));
eVI  = (Ve-meanVs(3));  iVI  = (Vi-meanVs(3));
% post/pre SOn            COn              I
ConnMat = [MatSS.*eVS,  MatSC.*eVS,  MatSI.*iVS;  % SOn
           MatCS.*eVC,  MatCC.*eVC,  MatCI.*iVC;  % COn
           MatIS.*eVI,  MatIC.*eVI,  MatII.*iVI];   % I
% Leak On/Off
LeakS  = gL_E * (0-meanVs(1)) * 1e3;
LeakC  = gL_E * (0-meanVs(2)) * 1e3;
LeakI =  gL_I * (0-meanVs(3)) * 1e3;
LeakV = [LeakS;
         LeakC;
         LeakI];
% Ext
ExtS = (lgn_S*NlgnS*S_Elgn + rS_amb*S_amb + rL6SU*S_EL6)*eVS *1e3 ...
     + (L4SIU*S_EI)*iVS + L4SEU*S_EE*eVS*(1-p_EEFail); 
ExtC = (lgn_C*NlgnC*S_Elgn + rC_amb*S_amb + rL6CU*S_EL6)*eVC *1e3 ...
     + (L4CIU*S_EI)*iVC + L4CEU*S_EE*eVC*(1-p_EEFail);
ExtI = (lgn_I*NlgnI*S_Ilgn + rI_amb*S_amb + rL6IU*S_IL6)*eVI *1e3...
     + (L4IIU*S_II)*iVI + L4IEU*S_IE*eVI;
ExtV = [ExtS;
        ExtC;
        ExtI];
    
Fr_MFinv = (eye(3)-RefM*ConnMat) \ (RefM * ( ExtV + LeakV)); 
%Fr_MFinv(Fr_MFinv<0) = 3;
% FrPreUseDim = Fr_MFinv<0; 
% FrPreUse = FrLIF(Fr_MFinv<0);
% if isempty(FrPreUse)
%     return
% else
%     Fr_MFinv(FrPreUseDim) = FrPreUse;
% %     FrPreVec = zeros(5,1); FrPreVec(FrPreUseDim) = FrPreUse;
% %     DimAll = 1:5; FrminusDim = DimAll(~ismember(DimAll,FrPreUseDim));
% %     FrPreReplace = RefM*ConnMat * FrPreVec;
% %     % Now compute in the subspace
% %     FrPreMinus = FrPreReplace(FrminusDim);
% %     RefMMinus  = RefM(FrminusDim, FrminusDim);
% %     ConnMatMinus = ConnMat(FrminusDim, FrminusDim);
% %     ExtVMinus = ExtV(FrminusDim);
% %     LeakVMinus = LeakV(FrminusDim);
% %     Fr_MFinvMinus = (eye(length(FrminusDim))-RefMMinus*ConnMatMinus) \ ...
% %         (RefMMinus * ( ExtVMinus + LeakVMinus) + FrPreMinus);
% %     % make up MF output
% %     Fr_MFinv = zeros(5,1);
% %     Fr_MFinv(FrminusDim) = Fr_MFinvMinus; Fr_MFinv(FrPreUseDim) = FrPreUse;
% end

end
%% single cell simulation: to collect mean V (and firing rates)


function [meanV,fr] = MEanFieldEst_SingleCell_L6(...
                      NeuronType, Fr_MFinv, ...
                      N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                      S_EE,S_EI,S_IE,S_II,p_EEFail,...
                      S_EL6,S_IL6,rS_L6,rC_L6,rI_L6,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                      lgn_S, lgn_C,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,...
                      tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                      rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                      gL_E,gL_I,Ve,Vi,LIFSimuT,dt, RecdThre, RecdDely)
%% attribute parameters for different types of neurons
if strcmpi(NeuronType,'s')
    N_S = N_PreSynPix(1,1);N_C = N_PreSynPix(2,1);N_I = N_PreSynPix(3,1);
    L4E = L4SEU; L4I = L4SIU;
    S_E = S_EE; S_I = S_EI;
    p_Fail = p_EEFail;
    lambda = lgn_S*NlgnS; S_lgn = S_Elgn; r_amb = rS_amb;
    S_L6 = S_EL6; r_L6 = rS_L6;
    gL = gL_E;
    rho_ampa = rhoE_ampa; rho_nmda = rhoE_nmda;
elseif strcmpi(NeuronType,'c')
    N_S = N_PreSynPix(1,2);N_C = N_PreSynPix(2,2);N_I = N_PreSynPix(3,2);
    L4E = L4CEU; L4I = L4CIU;
    S_E = S_EE; S_I = S_EI;
    p_Fail = p_EEFail;
    lambda = lgn_C*NlgnC; S_lgn = S_Elgn; r_amb = rC_amb;
    S_L6 = S_EL6; r_L6 = rC_L6;
    gL = gL_E;
    rho_ampa = rhoE_ampa; rho_nmda = rhoE_nmda;
elseif strcmpi(NeuronType,'i')
    N_S = N_PreSynPix(1,3);N_C = N_PreSynPix(2,3);N_I = N_PreSynPix(3,3);
    L4E = L4IEU; L4I = L4IIU;
    S_E = S_IE; S_I = S_II;
    p_Fail = 0;
    lambda = lgn_I*NlgnI; S_lgn = S_Ilgn; r_amb = rI_amb;
    S_L6 = S_IL6; r_L6 = rI_L6;
    gL = gL_I;
    rho_ampa = rhoI_ampa; rho_nmda = rhoI_nmda;
else 
    disp('***Unrecognized Neuron Type')    
end
if length(Fr_MFinv) == 3 % S C I
    Fr_MFinv(Fr_MFinv<0) = eps;
    rS = Fr_MFinv(1); 
    rC = Fr_MFinv(2); 
    rI = Fr_MFinv(3); % f_EnI in s^-1, but here we use ms^-1
else
    error('L4 Pix FRs dont match!')
end

TotalE = (L4E + rS*N_S + rC*N_C)*(1-p_Fail)/1000; 
TotalI = L4I + rI*N_I/1000; % f_EnI in s^-1, but here we use ms^-1
%% Evolve single neurons
T = LIFSimuT; % in ms. Default: 20*1e3 ms
t = 0:dt:T;
SampleProp = 9/10; % last half time for meanV

v = zeros(size(t)); 
G_gaba_D = 0; G_gaba_R = 0;
G_ampa_D = 0; G_ampa_R = 0;
G_nmda_D = 0; G_nmda_R = 0;
spike = [];

%rng(100)
% input determination: Assume all Poisson
%rng(100)
p_lgn = dt*lambda; Sp_lgn = double(rand(size(t))<=p_lgn);
%rng(101)
p_amb = dt*r_amb;  Sp_amb = double(rand(size(t))<=p_amb);
%rng(102)
p_EV1 = dt*TotalE;     Sp_EV1 = double(rand(size(t))<=p_EV1);
%rng(103)
p_IV1 = dt*TotalI;     Sp_IV1 = double(rand(size(t))<=p_IV1);
p_L6  = dt*r_L6;   Sp_L6  = double(rand(size(t))<=p_L6);

RefTimer = 0; 
for tInd = 1:length(t)-1
    % Firstly, refrectory neurons get out due to exponetial distributed time
     if isnan(v(tInd))
        RefTimer = RefTimer+dt; %RefTimer goes up
        if RefTimer>=tau_ref    %if timer reach tau_ref, kick v out of refrectory
            v(tInd+1) = 0;
            RefTimer = 0;
        else
            v(tInd+1) = nan;
        end
     else
         G_I = 1/(tau_gaba_D-tau_gaba_R) * (G_gaba_D - G_gaba_R); % S_EI is included in amplitude of GE_gaba
         G_E = 1/(tau_ampa_D-tau_ampa_R) * (G_ampa_D - G_ampa_R) ...
             + 1/(tau_nmda_D-tau_nmda_R) * (G_nmda_D - G_nmda_R); % 
         vv  = v(tInd) + dt*(-gL*v(tInd) - G_E.*(v(tInd)-Ve) - G_I.*(v(tInd)-Vi)); 
         if vv >= 1
             v(tInd+1) = nan;
             spike = [spike,t(tInd)];
         else
             v(tInd+1) = vv;
         end
     end

     % conductances
     G_gaba_Ro = (G_gaba_R +                                           S_I*Sp_IV1(tInd)                                     ) * exp(-dt/tau_gaba_R); 
     G_gaba_Do = (G_gaba_D +                                           S_I*Sp_IV1(tInd)                                     ) * exp(-dt/tau_gaba_D); 
     G_ampa_Ro = (G_ampa_R + S_lgn*Sp_lgn(tInd) + S_amb*Sp_amb(tInd) + S_E*Sp_EV1(tInd)*rho_ampa + S_L6*Sp_L6(tInd)*rho_ampa) * exp(-dt/tau_ampa_R);
     G_ampa_Do = (G_ampa_D + S_lgn*Sp_lgn(tInd) + S_amb*Sp_amb(tInd) + S_E*Sp_EV1(tInd)*rho_ampa + S_L6*Sp_L6(tInd)*rho_ampa) * exp(-dt/tau_ampa_D);
     G_nmda_Ro = (G_nmda_R +                                           S_E*Sp_EV1(tInd)*rho_nmda + S_L6*Sp_L6(tInd)*rho_nmda) * exp(-dt/tau_nmda_R);
     G_nmda_Do = (G_nmda_D +                                           S_E*Sp_EV1(tInd)*rho_nmda + S_L6*Sp_L6(tInd)*rho_nmda) * exp(-dt/tau_nmda_D);

     G_gaba_D = G_gaba_Do; G_gaba_R = G_gaba_Ro;
     G_ampa_D = G_ampa_Do; G_ampa_R = G_ampa_Ro;
     G_nmda_D = G_nmda_Do; G_nmda_R = G_nmda_Ro;
end

if RecdThre>0 && ~strcmpi(NeuronType,'s')
    v(v< RecdThre & v >0) = nan;
end
GridRef = find(isnan(v));
GridDly = floor(RecdDely/dt);
% or recording delay
if ~isempty(GridRef) && GridDly>0 && ~strcmpi(NeuronType,'s')% Only do more nan if spikes and nontrivial dly    
    v(unique(reshape(GridRef + (0: GridDly)',1,length(GridRef)*(GridDly+1)))) = nan;
end
meanV = mean(v(floor(end*(1-SampleProp)):end), 'omitnan');
fr = length(find(spike>(1-SampleProp)*T))/(T*SampleProp/1000);
end
