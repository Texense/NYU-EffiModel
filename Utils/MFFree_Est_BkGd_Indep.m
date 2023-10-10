%% mean-field stuff. We use parameters andto estimate firing rates
%% NOTE! mean Vs are now autonomus from single cell simulations
% Input: C_EE,C_EI,C_IE,C_II connectivity matrices
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rE_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
% Output:f_EnI               Estimation of firing rates, E and I; A sequences
%        meanVs              mean V of E and I
function [f_EnIOut,meanVs,loop] = MFFree_Est_BkGd_Indep(C_EE,C_EI,C_IE,C_II,... %4
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                                   lambda_E,S_Elgn,rE_amb,S_amb,... %4
                                   lambda_I,S_Ilgn,rI_amb,... %3
                                   tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,... %7
                                   rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                                   gL_E,gL_I,Ve,Vi,... %4
                                   N_HC,n_E_HC,n_I_HC,varargin) %3+x % if varagin non empty, we only export the last state
                                   %This line: we are taking spatial center 
% First define fr and mV
f_EnIOut = [];

% Get connectivity for the center HC
E_sideInd = floor(1*n_E_HC+1):2*n_E_HC;
[E_Ind_X,E_Ind_Y] = meshgrid(E_sideInd,E_sideInd);
E_Ind = (reshape(E_Ind_X,size(E_Ind_X,1)*size(E_Ind_X,2),1)-1)*n_E_HC*N_HC + reshape(E_Ind_Y,size(E_Ind_X,1)*size(E_Ind_X,2),1);

I_sideInd = floor(1*n_I_HC+1):2*n_I_HC;
[I_Ind_X,I_Ind_Y] = meshgrid(I_sideInd,I_sideInd);
I_Ind = (reshape(I_Ind_X,size(I_Ind_X,1)*size(I_Ind_X,2),1)-1)*n_I_HC*N_HC + reshape(I_Ind_Y,size(I_Ind_X,1)*size(I_Ind_X,2),1);

% Take averaged number of input neurons 
% NOTE! Only picking up the middle part
N_EE = mean(sum(C_EE(E_Ind,:),2)); N_EI = mean(sum(C_EI(E_Ind,:),2)); 
N_IE = mean(sum(C_IE(I_Ind,:),2)); N_II = mean(sum(C_II(I_Ind,:),2));

% initialize with a reasonable guess     
mVE = 0.57; mVI = 0.67; 
%mVE =  1; mVI = 1; 
meanVs = [mVE;mVI];
mVEpre  = -1; mVIpre = -1;
f_EnI0 = [3.5;12]; f_EnIpre = [0;0]; % start with an impossible value
N_Neuron = 40; % Number of Neurons we use
InitialStatesE = struct('v',zeros(1,N_Neuron),'G_gaba_D',zeros(1,N_Neuron),'G_gaba_R',zeros(1,N_Neuron),...
                                              'G_ampa_D',zeros(1,N_Neuron),'G_ampa_R',zeros(1,N_Neuron),...
                                              'G_nmda_D',zeros(1,N_Neuron),'G_nmda_R',zeros(1,N_Neuron));
InitialStatesI = struct('v',zeros(1,N_Neuron),'G_gaba_D',zeros(1,N_Neuron),'G_gaba_R',zeros(1,N_Neuron),...
                                              'G_ampa_D',zeros(1,N_Neuron),'G_ampa_R',zeros(1,N_Neuron),...
                                              'G_nmda_D',zeros(1,N_Neuron),'G_nmda_R',zeros(1,N_Neuron));

loop = 0;
while( norm([mVEpre;mVIpre] - [mVE;mVI]) > 0.01 || norm((f_EnIpre - f_EnI0)./f_EnI0)>0.1) %%% relative difference for firing rates!!
    tic
mVEpre = mVE; mVIpre = mVI;
f_EnINow = zeros(size(f_EnI0)); % from simulation
f_EnIpre = f_EnI0;
% simulate nultiple neurons   
[mVE,f_EnINow(1),EndStatesE] = MEanFieldEst_MultiCell('e', f_EnI0, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi,InitialStatesE,N_Neuron);
[mVI,f_EnINow(2),EndStatesI] = MEanFieldEst_MultiCell('i', f_EnI0, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi,InitialStatesI,N_Neuron);

loop = loop+1;

f_EnI0 = f_EnINow
f_EnIOut = [f_EnIOut,f_EnINow];
meanVs = [meanVs, [mVE;mVI]];

InitialStatesE = EndStatesE; InitialStatesI = EndStatesI; % Now if forms a loop
 toc
end

if nargin > 34 % Specify: I don't need the trajectory but the final ones
    f_EnIOut = f_EnIOut(:,end);
    meanVs = meanVs(:,end);
end
end


%% single cell simulation: to collect mean V (and firing rates)
% Input: NeuronType           'e' or 'i'
%        InitialStates        V and all Gs
%        NeuronNum            Number of NeuronSimulated. Generally we don't
%        want this number too large.

% Output: meanV, fr           Collected from all neurons, and time
%         EndStates           V and all Gs, for the next iteration

function [meanV,fr,EndStates] = MEanFieldEst_MultiCell(NeuronType, f_EnI, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi,InitialStates,NeuronNum)
%% attribute parameters for different types of neurons
if strcmpi(NeuronType,'e')
    N_E = N_EE; N_I = N_EI;
    S_E = S_EE; S_I = S_EI;
    p_Fail = p_EEFail;
    lambda = lambda_E; S_lgn = S_Elgn; r_amb = rE_amb;
    gL = gL_E;
    rho_ampa = rhoE_ampa; rho_nmda = rhoE_nmda;
elseif strcmpi(NeuronType,'i')
    N_E = N_IE; N_I = N_II;
    S_E = S_IE; S_I = S_II;
    p_Fail = 0;
    lambda = lambda_I; S_lgn = S_Ilgn; r_amb = rI_amb;
    gL = gL_I;
    rho_ampa = rhoI_ampa; rho_nmda = rhoI_nmda;
else 
    disp('***Unrecognized Neuron Type')    
end
rE = f_EnI(1)/1000; rI = f_EnI(2)/1000; % f_EnI in s^-1, but here we use ms^-1
%% Evolve single neurons
T = 15*1e3; % in ms
dt = 0.01; t = 0:dt:T;
SampleProp = 9/10; % last half time for meanV

v = zeros(length(t),NeuronNum); % time by rows and neuronInd by columns 
G_gaba_D = zeros(length(t),NeuronNum); G_gaba_R = zeros(length(t),NeuronNum);
G_ampa_D = zeros(length(t),NeuronNum); G_ampa_R = zeros(length(t),NeuronNum);
G_nmda_D = zeros(length(t),NeuronNum); G_nmda_R = zeros(length(t),NeuronNum);
spike = [];
G_I = zeros(length(t),NeuronNum); G_E = zeros(length(t),NeuronNum);

v(1,:) = InitialStates.v;
G_gaba_D(1,:)  = InitialStates.G_gaba_D; G_gaba_R(1,:)  = InitialStates.G_gaba_R;
G_ampa_D(1,:)  = InitialStates.G_ampa_D; G_ampa_R(1,:)  = InitialStates.G_ampa_R;
G_nmda_D(1,:)  = InitialStates.G_nmda_D; G_nmda_R(1,:)  = InitialStates.G_nmda_R;
%rng(100)
% input determination: Assume all Poisson
%rng(100)
p_lgn = 1-exp(-dt*lambda);            Sp_lgn = double(rand(size(v))<=p_lgn);
%rng(101)
p_amb = 1-exp(-dt*r_amb);             Sp_amb = double(rand(size(v))<=p_amb);
%rng(102)
p_EV1 = dt*rE*full(N_E)*(1-p_Fail);   Sp_EV1 = double(rand(size(v))<=p_EV1); % preserving expectation correct
%rng(103)
p_IV1 = dt*rI*full(N_I);              Sp_IV1 = double(rand(size(v))<=p_IV1); % preserving expectation correct

%% NOTE!! I am cheating a little here: RefTimer always restart from 0.
RefTimer = zeros(size(InitialStates.v));
for tInd = 1:length(t)-1
     % Taking all previous states, compute new v as is (defer ref decisions later)
     G_I(tInd,:) = 1/(tau_gaba_D-tau_gaba_R) * (G_gaba_D(tInd,:) - G_gaba_R(tInd,:)); % S_EI is included in amplitude of GE_gaba
     G_E(tInd,:) = 1/(tau_ampa_D-tau_ampa_R) * (G_ampa_D(tInd,:) - G_ampa_R(tInd,:)) ...
           + 1/(tau_nmda_D-tau_nmda_R) * (G_nmda_D(tInd,:) - G_nmda_R(tInd,:)); % 
     vv = v(tInd,:) + dt*(-gL*v(tInd,:) - G_E(tInd,:) .*(v(tInd,:)-Ve) - G_I(tInd,:).*(v(tInd,:)-Vi)); %Note: This ensures [nan v] gives [nan vv]
     spike = [spike,tInd*dt*ones(1,sum(vv>=1,'all'))];
     vv(vv>=1) = nan; 
     
     %Now, refrectory neurons get out due to exponetial distributed time    
     RefTimer(isnan(v(tInd,:))) = RefTimer(isnan(v(tInd,:)))+dt;%RefTimer goes up for those in ref period
     vv(RefTimer>=tau_ref) = 0;%if timer reach tau_ref, kick v out of refrectory
     RefTimer(RefTimer>=tau_ref) = 0; %Restart        
     v(tInd+1,:) = vv;     

     % conductances
     G_gaba_R(tInd+1,:) = (G_gaba_R(tInd,:)) * exp(-dt/tau_gaba_R) + S_I*Sp_IV1(tInd,:); 
     G_gaba_D(tInd+1,:) = (G_gaba_D(tInd,:)) * exp(-dt/tau_gaba_D) + S_I*Sp_IV1(tInd,:); 
     G_ampa_R(tInd+1,:) = (G_ampa_R(tInd,:)) * exp(-dt/tau_ampa_R) + S_lgn * Sp_lgn(tInd,:) + S_amb * Sp_amb(tInd,:) + S_E * Sp_EV1(tInd,:) * rho_ampa;
     G_ampa_D(tInd+1,:) = (G_ampa_D(tInd,:)) * exp(-dt/tau_ampa_D) + S_lgn * Sp_lgn(tInd,:) + S_amb * Sp_amb(tInd,:) + S_E * Sp_EV1(tInd,:) * rho_ampa;
     G_nmda_R(tInd+1,:) = (G_nmda_R(tInd,:)) * exp(-dt/tau_nmda_R) + S_E * Sp_EV1(tInd,:) * rho_nmda;
     G_nmda_D(tInd+1,:) = (G_nmda_D(tInd,:)) * exp(-dt/tau_nmda_D) + S_E * Sp_EV1(tInd,:) * rho_nmda;
end

meanV = nanmean(v(floor(end*(1-SampleProp)):end,:),'all');
fr = length(spike)/(T/1e3)/NeuronNum;
EndStates = struct('v',v(end,:),'G_gaba_D',G_gaba_D(end,:),'G_gaba_R',G_gaba_R(end,:),...
                                'G_ampa_D',G_ampa_D(end,:),'G_ampa_R',G_ampa_R(end,:),...
                                'G_nmda_D',G_nmda_D(end,:),'G_nmda_R',G_nmda_R(end,:));
end