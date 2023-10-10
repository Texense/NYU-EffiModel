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
function [f_EnIOut,meanVs,loop,SteadyIndicate,FailureIndicate] = MeanFieldEst_BkGd_Indep_StepSize_ref_testL6_NoInv(C_EE,C_EI,C_IE,C_II,... %4
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                                   lambda_E,S_Elgn,rE_amb,S_amb,... %4
                                   lambda_I,S_Ilgn,rI_amb,... %3
                                   S_EL6,S_IL6,rE_L6,rI_L6,...%4 L6 added
                                   tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,... %7
                                   rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                                   gL_E,gL_I,Ve,Vi,... %4
                                   N_HC,n_E_HC,n_I_HC,varargin) %3+x % if varagin non empty, we only export the last state
                                   %This line: we are taking spatial center 

if nargin > 35+4  % Specify: The number of loops I want, after stopping criteria met    
    AveLoop = varargin{2};
else 
    AveLoop = 100;
end

if nargin > 36+4 % Specify: Maximum nubmer of loops before stopping
        StopLoop = varargin{3};
else 
        StopLoop = AveLoop;
end

if nargin > 37+4 % Specify: stepsize h
    h_Step = varargin{4};
else
    h_Step = 1;        
end

if nargin > 38+4 % Specify: LIF simulation timef_pre
    LIFSimuT = varargin{5};
else
    LIFSimuT = 20*1e3;    % unit in ms    
end


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
%mVE =  0.539; mVI = 0.683; 
%mVE = 0.5; mVI = 0.5; 
meanVs = [mVE;mVI];
f_EnIIni = MeanFieldEst_BkGd_L6(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   S_EL6,S_IL6,rE_L6,rI_L6,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI);
f_EnIOut = f_EnIIni;

loop = 0;
SteadyCounter = 0; % Indicate the number of loops after steady condition
SteadyIndicate = false;
TestPoints = floor(15); % How many consecutive points we test
%while( norm([mVEpre;mVIpre] - [mVE;mVI]) > 0.01 || norm(f_EnIpre - f_EnI0)>0.1) %%% relative difference for firing rates!!

FailureIndicate = 0;
Suspicious = false;
while SteadyCounter<AveLoop %the formal ending condition
% simulate one neuron with input                             
[mVE,~] = MEanFieldEst_SingleCell_L6('e', f_EnIIni, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         S_EL6,S_IL6,rE_L6,rI_L6,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi,LIFSimuT);
[mVI,~] = MEanFieldEst_SingleCell_L6('i', f_EnIIni, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         S_EL6,S_IL6,rE_L6,rI_L6,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi,LIFSimuT);
                                     
%% NOW! consider the previous mVs if it already satisfies steady condition
if loop>100
    mVEIn = mean(meanVs(1,end-10+1:end)) * 0.9 + mVE*0.1;
    mVIIn = mean(meanVs(2,end-10+1:end)) * 0.9 + mVI*0.1;
else
     mVEIn = mVE;
     mVIIn = mVI;
end

% estimate with ref now
f_EnI0 = MeanFieldEst_BkGd_ref_L6_NoInv(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   S_EL6,S_IL6,rE_L6,rI_L6,...
                                   gL_E,gL_I,Ve,Vi,mVEIn,mVIIn,...
                                   tau_ref,f_EnIIni);

Suspicious = (min(f_EnI0)<0);
                               
%f_EnI0 = max([f_EnI0,[0;0]],[],2);
%f_EnI0 = abs(f_EnI0);                                     
%% The new input!
f_EnIIni = f_EnI0*h_Step + f_EnIIni*(1-h_Step);

loop = loop+1;

f_EnIOut = [f_EnIOut,f_EnI0];
meanVs = [meanVs, [mVE;mVI]];

if ~SteadyIndicate
  if loop >= TestPoints && max(std(f_EnIOut(:,end-TestPoints+1:end),0,2) ...
                            ./mean(f_EnIOut(:,end-TestPoints+1:end),2))<0.05 % std/mean for the last 10 samples
     SteadyIndicate = true;
  end  
else 
    SteadyCounter = SteadyCounter+1;
end

if (loop>100 | SteadyIndicate) & Suspicious
    FailureIndicate = 1;
end

% break out if not reaching convergence after 100 iterations. Tbis number
% should be larger than end condition of SteadyCounter
if (~SteadyIndicate && loop >= StopLoop) 
    disp('Firing rates unconverged')
    break
end
end

if nargin > 34+4 && strcmpi(varargin(1),'Mean') % Specify: I don't need the trajectory but the final ones    
    f_EnIOut = mean(f_EnIOut(:,end-50:end),2);
    meanVs   = mean(meanVs(:,end-50:end),2);
end
end



%% mean-field est With ref. We use parameters and mean Vs to estimate firing rates
% Input: C_EE,C_EI,C_IE,C_II connectivity matrices
%        N_EE,N_EI,N_IE,N_II number of presynaptic neurons
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        S_EL6,S_IL6,rE_L6,rI_L6 L6 Input
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rE_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
%        tau_ref             ref period, in ms
%        f_pre               from previous
% Output:f_EnI               Estimation of firing rates, E and I

function f_EnI = MeanFieldEst_BkGd_ref_L6_NoInv(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   S_EL6,S_IL6,rE_L6,rI_L6,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI,...
                                   tau_ref,f_pre)
                               
% E_sideInd = floor(1*n_E_HC+1):2*n_E_HC;
% [E_Ind_X,E_Ind_Y] = meshgrid(E_sideInd,E_sideInd);
% E_Ind = (reshape(E_Ind_X,size(E_Ind_X,1)*size(E_Ind_X,2),1)-1)*n_E_HC*N_HC + reshape(E_Ind_Y,size(E_Ind_X,1)*size(E_Ind_X,2),1);
% 
% I_sideInd = floor(1*n_I_HC+1):2*n_I_HC;
% [I_Ind_X,I_Ind_Y] = meshgrid(I_sideInd,I_sideInd);
% I_Ind = (reshape(I_Ind_X,size(I_Ind_X,1)*size(I_Ind_X,2),1)-1)*n_I_HC*N_HC + reshape(I_Ind_Y,size(I_Ind_X,1)*size(I_Ind_X,2),1);
% 
% % Take averaged number of input neurons 
% % NOTE! Only picking up the middle part
% N_EE = mean(sum(C_EE(E_Ind,:),2)); N_EI = mean(sum(C_EI(E_Ind,:),2)); 
% N_IE = mean(sum(C_IE(I_Ind,:),2)); N_II = mean(sum(C_II(I_Ind,:),2));

% Downplay current by a ref factor
ref_fac = 1-f_pre*tau_ref/1000;
% Matrix of firing rates
MatEI = [N_EE*S_EE*(Ve-mVE)*(1-p_EEFail)*ref_fac(1), N_EI*S_EI*(Vi-mVE)*ref_fac(1);
         N_IE*S_IE*(Ve-mVI)*ref_fac(2),              N_II*S_II*(Vi-mVI)*ref_fac(2)];

% External Drive and leaky current. Note that simulaton is by ms, so need *1000
VecInput = [(Ve-mVE) * (lambda_E*S_Elgn + rE_amb*S_amb + rE_L6*S_EL6);
            (Ve-mVI) * (lambda_I*S_Ilgn + rI_amb*S_amb + rI_L6*S_IL6)]*1000;
Leak = -[gL_E*mVE;
         gL_I*mVI]*1000; % leak has opposite signs of mean Vs

f_EnI = MatEI*f_pre + ((Leak+VecInput).*ref_fac);    
%f_EnI = -(ref_fac*MatEI-eye(2))^-1*ref_fac*(Leak+VecInput);

end

%% mean-field w/o ref. We use parameters and mean Vs to estimate firing rates
% Input: C_EE,C_EI,C_IE,C_II connectivity matrices
%        N_EE,N_EI,N_IE,N_II number of presynaptic neurons
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        S_EL6,S_IL6,rE_L6,rI_L6 L6 Input
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rE_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
% Output:f_EnI               Estimation of firing rates, E and I

function f_EnI = MeanFieldEst_BkGd_L6(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   S_EL6,S_IL6,rE_L6,rI_L6,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI)
                               
% E_sideInd = floor(1*n_E_HC+1):2*n_E_HC;
% [E_Ind_X,E_Ind_Y] = meshgrid(E_sideInd,E_sideInd);
% E_Ind = (reshape(E_Ind_X,size(E_Ind_X,1)*size(E_Ind_X,2),1)-1)*n_E_HC*N_HC + reshape(E_Ind_Y,size(E_Ind_X,1)*size(E_Ind_X,2),1);
% 
% I_sideInd = floor(1*n_I_HC+1):2*n_I_HC;
% [I_Ind_X,I_Ind_Y] = meshgrid(I_sideInd,I_sideInd);
% I_Ind = (reshape(I_Ind_X,size(I_Ind_X,1)*size(I_Ind_X,2),1)-1)*n_I_HC*N_HC + reshape(I_Ind_Y,size(I_Ind_X,1)*size(I_Ind_X,2),1);
% 
% % Take averaged number of input neurons 
% % NOTE! Only picking up the middle part
% N_EE = mean(sum(C_EE(E_Ind,:),2)); N_EI = mean(sum(C_EI(E_Ind,:),2)); 
% N_IE = mean(sum(C_IE(I_Ind,:),2)); N_II = mean(sum(C_II(I_Ind,:),2));

% Matrix of firing rates
MatEI = [N_EE*S_EE*(Ve-mVE)*(1-p_EEFail), N_EI*S_EI*(Vi-mVE);
         N_IE*S_IE*(Ve-mVI),              N_II*S_II*(Vi-mVI)];

% External Drive and leaky current. Note that simulaton is by ms, so need *1000
VecInput = [(Ve-mVE) * (lambda_E*S_Elgn + rE_amb*S_amb + rE_L6*S_EL6);
            (Ve-mVI) * (lambda_I*S_Ilgn + rI_amb*S_amb + rI_L6*S_IL6)]*1000;
Leak = -[gL_E*mVE;
        gL_I*mVI]*1000; % leak has opposite signs of mean Vs
    
f_EnI = -(MatEI-eye(2))^-1*(Leak+VecInput);

end
%% single cell simulation: to collect mean V (and firing rates)


function [meanV,fr] = MEanFieldEst_SingleCell_L6(NeuronType, f_EnI, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         S_EL6,S_IL6,rE_L6,rI_L6,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi,LIFSimuT)
%% attribute parameters for different types of neurons
if strcmpi(NeuronType,'e')
    N_E = N_EE; N_I = N_EI;
    S_E = S_EE; S_I = S_EI;
    p_Fail = p_EEFail;
    lambda = lambda_E; S_lgn = S_Elgn; r_amb = rE_amb;
    S_L6 = S_EL6; r_L6 = rE_L6;
    gL = gL_E;
    rho_ampa = rhoE_ampa; rho_nmda = rhoE_nmda;
elseif strcmpi(NeuronType,'i')
    N_E = N_IE; N_I = N_II;
    S_E = S_IE; S_I = S_II;
    p_Fail = 0;
    lambda = lambda_I; S_lgn = S_Ilgn; r_amb = rI_amb;
    S_L6 = S_IL6; r_L6 = rI_L6;
    gL = gL_I;
    rho_ampa = rhoI_ampa; rho_nmda = rhoI_nmda;
else 
    disp('***Unrecognized Neuron Type')    
end
rE = f_EnI(1)/1000; rI = f_EnI(2)/1000; % f_EnI in s^-1, but here we use ms^-1
%% Evolve single neurons
T = LIFSimuT; % in ms. Default: 20*1e3 ms
dt = 0.1; t = 0:dt:T;
SampleProp = 9/10; % last half time for meanV

v = zeros(size(t)); 
G_gaba_D = zeros(size(t)); G_gaba_R = zeros(size(t));
G_ampa_D = zeros(size(t)); G_ampa_R = zeros(size(t));
G_nmda_D = zeros(size(t)); G_nmda_R = zeros(size(t));
spike = [];

%rng(100)
% input determination: Assume all Poisson
%rng(100)
p_lgn = dt*lambda;                    Sp_lgn = double(rand(size(t))<=p_lgn);
%rng(101)
p_amb = dt*r_amb;                     Sp_amb = double(rand(size(t))<=p_amb);
%rng(102)
p_EV1 = dt*rE*full(N_E)*(1-p_Fail);   Sp_EV1 = double(rand(size(t))<=p_EV1);
%rng(103)
p_IV1 = dt*rI*full(N_I);              Sp_IV1 = double(rand(size(t))<=p_IV1);
p_L6  = dt*r_L6;                      Sp_L6  = double(rand(size(t))<=p_L6);

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
         G_I = 1/(tau_gaba_D-tau_gaba_R) * (G_gaba_D(tInd) - G_gaba_R(tInd)); % S_EI is included in amplitude of GE_gaba
         G_E = 1/(tau_ampa_D-tau_ampa_R) * (G_ampa_D(tInd) - G_ampa_R(tInd)) ...
             + 1/(tau_nmda_D-tau_nmda_R) * (G_nmda_D(tInd) - G_nmda_R(tInd)); % 
         vv  = v(tInd) + dt*(-gL*v(tInd) - G_E.*(v(tInd)-Ve) - G_I.*(v(tInd)-Vi)); 
         if vv >= 1
             v(tInd+1) = nan;
             spike = [spike,t(tInd)];
         else
             v(tInd+1) = vv;
         end
     end

     % conductances
     G_gaba_R(tInd+1) = (G_gaba_R(tInd) +                                           S_I*Sp_IV1(tInd)                                     ) * exp(-dt/tau_gaba_R); 
     G_gaba_D(tInd+1) = (G_gaba_D(tInd) +                                           S_I*Sp_IV1(tInd)                                     ) * exp(-dt/tau_gaba_D); 
     G_ampa_R(tInd+1) = (G_ampa_R(tInd) + S_lgn*Sp_lgn(tInd) + S_amb*Sp_amb(tInd) + S_E*Sp_EV1(tInd)*rho_ampa + S_L6*Sp_L6(tInd)*rho_ampa) * exp(-dt/tau_ampa_R);
     G_ampa_D(tInd+1) = (G_ampa_D(tInd) + S_lgn*Sp_lgn(tInd) + S_amb*Sp_amb(tInd) + S_E*Sp_EV1(tInd)*rho_ampa + S_L6*Sp_L6(tInd)*rho_ampa) * exp(-dt/tau_ampa_D);
     G_nmda_R(tInd+1) = (G_nmda_R(tInd) +                                           S_E*Sp_EV1(tInd)*rho_nmda + S_L6*Sp_L6(tInd)*rho_nmda) * exp(-dt/tau_nmda_R);
     G_nmda_D(tInd+1) = (G_nmda_D(tInd) +                                           S_E*Sp_EV1(tInd)*rho_nmda + S_L6*Sp_L6(tInd)*rho_nmda) * exp(-dt/tau_nmda_D);
end
meanV = mean(v(floor(end*(1-SampleProp)):end), 'omitnan');
fr = length(find(spike>(1-SampleProp)*T))/(T*SampleProp/1000);
end
