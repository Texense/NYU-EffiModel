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
% Ver Test3.1: %%% Maybe not now: We try LSY's idea about using previous mV to stablize iteration
% Ver Test3.2: Newton iteration to find steady point.
% 
function [f_EnIOut,meanVs,loop,SteadyIndicate] = MeanFieldEst_BkGd_Newton_StepSize_TestmV(C_EE,C_EI,C_IE,C_II,... %4
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                                   lambda_E,S_Elgn,rE_amb,S_amb,... %4
                                   lambda_I,S_Ilgn,rI_amb,... %3
                                   tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,... %7
                                   rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                                   gL_E,gL_I,Ve,Vi,... %4
                                   N_HC,n_E_HC,n_I_HC,varargin) %3+x % if varagin non empty, we only export the last state
                                   %This line: we are taking spatial center 
%% vara input process
if nargin > 35  % Specify: The number of loops I want, after stopping criteria met    
    AveLoop = varargin{2};
else 
    % *BOLD TEXT* NewtonL: Maybe not that large
    AveLoop = 100; 
end

if nargin > 36 % Specify: Maximum nubmer of loops before stopping
        StopLoop = varargin{3};
else 
        StopLoop = AveLoop;
end

if nargin > 37 % Specify: stepsize h
    h_Step = varargin{4};
else
    h_Step = 1;        
end

if nargin > 38 % Specify: LIF simulation time
    LIFSimuT = varargin{5};
else
    LIFSimuT = 20*1e3;    % unit in ms    
end

%% parameters and inital conditions
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
%mVE = 0.57; mVI = 0.67; 
%mVE =  0.539; mVI = 0.683; 
mVE = 0.57; mVI = 0.67; 
meanVs = [mVE;mVI];
[f_EnIIni,~] = MeanFieldEst_BkGd_Jacobi(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI);
f_EnIOut = f_EnIIni;

loop = 0;
SteadyCounter = 0; % Indicate the number of loops after steady condition
SteadyIndicate = false;
TestPoints = floor(15); % How many consecutive points we test

%% Now! All this section is about estimating derivatives! Hence, we need to fix seeds. 
%while( norm([mVEpre;mVIpre] - [mVE;mVI]) > 0.01 || norm(f_EnIpre - f_EnI0)>0.1) %%% relative difference for firing rates!!
while SteadyCounter<AveLoop %the formal ending condition                             
[mVE,ParmVE_f_EnI] = MEanFieldEst_SingleCell_Jacobi('e', f_EnIIni, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...NewtonL: Maybe not that large
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi,LIFSimuT);
[mVI,ParmVI_f_EnI] = MEanFieldEst_SingleCell_Jacobi('i', f_EnIIni, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi,LIFSimuT);

% The Jacobi of mV(fE,fI) is a 2*2 matrix                                     
JmV = [ParmVE_f_EnI;
       ParmVI_f_EnI];
% if SteadyIndicate, use previous mV
%     mVEIn = mean(meanVs(1,end-10+1:end)) * 0.9 + mVE*0.1;
%     mVIIn = mean(meanVs(2,end-10+1:end)) * 0.9 + mVI*0.1;
% else
%     mVEIn = mVE;
%     mVIIn = mVI;
% end
[f_EnIFunc, Jf_EnI0] = MeanFieldEst_BkGd_Jacobi(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI);

% But I am solving [fs] = F(fs), so it is G(fs) = [fs]-F(fs) = 0
J_F = Jf_EnI0 * JmV;
J_G = eye(2) - J_F;

% New firing rates based on Newton!
f_EnI0 = f_EnIIni - J_G^-1 * (f_EnIIni - f_EnIFunc);

f_EnI0 = max([f_EnI0,[0;0]],[],2); % If unreasonable: go zero
%f_EnI0 = abs(f_EnI0);                                     
%% The new input!
f_EnIIni = f_EnI0*h_Step + f_EnIIni*(1-h_Step);
%f_EnIIni = f_EnI0;
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

% break out if not reaching convergence after 100 iterations. Tbis number
% should be larger than end condition of SteadyCounter
if (~SteadyIndicate && loop >= StopLoop) 
    disp('Firing rates unconverged')
    break
end
end

if nargin > 34 && strcmpi(varargin(1),'Mean') % Specify: I don't need the trajectory but the final ones    
    f_EnIOut = mean(f_EnIOut(:,end-50:end),2);
    meanVs   = mean(meanVs(:,end-50:end),2);
end
end



%% firing rates and the Jacobi of mean field 
% Input: C_EE,C_EI,C_IE,C_II connectivity matrices
%        N_EE,N_EI,N_IE,N_II number of presynaptic neurons
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rE_amb rI_amb       Ambient inputParmV_f_EnI
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
% Output:f_EnI               Firing rates from the linear eqns
%        Jf_EnI              2*2 Jacobi mat: (fE,fI) of (mVE,mVI)

function [f_EnI, Jf_EnI] = MeanFieldEst_BkGd_Jacobi(N_EE,N_EI,N_IE,N_II,...
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                   lambda_E,S_Elgn,rE_amb,S_amb,...
                                   lambda_I,S_Ilgn,rI_amb,...
                                   gL_E,gL_I,Ve,Vi,mVE,mVI)
                               
% Matrix of firing rates
MatEI = [N_EE*S_EE*(Ve-mVE)*(1-p_EEFail), N_EI*S_EI*(Vi-mVE);
         N_IE*S_IE*(Ve-mVI),              N_II*S_II*(Vi-mVI)];

% External Drive and leaky current. Note that simulaton is by ms, so need *1000
VecInput = [(Ve-mVE) * (lambda_E*S_Elgn + rE_amb*S_amb);
            (Ve-mVI) * (lambda_I*S_Ilgn + rI_amb*S_amb)]*1000;
Leak = -[gL_E*mVE;
         gL_I*mVI]*1000; % leak has opposite signs of mean Vs

f_EnI = -(MatEI-eye(2))^-1*(Leak+VecInput);
%% Jacobian based on Formula: f_EnI = -(MatEI-eye(2))^-1*(Leak+VecInput)
% First, Matrix for effective connectivity
Con = -(MatEI-eye(2));
% Second, the current vector
Cur = (Leak+VecInput);
% Third, partial derivates of Con
ParCon_mVE = [N_EE*S_EE*(1-p_EEFail), N_EI*S_EI;
              0,                      0       ];
ParCon_mVI = [0,                      0;
              N_IE*S_IE,              N_II*S_II];
% Last, partials of Cur
ParCur_mVE = [-(lambda_E*S_Elgn + rE_amb*S_amb)-gL_E;
              0                                     ]*1000;
ParCur_mVI = [0;
              -(lambda_I*S_Ilgn + rI_amb*S_amb)-gL_I]*1000;          
%% Finally, get Jacobians of F(mVE,mVI)
ParF_mVE = Con^(-1) * (ParCur_mVE - ParCon_mVE*Con^(-1)*Cur);
ParF_mVI = Con^(-1) * (ParCur_mVI - ParCon_mVI*Con^(-1)*Cur);

Jf_EnI = [ParF_mVE,ParF_mVI];
end

%% single cell simulation: to collect mean V (and firing rates)
% New Task: Estimate partial derivatives on mV
% Output:   meanV          scalar     mVE or mVI
%           ParmV_f_EnI    1*2 vector [ParmV_fE, ParmV_fI]
function [meanV,ParmV_f_EnI] = MEanFieldEst_SingleCell_Jacobi(NeuronType, f_EnI, ...
                                         N_EE,N_EI,N_IE,N_II,...
                                         S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                         lambda_E,S_Elgn,rE_amb,S_amb,...
                                         lambda_I,S_Ilgn,rI_amb,...
                                         tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                         rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                         gL_E,gL_I,Ve,Vi,LIFSimuT)
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
%% firing rates: independent variables we want to do derivatives on
rE  = f_EnI(1)/1000; rI  = f_EnI(2)/1000; % f_EnI in s^-1, but here we use ms^-1
drE = rE/200;        drI = rI/200; % a small deviation for 1/200 of the firing rate
%% Evolve single neurons: Parameters
T = LIFSimuT; % in ms. Default: 20*1e3 ms
dt = 0.1; t = 0:dt:T;
SampleProp = 9/10; % last half time for meanV

v = zeros(size(t)); 
G_gaba_D = zeros(size(t)); G_gaba_R = zeros(size(t));
G_ampa_D = zeros(size(t)); G_ampa_R = zeros(size(t));
G_nmda_D = zeros(size(t)); G_nmda_R = zeros(size(t));
spike = [];

% input determination: Assume all Poisson
p_lgn = 1-exp(-dt*lambda);            Sp_lgn = double(rand(size(t))<=p_lgn);
p_amb = dt*r_amb;                     Sp_amb = double(rand(size(t))<=p_amb);
% rE spiking series, and small purturbations
SeedE = rand(size(t));
p_EV1 = dt*rE*full(N_E)*(1-p_Fail);       Sp_EV1 = double(SeedE<=p_EV1); % preserving expectation correct
p_EVd = dt*(rE+drE)*full(N_E)*(1-p_Fail); Sp_EVd = double(SeedE<=p_EVd);
% rI spiking series, and small purturbations
SeedI = rand(size(t));
p_IV1 = dt*rI*full(N_I);                  Sp_IV1 = double(SeedI<=p_IV1); % preserving expectation correct
p_IVd = dt*(rI+drI)*full(N_I);            Sp_IVd = double(SeedI<=p_IVd);

%% Simulate the original v seires
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
         vv = v(tInd) + dt*(-gL*v(tInd) - G_E.*(v(tInd)-Ve) - G_I.*(v(tInd)-Vi)); 
         if vv >= 1
             v(tInd+1) = nan;
             spike = [spike,t(tInd)];
         else
             v(tInd+1) = vv;
         end
     end

     % conductances
     G_gaba_R(tInd+1) = (G_gaba_R(tInd) + S_I*Sp_IV1(tInd)) * exp(-dt/tau_gaba_R); 
     G_gaba_D(tInd+1) = (G_gaba_D(tInd) + S_I*Sp_IV1(tInd)) * exp(-dt/tau_gaba_D); 
     G_ampa_R(tInd+1) = (G_ampa_R(tInd) + S_lgn * Sp_lgn(tInd) + S_amb * Sp_amb(tInd) + S_E * Sp_EV1(tInd) * rho_ampa) * exp(-dt/tau_ampa_R);
     G_ampa_D(tInd+1) = (G_ampa_D(tInd) + S_lgn * Sp_lgn(tInd) + S_amb * Sp_amb(tInd) + S_E * Sp_EV1(tInd) * rho_ampa) * exp(-dt/tau_ampa_D);
     G_nmda_R(tInd+1) = (G_nmda_R(tInd) + S_E * Sp_EV1(tInd) * rho_nmda) * exp(-dt/tau_nmda_R);
     G_nmda_D(tInd+1) = (G_nmda_D(tInd) + S_E * Sp_EV1(tInd) * rho_nmda) * exp(-dt/tau_nmda_D);
end
meanV = nanmean(v(floor(end*(1-SampleProp)):end));
%fr = length(find(spike>(1-SampleProp)*T))/(T*SampleProp/1000);

%% Simulate v for fE purturbation
v_drE = zeros(size(t)); 
G_gaba_D = zeros(size(t)); G_gaba_R = zeros(size(t));
G_ampa_D = zeros(size(t)); G_ampa_R = zeros(size(t));
G_nmda_D = zeros(size(t)); G_nmda_R = zeros(size(t));
spike = [];

RefTimer = 0;
for tInd = 1:length(t)-1
    % Firstly, refrectory neurons get out due to exponetial distributed time
     if isnan(v_drE(tInd))
        RefTimer = RefTimer+dt; %RefTimer goes up
        if RefTimer>=tau_ref    %if timer reach tau_ref, kick v out of refrectory
            v_drE(tInd+1) = 0;
            RefTimer = 0;
        else
            v_drE(tInd+1) = nan;
        end
     else
         G_I = 1/(tau_gaba_D-tau_gaba_R) * (G_gaba_D(tInd) - G_gaba_R(tInd)); % S_EI is included in amplitude of GE_gaba
         G_E = 1/(tau_ampa_D-tau_ampa_R) * (G_ampa_D(tInd) - G_ampa_R(tInd)) ...
             + 1/(tau_nmda_D-tau_nmda_R) * (G_nmda_D(tInd) - G_nmda_R(tInd)); % 
         vv = v_drE(tInd) + dt*(-gL*v_drE(tInd) - G_E.*(v_drE(tInd)-Ve) - G_I.*(v_drE(tInd)-Vi)); 
         if vv >= 1
             v_drE(tInd+1) = nan;
             spike = [spike,t(tInd)];
         else
             v_drE(tInd+1) = vv;
         end
     end

     % conductances
     G_gaba_R(tInd+1) = (G_gaba_R(tInd) + S_I*Sp_IV1(tInd)) * exp(-dt/tau_gaba_R); 
     G_gaba_D(tInd+1) = (G_gaba_D(tInd) + S_I*Sp_IV1(tInd)) * exp(-dt/tau_gaba_D); 
     G_ampa_R(tInd+1) = (G_ampa_R(tInd) + S_lgn * Sp_lgn(tInd) + S_amb * Sp_amb(tInd) + S_E * Sp_EVd(tInd) * rho_ampa) * exp(-dt/tau_ampa_R);
     G_ampa_D(tInd+1) = (G_ampa_D(tInd) + S_lgn * Sp_lgn(tInd) + S_amb * Sp_amb(tInd) + S_E * Sp_EVd(tInd) * rho_ampa) * exp(-dt/tau_ampa_D);
     G_nmda_R(tInd+1) = (G_nmda_R(tInd) + S_E * Sp_EVd(tInd) * rho_nmda) * exp(-dt/tau_nmda_R);
     G_nmda_D(tInd+1) = (G_nmda_D(tInd) + S_E * Sp_EVd(tInd) * rho_nmda) * exp(-dt/tau_nmda_D);
end
meanV_drE = nanmean(v_drE(floor(end*(1-SampleProp)):end));
%fr = length(find(spike>(1-SampleProp)*T))/(T*SampleProp/1000);

%% Simulate v for fI purturbation
v_drI = zeros(size(t)); 
G_gaba_D = zeros(size(t)); G_gaba_R = zeros(size(t));
G_ampa_D = zeros(size(t)); G_ampa_R = zeros(size(t));
G_nmda_D = zeros(size(t)); G_nmda_R = zeros(size(t));
spike = [];

RefTimer = 0;
for tInd = 1:length(t)-1
    % Firstly, refrectory neurons get out due to exponetial distributed time
     if isnan(v_drI(tInd))
        RefTimer = RefTimer+dt; %RefTimer goes up
        if RefTimer>=tau_ref    %if timer reach tau_ref, kick v out of refrectory
            v_drI(tInd+1) = 0;
            RefTimer = 0;
        else
            v_drI(tInd+1) = nan;
        end
     else
         G_I = 1/(tau_gaba_D-tau_gaba_R) * (G_gaba_D(tInd) - G_gaba_R(tInd)); % S_EI is included in amplitude of GE_gaba
         G_E = 1/(tau_ampa_D-tau_ampa_R) * (G_ampa_D(tInd) - G_ampa_R(tInd)) ...
             + 1/(tau_nmda_D-tau_nmda_R) * (G_nmda_D(tInd) - G_nmda_R(tInd)); % 
         vv = v_drI(tInd) + dt*(-gL*v_drI(tInd) - G_E.*(v_drI(tInd)-Ve) - G_I.*(v_drI(tInd)-Vi)); 
         if vv >= 1
             v_drI(tInd+1) = nan;
             spike = [spike,t(tInd)];
         else
             v_drI(tInd+1) = vv;
         end
     end

     % conductances
     G_gaba_R(tInd+1) = (G_gaba_R(tInd) + S_I*Sp_IVd(tInd)) * exp(-dt/tau_gaba_R); 
     G_gaba_D(tInd+1) = (G_gaba_D(tInd) + S_I*Sp_IVd(tInd)) * exp(-dt/tau_gaba_D); 
     G_ampa_R(tInd+1) = (G_ampa_R(tInd) + S_lgn * Sp_lgn(tInd) + S_amb * Sp_amb(tInd) + S_E * Sp_EV1(tInd) * rho_ampa) * exp(-dt/tau_ampa_R);
     G_ampa_D(tInd+1) = (G_ampa_D(tInd) + S_lgn * Sp_lgn(tInd) + S_amb * Sp_amb(tInd) + S_E * Sp_EV1(tInd) * rho_ampa) * exp(-dt/tau_ampa_D);
     G_nmda_R(tInd+1) = (G_nmda_R(tInd) + S_E * Sp_EV1(tInd) * rho_nmda) * exp(-dt/tau_nmda_R);
     G_nmda_D(tInd+1) = (G_nmda_D(tInd) + S_E * Sp_EV1(tInd) * rho_nmda) * exp(-dt/tau_nmda_D);
end
meanV_drI = nanmean(v_drI(floor(end*(1-SampleProp)):end));

%% mV derivatives on rE and rI
ParmV_f_EnI = [meanV_drE-meanV, meanV_drI-meanV]./[drE, drI];
end