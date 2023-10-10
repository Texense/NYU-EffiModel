%% single cell simulation: to collect mean V (and firing rates)


function [v, meanV,fr] = MEanFieldEst_SingleCell_L6(NeuronType, f_EnI, ...
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
    return
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
p_lgn = dt*lambda;                    Sp_lgn = poissrnd(p_lgn,size(t,1),size(t,2));
%rng(101)
p_amb = dt*r_amb;                     Sp_amb = poissrnd(p_amb,size(t,1),size(t,2));
%rng(102)
p_EV1 = dt*rE*full(N_E)*(1-p_Fail);   Sp_EV1 = poissrnd(p_EV1,size(t,1),size(t,2));
%rng(103)
p_IV1 = dt*rI*full(N_I);              Sp_IV1 = poissrnd(p_IV1,size(t,1),size(t,2));
p_L6  = dt*r_L6;                      Sp_L6  = poissrnd(p_L6, size(t,1),size(t,2));

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
