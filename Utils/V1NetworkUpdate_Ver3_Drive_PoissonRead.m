%% NOTE: Add ambient, add synaptic failure, maybe change ref setup
%% Great Change here:
%  Adding E-to-E synaptic delay for 1 ms
%  EEDlyRcd is an N*T sparse matrix, each column record the number of
%  pending kicks that is waiting, with T_i*dt ms remaining. Each raw record
%  all the pending E kicks for one E neuron

% New Adding 04/02/2021
% Read Poisson from existing event series, rather than generate every time.
% Update 04/14/2021: Actually, let's precompute input series, and feed it
% in. No more reading spike data.
% Update 04/21/2021: Now E-to-I also have delay
function [oRefTimeE,oVE,oSpE,oGE_ampa_R,oGE_nmda_R,oGE_gaba_R,... % Output
                             oGE_ampa_D,oGE_nmda_D,oGE_gaba_D,...
          oRefTimeI,oVI,oSpI,oGI_ampa_R,oGI_nmda_R,oGI_gaba_R,...
                             oGI_ampa_D,oGI_nmda_D,oGI_gaba_D,...
                             oEEDlyRcd,oIEDlyRcd] = ... % A updated N*T Mat recoding the time of kicks taking effect
          V1NetworkUpdate_Ver3_Drive_PoissonRead(RefTimeE,VE,SpE,GE_ampa_R,GE_nmda_R,GE_gaba_R,... % These are input kept updating
                                          GE_ampa_D,GE_nmda_D,GE_gaba_D,...
                          RefTimeI,VI,SpI,GI_ampa_R,GI_nmda_R,GI_gaba_R,...
                                          GI_ampa_D,GI_nmda_D,GI_gaba_D,...
                                          EEDlyRcd,IEDlyRcd,... % A N*T Mat recoding the time of kicks taking effect
                          C_EE,C_EI,C_IE,C_II,...
                          S_EE,S_EI,S_IE,S_II,...
                          tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,... % time unit is ms
                          dt,p_EEFail,...
                          gL_E,Ve,rhoE_ampa,rhoE_nmda,...
                          gL_I,Vi,rhoI_ampa,rhoI_nmda,...
                          EampaInp, IampaInp, EnmdaInp, InmdaInp) % L6 Parameters          % The lower/upper bounds for kick waiting time
%% Firstly, refrectory neurons get out due to timer. NaN stand for neuron in ref in Vs
RefTimeE(isnan(VE)) = RefTimeE(isnan(VE)) + dt; % timer for all ref neurons plus dt
VE(RefTimeE>=tau_ref) = 0; % For ref time up, kick ref neurons out
oRefTimeE = RefTimeE; oRefTimeE(RefTimeE>=tau_ref) = 0; % Those times for kicked-out neurons also reset
% same for I neurons
RefTimeI(isnan(VI)) = RefTimeI(isnan(VI)) + dt;
VI(RefTimeI>=tau_ref) = 0;
oRefTimeI = RefTimeI; oRefTimeI(RefTimeI>=tau_ref) = 0;
%% output of V and spikes, based on existing Gs                         + S_EL6 * SpE_L6 * rhoE_ampa        
GE_I = 1/(tau_gaba_D-tau_gaba_R) * (GE_gaba_D - GE_gaba_R); % S_EI is included in amplitude of GE_gaba
GE_E = 1/(tau_ampa_D-tau_ampa_R) * (GE_ampa_D - GE_ampa_R) + 1/(tau_nmda_D-tau_nmda_R) * (GE_nmda_D - GE_nmda_R); %                       
oVE = VE + dt*(-gL_E*VE - GE_E.*(VE-Ve) - GE_I.*(VE-Vi));
oSpE = sparse(double(oVE>=1)); 
oVE(oVE>=1) = nan; % nan represent refractory

GI_I = 1/(tau_gaba_D-tau_gaba_R) * (GI_gaba_D - GI_gaba_R); % S_EI is included in amplitude of GE_gaba
GI_E = 1/(tau_ampa_D-tau_ampa_R) * (GI_ampa_D - GI_ampa_R) + 1/(tau_nmda_D-tau_nmda_R) * (GI_nmda_D - GI_nmda_R); %                       
oVI = VI + dt*(-gL_I*VI - GI_E.*(VI-Ve) - GI_I.*(VI-Vi));
oSpI = sparse(double(oVI>=1)); 
oVI(oVI>=1) = nan; % nan represent refractory   

%% 'External' stimuli coin toss
% External stimuli. Can be
%%% NOTE!!: I use direct implementation of event generation now: generate
%%% poission random numbers
%p_Elgn = 1-exp(-dt*lambda_E); p_Ilgn = 1-exp(-dt*lambda_I); % Caution! lambda with unit ms^-1
% p_Elgn = dt*lambda_E; p_Ilgn = dt*lambda_I;
% %SpE_lgn = sparse(poissrnd(p_Elgn,size(VE))); SpI_lgn = sparse(poissrnd(p_Ilgn,size(VI)));
% SpE_lgn = sparse(double(rand(size(VE))<=p_Elgn)); SpI_lgn = sparse(double(rand(size(VI))<=p_Ilgn));
% % and Ambient
% %p_EAmb = 1-exp(-dt*rE_amb); p_IAmb = 1-exp(-dt*rI_amb); % Caution! r with unit ms^-1
% p_EAmb = dt*rE_amb; p_IAmb = dt*rI_amb;
% % SpE_amb = sparse(poissrnd(p_EAmb,size(VE))); SpI_amb = sparse(poissrnd(p_IAmb,size(VI)));
% SpE_amb = sparse(double(rand(size(VE))<=p_EAmb)); SpI_amb = sparse(double(rand(size(VI))<=p_IAmb));
% % L6
% p_EL6  = dt*rE_L6; p_IL6 = dt*rI_L6;
% % SpE_L6 = sparse(poissrnd(p_EL6,size(VE))); SpI_L6 = sparse(poissrnd(p_IL6,size(VI)));
% SpE_L6 = sparse(double(rand(size(VE))<=p_EL6)); SpI_L6 = sparse(double(rand(size(VI))<=p_IL6));

%% Instead, all input are already given now.
%% E-to-E has 1.Delay and 2.Synaptic Failure
RawEPSS_E = EEDlyRcd(:,1);%EPSS = e to e post synaptic spike effect
[nE1_E,nE2_E,Raw_NonZero_E]  = find(RawEPSS_E);
SynFail_E = binornd(Raw_NonZero_E,1-p_EEFail);
EPSS_E = sparse(nE1_E,nE2_E,SynFail_E,size(RawEPSS_E,1),size(RawEPSS_E,2));
% E-to-I has 1.Delay 
EPSS_I = IEDlyRcd(:,1);%EPSS = e to e post synaptic spike effect
% first respond to spikes
% gaba
GE_gaba_R_inter = GE_gaba_R + S_EI * (C_EI * SpI);  
GE_gaba_D_inter = GE_gaba_D + S_EI * (C_EI * SpI);
GI_gaba_R_inter = GI_gaba_R + S_II * (C_II * SpI);  
GI_gaba_D_inter = GI_gaba_D + S_II * (C_II * SpI);
% ampa
GE_ampa_R_inter = GE_ampa_R + full(S_EE * EPSS_E * rhoE_ampa) + EampaInp; % Only Recurrent and Ext input now!
GE_ampa_D_inter = GE_ampa_D + full(S_EE * EPSS_E * rhoE_ampa) + EampaInp;
GI_ampa_R_inter = GI_ampa_R + full(S_IE * EPSS_I * rhoI_ampa) + IampaInp;
GI_ampa_D_inter = GI_ampa_D + full(S_IE * EPSS_I * rhoI_ampa) + IampaInp;
% nmda
GE_nmda_R_inter = GE_nmda_R + full(S_EE * EPSS_E * rhoE_nmda) + EnmdaInp;
GE_nmda_D_inter = GE_nmda_D + full(S_EE * EPSS_E * rhoE_nmda) + EnmdaInp;
GI_nmda_R_inter = GI_nmda_R + full(S_IE * EPSS_I * rhoI_nmda) + InmdaInp;
GI_nmda_D_inter = GI_nmda_D + full(S_IE * EPSS_I * rhoI_nmda) + InmdaInp;

% then exponetial decay
% gaba
oGE_gaba_R = GE_gaba_R_inter * exp(-dt/tau_gaba_R);  
oGE_gaba_D = GE_gaba_D_inter * exp(-dt/tau_gaba_D);
oGI_gaba_R = GI_gaba_R_inter * exp(-dt/tau_gaba_R);  
oGI_gaba_D = GI_gaba_D_inter * exp(-dt/tau_gaba_D);
% ampa
oGE_ampa_R = GE_ampa_R_inter * exp(-dt/tau_ampa_R); 
oGE_ampa_D = GE_ampa_D_inter * exp(-dt/tau_ampa_D); 
oGI_ampa_R = GI_ampa_R_inter * exp(-dt/tau_ampa_R); 
oGI_ampa_D = GI_ampa_D_inter * exp(-dt/tau_ampa_D); 
% nmda
oGE_nmda_R = GE_nmda_R_inter * exp(-dt/tau_nmda_R); 
oGE_nmda_D = GE_nmda_D_inter * exp(-dt/tau_nmda_D); 
oGI_nmda_R = GI_nmda_R_inter * exp(-dt/tau_nmda_R); 
oGI_nmda_D = GI_nmda_D_inter * exp(-dt/tau_nmda_D); 

%% Update pending kick recording:
oEEDlyRcd = sparse(size(EEDlyRcd,1),size(EEDlyRcd,2));
oEEDlyRcd(:,1:end-1) = EEDlyRcd(:,2:end); % first, move all kicks one-dt forward
oEEDlyRcd(:,end) = C_EE * oSpE; % put the new generated kicks to the last column
oIEDlyRcd = sparse(size(IEDlyRcd,1),size(IEDlyRcd,2));
oIEDlyRcd(:,1:end-1) = IEDlyRcd(:,2:end); % first, move all kicks one-dt forward
oIEDlyRcd(:,end) = C_IE * oSpE;

end