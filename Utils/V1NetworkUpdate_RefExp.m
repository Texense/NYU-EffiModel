%% NOTE: Add ambient, add synaptic failure, maybe change ref setup
function [oVE,oSpE,oGE_ampa_R,oGE_nmda_R,oGE_gaba_R,...
                   oGE_ampa_D,oGE_nmda_D,oGE_gaba_D,...
          oVI,oSpI,oGI_ampa_R,oGI_nmda_R,oGI_gaba_R,...
                   oGI_ampa_D,oGI_nmda_D,oGI_gaba_D] = ...
          V1NetworkUpdate_RefExp(VE,SpE,GE_ampa_R,GE_nmda_R,GE_gaba_R,...
                                 GE_ampa_D,GE_nmda_D,GE_gaba_D,...
                          VI,SpI,GI_ampa_R,GI_nmda_R,GI_gaba_R,...
                                 GI_ampa_D,GI_nmda_D,GI_gaba_D,...
                          C_EE,C_EI,C_IE,C_II,...
                          S_EE,S_EI,S_IE,S_II,...
                          tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,... % time unit is ms
                          dt,p_EEFail,...
                          gL_E,Ve,S_Elgn,rhoE_ampa,rhoE_nmda,...
                          gL_I,Vi,S_Ilgn,rhoI_ampa,rhoI_nmda,...
                          S_amb,lambda_E,lambda_I,rE_amb,rI_amb)
%% Firstly, refrectory neurons get out due to exponetial distributed time
p_refout = 1-exp(-dt/tau_ref);  % the probability of getting out in such time interval
[NnE_ref,~] = find(isnan(VE));  % find ref neurons, which are indicated by nan
NnE_ref(rand(size(NnE_ref))>p_refout) = []; % only preserve neurons winning in coin-flip
VE(NnE_ref) = 0; % Kick the winning neurons out of ref
% same for I neurons
[NnI_ref,~] = find(isnan(VI));
NnI_ref(rand(size(NnI_ref))>p_refout) = []; 
VI(NnI_ref) = 0;
%% output of V and spikes, based on existing Gs                         
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

%% Conductances
% External stimuli . Can be
p_Elgn = 1-exp(-dt*lambda_E); p_Ilgn = 1-exp(-dt*lambda_I); % Caution! lambda with unit ms^-1
SpE_lgn = sparse(double(rand(size(VE))<=p_Elgn));
SpI_lgn = sparse(double(rand(size(VI))<=p_Ilgn));
% and Ambient
p_EAmb = 1-exp(-dt*rE_amb); p_IAmb = 1-exp(-dt*rI_amb); % Caution! r with unit ms^-1
SpE_amb = sparse(double(rand(size(VE))<=p_EAmb));
SpI_amb = sparse(double(rand(size(VI))<=p_IAmb));
% Synaptic Failure
RawEPSS = C_EE * SpE;%EPSS = e to e post synaptic spike effect
[nE1,nE2,Raw_NonZero]  = find(RawEPSS);
SynFail_ampa = binornd(Raw_NonZero,1-p_EEFail);
SynFail_nmda = binornd(Raw_NonZero,1-p_EEFail);
EPSS_ampa = sparse(nE1,nE2,SynFail_ampa,size(RawEPSS,1),size(RawEPSS,2));
EPSS_nmda = sparse(nE1,nE2,SynFail_nmda,size(RawEPSS,1),size(RawEPSS,2));
% first respond to spikes
% gaba
GE_gaba_R_inter = GE_gaba_R + S_EI * (C_EI * SpI);  
GE_gaba_D_inter = GE_gaba_D + S_EI * (C_EI * SpI);
GI_gaba_R_inter = GI_gaba_R + S_II * (C_II * SpI);  
GI_gaba_D_inter = GI_gaba_D + S_II * (C_II * SpI);
% ampa
GE_ampa_R_inter = GE_ampa_R + S_Elgn * SpE_lgn + S_amb * SpE_amb + S_EE * EPSS_ampa * rhoE_ampa;
GE_ampa_D_inter = GE_ampa_D + S_Elgn * SpE_lgn + S_amb * SpE_amb + S_EE * EPSS_ampa * rhoE_ampa;
GI_ampa_R_inter = GI_ampa_R + S_Ilgn * SpI_lgn + S_amb * SpI_amb + S_IE * (C_IE * SpE) * rhoI_ampa;
GI_ampa_D_inter = GI_ampa_D + S_Ilgn * SpI_lgn + S_amb * SpI_amb + S_IE * (C_IE * SpE) * rhoI_ampa;
% nmda
GE_nmda_R_inter = GE_nmda_R + S_EE * EPSS_nmda * rhoE_nmda;
GE_nmda_D_inter = GE_nmda_D + S_EE * EPSS_nmda * rhoE_nmda;
GI_nmda_R_inter = GI_nmda_R + S_IE * (C_IE * SpE) * rhoI_nmda;
GI_nmda_D_inter = GI_nmda_D + S_IE * (C_IE * SpE) * rhoI_nmda;

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
end