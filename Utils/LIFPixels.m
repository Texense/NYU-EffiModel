%% For each pixel, use 1E and 1I neurons to represent them
function [mVELIF,mVILIF] = LIFPixels(Fr_MFinv, L4Pix_EventsEIni, L4Pix_EventsIIni, L4EInputIni, L4IInputIni,...
                                     C_EE_Pixel_Us,C_EI_Pixel_Us,C_IE_Pixel_Us,C_II_Pixel_Us,...
                                     S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                     lgnEPix_Events,S_Elgn,AmbEPix_Events,S_amb,...
                                     lgnIPix_Events,S_Ilgn,AmbIPix_Events,...
                                     S_EL6,S_IL6,L6EPix_Events,L6IPix_Events,...
                                     tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                     rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                     gL_E,gL_I,Ve,Vi,LIFSimuT, dt, PixNum,...
                                     LGNCurInp, L6CurInp,...
                                     lambda_E_Pixel,lambda_I_Pixel,L6E_Pixel,L6I_Pixel)
                               % Last two lines in case we have used current input...
%% First check if L4 rates match neuron numbers
if length(Fr_MFinv) == 2*PixNum
rE = Fr_MFinv(1:PixNum); rI = Fr_MFinv(PixNum+1:2*PixNum); % f_EnI in s^-1, but here we use ms^-1
else
    error('L4 FRs dont match!')
end
%% Setup input and record
T = LIFSimuT; % in ms. Default: 20*1e3 ms
tt = 0:dt:T;
SampleProp = 9/10; % last half time for meanV
DroptInd = floor(length(tt)*(1-SampleProp));
% L4 input determination: Assume all Poisson
L4EInputNow = [C_EE_Pixel_Us*rE*(1-p_EEFail); 
               C_IE_Pixel_Us*rE             ]/1000;   
L4IInputNow = [C_EI_Pixel_Us*rI             ; 
               C_II_Pixel_Us*rI             ]/1000;   
L4EInputNow(L4EInputNow<0) = 0;
L4IInputNow(L4IInputNow<0) = 0;
L4Pix_EventsEU = L4Pix_EventsEIni.*repmat(L4EInputNow./L4EInputIni,1,size(L4Pix_EventsEIni,2));
L4Pix_EventsIU = L4Pix_EventsIIni.*repmat(L4IInputNow./L4IInputIni,1,size(L4Pix_EventsIIni,2));
% Incorporate L4 and other input
EampaInp = single(full(S_Elgn * (lgnEPix_Events + lambda_E_Pixel*dt*LGNCurInp) ...
                     + S_amb  *  AmbEPix_Events ...
                     + S_EL6  * (L6EPix_Events  + L6E_Pixel*dt*L6CurInp) * rhoE_ampa)); %  * rhoE_ampa
IampaInp = single(full(S_Ilgn * (lgnIPix_Events+ lambda_I_Pixel*dt*LGNCurInp) ...
                     + S_amb  *  AmbIPix_Events ...
                     + S_IL6  * (L6IPix_Events  + L6I_Pixel*dt*L6CurInp) * rhoI_ampa)); %  * rhoI_ampa

EnmdaInp = single(full(S_EL6  * (L6EPix_Events  + L6E_Pixel*dt*L6CurInp) * rhoE_nmda)); %
InmdaInp = single(full(S_IL6  * (L6IPix_Events  + L6I_Pixel*dt*L6CurInp) * rhoI_nmda));
L4InputfEampa = single(full( [S_EE * L4Pix_EventsEU(1:PixNum,:) * rhoE_ampa;
                              S_IE * L4Pix_EventsEU(PixNum+1:2*PixNum,:) * rhoI_ampa] ));
L4InputfEnmda = single(full( [S_EE * L4Pix_EventsEU(1:PixNum,:) * rhoE_nmda;
                              S_IE * L4Pix_EventsEU(PixNum+1:2*PixNum,:) * rhoI_nmda] ));
L4InputfI     = single(full( [S_EI * L4Pix_EventsIU(1:PixNum,:);
                              S_II * L4Pix_EventsIU(PixNum+1:2*PixNum,:)] ));

% Leakage
gL = [gL_E*ones(size(rE)); gL_I*ones(size(rI))];
% Initialize records
sumV = zeros(size(Fr_MFinv));
VNanCount = zeros(size(Fr_MFinv));
vt = zeros(size(Fr_MFinv));
%% Precompute All Conductances: Speed up by ifft/fft 
winAMPA = 0:dt:tau_ampa_D*5;
winNMDA = 0:dt:tau_nmda_D*5;
winGABA = 0:dt:tau_gaba_D*5;
KerAMPA = 1/(tau_ampa_D-tau_ampa_R) * (exp(-winAMPA/tau_ampa_D) - exp(-winAMPA/tau_ampa_R));
KerNMDA = 1/(tau_nmda_D-tau_nmda_R) * (exp(-winNMDA/tau_nmda_D) - exp(-winNMDA/tau_nmda_R));
KerGABA = 1/(tau_gaba_D-tau_gaba_R) * (exp(-winGABA/tau_gaba_D) - exp(-winGABA/tau_gaba_R));
KerAMPA = KerAMPA / (sum(KerAMPA)*dt);
KerNMDA = KerNMDA / (sum(KerNMDA)*dt);
KerGABA = KerGABA / (sum(KerGABA)*dt);
n = size(EampaInp,2); 
GAMPA = ifft(fft(([EampaInp; IampaInp] + L4InputfEampa)',n) .* repmat(fft(KerAMPA',n),1,2*PixNum))';
GNMDA = ifft(fft(([EnmdaInp; InmdaInp] + L4InputfEnmda)',n) .* repmat(fft(KerNMDA',n),1,2*PixNum))';
GGABA = ifft(fft(L4InputfI',n)                              .* repmat(fft(KerGABA',n),1,2*PixNum))'; 
GE = GAMPA + GNMDA; GI = GGABA;
%% Evolve single neurons
RefTimer = zeros(size(Fr_MFinv)); 
% G_ampa_R = zeros(size(Fr_MFinv)); 
% G_nmda_R = zeros(size(Fr_MFinv)); 
% G_gaba_R = zeros(size(Fr_MFinv)); 
% G_ampa_D = zeros(size(Fr_MFinv)); 
% G_nmda_D = zeros(size(Fr_MFinv)); 
% G_gaba_D = zeros(size(Fr_MFinv));                                          
%% Can precompute all Gs...
for tInd = 1:length(tt)-1
    % First, Get input matrices from series
    InpWin = 10; FrameNum = floor(InpWin/dt);
    FrameInd = mod(tInd, FrameNum);
    if FrameInd == 0
    FrameInd = FrameNum;
    end
    
    % When sample from Gs, GAMPA and GNMDA should use same time frame;
    if FrameInd == 1 % if the first frame, recompute input mats
        RandTimBin = randi([1 n-FrameNum],2,1);
        GEU = GE(:,RandTimBin(1):RandTimBin(1)+FrameNum);  
        GIU = GI(:,RandTimBin(2):RandTimBin(2)+FrameNum); 
    end
    
    % Firstly, refrectory neurons get out if timer reachs t_ref
    RefTimer(isnan(vt)) = RefTimer(isnan(vt)) + dt;
    vt(RefTimer>=tau_ref) = 0;
    RefTimer(RefTimer>=tau_ref) = 0;
    
%     G_I = 1/(tau_gaba_D-tau_gaba_R) * (G_gaba_D - G_gaba_R); % S_EI is included in amplitude of GE_gaba
%     G_E = 1/(tau_ampa_D-tau_ampa_R) * (G_ampa_D - G_ampa_R) ...
%         + 1/(tau_nmda_D-tau_nmda_R) * (G_nmda_D - G_nmda_R); %
    ov  = vt + dt*(-gL.*vt - GEU(:,FrameInd).*(vt-Ve) - GIU(:,FrameInd).*(vt-Vi));
    ov(ov>=1) = nan;

    % conductances
%     G_gaba_R = (G_gaba_R +                                               L4InputfI(:,FrameInd)) * exp(-dt/tau_gaba_R);
%     G_gaba_D = (G_gaba_D +                                               L4InputfI(:,FrameInd)) * exp(-dt/tau_gaba_D);
%     G_ampa_R = (G_ampa_R + [EampaInp(:,FrameInd);IampaInp(:,FrameInd)] + L4InputfEampa(:,FrameInd)) * exp(-dt/tau_ampa_R);
%     G_ampa_D = (G_ampa_D + [EampaInp(:,FrameInd);IampaInp(:,FrameInd)] + L4InputfEampa(:,FrameInd)) * exp(-dt/tau_ampa_D);
%     G_nmda_R = (G_nmda_R + [EnmdaInp(:,FrameInd);InmdaInp(:,FrameInd)] + L4InputfEnmda(:,FrameInd)) * exp(-dt/tau_nmda_R);
%     G_nmda_D = (G_nmda_D + [EnmdaInp(:,FrameInd);InmdaInp(:,FrameInd)] + L4InputfEnmda(:,FrameInd)) * exp(-dt/tau_nmda_D);
%     
    % Compute meanV
    if tInd > DroptInd
      % sumV = sum([sumV,ov],2,'omitnan');
       VNanCount = VNanCount + isnan(ov);
       osv = ov; osv(isnan(osv)) = 0;
       sumV = sumV + osv;
    end
    % Transfer variables
    vt = ov; 
end
mVLIF = sumV./(length(tt)-1 - DroptInd - VNanCount);
mVELIF = mVLIF(1:PixNum);
mVILIF = mVLIF(PixNum+1:2*PixNum);
% fr = length(find(spike>(1-SampleProp)*T))/(T*SampleProp/1000);
end