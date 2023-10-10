% single cell simulation: to collect mean V (and firing rates)

function [mVLIF,FrLIF] = LIF1Pixel(Fr_MFinv, N_PreSynPix, L4SE,L4SI, L4CE,L4CI, L4IE,L4II,... %8
                                   S_EE,S_EI,S_IE,S_II,p_EEFail,...                             %5
                                   S_EL6,S_IL6,rL6E,rL6I,S_amb,rE_amb,rI_amb,...                %7 L6 Amb                                   
                                   lgn_SOnOff,lgn_COnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,...%8
                                   tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...%7
                                   rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...                  % 4
                                   gL_E,gL_I,Ve,Vi,LIFSimuT, dt, varargin)                      % 6+x  %No more Freq
                                 % Last two lines in case we have used current input...
                               %% First check if L4 rates match neuron numbers
if length(Fr_MFinv) == 5 % Son Con Soff Coff, I
    Fr_MFinv(Fr_MFinv<0) = 0;
    rS = (Fr_MFinv(1) + Fr_MFinv(3))/2; 
    rC = (Fr_MFinv(2) + Fr_MFinv(4))/2; 
    rI = Fr_MFinv(5); % f_EnI in s^-1, but here we use ms^-1
else
    error('L4 Pix FRs dont match!')
end

if nargin > 45
    RecdThre = varargin{1};
else
    RecdThre = eps;
end
%% PreSynaptic neuron numbers. External and internal (pixel)
N_SSi = N_PreSynPix(1,1); N_CSi = N_PreSynPix(1,2); N_ISi = N_PreSynPix(1,3); %% REally?? Row/Col
N_SCi = N_PreSynPix(2,1); N_CCi = N_PreSynPix(2,2); N_ICi = N_PreSynPix(2,3); 
N_SIi = N_PreSynPix(3,1); N_CIi = N_PreSynPix(3,2); N_IIi = N_PreSynPix(3,3); 
%% Setup input and record
T = LIFSimuT; % in ms. Default: 20*1e3 ms
tt = 0:dt:T;
SampleProp = 19/20; % last half time for meanV
DroptInd = floor(length(tt)*(1-SampleProp));
TimeFrac = 1; dtMltp = 1;
% L4 input determination: Assume all Poisson
% Can be optimized in the future, but let's make var correct for now
L4SEAll = (L4SE + rS*N_SSi + rC*N_SCi)*(1-p_EEFail); L4SIAll = L4SI + rI*N_SIi; %fail multiply here. NOWHERE ELSE!
L4CEAll = (L4CE + rS*N_CSi + rC*N_CCi)*(1-p_EEFail); L4CIAll = L4CI + rI*N_CIi;
L4IEAll =  L4IE + rS*N_ISi + rC*N_ICi;               L4IIAll = L4II + rI*N_IIi;
L4EInputNow = [L4SEAll; L4CEAll; L4SEAll; L4CEAll; L4IEAll]/1000; 
L4IInputNow = [L4SIAll; L4CIAll; L4SIAll; L4CIAll; L4IIAll]/1000; 
L4Pix_EventsEU = PoissonInputForNetwork(5,L4EInputNow,LIFSimuT*TimeFrac,dt*dtMltp); % E to everyone
L4Pix_EventsIU = PoissonInputForNetwork(5,L4IInputNow,LIFSimuT*TimeFrac,dt*dtMltp); % I to everyone 
% LGN, amb, L6
LgnInputNow = [lgn_SOnOff(1)*NlgnS; 
               lgn_COnOff(1)*NlgnC;
               lgn_SOnOff(2)*NlgnS; 
               lgn_COnOff(2)*NlgnC;
               lgn_I        *NlgnI];
lgnPix_Events = PoissonInputForNetwork(5,LgnInputNow,LIFSimuT*TimeFrac,dt*dtMltp);
AmbPix_Events = PoissonInputForNetwork(5,[rE_amb*ones(4,1);rI_amb],LIFSimuT*TimeFrac,dt*dtMltp);
L6Pix_Events  = PoissonInputForNetwork(5,[rL6E*ones(4,1);    rL6I],LIFSimuT*TimeFrac,dt*dtMltp);
% Leakage
gL = [gL_E*ones(4,1); gL_I*ones(1)];
% Initialize records % sumV = zeros(5, 1);% VNanCount = zeros(5, 1);
%mVs = zeros(5,floor(T/T1s));
vt = zeros(5, 1);
vRecord = zeros(5, length(tt));
SpLIF = zeros(5, 1);
%% Can precompute all Gs...
winAMPA = 0:dt:tau_ampa_D*5;
winNMDA = 0:dt:tau_nmda_D*5;
winGABA = 0:dt:tau_gaba_D*5;
KerAMPA = 1/(tau_ampa_D-tau_ampa_R) * (exp(-winAMPA/tau_ampa_D) - exp(-winAMPA/tau_ampa_R));
KerNMDA = 1/(tau_nmda_D-tau_nmda_R) * (exp(-winNMDA/tau_nmda_D) - exp(-winNMDA/tau_nmda_R));
KerGABA = 1/(tau_gaba_D-tau_gaba_R) * (exp(-winGABA/tau_gaba_D) - exp(-winGABA/tau_gaba_R));
KerAMPA = KerAMPA / (sum(KerAMPA)*dt);
KerNMDA = KerNMDA / (sum(KerNMDA)*dt);
KerGABA = KerGABA / (sum(KerGABA)*dt);
% Incorporate L4 and other input
ampaInp = single([S_Elgn * lgnPix_Events(1:4,:);            S_Ilgn * lgnPix_Events(5,:)]...
               + [S_amb  * AmbPix_Events ]... % S_amb identical for E and I
               + [S_EL6  * L6Pix_Events(1:4,:) * rhoE_ampa; S_IL6  * L6Pix_Events(5,:) * rhoI_ampa]...
               + [S_EE * L4Pix_EventsEU(1:4,:) * rhoE_ampa; S_IE * L4Pix_EventsEU(5,:) * rhoI_ampa]);
nmdaInp = single([S_EL6  * L6Pix_Events(1:4,:) * rhoE_nmda; S_IL6  * L6Pix_Events(5,:) * rhoI_nmda]...
               + [S_EE * L4Pix_EventsEU(1:4,:) * rhoE_nmda; S_IE * L4Pix_EventsEU(5,:) * rhoI_nmda]);
gabaInp = single([S_EI * L4Pix_EventsIU(1:4,:);             S_II * L4Pix_EventsIU(5,:)] );
%% Precompute All Conductances: Speed up by ifft/fft
n = size(ampaInp,2);
GAMPA = ifft(fft(ampaInp',n) .* repmat(fft(KerAMPA',n),1,5))';
GNMDA = ifft(fft(nmdaInp',n) .* repmat(fft(KerNMDA',n),1,5))';
GGABA = ifft(fft(gabaInp',n) .* repmat(fft(KerGABA',n),1,5))';

% GAMPA = conv2(ampaInp,KerAMPA);
% GNMDA = conv2(nmdaInp,KerAMPA);
% GGABA = conv2(gabaInp,KerAMPA);
GE = GAMPA + GNMDA;  GI = GGABA;
%% Evolve single neurons
RefTimer = zeros(5, 1);
RecordInd = 0;
InpWin = 500;
mVs = zeros(5,floor(T/InpWin));
for tInd = 1:length(tt)-1 
    % First, Get input matrices from series
     FrameNum = floor(InpWin/dt);
    FrameInd = mod(tInd, FrameNum);
    if FrameInd == 0
    FrameInd = FrameNum;
    end
    
    % When sample from Gs, GAMPA and GNMDA should use same time frame;
    if FrameInd == 1 % if the first frame, recompute input mats
        GEU = GE(:,RecordInd*FrameNum+1:(RecordInd+1)*FrameNum);  
        GIU = GI(:,RecordInd*FrameNum+1:(RecordInd+1)*FrameNum); 
        vRecord = zeros(5,FrameNum);
        RecordInd = RecordInd +1;
    end
    % Firstly, refrectory neurons get out if timer reachs t_ref
    RefTimer(isnan(vt)) = RefTimer(isnan(vt)) + dt;
    vt(RefTimer>=tau_ref) = 0;
    RefTimer(RefTimer>=tau_ref) = 0;
    
    ov  = vt + dt*(-gL.*vt - GEU(:,FrameInd).*(vt-Ve) - GIU(:,FrameInd).*(vt-Vi));
    SpNow = single(ov>=1);
    ov(ov>=1) = nan;
    vRecord(:,FrameInd) = ov;
    if FrameInd == FrameNum
        vRecord(:,end) = [];
        vRecord(vRecord < RecdThre & vRecord >0) = nan; %% NOTE: Excluding too low voltages
        mVs(:,RecordInd) = nanmean(vRecord,2);
    end
    % Compute meanV
     if tInd > DroptInd
         SpLIF = SpLIF + SpNow;
     end
    % Transfer variables
    vt = ov; 
end
mVUseInd = (1:size(mVs,2))/size(mVs,2) > 1-SampleProp;
mVLIF = nanmean(mVs(:,mVUseInd),2);
FrLIF = SpLIF/(T*SampleProp/1000);

figure
subplot 221; plot(tt, GAMPA(3,1:length(tt))); xlim([T*9/10,T]); xlabel('T(ms)'); ylabel('GAMPA')
subplot 222; plot(tt, GNMDA(3,1:length(tt))); xlim([T*9/10,T]); xlabel('T(ms)'); ylabel('GNMDA')
subplot 223; plot(tt, GGABA(3,1:length(tt))); xlim([T*9/10,T]); xlabel('T(ms)'); ylabel('GGABA')
%subplot 224; plot(linspace(0,), vRecord(3,:)); xlim([T*9/10,T]); xlabel('T(ms)'); ylabel('V')
end