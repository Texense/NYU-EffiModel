%% MF+V for Pixels. E is divided into S and C cells
%% Input: 
%        N_PreSynPix   within Pixel->Pixel connectivity matrix 
%                           3*3 Post(R)/Pre(C): S, C, I
%                                                       S
%                                                       C
%                                                       I
%        L4E, L4I                 Total E/I input kick numbers
%        Assume L4SE/CE/IE are proportional to presynaptic neuron numbers
%        S_EE,S_EI,S_IE,S_II      synaptic strength
%        S_EL6,S_IL6,rE_L6,rI_L6, L6 Input, All Scalar
%        S_amb,rE_amb,rI_amb      amb Input, All Scalar
%        p_EEFail                 E-to-E Synaptic failure prob
%        lgn_EOnOff,lgn_I         LGN input rate. Phase OnOff for E (2*1 vec) 
%        S_Elgn S_Ilgn            LGN strength of drives
%        NlgnS,NlgnC,NlgnI        LGN neuron # to SCI cells
%        gL_E,gL_I                Leaky time constants
%        Ve,Vi                    Reversal potentials
%        HyperPara           Cell. Entries are simulation hyperparameters
%                            1st: 'Mean' or 'Traj', indicate the form of the output
%                            2nd: Sample Number after stopping criteria
%                            3rd: Max Iteration before converged
%                            4th: Stepsize h
%                            5th: LIF simulation time
% Output:f_EnI               Estimation of firing rates, E and I; A sequences
%        meanVs              mean V of E and I
%        loop                Number of loops
%        SteadyIndicate      logical value for convergence
%% Ver 5:       For single pixels
% iteration
% Zhuo-Cheng Xiao 07/27/2021
function [f_EnIOut,loop,SteadyIndicate,FailureIndicate]...
           = MF_free_LIF_SinglePixel(...
...% MF Parameters                     
                     N_PreSynPix, L4SE,L4SI, L4CE,L4CI, L4IE,L4II,... %3 
                     S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                     S_EL6,S_IL6,rL6E,rL6I,S_amb,rE_amb,rI_amb,...%7 L6 Amb                                   
                     lgn_EOnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... %7
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
dt = 0.1;
% initialize with a reasonable guess     
f_EnIOut = zeros(5,1);

loop = 0;
SteadyCounter = 0; % Indicate the number of loops after steady condition
SteadyIndicate = false;
TestPoints = floor(7); % How many consecutive points we test
%while( norm([mVEpre;mVIpre] - [mVE;mVI]) > 0.01 || norm(f_EnIpre - f_EnI0)>0.1) %%% relative difference for firing rates!!

FailureIndicate = 0;
Suspicious = false;
while SteadyCounter<AveLoop %the formal ending condition
% simulate one neuron with input                             
[~,FrLIF] = LIF1Pixel(f_EnIOut(:,end), N_PreSynPix, L4SE,L4SI, L4CE,L4CI, L4IE,L4II,...
                                 S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                 S_EL6,S_IL6,rL6E,rL6I,S_amb,rE_amb,rI_amb,...%7 L6 Amb                                   
                                 lgn_EOnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,...
                                 tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                 rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                 gL_E,gL_I,Ve,Vi,LIFSimuT, dt); % No more Freq                                   
%% NOW! consider the previous mVs if it already satisfies steady condition
if loop>50
    f_EnI0 = mean(f_EnIOut(:,end-10+1:end)) * 0.9 + FrLIF*0.1;
else
     f_EnI0 = FrLIF;
end

                               
loop = loop+1;

f_EnIOut = [f_EnIOut,f_EnI0];
%meanVs = [meanVs, mVIn];

if ~SteadyIndicate
  if loop >= TestPoints && max(std(f_EnIOut(:,end-TestPoints+1:end),0,2) ...
                            ./mean(f_EnIOut(:,end-TestPoints+1:end),2))<0.1 % std/mean for the last 10 samples
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

if N_Hyp > 0 && strcmpi(HyperPara(1),'Mean') % Specify: I don't need the trajectory but the final ones    
    f_EnIOut = mean(f_EnIOut(:,end-50:end),2);
    %meanVs   = mean(meanVs(:,end-50:end),2);
end
end


% single cell simulation: to collect mean V (and firing rates)

function [mVLIF,FrLIF] = LIF1Pixel(Fr_MFinv, N_PreSynPix, L4SE,L4SI, L4CE,L4CI, L4IE,L4II,...
                                 S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                 S_EL6,S_IL6,rL6E,rL6I,S_amb,rE_amb,rI_amb,...%7 L6 Amb                                   
                                 lgn_EOnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,...
                                 tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                 rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                 gL_E,gL_I,Ve,Vi,LIFSimuT, dt) % No more Freq
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
LgnInputNow = [lgn_EOnOff(1)*NlgnS; 
               lgn_EOnOff(1)*NlgnC;
               lgn_EOnOff(2)*NlgnS; 
               lgn_EOnOff(2)*NlgnC;
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

% figure
% subplot 221; plot(tt, GAMPA(3,1:length(tt))); xlim([T*9/10,T]); xlabel('T(ms)'); ylabel('GAMPA')
% subplot 222; plot(tt, GNMDA(3,1:length(tt))); xlim([T*9/10,T]); xlabel('T(ms)'); ylabel('GNMDA')
% subplot 223; plot(tt, GGABA(3,1:length(tt))); xlim([T*9/10,T]); xlabel('T(ms)'); ylabel('GGABA')
% subplot 224; plot(linspace(0,), vRecord(3,:)); xlim([T*9/10,T]); xlabel('T(ms)'); ylabel('V')
end