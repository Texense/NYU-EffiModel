%% LIFonly for Pixels. E is divided into S and C cells
%% Input: 
%        Th idea of this is similar to MF+v, but all L4 input 

%        LGNFreq                  on/off frequencies of LFN input
%        L4E, L4I                 Total E/I input kick numbers
%        Assume L4SE/CE/IE are proportional to presynaptic neuron numbers
%        S_EE,S_EI,S_IE,S_II      synaptic strength
%        S_EL6,S_IL6,rE_L6,rI_L6, L6 Input, All Scalar
%        S_amb,rS_amb,rC_amb,rI_amb      amb Input, All Scalar
%        p_EEFail                 E-to-E Synaptic failure prob
%        lgn_S/COnOff,lgn_I       LGN input rate. Phase OnOff for E (2*1 vec) 
%        S_Elgn S_Ilgn            LGN strength of drives
%        NlgnS,NlgnC,NlgnI        LGN neuron # to SCI cells
%        gL_E,gL_I                Leaky time constants
%        Ve,Vi                    Reversal potentials
%        HyperPara           Cell. Entries are simulation hyperparameters
%                            1th: LIF simulation time
%                            2th an 3th: testmode and test value: threshold
%                            or delay time
%                            
% Output:f_EnI               Estimation of firing rates, E and I; A sequences
%        meanVs              mean V of E and I

%% Ver 5:       For single pixels
%% Ver 6:       Bg features added: lgn_SU and lgn_CU are 2X1 vector (Then drive) or scalar (Then bg)?
% iteration
% Zhuo-Cheng Xiao 05/03/2022
function [f_EnIOut]...
           = LIFo_SinglePixel(...
...% MF Parameters                     
                     LGNFreq, L4SE,L4SI, L4CE,L4CI, L4IE,L4II,... %3 
                     S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                     S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                     lgn_SU, lgn_CU,lgn_I,N_Slgn,N_Clgn,N_Ilgn, S_Elgn,S_Ilgn,... %7
...% A big problem:  lgn_COnOff,lgn_I was unquoted - so using 45 Hz by defalut                    
                     gL_E,gL_I,Ve,Vi, tau_ref,... %5
...% Below are LIF details
                     tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,... %7
                     rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                     HyperPara) % if varagin non empty, we only export the last state
%% Hyperparameters
N_Hyp = length(HyperPara);

if N_Hyp > 0 % Specify: LIF simulation timef_pre
    LIFSimuT = HyperPara{1};
else
    LIFSimuT = 50*1e3;    % unit in ms    
end

if N_Hyp > 1 % Specify: Starting point of Recording
    if strcmpi(HyperPara{2},'thre')
        RecdThre = HyperPara{3};
        RecdDely = 0;
    elseif strcmpi(HyperPara{2},'delay')
        RecdThre = 0;
        RecdDely = HyperPara{3};
    else
        disp('***Wrong test mode')
        return
    end
else
    RecdThre = 0;    % unit in ms  
    RecdDely = 0;
end
dt = 0.1;

%% All input events
% Note: If rate is per second, need to be divided by 1e3 to rescale to ms
TimeFrac = 1/(LIFSimuT/1e3)/LGNFreq;

AmbPix_Events = PoissonInputForNetwork(3,[rS_amb;rC_amb;rI_amb],                          LIFSimuT*TimeFrac,dt,true); % rQ_amb are per ms
L6Pix_Events  = PoissonInputForNetwork(3,[rL6SU; rL6CU; rL6IU],                           LIFSimuT*TimeFrac,dt,true); % rL6 are per ms
L4Pix_EventsE = PoissonInputForNetwork(3,[L4SE*(1-p_EEFail); L4CE*(1-p_EEFail); L4IE]/1e3,LIFSimuT*TimeFrac,dt,true); % Need to count EE failure here 
L4Pix_EventsI = PoissonInputForNetwork(3,[L4SI;              L4CI;              L4II]/1e3,LIFSimuT*TimeFrac,dt,true); % L4E/I are per second
gL = [gL_E*ones(2,1); gL_I*ones(1)];

% LGN input pissons are handeled specifically:
if length(lgn_SU) == 2 % drive
    disp('Drive regime. Use Phase variant LGN input')
    F_LGN1ort = 45/1e3;
    lgn_Events = PoissonInputForNetwork(3,[N_Slgn;N_Clgn;N_Ilgn]*F_LGN1ort,             LIFSimuT*TimeFrac,dt,true);
    lgn_Dr_dflt = F_LGN1ort*ones(2,1);
    % Adjust S lgn according to phase. Do the same thing to C 
    lgnSAdj = lgn_SU./lgn_Dr_dflt;
    lgnS_AdjU = reshape(repmat(lgnSAdj,1,floor(size(L4Pix_EventsE,2)/2))', 1, 2*floor(size(L4Pix_EventsE,2)/2));
    
    lgnCAdj = lgn_CU./lgn_Dr_dflt;
    lgnC_AdjU = reshape(repmat(lgnCAdj,1,floor(size(L4Pix_EventsE,2)/2))', 1, 2*floor(size(L4Pix_EventsE,2)/2));
    % I will be same as ort, so actually lgnI input is redundant.
elseif length(lgn_SU) == 1 % bg
    disp('Background regime. Use Phase consistent LGN input')
    lgn_Events = PoissonInputForNetwork(3,[N_Slgn;N_Clgn;N_Ilgn].*[lgn_SU;lgn_CU;lgn_I],LIFSimuT*TimeFrac,dt,true); % lgn_Q should be per ms
else
    disp('***Illigal LGN input. Returning...')
end

%% Each on/off is a cycle for LIF. For each cyc, initate a new start.
TCyc = floor(1/TimeFrac);
VRcdCyc = cell(TCyc,1);
SpRcdCyc = cell(TCyc,1);

RefTimer = zeros(3, 1); % initiate the timer of ref
vt = zeros(3, 1); 

for TInt = 1:TCyc
    VRcd = zeros(size(L4Pix_EventsE));
    SpLIF = zeros(3, 1);

    % New Gs
    lgnPix_Events = lgn_Events(:,randperm(size(lgn_Events,2)));
    if length(lgn_SU) == 2 % only adjust in drive regime
        lgnPix_Events(1,1:length(lgnS_AdjU)) = lgnPix_Events(1,1:length(lgnS_AdjU)).*lgnS_AdjU; % first row: S
        lgnPix_Events(2,1:length(lgnC_AdjU)) = lgnPix_Events(2,1:length(lgnC_AdjU)).*lgnC_AdjU; % second row: C
    end

    winAMPA = 0:dt:tau_ampa_D*3;
    winNMDA = 0:dt:tau_nmda_D*3;
    winGABA = 0:dt:tau_gaba_D*3;
    KerAMPA = 1/(tau_ampa_D-tau_ampa_R) * (exp(-winAMPA/tau_ampa_D) - exp(-winAMPA/tau_ampa_R));
    KerNMDA = 1/(tau_nmda_D-tau_nmda_R) * (exp(-winNMDA/tau_nmda_D) - exp(-winNMDA/tau_nmda_R));
    KerGABA = 1/(tau_gaba_D-tau_gaba_R) * (exp(-winGABA/tau_gaba_D) - exp(-winGABA/tau_gaba_R));
    KerAMPA = KerAMPA / (sum(KerAMPA)*dt);
    KerNMDA = KerNMDA / (sum(KerNMDA)*dt);
    KerGABA = KerGABA / (sum(KerGABA)*dt);
    % Incorporate L4 and other input
    ampaInp = single([S_Elgn * lgnPix_Events(1:2,:);            S_Ilgn * lgnPix_Events(3,:)]...
        +  S_amb  * AmbPix_Events ... % S_amb identical for E and I
        + [S_EL6  * L6Pix_Events(1:2,:) * rhoE_ampa; S_IL6 * L6Pix_Events(3,:) * rhoI_ampa]...
        + [S_EE * L4Pix_EventsE(1:2,:) * rhoE_ampa; S_IE * L4Pix_EventsE(3,:) * rhoI_ampa]);
    nmdaInp = single([S_EL6  * L6Pix_Events(1:2,:) * rhoE_nmda; S_IL6 * L6Pix_Events(3,:) * rhoI_nmda]...
        + [S_EE * L4Pix_EventsE(1:2,:) * rhoE_nmda; S_IE * L4Pix_EventsE(3,:) * rhoI_nmda]);
    gabaInp = single([S_EI * L4Pix_EventsI(1:2,:);             S_II * L4Pix_EventsI(3,:)] );
    %% Precompute All Conductances: Speed up by ifft/fft
    n = size(ampaInp,2);
    GAMPA = ifft(fft(ampaInp',n) .* repmat(fft(KerAMPA',n),1,3))';
    GNMDA = ifft(fft(nmdaInp',n) .* repmat(fft(KerNMDA',n),1,3))';
    GGABA = ifft(fft(gabaInp',n) .* repmat(fft(KerGABA',n),1,3))';

    GE = GAMPA + GNMDA;  GI = GGABA;
    
    for tInd = 1:size(L4Pix_EventsE,2)
        % Firstly, refrectory neurons get out if timer reachs t_ref
        RefTimer(isnan(vt)) = RefTimer(isnan(vt)) + dt;
        vt(RefTimer>=tau_ref) = 0;
        RefTimer(RefTimer>=tau_ref) = 0;

        ov  = vt + dt*(-gL.*vt - GE(:,tInd).*(vt-Ve) - GI(:,tInd).*(vt-Vi));
        SpNow = single(ov>=1);
        ov(ov>=1) = nan;
        VRcd(:,tInd) = ov;
        % Compute meanV
            SpLIF = SpLIF + SpNow;
        % Transfer variables
        vt = ov;
    end
    %VRcdCyc{TInt} = VRcd;
    SpRcdCyc{TInt} = SpLIF;
end

%% I am not concerning the Vs now, but maybe need that in future
%VRcdMat = cell2mat(VRcdCyc')';
SpRcdMat = cell2mat(SpRcdCyc')';
f_EnIOut = SpRcdMat'*LGNFreq;
%FrRcdAve(LIFtestInd,:) = mean(SpRcdMat(10:end,:))*LGNFreq;


end





