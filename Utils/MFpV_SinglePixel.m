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
%        S_amb,rS_amb,rC_amb,rI_amb      amb Input, All Scalar
%        p_EEFail                 E-to-E Synaptic failure prob
%        lgn_S/COnOff,lgn_I       LGN input rate. Phase OnOff for E (2*1 vec) 
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
%                            6th an 7th: testmode and test value: threshold
%                            or delay time
%                            
% Output:f_EnI               Estimation of firing rates, E and I; A sequences
%        meanVs              mean V of E and I
%        loop                Number of loops
%        SteadyIndicate      logical value for convergence
%% Ver 5:       For single pixels
% iteration
% Zhuo-Cheng Xiao 07/27/2021
function [f_EnIOut,meanVs,loop,SteadyIndicate,FailureIndicate]...
           = MFpV_SinglePixel(...
...% MF Parameters                     
                     N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,... %3 
                     S_EE,S_EI,S_IE,S_II,p_EEFail,... %5
                     S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                     lgn_SOnOff, lgn_COnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... %7
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

if N_Hyp > 5 % Specify: Starting point of Recording
    if strcmpi(HyperPara{6},'thre')
        RecdThre = HyperPara{7};
        RecdDely = 0;
    elseif strcmpi(HyperPara{6},'delay')
        RecdThre = 0;
        RecdDely = HyperPara{7};
    else
        disp('***Wrong test mode')
        return
    end
else
    RecdThre = 0;    % unit in ms  
    RecdDely = 0;
end
dt = 0.1;
% initialize with a reasonable guess     
mVSOn = 0.67; mVSOff = 0.57; 
mVCOn = 0.67; mVCOff = 0.57;
mVI = 0.77; 
f_pre = ones(5,1);
meanVs = [mVSOn;mVSOff;mVCOn;mVCOff;mVI];
f_SCIIni = MF_SCI_1Pix(N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                       S_EE,S_EI,S_IE,S_II,p_EEFail,...
                       S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                       lgn_SOnOff, lgn_COnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... %7
                       gL_E,gL_I,Ve,Vi, tau_ref, meanVs, f_pre,f_pre);
f_EnIOut = f_SCIIni;

loop = 0;
SteadyCounter = 0; % Indicate the number of loops after steady condition
SteadyIndicate = false;
TestPoints = floor(7); % How many consecutive points we test
%while( norm([mVEpre;mVIpre] - [mVE;mVI]) > 0.01 || norm(f_EnIpre - f_EnI0)>0.1) %%% relative difference for firing rates!!

FailureIndicate = 0;
Suspicious = false;
while SteadyCounter<AveLoop %the formal ending condition
% simulate one neuron with input 
if max(f_EnIOut(:,end))>501
    disp('Large Fr >501Hz. Break')
    FailureIndicate = 1;
    break 
end
[mVLIF,FrLIF] = LIF1Pixel(f_EnIOut(:,end), N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                                 S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                 S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                                 lgn_SOnOff, lgn_COnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,...
                                 tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                 rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                 gL_E,gL_I,Ve,Vi,LIFSimuT, dt, RecdThre, RecdDely); % No more Freq                                   
Fail = sum(isnan(mVLIF))>0; 
if Fail
    disp('NaN appears. Break')
    FailureIndicate = 1;
    break    
end                             
%% NOW! consider the previous mVs if it already satisfies steady condition
if loop>100
    mVIn = mean(meanVs(:,end-10+1:end),2) * 0.9 + mVLIF*0.1; %%% WAS WRONG!!!
else
     mVIn = mVLIF;
end

% estimate with ref now
f_EnI0 = MF_SCI_1Pix(N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                     S_EE,S_EI,S_IE,S_II,p_EEFail,...
                     S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                     lgn_SOnOff, lgn_COnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... % lgn
                     gL_E,gL_I,Ve,Vi, tau_ref, mVIn, f_SCIIni, FrLIF); 
                               
                               
Suspicious = (min(f_EnI0)<-1);                              
%f_EnI0 = max([f_EnI0,[0;0]],[],2);
%f_EnI0 = abs(f_EnI0);                                     
%% The new input!
f_SCIIni = f_EnI0*h_Step + f_SCIIni*(1-h_Step);

loop = loop+1;

f_EnIOut = [f_EnIOut,f_EnI0];
meanVs = [meanVs, mVIn];

if ~SteadyIndicate
  if loop >= TestPoints && mean(std(f_EnIOut(:,end-TestPoints+1:end),0,2) ...
                            ./mean(f_EnIOut(:,end-TestPoints+1:end),2))<0.1 % std/mean for the last 10 samples
     SteadyIndicate = true;
  end  
else 
    SteadyCounter = SteadyCounter+1;
end

if ((loop>StopLoop || SteadyIndicate) && Suspicious) 
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
    meanVs   = mean(meanVs(:,end-50:end),2);
end
end



%% mean-field est With ref. We use parameters and mean Vs to estimate firing rates
% Input: N_PreSynPix         connectivity matrix for SCI within pixel. 3*3
%        S_EE,S_EI,S_IE,S_II synaptic strength
%        S_EL6,S_IL6,rE_L6,rI_L6 L6 Input
%        p_EEFail            Synaptic failure prob
%        lambda_E lambda_I   LGN input
%        rS_amb,rC_amb rI_amb       Ambient input
%        S_Elgn S_Ilgn S_amb Synaptic strength of drives
%        gL_E,gL_I           Leaky time constants
%        Ve,Vi               Reversak potentials
%        mVE,mVI             Mean Vs, collected from simulation
%        tau_ref             ref period, in ms
%        f_pre               from previous
% Output:Fr_MFinv            Estimation of firing rates, Son, Soff, Con, Coff, I

function Fr_MFinv = MF_SCI_1Pix(N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                             S_EE,S_EI,S_IE,S_II,p_EEFail,...
                             S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                             lgn_SOnOff, lgn_COnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,... % lgn
                             gL_E,gL_I,Ve,Vi, tau_ref, meanVs, f_pre, FrLIF)                             
%% PreSynaptic neuron numbers. External and internal (pixel)
N_SSi = N_PreSynPix(1,1); N_CSi = N_PreSynPix(1,2); N_ISi = N_PreSynPix(1,3); %% REally?? Row/Col
N_SCi = N_PreSynPix(2,1); N_CCi = N_PreSynPix(2,2); N_ICi = N_PreSynPix(2,3); 
N_SIi = N_PreSynPix(3,1); N_CIi = N_PreSynPix(3,2); N_IIi = N_PreSynPix(3,3); 

% Downplay current by a ref factor
% Use LIF firing rates if negative value emerges
f_preU = f_pre; 
f_preU(f_pre<0) = FrLIF(f_pre<0);
%f_preU(f_preU>200) = 200; % give a maximum
RefM = diag(1-f_preU*tau_ref/1000);
% Mats
MatSS = (S_EE*(1-p_EEFail))*N_SSi; MatCS = (S_EE*(1-p_EEFail))*N_CSi; MatIS = S_IE * N_ISi;
MatSC = (S_EE*(1-p_EEFail))*N_SCi; MatCC = (S_EE*(1-p_EEFail))*N_CCi; MatIC = S_IE * N_ICi;
MatSI = S_EI *              N_SIi; MatCI = S_EI *              N_CIi; MatII = S_II * N_IIi;

eVSOn  = (Ve-meanVs(1));  iVSOn  = (Vi-meanVs(1));
eVCOn  = (Ve-meanVs(2));  iVCOn  = (Vi-meanVs(2));
eVSOff = (Ve-meanVs(3));  iVSOff = (Vi-meanVs(3));
eVCOff = (Ve-meanVs(4));  iVCOff = (Vi-meanVs(4));
eVI    = (Ve-meanVs(5));  iVI    = (Vi-meanVs(5));
% post/pre SOn              COn              SOff             COff             I
ConnMat = [MatSS/2.*eVSOn,  MatSC/2.*eVSOn,  MatSS/2.*eVSOn,  MatSC/2.*eVSOn,  MatSI.*iVSOn;  % SOn
           MatCS/2.*eVCOn,  MatCC/2.*eVCOn,  MatCS/2.*eVCOn,  MatCC/2.*eVCOn,  MatCI.*iVCOn;  % COn
           MatSS/2.*eVSOff, MatSC/2.*eVSOff, MatSS/2.*eVSOff, MatSC/2.*eVSOff, MatSI.*iVSOff; % SOff
           MatCS/2.*eVCOff, MatCC/2.*eVCOff, MatCS/2.*eVCOff, MatCC/2.*eVCOff, MatCI.*iVCOff; % COff
           MatIS/2.*eVI,    MatIC/2.*eVI,    MatIS/2.*eVI,    MatIC/2.*eVI,    MatII.*iVI];   % I
% Leak On/Off
LeakSOn  = gL_E * (0-meanVs(1)) * 1e3;
LeakCOn  = gL_E * (0-meanVs(2)) * 1e3;
LeakSOff = gL_E * (0-meanVs(3)) * 1e3;
LeakCOff = gL_E * (0-meanVs(4)) * 1e3;
LeakI =    gL_I * (0-meanVs(5)) * 1e3;
LeakV = [LeakSOn;
         LeakCOn;
         LeakSOff;
         LeakCOff;
         LeakI];
% Ext
ExtSOn  = (lgn_SOnOff(1)*NlgnS*S_Elgn + rS_amb*S_amb + rL6SU*S_EL6)*eVSOn *1e3 ...
        + (L4SIU*S_EI)*iVSOn + L4SEU*S_EE*eVSOn*(1-p_EEFail); 
ExtCOn  = (lgn_COnOff(1)*NlgnC*S_Elgn + rC_amb*S_amb + rL6CU*S_EL6)*eVCOn *1e3 ...
        + (L4CIU*S_EI)*iVCOn + L4CEU*S_EE*eVCOn*(1-p_EEFail);
ExtSOff = (lgn_SOnOff(2)*NlgnS*S_Elgn + rS_amb*S_amb + rL6SU*S_EL6)*eVSOff *1e3...
        + (L4SIU*S_EI)*iVSOn + L4SEU*S_EE*eVSOff*(1-p_EEFail);
ExtCOff = (lgn_COnOff(2)*NlgnC*S_Elgn + rC_amb*S_amb + rL6CU*S_EL6)*eVCOff *1e3...
        + (L4CIU*S_EI)*iVCOn + L4CEU*S_EE*eVCOff*(1-p_EEFail);
ExtI =    (lgn_I        *NlgnI*S_Ilgn + rI_amb*S_amb + rL6IU*S_IL6)*eVI    *1e3...
        + (L4IIU*S_II)*iVI   + L4IEU*S_IE*eVI;
ExtV = [ExtSOn;
        ExtCOn;
        ExtSOff;
        ExtCOff
        ExtI];
    
Fr_MFinv = (eye(5)-RefM*ConnMat) \ (RefM * ( ExtV + LeakV)); 
%Fr_MFinv(Fr_MFinv<0) = 3;
% FrPreUseDim = Fr_MFinv<0; 
% FrPreUse = FrLIF(Fr_MFinv<0);
% if isempty(FrPreUse)
%     return
% else
%     Fr_MFinv(FrPreUseDim) = FrPreUse;
% %     FrPreVec = zeros(5,1); FrPreVec(FrPreUseDim) = FrPreUse;
% %     DimAll = 1:5; FrminusDim = DimAll(~ismember(DimAll,FrPreUseDim));
% %     FrPreReplace = RefM*ConnMat * FrPreVec;
% %     % Now compute in the subspace
% %     FrPreMinus = FrPreReplace(FrminusDim);
% %     RefMMinus  = RefM(FrminusDim, FrminusDim);
% %     ConnMatMinus = ConnMat(FrminusDim, FrminusDim);
% %     ExtVMinus = ExtV(FrminusDim);
% %     LeakVMinus = LeakV(FrminusDim);
% %     Fr_MFinvMinus = (eye(length(FrminusDim))-RefMMinus*ConnMatMinus) \ ...
% %         (RefMMinus * ( ExtVMinus + LeakVMinus) + FrPreMinus);
% %     % make up MF output
% %     Fr_MFinv = zeros(5,1);
% %     Fr_MFinv(FrminusDim) = Fr_MFinvMinus; Fr_MFinv(FrPreUseDim) = FrPreUse;
% end

end

% single cell simulation: to collect mean V (and firing rates)

function [mVLIF,FrLIF] = LIF1Pixel(Fr_MFinv, N_PreSynPix, L4SEU,L4SIU, L4CEU,L4CIU, L4IEU,L4IIU,...
                                 S_EE,S_EI,S_IE,S_II,p_EEFail,...
                                 S_EL6,S_IL6,rL6SU,rL6CU,rL6IU,S_amb,rS_amb,rC_amb,rI_amb,...%7 L6 Amb                                   
                                 lgn_SOnOff, lgn_COnOff,lgn_I,NlgnS,NlgnC,NlgnI, S_Elgn,S_Ilgn,...
                                 tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                                 rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                                 gL_E,gL_I,Ve,Vi,LIFSimuT, dt, RecdThre, RecdDely) % No more Freq
                               % Last two lines in case we have used current input...
                               %% First check if L4 rates match neuron numbers
if length(Fr_MFinv) == 5 % Son Con Soff Coff, I
    Fr_MFinv(Fr_MFinv<0) = eps;
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
L4SEUAll = (L4SEU + rS*N_SSi + rC*N_SCi)*(1-p_EEFail); L4SIUAll = L4SIU + rI*N_SIi; %fail multiply here. NOWHERE ELSE!
L4CEUAll = (L4CEU + rS*N_CSi + rC*N_CCi)*(1-p_EEFail); L4CIUAll = L4CIU + rI*N_CIi;
L4IEUAll =  L4IEU + rS*N_ISi + rC*N_ICi;               L4IIUAll = L4IIU + rI*N_IIi;
L4EInputNow = [L4SEUAll; L4CEUAll; L4SEUAll; L4CEUAll; L4IEUAll]/1000; 
L4IInputNow = [L4SIUAll; L4CIUAll; L4SIUAll; L4CIUAll; L4IIUAll]/1000; 
%Prevent too large input
if max([L4EInputNow;L4IInputNow])>40000
    mVLIF = nan;
    FrLIF = nan;
    return
end


L4Pix_EventsE = PoissonInputForNetwork(5,L4EInputNow,LIFSimuT*TimeFrac,dt*dtMltp); % E to everyone
L4Pix_EventsI = PoissonInputForNetwork(5,L4IInputNow,LIFSimuT*TimeFrac,dt*dtMltp); % I to everyone 
% LGN, amb, L6
LgnInputNow = [lgn_SOnOff(1)*NlgnS; 
               lgn_COnOff(1)*NlgnC;
               lgn_SOnOff(2)*NlgnS; 
               lgn_COnOff(2)*NlgnC;
               lgn_I        *NlgnI];
lgnPix_Events = PoissonInputForNetwork(5,LgnInputNow,                         LIFSimuT*TimeFrac,dt*dtMltp);
AmbPix_Events = PoissonInputForNetwork(5,[rS_amb;rC_amb;rS_amb;rC_amb;rI_amb],LIFSimuT*TimeFrac,dt*dtMltp);
L6Pix_Events  = PoissonInputForNetwork(5,[rL6SU; rL6CU; rL6SU; rL6CU; rL6IU], LIFSimuT*TimeFrac,dt*dtMltp);
% Leakage
gL = [gL_E*ones(4,1); gL_I*ones(1)];
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
               +  S_amb  * AmbPix_Events ... % S_amb identical for E and I
               + [S_EL6  * L6Pix_Events(1:4,:) * rhoE_ampa; S_IL6 * L6Pix_Events(5,:) * rhoI_ampa]...
               + [S_EE * L4Pix_EventsE(1:4,:) * rhoE_ampa; S_IE * L4Pix_EventsE(5,:) * rhoI_ampa]);
nmdaInp = single([S_EL6  * L6Pix_Events(1:4,:) * rhoE_nmda; S_IL6 * L6Pix_Events(5,:) * rhoI_nmda]...
               + [S_EE * L4Pix_EventsE(1:4,:) * rhoE_nmda; S_IE * L4Pix_EventsE(5,:) * rhoI_nmda]);
gabaInp = single([S_EI * L4Pix_EventsI(1:4,:);             S_II * L4Pix_EventsI(5,:)] );
%% Precompute All Conductances: Speed up by ifft/fft
n = size(ampaInp,2);
GAMPA = ifft(fft(ampaInp',n) .* repmat(fft(KerAMPA',n),1,5))';
GNMDA = ifft(fft(nmdaInp',n) .* repmat(fft(KerNMDA',n),1,5))';
GGABA = ifft(fft(gabaInp',n) .* repmat(fft(KerGABA',n),1,5))';

GE = GAMPA + GNMDA;  GI = GGABA;
%% Evolve single neurons
% Initialize records % sumV = zeros(5, 1);% VNanCount = zeros(5, 1);
vt = zeros(5, 1);
%vRecord = zeros(5, length(tt));
SpLIF = zeros(5, 1);
RefTimer = zeros(5, 1);
RecordInd = 0;
InpWin = 5;
mVs = zeros(5,floor(T/InpWin));
NeuNum = 5;
cellAdj = [2,4,5];
% If same for on/off., then cancel evolving off phases
if lgn_SOnOff(1)==lgn_SOnOff(2) && lgn_COnOff(1)==lgn_COnOff(2)
    NeuNum = 3; cellAdj = [2,3];
    GE(3:4,:) = []; GI(3:4,:) = []; gL(3:4,:) = [];
    mVs(3:4,:) = []; vt(3:4,:) = []; SpLIF(3:4,:) = [];
    RefTimer(3:4,:) = [];
end

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
        vRecord = zeros(NeuNum,FrameNum);
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
        GridDly = floor(RecdDely/dt);
        for cellInd = cellAdj
            % either recording thresold
            Vt = vRecord(cellInd,:);
            if RecdThre>0
            Vt(Vt< RecdThre & Vt >0) = nan; 
            end
            GridRef = find(isnan(Vt));
            % or recording delay
            if ~isempty(GridRef) && GridDly>0 % Only do more nan if spikes and nontrivial dly
                Vt(unique(reshape(GridRef + (0: GridDly)',1,length(GridRef)*(GridDly+1)))) = nan;
            end
            vRecord(cellInd,:) = Vt(1:length(vRecord(cellInd,:)));
        end
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
if lgn_SOnOff(1)==lgn_SOnOff(2) && lgn_COnOff(1)==lgn_COnOff(2)
    mVLIF(5) = mVLIF(3);mVLIF(3) = mVLIF(1); mVLIF(4) = mVLIF(2);
    FrLIF(5) = FrLIF(3);FrLIF(3) = FrLIF(1); FrLIF(4) = FrLIF(2);
end

%mVLIF(isnan(mVLIF)) = 0;
end