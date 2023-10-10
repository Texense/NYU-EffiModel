%% MF+v for the field. Ver1: LGN on/off not considered so far.
% Input: NW parameters, 
%        connectivity matrices,
%        External input rates
%        Pixel setups
% Output: FrE, FrI, mVE, mVI

% Note: Let's not use anything to speed up linear inverse for now.
%% Functions quoting outside this file:
% 1. function [SpaPixMat,NnEPixel] = NeuVec2Pixel(FrE,NnE,NPxX,NPxY)
function [FrEPixVecOut, FrIPixVecOut, ...
          mVEPixVecOut, mVIPixVecOut] = MFvField_NoCorrection(C_EE_Fix,C_EI_Fix,C_IE_Fix,C_II_Fix,...     %4
                                                        S_EE,S_EI,S_IE,S_II,p_EEFail,...            %5
                                                        lambda_EOn_drive, lambda_EOff_drive, LGNFreq, S_Elgn,rE_amb,S_amb,...      %6 lambda_E drive is vector for each cell
                                                        lambda_I_drive,                               S_Ilgn,rI_amb,...            %3
                                                        S_EL6,S_IL6,rE_L6_Drive,rI_L6_Drive,...     %4 
                                                        ... % L6 drives are vectors for each cell. Vital since length(rQL6) for NE and NI
                                                        tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,... %7
                                                        rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,... %4
                                                        gL_E,gL_I,Ve,Vi,...                         %4
                                                        N_HC,NnE,NnI,NPixX, NPixY, varargin)        %5+x % if varagin non empty, we only export the last state
                                                        %This line: spatial indexes of E and I neurons. Numbers of pixels
%% Hyperparameters for iterations
if nargin > 41+1  % Specify: The number of loops I want, after stopping criteria met    
    AveLoop = varargin{2};
else 
    AveLoop = 100;
end

if nargin > 42+1 % Specify: Maximum nubmer of loops before stopping
        StopLoop = varargin{3};
else 
        StopLoop = AveLoop;
end

if nargin > 43+1 % Specify: stepsize h
    h_Step = varargin{4};
else
    h_Step = 1;        
end

if nargin > 44+1 % Specify: LIF simulation timef_pre
    LIFSimuT = varargin{5};
else
    LIFSimuT = 20*1e3;    % unit in ms    
end
%% Get matrices for pixels
% First, get connectivity matrix between Pixels
PixNum = NPixX*N_HC*NPixY*N_HC; % Number of pixels
C_EE_Pixel = zeros(PixNum);
C_EI_Pixel = zeros(PixNum);
C_IE_Pixel = zeros(PixNum);
C_II_Pixel = zeros(PixNum);
[~,NnEPixel] = NeuVec2Pixel(zeros(length(NnE.X),1),NnE,NPixX*N_HC,NPixY*N_HC);
[~,NnIPixel] = NeuVec2Pixel(zeros(length(NnI.X),1),NnI,NPixX*N_HC,NPixY*N_HC);
tic
parfor PixIndPost = 1:PixNum
    for PixIndPre = 1:PixNum
        C_EE_Pixel(PixIndPost,PixIndPre) = mean(sum(C_EE_Fix(NnEPixel.Vec == PixIndPost, NnEPixel.Vec == PixIndPre),2));
        C_EI_Pixel(PixIndPost,PixIndPre) = mean(sum(C_EI_Fix(NnEPixel.Vec == PixIndPost, NnIPixel.Vec == PixIndPre),2));
        C_IE_Pixel(PixIndPost,PixIndPre) = mean(sum(C_IE_Fix(NnIPixel.Vec == PixIndPost, NnEPixel.Vec == PixIndPre),2));
        C_II_Pixel(PixIndPost,PixIndPre) = mean(sum(C_II_Fix(NnIPixel.Vec == PixIndPost, NnIPixel.Vec == PixIndPre),2));
    end    
end
toc
C_EE_Pixel_Us = full(C_EE_Pixel);
C_IE_Pixel_Us = full(C_IE_Pixel);nonlinear dynamics
C_EI_Pixel_Us = full(C_EI_Pixel);
C_II_Pixel_Us = full(C_II_Pixel);
% Then get input vec for each pixel
lambda_EOn_Pixel = zeros(size(mVEPixVec));
lambda_EOff_Pixel = zeros(size(mVEPixVec));
L6E_Pixel = zeros(size(mVEPixVec));
L6I_Pixel = zeros(size(mVEPixVec));
for PixInd = 1:PixNum
    lambda_EOn_Pixel(PixInd)  = mean(lambda_EOn_drive(NnEPixel.Vec == PixInd));
    lambda_EOff_Pixel(PixInd) = mean(lambda_EOff_drive(NnEPixel.Vec == PixInd));
    L6E_Pixel(PixInd) = mean(rE_L6_Drive(NnEPixel.Vec == PixInd));
    L6I_Pixel(PixInd) = mean(rI_L6_Drive(NnIPixel.Vec == PixInd));
end
lambda_I_Pixel = lambda_I_drive;
% Then initialize FrPix and mVPix vectors
% Start from uniform guessed values
FrEIni = 16; FrIIni = 64;
FrEPixIni = FrEIni *ones(PixNum,1);
FrIPixIni = FrIIni *ones(PixNum,1);
mVEPixVec = 0.7 *ones(PixNum,1);
mVIPixVec = 0.8 *ones(PixNum,1);
% L4EInputIni = [C_EE_Pixel_Us*FrEPixIni*(1-p_EEFail);C_IE_Pixel_Us*FrEPixIni]/1e3;   
% L4IInputIni = [C_EI_Pixel_Us*FrIPixIni           ;  C_II_Pixel_Us*FrIPixIni]/1e3;   
%% PreSet Input events to Pixel-LIF neurons
LGNCurInp = 0; L6CurInp = 0;
dt = 0.1; TimeFrac = 1e3/LIFSimuT;
lambda_E_Pixel = (lambda_EOn_Pixel + lambda_EOff_Pixel)/2;
tic
lgnEPix_Events = PoissonInputForNetwork(PixNum,lambda_E_Pixel*(1-LGNCurInp),LIFSimuT*TimeFrac,dt);
lgnIPix_Events = PoissonInputForNetwork(PixNum,lambda_I_Pixel*(1-LGNCurInp),LIFSimuT*TimeFrac,dt);
AmbEPix_Events = PoissonInputForNetwork(PixNum,rE_amb,LIFSimuT*TimeFrac,dt);
AmbIPix_Events = PoissonInputForNetwork(PixNum,rI_amb,LIFSimuT*TimeFrac,dt);
L6EPix_Events  = PoissonInputForNetwork(PixNum,L6E_Pixel*(1-L6CurInp),LIFSimuT*TimeFrac,dt);
L6IPix_Events  = PoissonInputForNetwork(PixNum,L6I_Pixel*(1-L6CurInp),LIFSimuT*TimeFrac,dt);
% Idea for L4: changing every loop, but we can always rescale from the
% original event counts in each bin
toc

%% Do the following iterations recursively
for Epoch = 1:10
    tic
% PH for MF: Need variable modification
    Fr_MFinv = MFgivV(S_EE,S_EI,S_IE,S_II,p_EEFail,...
                     S_Elgn,S_Ilgn, rE_amb,rI_amb,S_amb, S_EL6,S_IL6, ...
                     gL_E,gL_I, Ve,Vi,...    
                     C_EE_Pixel_Us, C_EI_Pixel_Us, C_IE_Pixel_Us, C_II_Pixel_Us,...
                     lambda_EOn_Pixel,lambda_EOff_Pixel,lambda_I_Pixel, L6E_Pixel, L6I_Pixel,...
                     mVEOnPixVec, mVEOffPixVec, mVIPixVec, PixNum, ...
                     FrEOnPixVec, FrEOffPixVec, FrIPixVec, tau_ref);
    
    if Epoch < 10
       [mVEOnPixNew,mVEOffPixNew,mVIPixNew] =  ...
           LIFPixels(Fr_MFinv, L4Pix_EventsEIni, ...
                     C_EE_Pixel_Us,C_EI_Pixel_Us,C_IE_Pixel_Us,C_II_Pixel_Us,...
                     S_EE,S_EI,S_IE,S_II,p_EEFail,...
                     lgnEPix_Events,S_Elgn,AmbEPix_Events,S_amb,...
                     lgnIPix_Events,S_Ilgn,AmbIPix_Events,...
                     S_EL6,S_IL6,L6EPix_Events,L6IPix_Events,...
                     tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,...
                     rhoE_ampa,rhoE_nmda,rhoI_ampa,rhoI_nmda,...
                     gL_E,gL_I,Ve,Vi,LIFSimuT, dt, PixNum,...
                     LGNCurInp, L6CurInp,...
                     lambda_E_Pixel,lambda_I_Pixel,L6E_Pixel,L6I_Pixel,...
                     lambda_EOn_Pixel, lambda_EOff_Pixel, LGNFreq);
       mVEOnPixVec = mVEOnPixVec*(1-h_Step) + mVEOnPixNew * h_Step;
       mVEOffPixVec = mVEOffPixVec*(1-h_Step) + mVEOffPixNew * h_Step;
       mVIPixVec = mVIPixVec*(1-h_Step) + mVIPixNew * h_Step;
   
    end
   toc
   figure(1)
   subplot 321
   ShowField(Fr_MFinv, 1:900, 30, 30 )
   caxis([0 40])
   subplot 322
   ShowField(mVEOnPixVec, 1:900, 30, 30 )
   caxis([0.64 0.77])
   subplot 323
   ShowField(Fr_MFinv, 901:1800, 30, 30 )
   caxis([0 40])
   subplot 324
   ShowField(mVEOffPixVec, 1:900, 30, 30 )
   caxis([0.64 0.77])
      subplot 323
   ShowField(Fr_MFinv, 1801:2700, 30, 30 )
   caxis([20 120])
   subplot 324
   ShowField(mVIPixVec, 1:900, 30, 30 )
   caxis([0.74 0.83])
   drawnow
end
end

%% MF Computing: giving mVs
function [Fr_MFinv] = MFgivV(S_EE,S_EI,S_IE,S_II,p_EEFail,...
                     S_Elgn,S_Ilgn, rE_amb,rI_amb,S_amb, S_EL6,S_IL6, ...
                     gL_E,gL_I, Ve,Vi,...    
                     C_EE_Pixel_Us, C_EI_Pixel_Us, C_IE_Pixel_Us, C_II_Pixel_Us,...
                     lambda_EOn_Pixel,lambda_EOff_Pixel,lambda_I_Pixel, L6E_Pixel, L6I_Pixel,...
                     mVEOnPixVec, mVEOffPixVec, mVIPixVec, PixNum, ...
                     FrEOnPixVec, FrEOffPixVec, FrIPixVec, tau_ref,...
                     FrPreUse, FrPreUseDim) % The last 2: Use preset Frs, and their entry location
    % Mats
    MatEE = (S_EE*(1-p_EEFail))*C_EE_Pixel_Us; % Need Ref here??
    MatEI = S_EI *              C_EI_Pixel_Us;
    MatIE = S_IE *              C_IE_Pixel_Us.*(Ve-repmat(mVIPixVec,1,PixNum));
    MatII = S_II *              C_II_Pixel_Us.*(Vi-repmat(mVIPixVec,1,PixNum));
    ConnMat = [MatEE.*(Ve-repmat(mVEOnPixVec,1,PixNum))/2 , MatEE.*(Ve-repmat(mVEOnPixVec,1,PixNum))/2 , MatEI.*(Vi-repmat(mVEOnPixVec,1,PixNum));
               MatEE.*(Ve-repmat(mVEOffPixVec,1,PixNum))/2, MatEE.*(Ve-repmat(mVEOffPixVec,1,PixNum))/2, MatEI.*(Vi-repmat(mVEOffPixVec,1,PixNum));
               MatIE/2                                    , MatIE/2                                    , MatII];
    % Leak On/Off
    LeakEOn  = gL_E * (0-mVEOnPixVec)  * 1e3;
    LeakEOff = gL_E * (0-mVEOffPixVec) * 1e3;
    LeakI =    gL_I * (0-mVIPixVec)    * 1e3;
    LeakV = [LeakEOn;LeakEOff;LeakI];
    
    % Ext
    ExtEOn  = (lambda_EOn_Pixel*S_Elgn  + rE_amb*S_amb + L6E_Pixel*S_EL6).*(Ve-mVEOnPixVec ) * 1e3;
    ExtEOff = (lambda_EOff_Pixel*S_Elgn + rE_amb*S_amb + L6E_Pixel*S_EL6).*(Ve-mVEOffPixVec) * 1e3;
    ExtI =    (lambda_I_Pixel*S_Ilgn    + rI_amb*S_amb + L6I_Pixel*S_IL6).*(Ve-mVIPixVec)    * 1e3;
    ExtV = [ExtEOn;ExtEOff;ExtI];
    % Ref Vec
    RefEOn = 1-FrEOnPixVec*tau_ref/1e3;
    RefEOff = 1-FrEOffPixVec*tau_ref/1e3;
    RefI = 1-FrIPixVec*tau_ref/1e3;
    RefM = sparse(diag([RefEOn;RefEOff;RefI]));
    
    % MF Equations:
    if isempty(FrPreUse)
    Fr_MFinv = (sparse(eye(3*PixNum))-RefM*ConnMat) \ (RefM * ( ExtV + LeakV));
    else
        FrPreVec = zeros(PixNum*3,1); FrPreVec(FrPreUseDim) = FrPreUse;
        DimAll = 1:PixNum*3; FrminusDim = DimAll(~ismember(DimAll,FrPreUse));
        FrPreReplace = RefM*ConnMat * FrPreVec; 
        % Now compute in the subspace
        FrPreMinus = FrPreReplace(FrminusDim);
        RefMMinus  = RefM(FrminusDim, FrminusDim);
        ConnMatMinus = ConnMat(FrminusDim, FrminusDim);
        ExtVMinus = ExtV(FrminusDim);
        LeakVMinus = LeakV(FrminusDim);
        Fr_MFinvMinus = (sparse(eye(length(FrminusDim)))-RefMMinus*ConnMatMinus) \ ...
                        (RefMMinus * ( ExtVMinus + LeakVMinus) + FrPreMinus);
        % make up MF output            
        Fr_MFinv = zeros(PixNum*3,1);   
        Fr_MFinv(FrminusDim) = Fr_MFinvMinus; Fr_MFinv(FrPreUseDim) = FrPreUse; 
    end


end

