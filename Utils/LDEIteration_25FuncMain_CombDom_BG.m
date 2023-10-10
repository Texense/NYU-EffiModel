%% LDE iteration w/o figures
% L4EmeshX,L4ImeshY,LDEFrfunc: compsed by several different response
% functions, and we choose the best from them for every step
%% Problem for the BG parts: maybe not fine enough! 
% Let's involve the previously computed local responses...
%% Starting from the finest function!
% Output: Iteration firing maps, and L1 diff norm

% New feature 051823: 
function [LDERepFinal,LDEIL2Diff,L2DiffNormNeib,L2Diameter,LDEequv,FuncUseAll,NANFinal] = ...
    LDEIteration_25FuncMain_CombDom_BG(...
    PixInptCtgrUse,LDEIni,p,Epoc,...
    C_SS_mean,C_CS_mean,C_IS_mean,...
    C_SC_mean,C_CC_mean,C_IC_mean,...
    C_SI_mean,C_CI_mean,C_II_mean,...
    L4SEp, L4SIp, ...
    L4CEp, L4CIp, ...
    L4IEp, L4IIp, ...
    L4EmeshXAll,L4ImeshYAll,LDEFrfuncAll, ...
    L4EmeshXBG, L4ImeshYBG, LDEFrfuncBG, ...
    varargin)
NANFinal = false;

% specify NHCout
if ~isempty(varargin)
    N_HCOut = varargin{1};
    NPixX = varargin{2};
    NPixY = varargin{3};
else
    N_HCOut = 4; NPixX = 10; NPixY = 10;
end

if length(varargin)>3
    InhKillFlag = varargin{4};
else
    InhKillFlag = true;
end

if length(varargin)>4
   Outflag = varargin{5};
else
   Outflag = 'xn';
end

if length(varargin)>5
    EKp = varargin{6}; % excitation suppresion parameters
    IKp = varargin{7};
else
    EKp.Thrsld = 50; EKp.Highist = [200,150]; EKp.HardBound = 200; EKp.Slope = 0;
    
    IKp.Thrsld = 70; IKp.Highist = [100,95];  IKp.HardBound = 120; IKp.Slope = 0.95;
end

LDEItr = cell(Epoc+1,1);LDEItr{1} = LDEIni;
LDEEpoOut = cell(Epoc+1,1); LDEEpoOut{1} = LDEIni;
LDEoutVec = zeros(size(LDEIni.I,1)*3,Epoc+1);
LDEoutVec(:,1) = [LDEIni.S;LDEIni.C;LDEIni.I];
BGFlagall = false(3, Epoc+1);

% record the function quoted every step
FuncUseAll = zeros(Epoc,1);

%InhKillFlag = true; Thrsld = 70; Highist = [100,97]; % pars for inhibition suppresion
NANGlobalFlag = false;
for Epc = 1:Epoc
    if NANGlobalFlag
        continue
    end
    
    LDEInpt = LDEItr{Epc};
    LDEUse = LDEInpt;
    if InhKillFlag % apply inhibition suppresion
        LDEUse.E = LDEUse.S * (1-0.3077) + LDEUse.C * 0.3077;
        %ThrsldE = 50; HighistE = [200,150]; HardBoundE = 200; SlopeE = 0;
        LDEUseEAdj = InhKill(LDEUse.E, EKp.Thrsld, EKp.Highist, EKp.HardBound, EKp.Slope)./LDEUse.E;
        LDEUse.S = LDEUse.S .* LDEUseEAdj;
        LDEUse.C = LDEUse.C .* LDEUseEAdj;
        
        %ThrsldI = 70; HighistI = [100,95];  HardBoundI = 120; SlopeI = 0.95;
        LDEUse.I   = InhKill(LDEUse.I, IKp.Thrsld, IKp.Highist, IKp.HardBound, IKp.Slope);
        LDEUse.I   = InhKill(LDEUse.I, IKp.Thrsld, IKp.Highist, IKp.HardBound, IKp.Slope);%*0.8+LDEInpt.I*0.2 % compensate for the original recursive
    end
%     L4EUse = C_SS_mean*LDEUse.S + ...%- diag(C_SS_mean).*LDEUse.S + ...
%              C_CS_mean*LDEUse.S + ...%- diag(C_CS_mean).*LDEUse.S + ...
%              C_IS_mean*LDEUse.S + ...%- diag(C_IS_mean).*LDEUse.S + ...
%              C_SC_mean*LDEUse.C + ...%- diag(C_SC_mean).*LDEUse.C + ...
%              C_CC_mean*LDEUse.C + ...%- diag(C_CC_mean).*LDEUse.C + ...
%              C_IC_mean*LDEUse.C  ;%- diag(C_IC_mean).*LDEUse.C ;
%     L4IUse = C_SI_mean*LDEUse.I + ...%- diag(C_SI_mean).*LDEUse.I + ...
%              C_CI_mean*LDEUse.I + ...%- diag(C_CI_mean).*LDEUse.I + ...
%              C_II_mean*LDEUse.I  ;%- diag(C_II_mean).*LDEUse.I ;
    L4EUse_S = (C_SS_mean*LDEUse.S + C_SC_mean*LDEUse.C)/L4SEp;
    L4EUse_C = (C_CS_mean*LDEUse.S + C_CC_mean*LDEUse.C)/L4CEp;
    L4EUse_I = (C_IS_mean*LDEUse.S + C_IC_mean*LDEUse.C)/L4IEp;

    L4IUse_S = (C_SI_mean*LDEUse.I)/L4SIp;
    L4IUse_C = (C_CI_mean*LDEUse.I)/L4CIp;
    L4IUse_I = (C_II_mean*LDEUse.I)/L4IIp;
            
    % No Need to resymmetrize!!...
%     L4EUse = symmHCs(L4EUse,N_HCOut,NPixX,NPixY);
%     L4IUse = symmHCs(L4IUse,N_HCOut,NPixX,NPixY);
    
    % Choose the best option of domain: Large by defalt; if better move to
    % smalle and smaller
    FuncN = length(L4EmeshXAll); % should be consistent for all three vars
        
    % Refer to functions
    LDEOut = struct('S',[],'C',[],'I',[]);    
    % adjust for different functions
    nanFlag = true;
    FuncUse = 1; % using the finest by defalt, then go coarser if not good enough
    while nanFlag && FuncUse<=FuncN
        L4EmeshX = L4EmeshXAll{FuncUse};
        L4ImeshY = L4ImeshYAll{FuncUse};
        LDEFrfunc = LDEFrfuncAll{FuncUse};
        LDEOut.S = LDEIterFunc_Grating_25Func_BG(...
            L4EmeshX, L4ImeshY,...
            LDEFrfunc.S,...
            L4EmeshXBG, L4ImeshYBG,...
            LDEFrfuncBG.S,...
            L4EUse_S,L4IUse_S,...
            PixInptCtgrUse);
        %ORTPix,OBLPix,OPTPix,ORTOBLBd_Pix,OPTOBLBd_Pix);
        LDEOut.C = LDEIterFunc_Grating_25Func_BG(...
            L4EmeshX, L4ImeshY,...
            LDEFrfunc.C,...
            L4EmeshXBG, L4ImeshYBG,...
            LDEFrfuncBG.C,...
            L4EUse_C,L4IUse_C,...
            PixInptCtgrUse);
        % ORTPix,OBLPix,OPTPix,ORTOBLBd_Pix,OPTOBLBd_Pix);
        LDEOut.I = LDEIterFunc_Grating_25Func_BG(...
            L4EmeshX, L4ImeshY,...
            LDEFrfunc.I,...
            L4EmeshXBG, L4ImeshYBG,...
            LDEFrfuncBG.I,...
            L4EUse_I,L4IUse_I,...
            PixInptCtgrUse);
        %ORTPix,OBLPix,OPTPix,ORTOBLBd_Pix,OPTOBLBd_Pix);
        
        % Now check if output contains nan
        nanFlag = any(isnan([LDEOut.S,LDEOut.C,LDEOut.I]),'all');
        FuncUse = FuncUse+1;
    end
    
    
    % record the function used.
    FuncUseAll(Epc) = FuncUse;
    
    LDENext = struct(...
        'S',LDEOut.S*p + LDEInpt.S*(1-p),...
        'C',LDEOut.C*p + LDEInpt.C*(1-p),...
        'I',LDEOut.I*p + LDEInpt.I*(1-p));
    LDEItr{Epc+1} = LDENext; LDEEpoOut{Epc+1} = LDEOut;
    LDEoutVec(:,Epc+1) = [LDEOut.S; LDEOut.C; LDEOut.I];
    %BGFlagall(:,Epc+1) = [BGflagS;  BGflagC;  BGflagI];

    if FuncUse>FuncN && nanFlag % return if used up all functions but still getting nans
        fprintf('***Warning! NAN results in the %d epoch. Returning...\n',Epc)
        NANGlobalFlag = true;
    end
end

if strcmpi(Outflag,'f(xn)')
   %disp('Showing f(xn)')
   LDERepFinal = LDEoutVec;
elseif strcmpi(Outflag,'xn')
   %disp('Showing xn')
   LDERepFinal = LDEItr;
else
   disp('Unknown export flag, using xn')
   LDERepFinal = LDEItr;
end

%L2DiffNorm = zeros(Epoc,1);
L2DiffNormNeib = zeros(Epoc,1);
L2Diameter = zeros(Epoc,1);
LDEequv = [];
LDEIL2Diff = [];
% get an "equilibrium" of I, and assume >150 is fine
NANFinal = NANGlobalFlag;
if ~NANGlobalFlag
    LDEequv = mean(LDEoutVec(:,floor((Epc+1)*2/3):end),2);
    LDEIL2Diff = sqrt(sum((LDEoutVec - repmat(LDEequv,1,Epoc+1)).^2, 1));
    for Epc = 1:Epoc
        %     L2DiffNorm(Epc) = ...
        %         norm([LDEEpoOut{Epc}.S;LDEEpoOut{Epc}.C;LDEEpoOut{Epc}.I] - ...
        %         [LDEEpoOut{end}.S;LDEEpoOut{end}.C;LDEEpoOut{end}.I], 2);
        
        L2DiffNormNeib(Epc) = ...
            norm([LDEEpoOut{Epc}.S;LDEEpoOut{Epc}.C;LDEEpoOut{Epc}.I] - ...
            [LDEEpoOut{Epc+1}.S;LDEEpoOut{Epc+1}.C;LDEEpoOut{Epc+1}.I], 2);
        L2Diameter(Epc) = max(LDEIL2Diff(Epc:end));
    end
end
end