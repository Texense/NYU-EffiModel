%% LDE iteration w/o figures
% Output: Iteration firing maps, and L1 diff norm
% Outflag: 'f(xn)' or 'xn'
function [LDERepFinal,L2DiffNorm,L2DiffNormNeib,L2Diameter,LDEequv] = ...
    LDEIteration_16FuncMain(...
    PixInptCtgrUse,LDEIni,p,Epoc,...
    C_SS_mean,C_CS_mean,C_IS_mean,...
    C_SC_mean,C_CC_mean,C_IC_mean,...
    C_SI_mean,C_CI_mean,C_II_mean,...
    L4EmeshX,L4ImeshY,LDEFrfunc,varargin)
% specify NHCout
if length(varargin)>0
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

LDEItr = cell(Epoc+1,1);LDEItr{1} = LDEIni;
LDEEpoOut = cell(Epoc+1,1); LDEEpoOut{1} = LDEIni;
LDEoutVec = zeros(size(LDEIni.I,1)*3,Epoc+1);
LDEoutVec(:,1) = [LDEIni.S;LDEIni.C;LDEIni.I];


%InhKillFlag = true;  % pars for inhibition suppresion
for Epc = 1:Epoc
    % First do convolution
    LDEInpt = LDEItr{Epc};
    LDEUse = LDEInpt;

    if InhKillFlag % apply inhibition suppresion
        Thrsld = 70; Highist = [100,95];
        LDEUse.I = InhKill(LDEUse.I, Thrsld, Highist);
        LDEUse.I = InhKill(LDEUse.I, Thrsld, Highist);
    end
    L4EUse = C_SS_mean*LDEUse.S + ...%- diag(C_SS_mean).*LDEUse.S + ...
        C_CS_mean*LDEUse.S + ...%- diag(C_CS_mean).*LDEUse.S + ...
        C_IS_mean*LDEUse.S + ...%- diag(C_IS_mean).*LDEUse.S + ...
        C_SC_mean*LDEUse.C + ...%- diag(C_SC_mean).*LDEUse.C + ...
        C_CC_mean*LDEUse.C + ...%- diag(C_CC_mean).*LDEUse.C + ...
        C_IC_mean*LDEUse.C  ;%- diag(C_IC_mean).*LDEUse.C ;
    L4IUse = C_SI_mean*LDEUse.I + ...%- diag(C_SI_mean).*LDEUse.I + ...
        C_CI_mean*LDEUse.I + ...%- diag(C_CI_mean).*LDEUse.I + ...
        C_II_mean*LDEUse.I  ;%- diag(C_II_mean).*LDEUse.I ;
    
    % Need to resymmetrize!!...
    L4EUse = symmHCs(L4EUse,N_HCOut,NPixX,NPixY);
    L4IUse = symmHCs(L4IUse,N_HCOut,NPixX,NPixY);
    
    % Refer to functions
    LDEOut = struct('S',[],'C',[],'I',[]);
    LDEOut.S = LDEIterFunc_Grating_16Func(...
        L4EmeshX, L4ImeshY,...
        LDEFrfunc.S,...
        L4EUse,L4IUse,...
        PixInptCtgrUse);
    %ORTPix,OBLPix,OPTPix,ORTOBLBd_Pix,OPTOBLBd_Pix);
    LDEOut.C = LDEIterFunc_Grating_16Func(...
        L4EmeshX, L4ImeshY,...
        LDEFrfunc.C,...
        L4EUse,L4IUse,...
        PixInptCtgrUse);
    % ORTPix,OBLPix,OPTPix,ORTOBLBd_Pix,OPTOBLBd_Pix);
    LDEOut.I = LDEIterFunc_Grating_16Func(...
        L4EmeshX, L4ImeshY,...
        LDEFrfunc.I,...
        L4EUse,L4IUse,...
        PixInptCtgrUse);
    %ORTPix,OBLPix,OPTPix,ORTOBLBd_Pix,OPTOBLBd_Pix);
    LDENext = struct(...
        'S',LDEOut.S*p + LDEInpt.S*(1-p),...
        'C',LDEOut.C*p + LDEInpt.C*(1-p),...
        'I',LDEOut.I*p + LDEInpt.I*(1-p));
    LDEItr{Epc+1} = LDENext; LDEEpoOut{Epc+1} = LDEOut;
    LDEoutVec(:,Epc+1) = [LDEOut.S; LDEOut.C; LDEOut.I];
end

if strcmpi(Outflag,'f(xn)')
   disp('Showing f(xn)')
   LDERepFinal = LDEEpoOut;
elseif strcmpi(Outflag,'xn')
   disp('Showing xn')
   LDERepFinal = LDEItr;
else
   disp('Unknown export flag, using xn')
   LDERepFinal = LDEItr;
end

L2DiffNorm = zeros(Epoc,1);
L2DiffNormNeib = zeros(Epoc,1);
L2Diameter = zeros(Epoc,1);
% get an "equilibrium" of I, and assume >150 is fine
LDEequv = mean(LDEoutVec(:,floor((Epc+1)*2/3):end),2);
LDEIL2Diff = sqrt(sum((LDEoutVec - repmat(LDEequv,1,Epoc+1)).^2, 1));
for Epc = 1:Epoc
    L2DiffNorm(Epc) = ...
        norm([LDEEpoOut{Epc}.S;LDEEpoOut{Epc}.C;LDEEpoOut{Epc}.I] - ...
        [LDEEpoOut{end}.S;LDEEpoOut{end}.C;LDEEpoOut{end}.I], 2);
    
    L2DiffNormNeib(Epc) = ...
        norm([LDEEpoOut{Epc}.S;LDEEpoOut{Epc}.C;LDEEpoOut{Epc}.I] - ...
        [LDEEpoOut{Epc+1}.S;LDEEpoOut{Epc+1}.C;LDEEpoOut{Epc+1}.I], 2);
    L2Diameter(Epc) = max(LDEIL2Diff(Epc:end));
end
end