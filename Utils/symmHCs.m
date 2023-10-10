% Symmetrize 
% Mode 1: Symmetrize a Field with some noise
% Mode 2: add symmetric perturbations

% PertMode: 'add' add the Pert to one HC
%           'prod' multiply to one HC
function L4EUseOut = symmHCs(L4EUse,N_HCOut,NPixX,NPixY,varargin)
if isempty(varargin)
    FlagPerturb = false;
else
    FlagPerturb = true;
    Pert = varargin{1};
    if length(varargin)>1
        PertMode = varargin{2};
    else
        PertMode = 'add'; % additive noise by defalut
    end
end


L4EMap = reshape(L4EUse,N_HCOut*NPixY,N_HCOut*NPixX);
[HCX, HCY] = meshgrid(1:N_HCOut,1:N_HCOut);
HCX = mod(HCX,2); HCY = mod(HCY,2); 

% pack all HCs together then get symmetry
L4EUseOutAllHC = zeros(NPixY,NPixX,N_HCOut^2);
HCid = 1;
for xid = 1:N_HCOut
    for yid = 1:N_HCOut
        xmod = HCX(yid, xid); ymod = HCY(yid, xid); 
        HCOrig = L4EMap((yid-1)*NPixY+1:yid*NPixY, (xid-1)*NPixX+1:xid*NPixX);
        if xmod == 0
            HCOrig = HCOrig(:,end:-1:1);
        end
        
        if ymod == 0
            HCOrig = HCOrig(end:-1:1,:);
        end
        L4EUseOutAllHC(:,:,HCid) = HCOrig;
        HCid = HCid +1;
    end    
end
L4EUseOutOneHC = mean(L4EUseOutAllHC,3);

% add the perturbation
if FlagPerturb
    if strcmpi(PertMode,'prod')
        L4EUseOutOneHC = L4EUseOutOneHC .* reshape(Pert,NPixY,NPixX);
    elseif strcmpi(PertMode,'add')
    L4EUseOutOneHC = L4EUseOutOneHC + reshape(Pert,NPixY,NPixX);    
    else
        diap('***Illigal pert mode. Quiting...')
    end
end

% Map the HC back
L4EUseOut = zeros(size(L4EMap));
for xid = 1:N_HCOut
    for yid = 1:N_HCOut
        xmod = HCX(yid, xid); ymod = HCY(yid, xid); 
        HCHold = L4EUseOutOneHC;
        if xmod == 0
            HCHold = HCHold(:,end:-1:1);
        end
        
        if ymod == 0
            HCHold = HCHold(end:-1:1,:);
        end
        L4EUseOut((yid-1)*NPixY+1:yid*NPixY, (xid-1)*NPixX+1:xid*NPixX) = HCHold;
    end
end

L4EUseOut = reshape(L4EUseOut,N_HCOut^2*NPixX*NPixY,1);
end