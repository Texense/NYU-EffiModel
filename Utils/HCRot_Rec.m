% rotate (and flip ICs)
% varargin: Mirind: true for flip, false for not flip
% Abanden the move for every field: just do one field at a time.
function [LDEequvUse] = HCRot_Rec(LDEequv,rotID,N_HCOutX,N_HCOutY,NPixX,NPixY, varargin)
if ~isempty(varargin)
    Mirind = varargin{1};
else
    Mirind = false;
end

if Mirind
    flipDim = mod(rotID,2); %1 for col 2 ofr row
    if flipDim == 0
        flipDim = 2;
    end
    rotID = mod(rotID + 1, 4);% 0-3
end


    [HCX, HCY] = meshgrid(1:N_HCOutX,1:N_HCOutY);
    HCX = mod(HCX,2); HCY = mod(HCY,2);
    CurrFieldVec = LDEequv;
    CurrFieldMap = reshape(CurrFieldVec,N_HCOutY*NPixY,N_HCOutX*NPixX);
    OneHC = rot90(CurrFieldMap(1:NPixY,1:NPixX),rotID); % get rotated one HC
    % OneHC: flip or not?
    if Mirind
        OneHC = flip(OneHC,flipDim);
    end
    
    MapOut = zeros(size(CurrFieldMap));
    for xid = 1:N_HCOut
        for yid = 1:N_HCOut
            xmod = HCX(yid, xid); ymod = HCY(yid, xid);
            HCHold = OneHC;
            if xmod == 0
                HCHold = HCHold(:,end:-1:1);
            end
            
            if ymod == 0
                HCHold = HCHold(end:-1:1,:);
            end
            MapOut((yid-1)*NPixY+1:yid*NPixY, (xid-1)*NPixX+1:xid*NPixX) = HCHold;
        end
    end
    LDEequvUse= reshape(MapOut,...
        length(CurrFieldVec),1);

end