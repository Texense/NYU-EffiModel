% rotate (and flip ICs)
% varargin: Mirind: true for flip, false for not flip
% LDEequv contains SCI
function [LDEequvUse] = HCRot(LDEequv,rotID,N_HCOut,NPixX,NPixY, varargin)
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

if isstruct(LDEequv)
    
    fields = fieldnames(LDEequv);
    LDEequvUse = LDEequv;
    
    [HCX, HCY] = meshgrid(1:N_HCOut,1:N_HCOut);
    HCX = mod(HCX,2); HCY = mod(HCY,2);
    for FInd = 1:length(fields)
        CurrFieldVec = LDEequv.(fields{FInd});
        CurrFieldMap = reshape(CurrFieldVec,N_HCOut*NPixY,N_HCOut*NPixX);
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
        LDEequvUse.(fields{FInd})= reshape(MapOut,...
            length(CurrFieldVec),1);
    end
    
else
    [HCX, HCY] = meshgrid(1:N_HCOut,1:N_HCOut);
    HCX = mod(HCX,2); HCY = mod(HCY,2);
    CurrFieldVec = LDEequv;
    CurrFieldMap = reshape(CurrFieldVec,N_HCOut*NPixY,N_HCOut*NPixX);
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

end