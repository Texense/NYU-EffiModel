% The goal of this functions is to take an input firing rate map with
% arbitrary rectangular shape, then take the first HC, then output another
% firing rate map with arbitrary rectangular shape based on rotating the HC 
% periodicly. 
%
% The used HC can be rotated or flipped to implement ICs based on NESSs of 
% different angles.
% varargin: 
%           Mirind: true for flip, false for not flip
% LDEequv can contains SCI (3 fields of maps) or just one field of map
%% New version: can tolerate nonsquare patch in output now
function [LDEequvUse] = HCRotXY(LDEequv,rotID,N_HCinX,N_HCinY,N_HCOutX,N_HCOutY,NPixX,NPixY, varargin)
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
    
    [HCX, HCY] = meshgrid(1:N_HCOutX,1:N_HCOutY);
    HCX = mod(HCX,2); HCY = mod(HCY,2);
    for FInd = 1:length(fields)
        CurrFieldVec = LDEequv.(fields{FInd});
        CurrFieldMap = reshape(CurrFieldVec,N_HCinY*NPixY,N_HCinX*NPixX);
        OneHC = rot90(CurrFieldMap(1:NPixY,1:NPixX),rotID); % get rotated one HC
        % OneHC: flip or not?
        if Mirind
            OneHC = flip(OneHC,flipDim);
        end
        
        MapOut = zeros(N_HCOutY*NPixY,N_HCOutX*NPixX);
        for xid = 1:N_HCOutX
            for yid = 1:N_HCOutY
                xmod = HCX(yid, xid); ymod = HCY(yid, xid);
                HCHold = OneHC;
                if xmod == 0
                    HCHold = HCHold(:,end:-1:1);
                end
                
                if ymod == 0
                    HCHold = HCHold(end:-1:1,:);
                end
                MapOut((yid-1)*NPixY+1:yid*NPixY, (xid-1)*NPixX+1:xid*NPixX)...
                    = HCHold;
            end
        end
        LDEequvUse.(fields{FInd})= reshape(MapOut,...
            N_HCOutY*NPixY*N_HCOutX*NPixX,1);
    end
    
else
    [HCX, HCY] = meshgrid(1:N_HCOutX,1:N_HCOutY);
    HCX = mod(HCX,2); HCY = mod(HCY,2);
    CurrFieldVec = LDEequv;
    CurrFieldMap = reshape(CurrFieldVec,N_HCinY*NPixY,N_HCinX*NPixX);
    OneHC = rot90(CurrFieldMap(1:NPixY,1:NPixX),rotID); % get rotated one HC
    % OneHC: flip or not?
    if Mirind
        OneHC = flip(OneHC,flipDim);
    end
    
    MapOut = zeros(N_HCOutY*NPixY,N_HCOutX*NPixX);
    for xid = 1:N_HCOutX
        for yid = 1:N_HCOutY
            xmod = HCX(yid, xid); ymod = HCY(yid, xid);
            HCHold = OneHC;
            if xmod == 0
                HCHold = HCHold(:,end:-1:1);
            end
            
            if ymod == 0
                HCHold = HCHold(end:-1:1,:);
            end
            MapOut((yid-1)*NPixY+1:yid*NPixY, (xid-1)*NPixX+1:xid*NPixX)...
                = HCHold;
        end
    end
    LDEequvUse= reshape(MapOut,...
        N_HCOutY*NPixY*N_HCOutX*NPixX,1);
end

end