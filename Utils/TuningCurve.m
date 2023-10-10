% Compute tuning curve from only 4 angles



function [AnglePlot,TuningCurveTrunc,TunCurColor,PreferAng] = TuningCurve(...
    TruncList,AngleList,FrTunCurveRawDataAll,...
    x,y,AngleAll,NPixX,NPixY,cmap)

FrPixsAll.S = zeros(8,length(TruncList),length(AngleList));
FrPixsAll.C = zeros(8,length(TruncList),length(AngleList));
for AngCtgr = 1: length(AngleList)
    FrPixsAll.S(:,:,AngCtgr) = FrTunCurveRawDataAll{AngCtgr}.FrTunCurAllPix{y,x}.S;
    FrPixsAll.C(:,:,AngCtgr) = FrTunCurveRawDataAll{AngCtgr}.FrTunCurAllPix{y,x}.C;
end

% Get the tuning curves for different truncations
TuningCurveTrunc.S = zeros(8*length(AngleList),length(TruncList));
TuningCurveTrunc.C = zeros(8*length(AngleList),length(TruncList));
PixOrder = [1,2,7,4,5,6,3,8];% Group order of pixels
%PixOrder = 1:8;% Group order of pixels
for TruncInd = 1:length(TruncList)
    FrData1Trunc.S = squeeze(FrPixsAll.S(:,TruncInd,:));
    FrData1Trunc.C = squeeze(FrPixsAll.C(:,TruncInd,:));
    
    for PixUseInd = 1:8
        PixU = PixOrder(PixUseInd);
        if mod(PixUseInd,2) == 1
            TuningCurveTrunc.C(...
                (PixUseInd-1)*length(AngleList)+1:PixUseInd*length(AngleList),...
                TruncInd) = FrData1Trunc.C(PixU,:)';
            TuningCurveTrunc.S(...
                (PixUseInd-1)*length(AngleList)+1:PixUseInd*length(AngleList),...
                TruncInd) = FrData1Trunc.S(PixU,:)';
        else
            TuningCurveTrunc.C(...
                (PixUseInd-1)*length(AngleList)+1:PixUseInd*length(AngleList),...
                TruncInd) = FrData1Trunc.C(PixU,end:-1:1)';
            TuningCurveTrunc.S(...
                (PixUseInd-1)*length(AngleList)+1:PixUseInd*length(AngleList),...
                TruncInd) = FrData1Trunc.S(PixU,end:-1:1)';
        end
    end
end

AnglePlot = AngleAll;
for ODInd = 1:7
    AnglePlot = [AnglePlot,AngleAll+(ODInd)*22.5];
end
%% Theoretically the preferred angle
xC = NPixX*1.5+0.5; yC = NPixY*1.5+0.5;
xx = x+NPixX-xC; yy = y+NPixY-yC;
PreferAng = mod(angle(-xx + yy*1i),2*pi)/(2*pi);
ColorId = ceil(PreferAng*size(cmap,1));
TunCurColor = cmap(ColorId,:);


end