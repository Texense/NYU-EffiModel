% Construct a symmetrical input map!

function [PixInptCtgr,PixInptCtgrUse]  =  SymInptCtgr...
    (N_HC,N_HCOut,n_S_HCsur,n_I_HC,NPixX,NPixY,ODNum,...
    L6smear,LGNsmear,Trunc,FigOn,PixNumOut,LGNlist,L6list)


[NnSsur.X,NnSsur.Y] = ...
    V1Field_Generation(N_HC,1:n_S_HCsur^2*N_HC^2,'e',n_S_HCsur,n_I_HC);
NnSPixelsur.X = ceil(NnSsur.X/(max(NnSsur.X)/(NPixX*N_HC)));
NnSPixelsur.Y = ceil(NnSsur.Y/(max(NnSsur.Y)/(NPixY*N_HC)));
NnSPixelsur.Vec = NnSPixelsur.Y + (NnSPixelsur.X-1)*(NPixY*N_HC);

OD_SMapsur = zeros(n_S_HCsur*N_HC,'single');
for NeuInd = 1:length(NnSsur.X)
    OD_SMapsur(NnSsur.Y(NeuInd),NnSsur.X(NeuInd)) = ...
        OrientDom(ODNum,NnSsur.X(NeuInd),NnSsur.Y(NeuInd),n_S_HCsur);
end


L6SFilt_Grating = SymDiag(OD_SMapsur,...
    N_HC,n_S_HCsur,L6smear,Trunc,FigOn);

LGNFilt_Grating = SymDiag(OD_SMapsur,...
    N_HC,n_S_HCsur,LGNsmear,Trunc,FigOn);

PixL6Ctgr = LGNIndSpat(L6SFilt_Grating,1:4,NnSPixelsur,N_HC,N_HCOut,NPixX,NPixY);
PixLGNCtgr = LGNIndSpat(LGNFilt_Grating,1:4,NnSPixelsur,N_HC,N_HCOut,NPixX,NPixY);
% Just note: This doesn't mean anything now...
PixInptCtgr = PixL6Ctgr*0.7 + PixLGNCtgr*0.3;
PixInptCtgrUse  = zeros(PixNumOut,LGNlist,L6list);
for PixInd = 1:PixNumOut
    PixInptCtgrUse(PixInd,:,:) = PixLGNCtgr(PixInd,:)' *PixL6Ctgr(PixInd,:);% by multiplying both indexes
end

end

% RotId: 1-4, 1: same; 2: diagonal; 3: antidiag; 4: both
function L6SFilt_Grating = SymDiag(OD_SMapsur,...
    N_HC,n_S_HCsur,L6smear,Trunc,FigOn)

L6SFilt_GratingSym(:,:,1) = SpatialGaussianFilt_my(... % original
    ODSwap(OD_SMapsur,[1,2,3,4;1,2,3,4]),...
    N_HC,n_S_HCsur,L6smear,Trunc,FigOn);
L6SFilt_GratingSym(:,:,2) = SpatialGaussianFilt_my(... % diagonal
    ODSwap(OD_SMapsur',[1,2,3,4;2,1,4,3]),...
    N_HC,n_S_HCsur,L6smear,Trunc,FigOn);
L6SFilt_GratingSym(:,:,3) = SpatialGaussianFilt_my(... % antidiagonal
    ODSwap(rot90(OD_SMapsur,2)',[1,2,3,4;4,3,2,1]),...
    N_HC,n_S_HCsur,L6smear,Trunc,FigOn);
L6SFilt_GratingSym(:,:,4) = SpatialGaussianFilt_my(... % both
    ODSwap(rot90(OD_SMapsur,2), [1,2,3,4;3,4,1,2]),...
    N_HC,n_S_HCsur,L6smear,Trunc,FigOn);

L6SFilt_Grating = mean(L6SFilt_GratingSym,3);
end

% Mapping: [1,2,3,4;a,b,c,d]
function OD_SMapSwap = ODSwap(OD_SMapsur,Mapping)
OD_SMapSwap = zeros(size(OD_SMapsur));
for MapInd = 1:size(Mapping,2)
    OD_SMapSwap(OD_SMapsur==Mapping(1,MapInd)) = Mapping(2,MapInd);
end
end