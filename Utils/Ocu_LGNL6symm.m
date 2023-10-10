% Has to find a way to symmetrize the LGN and L6 inputs to pixels. The
% insymmetry largely comes from the inhomogeneous distributions of cells in
% different pixels.

% PixL6Ctgr: N*5 (or other number) matrix. Each column represents one
% domain (orientation * ocular).

%% NOTE:!! The code here largely depending on that N_HCs are even numbers, and the periodic pattern is 2*2.

function PixL6CtgrOut = Ocu_LGNL6symm(PixL6Ctgr,N_HCOutX,N_HCOutY,NPixX,NPixY)
% every domain admits a horizontal mirror symmetry
PixL6CtgrOut = zeros(size(PixL6Ctgr));
PixCtgr24Hold = zeros(N_HCOutY*NPixY,N_HCOutX*NPixX,2);
HoldId = 1;
for ii = 1:size(PixL6Ctgr,2)
    PixL6Now = reshape(PixL6Ctgr(:,ii),N_HCOutY*NPixY,N_HCOutX*NPixX);
    PixL6NowHorzSymm = PixL6Now(:,end:-1:1);
    PixL6Out = (PixL6Now+PixL6NowHorzSymm)/2;
    % for domain 1, 3, 5: vertical mirror symmetry (mode every two rows) as
    % well
    if mod(ii,2) == 1 
        PixL6OutVerSymm = PixL6Out(end:-1:1,:);
        PixL6OutVerSymm = PixL6OutVerSymm([NPixY+1:end,1:NPixY],:); % adjust for the row deviation
        PixL6Out = (PixL6Out + PixL6OutVerSymm)/2;
    end
    PixL6CtgrOut(:,ii) = reshape(PixL6Out,N_HCOutY*NPixY*N_HCOutX*NPixX,1);
    if (ii == 2 || ii == 4)
        PixCtgr24Hold(:,:,HoldId) = PixL6Out;
        HoldId = HoldId + 1;
    end
end
PixCtgr24symm = squeeze(PixCtgr24Hold(end:-1:1,:,2));
PixCtgr24symm = PixCtgr24symm([NPixY+1:end,1:NPixY],:);
PixCtgr2ave = (squeeze(PixCtgr24Hold(:,:,1)) + PixCtgr24symm)/2;

PixL6CtgrOut(:,2) = reshape(PixCtgr2ave,N_HCOutY*NPixY*N_HCOutX*NPixX,1);

%finally do Ctgry 4;
PixCtgr4ave = PixCtgr2ave(end:-1:1,:);
PixCtgr4ave = PixCtgr4ave([NPixY+1:end,1:NPixY],:);
PixL6CtgrOut(:,4) = reshape(PixCtgr4ave,N_HCOutY*NPixY*N_HCOutX*NPixX,1);
end