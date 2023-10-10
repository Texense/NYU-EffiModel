%% This function generates arbitrary large field with lgn indexes, given a 3*3HC input
% Input:  lambda_SOn_drive: LGNon input for all S cells. There should be
%                           three unique values
%         temp: A template 1*N (same length of number of functions) that composes lambda_SOn_drive 
%               acsending, should be a row vector
%         NnSPixel: location of every S cell. Should be composed by field X
%                   and Y
%         N_HC,N_HCout: Number of HC for the whole map. n*n
%         NPixX,NPixY: Number of pixel per side
% Output: PixLGNCtgr: (N_HCout^2*NPixX*NPixY)-by-lgn matrix. lgn stand for
%                     different lgn types. Entries are percentage of each 
%                     type, sum of rows=1
%         LGNon: 1*lgn vecter, entries are lgn values of each type. from
%                 small to large

% A important variant: If InptS_drive is already ctgr type, then no need to
% assign anymore
function [PixLGNCtgr] = LGNIndSpat(InptS_drive, temp, NnSPixel,N_HC,N_HCout,NPixX,NPixY)
PixNum = NPixX*N_HC * NPixY*N_HC;
%% get LGN type of each cell
if size(InptS_drive,2) == 1
    InptTypeNeu = zeros(length(InptS_drive),length(temp));
    for NeuInd = 1:length(InptS_drive)
        CurrInpt = InptS_drive(NeuInd);
        if CurrInpt>temp(end) % larger than the maximum
            InptTypeNeu(NeuInd,end) = 1;
        elseif CurrInpt<temp(1) % smaller than the minimum
            InptTypeNeu(NeuInd,1) = 1;
        elseif ismember(CurrInpt,temp) % if exactly equals to certain template
            InptTypeNeu(NeuInd,:) = double(temp==CurrInpt);
        else % between two templates
            loc = find(sort([temp,CurrInpt])==CurrInpt);
            InptTypeNeu(NeuInd,loc-1) = (CurrInpt-temp(loc))/(temp(loc-1)-temp(loc));
            InptTypeNeu(NeuInd,loc)   = 1-InptTypeNeu(NeuInd,loc-1);
        end
    end
elseif size(InptS_drive,2) == length(temp)
    InptTypeNeu = InptS_drive;
else
    disp('illigal Inpt information. Check')
    return
end

InptTypePix = zeros(PixNum,length(temp));

for PixInd = 1:PixNum % for each pixel:
    % First get all neruons in this pixel
    CurrentNeu = InptTypeNeu(NnSPixel.Vec == PixInd,:);
    % get their lgn type and distribute in matrix
    %[a,~,c] = find(sum(ind2vec(CurrentNeu'),2)/length(CurrentNeu));
    InptTypePix(PixInd,:)  = mean(CurrentNeu,1);
end

%% From here doesn't need to change
% Reshape the lgn index vectors to maps
Ind_SOn_Map = zeros(NPixY*N_HC,NPixX*N_HC,length(temp));
for lgnTy = 1:length(temp)
    Ind_SOn_Map(:,:,lgnTy) = reshape(InptTypePix(:,lgnTy),NPixY*N_HC,NPixX*N_HC);
end

% Now extract a basic periodic 2*2HC map, then expand to a new N_HCout^2 map
Ind_SOn_Map2by2 = Ind_SOn_Map(1:2*NPixY , 1:2*NPixX, :);
% Do a symmetry here! First horizontal and vertical
Ind_SOn_Map2by2Sym = zeros(2*NPixY,2*NPixX,length(temp),4);
Ind_SOn_Map2by2Sym(:,:,:,1) = Ind_SOn_Map2by2;
Ind_SOn_Map2by2Sym(:,:,:,2) = Ind_SOn_Map2by2(end:-1:1,:,:);
Ind_SOn_Map2by2Sym(:,:,:,3) = Ind_SOn_Map2by2(:,end:-1:1,:);
Ind_SOn_Map2by2Sym(:,:,:,4) = Ind_SOn_Map2by2(end:-1:1,end:-1:1,:);
Ind_SOn_Map2by2 = mean(Ind_SOn_Map2by2Sym,4);
% Do a fold symmetry here WITHIN 2by2
Ind_SOn_Map1by1 = Ind_SOn_Map2by2(1:NPixY,1:NPixX,:);
Ind_SOn_Map2by2(1:NPixY,NPixX+1:end,:)     = Ind_SOn_Map1by1(:,end:-1:1,:);
Ind_SOn_Map2by2(NPixY+1:end,1:NPixX,:)     = Ind_SOn_Map1by1(end:-1:1,:,:);
Ind_SOn_Map2by2(NPixY+1:end,NPixX+1:end,:) = Ind_SOn_Map1by1(end:-1:1,end:-1:1,:);
% Ind_SOn_Map2by2 = mean(Ind_SOn_Map2by2Dia,4);

PixIndX = mod(1:N_HCout*NPixX,2*NPixX); PixIndX(PixIndX==0) = 2*NPixX;
PixIndY = mod(1:N_HCout*NPixY,2*NPixY); PixIndY(PixIndY==0) = 2*NPixY;
Ind_SOn_MapOut = Ind_SOn_Map2by2(PixIndY,PixIndX,:);

PixLGNCtgr = zeros(NPixY*N_HCout*NPixX*N_HCout,length(temp));
for lgnTy = 1:length(temp)
    PixLGNCtgr(:,lgnTy) = reshape(Ind_SOn_MapOut(:,:,lgnTy),NPixY*N_HCout*NPixX*N_HCout,1);
end
end