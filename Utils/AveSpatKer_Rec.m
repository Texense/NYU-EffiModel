%% function C_SS_mean = AveSpatKer_Rec(C_SS_Pixel_Us,N_HC,NPixX,NPixY)
% Input: C_SS_Pixel_Us: A Pixel connectivity matrix (from network).
% Pre(col)*Post(row)
%        N_HC,NPixX,NPixY: # of hypercolumns, # of pixels per HC
% Output: C_SS_mean: Averaged spatial kernal of firing rats. Same size as
% C_SS_Pixel_Us

% Zhuo-Cheng Xiao 04/21/2023
% New Version: Want to expand from 3HC to arbitrary # of HC
% N_HCout should be larger than 3

%% NOTE: I should implement the ocular modulation of sptial kernels in this function
% Hence it should be further modified. Now I am just adapting for 4*8 HC
% outputs...
% Parax, Paray stands for the two parameters for ocular modulation:
%       Parax: how much is killed for cross ocular column connections. 
%               0: all killed; 1: all preserved
%       Paray: how much of the killed connections is returned to the mirror 
%       pixel in the same ocular column
%               0: no returned; 1: all returned

function C_PQ_mean = AveSpatKer_Rec(...
    C_PQ_Pixel_Us,N_HCin,N_HCoutX,N_HCoutY,NPixX,NPixY, ...
    varargin)
%% Check if C_SS is an eligible matrix
if isempty(varargin)
    Parax = 1; Paray = 0;
else
    Parax = varargin{1}; 
    Paray = varargin{2};
end

if size(C_PQ_Pixel_Us,1) ~= size(C_PQ_Pixel_Us,2)
    error('Input Matrix not square!')
end
if size(C_PQ_Pixel_Us,2) ~= N_HCin*NPixX*N_HCin*NPixY
    error('# of pixels doesnt match!')
end

if (N_HCin<3)
    disp('Warning!: Reverse averaged kernal may have duplicate pixel infos')  
end
% if N_HCout<3
%     disp('Not enough HC for output. Return...') 
%     return
% end
%% For each pixel: Extract a density map around it. Radius at most one HC,
% so 2HC-by-2HC
PrySynDist_all = zeros(2*NPixY,2*NPixX,size(C_PQ_Pixel_Us,1));
for PixInd = 1:size(C_PQ_Pixel_Us,1)
    % first get the spatial coord of the pixel
    PX = ceil(PixInd/(N_HCin*NPixY)); 
    PY = mod(PixInd,N_HCin*NPixY);
    if PY == 0
        PY = N_HCin*NPixY;
    end
    
    % Make an extended map
    %% We are having big problem here. Actually this only works for 3*3 inputs 
    % -- let's modify later for other inputs
    ConnVec = C_PQ_Pixel_Us(PixInd,:)';
    ConnMap = reshape(ConnVec,N_HCin*NPixY,N_HCin*NPixX);
    Cen = floor((min(N_HCin,N_HCin)-1)/2);
    HC_Center = ConnMap(Cen*NPixY+1:(Cen+1)*NPixY,Cen*NPixX+1:(Cen+1)*NPixX);
    Bar_Center_hor = ConnMap(Cen*NPixY+1:(Cen+1)*NPixY,:);
    Bar_Center_ver = ConnMap(:,Cen*NPixX+1:(Cen+1)*NPixX);
    
    ConnMap_Ext = zeros((N_HCin+2)*NPixY,(N_HCin+2)*NPixX); %extend each side by 1 HC
    ConnMap_Ext(NPixY+1:(N_HCin+1)*NPixY,...
                NPixX+1:(N_HCin+1)*NPixX) = ConnMap;
    ConnMap_Ext(1:NPixY,...
                1:NPixX) = HC_Center; % 4 corners
    ConnMap_Ext((N_HCin+1)*NPixY+1:(N_HCin+2)*NPixY,...
                1:NPixX) = HC_Center;
    ConnMap_Ext(1:NPixY,...
                (N_HCin+1)*NPixX+1:(N_HCin+2)*NPixX) = HC_Center; 
    ConnMap_Ext((N_HCin+1)*NPixY+1:(N_HCin+2)*NPixY,...
                (N_HCin+1)*NPixX+1:(N_HCin+2)*NPixX) = HC_Center;
    ConnMap_Ext(1:NPixY,...
                NPixX+1:(N_HCin+1)*NPixX) = Bar_Center_hor; % 4 sides
    ConnMap_Ext((N_HCin+1)*NPixY+1:(N_HCin+2)*NPixY,...
                NPixX+1:(N_HCin+1)*NPixX) = Bar_Center_hor; 
    ConnMap_Ext(NPixY+1:(N_HCin+1)*NPixY,...
                1:NPixX) = Bar_Center_ver; 
    ConnMap_Ext(NPixY+1:(N_HCin+1)*NPixY,...
                (N_HCin+1)*NPixX+1:(N_HCin+2)*NPixX) = Bar_Center_ver; 
    
    % Extract a 2HC-by-2HC Presynaptic distribution
    PX_New = PX+NPixX; PY_New = PY+NPixY;
    PrySynDist_all(:,:,PixInd) = ConnMap_Ext(PY_New-NPixY:PY_New+NPixY-1,...
                                             PX_New-NPixX:PX_New+NPixX-1);
end
Ker_PQ_mean = mean(PrySynDist_all,3);
%% We then symmetricalize. The center should always be NPixY+1, NPixX+1
Ker_PQ_meanSym = zeros(2*NPixY,2*NPixX,2);
Ker_PQ_meanSym(2:end,2:end,1) = Ker_PQ_mean(2:end,2:end);
Ker_PQ_meanSym(2:end,2:end,2) = Ker_PQ_mean(end:-1:2,2:end);
Ker_PQ_mean = mean(Ker_PQ_meanSym,3);
Ker_PQ_meanSym = zeros(2*NPixY,2*NPixX,2);
Ker_PQ_meanSym(2:end,2:end,1) = Ker_PQ_mean(2:end,end:-1:2);
Ker_PQ_meanSym(2:end,2:end,2) = Ker_PQ_mean(2:end,2:end);
Ker_PQ_mean = mean(Ker_PQ_meanSym,3);
% Do a symmetry here! Then diagonal and antidiagonal
Ker_PQ_meanDia = zeros(2*NPixY,2*NPixX,2);
Ker_PQ_meanDia(2:end,2:end,1) = Ker_PQ_mean(2:end,2:end);
Ker_PQ_meanDia(2:end,2:end,2) = Ker_PQ_mean(2:end,2:end)';
Ker_PQ_mean = mean(Ker_PQ_meanDia,3);

Ker_PQ_meanDia(2:end,2:end,3) = rot90(Ker_PQ_mean(2:end,2:end),2);
Ker_PQ_meanDia(2:end,2:end,4) = Ker_PQ_mean(2:end,2:end);
Ker_PQ_mean = mean(Ker_PQ_meanDia,3);

%% Put Kernel back to Matrix
C_PQ_mean  = zeros(N_HCoutX*NPixX*N_HCoutY*NPixY);

for PixInd = 1:N_HCoutX*NPixX*N_HCoutY*NPixY
    PX = ceil(PixInd/(N_HCoutY*NPixY));
    PY = mod(PixInd,N_HCoutY*NPixY);
    if PY == 0
        PY = N_HCoutY*NPixY;
    end
    
    %% Below for odd N_HC is also problematic -- leave for now as well.
    % If I really want to implement for odd HCs, boundary conditions &
    % different # rows/cols will be problems
    if mod(N_HCoutY,2) == 1
        N_HCout = N_HCoutX; %% REDUANDENT! This is only for not making errors...
        % If odd: Note, this will only work for 4n-1. 4n+1 will introduce error on
        % bdry
        % Create an extended map and center the smoothed kernel around a pixel
        PX_New = PX+NPixX; PY_New = PY+NPixY;
        ConnMap_Rev = zeros((N_HCout+2)*NPixY,(N_HCout+2)*NPixX);
        ConnMap_Rev(PY_New-NPixY:PY_New+NPixY-1,...
                    PX_New-NPixX:PX_New+NPixX-1) = Ker_PQ_mean;
        %Add cancidates up, since the 2HC-by-2HC mat can't extend to three consecutive HCs...
        CenExt = floor((N_HCout+2-1)/2);
        %HC_Cen_Ext = ConnMap_Rev(CenExt*NPixY+1:(CenExt+1)*NPixY, CenExt*NPixX+1:(CenExt+1)*NPixX);        
        % Fisrt Pile back four sides (get rid of corners)
              ConnMap_Rev(     CenExt*NPixY+1: (CenExt+1)*NPixY, NPixX+1:(N_HCout+1)*NPixX) ...% Center Row...
            = ConnMap_Rev(     CenExt*NPixY+1: (CenExt+1)*NPixY, NPixX+1:(N_HCout+1)*NPixX) ...
            + ConnMap_Rev(                  1:            NPixY, NPixX+1:(N_HCout+1)*NPixX) ...
            + ConnMap_Rev((N_HCout+1)*NPixY+1:(N_HCout+2)*NPixY, NPixX+1:(N_HCout+1)*NPixX);
        
              ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY,      CenExt*NPixX+1: (CenExt+1)*NPixX) ...% center column,
            = ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY,      CenExt*NPixX+1: (CenExt+1)*NPixX) ...
            + ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY,                   1:            NPixX) ...
            + ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY, (N_HCout+1)*NPixX+1:(N_HCout+2)*NPixX);
        
              ConnMap_Rev(     CenExt*NPixY+1: (CenExt+1)*NPixY,      CenExt*NPixX+1: (CenExt+1)*NPixX) ... % center, 4 corners, 4 midHC come together
            = ConnMap_Rev(     CenExt*NPixY+1: (CenExt+1)*NPixY,      CenExt*NPixX+1: (CenExt+1)*NPixX) ...
            + ConnMap_Rev(                  1:            NPixY,                   1:            NPixX) ...
            + ConnMap_Rev((N_HCout+1)*NPixY+1:(N_HCout+2)*NPixY,                   1:            NPixX) ...
            + ConnMap_Rev(                  1:            NPixY, (N_HCout+1)*NPixX+1:(N_HCout+2)*NPixX) ...
            + ConnMap_Rev((N_HCout+1)*NPixY+1:(N_HCout+2)*NPixY, (N_HCout+1)*NPixX+1:(N_HCout+2)*NPixX);
        ConnMapPix = ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY, NPixX+1:(N_HCout+1)*NPixX);   
    % We are actually only using the even case  % If even: Periodic boundary conditions  
    else 
        ConnMap_Rev = zeros(N_HCoutY*NPixY,N_HCoutX*NPixX);
        % Periodic boundarys for projections
        YRange = mod(PY-NPixY:PY+NPixY-1,N_HCoutY*NPixY); YRange(YRange==0) = N_HCoutY*NPixY;
        XRange = mod(PX-NPixX:PX+NPixX-1,N_HCoutX*NPixX); XRange(XRange==0) = N_HCoutX*NPixX;
        ConnMap_Rev(YRange,XRange) = Ker_PQ_mean;
        
        %% Now I want to implement ocular modulation to connectivity.
        % The big problema in the last round: I forgot that some kernals
        % could extend to as many as 3 rows! Need to mirror back both ends!
        OcuIndi = ones(size(ConnMap_Rev)); % ocular column indicator
        OcuIndi(mod(ceil((1:N_HCoutY*NPixY)/NPixY),2)==1,:) = 0;
        % first, determine the ocular column of center of the kernel
        CenterOcuFlag = OcuIndi(PY,PX);
        % Then determine the same-ocu and diff-ocu connections
        ConnMap_SameOcu = ConnMap_Rev; ConnMap_DiffOcu = ConnMap_Rev; ConnMap_MirrOcu = zeros(size(ConnMap_Rev));
        ConnMap_SameOcu(OcuIndi~=CenterOcuFlag) = 0;
        ConnMap_DiffOcu(OcuIndi==CenterOcuFlag) = 0;
        % Then kill the diff-ocu connections using Parax
        ConnMap_DiffOcuU = ConnMap_DiffOcu*Parax;
        % Then, map the killed part to the same-ocu            
            % Dscarded code: [DiffRowInd,~] = find(ConnMap_DiffOcu); % First, which row is the DiffKernel?
            %               DiffRow_HC = mode(ceil(DiffRowInd/NPixY)); %% NOTE!!: 
            % here I use the most frequent row -- but this will FAIL if the kernel is too large and covers multiple rows...
        %% Basic assumption: any L4 kernel does not extend more than 3 rows!
        Mirror_HC = ceil(PY/NPixY);             % Then, which row does the center lie in?
        DiffRow_HC = mod([Mirror_HC+1,Mirror_HC-1],N_HCoutY);
        DiffRow_HC(DiffRow_HC == 0) = N_HCoutY;
        ConnMap_MirrOcu(Mirror_HC*NPixY:-1:(Mirror_HC-1)*NPixY+1,:) = ...
            ConnMap_DiffOcu((DiffRow_HC(1)-1)*NPixY+1:DiffRow_HC(1)*NPixY,:) + ...
            ConnMap_DiffOcu((DiffRow_HC(2)-1)*NPixY+1:DiffRow_HC(2)*NPixY,:); % map to the mirror
        ConnMap_MirrOcuU = Paray*(1-Parax)*ConnMap_MirrOcu;
        % lastly, put up everthing together
        ConnMapPix = ConnMap_SameOcu + ConnMap_DiffOcuU + ConnMap_MirrOcuU;      
    end
    % Put back to Build a mat with the averaged kernal 
    C_PQ_mean(PixInd,:) = reshape(ConnMapPix,N_HCoutX*NPixX*N_HCoutY*NPixY,1)';
end
end