%% function C_SS_mean = AveSpatKer(C_SS_Pixel_Us,N_HC,NPixX,NPixY)
% Input: C_SS_Pixel_Us: A Pixel connectivity matrix (from network).
% Pre(col)*Post(row)
%        N_HC,NPixX,NPixY: # of hypercolumns, # of pixels per HC
% Output: C_SS_mean: Averaged spatial kernal of firing rats. Same size as
% C_SS_Pixel_Us

% Zhuo-Cheng Xiao 09/16/2021
% New Version: Want to expand from 3HC to arbitrary # of HC
% N_HCout should be larger than 3
function C_SS_mean = AveSpatKer(C_PQ_Pixel_Us,N_HCin,N_HCout,NPixX,NPixY)
%% Check if C_SS is an eligible matrix
if size(C_PQ_Pixel_Us,1) ~= size(C_PQ_Pixel_Us,2)
    error('Input Matrix not square!')
end
if size(C_PQ_Pixel_Us,2) ~= N_HCin*NPixX*N_HCin*NPixY
    error('# of pixels doesnt match!')
end

if N_HCin<3
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
    ConnVec = C_PQ_Pixel_Us(PixInd,:)';
    ConnMap = reshape(ConnVec,N_HCin*NPixY,N_HCin*NPixX);
    Cen = floor((N_HCin-1)/2);
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
C_SS_mean  = zeros(N_HCout*NPixX*N_HCout*NPixY);

for PixInd = 1:N_HCout*NPixX*N_HCout*NPixY
    PX = ceil(PixInd/(N_HCout*NPixY));
    PY = mod(PixInd,N_HCout*NPixY);
    if PY == 0
        PY = N_HCout*NPixY;
    end
    
    if mod(N_HCout,2) == 1
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
    else % If even: Periodic boundary conditions
        ConnMap_Rev = zeros(N_HCout*NPixY,N_HCout*NPixX);
        YRange = mod(PY-NPixY:PY+NPixY-1,N_HCout*NPixY); YRange(YRange==0) = N_HCout*NPixY;
        XRange = mod(PX-NPixX:PX+NPixX-1,N_HCout*NPixX); XRange(XRange==0) = N_HCout*NPixX;
        ConnMap_Rev(YRange,XRange) = Ker_PQ_mean;
        
        ConnMapPix = ConnMap_Rev;
        
%         ConnMap_Rev(PY_New-NPixY:PY_New+NPixY-1,...
%                     PX_New-NPixX:PX_New+NPixX-1) = Ker_SS_mean;
%               ConnMap_Rev(            NPixY+1:          2*NPixY, NPixX+1:(N_HCout+1)*NPixX) ...% First row
%             = ConnMap_Rev(            NPixY+1:          2*NPixY, NPixX+1:(N_HCout+1)*NPixX) ...
%             + ConnMap_Rev((N_HCout+1)*NPixY+1:(N_HCout+2)*NPixY, NPixX+1:(N_HCout+1)*NPixX);
%          
%               ConnMap_Rev(N_HCout*NPixY+1:(N_HCout+1)*NPixY, NPixX+1:(N_HCout+1)*NPixX) ...% last row
%             = ConnMap_Rev(N_HCout*NPixY+1:(N_HCout+1)*NPixY, NPixX+1:(N_HCout+1)*NPixX) ...
%             + ConnMap_Rev(              1:            NPixY, NPixX+1:(N_HCout+1)*NPixX);
%          
%               ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY,             NPixX+1:          2*NPixX) ... % First Col
%             = ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY,             NPixX+1:          2*NPixX) ...
%             + ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY, (N_HCout+1)*NPixX+1:(N_HCout+2)*NPixX);
%          
%               ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY, N_HCout*NPixX+1:(N_HCout+1)*NPixX) ... % last Col
%             = ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY, N_HCout*NPixX+1:(N_HCout+1)*NPixX) ...
%             + ConnMap_Rev(NPixY+1:(N_HCout+1)*NPixY,               1:            NPixX);
    end
    % Put back to Build a mat with the averaged kernal 
    C_SS_mean(PixInd,:) = reshape(ConnMapPix,N_HCout*NPixX*N_HCout*NPixY,1)';
end
end