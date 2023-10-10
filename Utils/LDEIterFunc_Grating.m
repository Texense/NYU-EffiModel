%% Iterations of LDE: Use precomputed function to determine output of a cell type
% Input: L4EmeshX, L4ImeshY  Domain of precomputed functions
%        LDEFrfunc_Subf      Functions of all input type for certain cell
%        population
%        L4EUse,L4IUse       All input L4E L4I
%        LDEUse              Just for output formality
%        varargin:
%        ORTPix,OBLPix,OPTPix,ORTOBLBd_Pix,OPTOBLBd_Pix   Pixel index of different ODs
%        or
%        PixLGNCtgr: n*3

% Zhuo-Cheng Xiao 10/01/2021

function [LDEOutS] = LDEIterFunc_Grating(L4EmeshX, L4ImeshY,...
                                 LDEFrfunc_Subf,...
                                 L4EUse,L4IUse,...
                                 PixLGNCtgr)
    % first look up the library
    LDEOutLIBy = cell(size(LDEFrfunc_Subf));
    LibyAll = zeros(length(L4EUse),length(LDEFrfunc_Subf));
    for FuncInd = 1:length(LDEFrfunc_Subf)
        LDEOutLIBy{FuncInd} = ...
            interp2(L4EmeshX, L4ImeshY,LDEFrfunc_Subf{FuncInd}, L4EUse,L4IUse);
        LibyAll(:,FuncInd) = LDEOutLIBy{FuncInd};
    end
    % Then compose
    LDEOutS = sum(PixLGNCtgr.*LibyAll,2);




end