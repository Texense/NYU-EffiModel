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

function [LDEOutS,LibyAll] = LDEIterFunc_Grating_16Func(L4EmeshX, L4ImeshY,...
                                 LDEFrfunc_Subf,...
                                 L4EUse,L4IUse,...
                                 PixInptCtgrUse)
    % first look up the library
    LDEOutLIBy = cell(size(LDEFrfunc_Subf));
    LibyAll = zeros(length(L4EUse),size(LDEFrfunc_Subf,1),size(LDEFrfunc_Subf,2));
    for LGNInd = 1:size(LDEFrfunc_Subf,1)
        for L6Ind = 1:size(LDEFrfunc_Subf,2)
        LDEOutLIBy{LGNInd,L6Ind} = ...
            interp2(L4EmeshX, L4ImeshY,LDEFrfunc_Subf{LGNInd,L6Ind}, L4EUse,L4IUse);%,'makima'
        LibyAll(:,LGNInd,L6Ind) = LDEOutLIBy{LGNInd,L6Ind};
        end
    end
    % Then compose
    LDEOutS = sum(PixInptCtgrUse.*LibyAll,[2,3]);




end