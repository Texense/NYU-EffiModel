%% Just make sure 22.5 is using the same function at all.....
% LDEFrfuncS 4*4 cell
% mapping: [1,2,3,4;a,b,c,d]; generally 4,3,2,1
function LDEFrfuncSOut = FuncSwap(LDEFrfuncS,mapping)
LDEFrfuncSOut = LDEFrfuncS;
for ii = 1:4
    for jj = 1:4
        LDEFrfuncSOut{ii,jj} = ...
            (LDEFrfuncS{ii,jj}           + LDEFrfuncS{mapping(2,ii),jj} + ...
            LDEFrfuncS{ii,mapping(2,jj)} + LDEFrfuncS{mapping(2,ii),mapping(2,jj)})/4;
    end
    
end
%ii,jj,mapping(2,ii),jj,ii,mapping(2,jj),mapping(2,ii),mapping(2,jj)

end