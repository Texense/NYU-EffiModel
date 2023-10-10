% generate random binary matrix
% p row q col
% plus and minus randomly selected from n entries
function [PlusInd, MinusInd] = RandPermIndMat(p,q,n)
a = zeros(p,q);
perturbIndBoth = randperm(numel(a), n);
PMInd = logical(randi(2,size(perturbIndBoth,1),size(perturbIndBoth,2))-1);

PlusInd = a; MinusInd = a;
PlusInd(perturbIndBoth(PMInd)) = 1;
MinusInd(perturbIndBoth(~PMInd)) = 1;

PlusInd = logical(PlusInd);
MinusInd = logical(MinusInd);
end