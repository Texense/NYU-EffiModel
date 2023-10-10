%% Modify rate map
% CurrentFrE: The rough ratemap
% ConvMat: Convolution kernel
% CurrentBlowupModi: Blowup mat. Should be same dim of CurrentFrE. We get
% rid of where blowup is true
% CurrentsteadyModi: Steady mat. Should be same dim of CurrentFrE. We get
% rid of where steady is false.

% Output: FringMap
function FringMap = FrRateModify(CurrentFrE, ConvMat,CurrentBlowupModi,CurrentsteadyModi)
%CurrentFrE(CurrentFrE<eps) = nan; 
CurrentFrE(CurrentBlowupModi) = nan;
CurrentFrE(~CurrentsteadyModi) = nan;

CurrentFrE = conv2(CurrentFrE,ConvMat);
[a,b] = size(ConvMat);

% Cut out edges after convolution
% Cutting out 'a-1' rows in all
if a>1
    a_half1 = floor((a-1)/2);
    a_half2 = (a-1) - a_half1;
    CurrentFrE(end-a+2:end-a_half1,:) = nan;
    CurrentFrE(end-a_half1+1:end,:) = [];
    CurrentFrE(a_half2+1:a-1,  :) = nan;
    CurrentFrE(1:a_half2,      :) = [];
end

if b>1
    b_half1 = floor((b-1)/2);
    b_half2 = (b-1) - a_half1;    
    CurrentFrE(:,end-b+2:end-b_half1) = nan;
    CurrentFrE(:,end-b_half1+1:end) = [];
    CurrentFrE(:,b_half2+1:b-1) = nan;
    CurrentFrE(:,1:b_half2      ) = [];
end

FringMap = CurrentFrE;
end