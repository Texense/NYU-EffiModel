%% Inhibition suppresion
% input:
%        LDEUseI  I firing rates
%        Thrsld    Threshold we start suppresion, i.e., below Thrsld should
%        be identity map
%        Highist   [largestinput, largestoutput]
%        varargin1: hard bound
%        varargin2: Linear cap or constant cap
% Using second order
function ...
    LDEUseIKill = InhKill(LDEUseI, Thrsld, Highist,varargin)

%% First make up a substitute function
x0 = Thrsld; x1 = Highist(1); y1 = Highist(2);
if length(varargin)<1
    hardbound = 120; % back again to linear
else
    hardbound = varargin{1};
end
if length(varargin)<2
    LinearSlop = 0; 
else
    LinearSlop = varargin{2};
end

a = (y1-x1)/(x1-x0)^2;
b = 1 - 2*a*x0;
c = x0 - a*x0^2 - b*x0;

f = @(x) a*x.^2 + b*x + c;

% back again to linear
d = f(hardbound);
f1 = @(x) LinearSlop*(x-hardbound)+d;

% Manupulate I frs
LDEUseIKill = LDEUseI;
KillInd = LDEUseIKill>Thrsld & LDEUseIKill<=hardbound;
LDEUseIKill(KillInd) = f(LDEUseIKill(KillInd));

LinearInd = LDEUseIKill>hardbound;

LDEUseIKill(LinearInd) = f1(LDEUseIKill(LinearInd));

end