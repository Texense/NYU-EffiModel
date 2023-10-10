%% Generale REAL Poisson processes for each cell;
% Input: N: Number of Neurons
%        lambda: Rates of Poisson; Can be N*1 vec or 1 number
%        T: Total time (Unit: ms)
%        Normal: boolean: normalization or not

% Output: eventsN: N*X Mat containing the series time for each neuron.
%         EventTimeMat: N*timebins, entries are number of events of
%         corresponding neuron in bin
% single by default
function EventTimeMat = PoissonInputForNetwork(N,lambda,T,dt, varargin)
switch length(lambda)
    case 1
        lambda_use = lambda*ones(N,1);
    case N
        lambda_use = lambda;
    otherwise
        disp('Illigal rate form')
        return
end
% Use a two-time 
NStoreEvent = floor(T*max(lambda)*1.2);
eventsN = zeros(N,NStoreEvent);

for NeuInd = 1:N
    OneNeuEvent = Poisson_Process_2(lambda_use(NeuInd),T);
    eventsN(NeuInd,1:length(OneNeuEvent)) = OneNeuEvent';
end

[A,~,Times] = find(eventsN);
EventTimeMat = full(sparse(A,floor(Times/dt)+1,ones(size(A)),N,floor(T/dt)+1));

if ~isempty(varargin)
    normal = varargin{1};
    if normal
        norm4sum = sum(EventTimeMat,2)/T./lambda;
        EventTimeMat = EventTimeMat./repmat(norm4sum,1,size(EventTimeMat,2));
    end
    
end

end


function arrt = Poisson_Process_2(lambda,T)
npoints = poissrnd(lambda*T);
% When event number N is fixed, the event time are distributed averagely.
if (npoints>0)
  arrt = [0; sort(rand(npoints,1)*T)];
else
  arrt = [];
end
end