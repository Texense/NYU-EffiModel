%% Get Input for the current time bin...
% Input: Events: Event matrix N*X. Entries are the event time for each
% neuron (indexed by row number)
%        TimeN:  Number of time bin
%        dt:     timestep
%        N:      Number of E or I neurons
% STASHED!!!
function InputVec = GetInputFromSeries(Events,TimeN,dt,N)

[A,~,~] = find(Events>(TimeN-1)*dt & Events<(TimeN)*dt); % How the neuron indexes of events in the current time bin
% [GC,GR] = groupcounts(A); % Get Unique neuind and occurence
% InputVec = zeros(N,1);
% InputVec(GR) = GC;
GR = unique(A);
n  = histcounts(A,GR);
