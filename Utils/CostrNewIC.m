%% Use old terminate condition to get the new initial condition for NW simulation
% If Neuron Num not identical, then merge or create new ones
% If Delay Mats not identical, also adjust.

function [RefTimeE, VE, SpE, GE_ampa_R, GE_nmda_R, GE_gaba_R, GE_ampa_D, GE_nmda_D, GE_gaba_D,...
          RefTimeI, VI, SpI, GI_ampa_R, GI_nmda_R, GI_gaba_R, GI_ampa_D, GI_nmda_D, GI_gaba_D,...
          EEDlyRcd, IEDlyRcd, EIDlyRcd]...
         = CostrNewIC(EndState,N_E,N_I,N_EEDly,N_IEDly,N_EIDly)
%% First, get old neuron number
OldNE = length(EndState{1});
OldNI = length(EndState{10});
% predefine Outputs
EndStateU = EndState;
if length(EndStateU)<21 % 18 fields for EnI, 3 for delay matrices
    for FieldInd = length(EndStateU)+1:21
        EndStateU{FieldInd} = [];
    end
end
% if Old simulation has more neurons, then randomly pick a subset
if OldNE>N_E
    NewSubE = randsample(OldNE,N_E);
    NewSubI = randsample(OldNI,N_I);
    for FieldInd = 1:9
        EndStateU{FieldInd} = MatShrink(EndState{FieldInd},NewSubE);
    end
    for FieldInd = 10:18
        EndStateU{FieldInd} = MatShrink(EndState{FieldInd},NewSubI);
    end
    for FieldInd = [19,21]
        if ~isempty(EndStateU{FieldInd})
        EndStateU{FieldInd} = MatShrink(EndState{FieldInd},NewSubE);
        end
    end
        if ~isempty(EndStateU{20})
        EndStateU{20} = MatShrink(EndState{20},NewSubI);
        end
    
elseif OldNE<N_E
    NewSubE = datasample(1:OldNE,N_E-OldNE);
    NewSubI = datasample(1:OldNI,N_I-OldNI);
    for FieldInd = 1:9
        EndStateU{FieldInd} = MatExtend(EndState{FieldInd},NewSubE);
    end
    for FieldInd = 10:18
        EndStateU{FieldInd} = MatExtend(EndState{FieldInd},NewSubI);
    end
    for FieldInd = [19,21]
        if ~isempty(EndStateU{FieldInd})
        EndStateU{FieldInd} = MatExtend(EndState{FieldInd},NewSubE);
        end
    end
    if ~isempty(EndStateU{20})
        EndStateU{20} = MatExtend(EndState{20},NewSubI);
    end
    
end     
     
     
RefTimeE = EndStateU{1}; VE = EndStateU{2}; SpE = EndStateU{3}; GE_ampa_R = EndStateU{4}; GE_nmda_R = EndStateU{5}; GE_gaba_R = EndStateU{6};
                                                             GE_ampa_D = EndStateU{7}; GE_nmda_D = EndStateU{8}; GE_gaba_D = EndStateU{9};
RefTimeI = EndStateU{10}; VI = EndStateU{11}; SpI = EndStateU{12}; GI_ampa_R = EndStateU{13}; GI_nmda_R = EndStateU{14}; GI_gaba_R = EndStateU{15};
                                                                GI_ampa_D = EndStateU{16}; GI_nmda_D = EndStateU{17}; GI_gaba_D = EndStateU{18};

%% Second, adjust delay matrices
if length(EndStateU)<21 % 18 fields for EnI, 3 for delay matrices
    for FieldInd = length(EndStateU)+1:21
        EndStateU{FieldInd} = [];
    end
end

if N_EEDly>=size(EndStateU{19},2)
EEDlyRcd = [EndStateU{19},zeros(N_E,N_EEDly-size(EndStateU{19},2))]; 
else
    EEDlyRcd = EndStateU{19}(:,1:N_EEDly);
end
if N_IEDly>=size(EndStateU{20},2)
IEDlyRcd = [EndStateU{20},zeros(N_I,N_IEDly-size(EndStateU{20},2))]; 
else
    IEDlyRcd = EndStateU{20}(:,1:N_IEDly);
end
if N_EIDly>=size(EndStateU{21},2)
EIDlyRcd = [EndStateU{21},zeros(N_E,N_EIDly-size(EndStateU{20},2))]; 
else
    EIDlyRcd = EndStateU{21}(:,1:N_EIDly);
end

end


%% Select from matrix
function MatPost = MatShrink(Mat,NewSub)
MatPost = Mat(NewSub,:);
end

%% Extend from matrix
function MatPost = MatExtend(Mat,NewSub)
MatPost = [Mat;Mat(NewSub,:)];
end