%% A function from Neuron Ind and spatial scales to Orientation Domains
% Input: ODNum: 4,6,8,...Number of orientation domain per HC %% Should be
% even!!!!
%        NXIND, NYIND: x and y ind of the neuron
%        N_HC : Number of E/I neurons Per HyperColumn
% Output: The index of the orientation domain,1-ODNum
%% Th core idea is 1. divide each HC into ODNum angles from its center
%                  2. label angle domains 1-ODNum. 
%                  3. Based on the HC index, the order of angle domains to
%                  OD Domains
%                  4. First HC clockwise, and any neibours go for x/y-axis
%                  symmetry. So HC OD index becomes periodic for every 2*2
%                  HCs.
function ODInd = OrientDom(ODNum,NXInd,NYInd,Neu_HC)
%% First, determine HC Index
HC_X = ceil(NXInd/Neu_HC);
HC_Y = ceil(NYInd/Neu_HC);
HCMode = mod(HC_X,2)+ 2*(mod(HC_Y,2)-1);

% Depending on the HCMode, assign the order of OD. 
switch HCMode
    % Original: Starting from the angle at 3 o'clock, go counter clockwise
    case 1
    ODorder = 1:ODNum;
    % Y-axis symmetry    
    case 0
    ODorder = [floor(ODNum/2)+1:-1:1, ODNum:-1:floor(ODNum/2)+2];
    % X-axis symmetry    
    case -1
    ODorder = [1, ODNum:-1:2];
    % X and Y    
    case -2    
    ODorder = [floor(ODNum/2)+1:ODNum, 1:floor(ODNum/2)]; 
    otherwise
        disp('Unrecognized HC Mode')
        return
end

%% Second, check the angle phase of the neuron in OD
MidPoint = Neu_HC/2; %x and y ind of the middle point
NX_HC = NXInd - (HC_X-1)*Neu_HC;
NY_HC = NYInd - (HC_Y-1)*Neu_HC; % x and y ind of the neuron inside the HC

% Get the angle, and assign angle in phases
if isempty(find([NY_HC - MidPoint, NX_HC - MidPoint], 1)) 
    % If the neuron lies exactly at the center, then...
    ODInd = ODorder(1);
else
    Angle_HC = angle(NX_HC - MidPoint + 1i*(NY_HC - MidPoint)); % get angle by atan
    % First, set range in 0-2pi and rotate conterclockwise for half phase
    % Then get the angle number
    Angle_Phase = ceil(mod(Angle_HC+pi/ODNum, 2*pi) / (2*pi/ODNum)); 
    if Angle_Phase == 0
       Angle_Phase = ODNum;
    end
    ODInd = ODorder(Angle_Phase);
end

end