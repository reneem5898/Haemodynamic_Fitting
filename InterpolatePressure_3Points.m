function [LV_Pressure, MRI_LVP_ed2es, MRI_LVP_es2ds, MRI_LVP_ds2ed] = InterpolatePressure_3Points(timepoints_MRI, LV_cycles, AO_cycles, ED, ES, DS, maxLVP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function takes in the timepoints from MRI cine images and timepoints from
% pressure traces and returns an average pressure trace
%
% Note: function requires THREE identified cardiac time points
%
% Inputs: 1) timepoints_MRI: vector with cardiac MRI timepoints and number
% of frames (1 x 4)
%         2) LV_cycles: all LV pressure traces
%         3) AO_cycles: all AO pressure traces
%         4) ED: ED timepoints from pressure traces
%         5) ES: ES timepoints from pressure traces
%         6) DS: DS timepoints from pressure traces
%         7) maxLVP: timepoints of maximum LV pressure from traces
%
% Output: 1) LV_Pressure: vector of interpolated LV pressure trace
%
% Written by: Renee Miller (renee.miller@auckland.ac.nz)
% Date: 22 June 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MRI cardiac time points
edMRI = timepoints_MRI(1);
esMRI = timepoints_MRI(2);
dsMRI = timepoints_MRI(3);

% Get number of frames for each segment of the pressure curve
num_MRI_frames_ed2es = esMRI - edMRI + 1; % Number of frames in images from ED to end ES 
num_MRI_frames_es2ds = dsMRI - esMRI + 1; % Number of frames in images from ES to end DS 
num_MRI_frames_ds2ed = numFrames - dsMRI + edMRI; % Number of frames in images from end DS to ED 

c = 1; % counter
ES_LV = zeros(size(LV_cycles,2)*size(AO_cycles,2),1); % Initialise variable to save all possible ES indices

% Initialise variables to save pressures at MRI frames
% Sizes differ depending on options
MRI_LVP_ed2es = zeros(size(LV_cycles,2)*size(AO_cycles,2), num_MRI_frames_ed2es); 
MRI_LVP_es2ds = zeros(size(LV_cycles,2)*size(AO_cycles,2), num_MRI_frames_es2ds-1); % -1 because data will be exclusive of eivc
MRI_LVP_ds2ed = zeros(size(LV_cycles,2), num_MRI_frames_ds2ed); 

for i = 1:size(LV_cycles,2)
    
    LVP_c = LV_cycles{2,i}; % Current LV pressure cycle
    
%     LVP_ed2MaxLVP = LVP_c(ED(i):maxLVP(i)); % Section from ED to maximum LVP
    LVP_maxLVP2ds = LVP_c(maxLVP(i):DS(i)); % Section from maximum LVP to DS
    
    
    %% DS to ED
    
    % Interpolate pressure at each MRI frame between ds and ed
    pressure_ds2ed = LVP_c(DS(i):end); % LV pressure between DS and ED
    interpPressure = LinearInterpolatePressure(pressure_ds2ed, length(pressure_ds2ed), num_MRI_frames_ds2ed);
    MRI_LVP_ds2ed(i,:) = interpPressure;
    
    
    for j = 1:size(AO_cycles,2)
        
        AOP_c = AO_cycles{2,j}; % Current AO pressure cycle
        
        %% ED to ES       
        
        % Get point in LVP when LVP == AOP(ES)
        AOP_ES = AOP_c(ES(j)); % Aortic pressure at ES
        ES_LV(c) = find(abs(LVP_maxLVP2ds - AOP_ES) == min(abs(LVP_maxLVP2ds - AOP_ES))) + maxLVP(i);
         
        % Interpolate pressure at MRI image frames from ED to ES
        pressure_ed2es = LVP_c(ED(c):ES_LV(c));
        interpPressure = LinearInterpolatePressure(pressure_ed2es, length(pressure_ed2es), num_MRI_frames_ed2es);
        MRI_LVP_ed2es(c,:) = interpPressure(2:end); % Exclude eIVC
        
        %% ES to DS
        
        pressure_es2ds = LVP_c(ES_LV(c):DS(i));
        interpPressure = LinearInterpolatePressure(pressure_es2ds, length(pressure_es2ds), num_MRI_frames_es2ds);
        MRI_LVP_es2ds(c,:) = interpPressure(2:end); % Exclude ES
        
        % Count up
        c = c + 1;
        
    end
end

% Average pressures and put together segments
%                   ED to ES                 ES to DS             DS to ED                 
LV_Pressure = [mean(MRI_LVP_ed2es,1), mean(MRI_LVP_es2ds,1), mean(MRI_LVP_ds2ed,1)]';
