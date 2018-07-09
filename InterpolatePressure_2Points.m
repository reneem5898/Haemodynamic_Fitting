function [LV_Pressure, MRI_LVP_ed2ds, MRI_LVP_ds2ed] = InterpolatePressure_2Points(timepoints_MRI, LV_cycles, DS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function takes in the timepoints from MRI cine images and timepoints from
% pressure traces and returns an average pressure trace
%
% Note: function requires THREE identified cardiac time points
%
% Inputs: 1) timepoints_MRI: vector with cardiac MRI timepoints and number
% of frames (1 x 4)
%         2) LV_cycles: all LV pressure traces
%         4) DS: DS timepoints from pressure traces
%
% Output: 1) LV_Pressure: vector of interpolated LV pressure trace
%
% NOTE: Function assumes that ED is frame #1
%
% Written by: Renee Miller (renee.miller@auckland.ac.nz)
% Date: 26 June 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MRI cardiac time points
edMRI = timepoints_MRI(1);
dsMRI = timepoints_MRI(2);

% Number of frames
numFrames = timepoints_MRI(3);

% Get number of frames for each segment of the pressure curve
num_MRI_frames_ed2ds = dsMRI - edMRI + 1; % Number of frames in images from ED to end ES 
num_MRI_frames_ds2ed = numFrames - dsMRI + edMRI; % Number of frames in images from end DS to ED 

% Initialise variables to save pressures at MRI frames
% Sizes differ depending on options
MRI_LVP_ed2ds = zeros(size(LV_cycles,2), num_MRI_frames_ed2ds - 1); % Will exclude ED and DS
MRI_LVP_ds2ed = zeros(size(LV_cycles,2), num_MRI_frames_ds2ed); 

for i = 1:size(LV_cycles,2)
    
    LVP_c = LV_cycles{2,i}; % Current LV pressure cycle
    
    %% DS to ED
    
    % Interpolate pressure at each MRI frame between ds and ed
    pressure_ds2ed = LVP_c(DS(i):end); % LV pressure between DS and ED
    interpPressure = LinearInterpolatePressure(pressure_ds2ed, length(pressure_ds2ed), num_MRI_frames_ds2ed);
    MRI_LVP_ds2ed(i,:) = interpPressure;
    
    %% ED to DS
    
    % Interpolate pressure at each MRI frame between ds and ed
    pressure_ed2ds = LVP_c(1:DS(i)); % LV pressure between DS and ED
    interpPressure = LinearInterpolatePressure(pressure_ed2ds, length(pressure_ed2ds), num_MRI_frames_ed2ds);
    MRI_LVP_ed2ds(i,:) = interpPressure(1:end-1);    

end

% Average pressures and put together segments
%                   ED to DS                 DS to ED                 
LV_Pressure = [mean(MRI_LVP_ed2ds,1), mean(MRI_LVP_ds2ed,1)]';
