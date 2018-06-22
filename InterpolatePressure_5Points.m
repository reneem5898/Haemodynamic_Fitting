function [LV_Pressure, MRI_LVP_ed2eivc, MRI_LVP_eivc2es, MRI_LVP_es2eivr, MRI_LVP_eivr2ds, MRI_LVP_ds2ed] = InterpolatePressure_5Points(timepoints_MRI, LV_cycles, AO_cycles, ED, eIVC, ES, eIVR, DS, maxLVP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function takes in the timepoints from MRI cine images and timepoints from
% pressure traces and returns an average pressure trace
%
% Note: function requires FIVE identified cardiac time points
%
% Inputs: 1) timepoints_MRI: vector with cardiac MRI timepoints and number of frames (1 x 6)
%         2) LV_cycles: all LV pressure traces
%         3) AO_cycles: all AO pressure traces
%         4) ED: ED timepoints from pressure traces
%         5) eIVC: eIVC timepoints from pressure traces
%         6) ES: ES timepoints from pressure traces
%         7) eIVR: eIVR timepoints from pressure traces
%         8) DS: DS timepoints from pressure traces
%         9) maxLVP: timepoints of maximum LV pressure
%
% Output: 1) LV_Pressure: vector of interpolated LV pressure trace
%
% Written by: Renee Miller (renee.miller@auckland.ac.nz)
% Date: 22 June 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MRI cardiac time points
edMRI = timepoints_MRI(1);
eivcMRI = timepoints_MRI(2);
esMRI = timepoints_MRI(3);
eivrMRI = timepoints_MRI(4);
dsMRI = timepoints_MRI(5);

% Total number of frames
numFrames = timepoints_MRI(6);

% Get number of frames for each segment of the pressure curve
num_MRI_frames_ed2eivc = eivcMRI - edMRI + 1; % Number of frames in images from ED to end IVC 
num_MRI_frames_eivc2es = esMRI - eivcMRI + 1; % Number of frames in images from end IVC to ES 
num_MRI_frames_es2eivr = eivrMRI - esMRI + 1; % Number of frames in images from ES to end IVR 
num_MRI_frames_eivr2ds = dsMRI - eivrMRI + 1; % Number of frames in images from end IVR to DS 
num_MRI_frames_ds2ed = numFrames - dsMRI + edMRI; % Number of frames in images from DS to ED (excluding last point)

c = 1; % counter
eIVC_LV = zeros(size(LV_cycles,2)*size(AO_cycles,2),1); % Initialise variable to save all possible eIVC points
ES_LV = zeros(size(LV_cycles,2)*size(AO_cycles,2),1); % Initialise variable to save all possible ES indices

% Initialise variables to save pressures at MRI frames
% Sizes differ depending on options
MRI_LVP_ed2eivc = zeros(size(LV_cycles,2)*size(AO_cycles,2), num_MRI_frames_ed2eivc); 
MRI_LVP_eivc2es = zeros(size(LV_cycles,2)*size(AO_cycles,2), num_MRI_frames_eivc2es-1); % -1 because data will be exclusive of eivc
MRI_LVP_es2eivr = zeros(size(LV_cycles,2)*size(AO_cycles,2), num_MRI_frames_es2eivr-1); % -1 because data will be exclusive of es
MRI_LVP_eivr2ds = zeros(size(LV_cycles,2), num_MRI_frames_eivr2ds-2); % -1 because data will be exclusive of eivr AND ds
MRI_LVP_ds2ed = zeros(size(LV_cycles,2), num_MRI_frames_ds2ed); 

for i = 1:size(LV_cycles,2)
    
    LVP_c = LV_cycles{2,i}; % Current LV pressure cycle
    
    LVP_ed2MaxLVP = LVP_c(ED(i):maxLVP(i)); % Section from ED to maximum LVP
    LVP_maxLVP2ds = LVP_c(maxLVP(i):DS(i)); % Section from maximum LVP to DS
    
    %% eIVR to DS
    
    % Interpolate pressure at each MRI frame between ds and ed
    pressure_eivr2ds = LVP_c(eIVR(i):DS(i)); % LV pressure between DS and ED
    interpPressure = LinearInterpolatePressure(pressure_eivr2ds, length(pressure_eivr2ds), num_MRI_frames_eivr2ds);
    MRI_LVP_eivr2ds(i,:) = interpPressure(2:end-1); % Exclude end IVR
    
    %% DS to ED
    
    % Interpolate pressure at each MRI frame between ds and ed
    pressure_ds2ed = LVP_c(DS(i):end); % LV pressure between DS and ED
    interpPressure = LinearInterpolatePressure(pressure_ds2ed, length(pressure_ds2ed), num_MRI_frames_ds2ed);
    MRI_LVP_ds2ed(i,:) = interpPressure;
    
    for j = 1:size(AO_cycles,2)
        
        AOP_c = AO_cycles{2,j}; % Current AO pressure cycle
        
        %% ED to eIVC
        
        % Get point in LVP when LVP == AOP(eIVC)
        AOP_eIVC = AOP_c(eIVC(j)); % Aortic pressure at endIVC
        eIVC_LV(c) = find(abs(LVP_ed2MaxLVP - AOP_eIVC) == min(abs(LVP_ed2MaxLVP - AOP_eIVC)));
        
        % Interpolate pressure at MRI image frames from ED to eIVC
        pressure_ed2eivc = LVP_c(ED(i):eIVC_LV(c));
        interpPressure = LinearInterpolatePressure(pressure_ed2eivc, length(pressure_ed2eivc), num_MRI_frames_ed2eivc);
        MRI_LVP_ed2eivc(c,:) = interpPressure;
        
        %% eIVC to ES
        
        % Get point in LVP when LVP == AOP(ES)
        AOP_ES = AOP_c(ES(j)); % Aortic pressure at ES
        ES_LV(c) = find(abs(LVP_maxLVP2ds - AOP_ES) == min(abs(LVP_maxLVP2ds - AOP_ES))) + maxLVP(i);
        %ESt_LV(c) = LV_cycles{1,i}(ES_LV(c)) * LV_cycles{1,i}(DS(i))/mean(LV_cycles{1,i}(DS(:))); % DON'T UNDERSTAND THIS LINE - WHY SCALE BY DS TIMING????
        
        % Interpolate pressure at MRI image frames from eIVC to ES
        pressure_eivc2es = LVP_c(eIVC_LV(c):ES_LV(c));
        interpPressure = LinearInterpolatePressure(pressure_eivc2es, length(pressure_eivc2es), num_MRI_frames_eivc2es);
        MRI_LVP_eivc2es(c,:) = interpPressure(2:end); % Exclude eIVC
        
        %% ES to eIVR
        
        pressure_es2eivr = LVP_c(ES_LV(c):eIVR(i));
        interpPressure = LinearInterpolatePressure(pressure_es2eivr, length(pressure_es2eivr), num_MRI_frames_es2eivr);
        MRI_LVP_es2eivr(c,:) = interpPressure(2:end); % Exclude ES
        
        % Count up
        c = c + 1;
        
    end
end

% Average pressures and put together segments
%                   ED to eIVC                 eIVC to ES             ES to eIVR                 eIVR to DS               DS to ED
LV_Pressure = [mean(MRI_LVP_ed2eivc,1), mean(MRI_LVP_eivc2es,1), mean(MRI_LVP_es2eivr,1), mean(MRI_LVP_eivr2ds,1), mean(MRI_LVP_ds2ed,1)]';
