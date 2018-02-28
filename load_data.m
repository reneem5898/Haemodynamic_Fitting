function [time, pressure, ECG, samplerate] = load_data_LVP(directory, filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to load pressure data from mat file - to use with Ohio State 
% University data
% Originally written by: Angus Cheng
% 
% Adapted by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 11 December 2017
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load .mat file
load(sprintf('%s/%s', directory, filename));

% Get time axis
time = (datastart(1):dataend(1))/samplerate(1);

% Get LV pressure data
pressure = data(datastart(1):dataend(1)); 



% There was an offset in SOME of the pressure acquisitions. Check mean value.
% If mean(pressure) > 700, subtract 706 mmHg from all pressure values. 
if mean(pressure) > 700
    
    % Print to log
    disp('Data was offset by 706 mmHg.')
    
    % Offset pressure
    pressure = pressure - 706; %706 mmHg offset
end



% Save ECG data - not always in the same mat file as with AOP and LVP
if length(datastart) == 2

    % Save ECG data
    ECG = data(datastart(2):dataend(2));

else
    
    % Get ECG file name
    file_spec = sprintf('%s/*ECG*L*.mat', directory);
    if length(dir(file_spec)) == 1
        ECG_file = dir(file_spec).name;
    else
        ECG_file = uigetfile(file_spec, 'Choose ECG file corresponding with LVP');
    end
    
    % Load ECG mat file
    load(sprintf('%s/%s', directory, ECG_file));
    
    % Print to log
    disp('Loading ECG data from separate file:');
    disp(sprintf('%s/%s', directory, ECG_file));
    
    ECG = data(datastart(1):dataend(1));
 
end

