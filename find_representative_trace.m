function [representative_times, representative_pressure, representative_ECG] = find_representative_trace(RR, cycles, samplerate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function finds representative pressure trace by normalising them over time
% and averaging pressure traces
% 
% Inputs: 1) RR - R-R intervals
%         2) cycles - cell array (output from extract_individual_cycles.m)
%         3) samplerate - sampling frequency
%
% Outputs: 1) representative_times - time vector
%          2) representative_pressure - pressure vector
%          3) representative_ECG - ecg vector
%
% Originally written by: Angus Cheng
% Adapted by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 12 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define time scale for representative trace

% Mean RR interval
average_LV_time = mean(RR);
representative_times = 0:1/samplerate:average_LV_time;

% Normalise times for each individual cycle
for i = 1:length(RR)
    for j = 1:length(cycles{1,i})
        cycles{1,i}(j) = (cycles{1,i}(j) - min(cycles{1,i}(:)))/(max(cycles{1,i}(:))-min(cycles{1,i}(:)))*average_LV_time;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract closest time information for each trace

% Initialise cell array to store cycle pressures at normalised time point
representative_cycles = cell(3,length(RR));

for i = 1:length(RR) % Loop through number of cycles
    for j = 1:length(representative_times) % Loop through each time point in the representative trace
        
        representative_cycles{1,i} = representative_times; % time

        
        difference = abs(cycles{1,i}-representative_times(j));
        index = find(difference == min(difference));
        representative_cycles{2,i}(j) = cycles{2,i}(index(1));
        representative_cycles{3,i}(j) = cycles{3,i}(index(1));
    end
end

% average by extracting nearest time information
representative_pressure = zeros(1,length(representative_times));
representative_ECG = zeros(1,length(representative_times));
for i = 1:length(RR)
    for j = 1:length(representative_times)
        representative_pressure(j) = representative_pressure(j)+representative_cycles{2,i}(j);
        representative_ECG(j) = representative_ECG(j)+representative_cycles{3,i}(j);
    end
end

representative_pressure = (representative_pressure)/(length(RR));
representative_ECG = (representative_ECG)/(length(RR));



end

