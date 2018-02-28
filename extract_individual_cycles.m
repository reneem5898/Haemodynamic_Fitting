function [cycles] = extract_individual_cycles(time, R_positions, cyclesToKeep, Pressure, ECG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses the RR intervals to extract each individual cycle from 
% a pressure trace 
%
% Inputs: 1) time - vector over entire trace
%         2) R_positions - indices of R peaks
%         3) cyclesToKeep - list of cycles to keep in pressure waveform
%         4) Pressure - vector of entire pressure trace
%         5) ECG - vector of ECG data over entire trace
%
% Output: 1) cycles - cell array: row 1 = time vectors
%                                 row 2 = pressure vectors
%                                 row 3 = ecg vectors
%
% Originally written by: Angus Cheng
% Adapted by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 11 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise variable to save each pressure cycle
cycles = cell(3, length(cyclesToKeep));

% Loop through number of cycles
for i = 1:length(cyclesToKeep)
    
    % Cycle number
    cycleNum = cyclesToKeep(i);
    
    % Time
    cycle_time = time(R_positions(cycleNum+1)) - time(R_positions(cycleNum)); % Length of cycle in time
    cycle_points = (R_positions(cycleNum+1)-R_positions(cycleNum))+1; % Number of points in cycle
    cycles{1,i} = linspace(0, cycle_time, cycle_points); % time
    
    % Pressure
    cycles{2,i} = Pressure(R_positions(cycleNum):R_positions(cycleNum+1)); 
    
    % ECG
    cycles{3,i} = ECG(R_positions(cycleNum):R_positions(cycleNum+1));
end

end

