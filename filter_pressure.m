function [Pressure_Filtered] = filter_pressure(pressure, samplerate, cut_off_freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function filters pressure using a 4th order low pass filter 
% using a specified cutoff frequency
%
% Inputs: 1) pressure (vector)
%         2) samplerate (single value)
%         3) cut_off_freq (single value, e.g. 25)
%
% Output: 1) Pressure_Filtered = filtered pressure (vector)
%
% Originally written by: Angus Cheng
% Adapted by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 11 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create butterworth filter
Wnl = cut_off_freq/(samplerate/2);
[D, C] = butter(4, Wnl, 'low');

% Filter data
Pressure_Filtered = filtfilt(D, C, pressure);

end

