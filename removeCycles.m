function cycles_new = removeCycles(cycles, cycles2keep)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function removes cycles from the pressure data which the user wants
% to omit
%
% Inputs: 1) cycles = cell array of pressure data (output from extract_individual_cycles.m)
%         2) cycles2keep = list of cycle numbers to keep
%
% Output: 1) cycles_new = cell array of pressure data only containing cycles to keep
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 19 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totalCycleList = size(cycles,2);