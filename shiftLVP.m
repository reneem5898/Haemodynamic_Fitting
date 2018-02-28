function shifted_LVP = shiftLVP(LVP, DS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes a single cycle for an LV pressure trace as well as 
% the diastasis index and then shifts the entire LV pressure trace so that
% LVP(DS) = 0
%
% Inputs: 1) LVP = 1 cycle of LV pressure trace
%         2) DS = index of diastasis time point
%
% Output: shifted_LVP = LV pressure trace shifted up or down
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 18 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DS_pressure = LVP(DS);
shifted_LVP = LVP - DS_pressure;