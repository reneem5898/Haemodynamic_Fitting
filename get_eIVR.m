function eIVR = get_eIVR(LVP, d2LVP_dt2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function finds eIVR - end isovolumetric relaxation
%
% Input: 1) LVP = LV pressure (vector, 1 cycle)
%        2) d2LVP_dt2 = second derivative of LV pressure trace
%
% Output: index of eIVR in LVP
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 20 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find peaks in second derivative pressure trace
[peaks, locs] = findpeaks(d2LVP_dt2); 

% Get rid of first and last peak
peaks = peaks(2:end-1); 
locs = locs(2:end-1); 

% Get biggest peak
eIVR = locs(peaks == max(peaks)); 

% Plot second derivative
t = linspace(1, length(LVP), length(LVP));
f = figure;
subplot(2,1,1)
plot(t, d2LVP_dt2)
hold on
scatter(eIVR, d2LVP_dt2(eIVR), 'ro', 'filled')
xlabel('Time', 'FontSize', 12)
ylabel('d2P/dt2', 'FontSize', 12)

% Plot LVP with boundaries and eIVR
subplot(2,1,2)
plot(t, LVP)
hold on
scatter(eIVR, LVP(eIVR), 'ro', 'filled')
title('eIVR Selection')
ylabel('Pressure (mmHG)', 'FontSize', 12)
xlabel('Time', 'FontSize', 12)
text(eIVR+3, LVP(eIVR), 'eIVR', 'FontSize', 12)
pause(1)
close(f)