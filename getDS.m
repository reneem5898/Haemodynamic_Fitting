function DS = getDS(LVP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function asks user to choose DS
%
% Input: LVP = LV pressure (vector, 1 cycle)
%
% Output: DS = index of diastasis in LVP
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 16 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE TO USER: In the pig data - have NO CLUE how to find diastasis. Good luck.

% X axis
pt = linspace(1, length(LVP), length(LVP));

% Plot cycle with minima identified
FH = figure('position', [100 100 600 600]);
plot(pt, LVP, '-', 'markersize', 2);
title('LVP')
xlabel('Time')
ylabel('Pressure (mmHg)')

% Get user to choose diastasis point
waitfor(msgbox('Select the approximate location of diastasis and press enter.', 'Dicrotic notch'));
[x,~] = getpts(FH); % Get user to choose approximate location of dicrotic notch
DS = round(x);

% Close figure
close(FH);