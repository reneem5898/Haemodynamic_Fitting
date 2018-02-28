function ES = getES_dicroticNotch(AOP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function finds ES - dicrotic notch in AOP
%
% Input: AOP = aortic pressure (vector, 1 cycle)
%
% Output: ES = index of dicrotic notch in AOP
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 14 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE TO USER: This may need to be altered depending on data. In the Ohio 
% pig data, the dicrotic notch is quite prominent and could be found using 
% findpeaks.m However, in human data, the dicrotic notch may not be a 
% minimum and must be found another way. Feel free to make another version
% or edit this function to find the dicrotic notch in your data.

% Find minima in AOP
[~, l] = findpeaks(-AOP);

% Make a pseudo-time vector - arbitrary
pt = linspace(1, length(AOP), length(AOP));
dataPointLabels = num2cell(linspace(1,length(l), length(l)));

% Plot cycle with minima identified
FH = figure('position', [100 100 600 600]);
plot(pt, AOP, '-', 'markersize', 2);
hold on
scatter(l, AOP(l), 'ro', 'filled')
title('AOP')
xlabel('Time')
ylabel('Pressure (mmHg)')
for i = 1:length(l)
    text(l(i), AOP(l(i)) - 0.2, dataPointLabels(i), 'FontSize', 14)
end

% Create a list of minima found by findpeaks.m
minima = [linspace(1,length(l),length(l)), NaN];
minima = strtrim(cellstr(num2str(minima'))');

% Use dialogue box for selection of dicrotic notch - BETTER WAY TO DO THIS?? I NEED COFFEE...
tmp = listdlg('PromptString','Choose the dicrotic notch',...
    'SelectionMode','single','ListString',minima);

% If user chose NaN - Dicrotic notch was not identified as a local minimum
if tmp == length(minima)
    waitfor(msgbox('Select the approximate location of the dicrotic notch and press enter.', 'Dicrotic notch'));
    [x,~] = getpts(FH); % Get user to choose approximate location of dicrotic notch
    ES = round(x);
else
    ES = l(tmp); % Else, if user chose a minimum on the list, save that location
end

% Close figure
close(FH);