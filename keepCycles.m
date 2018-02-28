function cyclesToKeep = keepCycles(pressure, R_positions, fig_filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots the pressure trace with cycles identified. The user 
% must select which cycles SHOULD BE KEPT. The remaining cycles will omitted
% in subsequent analyses.
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 19 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get max LVP for putting text on graph
maxLVP = max(pressure);
yText = maxLVP + 0.2 * maxLVP;

% Get x text positions - approximate centre of each cycle
xText = zeros(length(R_positions)-1,1);
for i = 1: length(R_positions)-1
    xText(i) = mean(R_positions(i:i+1));
end

% Plot filtered LVP with cycle numbers shown
FH15 = figure('position', [100, 100, 1800, 600]);
plot(pressure, '-', 'LineWidth', 3);
hold on
x = xlim;
plot(linspace(x(1), x(2), 100), zeros(100,1), 'k--', 'LineWidth', 2)
ylabel('Pressure (kPa)', 'FontSize', 14);
xlabel('Time', 'FontSize', 14);
y = ylim;
for i = 1:length(R_positions)
    plot(R_positions(i)*ones(100,1), linspace(y(1), y(2), 100), 'r--')
end
for i = 1:length(xText)
    text(xText(i), yText, num2str(i), 'FontSize', 12)
end
ylim([y(1) y(2) + 0.2*y(2)]) % hack to get the ylimits to show the labels % I need coffee

% Save filename
saveas(FH15, fig_filename);

% Directions for the user
waitfor(msgbox('Select cycles to keep in subsequent analysis.', 'LVP'));

% Create a list of cycles
cycleList = linspace(1,length(R_positions)-1,length(R_positions)-1);
cycleListStr = strtrim(cellstr(num2str(cycleList'))');

% Use dialogue box for selection of cycles to omit
cyclesToKeep = listdlg('PromptString','Cycles to keep',...
        'SelectionMode','multiple','ListString',cycleListStr);

% Close figure
close(FH15);