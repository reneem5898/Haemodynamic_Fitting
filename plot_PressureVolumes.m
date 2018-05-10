function plot_PressureVolumes(p, v, parentDir, descriptor, omit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots a) pressure-volume loops, b) plots from diastasis (DS)
% to end-diastole (ED) and c) plots of only DS and ED
% 
% Inputs: 1) p - pressure (kPa)
%         2) v - volume (mL)
%
% Renee Miller (renee.miller@auckland.ac.nz)
% Date Modified: 10 May 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of cases
numCases = size(p, 1);

% Plot all pressure-volume loops together
FH = figure;
for i = 1:3:numCases
    if ~any(omit == i)
        h1 = plot(v(i,:), p(i,:), 'r-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        hold on
    end
end
for i = 2:3:numCases
    if ~any(omit == i)
        h2 = plot(v(i,:), p(i,:), 'b-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
        hold on
    end
end
for i = 3:3:numCases
    if ~any(omit == i)
        h3 = plot(v(i,:), p(i,:), 'g-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
        hold on
    end
end
xlabel('Volume (mL)', 'FontSize', 12)
ylabel('Pressure (kPa)', 'FontSize', 12)
legend([h1, h2, h3], {'Bx', 'M1', 'M2'}, 'Location', 'Best')
saveas(FH, sprintf('%s/PV-Loops%s.png', parentDir, descriptor));


% Plot all diastolic portions of the PV loops together
FH2 = figure;
for i = 1:3:numCases
    if ~any(omit == i)
        ds = find(p(i,:) == min(p(i,:)));
        h1 = plot(v(i,ds:end), p(i,ds:end), 'r-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        hold on
    end
end
for i = 2:3:numCases
    if ~any(omit == i)
        ds = find(p(i,:) == min(p(i,:)));
        h2 = plot(v(i,ds:end), p(i,ds:end), 'b-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
        hold on
    end
end
for i = 3:3:numCases
    if ~any(omit == i)
        ds = find(p(i,:) == min(p(i,:)));
        h3 = plot(v(i,ds:end), p(i,ds:end), 'g-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
        hold on
    end
end
xlabel('Volume (mL)', 'FontSize', 12)
ylabel('Pressure (kPa)', 'FontSize', 12)
legend([h1, h2, h3], {'Bx', 'M1', 'M2'}, 'Location', 'Best')
saveas(FH2, sprintf('%s/PV-Diastole%s.png', parentDir, descriptor));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot DS to ED
FH4 = figure;
for i = 1:3:numCases
    if ~any(omit == i)
        ds = find(p(i,:) == min(p(i,:)));
        h1 = plot(v(i,[ds,end]), p(i,[ds,end]), 'r-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        hold on
    end
end
for i = 2:3:numCases
    if ~any(omit == i)
        ds = find(p(i,:) == min(p(i,:)));
        h2 = plot(v(i,[ds,end]), p(i,[ds,end]), 'b-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
        hold on
    end
end
for i = 3:3:numCases
    if ~any(omit == i)
        ds = find(p(i,:) == min(p(i,:)));
        h3 = plot(v(i,[ds,end]), p(i,[ds,end]), 'g-o', 'MarkerSize', 3, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
        hold on
    end
end
xlabel('Volume (mL)', 'FontSize', 12)
ylabel('Pressure (kPa)', 'FontSize', 12)
legend([h1, h2, h3], {'Bx', 'M1', 'M2'}, 'Location', 'Best')
saveas(FH4, sprintf('%s/PV-Diastole-DS-ED%s.png', parentDir, descriptor));

close all