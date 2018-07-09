function plots_2points(timepoints_MRI, mri_time, LVP, ed2ds, ds2ed, V, directory)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates plots of resulting interpolated pressures
%
% Inputs: 1) timepoints_MRI: vector of cardiac timepoints from MRI
%         2) mri_time: vector of time of one cardiac cycle
%         3) LVP: final averaged left ventricular pressure trace
%         4) ed2ds: LV presure between ED and ES
%         6) ds2ed: LV pressure between DS and ED
%         7) V: volumes (from CIM model)
%         8) directory: directory where to save plots
%
% Plots: 1) beat averaged left ventricular pressure
%        2) pressure volume loop
%        3) volume curve with highlighted points
%
% Written by: Renee Miller (renee.miller@auckland.ac.nz)
% Date: 22 June 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MRI cardiac time points
edMRI = timepoints_MRI(1);
dsMRI = timepoints_MRI(2);


%%%%%%%%%%%%%%%% Beat-averaged left ventricular pressure %%%%%%%%%%%%%%%%%%
FH8 = figure;

% Plot ED to DS
t = repmat(mri_time(edMRI:dsMRI-1),size(ed2ds,1),1);
disp(length(t(:)))
disp(length(ed2ds(:)))
scatter(t(:), ed2ds(:), 'g*')
hold on

% Plot DS to ED
t = repmat(mri_time(dsMRI:end),size(ds2ed,1),1);
scatter(t(:), ds2ed(:), 'b*')

% Plot average
plot(mri_time, LVP, 'k*-')
xlabel('Time (s)', 'FontSize', 16)
ylabel('LV Pressure (kPa)', 'FontSize', 16)
title('Beat-averaged LV Pressure', 'FontSize', 16)
set(gca, 'FontSize', 12)
saveas(FH8, sprintf('%s/beat-averaged-lv-pressure.png', directory));


%%%%%%%%%%%%%%%%%%% Pressure-volume loop %%%%%%%%%%%%%%%%%%
FH9 = figure;

% Plot ED to DS
t = repmat(V(edMRI:dsMRI-1)',size(ed2ds,1),1);
scatter(t(:), ed2ds(:), 'g*')
hold on

% Plot DS to ED
t = repmat(V(dsMRI:end)',size(ds2ed,1),1);
scatter(t(:), ds2ed(:), 'b*')

% Plot Pressure-Volume
plot(V, LVP, 'k*-')
xlabel('Volume (mL)', 'FontSize', 16)
ylabel('LV Pressure (kPa)', 'FontSize', 16)
title('Pressure-Volume Loop', 'FontSize', 16)
set(gca, 'FontSize', 12)
saveas(FH9, sprintf('%s/pv-loop.png', directory));


%%%%%%%%%%%%%%%%%%% Volume with cardiac events %%%%%%%%%%%%%%%%%%
FH10 = figure;
plot(mri_time, V, 'ksq-', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
hold on
scatter(mri_time(dsMRI), V(dsMRI), 50, 'ro', 'filled')
scatter(mri_time(edMRI), V(edMRI), 50, 'ro', 'filled')
xlabel('Time (s)', 'FontSize', 16)
ylabel('Volume (mL)', 'FontSize', 16)
set(gca, 'FontSize', 12)
saveas(FH10, sprintf('%s/volume-mri-cardiac-events.png', directory));