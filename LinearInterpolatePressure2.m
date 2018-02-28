function [mri_pressure] = LinearInterpolatePressure2(raw_pressure, numMRIframes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate pressure points at each image frame assuming linear 
% interpolation between pressure data points.
% Author: ZJW
% Date: 20 Feb 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Initialise interpolated pressure and set first and last values
% mri_pressure = zeros(1, numMRIframes);
% mri_pressure(1) = raw_pressure(1);
% mri_pressure(end) = raw_pressure(end);

% xi and xq values 
xi = linspace(1,numMRIframes,length(raw_pressure));
xq = linspace(1,numMRIframes,numMRIframes);

% Linearly interpolate
mri_pressure = interp1(xi, raw_ressure, xq, 'linear');


% p_m_ratio = (numPressurePoints-1)/(numMRIframes-1); % RM - why - 1???
% 
% for m = 2:numMRIframes-1
%     
%     % What does this do?? RM
%     idx_interp = (m-1) * p_m_ratio + 1; % RM - why -1 and +1 ???
%     idx_before = floor(idx_interp);
%     idx_after = ceil(idx_interp);
%     fraction = idx_interp - idx_before;
%     
%     % Weight the contribution of the pressure data point before and
%     % after the interpolation point using linear interpolation.
%     mri_pressure(m) = raw_pressure(idx_before)*(1-fraction) + raw_pressure(idx_after)*fraction;
% end
% 
% end