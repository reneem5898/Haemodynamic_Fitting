function [RR, numFrames] = getRR_CIM(directory)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads the RR interval from the CIM spec file
%
% Input: directory - CIM model directory
%
% Output: RR interval
%
% Written by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 14 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get model name
tmp = strsplit(directory, '/');

% Load relevant file
f = sprintf('%s/system/%s.spec', directory, tmp{end});

% Open file
fid = fopen(f);

% Get lines until 'RR interval' is found
while ~feof(fid)
   
    tline = fgetl(fid); % Get next line
   
   if strfind(tline, 'RR') % When 'RR' is found in line
       tmp = strsplit(tline, ':'); % Split line at :
       RR = str2double(tmp{end}); % Get the RR interval
   
   elseif strfind(tline, 'Series 0 number of times') 
       tmp = strsplit(tline, ':'); % Split line at :
       numFrames = str2num(tmp{end}); % Get the number of frames in the short-axis series
   end
   
end


 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% File to read should look similar to the text below (without the %'s): 

% #CIM V8.1.6 study spec file
% 
% Study name        : Pig1_ET38_HFPEF_Base
% Problem type      : LV_GEOM_BP
% Number of series  : 2
% Apex is slice one : true
% Number of models  : 2
% Focal length      :   29.26407
% RR interval       :  780.00000                     <-- RELEVANT LINE
% Max frame time    :  754.96667
% User defined heart rate : -1
% User defined height     : -1.0000
% User defined weight     : -1.0000
% 
% Series 0 number of slices : 12
% Series 0 number of times  : 30                     <-- RELEVANT LINE
% Series 0 label            : Short Axis            
% 
% Series 1 number of slices : 4
% Series 1 number of times  : 30
% Series 1 label            : Long Axis

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

