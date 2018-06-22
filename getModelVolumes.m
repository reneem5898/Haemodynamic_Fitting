function V = getModelVolumes(directory, numFrames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function gets the endocardial volumes from the CIM model
% 
% Input: directory - model directory
% Output: Vector of volumes (over cardiac cycle)
%
% Written by: Jenny Wang
% Modified by: Renee Miller (rmil520@aucklanduni.ac.nz)
% Last modified: 14 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volume directory
volume_dir = sprintf('%s/volumes/info/', directory);

% Get volume file
% Choose *_ca* file from the CIM model that you want to use
file_spec = sprintf('%s/*_ca*', volume_dir);
fname = getlatestfile(file_spec);

% Print to log file
disp('CIM model volume file:')
fprintf('%s/%s', volume_dir, fname);

% Read file
fid = fopen(sprintf('%s/%s', volume_dir, fname), 'r');
i = 1; % counter

% Read in junk lines - HARD CODED FOR THE FORMAT OF CIM MODEL FILE: <patient name>_model_ca.model_<model name>
junk = fgets(fid);
while ~contains(junk, 'TOTAL') % loop until the word 'TOTAL' is found
    junk = fgets(fid);
    i = i + 1;
end
for i = 1:2 % get two more lines
    junk = fgets(fid);
end
junk = fscanf(fid, '%s', 2);

% Finally, read in volumes
V = fscanf(fid, '%f', numFrames);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Section of file to read should look similar to the text below (without the %'s): 
% Lines 169 - 177 are shown below.
% V = vector of endocardial volumes (in LV chamber) at each cardiac time point 
% 
%  TOTAL
%     {     /* Cross are phases starting from phase 1 (time 0) */
%       Epi  :  142.5387  142.2661  140.9335  137.0897  132.9766  129.0236  124.9715  120.7311  116.6967  113.1890  109.8485  107.3701  106.9944  108.5053  109.8687  112.8330  121.2663  129.3819  135.1139  136.8398  136.9601  136.8427  136.7795  136.1650  136.1787  136.0120  136.1216  137.1403  138.9960  141.3227
%       Endo :   73.4868   72.9368   71.5210   68.3716   64.7001   60.8805   56.9549   52.9616   49.1178   45.5650   42.2919   39.8493   39.0222   39.9971   41.9877   45.6593   52.6131   59.9127   65.7352   68.4765   69.2054   68.9980   68.5213   67.8736   67.8660   68.0508   68.5157   69.4444   70.7671   72.3207
%       Wall :   69.0518   69.3292   69.4125   68.7181   68.2765   68.1431   68.0166   67.7695   67.5789   67.6240   67.5566   67.5208   67.9722   68.5082   67.8809   67.1736   68.6533   69.4692   69.3787   68.3633   67.7547   67.8447   68.2582   68.2914   68.3127   67.9612   67.6059   67.6959   68.2289   69.0020
%       WT   :    6.6168    6.6608    6.7340    6.8479    7.0137    7.2230    7.4587    7.7198    8.0021    8.3075    8.5896    8.8133    8.9285    8.8784    8.6495    8.2854    7.8761    7.4467    7.0928    6.8738    6.7944    6.8084    6.8538    6.8778    6.8667    6.8204    6.7624    6.7112    6.6736    6.6436
%       MinWT:    3.8012    3.8796    3.9525    4.0442    4.1666    4.3416    4.5784    4.8834    5.2174    5.5658    5.8516    6.0968    6.2870    6.3732    6.3087    6.1360    5.8257    5.4776    5.1285    4.8604    4.6549    4.5197    4.4419    4.3959    4.3380    4.2940    4.2455    4.1861    4.1282    4.0611
%       MaxWT:    8.5257    8.5085    8.5337    8.6424    8.8231    9.0383    9.2703    9.5108    9.7696   10.0654   10.3341   10.5424   10.6270   10.5127   10.1950    9.7392    9.2418    8.8123    8.4809    8.2943    8.2576    8.3191    8.4175    8.4986    8.5420    8.5464    8.5294    8.5087    8.5038    8.4943
%     }
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %