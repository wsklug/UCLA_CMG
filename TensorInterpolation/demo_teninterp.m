% DT-MRI tensor field interpolation script example for using the tensor 
% interpolation functions
%
% SYNTAX:  demo_interp;
%
% INPUTS:  None.
%
% OUTPUTS: None
%
% Written by JK Gahm @ UCLA on 03/18/2014.
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

%% load and validate data
clear

% load mesh nodes (x,y,z): nx3 (spatial resolution: 1x1x1 cm)
load ./mesh_nodes

% load DT-MRI data (spatial resolution: 0.5x0.5x0.5 mm)
load ./tensor_data % normal rabbit heart

% match mesh nodes and tensor data (resolution and orientation)
data1 = 10/0.5 * data; % ***match resolutions

% % validate data
% m = false(size(MASK1));
% for i=1:size(data1,1)edit 
%   m(round(data1(i,1)),round(data1(i,2)),round(data1(i,3))) = true;
% end
% 
% % roughly make sure they match
% figure, imagescn(m)
% figure, imagescn(MASK1)

%% tensor field interpolation
Y = teninterp(Dspd,data1,'DY'); % Dyadic tensor-based interpolation
