function [ rawPoints ] = readCloudCsv( filename, subSampleRatio )
%READCLOUDCSV read csv files of ETH data set
% INPUTS:
% - filename: CSV file name of the point cloud to be ploted. It can be
%                   local point cloud (ex.: Hokuyo_0.csv) or global map (ex.:
%                   PointCloud0.csv).
% - subSampleRatio: Ratio of point to be removed (used for performance
%                               reasons)
%

% Load data
cloud = importdata(filename);

% Extract coordinates
x = cloud.data(:, strcmp('x', cloud.colheaders));
y = cloud.data(:, strcmp('y', cloud.colheaders));
z = cloud.data(:, strcmp('z', cloud.colheaders));

% Subsample
f = rand(1, length(x));

x = x(f > subSampleRatio);
y = y(f > subSampleRatio);
z = z(f > subSampleRatio);

rawPoints = [x y z];
end

