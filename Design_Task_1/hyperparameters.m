function [rng_seed, N_data_points, data_dimension, layer_0_dimension, output_dimension] = hyperparameters(group_number)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [rng_seed, N_data_points, data_dimension, layer_0_dimension, output_dimension] = hyperparameters(group_number)
%
% Parameters
%   group_number      Your group number
%
% Return values
%   rng_seed          Use this value to initialize the random number
%                     generator
%   N_data_points     Number of data points in the dataset.
%   data_dimension    Number of dimensions in a data sample.
%   layer_0_dimension Number of neurons in the first fully connected layer.
%   output_dimension  Number of neurons in the output layer.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rng_seed, N_data_points, data_dimension, layer_0_dimension, output_dimension] = get_hyperparameters(group_number);
