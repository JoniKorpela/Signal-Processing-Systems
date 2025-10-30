function [input_data, calibration_data, layer_0_weights, layer_0_bias, output_weights, output_bias] = generate_data_and_network_parameters(data_dimension, N_data_points, layer_0_dimension, output_dimension)

data_multiplier = 5;
data_offset = -2.2;
 
input_data = rand(data_dimension,N_data_points)*data_multiplier + data_offset;
calibration_data = rand(data_dimension,N_data_points)*data_multiplier + data_offset;

layer_0_weights = rand(layer_0_dimension, data_dimension)* 2 - 1.0;
layer_0_bias = rand(layer_0_dimension, 1) - 0.6;

output_weights = rand(output_dimension, layer_0_dimension) * 3 - 2;
output_bias = rand(output_dimension, 1) * 2 - 1.0;
