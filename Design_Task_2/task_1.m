clc;
clear;
% Group 2 Parameters for task1.
c = -2.851;
d = 3.402;
p_c = 7;
n_c = 2;
p_d = 8;
n_d = 3;
p_o = 12;
n_o = 3;

F = fimath('RoundingMethod', 'zero');
F.ProductMode = 'SpecifyPrecision';
F.ProductWordLength = p_o;
F.ProductFractionLength = n_o;
F.SumMode = 'SpecifyPrecision';
F.SumWordLength = p_o;
F.SumFractionLength = n_o;

% fi uses round to nearest by default for inputs
c_fixed = fi(c, 1, p_c, n_c, F);
d_fixed = fi(d, 1, p_d, n_d, F);

intermediate_1 = c_fixed * d_fixed;
intermediate_1
intermediate_2 = intermediate_1 + c_fixed;
intermediate_2
final_result = intermediate_2 - d_fixed;
final_result
bin(final_result)


