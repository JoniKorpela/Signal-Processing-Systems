
% Run the simulink model with params from T2
D = 39;
L = 37;
step_size = 0.04867032;
[~, ~, ~, ~] = exsim_run(L, D, step_size);

bits = string_data;
bits = bits(:)';

% Remove training pattern (for some reason there is a 4 element offset,
% maybe filtering delay?)
msg_bits = bits(10005:end);

% reshape into 7-bit chunks
n_chars = floor(length(msg_bits)/7);
msg_bits = msg_bits(1:n_chars*7); 
msg_matrix = reshape(msg_bits, 7, n_chars)'; % each row = 1 character

% Convert to decimal ASCII codes
ascii_codes = bi2de(msg_matrix, 'left-msb');

% Convert to string
decoded_msg = char(ascii_codes);
disp(decoded_msg')

% On the first row, there should read "correlate"
