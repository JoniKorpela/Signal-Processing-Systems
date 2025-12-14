% Assume bitstream is a vector of 0/1 called 'bits'
bits = string_data;  % from Simulink output
bits = bits(:)'; % make it row vector

% Remove training pattern
msg_bits = bits(10005:end);

% Reshape into 7-bit chunks
n_chars = floor(length(msg_bits)/7);
msg_bits = msg_bits(1:n_chars*7);  % truncate excess bits
msg_matrix = reshape(msg_bits, 7, n_chars)'; % each row = 1 character

% Convert to decimal ASCII codes
ascii_codes = bi2de(msg_matrix, 'left-msb'); % left-msb matches standard ASCII

% Convert to string
decoded_msg = char(ascii_codes);
disp(decoded_msg')