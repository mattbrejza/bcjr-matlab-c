function test_ir()

addpath('testing')

% fractions
frac_urc = [0 0 0 0 0.6094 0.0725 0.0226 0 0.2955 0];
frac_uec = [0.25 0 0 0.50 0 0.25 0]; %[0.9191 0 0.0415 0.0118 0 0 0.0276];

% indexs of codewords in codebooks
codes_ind = find(frac_uec);
codebooks = cell(1,7);
codebooks{1, 1} = [0 0; 0 0]%; 0 0; 0 0];%;0 0];%
codebooks{1, 2} = [0 0; 0 1; 0 1; 0 1; 0 1];
codebooks{1, 3} = [0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0];%
codebooks{1, 4} = [0 0 0; 0 1 1; 0 1 1; 0 1 1; 0 1 1];%
codebooks{1, 5} = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
codebooks{1, 6} = [0 0 0 0; 0 1 1 1; 0 1 1 1; 0 1 1 1];% 0 1 1 1];%
codebooks{1, 7} = [0; 1; 1; 1; 1];
codewords_cell = cell(1, length(codes_ind));
for i = 1:1:length(codes_ind)
    codewords_cell{1, i} = codebooks{1, codes_ind(i)};
end

frac_uec(frac_uec == 0) = [];
% p1 = 0.797
s = 2.7705;
l = 1.5368;
a = 100;
distribution = 'distribution_zeta';
p_cell = cell(1, length(frac_uec));
for i_trellis = 1:1:length(frac_uec)
    p = feval(distribution,1:size(codewords_cell{1, i_trellis},1)-1,s);
    p_cell{1, i_trellis} = p;
end

while 1
% Generate some random symbols
x = generate_random_symbols(a, distribution, s);
%x = [1 1 2 1 1 4 1 1 2 1];
y = unary_encoder(x);
%y=randn([1 10])>0;

% Irregular trellis encoder
[z_cell, z_num_vec] = ir_trellis_encoder(y, codewords_cell, frac_uec);
z_vec = cell2vec(z_cell, codewords_cell);
z_tilde_a = zeros(size(z_vec));
z_tilde_a(z_vec>0) = 1.5;
z_tilde_a(z_vec<=0) = -1.5;
z_tilde_a_cell = vec2cell(z_tilde_a, codewords_cell, z_num_vec);


cell_trans = cell(size(codewords_cell));
for ct = 1:numel(codewords_cell)
   cell_trans{ct} =  gen_uec_transitions(codewords_cell{ct}, p_cell{ct}, l );
end
% Irregular trellis decoder
z_tilde_a_v = bits_cell_to_vec(z_tilde_a_cell);
frame_len = size(z_tilde_a_v,2);
tic;
[z_tilde_e_cell, ytildep] = ir_trellis_decoder_bcjr(z_tilde_a_cell, codewords_cell, p_cell, l, sum(y==0), frac_uec);
matlab_run_time = toc
tic
[y_tilde_p_c, ztilde_e_c] = bcjr_decoder_iruec(zeros([1 frame_len]), z_tilde_a_v, 'iruec', cell_trans, frac_uec, mod(a,2)+1);
c_run_time = toc
speedup = matlab_run_time / c_run_time
%%%% compare c and matlab

y_tilde_p_c(y_tilde_p_c < -7000) = -inf;
ztilde_e_c(ztilde_e_c < -7000) = -inf;
ztilde_e_c(ztilde_e_c > 7000) = inf;
ztilde_e_c_cell = bits_vec_to_cell(ztilde_e_c, z_tilde_e_cell);


if ~cell_equal(ztilde_e_c_cell,z_tilde_e_cell)
    error('not equal')
end
if sum(abs(ytildep(1:end-1)-y_tilde_p_c(1:end-1))) > 0.001
   error('not equal')
end


% %%
ytildep;
y_hat = ytildep>0;
y_error = sum(y ~= y_hat);
    x_hat = unary_decoder_soft(ytildep, a);
x_error = sum(x ~= x_hat);
x_error;
y_error;




y_hat_c = y_tilde_p_c>0;
y_error_c = sum(y ~= y_hat_c);
    x_hat_c = unary_decoder_soft(y_tilde_p_c, a);
x_error_c = sum(x ~= x_hat_c);
x_error_c;
y_error_c;
end

end

function c_vec = cell2vec(c_cell, codewords_cell)

if(length(c_cell) ~= length(codewords_cell))
    error('Must have same number of c_cell and codewords_cell!');
end

c_vec = [];

for i_trellis = 1:1:length(c_cell)
    
    c = reshape(c_cell{1, i_trellis}, 1, []);
    
    c_vec = [c_vec c];
    
end
end

function c_cell = vec2cell(c_vec, codewords_cell, z_num_vec)

% initialize the output c_cell
c_cell = cell(1, length(codewords_cell));

% identify where to start and end
c_start = 0;
c_end   = z_num_vec(1);

% loop for each trellis
for i_trellis = 1:1:length(z_num_vec)
    
    %
    c = reshape(c_vec((c_start+1):c_end), size(codewords_cell{1, i_trellis}, 2), []);
    %
    c_cell{1, i_trellis} = c;
        
    % only need to move the positions of c_end & c_start when i_trellis < length(z_num_vec)
    if i_trellis < length(z_num_vec)
        c_start = c_end;
        c_end   = c_end + z_num_vec(i_trellis+1);
    end
    
end
end

function out = bits_cell_to_vec(cell_in)

M = 0;
N = 0;
for i = 1:numel(cell_in)
    M = max(M,size(cell_in{i},1));
    N = N + size(cell_in{i},2);
end

out = zeros(M,N);

bc = 1;
for i = 1:numel(cell_in)
    M = size(cell_in{i},1);
    N = size(cell_in{i},2);
   out(1:M,bc:bc+N-1) = cell_in{i};
   bc = bc+N;
end

end

function out = bits_vec_to_cell(bits_in, ref_cell)


out = cell(size(ref_cell));
bc = 1;
for i = 1:numel(ref_cell)
    M = size(ref_cell{i},1);
    N = size(ref_cell{i},2);
   out{i} = bits_in(1:M,bc:bc+N-1);
   bc = bc + N;
end


end

function out = cell_equal(cell1,cell2)
    out = 1;
    if ~isequal(size(cell1),size(cell2))
        out = 0;
        return
    end

    
    for i = 1:numel(cell1)
        cell1{i}(isinf(cell1{i})) = 900000;
        cell2{i}(isinf(cell2{i})) = 900000;
        if sum(sum(abs(cell1{i}-cell2{i}))) > 0.001
            out = 0;
        end
    end
    
end
