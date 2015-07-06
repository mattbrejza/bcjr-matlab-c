% Wenbo, 2013.08.2
% This encoder is based on Rob's new trellis_encoder that on his personal website


function [z_cell, elmt_vec] = ir_trellis_encoder(y,C_cell, frac)

% calculate the index of each fractions
y_frac = round(cumsum(length(y) * frac));
y_frac = [0 y_frac];

% Always start from a previous state of 1
mprime = 1;

z_cell = cell(1, length(frac));
elmt_vec = zeros(1, length(frac));
for i_trellis = 1:1:length(y_frac) - 1
    
    % select the current trellis
    C = C_cell{1, i_trellis};
    
    % Determine the number of states
    r = 2*size(C,1);
    
    % Perform the trellis encoding
    z = inf(size(C,2), y_frac(i_trellis+1)-y_frac(i_trellis));
    
    for j = (y_frac(i_trellis) + 1):1:y_frac(i_trellis + 1)
        
        odd = mod(mprime,2);
        
        % in the trellis
        if j < y_frac(i_trellis + 1)
            
            % Determine the next state
            if y(j) == 0
                m = 1 + odd;
            else
                m = min(mprime+2,r - odd);
            end
            
            % on the edge of the joint trellis
        else
            if y(j) == 0
                m = 1 + odd;
            else
                if i_trellis < length(y_frac) - 1
                    m = min(mprime+2,2*size(C_cell{1, i_trellis+1}, 1) - odd);
                else
                    m = min(mprime+2,r - odd);
                end
            end
        end
        
        % Select the corresponding codeword
        if y(j) == odd
            z(:,j - y_frac(i_trellis)) = ~C(ceil(mprime/2), :);
        else
            z(:,j - y_frac(i_trellis)) = C(ceil(mprime/2), :);
        end
        
        % Enter the next state
        mprime = m;
    end
    z_cell{i_trellis} = z;
    elmt_vec(i_trellis) = numel(z);
end
end


