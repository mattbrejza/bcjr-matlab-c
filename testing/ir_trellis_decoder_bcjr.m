% Copyright (C) 2013  Robert G. Maunder

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

% The GNU General Public License can be seen at
% http://www.gnu.org/licenses/.



% The soft-input soft-output trellis decoder of Section IV-A in (Maunder et al., 2013)
% ztildea is the LLR matrix comprising n times b a priori LLRs pertaining to z, where LLR=ln[Pr(bit=1)/Pr(bit=0)]
% C is the set comprising r/2 codewords, each having n bits - e.g. C = [0 1; 1 1] in Fig 3(b) of (Maunder et al., 2013)
% Px_vector is a vector comprising the probabilities of the first r/2-1 symbol values
% l is the average unary codeword length
% a is the number of symbols to decode
% ztildee is the LLR matrix comprising n times b extrinsic LLRs pertaining to z, where LLR=ln[Pr(bit=1)/Pr(bit=0)]
% ytildep is the LLR vector comprising b a posteroiri LLRs pertaining to y, where LLR=ln[Pr(bit=1)/Pr(bit=0)]
function [ztildee_cell, ytildep] = ir_trellis_decoder_bcjr(ztildea_cell, codewords_cell, p_cell, l, a, frac)

% All calculations are performed in the logarithmic domain in order to
% avoid numerical issues. These occur in the normal domain, because some of
% the confidences can get smaller than the smallest number the computer can
% store.
%
% A multiplication of two confidences is achieved using the addition of the
% corresponding log-confidences. If A = log(a) and B = log(b), then
% log(a*b) = A+B.
%
% An addition of two confidences is achieved using the Jacobian logarithm
% of the corresponding log-confidences. The Jacobian logarithm is defined
% in the jac.m file. If A = log(a) and B = log(b), then
% log(a+b) = max(A,B) + log(1+exp(-abs(A-B))).

% maximum of the states, number of rows in the big trellis
r_max = 0;
% number of colmuns in the big trellis
b = 0;

for i_trellis = 1:1:length(frac-1)
    
    % if length(Px_vector) ~= size(C,1)-1
    %     error('length(Px_vector) should equal size(C,1)-1');
    % end
    
    if size(ztildea_cell{1, i_trellis},1) ~= size(codewords_cell{1, i_trellis},2)
        error('size(ztildea,1) should equal size(C,2)');
    end
    
    % Determine number of states
    r_max = max(r_max, 2*size(codewords_cell{1, i_trellis},1));
    % Determine length of unary-encoded bit sequence
    b = b + size(ztildea_cell{1, i_trellis},2);    
end

% % Determine number of states
% r = 2*size(C,1);
%
% % Determine number of codeword bits
% n = size(C,2);

% Determine length of unary-encoded bit sequence
% b = size(ztildea,2);

% Build the trellis and calculate conditional transition probabilities
% transitions = zeros(2*r, n+3);
gammas = zeros(2*r_max, b);

if sum(frac) < 0.99999
    error('sum(frac) should be one')
end

% calculate the index of each fractions
b_frac = round(b * cumsum(length(b) * frac));
b_frac = [0 b_frac];

for i_trellis = 1:1:(length(b_frac) - 1)
    % number of states of current codewords
    r = 2*size(codewords_cell{1, i_trellis}, 1);
    %
    Px_vector = p_cell{1, i_trellis};
    
    for mprime = 1:r
        
        odd = mod(mprime,2);
        
        for y = 0:1
            
            if y == 0
                m = 1 + odd;
            else
                m = min(mprime+2,r - odd);
            end
            
            if mprime <= r-2
                if m == mprime + 2
                    Pmmprime = (1-sum(Px_vector(1:ceil(mprime/2))))/(1-sum(Px_vector(1:ceil(mprime/2)-1)));
                elseif m == 1+odd
                    Pmmprime = Px_vector(ceil(mprime/2))/(1-sum(Px_vector(1:ceil(mprime/2)-1)));
                else
                    Pmmprime = 0;
                end
            else
                if m == 1 + odd
                    Pmmprime = (1-sum(Px_vector(1:r/2-1)))/(1+l-r/2-sum(Px_vector(1:r/2-1).*(1+(1:r/2-1)-r/2)));
                elseif m == mprime
                    Pmmprime = (l-r/2-sum(Px_vector(1:r/2-1).*((1:r/2-1)-r/2)))/(1+l-r/2-sum(Px_vector(1:r/2-1).*(1+(1:r/2-1)-r/2)));
                else
                    Pmmprime = 0;
                end
            end
            %             gammas(y*r+mprime,:) = log(Pmmprime);
            gammas(y*r+mprime,(b_frac(i_trellis)+1:b_frac(i_trellis+1))) = log(Pmmprime);
        end
    end
end

% Calculate a priori transition log-probabilities
for i_trellis = 1:1:(length(b_frac) - 1)
    transitions = trans_cur(codewords_cell{1, i_trellis});
    n = size(codewords_cell{1, i_trellis},2);
    for codebit_index = 1:n
        gammas(transitions(:,3+codebit_index)==1,(b_frac(i_trellis)+1:b_frac(i_trellis+1))) = gammas(transitions(:,3+codebit_index)==1,(b_frac(i_trellis)+1:b_frac(i_trellis+1))) + repmat(ztildea_cell{1, i_trellis}(codebit_index, :), sum(transitions(:,3+codebit_index)==1),1);
    end
end

% Forward recursion to calculate state log-probabilities
b_frac(1) = 1;
alphas=-inf(r_max,b);
alphas(1,1)=0; % We know that this is the first state
for i_trellis = 1:1:(length(b_frac) - 1)
    transitions = trans_cur(codewords_cell{1, i_trellis});
    sb = (b_frac(i_trellis) + 1);
    eb = b_frac(i_trellis + 1);
    for j = sb:1:eb
        if (j == sb) && (i_trellis > 1) && (j>1)
            transitions_j = trans_jnt(codewords_cell{1, i_trellis-1}, codewords_cell{1, i_trellis});
            temp = alphas(transitions_j(:,1),j-1)+gammas((1:size(transitions_j,1)),j-1);
            for m = 1:max(transitions_j(:,2)) %1:2*size(codewords_cell{1, i_trellis+1},1)   %not sure the +1 should really be here
                alphas(m,j) = maxstar_rob(temp(transitions_j(:,2) == m));
            end

        else
            temp = alphas(transitions(:,1),j-1)+gammas((1:size(transitions,1)),j-1);
            for m = 1:2*size(codewords_cell{1, i_trellis},1)
                alphas(m,j) = maxstar_rob(temp(transitions(:,2) == m));
            end
        end
    end
end

% Backward recursion to calculate state log-probabilities
betas=-inf(r_max,b);
betas(1 + mod(a,2),end)=0; % We know that this is the final state
for i_trellis = (length(b_frac) - 1):-1:1
    transitions = trans_cur(codewords_cell{1, i_trellis});
    sb = b_frac(i_trellis);
    eb = (b_frac(i_trellis+1) - 1);
    for j = eb:-1:sb
        if (j == eb) && (i_trellis < (length(b_frac) - 1)) && (j < b)           
            transitions_j = trans_jnt(codewords_cell{1, i_trellis}, codewords_cell{1, i_trellis+1});
            temp = betas(transitions_j(:,2),j+1)+gammas((1:size(transitions_j,1)),j+1);
            for mprime = 1:max(transitions_j(:,1)) %1:2*size(codewords_cell{1, i_trellis-1},1)
                betas(mprime,j) = maxstar_rob(temp(transitions_j(:,1) == mprime));
            end

        else
            temp = betas(transitions(:,2),j+1)+gammas((1:size(transitions,1)),j+1);
            for mprime = 1:2*size(codewords_cell{1, i_trellis},1)
                betas(mprime,j) = maxstar_rob(temp(transitions(:,1) == mprime));
            end
        end
    end
end

% Calculate a posteriori transition log-probabilities
b_frac(1) = 0;
deltas = zeros(size(gammas));
for i_trellis = 1:1:(length(b_frac)-1)
    transitions = trans_cur(codewords_cell{1, i_trellis});
    deltas((1:size(transitions,1)),(b_frac(i_trellis)+1:(b_frac(i_trellis+1)-1))) = alphas(transitions(:,1),(b_frac(i_trellis)+1:(b_frac(i_trellis+1)-1))) + betas(transitions(:,2),(b_frac(i_trellis)+1:(b_frac(i_trellis+1)-1))) + gammas((1:size(transitions,1)),(b_frac(i_trellis)+1:(b_frac(i_trellis+1)-1)));
    % for the joint part
    if i_trellis < (length(b_frac)-1)
        transitions = trans_jnt(codewords_cell{1, i_trellis}, codewords_cell{1, i_trellis+1});
        deltas((1:size(transitions,1)),b_frac(i_trellis+1)) = alphas(transitions(:,1),b_frac(i_trellis + 1)) + betas(transitions(:,2),b_frac(i_trellis + 1)) + gammas((1:size(transitions,1)),b_frac(i_trellis + 1));
    else
        deltas((1:size(transitions,1)),end) = alphas(transitions(:,1),end) + betas(transitions(:,2),end) + gammas((1:size(transitions,1)),end);
    end    
end

% Calculate a posteriori LLRs pertaining to y and z
ztildee_cell = cell(1, (length(b_frac)-1));
ytildep = [];
for i_trellis = 1:1:(length(b_frac)-1)
    transitions = trans_cur(codewords_cell{1, i_trellis});
    yztildep = zeros(size(codewords_cell{1, i_trellis},2)+1, (b_frac(i_trellis+1) - b_frac(i_trellis)));
    for codebit_index = 1:(size(codewords_cell{1, i_trellis},2) + 1)
        log_p1=maxstar_rob(deltas(transitions(:,2+codebit_index) == 1,(b_frac(i_trellis)+1):(b_frac(i_trellis+1))));
        log_p0=maxstar_rob(deltas(transitions(:,2+codebit_index) == 0,(b_frac(i_trellis)+1):(b_frac(i_trellis+1))));
        yztildep(codebit_index,:) = log_p1-log_p0;
    end
    ytildep = [ytildep yztildep(1, :)];
    ztildee_cell{1, i_trellis} = yztildep(2:end,:) - ztildea_cell{1, i_trellis};
end

end

function transitions = trans_cur(C)

% C = [0 1; 1 1; 1 1]
% C_nxt = [1; 1; 0]

% Determine number of states
r = 2 * size(C,1);

% Determine number of codeword bits
n = size(C,2);

% Build the trellis and calculate conditional transition probabilities
transitions = zeros(2*r, n+3);

for mprime = 1:r
    odd = mod(mprime,2);
    
    for y = 0:1
        
        % Set the from state
        transitions(y*r+mprime,1) = mprime;
        
        % Set the to state
        if y == 0
            m = 1 + odd;
        else
            m = min(mprime+2,r - odd);
        end
        transitions(y*r+mprime,2) = m;
        
        % Set the uncoded bit
        transitions(y*r+mprime,3) = y;
        
        % Set the codeword
        if y == odd
            transitions(y*r+mprime,4:end) = ~C(ceil(mprime/2), :);
        else
            transitions(y*r+mprime,4:end) = C(ceil(mprime/2), :);
        end
    end
end
end


function transitions_joint = trans_jnt(C_cur, C_nxt)

% C_cur = [0 1; 1 1; 1 1]
% C_nxt = [1]

% Determine number of states
r_cur = 2 * size(C_cur,1);
r_nxt = 2 * size(C_nxt,1);

% Determine number of codeword bits
n_cur = size(C_cur,2);
% n_nxt = size(C_nxt,2);

% Build the trellis and calculate conditional transition probabilities
transitions_joint = zeros(2*r_cur, n_cur+3);

for mprime = 1:r_cur
    odd = mod(mprime,2);
    
    for y = 0:1
        
        % Set the from state
        transitions_joint(y*r_cur+mprime,1) = mprime;
        
        % Set the to state
        if y == 0
            m = 1 + odd;
        else
            m = min(mprime+2,r_nxt - odd);
        end
        transitions_joint(y*r_cur+mprime,2) = m;
        
        % Set the uncoded bit
        transitions_joint(y*r_cur+mprime,3) = y;
        
        % Set the codeword
        if y == odd
            transitions_joint(y*r_cur+mprime,4:end) = ~C_cur(ceil(mprime/2), :);
        else
            transitions_joint(y*r_cur+mprime,4:end) = C_cur(ceil(mprime/2), :);
        end
    end
end
%transitions_joint = C_cur; %remove
end




function out = maxstar_rob(in)


    if size(in,1) == 1
        in = in';
    end
    out = in(1,:);
    for index = 2:size(in,1)        
        difference = out-in(index,:);
        difference(isnan(difference)) = 0;
        out = max(out,in(index,:)) + log(1+exp(-abs(difference))); % Log-MAP
    end
end
   

    
