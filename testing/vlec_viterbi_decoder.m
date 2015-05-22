


% c_tilde_a is a matrix of a priori encoded LLRs
% b_tilde_a is a vector of a priori uncoded LLRs
% c_tilde_e is a matrix of extrinsic encoded LLRs
% b_tilde_p is a vector of a posteriori uncoded LLRs
function b_hat = vlec_viterbi_decoder(b_tilde_a, transitions)

bit_count = length(b_tilde_a);



% All calculations are performed in the logarithmic domain in order to
% avoid numerical issues. These occur in the normal domain, because some of
% the confidences can get smaller than the smallest number the computer can
% store. See Section 2.2.2 of Liang Li's nine month report for more information
% on this.
%
% A multiplication of two confidences is achieved using the addition of the
% corresponding log-confidences. If A = log(a) and B = log(b), then
% log(a*b) = A+B (Equation 2.12 in Liang Li's nine month report).
%
% An addition of two confidences is achieved using the Jacobian logarithm
% of the corresponding log-confidences. The Jacobian logarithm is defined
% in the jac.m file. If A = log(a) and B = log(b), then
% log(a+b) = max(A,B) + log(1+exp(-abs(A-B))) (Equation 2.12 in Liang Li's
% nine month report).

% Matrix to describe the trellis
% Each row describes one transition in the trellis
% Each state is allocated an index 1,2,3,... Note that this list starts
% from 1 rather than 0.
%%               FromState,  ToState,    UncodedBit, EncodedBit1
%%transitions =  [1,          1,          0,          0;
%%                1,          2,          1,          1;
%%                2,          2,          0,          1;
%%                2,          1,          1,          0];
%% Find the largest state index in the transitions matrix
%% In this example, we have eight states since the code has three memory elements
state_count = max(max(transitions(:,1)),max(transitions(:,2)));
transition_count = size(transitions,1);
codeword_length = size(b_tilde_a,1);

gammas = zeros(transition_count, bit_count);
for transition_index = 1:transition_count
    if transitions(transition_index, 3)
        gammas(transition_index, b_tilde_a ~= inf) = b_tilde_a(b_tilde_a ~= inf); %repmat(b_tilde_a(b_tilde_a ~= inf),[1,1,codeword_length+1]);
    else
        gammas(transition_index, b_tilde_a == inf) = -inf;
    end
end

gammas = gammas + repmat(transitions(:,5), [1 numel(b_tilde_a)]);
%for codebit_index = 1:codeword_length
%    for transition_index = 1:transition_count
%        if transitions(transition_index, 3+codebit_index)
%            gammas(transition_index, c_tilde_a(codebit_index,:) ~= inf,[1:codebit_index, codebit_index+2:codeword_length+1]) = add(gammas(transition_index, c_tilde_a(codebit_index,:) ~= inf,[1:codebit_index, codebit_index+2:codeword_length+1]),repmat(c_tilde_a(codebit_index, c_tilde_a(codebit_index,:) ~= inf),[1,1,codeword_length]));
%        else
%            gammas(transition_index, c_tilde_a(codebit_index,:) == inf,[1:codebit_index, codebit_index+2:codeword_length+1]) = -inf;
%        end
%    end
%end

%gammas

% Forward recursion to calculate state log-confidences. This is similar to
% Equation 1.13 in Rob's thesis or Equations 5 and 6 in the BCJR paper.
alphas=-inf*ones(state_count,bit_count);
alphas(1,1)=0; % We know that this is the first state
for bit_index = 2:bit_count
    set = zeros(1,state_count);
    for transition_index = 1:transition_count
        if(set(transitions(transition_index,2)))
            alphas(transitions(transition_index,2),bit_index) = maxstar(alphas(transitions(transition_index,2),bit_index), add(alphas(transitions(transition_index,1),bit_index-1), gammas(transition_index, bit_index-1,1)));
        else
            alphas(transitions(transition_index,2),bit_index) = add(alphas(transitions(transition_index,1),bit_index-1), gammas(transition_index, bit_index-1));
            set(transitions(transition_index,2)) = 1;
        end
    end
end

%alphas
% Viterbi
state = 1;   %VLEC starts/ends at state 1  %[large_state, state] = max(alphas(:,end));
%best_state = 1;  %why was this commented out?
for bit_index = bit_count:-1:1
    best_alpha = -inf;    
    for transition_index = 1:size(transitions,1)
        if transitions(transition_index, 2) == state && alphas(transitions(transition_index,1),bit_index) > best_alpha
            best_alpha = alphas(transitions(transition_index,1), bit_index);
            b_hat(bit_index) = transitions(transition_index,3);
            best_state = transitions(transition_index,1);
        end       
    end
    state = best_state;
end

end

function s= add (a,b)
 s = a + b;
end

function m = maxstar(a,  b)
if (a==-inf) && (b==-inf)
    m = -inf;
else
    
    m =  max(a,b) + log(1+exp(-abs(a-b)));
end
end