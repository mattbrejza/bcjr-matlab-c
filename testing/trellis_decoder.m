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
function [ztildee, ytildep] = trellis_decoder(ztildea, C, Px_vector, l)

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

% C = settings.uec_c;
% Px_vector = settings.uec_Px_vector(size(C,1),:);
% Px_vector = Px_vector(1:size(C,1)-1);
% l = settings.uec_l;

if length(Px_vector) ~= size(C,1)-1
    error('length(Px_vector) should equal size(C,1)-1');
end

if size(ztildea,1) ~= size(C,2)
    error('size(ztildea,1) should equal size(C,2)');
end

% Determine number of states
r = 2*size(C,1);

% Determine number of codeword bits
n = size(C,2);

% Determine length of unary-encoded bit sequence
b = size(ztildea,2);

% Build the trellis and calculate conditional transition probabilities
transitions = zeros(2*r, n+3);
gammas = zeros(2*r, b);
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
        
        % See Eq (6) in Derivations for (Maunder et al., 2013)
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
        gammas(y*r+mprime,:) = log(Pmmprime);
        
    end
end


% Calculate a priori transition log-probabilities
for codebit_index = 1:n
    gammas(transitions(:,3+codebit_index)==1,:) = gammas(transitions(:,3+codebit_index)==1,:) + repmat(ztildea(codebit_index, :), sum(transitions(:,3+codebit_index)==1),1);
end

% Forward recursion to calculate state log-probabilities
alphas=-inf(r,b);
%alphas(1,1)=0;
%alphas(:,1)=0; % We know that this is the first state

%if settings.generate_unary_mode == 4
%	alphas(settings.start_trellis_state,1) = 0;
%elseif settings.trellis_knows_first>0
    alphas(1,1)=0; % We know that this is the first state
%else
%    alphas(:,1)=0; % We dont know that this is the first state
%end

for j = 2:b
    temp = alphas(transitions(:,1),j-1)+gammas(:,j-1);
    for m = 1:r
        alphas(m,j) = maxstar(temp(transitions(:,2) == m));
    end
end


  

% % Backward recursion to calculate state log-probabilities
betas=ones(r,b)*-4000;
% %betas(:,end)=0; % We know that this is the final state
% if settings.trellis_knows_last>0
%     betas(settings.end_trellis_state,end) = 0;
% else
    betas(:,end)=0; % We dont know that this is the first state
% end
for j = b-1:-1:1
    temp = betas(transitions(:,2),j+1)+gammas(:,j+1);
    for mprime = 1:r
        betas(mprime,j) = maxstar(temp(transitions(:,1) == mprime));
    end
end

% Calculate a posteriori transition log-probabilities
deltas = alphas(transitions(:,1),:) + betas(transitions(:,2),:) + gammas;

% % Calculate a posteriori LLRs pertaining to y and z
% yztildep = zeros(n+1,b);
% for codebit_index = 1:n+1
%     for bit = 1:b
%         log_p1=max(deltas(transitions(:,2+codebit_index) == 1,b));
%         log_p0=max(deltas(transitions(:,2+codebit_index) == 0,b));
%         yztildep(codebit_index,b) = log_p1-log_p0;
%     end    
%     
%     
% end
yztildep = zeros(n+1,b);
for codebit_index = 1:n+1
    log_p1=maxstar(deltas(transitions(:,2+codebit_index) == 1,:));
    log_p0=maxstar(deltas(transitions(:,2+codebit_index) == 0,:));
    yztildep(codebit_index,:) = log_p1-log_p0;
end

% Extract a posteriori LLRs pertaining to y
ytildep = yztildep(1,:);

% Calculate extrinsic LLRs pertaining to z
ztildee = yztildep(2:end,:) -ztildea;

end