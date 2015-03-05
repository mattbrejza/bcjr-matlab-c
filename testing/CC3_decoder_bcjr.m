% BCJR algorithm for an accumulator.
% Copyright (C) 2013  Robert G. Maunder

% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

% The GNU General Public License can be seen at http://www.gnu.org/licenses/.

% apriori_uncoded_llrs is a 1xN vector of a priori uncoded LLRs
% apriori_encoded_llrs is a KxN vector of a priori encoded LLRs
% extrinsic_uncoded_llrs is a 1xN vector of extrinsic encoded LLRs
% extrinsic_encoded_llrs is a KxN vector of extrinsic encoded LLRs
function [extrinsic_uncoded, extrinsic_encoded] = CC3_decoder_bcjr(apriori_uncoded, apriori_encoded)

% Matrix to describe the trellis
% Each row describes one transition in the trellis
% Each state is allocated an index 1,2,3,... Note that this list starts
% from 1 rather than 0.
%               FromState,  ToState,    UncodedBit, EncodedBit1,EncodedBit2,EncodedBit3
transitions =  [1,          1,          0,          0,          0,           0;
                1,          3,          1,          1,          1,           1;
                2,          1,          0,          0,          1,           1;
                2,          3,          1,          1,          0,           0;
                3,          4,          0,          1,          0,           0;
                3,          2,          1,          0,          1,           1;
                4,          4,          0,          1,          1,           1;
                4,          2,          1,          0,          0,           0];
            
    if(size(apriori_uncoded,2) ~= size(apriori_encoded,2))
        error('LLR sequences must have the same length');
    end

    if(size(apriori_encoded,1) ~= size(transitions,2)-3)
        error('LLR sequences must have the correct dimensions');
    end
            
            
    % Find the largest state index in the transitions matrix           
    % In this example, we have eight states since the code has three memory elements
    state_count = max(max(transitions(:,1)),max(transitions(:,2)));

    
    
    % Calculate the a priori transition log-probabilities
    gammas = zeros(size(transitions,1),length(apriori_uncoded));
    gammas(transitions(:,3)==1,:) = repmat(apriori_uncoded, sum(transitions(:,3)==1),1);    
    for codebit_index = 1:size(apriori_encoded,1)
        gammas(transitions(:,3+codebit_index)==1,:) = gammas(transitions(:,3+codebit_index)==1,:) + repmat(apriori_encoded(codebit_index,:), sum(transitions(:,3+codebit_index)==1),1);
    end
    
    % Recursion to calculate forward state log-probabilities
    alphas=zeros(state_count,length(apriori_uncoded));
    alphas(2:end,1)=-inf; % We know that these are not the first state
    for bit_index = 2:length(apriori_uncoded)        
        temp = alphas(transitions(:,1),bit_index-1)+gammas(:,bit_index-1);
        for state_index = 1:state_count
            alphas(state_index,bit_index) = maxstar(temp(transitions(:,2) == state_index));
        end
    end
    
    % Recursion to calculate backward state log-probabilities
    betas=zeros(state_count,length(apriori_uncoded));
    for bit_index = length(apriori_uncoded)-1:-1:1
        temp = betas(transitions(:,2),bit_index+1)+gammas(:,bit_index+1);
        for state_index = 1:state_count
            betas(state_index,bit_index) = maxstar(temp(transitions(:,1) == state_index));
        end
    end

    % Calculate a posteriori transition log-probabilities
    deltas = alphas(transitions(:,1),:) + betas(transitions(:,2),:) + gammas;
    
    % Calculate the uncoded extrinsic LLRs
    log_p0=maxstar(deltas(transitions(:,3) == 0,:));
    log_p1=maxstar(deltas(transitions(:,3) == 1,:));
    extrinsic_uncoded = log_p1-log_p0-apriori_uncoded;
    
    extrinsic_encoded = zeros(size(apriori_encoded));
    for codebit_index = 1:size(apriori_encoded,1)
        log_p0=maxstar(deltas(transitions(:,3+codebit_index) == 0,:));
        log_p1=maxstar(deltas(transitions(:,3+codebit_index) == 1,:));
        extrinsic_encoded(codebit_index,:) = log_p1-log_p0-apriori_encoded(codebit_index,:);
    end
   
