function [ transitions ] = gen_uec_transitions( C, Px_vector, len )
%GEN_UEC_TRANSITIONS Generates a Transitions X 5 matrix with the UEC
%transitions and the conditional transition probablity
%   Detailed explanation goes here

r = 2*size(C,1);

% Determine number of codeword bits
n = size(C,2);

if (n > 2)
    error('only works for 2 codeword bit atm');
end

l=len;

transitions = zeros(2*r, n+3);
gammas = zeros(2*r, 1);
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

transitions = [transitions , gammas];

end

