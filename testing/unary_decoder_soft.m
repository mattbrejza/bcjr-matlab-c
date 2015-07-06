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

% The soft-input hard-output unary decoder of Section IV-C in (Maunder et al., 2013)
% ytildep is the LLR vector comprising b a posteriori LLRs, where LLR=ln[Pr(bit=1)/Pr(bit=0)]
% a is the number of symbols to decode
% xhat is the symbol vector comprising a symbols
function xhat = unary_decoder_soft(ytildep, a)

yhat = ones(size(ytildep));
[~,indices] = sort(ytildep);
yhat(indices(1:a)) = 0;
xhat = unary_decoder_hard(yhat);

end

function xhat = unary_decoder_hard(yhat)
indices = find(yhat==0);
xhat = [indices(1), diff(indices)];
end