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

% The unary encoder of Section III-A in (Maunder et al., 2013)
% x is the symbol vector comprising a symbols
% y is the bit vector comprising b bits, formed of concatenated unary codewords
function y = unary_encoder(x)
    cumulative = cumsum(x);
    y = ones(1,cumulative(end));
    y(cumulative) = 0;
end