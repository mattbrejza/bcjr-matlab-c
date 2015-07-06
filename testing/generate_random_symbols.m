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




% Generate a vector of random symbols
% a is the length of the symbol vector to generate
% distribution is a string specifying the name of a distribution function, e.g. 'distribution_zeta'
% parameter is the parameter of that distribution
% x is the vector of random symbols
function x = generate_random_symbols(a, distribution, parameter)

    x = zeros(1,a);
    random_number = rand(size(x));
    sum_prob = 0;
    symbol_value = 0;

    while sum(random_number > sum_prob) > 0
        symbol_value = symbol_value + 1;
        x(random_number > sum_prob) = symbol_value;
        sum_prob = sum_prob + feval(distribution, symbol_value, parameter);
    end
end