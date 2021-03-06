function E = KronigPenney(k) %actually, k *(a+b) is inserted, from -pi to pi
% constants
eV = 1.602176634e-19;

global h k_n m a b U_0;
h = 1.054571817 * 10^(-34);
m = 9.109383701528 * 10^(-31);
a = 0.5e-9;
b = 2e-9;
U_0 = -4 * eV;

k = k / (a+b);
E = zeros(size(k));

%ranges of E for which F(E) is monotonous, abs(F(E)) <= 1
energy_range = [-3.9222, -3.9216; -3.69, -3.687; -3.305, -3.295; -2.77, -2.758;...
    -2.1, -2.06; -1.32, -1.24; -0.45, -0.29; 0.36, 0.78; 1.04, 1.784;...
    1.904, 2.9; 3, 4.1; 4.3, 5.54; 5.56, 7; 7.05, 8.6; 8.7, 10.35; ...
    10.45, 12.24; 12.28, 14.22; 14.24, 16.28; 16.38, 18.5; 18.6, 20.85;...
    20.92, 23.35; 23.366, 25.92];
 
figure; grid on; hold on; xlabel('k*(a+b)'); ylabel('Energy, eV')
for n = 1 : length(k) %cycle 
    k_n = k(n);
    column = []; %column for solutions of F(E) = cos(k(a+b)) for the given k
    %cycle that finds all solutions 
    for ii = 1 : length(energy_range)
        solution = fzero(@difference, energy_range(ii, :));
        if isempty(column(abs(column - solution) < 1e-10)) %checks if the solution is unique
            column = [column; 0];
            column(length(column)) = solution; %adds solution to column
        end
    end
    %plots the results upon the previous plots
    plot(k_n * (a+b), column, 'linestyle', 'none', 'marker', '.')
    E(1:length(column), n) = column; %inserts the column into E matrix
end

end

function x = F(E)
global h m a b U_0;
    mu = sqrt(2 * m * E) / h;
    lambda = sqrt(2 * m * (E - U_0)) / h;
    x = cos(mu*a) * cos(lambda*b) - (lambda^2 + mu^2)/(2*mu*lambda)...
        * sin(mu*a) * sin(lambda*b);
end

function x = difference(E)
global k_n a b;
E = E * 1.602176634e-19;
    x = real(F(E)) - cos(k_n * (a+b));
end
