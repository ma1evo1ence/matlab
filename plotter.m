eV = 1.602176634e-19;
E = [-3.99*eV : 0.005*eV : 34*eV];
f = zeros(size(E));
for ii = 1 : length(E)
    f(ii) = F(E(ii)/eV);
end

figure; hold on; grid on; plot(E/eV, f)


function y = F(curE)
curE = curE*1.602176634e-19;
m = 9.109383701528 * 10^(-31);
U_0 = -4*1.602176634*10^(-19);
a = 0.5*10^(-9);
b = 2*10^(-9);
h = 1.054571817 * 10^(-34);
mu = sqrt(2*m*curE/(h*h));
lambda = sqrt(2*m*(curE - U_0)/(h*h));
y = cos(mu*a)*cos(lambda*b)-(lambda^2+mu^2)/(2*mu*lambda)*sin(mu*a)*sin(lambda*b);
end