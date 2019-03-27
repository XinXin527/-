% Rayleigh scattering
% Use SI units
function Qs=Rayleigh(x)
lambda=400*1e-9;
n_sphere=1.40
n_background=1.33;
w_sphere=1.05*1e-3/(1e-2)^3;
w_background=1.0*1e-3/(1e-2)^3;
concentration=1e-5;
n_rel=n_sphere/n_background;

Qs = 8*x^4/3*abs((n_rel^2 - 1)/(n_rel^2 + 2))^2;