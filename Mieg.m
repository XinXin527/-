% Mie theory
% Use SI units
function g=Mie(x)
lambda =400*1e-9;
n_s = 1.40;
n_b = 1.33;
w_s =1.05*1e3;
w_b = 1.0*1e3;
concentration =0.02;
n_rel = n_s/n_b
y = n_rel*x

% Calculate the summations
err = 1e-8;
Qs = 0;
gQs = 0;
for n = 1:1000
    Snx = sqrt(pi*x/2)*besselj(n+0.5,x);
    Sny = sqrt(pi*y/2)*besselj(n+0.5,y);
    Cnx = -sqrt(pi*x/2)*bessely(n+0.5,x);
    Zetax = Snx+i*Cnx;

    % Calculate the first-order derivatives
    Snx_prime = - (n/x)*Snx+sqrt(pi*x/2)*besselj(n-0.5,x);
    Sny_prime = - (n/y)*Sny+sqrt(pi*y/2)*besselj(n-0.5,y);
    Cnx_prime = - (n/x)*Cnx-sqrt(pi*x/2)*bessely(n-0.5,x);
    Zetax_prime = Snx_prime + i*Cnx_prime;

    an_num = Sny_prime*Snx-n_rel*Sny*Snx_prime;
    an_den = Sny_prime*Zetax-n_rel*Sny*Zetax_prime;
    an = an_num/an_den;

    bn_num = n_rel*Sny_prime*Snx-Sny*Snx_prime;
    bn_den = n_rel*Sny_prime*Zetax-Sny*Zetax_prime;
    bn = bn_num/bn_den;

    Qs1 = (2*n+1)*(abs(an)^2+abs(bn)^2);
    Qs = Qs+Qs1;

    if n > 1
        gQs1 = (n-1)*(n+1)/n*real(an_1*conj(an)+bn_1*conj(bn))...
            +(2*n-1)/((n-1)*n)*real(an_1*conj(bn_1));
        gQs = gQs+gQs1;
    end

    an_1 = an;
    bn_1 = bn;

    if abs(Qs1)<(err*Qs) & abs(gQs1)<(err*gQs)
        break;
    end
end

Qs = (2/x^2)*Qs;
gQs = (4/x^2)*gQs;
g = gQs/Qs;
end