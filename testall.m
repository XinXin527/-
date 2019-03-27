a=[-1:0.01:3];
for i=1:length(a)
    x(i)=10^a(i);
    Qs(i)=Rayleigh(x(i));
end
loglog(x,Qs);
hold on
for i=1:length(a)
    x(i)=10^a(i);
    Qs(i)=Mie(x(i));
end
loglog(x,Qs);
xlabel('x=ka');
ylabel('Scattering efficiency Qs');
legend('Rayleigh approximation');