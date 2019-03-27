a=[-1:0.01:3];
for i=1:length(a)
    x(i)=10^a(i);
    g(i)=Mie(x(i));
end
loglog(x,g);
xlabel('x=ka');
ylabel('Anisotropy g');