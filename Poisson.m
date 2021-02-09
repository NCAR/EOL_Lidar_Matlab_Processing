clear all
alpha = 50;
min = 0;
max = 100;
s=min:1:max;
P = (alpha.^s*exp(-alpha))./factorial(s);
figure(1)
bar(s,P)
area_P = trapz(P)

%Y = poisspdf(s,alpha);
%figure(2)
%bar(s,Y)
%area_Y = trapz(Y)

%v = var(P,)
%pd = fitdist(P', 'Poisson');
%v = var(pd)
%P = @(s,alpha) (alpha.^s*exp(-alpha))./factorial(s);
%q = integral(@(s) P(s,alpha),min,max)




