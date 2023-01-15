function [x,f] = dichte(parsp)

%  parsp ... Spaltenvektor

par = parsp';

%  par ... zeilenvektor!!!!

par = squeeze(par);
pars = sort(par);
nd = 100;
smooth = sqrt(10);
m = size(pars,2);
alpha = 0.;
n0 = fix(alpha*m)+1;
n1 = fix((1-alpha)*m);
x = 1:nd;
x = pars(n0) + x'*(pars(n1)-pars(n0))/nd;
sd = std(par)/smooth;
const = 1./((2*pi)^0.5*sd*m);
par = par(ones(100,1),:);
xp = x(:,ones(1,m));
std22 = -2*sd^2; 
fm = exp((xp-par).^2/std22);
f = const*sum(fm,2);