function [ac, hm, sd, med, ki, ineff] = autocovneu(h,r)
%
% h ... anzahl der spalten = anzahl der variablen
% r ... maximaler lag
% output:
% ac(1:r,k) 
%   k .. anzahl der variablen
%   ac(s,k) autocorrelation at lag s
% hm(1:k)
%   hhat(j) ... MW der j.ten variable
% sd(1:k)  ... SD der j.ten Variablen

% varhhat(1:k,1:k)
% ineff(1:k) .. inefficiency compared to i.i.d mean
%   
h1=h;
h=h';
n = size(h,2);
k = size(h,1);
hm = mean(h,2);
med = median(h,2);
omega = zeros(k,k,r+1);
omegat = omega;
svec = zeros(k,k,r+1);
ac = zeros(r,k);
ineff = zeros(k,1);

for j=1:k
hhat(j,1:n) = hm(j);
end

for s = 1:r+1
omegas(1:k,1:k) = (h(:,s:n)-hhat(:,1:n-s+1))*(h(:,1:n-s+1)-hhat(:,1:n-s+1))'./n;
omega(:,:,s) = omegas(:,:);
omegat(:,:,s) = omegas(:,:)';
svec(1:k,1:k,s) = (1-(s-1)/(r+1))/n;
end
svec(:,:,1) = svec(:,:,1)*0.5;

sd=(diag(omega(:,:,1))).^.5;

% neue Berechnung für KI

[ssort] = sort(h1);
ssort = ssort';
F = floor(n*0.05);
kiu = zeros(k,F+1);
kio = zeros(k,F+1);
kid = zeros(k,F+1);
ki = zeros(k,2);
for i = 0:F,
	kiu(:,i + 1) = ssort(:,i + 1);
	kio(:,i + 1) = ssort(:,(n -(F - i)));
	kid(:,i + 1) = kio(:,i+1) - kiu(:,i+1);
end
[y I] = min(kid,[],2);
for i = 1:k,
   ki(i,:) = [kiu(i,I(i)) kio(i,I(i))];
end   


for j=1:k
ac(:,j) = squeeze(omega(j,j,2:end))./omega(j,j,1);
end

varhhat = squeeze(sum(svec.*(omega+omegat),3));
for j=1:k
ineff(j) = varhhat(j,j)/omega(j,j,1)*n;
end


