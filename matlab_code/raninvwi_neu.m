function [Q, Qinv, detQ] = raninvwi(nu,S)

%  sampling from an Inverted Wishart distribution

% input: nu ... shape parameter
%        S .... scale parameter

n=2*nu+1;
r=size(S,1);
%[C ierror]=chol(S)';

[us t] = schur(S);

Tdiag=max(diag(real(t)),0);
detS = sum(log(Tdiag));

th=diag(Tdiag.^.5);
u=real(us)*th;

thinv=diag(Tdiag.^(-.5));
uinv=real(us)*thinv;


z=randn(n,r);
E=nu*cov(z,1);

Q = u*inv(E)*u';
Qinv = uinv*E*uinv';
detQ=log(det(E))-detS;
