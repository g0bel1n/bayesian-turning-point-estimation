function [mr,sr,r]=draw_comps_mix(res);
% drawing the components of the mixture for the residuals of the regression on ln tau  

p=       [0.2924  0.2599  0.2480  0.1525  0.0472]';
mu=      [0.0982 -1.5320 -0.7433  0.8303 -3.1428]';
sig=sqrt([0.2401  1.1872  0.3782  0.1920  3.2375]');

r=ones(size(res)); % component indicator   
mr=mu(1)*ones(size(res)); % component indicator   
sr=sig(1)*ones(size(res)); % component indicator   

lpt=zeros([size(res) 5]); 

for s=1:5;            
    % log-likelihood for each component
    lpt(:,:,s) =-log(sig(s))-0.5*(((res-mu(s))./sig(s)).^2)+log(p(s));
end      
cp=cumsum(exp(lpt),3);
cp=cp./cp(:,:,5*ones(1,5));

u=rand(size(res)); 
for j=1:4; % 4= number of components-1 %
    is=(u>cp(:,:,j)); % is is equal to 1, if the component falls into this category
    r=r+is;
    mr= mr.*(1-is)+mu(j+1).*is;
    sr= sr.*(1-is)+sig(j+1).*is;
end;