function [IMS,IMS_enc] = sim_ecompIMS(etaMS_enc,lik);

%function to simulate two state indicator, whereby state ind. 2 is leading state ind. 1
ngrp=size(lik,2);
% if size(etaMS,3)==1;
%     etaMS=etaMS(:,:,ones(1,ngrp));
% end

% xi=[xi^2_11   xi^2_12          0            0
%       0       xia(1,1)         0          xia(1,2)
%     xia(2,1)    0          xia(2,2)          0
%       0         0          xi^2_21        xi^2_22 ];

% etaMS_enc=[etaMS(1,:,1) 0 0; ...
%            0  xia(1,1) 0 xia(1,2); ...
%            xia(2,1) 0 xia(2,2) 0; ...
%            0 0  etaMS(2,:,1)];

%state 1 = I_1t=1,I_2t=1 ; state 2 = I_1t=1,I_2t=2 ; state 3 = I_1t=2,I_2t=1  ; state 4 = I_1t=2,I_2t=2 
lik_enc = [lik(1,1,:)+lik(1,2,:);lik(1,1,:)+lik(2,2,:);lik(2,1,:)+lik(1,2,:);lik(2,1,:)+lik(2,2,:)];

IMS_enc = simstate_ms_enc(etaMS_enc,lik_enc)';

IMS(1:2,find(IMS_enc==1))=kron([1;1],ones(1,sum(IMS_enc==1)))-1;
IMS(1:2,find(IMS_enc==2))=kron([1;2],ones(1,sum(IMS_enc==2)))-1;
IMS(1:2,find(IMS_enc==3))=kron([2;1],ones(1,sum(IMS_enc==3)))-1;
IMS(1:2,find(IMS_enc==4))=kron([2;2],ones(1,sum(IMS_enc==4)))-1;
