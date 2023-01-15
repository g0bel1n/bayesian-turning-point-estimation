%K = size(etamc,2);
M = size(alphamc,1);

%Mb=size(bmc,1);
%N=size(bmc,2);
%dd=size(Q,1);

figure(double(gcf)+1);

if M< 500
    indexmc=[1:M];
else
    ispace=fix(M/500);indexmc=[1:ispace:M];
end   

for i1=1:dd
    subplot(nv,nh,i1);
    for k=1:K-(K>2)
        [xa,fa]=dichte(alphamc(indexmc,(k-1)*dd+i1));
        plot(xa,fa,str(k));hold on;
        xlabel(['\beta_{' num2str(i1) '}']);
        if i1==dd
            axis([-2 2 min(fa) max(fa)])
        end
    end
    if i1==1 title('group specific parameters (I_t=1)');
        legend(leg(1:k,:)) ;end
end

