% checa velocidad y taza de reproduccion
ncr=0;
aux=unique(monimeris(:,2));
for i=2:cuentas-2
    ncr=ncr+xoric(i);
    [ab abb]=sort(YT(1:ncr,i));
    veloy(abb,i)=(YT(abb,i)-YT(abb,i-1))/((aux(i+1)-aux(i)*dt*beta));
    velox(abb,i)=(XT(abb,i)-XT(abb,i-1))/((aux(i+1)-aux(i)*dt*beta));
    Velo(abb,i)=sqrt(VXT(abb,i).^2+VYT(abb,i).^2);
    
    plot(Velo(abb,i)*100);
    %axis([0 nc 0 0.001]);
    F2r=getframe;
F2Nr(i)=F2r;
end
% 
% 
% for mir=1:nc
%     
%     if all(listmonitor==mir)
% 
%         
        
        
for ir=2:length(listmonitor)
    timerate=0;
    arr=find(moni(listmonitor(ir),:)~=0);
    if length(arr)<2
        %timerate=moni(ir,arr(1));
    else
        for irr=2:length(arr)
            timerate=(timerate+abs(moni(ir,arr(irr))-moni(ir,arr(irr-1))))/(irr-1);
        end
    end
    raterepro(ir)=timerate;
end