%function CCin=concet(x,y,vert,celdas,CCp,tria,ind,indm,ncc)
xb=[];
yb=[];
JJp=[];
xf=x;%(xori+length(8xv2xv2)+1:end);
yf=y;%(xori+length(xv2)+1:end);
for j=1:xori+length(xv2)

    neibr=neib{j}';
    %neibr=neibr(1:find(neibr<=xori+length(xv2), 1, 'last' ));

  
    for nn=1:length(neibr)  
        if isempty(find((xf(neibr(nn))-xv3f1==0) & (yf(neibr(nn))-yv3f1==0), 1))==0
            
            xb=[xb xf(j)];
            yb=[yb yf(j)];
            %if isempty(j~=JJ(1:length(JJ)))
            JJp=[JJp j];
            %end 
        end  
    end  
end
JJp=unique(JJp);
%JJp(find(JJp>xori))=[];
