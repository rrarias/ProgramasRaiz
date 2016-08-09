%function CCin=concet(x,y,vert,celdas,CCp,tria,ind,indm,ncc)
xb=[];
yb=[];
JJpf=[];
xf=x;%(xori+length(xv2)+1:end);
yf=y;%(xori+length(xv2)+1:end);
for j=1:xori

    neibr=neib{j}';
    %neibr=neibr(1:find(neibr<=xori+length(xv2), 1, 'last' ));

  
    for nn=1:length(neibr)  
        if isempty(find((xf(neibr(nn))-xv3f==0) & (yf(neibr(nn))-yv3f==0)))==0
            
            xb=[xb xf(j)];
            yb=[yb yf(j)];
            %if isempty(j~=JJ(1:length(JJ)))
            JJpf=[JJpf j];
            %end 
        end  
    end  
end
JJpf=unique(JJpf);
%JJp(find(JJp>xori))=[];
