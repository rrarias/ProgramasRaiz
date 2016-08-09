%function CCin=concet(x,y,vert,celdas,CCp,tria,ind,indm,ncc)


neib=[];

for q=1:length(tria)
    for jj=1:3
        if x(j)==x(tria(q,jj))
            neib=[neib;tria(q,jj+1)];
        end
    end
end
neib=unique(neib);

rx=[];
ry=[];
%bn=0;
an=0;

am=0;
ams=0;


for nn=1:length(neib)
   
    rx=intersect(vert(celdas{j},1),vert(celdas{neib(nn)},1));
    ry=intersect(vert(celdas{j},2),vert(celdas{neib(nn)},2));
    %bn(nn)=distance([x(j),y(j)],[x(neib(nn)),y(neib(nn))]);
    an(nn)=distance([rx(1),ry(1)],[rx(2),ry(2)]);
    
    %flujo de pines hacia arriba
    if isempty(find(j==JJ, 1))==0 & isempty(find(neib(nn)==JJ, 1))==0 & y(j)>y(neib(nn))    
        am=am+an(nn)*(c(j)-c(neib(nn)))*(abs(Ep(j)-Ep(neib(nn))));
        %ams=ams+an(nn)*(abs(Ep(j)-Ep(neib(nn))));
    end
    %flujo de pines hacia abajo
    if isempty(find(j~=JJ, 1))==0 & isempty(find(neib(nn)~=JJ, 1))==0 & y(j)<y(neib(nn))    
        am=am+an(nn)*(c(j)-c(neib(nn)))*(abs(Ep(j)-Ep(neib(nn))));
        %ams=ams+an(nn)*(abs(Ep(j)-Ep(neib(nn))));
    end
    %flujo de pines hacia los lados-centro
    if isempty(find(j~=JJ, 1))==0 & isempty(find(neib(nn)==JJ, 1))==0 & y(j)>y(neib(nn))    
        am=am+an(nn)*(c(j)-c(neib(nn)))*(abs(Ep(j)-Ep(neib(nn))));
        %ams=ams+an(nn)*(abs(Ep(j)-Ep(neib(nn))));
    end
    %am=am+an(nn)*(c(j)-c(neib(nn)))*(abs(Ep(j)-Ep(neib(nn))));
        %ams=ams+an(nn)*((Ep(j)-Ep(neib(nn))));
        
        
        %          if (x(j)~=xb(1:length(xb))) & y(j)~=yb(1:length(yb)) & x(neib(nn))==xb(1:length(xb)) & y(neib(nn))==yb(1:length(yb))
        %              s=0;
        %          end
        
%                 if isempty(find(x(neib(nn))~=x(nc+1:length(x)), 1))==0 & isempty(find(y(neib(nn))~=y(nc+1:length(y)), 1))==0
%                 
%                 end

    
end




%Diff(x,y,vert,celdas,CC,tria,j,m,nc)