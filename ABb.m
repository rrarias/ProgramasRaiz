%function CCin=concet(x,y,vert,celdas,CCp,tria,ind,indm,ncc)
xb=[];
yb=[];
JJ=[];
for j=1:nc
    
    neib=[];
    
    for q=1:length(tria)
        for jj=1:3
            if x(j)==x(tria(q,jj))
                neib=[neib;tria(q,jj+1)];
            end
        end
    end
    neib=unique(neib);
    

    
    for nn=1:length(neib)
        
        sc=0;
         if isempty(find((x(neib(nn))-x(nc+1:length(x))==0) & (y(neib(nn))-y(nc+1:length(x))==0), 1))==0
           sc=sc+1;  
         end
       
         end
         
         if sc > 0
             xb=[xb x(j)];
             yb=[yb y(j)];
             JJ=[JJ j]; 
         end
             
        
        
        
    
    
end

