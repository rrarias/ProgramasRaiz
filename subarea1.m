function [amasn,amasnnx,amasnny]=subarea1(x,y,vert,celdas,tria,ind)


%[vert,celdas]=voronoin([x,y]);
%tri=delaunay(x,y);
%tria=[tri tri(:,1)];
	%ind=5;
	neib=[];
	for q=1:length(tria)
		for jj=1:3
			if x(ind)==x(tria(q,jj)) 
			neib=[neib;tria(q,jj+1)];
			end
		end
    end
   
	neib=unique(neib);
	an=[];
    bn=[];
	rx=[];
	ry=[];
	amasn=0;
	amasnnx=0;
	amasnny=0;
 
	for nn=1:length(neib)
		%if x(neib(nn))<all(x(K))		
		rx=intersect(vert(celdas{ind},1),vert(celdas{neib(nn)},1));
		ry=intersect(vert(celdas{ind},2),vert(celdas{neib(nn)},2));		
		bn=[bn;distance([x(ind),y(ind)],[x(neib(nn)),y(neib(nn))])];
		an=[an;distance([rx(1),ry(1)],[rx(2),ry(2)])];
        
        %ele=
		amasn=amasn+0.25*an(nn)*bn(nn);
		amasnnx=amasnnx+0.25*an(nn)*((x(ind)-x(neib(nn)))/bn(nn));		
		amasnny=amasnny+0.25*an(nn)*((y(ind)-y(neib(nn)))/bn(nn));
        
        


		%end		
	end      
	


