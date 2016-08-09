function [amasn,amasnnx,amasnny,amasnnz]=subarea1cil(x,y,vert,celdas,tria,ind,tama,neib,JJ1)


    neibr=neib{ind}';
    neibr=neibr(1:find(neibr>tama, 1, 'last' ));
    %neibr=intersect(neibr,JJ1);
	an=[];
    bn=[];
	cn=[];
	rx=[];
	ry=[];
	amasn=0;
	amasnnx=0;
	amasnny=0;
	amasnnz=0;
   
    
	for nn=1:length(neibr)
				
		rx=intersect(vert(celdas{ind},1),vert(celdas{neibr(nn)},1));
		ry=intersect(vert(celdas{ind},2),vert(celdas{neibr(nn)},2));
        numel(rx);
        if numel(rx)==0
            rx(1)=0;
            rx(2)=0;
            ry(1)=0;
            ry(2)=0;
        end
        if numel(ry)==0
            rx(1)=0;
            rx(2)=0;
            ry(1)=0;
            ry(2)=0;
        end
        if length(ry)<2 || length(rx)<2
            rx(2)=rx(1);
            ry(2)=ry(1);
        end
		bn=[bn;distance([x(ind),y(ind)],[x(neibr(nn)),y(neibr(nn))])];
        if bn(nn)==0
            bn(nn)=1;
        end
		an=[an;distance([rx(1),ry(1)],[rx(2),ry(2)])];
		amasn=amasn+0.25*an(nn)*(abs(y(ind)-y(neibr(nn)))-bn(nn));
		amasnnx=amasnnx-0.25*an(nn)*((x(ind)-x(neibr(nn)))/bn(nn));		
		amasnny=amasnny-0.25*an(nn)*(sign(y(ind)-y(neibr(nn)))+((y(ind)-y(neibr(nn)))/bn(nn)));
		amasnnz=amasnnz+(y(neibr(nn))-y(ind)-.9);
        
    end
	



