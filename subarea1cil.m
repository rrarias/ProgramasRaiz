function [amasn,amasnnx,amasnny,rpfy,rpfx]=subarea1cil(x,y,vert,celdas,tria,ind,tama,neib)

% 	neib=[];
% 	for q=1:length(tria)
% 		for jj=1:3
% 			if x(ind)==x(tria(q,jj)) 
% 			neib=[neib;tria(q,jj+1)];
% 			end
% 		end
% 	end
% 	neib=unique(neib);

    neibr=neib{ind}';
    neibr=neibr(1:find(neibr<=tama, 1, 'last' ));
	an=[];
    bn=[];
	rx=[];
	ry=[];
	amasn=0;
	amasnnx=0;
	amasnny=0;
    rpfx=0;
    rpfy=0;
   
    
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
		amasn=amasn+0.25*an(nn)*bn(nn);
		amasnnx=amasnnx+0.25*an(nn)*((x(ind)-x(neibr(nn)))/bn(nn));		
		amasnny=amasnny+0.25*an(nn)*((y(ind)-y(neibr(nn)))/bn(nn));
		rpfy=((y(ind)-y(neibr(nn)))-0.7);
        
    end
	



