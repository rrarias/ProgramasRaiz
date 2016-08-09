

fluxin=0;

neibr=neib{j}';


rx=[];
ry=[];
rxaux=[];
ryaux=[];
bn=[];
an=[];

am=0;
JJpfa=JJpf(1:length(JJpf)-2);

for nn=1:length(neibr)
 
        rx=intersect(vert(celdas{j},1),vert(celdas{neibr(nn)},1));
        ry=intersect(vert(celdas{j},2),vert(celdas{neibr(nn)},2));
       
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
       
        an=[an;distance([rx(1),ry(1)],[rx(2),ry(2)])];
         
door=an(nn)*heaviside((Ep(neibr(nn))-Ep(j))*(c(neibr(nn))-c(j)))*sign(Ep(neibr(nn))-sum(Ep(1:xori))/(xori));
        
am=am+door*(gradcx(j)*gradvx(j)+gradcy(j)*gradvy(j))+Do*lapc(j)/alfa;

fluxin=fluxin+alfa*door;

end

fluxt=[fluxt; x(j) y(j) fluxin*gradcx(j)*Ep(j)+Do*gradcx(j) fluxin*gradcy(j)*Ep(j)+Do*gradcy(j)]; 


