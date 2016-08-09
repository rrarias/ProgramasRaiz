%%
% %******** funcion para calcular la nueva frontera, [xv yv]

%if coa>0 %& sobre el potencial  

%% crecimiento de la base
nAo=coa;

addcoa1=addcoa1+nAo;

nsourup=(pi/(2*beta*dt*20))^2;
addcoa=(-1+addcoa1)*nsourup;
dfup=dfup+3*dsourup*dt;
dsourup=dsourup+3*(-nsourup*dfup+addcoa)*dt;

addcoa1=addcoa1-dfup;

if addcoa1<0
    addcoa1=0;
    dsourup=0;
    dfup=0;
end
 if dsourup<0
 dsourup=0;
 dfup=0;
 end



up=up+(dfup/(2*parax2)/deltaup);
%up=up+(coa/(2*parax2)); % el cilindro exterior donde todas las celulas son cerradas

coar=coar+(coa/(2*parax2));
sidedivaux=sidediv;
lxv2=length(xv2);
lxv2paracil1=length(a12)+length(a22);
lxv2paracil2=length(a42);
ib2f=0.7*ib2;               
coa=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%% Tercer fila en forma parabolica,,aumento de
%%%%%%%%%%%%%%%%%%%%%%%%%% dominio

b23=(ppara3*parax3^2)-1.2:ib3:up;
a23=(parax3-.14)*ones(1,length(b23));
a33=(parax3-.21):-ia3:-(parax3-.14);
b33=(up+.15)*ones(1,length(a33));
%b43=up-.05:-ib3:(ppara3*parax3^2)-1.3;
b43=b23;
a43=-(parax3-.14)*ones(1,length(b43));

xv3=[a13';a23';a33';a43'];
yv3=[b13';b23';b33';b43'];
xv3f1=[a13';a23';a43'];
yv3f1=[b13';b23';b43'];
xv3f=[a33'];
yv3f=[b33'];

%%

    
%%  Constantes el??sticas en la frontera    
    
%%%%%%%%%%%%%%%%%%% Segunda fila 
%relajacion para que todos los puntos de la frontera tengan la misma
%distancia
%ddt=dt;  %tiempo para relajacion de frontera


xv=[xv2;xv3];
yv=[yv2;yv3];

x=x(1:xori);
y=y(1:xori);
vx=vx(1:xori+length(yv2));
vy=vy(1:xori+length(yv2));
r=r(1:xori+length(yv2));
Ep=Ep(1:xori+length(xv2));
Ec=Ec(1:xori+length(xv2));
c=c(1:xori+length(xv2));
x=[x;xv];
y=[y;yv];
vx=[vx;zeros((length(xv3)),1)];
vy=[vy;zeros((length(yv3)),1)];
Ep=[Ep;zeros((length(xv3)),1)];
Ec=[Ec;zeros((length(xv3)),1)];
r=[r;zeros((length(yv3)),1)];
c=[c;zeros((length(yv3)),1)];

%%
% if prr<4
%         x(xori+1)=centroxx(xori+1);
%         x(xori+8)=centroxx(xori+8);
%         y(xori+1)=centroyy(xori+1);
%         y(xori+8)=centroyy(xori+8);
% %         
% %         x(xori+15)=centroxx(xori+15);
% %         x(xori+18)=centroxx(xori+18);
% %         y(xori+15)=centroyy(xori+15);
% %         y(xori+18)=centroyy(xori+18);
% end


       

                        
 


 


