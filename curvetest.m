%%
% %******** funcion para calcular la nueva frontera, [xv yv]
% nAo=coa;
% 
% addcoa1=addcoa1+nAo;
% 
% nsourup=(pi/(2*beta*dt*20))^2;
% addcoa=(-1+addcoa1)*nsourup;
% dfup=dfup+3*dsourup*dt;
% dsourup=dsourup+3*(-nsourup*dfup+addcoa)*dt;
% 
% addcoa1=addcoa1-dfup;
% 
% if addcoa1<0
%     addcoa1=0;
%     dsourup=0;
%     dfup=0;
% end
%  if dsourup<0
%  dsourup=0;
%  dfup=0;
%  end
% 
% 
% 
% updw=updw-(coa/(2*parax3)/deltaup);
coa=0;






    b43=b33(end)-ib2:-ib3:updw+.3;
    a43=a33(end)*ones(1,length(b43))+0.001*rand(1,length(b43));
    
    b53=b43(end):ib3:b33(end);
    a53=(3*parax3+2*raddw3)*ones(1,length(b53))+0.001*rand(1,length(b53));
    

 
    a83=a33(end):ia3:(3*parax3+2*raddw3);
    b83=(updw-.4)*ones(1,length(a83));


xv3=[a13';a23';a33';a43';a83';a53';a63';a73'];
yv3=[b13';b23';b33';b43';b83';b53';b63';b73'];


xv3f1=[a13';a23';a33';a43';a53';a63';a73'];
yv3f1=[b13';b23';b33';b43';b53';b63';b73'];

xv3f=[a83'];
yv3f=[b83'];

%%
    
    
    
    
    
    
    
%%%%%%%%%%%%%%%%%%% Segunda fila 
%relajacion para que todos los puntos de la frontera tengan la misma
%distancia
ddt=dt;  %tiempo para relajacion de frontera
kbs=120;      %constante de armortiguacion frontera
kbsr=4;     %constante de disipacion frontera

%b42=sort(b42);





    
%xva=[xv2;xv3];
%yva=[yv2;yv3];

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






       


