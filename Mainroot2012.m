%Programa para relajar y hacer crecimiento en funcion del area
%****************************************************************
%********* condiciones iniciales ********************************

clear all
% XXXXXXXXXXXXXX    DATOS   XXXXXXXXXXXXXXXXXXXXXXXXX
%%%%%%%%%% constantes de fuerza
kx=3;
ky=kx;
kvol=80;
kcen=80;

cmin=0.005;% .0001;   % el cutoff de la concentracion de auxinas.
dt=0.01;  %3 minutos
mo=20;    % numero de iteraciones para la relajacion
gamma=200;%200       %tiempo de reproduccion (gamma*dt=10 hrs)

alfa=-50;  %     el treshhold de la c para la rep.
%beta=gamma/(20*mo*dt);   %deltat>alfa*dt
beta=35; %50
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc=15;     %     numero inicial de celulas, da la densidad inicial
parax=1.2;%1.5;   %    ancho del dominio
ppara=1;%1.;     %  coefciente de la parabola
d=2;%5;
up=ppara*parax^2+d;    % altura inicial del cilindro
bo=(up-d)/2; %    posicion del nicho

up1=up;           %altura mas incremento
ua=0.1;                        % incremento de altura
boc=2;          %alturade reproduccion vertical
delta=1;%.65;%     division vertical
deltaup=delta;   % tolerancia del area
deltadown=delta;
ia=0.1;%0.1;  % intervalo puntos parabola
ib=0.02;%0.2;    % intervalo puntos cilindro

%indice=0;  % indice del nicho (cambia)
invc=0;   % contador del incho
timequies=10; %     tasa de reprod del nicho
xad=[];
cuentas=1;   %  contador del numero de reproducciones por ciclo
pp=0.;  % porcentaje de area de la parabola creciendo pp=0 solo el rectangulo crece

%*********************************************************
moni=zeros(600,1); % variable que me ayuda a monitoriar las celulas
monipy=moni;
monimeris=[];
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%% **************** condiciones iniciales************************************************

x=1.5*(parax-ib)*(rand(nc,1)-0.5);%0.5==1.5

y=ppara*x.^2+((up-5*ia)-(ppara*x.^2+5*ia)).*rand(nc,1)+10*ia;
%XXXXX   posicion de las 2 troncales XXXXXXXXXXXXX
x(1)=-0.051;
x(2)=0.05;
y(1)=0.7*ppara+.0001; %0.7=0.3
y(2)=0.7*ppara-.0001;

%load xycinicial
%x=xin;
%y=yin;
xo1=x(1);
xo2=x(2);
yo1=y(1);
yo2=y(2);
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%y=(0.5*up + (0.8)*up*(rand(nc,1)-0.5));
% load pinicial
% Vth=.1;
% x=xin;
% y=yin;
vx=zeros((length(x)),1);
vy=vx;
Ep=vx;
Ec=vx;
c=vx;
r=vx;

%********************Boundaries**************************************
%%
a1(1)=-(parax-0.01);
    ii=1;
    while  a1(ii)<=(parax-0.01)
        ii=ii+1;
        a1(ii)=a1(ii-1) + .01/sqrt(4*ppara*a1(ii-1)^2+.05);
    end
    b1=ppara*a1.^2;
    b2=(ppara*parax^2):ib:up-.1;
    a3=-parax*ones(1,length(b2));
    a2=-a3;
    b3=b2;
    a4=-parax:ia:parax;
    b4=(up)*ones(1,length(a4));



xv=[a1';a2';a3';a4'];
yv=[b1';b2';b3';b4'];
%*****************************************************************
xori=length(x);
xoric(cuentas)=xori;
x=[x;xv];
y=[y;yv];
vx=[vx;zeros((length(xv)),1)];
vy=[vy;zeros((length(yv)),1)];
Ep=[Ep;zeros((length(xv)),1)];
Ec=[Ec;zeros((length(xv)),1)];
c=[c;zeros((length(xv)),1)];
areacero=((up-ppara*parax^2)*2*parax+(parax^3)/3)/nc;
ain=areacero;
rho=1/areacero;


XT(:,1)=zeros(1500,1);
YT=XT;
VXT=XT;
VYT=XT;
EPT=XT;
CT=XT;

tf=100;%    numero final de iteraciones
counttotal=1; %contador de iteraciones totales


%load runsept23    % uncomment si se dispone de un calculo anterior

 %prr=1;
 tri= delaunay(x,y);
        %tria=[tri tri(:,1)];
%encuentra las celdas de la orilla
AB
c=0*x;
c(1:nc)=.05*rand(nc,1);
for i=1:length(JJ)
c(JJ(i))=.2*rand;
end

%c(1:nc)=cin(1:nc);

c(1)=max(c);
c(2)=c(1);
c=c/sum(c);

% fin=length(x)-length(xv);
% x=x(1:fin);
% y=y(1:fin);
% 
% x=[x;xv];
% y=[y;yv];

XT(1:length(x),1)=x(1:end);
YT(1:length(x),1)=y(1:end);
VXT(1:length(vx),1)=vx(1:end);
VYT(1:length(vy),1)=vy(1:end);
EPT(1:length(Ep),1)=Ep(1:end);
CT(1:length(c),1)=c(1:end);

% for i=1:size(x)
%     if x(i)==xoo && y(i)==yoo
%         nii=i;
%     end
% end

%%
for prr=1:240  %para ciclos de visualizacion
    

    %**********************************************************
     for m=1:mo %tiempos
           
        [num2str(prr) '.' num2str(m) '  invc = ' num2str(invc)]
        %x=[x;xv];
        %y=[y;yv];
        counttotal=counttotal+1;
        tri= delaunay(x,y);
        tria=[tri tri(:,1)];

        [vert celdas]=voronoin([x,y],{'Qbb','Qz'});
        K=convhull(x,y);
        %areacero=polyarea(x(K),y(K))/nc;    %Area de la frontera
        AB
        
        %xold=x;
        %yold=y;
        AREA=[];
        %BN=[];
        %AREANX=[];
        %AREANY=[];
        centroxx=zeros(length(x),1);
        centroyy=centroxx;
        for j=1:nc                 %for particule
       
%                 
                [centrox centroy]=centroid(vert(celdas{j},1),vert(celdas{j},2));				
                [amasnj,amasnnxj,amasnnyj]=subarea1(x,y,vert,celdas,tria,j);
                AREA=[AREA;amasnj];
                centroxx(j)=centrox;
                centroyy(j)=centroy;
                xj1=(x(j)-centrox);
                yj1=(y(j)-centroy);
                %rjo=distance([x(j),y(j)],[centrox,centroy]);
                %rj=sqrt((x(j))^2 + (y(j))^2)+ 0.0001;
                %[area,AREANX,AREANY]=subarea(x,y,vert,celdas,tria,j);
                
                x(j)=x(j)+ vx(j)*dt;
                y(j)=y(j)+ vy(j)*dt;
                vx(j)=vx(j)+(kvol*(amasnj-areacero)*amasnnxj-kcen*xj1-kx*vx(j))*dt;
                vy(j)=vy(j)+(kvol*(amasnj-areacero)*amasnnyj-kcen*yj1-ky*vy(j))*dt;
                Ep(j)=kvol*(amasnj-areacero)^2/2+kcen*(xj1^2+yj1^2)/2;
                Ec(j)=((vx(j).^2+vy(j).^2))/2;
            
            simpdir
            %alfa=beta/(ams*gamma);
            %alfa=-20;%*mo;
            c(j)=c(j)+dt*(am)*alfa;
             r(j)=r(j)+dt*beta/sqrt(c(j)-cmin);

        end
                x(1)=xo1;
                y(1)=yo1;
                vx(1)=0;
                vy(1)=0;   %las nicho no se mueven
                x(2)=xo2;
                y(2)=yo2;
                vx(2)=0;
                vy(2)=0;         
     
       areacero=(sum(AREA)/nc);
        
 %=======Fuentes=================               
        c(1)=max(c); 
        c(2)=c(1);
 %===============================
        c=c/sum(c);
        %r=r/sum(r);
   
    
%%    reproduction
        
        ar=find(r(1:nc).*sqrt(c(1:nc))/gamma>=1);
        %bocn=find((c(1:nc)/max(c(1:nc)))<=0.8); %sustituye a boc en funcion de las auxinas
        if length(ar)>0
        reptot
        r(ar)=0;
        end

        %c(1:nc)=sqrt(vx(1:nc).^2+vy(1:nc).^2)./(AREA(1:nc).*(Ep(1:nc)+Ec(1:nc)'));
        %sum(c)
        %c=c/sum(c);
        
%         fin=length(x)-length(xv);
%         x=x(1:fin);
%         y=y(1:fin);
% 
% %%
%         figure(1)
%         clf
%         hold on
%         %voronoi(x,y)
%         plot(x,y,'r.');
%         %plot(centroxx,centroyy,'g.');
%         plot(xb,yb,'*k');
%         %axis([45 55 255 270]);
%         axis equal
%         hold off
%         F1=getframe;
%         F1N(prr)=F1;
%         %xp=x(1:length(xold));
%         %yp=y(1:length(xold));
%         %x=[];
%         %y=[];
%         %x=xp;
%         %y=yp;
%         
        
    end  % tiempos
    
    
%**********************************************************



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %figure(5)
    %clf
    [vv,cc] = voronoin([x,y]);
    x5=[a3 a1 a2 a4];
    y5=[b3 b1 b2 b4];
    [IN ON] = inpolygon(x,y,x5,y5);
    for i = 1:length(cc)
        if ON(i)==0;
            nuevovoro{i}=cc{i};
            nuevovoro{i}=horzcat(nuevovoro{i},nuevovoro{i}(1));
        end
    end
%     hold on
%     for i=1:length(nuevovoro)
%         line(vv(nuevovoro{i},1),vv(nuevovoro{i},2),'Color','k')
%     end
%          plot(x(1:2),y(1:2),'*b')
%          hold off
%     axis equal
%     axis off
%     
%     M=getframe(gcf);
%     MN(prr)=M;

 %%    %**********************************************************
figure(3)
clf
[aa bb]=sort(y(1:nc));

%subplot(2,1,1)
axis([0 nc 0 1]);
em=max(Ep(bb));
hold on
%bar((Ep(bb))/(max(em)),.1)
plot((Ep(bb))/em);
plot(c(bb)/max(c(bb)),'r')
plot(real(r(bb).*sqrt(c(bb))/gamma),'k')
%pause(.01)
%plot(y(bb==1),0,'*r')
plot((find(bb==1)),0,'*r')
hold off
F2=getframe;
F2N(prr)=F2;

%%********************

XT(bb,prr+1)=x(bb);
YT(bb,prr+1)=y(bb);
VXT(bb,prr+1)=vx(bb);
VYT(bb,prr+1)=vy(bb);
EPT(bb,prr+1)=Ep(bb);
CT(bb,prr+1)=c(bb);




%%    
% figure(3)
% 
% trisurf(tri,x,y,Ep*2/max(Ep))
% shading interp
% axis equal
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     figure(4)
clf

    %subplot(2,1,2)

     fEc=10*c/(max(c));
     %fE(K)=10;
     
    Z=nuevovoro;
    for i=1:length(nuevovoro) 
        tt=Z{i};
        tt(:)=1;
        Z{i}=fEc(i)*tt;

    end

    hold on
    axis xy ,colormap 'jet'
    axis equal, axis off,axis([-2 2 0 max(y) 0 max(fEc)]),view(0,90);% view(26,84);
    for i=1:length(nuevovoro)
            line(vv(nuevovoro{i},1),vv(nuevovoro{i},2),'Color','k')
            patch(vv(nuevovoro{i},1),vv(nuevovoro{i},2),Z{i}',Z{i}')

                A=vv(nuevovoro{i},1);
                B=vv(nuevovoro{i},2);
                E=Z{i}';
                Xi=[A,B,E];
                Xf=[A,B,0*E];
                hold on
               for j=1:length(A)         
                line('xdata',[Xi(j,1),Xf(j,1)],'ydata',[Xi(j,2),Xf(j,2)],'zdata',[Xi(j,3),Xf(j,3)],'linestyle',':')
               end                   
    end
     plot3(x(1),y(1),fEc(1),'dw','markersize',4)
     plot3(x(2),y(2),fEc(2),'dw','markersize',4)
          hold off    
   

     Q=getframe(gca);
    QN(prr)=Q;
    
%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     figure(2)
clf

    %subplot(2,1,2)

     fE=4*Ep/(max(Ep));
     %fE(K)=10;
     
    Z=nuevovoro;
    for i=1:length(nuevovoro) 
        tt=Z{i};
        tt(:)=1;
        Z{i}=fE(i)*tt;

    end

    hold on
    axis xy ,colormap 'jet'
    axis equal, axis off,axis([-2 2 0 max(y) 0 max(fE)]), view(0,90);%view(26,84);
    for i=1:length(nuevovoro)
            line(vv(nuevovoro{i},1),vv(nuevovoro{i},2),'Color','k')
            patch(vv(nuevovoro{i},1),vv(nuevovoro{i},2),Z{i}',Z{i}')

                A=vv(nuevovoro{i},1);
                B=vv(nuevovoro{i},2);
                E=Z{i}';
                Xi=[A,B,E];
                Xf=[A,B,0*E];
                hold on
               for j=1:length(A)         
                line('xdata',[Xi(j,1),Xf(j,1)],'ydata',[Xi(j,2),Xf(j,2)],'zdata',[Xi(j,3),Xf(j,3)],'linestyle',':')
               end                   
    end
     plot3(x(1),y(1),fE(1),'dw','markersize',2)
     plot3(x(2),y(2),fE(2),'dw','markersize',2)
          hold off    
   

     Q1=getframe(gca);
    QN1(prr)=Q1;
    

   %%
%%
% figure(6)
% 
%      [xima yima]=voronoi(x,y,'tri');
%      x4=[a1';a2';-a4';a3'];
%     y4=[b1';b2';b4';b3'];
%                
%     clf
%     hold on
%     patch(x4,y4,'y','facealpha',.2,'edgecolor','none')
%     patch([-parax-.2,a1,parax+.2,parax+.2,-parax-.2],[ppara*(parax)^2+.2,b1-.5,ppara*(parax)^2+.2,-.5,-.5],'k','facecolor','k','edgecolor','none')
%                 patch(vv(nuevovoro{1},1),vv(nuevovoro{1},2),'g','facealpha',.5)
%                     patch(vv(nuevovoro{2},1),vv(nuevovoro{2},2),'c','facealpha',.5) 
%     h = plot(xima,yima,'k-',x(1:nc),y(1:nc),'r.',xv,yv,'.g','markersize',3);
%     set(h(1:end-1),'xliminclude','off','yliminclude','off')
%     axis equal,axis([-parax-.2 parax+.2 -.5 up+.2]), box 'on'
%     %axis equal, axis([-parax-.2 parax+.2 -.5 8]), box 'on'
%     hold off
%         S=getframe(gcf);
%         SN(prr)=S;
%     
    %******************************************************

    
end   %*** ciclos


