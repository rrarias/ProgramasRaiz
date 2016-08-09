%Programa para relajar y hacer crecimiento en funcion del area
%****************************************************************
%********* condiciones iniciales ********************************

clear all
% XXXXXXXXXXXXXX    DATOS   XXXXXXXXXXXXXXXXXXXXXXXXX
%%%%%%%%%% constantes de fuerza
kx=3;      %constante de armortiguacion
ky=3;
kbsrx=3;     %constante de disipacion frontera
kbsry=3;     %constante de disipacion frontera
kvol=80;     %constante elastica celulas internas- area
kcen=80;     %constante elastica celulas internas- centroide
kbs=120;     %constante elastica celulas externas- area
kecen=120;   %constante elastica celulas externa- centroide

rep=0;       % reproduciones
edges=0:.2:30; %para generar histogramas de reproduccion

Nf=1500;      % pasos de tiempo

cmin=.0005;   % el cutoff de la concentracion de auxinas.
alfa=10;      % escala de tiempo para las auxinas
alfac=alfa;
Do=0.3;       % Coef difusion isotropico
beta=22;      % escala de tiempo de las ciclinas- celulas internas
betae=5;      % escala de tiempo de las ciclinas- celulas externas
gamma=200;    % ajusta el tiempo del ciclo de cyclinas
rcon=1;
rconMN=0;
sigma=0.0;   % se usa para ajustar la concentracion de auxina

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.01;   % escala de tiempo de la dinamica
nc=2;     %     numero inicial de celulas, da la densidad inicial
parax=2;   %    ancho del dominio
ppara=1;     %  coefciente de la parabola
d=2;        % altura del dominio
up=ppara*parax^2+d;    % altura inicial del cilindro
%up=4;
bo=(up-d)/2; %    posicion del nicho
%Ez=20;
up1=0;          %altura mas incremento
uax=d;          %incremento de altura
delta=2;        %division vertical
deltaup=delta;   % tolerancia del area
deltadown=delta;
ia=.2;  % intervalo puntos parabola
ib=.3;    % intervalo puntos cilindro

mo=5;    % numero de iteraciones por imagen

%indice=0;  % indice del nicho (cambia)
invc=0;   % contador del incho
timequies=1; %     tasa de reprod del nicho
xad=[];
cuentas=1;   %  contador del numero de reproducciones por ciclo
pp=0.;  % porcentaje de area de la parabola creciendo pp=0 solo el rectangulo crece
rxad=[];  %contar las coordenadas de las reproducciones
ryad=[];
sidediv=[];
coa=0;
coar=0;   %medidor de area
coaux=0;
dfup=0;   %condiciones iniciales para aumentar up, que es del dominio
dsourup=0;
addcoa1=0;
% adddfup=0;
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% **************** condiciones iniciales************************************************

x=.9*(parax-ib)*(rand(nc,1)-0.5);

y=ppara*x.^2+((up-5*ia)-(ppara*x.^2+5*ia)).*rand(nc,1)+5*ia;

%% ********************Boundaries**************************************
%condiciones de la frontera del dominio

%*********************Primera fila en forma parabolica
ia1=0.2;
ib1=0.7;
parax1=2.2;
ppara1=4/4;
% %%%%%%%%%%%%%%%%%%
a11(1)=-(parax1-0.35);
ii=1;
while  a11(ii)<(parax1-0.35)
    ii=ii+1;
    a11(ii)=a11(ii-1) + .7/sqrt(4*ppara1*a11(ii-1)^2+1);
end
b11=ppara1*a11.^2+.6;
b21=(ppara1*parax1^2)+.5:ib1:up-.3;
a21=parax1*ones(1,length(b21))+0.001*rand(1,length(b21));
%a31=parax1:-ia1:-parax1;
%b31=(up)*ones(1,length(a31));
b41=b21(end):-ib1:(ppara1*parax1^2)+.5;
a41=-parax1*ones(1,length(b41))+0.001*rand(1,length(b41));

xv1=[a11';a21';a41'];
yv1=[b11';b21';b41'];
%%%%%%%%%%%%%%%%%%% Segunda fila en forma parabolica
ia2=0.2;
ib2=0.9;
parax2=3;
ppara2=4/9;
a12(1)=-(parax2-0.11);
a12(1)=0.09;
ii=1;
while  a12(ii)<(parax2-.3)
    ii=ii+1;
    a12(ii)=a12(ii-1) + 0.9/sqrt(4*ppara2*a12(ii-1)^2+1);
end
b12=ppara2*a12.^2;
b22=(ppara2*parax2^2)+.1:ib2:up-.3;
a22=parax2*ones(1,length(b22));
%a32=parax2:-ia2:-parax2;
%b32=(up)*ones(1,length(a32));
b42=b22(end):-ib2:(ppara2*parax2^2)+.1;
a42=-parax2*ones(1,length(b42));

a12r=[-a12';a12'];
b12r=[b12';b12'];
a12=a12r';
b12=b12r';
xv2=[a12';a22';a42'];
yv2=[b12';b22';b42'];
%%%%%%%%%%%%%%%%%%%%%%%%%% Tercer fila en forma parabolica,,,esta ya no
%%%%%%%%%%%%%%%%%%%%%%%%%% presenta forma cerrada en las celdas, solo es
%%%%%%%%%%%%%%%%%%%%%%%%%% para contener. frontera ficticia.
ia3=0.2;
ib3=0.4;
parax3=4;
ppara3=5/16;
a13(1)=-(parax3-0.17);
ii=1;
while  a13(ii)<(parax3-.2)
    ii=ii+1;
    a13(ii)=a13(ii-1) + 0.4/sqrt(4*ppara3*a13(ii-1)^2+1);
end
b13=ppara3*a13.^2-1;
b23=(ppara3*parax3^2)-.8:ib3:up-.1;
a23=parax3*ones(1,length(b23));
a33=parax3:-ia3:-parax3;
b33=(up+.3)*ones(1,length(a33));
%b43=b23(end):-ib3:(ppara3*parax3^2)-.8;
b43=b23;
a43=-parax3*ones(1,length(b43));

xv3=[a13';a23';a33';a43'];
yv3=[b13';b23';b33';b43'];
xv3f1=[a13';a23';a43'];
yv3f1=[b13';b23';b43'];

xv3f=[a33'];
yv3f=[b33'];
%%%%%%%%%%%%%%%%%%%%%%%%%

vb22=zeros(1,length(b22));
vb22=zeros(1,length(b22));
vb42=zeros(1,length(b42));
vb42=zeros(1,length(b42));

xv=[xv2;xv3];%a33'];
yv=[yv2;yv3];%b33'];

x(1)=-0.15;
x(2)=0.15;
y(1)=1.5*ppara+.0001;
y(2)=1.5*ppara-.0001;
xo1=x(1);
xo2=x(2);
yo1=y(1);
yo2=y(2);

x=[x;xv1];  % se unen todos los puntos, tanto los de la frontera como los de las celulas
y=[y;yv1];

xori=length(x);
xorin=xori+length(xv2); % condicion sobre los puntos de internos del dominio
xoric(cuentas)=xori;
vx=zeros(xori,1);
vy=vx;
Ep=vx;
Ec=vx;
c=vx;
r=vx;

x=[x;xv];  % se unen todos los puntos, tanto los de la frontera como los de las celulas
y=[y;yv];

%condiciones en la frontera
vx=[vx;zeros((length(xv)),1)];  % velocidad
vy=[vy;zeros((length(yv)),1)];
Ep=[Ep;zeros((length(xv)),1)];  %potencial
Ec=[Ec;zeros((length(xv)),1)];  %energia cinetica
c=[c;zeros((length(xv)),1)];   % auxina
r=[r;zeros((length(xv)),1)];    %reloj ciclinas
%%


tri= delaunay(x,y);       % se genera una triangulacion de Delaunay para el voronoi
tria=[tri tri(:,1)];
neib=neibs(x,y,length(x),tria);  % funcion para encontrar los indices y luego los vecinos
ABpinnorb                 % Programa para saber quienes son los puntos de la frontera despues de la triangulacion
[vert celdas]=voronoin([x,y],{'Qbb','Qz'});
%%    Condiciones sobre las auxinas y el reloj del ciclo celular
c=0*x;
r=0*x;
r(1:xori+length(xv2))=15*rand(xori+length(xv2),1);
c(1:xori)=.1*rand(xori,1);
c(xori+1:xori+1*length(xv2))=0.2*max(c);
c=c/sum(c);
c(1)=max(c);  % estas celulas representan el centro quiescente
c(2)=c(1);
%c=c/sum(c);


lapmat=lapmr(x,y,xori+length(xv2),tria,neib); %calcula el laplaciana en la red
[gradcx,gradcy]=trigradient(x,y,c,tri);       %calcula el gradiente de la auxina
[gradvx,gradvy]=trigradient(x,y,Ep,tri);      %calcula el gradiente del potencial
lapc=lapmat*c(1:xori+length(xv2));

gradcx(JJpf)=0;
gradcy(JJpf)=0;
gradvx(JJpf)=0;
gradvy(JJpf)=0;

areacero=polyarea(xv3,yv3)/(xori+length(xv2));    %area cero promedio
areacero1=0.4;          %area cero promedio alternativa
prr=0;            % contador de ciclos


contfig=0;     %contador de inicio de las figuras
ABpinnora      % encuentra los puntos que estan en la epidermis
JJ=JJp;        % indices de los puntos de la epidermis
JJs=xori+1:1:xori+18;


for prr=1:Nf %(240==10 dias si mo=20)  %para ciclos de visualizacion
    
    flux=[];     % guardan el flujo de auxina
    fluxt=[];
    %**********************************************************
    for m=1:mo %tiempos
        % [num2str(prr) '.' num2str(m) ' reprod = ' num2str(rep) ] % para
        % visualizar la variable
        
        
        %%    reproduccion celulas internas
        
        ar=find(r(1:xori).*sqrt(c(1:xori))/gamma>=1);  % discriminate para saber si hay reproduccion
        
        if ~isempty(ar)  %&& dfup==0
            reptotroot       % programa para efectuar la reproduccion
            rcon=rcon+1;
            rxad=[rxad;xad];   % posiciones ageagadas
            ryad=[ryad;yad];
            %       %%
            
            %figure(11) %%% reproducciones
            %clf
            %W=histc(ryad,edges);
            %hold on
            %%hist(ryad,100)
            %bar(edges,W,'histc')
            %axis([0 max(y) 0 12])
            %plot(y(1:xori),Ep(1:xori)/max(Ep)*8,'.r')
            %%axis equal
            %hold off
            %R=getframe;
            %RN(rcon)=R;
            vxp=vx;
            vyp=vy;
            tri= delaunay(x,y);
            tria=[tri tri(:,1)];
            neib=neibs(x,y,length(x),tria);
            
            [vert celdas]=voronoin([x,y],{'Qbb','Qz'});
            
            
            ABpinnora
            JJ=JJp;
            JJs=xori+1:1:xori+18;
            coaux=coa;
        else
            % r(ar)=0;
            
        end
        
        
        %%    reproduccion de celulas epidermicas
        
        
        ar=find(r(xori+1:xori+length(xv2)).*sqrt(c(xori+1:xori+length(xv2)))/gamma>=1);
        %r(ar)=0;
        if ~isempty(ar) && dfup<=0.0001
            r(ar+xori)=0.5;
        end
        
        
        AREA=zeros(xori+length(xv2),1);
        
        centroxx=zeros(length(x),1);
        centroyy=centroxx;
        %Ep=zeros((xori+length(xv2)+length(xv3)),1);
        
        
        testjun   % funcion para aumentar el dominio
        %%
        tri= delaunay(x,y);
        tria=[tri tri(:,1)];
        neib=neibs(x,y,length(x),tria);
        
        [vert celdas]=voronoin([x,y],{'Qbb','Qz'});
        
        
        ABpinnora
        ABpinnorb
        JJ=JJp;
        
        for j=1:xori+length(xv2)                 %para celulas
            [centrox centroy]=centroid(vert(celdas{j},1),vert(celdas{j},2));  %calcula centroides
            [amasnj,amasnnxj,amasnnyj]=subarea1cil(x,y,vert,celdas,tria,j,xori,neib);  % calcula las ctes de fuerza derivadas del potencial
            AREA(j)=polyarea(vert(celdas{j},1),vert(celdas{j},2));
            %amasnj=AREA(j);
            centroxx(j)=centrox;
            centroyy(j)=centroy;
            xj1=(x(j)-centrox);
            yj1=(y(j)-centroy);
            
            if isempty(find(j==JJ(:), 1))==1  % puntos de interior
                
                vx(j)=vx(j)+(kvol*(AREA(j)-areacero)*amasnnxj-kcen*xj1-kx*vx(j))*dt;
                vy(j)=vy(j)+(kvol*(AREA(j)-areacero)*amasnnyj-kcen*yj1-ky*vy(j))*dt;
                
                x(j)=x(j)+ vx(j)*dt;
                y(j)=y(j)+ vy(j)*dt;
                
                Ep(j)=kvol*(AREA(j)-areacero)^2/2+kcen*(xj1^2+yj1^2)/2;
                Ec(j)=(vx(j).^2+vy(j).^2)/2;
            end
            
            if isempty(find(j==JJ(:), 1))==0 % puntos de la epidermis
                
                [amasnjf,amasnnxjf,amasnnyjf,amasnnzjf]=subarea1cil2(x,y,vert,celdas,tria,j,xori,neib,JJ);
                
                vx(j)=vx(j)+(kbs*(AREA(j)-areacero1)*(amasnnxj-amasnnxjf)-kecen*xj1)*dt-kbsrx*vx(j)*dt;
                vy(j)=vy(j)+(kbs*(AREA(j)-areacero1)*(amasnnyj-amasnnyjf)-kecen*yj1)*dt-kbsry*vy(j)*dt;
                
                x(j)=x(j)+vx(j)*dt;
                y(j)=y(j)+vy(j)*dt;
                
                Ep(j)=kbs*(AREA(j)-areacero1)^2/2+kecen*(xj1^2+yj1^2)/2;
                Ec(j)=(vx(j).^2+vy(j).^2)/2;
            end
            
            if isempty(find(j==JJs(:), 1))==0   % para los puntos restantes del dominio
                vx(j)=0;
                vy(j)=0;
            end
            
        end
        
        %           x(1)=xo1;
        %           x(2)=xo2;      %para fijar el centro quiescente.
        %           y(1)=yo1;      %si se activan, comentar lineas 475-476 y 491
        %           y(2)=yo2;      % y descomentar lineas 477-490
        x(JJ)=centroxx(JJ); % fija el centroide de las celdas epidermales
        y(JJ)=centroyy(JJ);
        %testjun
        %%
        tri= delaunay(x,y);
        tria=[tri tri(:,1)];
        neib=neibs(x,y,length(x),tria);
        
        [vert celdas]=voronoin([x,y],{'Qbb','Qz'});
        
        
        
        lapmat=lapmr(x,y,xori+length(xv2),tria,neib); %calcula el laplaciana en la red
        lapc=lapmat*c(1:xori+length(xv2));
        
        [gradvx,gradvy]=trigradient(x,y,Ep,tri);
        [gradcx,gradcy]=trigradient(x,y,c,tri);
        
        
        %%*****************************************************************************
        %******************************************************************************
        %        Condiciones de frontera para el flujo de auxina
        
        ABpinnora
        ABpinnorb
        JJ=JJp;
        
        lastpr=find(y(JJ) & x(JJ)>0);
        lastpr1=find(y(JJ(lastpr))==max(y(JJ(lastpr))) & x(JJ(lastpr))>0);
        lastpl=find(y(JJ) & x(JJ)<0);
        lastpl1=find(y(JJ(lastpl))==max(y(JJ(lastpl))) & x(JJ(lastpl))<0);
        
        %JJpf=[JJpf JJ(lastpr(lastpr1)) JJ(lastpl(lastpl1))];
        
        cp=sum(c(1:xori+length(xv2)))/(xori+length(xv2));
        c(JJpf)=0*cp;
        for jir=1:length(JJpf)    %Celulas de arriba
            nbsup=neib(JJpf(jir));
            ct=cell2mat(nbsup);
            c1=ct(ct(:) > xori+length(xv2));
            Ep(c1)=Ep(JJpf(jir));
        end
        
        a1=JJ(lastpl(lastpl1));   % celuas esquina izquierda
        %Ep(a1)=0*Ep(a1);
        c(a1)=0;
        a2=JJ(lastpr(lastpr1));    % celulas esquina derecha
        %Ep(a2)=0*Ep(a2);
        c(a2)=0;
        %crr=c;
        
        
        for jir1=1:length(JJ) % celulas epidermicas
            nbsup=neib(JJ(jir1));
            ct=cell2mat(nbsup);
            c1=ct(ct(:) > xori+length(xv2));
            rt=[sum(x(c1))-length(c1)*x(JJ(jir1)),sum(y(c1))-length(c1)*y(JJ(jir1))];
            nt=-rt/sqrt(rt(1)^2+rt(2)^2);
            
            gr=[gradcx(JJ(jir1)),gradcy(JJ(jir1))];
            png=nt(1)*(gr(1))+nt(2)*(gr(2));
            grv=gr-[png*(nt(1)),png*(nt(2))];
            gradcx(JJ(jir1))=grv(1);
            gradcy(JJ(jir1))=grv(2);
            
            
        end
  
        %********************************************************************************
        %********************************************************************************
        
        %% Para el transporte de auxina
        
        Auxin=c;  %funcion auxiliar
        crr=0;
        for j=1:xori+length(xv2)                
            simpflujo    % funcion que calcula el flujo de auxina por vecino
            Auxin(j)=Auxin(j)+dt*(am)*alfa-sigma*Auxin(j)*dt;
            crr=crr+(am)*alfa;
        end 
        
        c=Auxin;
        
        c(xori+length(xv2)+1:end)=zeros(length(xv3),1);
        Ep(xori+length(xv2)+1:end)=zeros(length(xv3),1);
 
        %% Para las ciclinas
        for j=1:xori+length(xv2)
            if c(j)>cmin && isempty(find(j==JJ(:), 1))==1
                r(j)=r(j)+dt*beta/sqrt(abs(c(j)-cmin));
            elseif c(j)>cmin && isempty(find(j==JJ(:), 1))==0
                r(j)=r(j)+dt*betae/sqrt(abs(c(j)-cmin));
            else
                r(j)=r(j);
            end
        end
        
    
        
        Areazero(prr)=areacero;
        areacero=(sum(AREA)-sum(AREA(JJ)))/(length(AREA)-length(JJ));
        areacero1=.4;
        mtcp(prr)=sum(c(1:xori+length(xv2)));
        mtup(prr)=up;
        
        
        c(find(c<0))=0;
        c(find(c>1))=1;
        %      c(1)=max(c);
        %      c(2)=c(1);
        
        %      c=(c-min(c))/(max(c)-min(c));
        %      c=c/sum(c);
        %c=(c-min(c))/(max(c)-min(c));
        
        %  c=c/sum(c);
        
        %        cpre=(max(c)-c(1));
        %        c(1)=max(c)+cpre;
        %        c(1:(xori+length(xv2)))=(c(1:(xori+length(xv2)))-cpre/(max(c)-min(c)));
        %        c=c/(sum(c)-cpre);
        %      sigma=crr/sum(c);
        sigma=0;
        
        
    end  % tiempos
    
    x5=xv3;
    y5=yv3;
    
    %para salvar
    mtcp(prr)=sum(c(1:xori+length(xv2)));
    mtEc(prr)=sigma;
    mtEp(prr)=sum(Ep(1:xori+length(xv2)))/(xori+length(xv2));
    contfig=contfig+1;
    %save(['z1save-' num2str(prr)])
    
    
    FigRBc   % figuras 
end   %*** ciclos

