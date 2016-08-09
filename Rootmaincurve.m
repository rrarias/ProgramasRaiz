%Programa para relajar y hacer crecimiento en funcion del area
%****************************************************************
%********* condiciones iniciales ********************************
 %RAIZ CURVA 
clear all
% XXXXXXXXXXXXXX    DATOS   XXXXXXXXXXXXXXXXXXXXXXXXX
%%%%%%%%%% constantes de fuerza
kx=3;
ky=kx;
kbsrx=3;     %constante de disipacion frontera
kbsry=3;
kvol=80;
kcen=80;
kbs=120;
kecen=120;
rep=0;
edges=0:.12:30;

cmin=.0005;   % el cutoff de la concentracion de auxinas.
alfa=8;%     escala de tiempo para las auxinas
alfac=alfa;
Do=10.00;
beta=15; %     escala de tiempo de las ciclinas
betae=5;
gamma=220;%     ajusta el tiempo del ciclo de ciclinas
rcon=1;
sigma=0.3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=0.01;
nc=10;     %     numero inicial de celulas, da la densidad inicial
parax=2;   %    ancho del dominio
ppara=1;     %  coefciente de la parabola
d=8;      %el valor min es 2
up=ppara*parax^2+d;    % altura inicial del cilindro
%up=4;
bo=(up-d)/2; %    posicion del nicho
Nf=200;

up1=0;          %altura mas incremento
uax=d;          %incremento de altura
delta=3;        %division vertical
deltaup=delta;   % tolerancia del area
deltadown=delta;
ia=.2;  % intervalo puntos parabola
ib=.3;    % intervalo puntos cilindro

mo=5;    % numero de iteraciones por imagen

%indice=0;  % indice del nicho (cambia)
invc=0;   % contador del incho
timequies=10; %     tasa de reprod del nicho
xad=[];
cuentas=1;   %  contador del numero de reproducciones por ciclo
pp=0.;  % porcentaje de area de la parabola creciendo pp=0 solo el rectangulo crece
rxad=[];  %contar las coordenadas de las reproducciones
ryad=[];
sidediv=[];
coa=0;
coar=0;   %medidor de area
dfup=0;   %condiciones iniciales para aumentar up
dsourup=0;
addcoa1=0;
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ff=1;


%% ********************Boundaries**************************************
%boundary
ia1=0.2;
ib1=0.9;%0.7
parax1=2.2;
ppara1=4/4;

ia2=0.7;
ib2=1.2;%0.9
parax2=3;
ppara2=4/9;


ia3=0.3;
ib3=0.7;%0.5
parax3=4;
ppara3=5/16;

raddw1=2.3; %radio de los semicirculos
updw=-1.5; %altura hacia abajo de la raiz
raddw2=raddw1-(parax2-parax1); %radio del circulo inferio
raddw3=raddw1-(parax3-parax1); %radio del circulo superior
npoinu=15;    %numero de puntos circulo superior
npoind=6;     %numero de punto circulo inferior
npoinu3=40;   %numero de puntos circulo superior tercera fila


%save valrootinit
%*********************Primera fila en forma parabolica

% %%%%%%%%%%%%%%%%%%
a11(1)=-(parax1-0.25);
ii=1;
while  a11(ii)<(parax1-0.25)
    ii=ii+1;
    a11(ii)=a11(ii-1) + 0.9/sqrt(4*ppara1*a11(ii-1)^2+1); %0.7
end
b11=ppara1*a11.^2+.6;
b21=(ppara1*parax1^2)+.5:ib1:up-.3;
a21=parax1*ones(1,length(b21))+0.001*rand(1,length(b21));


iil=0;
for ii=1:npoind
    iil=1+iil;
    radians=pi*(ii)/npoind+0.01*rand();
    a31(npoind+1-iil)= parax1+raddw1 +raddw1*cos(radians);
    b31(npoind+1-iil)= b21(end)+ib1 +raddw1*sin(radians);
end


b41=b31(end)-ib1-.001:-ib1:updw+.3;
a41=a31(end)*ones(1,length(b41))+0.001*rand(1,length(b41));

b51=b41(end):ib1:b31(end)-ib1-0.001;
a51=(3*parax1+2*raddw1)*ones(1,length(b51))+0.001*rand(1,length(b51));


iil=0;
for ii=1:npoinu
    iil=1+iil;
    radians=pi*(ii)/npoinu+0.01*rand();
    a61(iil)= parax1+raddw1 +(2*parax1+raddw1)*cos(radians);
    b61(iil)= b21(end)+ib1 +(2*parax1+raddw1)*sin(radians);
end

b71=b61(end)-ib1:-ib1:(ppara1*parax1^2)+.5;
a71=-parax1*ones(1,length(b71))+0.001*rand(1,length(b71));


xv1=[a11';a21';a31';a41';a51';a61';a71'];
yv1=[b11';b21';b31';b41';b51';b61';b71'];



%%%%%%%%%%%%%%%%%%% Segunda fila en forma parabolica

a12(1)=-(parax2-0.11);
%a12(1)=0.09;
ii=1;
while  a12(ii)<(parax2-.33)
    ii=ii+1;
    a12(ii)=a12(ii-1) + 1.1/sqrt(4*ppara2*a12(ii-1)^2+1);%0.8
end
b12=ppara2*a12.^2;
b22=(ppara2*parax2^2)+.43:ib2:up-.1;
a22=parax2*ones(1,length(b22));

%a12r=[-a12';a12'];
%b12r=[b12';b12'];
%a12=a12r';
%b12=b12r';


iil=0;
for ii=1:npoind
    iil=1+iil;
    radians=pi*(ii)/npoind +0.01*rand();
    a32(npoind+1-iil)= parax2+raddw2 +raddw2*cos(radians);
    b32(npoind+1-iil)= b22(end)+ib2 +raddw2*sin(radians);
end


b42=b32(end)-ib2:-ib2:updw+.3;
a42=a32(end)*ones(1,length(b42))+0.001*rand(1,length(b42));

b52=b42(end):ib2:b32(end)-ib2;
a52=(3*parax2+2*raddw2)*ones(1,length(b52))+0.001*rand(1,length(b52));


iil=0;
for ii=1:npoinu
    iil=1+iil;
    radians=pi*(ii)/npoinu-0.01*rand();
    a62(iil)= parax2+raddw2 +(2*parax2+raddw2)*cos(radians);
    b62(iil)= b22(end)+ib2 +(2*parax2+raddw2)*sin(radians);
end

b72=b62(end)-ib2:-ib2:(ppara2*parax2^2)+.43;
a72=-parax2*ones(1,length(b72))+0.001*rand(1,length(b72));

a82=a32(end):ia2:(3*parax2+2*raddw2);
b82=(updw-.3)*ones(1,length(a82));

%a32=parax2:-ia2:-parax2;
%b32=(up)*ones(1,length(a32));
%b42=b22(end):-ib2:(ppara2*parax2^2)+.1;
%a42=-parax2*ones(1,length(b42));
xv2=[a12';a22';a32';a42';a52';a62';a72'];
yv2=[b12';b22';b32';b42';b52';b62';b72'];

%xv2=[a12';a22';a42'];
%yv2=[b12';b22';b42'];
%%%%%%%%%%%%%%%%%%%%%%%%%% Tercer fila en forma parabolica,,,esta ya no
%%%%%%%%%%%%%%%%%%%%%%%%%% presenta forma cerrada en las celdas, solo es
%%%%%%%%%%%%%%%%%%%%%%%%%% para contener. frontera ficticia.

a13(1)=-(parax3-0.1);
ii=1;
while  a13(ii)<(parax3-.16)
    ii=ii+1;
    a13(ii)=a13(ii-1) + 0.7/sqrt(4*ppara3*a13(ii-1)^2+1);%0.4
end
b13=ppara3*a13.^2-1;
b23=(ppara3*parax3^2)-.8:ib3:up-.1;
a23=parax3*ones(1,length(b23));


iil=0;
for ii=1:npoind
    iil=1+iil;
    radians=pi*(ii)/npoind+0.00001*rand();
    a33(npoind+1-iil)= parax3+raddw3 +raddw3*cos(radians);
    b33(npoind+1-iil)= b23(end)+ib3 +raddw3*sin(radians);
end


b43=b33(end)-ib2:-ib3:updw+.3;
a43=a33(end)*ones(1,length(b43))+0.001*rand(1,length(b43));

b53=b43(end):ib3:b33(end);
a53=(3*parax3+2*raddw3)*ones(1,length(b53))+0.001*rand(1,length(b53));


iil=0;
for ii=1:npoinu3
    iil=1+iil;
    radians=pi*(ii)/npoinu3+0.00001*rand();
    a63(iil)= parax3+raddw3 +(2*parax3+raddw3)*cos(radians);
    b63(iil)= b23(end)+ib3 +(2*parax3+raddw3)*sin(radians);
end

b73=b63(end)-ib2:-ib3:(ppara3*parax3^2)-.7;
a73=-parax3*ones(1,length(b73))+0.001*rand(1,length(b73));


a83=a33(end):ia3:(3*parax3+2*raddw3);
b83=(updw-.3)*ones(1,length(a83));


xv3=[a13';a23';a33';a43';a83';a53';a63';a73'];
yv3=[b13';b23';b33';b43';b83';b53';b63';b73'];

xv3f1=[a13';a23';a33';a43';a53';a63';a73'];
yv3f1=[b13';b23';b33';b43';b53';b63';b73'];

xv3f=[a83'];
yv3f=[b83'];

vb42=zeros(1,length(b42));
vb52=zeros(1,length(b52));



xv=[xv2;xv3];
yv=[yv2;yv3];




%% **************** condiciones iniciales celdas interiores************************************************

xa1=.8*(parax-ib)*(rand(3.5*nc,1)-0.5);

ya1=ppara*xa1.^2+((up)-(ppara*xa1.^2+5*ia)).*rand(3.5*nc,1)+5*ia;


for ii=1:3.0*nc
    radians=pi*rand();
    ax1(ii)= parax3+raddw3 +(((3*parax1+raddw1-4*ib1))*(rand()))*cos(radians);
    bx1(ii)= b21(end)+raddw1+1.5*ib1 + (((3*parax1+raddw1-6*ib1))*(rand()))*sin(radians);
end

ax2=(3*parax1+2*raddw1-5*ia1)-(a13(end)-5*ia1)*rand(3.0*nc,1);
bx2=(up-updw)*rand(3.0*nc,1)+updw+0.5;


x=[xa1;ax1';ax2];
y=[ya1;bx1';bx2];

%d=8 xa1=1.2*nc; ax1=0.7*nc ax2=1.6*nc
%d=4 xa1=1.0*nc; ax1=0.7*nc ax2=1.2*nc
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%XXXXX   posicion de las 2 troncales XXXXXXXXXXXXX
x(1)=-0.25;
x(2)=0.25;
y(1)=2.5*ppara+.0001;
y(2)=2.5*ppara-.0001;
xo1=x(1);
xo2=x(2);
yo1=y(1);
yo2=y(2);

x=[x;xv1];
y=[y;yv1];

xori=length(x);
xorin=xori+length(xv2);
xoric(cuentas)=xori;
vx=zeros(xori,1);
vy=vx;
Ep=vx;
Ec=vx;
c=vx;
%r=ones(xori,1);

x=[x;xv];
y=[y;yv];


vx=[vx;zeros((length(xv)),1)];
vy=[vy;zeros((length(yv)),1)];
Ep=[Ep;zeros((length(xv)),1)]; %potencial en la frontera
Ec=[Ec;zeros((length(xv)),1)];
c=[c;zeros((length(xv)),1)];

figure(1)
voronoi(x,y,'g')
hold on
plot(x,y,'b.')
plot(ax1,bx1,'ro')
hold off
%axis([-6 14 -11 18])
axis equal


%%
prr=1;
tri= delaunay(x,y);
tria=[tri tri(:,1)];
neib=neibs(x,y,length(x),tria);

[vert celdas]=voronoin([x,y],{'Qbb','Qz'});

%AB     %      encuentra las celdas de la orilla
c=0*x;
r=0*x;
%c(1:xori+length(xv2))=1*rand(xori+length(xv2),1);
r(1:xori)=.2*rand(xori,1);
c(1:xori)=1*rand(xori,1);
c(xori+1:xori+1*length(xv2))=0.1*max(c);


c(1)=max(c);
c(2)=c(1);
%c=c/sum(c);

%load perl2
r=0*x;
r(1:xori)=2*rand(xori,1);

lapmat=lapmr(x,y,xori+length(xv2),tria,neib); %calcula el laplaciana en la red
[gradcx,gradcy]=trigradient(x,y,c,tri);
[gradvx,gradvy]=trigradient(x,y,Ep,tri);
lapc=lapmat*c(1:xori+length(xv2));

areacerot=polyarea(xv2,yv2)/(xori);
areacero=areacerot;
areacero1=0.4;

ABpinnora
ABpinnorb
JJ=JJp;
%%
for prr=1:Nf %(8 dias)  %para ciclos de visualizacion
    
    flux=[];
    fluxt=[];
    %**********************************************************
    for m=1:mo %tiempos
        [num2str(prr) '.' num2str(m) '  invc = ' num2str(invc)]
        vipol=[];
        %    reproduction  
        
        %ar=find(r(1:xori).*sqrt(c(1:xori))/gamma>=1);
%         ar=find(r(1:xori)./sqrt(gamma)>=1);
%         if ~isempty(ar)  %&& dfup==0
%             reptotroot
%             rcon=rcon+1;
%             rxad=[rxad;xad];
%             ryad=[ryad;yad];
%             %       %%
%             
%             %figure(11) %%% reproducciones
%             %clf
%             %W=histc(ryad,edges);
%             %hold on
%             %%hist(ryad,100)
%             %bar(edges,W,'histc')
%             %axis([0 max(y) 0 12])
%             %plot(y(1:xori),Ep(1:xori)/max(Ep)*8,'.r')
%             %%axis equal
%             %hold off
%             %R=getframe;
%             %RN(rcon)=R;
%             vxp=vx;
%             vyp=vy;
%             tri= delaunay(x,y);
%             tria=[tri tri(:,1)];
%             neib=neibs(x,y,length(x),tria);
%             
%             [vert celdas]=voronoin([x,y],{'Qbb','Qz'});
%             
%             
%             ABpinnora
%             JJ=JJp;
%             JJs=xori+1:1:xori+18;
%             coaux=coa;
%             r(ar)=0;
%             
%         end
%         
%         
        %%    reproduction epidermal cells
        
        
        ar=find(r(xori+1:xori+length(xv2)).*sqrt(c(xori+1:xori+length(xv2)))/gamma>=1);
        %r(ar)=0;
        if ~isempty(ar) && dfup<=0.0001
            r(ar+xori)=0.5;
        end
        
        AREA=zeros(xori+length(xv2),1);
        
        centroxx=zeros(length(x),1);
        centroyy=centroxx;
        
        curvetest   % mide el tama?o de la frontera
        %%
        tri= delaunay(x,y);
        tria=[tri tri(:,1)];
        neib=neibs(x,y,length(x),tria);
        
        [vert celdas]=voronoin([x,y],{'Qbb','Qz'});
        
        
        ABpinnora
        ABpinnorb
        JJ=JJp;
        
        
        for j=1:xori+length(xv2)                 %for particule
            [centrox centroy]=centroid(vert(celdas{j},1),vert(celdas{j},2));
            [amasnj,amasnnxj,amasnnyj]=subarea1cil(x,y,vert,celdas,tria,j,xori,neib);
            AREA(j)=polyarea(vert(celdas{j},1),vert(celdas{j},2));
            %amasnj=AREA(j);
            centroxx(j)=centrox;
            centroyy(j)=centroy;
            xj1=(x(j)-centrox);
            yj1=(y(j)-centroy);
            
            if isempty(find(j==JJ(:), 1))==1
                
                
                vx(j)=vx(j)+(kvol*abs(AREA(j)-areacero)*amasnnxj-kcen*xj1-kx*vx(j))*dt;
                vy(j)=vy(j)+(kvol*abs(AREA(j)-areacero)*amasnnyj-kcen*yj1-ky*vy(j))*dt;
                
                x(j)=x(j)+ vx(j)*dt;
                y(j)=y(j)+ vy(j)*dt;
                
                Ep(j)=kvol*(AREA(j)-areacero)^2/2+kcen*(xj1^2+yj1^2)/2;
                Ec(j)=(vx(j).^2+vy(j).^2)/2;
            end
            if isempty(find(j==JJ(:), 1))==0
                
                [amasnjf,amasnnxjf,amasnnyjf,amasnnzjf]=subarea1cil2(x,y,vert,celdas,tria,j,xori,neib,JJ);
                
                vx(j)=vx(j)+(kbs*abs(AREA(j)-areacero1)*(amasnnxj-amasnnxjf)-kecen*xj1)*dt-kbsrx*vx(j)*dt;
                vy(j)=vy(j)+(kbs*abs(AREA(j)-areacero1)*(amasnnyj-amasnnyjf)-kecen*yj1)*dt-kbsry*vy(j)*dt;
                
                x(j)=x(j)+vx(j)*dt;
                y(j)=y(j)+vy(j)*dt;
                
                Ep(j)=kbs*(AREA(j)-areacero1)^2/2+kecen*(xj1^2+yj1^2)/2;
                Ec(j)=(vx(j).^2+vy(j).^2)/2;
            end
            
            %if isempty(find(j==JJs(:), 1))==0
            %       %x(j)=centroxx(j);
            %       %y(j)=centroyy(j);
            %       vx(j)=0;
            %       vy(j)=0;
            %end
            
        end
        
        %x(1)=xo1;
        %x(2)=xo2;
        %y(1)=yo1;
        %y(2)=yo2;
        x(JJ)=centroxx(JJ); % fija el centroide de las celdas epidermales
        y(JJ)=centroyy(JJ);
        
        
        
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
        %        Boundary conditions for flux
        
        ABpinnora
        ABpinnorb
        JJ=JJp;
        
        lastpr=find(y(JJ) & x(JJ)>9);
        lastpr1=find(y(JJ(lastpr))==min(y(JJ(lastpr))) & x(JJ(lastpr))>9);
        lastpl=find(y(JJ) & x(JJ)<9 & x(JJ)>4);
        lastpl1=find(y(JJ(lastpl))==min(y(JJ(lastpl))) & x(JJ(lastpl))<9 & x(JJ(lastpl))>4);
        
        %JJpf=[JJpf JJ(lastpr(lastpr1)) JJ(lastpl(lastpl1))];
        
        cp=sum(c(1:xori+length(xv2)))/(xori+length(xv2));
        c(JJpf)=0*cp;
        for jir=1:length(JJpf)    %top cells
            nbsup=neib(JJpf(jir));
            ct=cell2mat(nbsup);
            c1=ct(ct(:) > xori+length(xv2));
            Ep(c1)=Ep(JJpf(jir));
        end
        
        a1=JJ(lastpl(lastpl1));   % corner cells left
        %Ep(a1)=0*Ep(a1);
        c(a1)=0;
        a2=JJ(lastpr(lastpr1));    % corner cells right
        %Ep(a2)=0*Ep(a2);
        c(a2)=0;
        %crr=c;
        
        
        for jir1=1:length(JJ) % epidermal cells
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
        %%
        %Ep=Ep/max(Ep);
        
        
        %% For auxin transport
        
        Auxin=c;
        crr=0;
        %flujo
        for j=1:xori+length(xv2)                 %for particule
            simpflujo
            Auxin(j)=Auxin(j)+dt*(am)*alfa-sigma*Auxin(j)*dt;
            crr=crr+(am)*alfa;
        end
        
        
        c=Auxin;
       
        
        c(xori+length(xv2)+1:end)=zeros(length(xv3),1);
        Ep(xori+length(xv2)+1:end)=zeros(length(xv3),1);
       
        %% For cyclins
        for j=1:xori+length(xv2)
            if c(j)>cmin && isempty(find(j==JJ(:), 1))==1
                %r(j)=r(j)+dt*beta/sqrt(abs(c(j)-cmin));
                r(j)=r(j)+dt*beta*sqrt(abs(c(j)-cmin));
            elseif c(j)>cmin && isempty(find(j==JJ(:), 1))==0
                %r(j)=r(j)+dt*betae/sqrt(abs(c(j)-cmin));
                r(j)=r(j)+dt*betae*sqrt(abs(c(j)-cmin));
            else
                r(j)=r(j);
            end
        end
        
        %r(xori+1:xori+15)=0;
        %r(18)=0;
        
        Areazero(prr)=areacero;
        %areacero=(sum(AREA(:))-sum(AREA(JJ)))/(xori+length(xv2)-length(JJ));
        areacero=(sum(AREA)-sum(AREA(JJ)))/(length(AREA)-length(JJ));
        areacero1=sum(AREA(JJ))/length(AREA(JJ));%0.4
        mtcp(prr)=sum(c(1:xori+length(xv2)));
        mtup(prr)=up;
        
        
        c(find(c<0))=0;
        c(find(c>1))=1;
        %c(1)=max(c);
        %c(2)=c(1);
        
        %c=(c-min(c))/(max(c)-min(c));
        %c=c/sum(c);
        %c=(c-min(c))/(max(c)-min(c));
        
        %  c=c/sum(c);
        
        %   cpre=(max(c)-c(1));
        %   c(1)=max(c)+cpre;
        %   c(1:(xori+length(xv2)))=(c(1:(xori+length(xv2)))-cpre/(max(c)-min(c)));
        %   c=c/(sum(c)-cpre);
        sigma=crr/sum(c);
        sigma=0;
        
        
        

        
    end  % tiempos
    
    x5=xv3;
    y5=yv3;
    mtcp(prr)=sum(c(1:xori+length(xv2)));
    mtEc(prr)=sigma;
    mtEp(prr)=sum(Ep(1:xori+length(xv2)))/(xori+length(xv2));
    
    FigRBc  
    
    
end   %*** ciclos


