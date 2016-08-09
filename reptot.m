%**********para reproduccion cuando se pase del tamano**********


xad=[];
yad=[];
vxad=[];
vyad=[];
Epad=[];
Ecad=[];
cad=[];
rad=[];
coa=0;         %aumento de altura
mnc=0;         %contador de celdas
%%%%%%%%%%%%reproduccion del nicho %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if    mod(invc,timequies)==0
if  invc>timequies
    invc=0;
    rpro=sqrt(AREA(1)/pi)/2;
    %theta=pi*(rand(1,1)-.5)/2;
    xad=[xad;centroxx(1)-rpro];%*cos(theta)];
    yad=[yad;centroyy(1)];%+rpro*sin(theta)];
    
    vxad=[vxad;0];
    vyad=[vyad;0];
    Epad=[Epad;0.5*Ep(1)];
    Ecad=[Ecad;Ec(1)];
    cad=[cad;1.0*c(1)];
    c(1)=1.0*c(1);
    %rad=[rad;r(1)];
    Ep(1)=0.5*Ep(1);
    Ec(1)=0;
    coa=coa+AREA(1)/(2*parax)/2;
    
    rpro=sqrt(AREA(2)/pi)/2;
    %theta=pi*(rand(1,1)-.5);
    xad=[xad;centroxx(2)+rpro];%*cos(theta)];
    yad=[yad;centroyy(2)];%-rpro*sin(theta)];
    
    vxad=[vxad;0];
    vyad=[vyad;0];
    Epad=[Epad;0.5*Ep(2)];
    Ecad=[Ecad;0.5*Ec(2)];
    cad=[cad;1.0*c(2)];
    c(2)=1.0*c(2);
    %rad=[rad;r(2)];
    Ep(2)=0.5*Ep(2);
    Ec(2)=0;
    coa=coa+AREA(2)/(2*parax)/2;
   
    mnc=mnc+1;
    ip=1;
    monitor
    mnc=mnc+1;
    ip=2;
    monitor
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for jj = 1:length(ar)
    ip=ar(jj);
    
    if (ip~=1 && ip~=2) %& AREA(ip)>0.75*areacero
        
        rpro=sqrt(AREA(ip)/3.1416)/2;
        
        if  ip>boc %all(bocn==ip)            %reproduccion vertical arriba de boc
            xad=[xad; centroxx(ip)];
            yad=[yad; (centroyy(ip)+rpro)];
            y(ip)=centroyy(ip);
            vxad=[vxad;0.*vx(ip)];
            vyad=[vyad;sqrt(vx(ip)^2+vy(ip)^2)];
            Epad=[Epad;0.5*Ep(ip)];
            Ecad=[Ecad;Ec(ip)];
            cad=[cad;1.0*c(ip)];
            %rad=[rad;0];
            vx(ip)=0.*vx(ip);
            vy(ip)=0.*vy(ip);
            Ep(ip)=0.5*Ep(ip);
            Ec(ip)=0*Ec(ip);
            c(ip)=1.0*c(ip);
            %r(ip)=0;
            coa=coa+AREA(ip)/(2*parax)/2/deltaup;
            mnc=mnc+1;
            monitor
        end
        
        
        if ip<=boc%all(bocn~=ip)       %      reproduccion isotropica excluyendo el nicho
            theta=2*pi*rand;
            xad=[xad; (centroxx(ip)+2*rpro*cos(theta))];
            yad=[yad; (centroyy(ip)+2*rpro*sin(theta))];
            x(ip)=centroxx(ip)-0*rpro*cos(theta);
            y(ip)=centroyy(ip)-0*rpro*sin(theta);
            vxad=[vxad;-vx(ip)];
            vyad=[vyad;-vy(ip)];
            Epad=[Epad;0.5*Ep(ip)];
            Ecad=[Ecad;Ec(ip)];
            cad=[cad;1.0*c(ip)];
            %rad=[rad;0];
            vx(ip)=0.*vx(ip);
            vy(ip)=0.*vy(ip);
            Ep(ip)=0.5*Ep(ip);
            Ec(ip)=0*Ec(ip);
            c(ip)=1.0*c(ip);
            %r(ip)=0;
            coa=coa+AREA(ip)/(2*parax)/2/deltadown;
            mnc=mnc+1;
            monitor
        end
    end
    
end






%%%%%%%%conteo de reproduccion %%%%%%%%%%%%%%%
if ~isempty(xad)
    %xori=length(x);
    
    x=x(1:xori);
    y=y(1:xori);
    vx=vx(1:xori);
    vy=vy(1:xori);
    Ep=Ep(1:xori);
    Ec=Ec(1:xori);
    c=c(1:xori);
    r=r(1:xori);
    x=[x;xad];
    y=[y;yad];
    vx=[vx;vxad];
    vy=[vy;vyad];
    Ep=[Ep;Epad];
    Ec=[Ec;Ecad];
    c=[c;cad];
    r=[r;rad];
    %**********************************************************
    
    xori=length(x);
    cuentas=cuentas+1;
    xoric(cuentas)=xori-nc;
    nc=xori;
    x=[x;xv];
    y=[y;yv];
    vx=[vx;zeros((length(xv)),1)];
    vy=[vy;zeros((length(yv)),1)];
    Ep=[Ep;zeros((length(xv)),1)];
    Ec=[Ec;zeros((length(xv)),1)];
    c=[c;zeros((length(xv)),1)];
    r=[r;zeros((length(vx)),1)];
    
    invc=invc+2*up*parax*xoric(cuentas)/nc;
    
    xad=[];
    %crecio=crecio+1;
    
    %*******************************************************************
    %***********boundary con incrementos*********************
    
    
    %up1=xoric(cuentas)*sum(AREA(1:nc-xoric(cuentas)))/(nc-xoric(cuentas))/(2*parax)/delta;%+sqrt((abs(1/rho-1/rhoo(prr-1))*(2*parax)));
    
    %up1=coa;
    up1=coa/delta;
    
    
    
    ancho=pp*up1;
    up1=up1*(1-pp);
    
    parax=(parax^3+3*ancho)^(1/3);
    

    up=ppara*parax^2+d+up1;
    
        %up1=ua+up1;
    
    
    
    a1(1)=-(parax-0.05);
    ii=1;
    while  a1(ii)<=(parax-0.05)
        ii=ii+1;
        a1(ii)=a1(ii-1) + .05/sqrt(4*ppara*a1(ii-1)^2+.05);
    end
    b1=ppara*a1.^2;
    b2=(ppara*parax^2):ib:up-.1;
    a3=-parax*ones(1,length(b2));
    a2=-a3;
    b3=b2;
    a4=-parax:ia:parax;
    b4=(up)*ones(1,length(a4));
    %b4=up*ones(1,sqrt(ceil(length(a4))));
    xv=[a1';a2';a3';a4'];
    yv=[b1';b2';b3';b4'];
    x=x(1:nc);
    y=y(1:nc);
    vx=vx(1:nc);
    vy=vy(1:nc);
    Ep=Ep(1:nc);
    Ec=Ec(1:nc);
    c=c(1:nc);
    x=[x;xv];
    y=[y;yv];
    vx=[vx;zeros((length(xv)),1)];
    vy=[vy;zeros((length(yv)),1)];
    Ep=[Ep;zeros((length(xv)),1)];
    Ec=[Ec;zeros((length(xv)),1)];
    c=[c;zeros((length(xv)),1)];
end
ua=ua+up1;


%monitor de posiciones

% XT(1:length(x),cuentas)=x(1:end);
% YT(1:length(x),cuentas)=y(1:end);
% VXT(1:length(vx),cuentas)=vx(1:end);
% VYT(1:length(vy),cuentas)=vy(1:end);
% EPT(1:length(Ep),cuentas)=Ep(1:end);
% CT(1:length(c),cuentas)=c(1:end);
