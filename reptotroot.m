%**********para reproduccion cuando se pase del tamano**********

xad=[];
yad=[];
vxad=[];
vyad=[];
Epad=[];
Ecad=[];
cad=[];
rad=[];
AREAad=[];
coa=0;
%%%%%%%%%%%%reproduccion del nicho %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%if    mod(invc,timequies)==0
% if  invc>timequies
%     invc=0;
%     rpro=sqrt(AREA(1)/pi)/2;
%     theta=2*pi*rand;
%     xad=[xad;centroxx(1)+2*rpro*cos(theta)];
%     yad=[yad;centroyy(1)+2*rpro*sin(theta)];
%     
%     vxad=[vxad;0.5*vx(1)];
%     vyad=[vyad;0.5*vy(1)];
%     Epad=[Epad;0.5*Ep(1)];
%     Ecad=[Ecad;Ec(1)];
%     cad=[cad;c(1)];
%     rad=[rad;0];
%     AREAad=[AREAad;AREA(ip)];
%     Ep(1)=0.5*Ep(1);
%     Ec(1)=0;
%     vx(1)=0.5*vx(1);
%     vy(1)=0.5*vy(1);
%     coa=coa+AREA(1);
%     sidediv=[sidediv ip];
%     %---------------------------------------
%     rpro=sqrt(AREA(2)/pi)/2;
%     theta=2*pi*rand;
%     xad=[xad;centroxx(2)+2*rpro*cos(theta)];
%     yad=[yad;centroyy(2)+2*rpro*sin(theta)];
%     
%     vxad=[vxad;0.5*vx(2)];
%     vyad=[vyad;0.5*vy(2)];
%     Epad=[Epad;0.5*Ep(2)];
%     Ecad=[Ecad;0.5*Ec(2)];
%     cad=[cad;c(2)];
%     rad=[rad;0];
%     AREAad=[AREAad;AREA(ip)];
%     Ep(2)=0.5*Ep(2);
%     Ec(2)=0;
%     vx(2)=0.5*vx(2);
%     vy(2)=0.5*vy(2);
%     coa=coa+AREA(2);
%     sidediv=[sidediv ip];
%     r(1)=0;
%     r(2)=0;
%     
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for jj = 1:length(ar)
    ip=ar(jj);
    r(ip)=0;
    
    if isempty(find(ip==JJ(:), 1))==1% && AREA(ip)>0.75*areacero
      
    
         rpro=sqrt(AREA(ip)/3.1416)/2;

            theta=2*pi*rand;
            xad=[xad; (centroxx(ip)+rpro*cos(theta))];
            yad=[yad; (centroyy(ip)+rpro*sin(theta))];
            x(ip)=centroxx(ip)-rpro*cos(theta);
            y(ip)=centroyy(ip)-rpro*sin(theta);
            vxad=[vxad;0.5*vx(ip)];
            vyad=[vyad;0.5*vy(ip)];
            Epad=[Epad;0.5*Ep(ip)];
            Ecad=[Ecad;Ec(ip)];
            cad=[cad;c(ip)];
            rad=[rad;0];
            AREAad=[AREAad;AREA(ip)];
            vx(ip)=0.5*vx(ip);
            vy(ip)=0.5*vy(ip);
            Ep(ip)=0.5*Ep(ip);
            Ec(ip)=0*Ec(ip);
            c(ip)=c(ip);
            coa=coa+AREA(ip);
            sidediv=[sidediv ip];
    else
        
         rpro=sqrt(AREA(ip)/3.1416)/2;

            
            xad=[xad; (centroxx(ip))];
            yad=[yad; (centroyy(ip)+rpro)];
            x(ip)=centroxx(ip);
            y(ip)=centroyy(ip)-rpro;
            vxad=[vxad;0.0*vx(ip)];
            vyad=[vyad;1.0*vy(ip)];
            Epad=[Epad;0.5*Ep(ip)];
            Ecad=[Ecad;Ec(ip)];
            cad=[cad;c(ip)];
            rad=[rad;0];
            AREAad=[AREAad;AREA(ip)];
            vx(ip)=0.0*vx(ip);
            vy(ip)=0.0*vy(ip);
            Ep(ip)=0.5*Ep(ip);
            Ec(ip)=0*Ec(ip);
            c(ip)=c(ip);
            coa=coa+AREA(ip);
            sidediv=[sidediv ip];
        
    end

end
    
%end


%%%%%%%%conteo de reproduccion %%%%%%%%%%%%%%%
if ~isempty(xad)
    %xori=length(x);
   
    xs=x(xori+1:xori+length(xv2));
    ys=y(xori+1:xori+length(xv2));
    Eps=Ep(xori+1:xori+length(xv2));
    vxs=vx(xori+1:xori+length(xv2));
    vys=vy(xori+1:xori+length(xv2));
    cs=c(xori+1:xori+length(xv2));
    rs=r(xori+1:xori+length(xv2));
    AREAs=AREA(xori+1:xori+length(xv2));
	
    x=x(1:xori);
    y=y(1:xori);
    vx=vx(1:xori);
    vy=vy(1:xori);	
    Ep=Ep(1:xori);
    Ec=Ec(1:xori);
    c=c(1:xori);
    r=r(1:xori);
    AREA=AREA(1:xori);
    x=[x;xad];
    y=[y;yad];
    vx=[vx;vxad];
    vy=[vy;vyad];
    Ep=[Ep;Epad];
    Ec=[Ec;Ecad];
    c=[c;cad];
    r=[r;rad];
    AREA=[AREA;AREAad];
    %**********************************************************
    
    xori=length(x);
    cuentas=cuentas+1;
    xoric(cuentas)=xori-nc;
    nc=xori;
    x=[x;xs];
    y=[y;ys];
    vx=[vx;vxs];
    vy=[vy;vys];
    Ep=[Ep;Eps];
    c=[c;cs];
    r=[r;rs];
    AREA=[AREA;AREAs];
    x=[x;xv3];
    y=[y;yv3];
    vx=[vx;zeros((length(xv3)),1)];
    vy=[vy;zeros((length(xv3)),1)];
    Ep=[Ep;zeros((length(xv3)),1)];
    Ec=[Ec;zeros((length(xv3)),1)];
    c=[c;zeros((length(xv3)),1)];
    r=[r;zeros((length(xv3)),1)];
    
    invc=invc+1*(coa/(2*parax2))*length(xad);
    rep=length(xad)+rep;
    %coa=coa/length(xad);
end

