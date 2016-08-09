%programa para monitoriar la reproduccion celular

varip=ip;

moni(varip,(length(moni(varip,:))+1))=counttotal;
monipy(varip,(length(monipy(varip,:))+1))=y(varip);

    [copt coptt]=max(monipy);
optt=[cuentas counttotal ip (nc+mnc) copt(end) ];%sum(y(bocn))/length(bocn)];
%mide la reprod,el tiempo de reprod,hija,madre,tama?o meristemo,boc
monimeris=[monimeris;optt];


%seguir hija, boc, meristemo, posicion

