function [errorparams logmlh1 bic]=mlhDiffNbkgerrorC(initparams,frburstdata,cumindex,indexone,Diffparams,iscolor,drvintv,drvpnt)
% Beacuse evaluating mlhratesub takes long time, substitute this
% sub-routine to a C-translated version.

if size(initparams,1) < size(initparams,2) % It should be a column vector
    initparams=initparams';
end

drvintvs=1+drvintv*(-drvpnt:drvpnt);
if iscolor
    nstate = (length(initparams)-1)/4;
else
    nstate = length(initparams)/3;
end

initparammat=[];
for kk=1:length(initparams)
    for nn=1:length(drvintvs)
        initparamone=initparams;
        initparamone(kk)=initparams(kk)*drvintvs(nn);
        if nstate > 1
            peqtemp=initparamone(2*nstate+1:3*nstate-1);
            frn=zeros(size(peqtemp));
            frn(1)=peqtemp(1);
            for jj=1:nstate-2;
                fdiv=prod(1-frn(1:jj));
                frn(jj+1)=peqtemp(jj+1)/fdiv;
            end
            initparamone(2*nstate+1:3*nstate-1)=frn;
        end
        initparammat=[initparammat initparamone];
    end
end

for kk=1:length(initparams)
    for mm=kk+1:length(initparams)
        for nn=1:length(drvintvs)
            initparamone=initparams;
            initparamone(kk)=initparams(kk)*drvintvs(nn);initparamone(mm)=initparams(mm)*drvintvs(nn);
            if nstate > 1
                peqtemp=initparamone(2*nstate+1:3*nstate-1);
                frn=zeros(size(peqtemp));
                frn(1)=peqtemp(1);
                for jj=1:nstate-2;
                    fdiv=prod(1-frn(1:jj));
                    frn(jj+1)=peqtemp(jj+1)/fdiv;
                end
                initparamone(2*nstate+1:3*nstate-1)=frn;
            end
            initparammat=[initparammat initparamone];
        end
    end
end

LUbounds=[initparammat(:,(length(drvintvs)+1)/2)*0.8 initparammat(:,(length(drvintvs)+1)/2)*1.2]; % dummy LUbounds
logmlh=mlhDiffNbkg_MT(initparammat,LUbounds,frburstdata,cumindex,indexone,Diffparams{1},iscolor);

Hessianmat=zeros(length(initparams),length(initparams));
logmlhdiag=zeros(1,length(drvintvs));

mlhIndex=1;
for kk=1:length(initparams)
    logmlhdiag=logmlhdiag*0;
    for nn=1:length(drvintvs)
        logmlhdiag(nn)=logmlh(mlhIndex);
        mlhIndex=mlhIndex+1;
    end
    diagcoeff=polyfit(drvintvs,logmlhdiag,2);
    Hessianmat(kk,kk)=diagcoeff(1)/initparams(kk)^2;
end

logmlhoffdiag=zeros(1,length(drvintvs));
for kk=1:length(initparams)
    for mm=kk+1:length(initparams)
        logmlhoffdiag=logmlhoffdiag*0;
        for nn=1:length(drvintvs)
            logmlhoffdiag(nn)=logmlh(mlhIndex);
            mlhIndex=mlhIndex+1;
        end
        offdiagcoeff=polyfit(drvintvs,logmlhoffdiag,2);
        Hessianmat(kk,mm)=(offdiagcoeff(1)-Hessianmat(kk,kk)*initparams(kk)^2-Hessianmat(mm,mm)*initparams(mm)^2)/2/initparams(kk)/initparams(mm);
        Hessianmat(mm,kk)=Hessianmat(kk,mm);
    end
end

covmat=inv(2*Hessianmat);
errorparams=sqrt(diag(covmat));

logmlh1=-logmlhdiag((length(drvintvs)+1)/2);
numphotons=sum(cumindex(indexone+1)-cumindex(indexone));
bic=-2*logmlh1 + length(initparams)*log(numphotons);

end