clear all

channels=1:5;
numchans=length(channels);
donchan=3;
accchan=1;
accchan2=2;
accchanAex=4;
accchan2Aex=5;

numrep=1;

trjT = 1*1e7;  % in 100 ns
syncrate = 1; % 1: 100 ns, 2: 50 ns
Tsep = 5000; % 5 ms, time separation between trajectories of trjT
tick = 100/syncrate;    % 20 MHz, syncrate = 2;
focalvol = [0.5 1.5]; % confocal volume shape [xy z]*(1 um)

simid=1;
switch simid
    case 1
% two components, simulation was performed with pulsed interleaved exictation (PIE) or alternating laser excitation (ALEX)
n0 = [180 240 280 300];     % (ms-1) [F_3c, U_3c, F_DA1, U_DA1]
D_diff = focalvol(1)^2*1.608./[1.6 2.1 1.5 2]/1000;     % um^2/ms, [F_3c, U_3c, F_DA1, U_DA1]
eff1 = [0.16 0.28];     % acceptor 1 (A1) fraction, 3-color component, [Folded (F), Unfolded (U)]
eff2 = [0.7 0.5];     % acceptor 2 (A2) fraction, 3-color component, [F, U]
leak12 = 0.25;  % Leak of A1 photons into A2 channel
eff = [0.9 0.78];   % DA1, 2-color FRET efficiency (D and A1), [F, U]
dleak=0.05;     % Leak of donor (D) photons into A1 channel
eff1 = [eff1 eff*(1-leak12)+(1-eff)*dleak];
eff2 = [eff2 eff*leak12];   % ignore dleak to A2 channel
eff12 = [0.84 0.64];    % A1A2, 2-color FRET efficeincy between A1 and A2
stoich = [0.5 0.5 0.5 0.5];   % Stoichiometry, Fraction of photons by acceptor 1 excitation (Aex), [F_3c, U_3c, F_DA1, U_DA1]
peq = [0.5 0.5];    % equilibrium population, [F, U]
ratesum = 1;  % relaxation rate, (ms-1)
rates=ratesum*([peq(2) peq(1)]);
Ntrj = [50 50];    % number of trajectories, [3c, DA1]
ncolor = length(Ntrj);
dt = 0.1;     % in us, simulated photon arrival times in 100 ns
rout = 4;
bkgcnt = [1 0.5 1.5 0.8 0.7];     % ms-1, [A1, A2, D, A1Aex, A2Aex]
filetext='A3D_pF0.5_elong';
    case 2
end

%% burst selection CW
%$$$$$$$$$$$ CW analysis (i.e., Dex) $$$$$$$$$$$$$
phintvcrit = 300:-50:200; 	% maximum photon interval in microsec to be merged into the same burst
phintv = phintvcrit(1);
phthreshold = [30 40 50 60];
% find bursts ------------
load(['trjkin3c_' filetext '.mat'],'bint3r','-mat');
clear burstbint3r
for kk=1:length(phthreshold)
    for nn=1:numrep
        t3r=zeros(size(bint3r{nn},1),3);
        t3r(:,1)=bint3r{nn}(:,1);
        t3r(:,3)=bint3r{nn}(:,2);

        t3r(find(t3r(:,end) == accchanAex | t3r(:,end) == accchan2Aex),:)=[];   % Delete photons by A1 excitation, i.e., mimicking CW donor excitation

        photoninterval = diff(t3r(:,1));
        burstsepid = [0; find(photoninterval*tick*1e-3 > phintv); size(t3r,1)];
        phthreshid = find(diff(burstsepid) >= phthreshold(kk));  % number of photons >= threshold

        burstids = zeros(size(t3r,1),1);
        burstindex = 1;
        for j = 1:length(phthreshid)
            burstids(burstsepid(phthreshid(j))+1:burstsepid(phthreshid(j)+1)) = burstindex;
            burstindex = burstindex + 1;
        end

        sigphotonids = find(burstids > 0);
        burstbint3r{kk}{nn}=zeros(length(sigphotonids),5);
        burstbint3r{kk}{nn}(:,1)=burstids(sigphotonids);
        burstbint3r{kk}{nn}(:,3:end) = t3r(sigphotonids,:);
    end
    save(['trjkin3cCWBurst_' filetext '.mat'] ,'burstbint3r','-mat');
end

% 2D plot FRET histograms for check
load(['trjkin3cCWBurst_' filetext '.mat']);

xlabeltext='E1';
ylabeltext='E2';
numrow=length(phthreshold);
islog = true;

clear frhist2
figure;
for kk=1:length(phthreshold)
    for nn=1:numrep
        [frhist xhist ev frburst] = calcfret3cCW(burstbint3r{kk}{nn}, donchan, accchan, accchan2);
        frhist2=histcounts2(frburst(:,3),frburst(:,2),xhist,xhist);
        frhist2one=frhist2;
        subplot(numrow,numrep,nn+(kk-1)*numrep);
        if islog
            tempval=log10(frhist2one);tempval(find(frhist2one == 0))=0;
            [xyz phandle]=contourf(xhist(1:end-1)+diff(xhist(1:2))/2,xhist(1:end-1)+diff(xhist(1:2))/2,tempval,max(max(tempval))*(0:0.1:1));ylim([0 1]);xlim([0 1]);axis square
        else
            [xyz phandle]=contourf(xhist(1:end-1)+diff(xhist(1:2))/2,xhist(1:end-1)+diff(xhist(1:2))/2,frhist2one,max(max(frhist2one))*(0:0.1:1));ylim([0 1]);xlim([0 1]);axis square
        end
        colormap([[1 1 1]; colormap2d1(8)]);
        set(phandle,'EdgeColor','none');
        text(0.7,0.9,sprintf('T=%d', phthreshold(kk)),'fontsize',12);
        colorbar;
        xlabel(xlabeltext);ylabel(ylabeltext);
    end
end

%% burstML analysis CW
% extract rates, FRET efficiencies, and equilibrium constant ---------------
drvintv=0.005;
drvpnt=5;
iscolor = true;
jmax = 30;
qmax = 4;
t_th = phintvcrit(1)/1000;

ww=focalvol(2)/focalvol(1);
y0=2*ww*log(ww+sqrt(ww^2-1))/sqrt(ww^2-1);
vv=2*log(ww+sqrt(ww^2-1))/log((ww+sqrt(ww^2-1))*(sqrt(ww^2+y0)-sqrt(ww^2-1))/sqrt(1+y0));
elongfactor=2*y0/(2*vv-1-sqrt(4*vv+1));     % Conversion factor of extracted diffusion time in isotropic approximation into that of elongated focal shape

trjcolortype = 1;     % 1: 3color + DA1, 2: 3-color
clear resparam errorparam logmlh bic
for kk=1:length(phthreshold)
    kk
    load(['trjkin3cCWBurst_' filetext '.mat']);
    N_th = phthreshold(kk);
    Diffparams={[jmax; qmax; t_th; N_th],bkgcnt};

    burstbint3rtmp=burstbint3r{kk};
    for nn=1:numrep
        nn
        burstbint3r=[];burstphotons=[];fileidtrack=[];
        burstbint3r=burstbint3rtmp{nn};

        burstfridone=unique(burstbint3r(:,1));
        histone=histc(burstbint3r(:,1),burstfridone);
        cumindex=[0; cumsum(histone)];
        indexone=1:length(burstfridone);

        switch trjcolortype
            case 1
                % 2-state, burstML1
                numC=3;
                initparams=[100 120 150 150 ... % (1:2n): n0(1) ~ n0(n) for 3c, DA1 
                            1.5 2.5 1.5 2.5 ... % (2n+1:4n): tau_d(1) ~ tau_d(n) for 3c, DA1, 
                            1.2 ...   % (4n+1:5n-1): rate, k(1) ~ k(n-1)
                            0.45 ...   % (5n:6n-2): Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
                            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
                            0.4 ...  % (6n-1): fraction parameters of DA1A2 (same format as p above), fsp(1)
                            0.16 0.26 0.7 0.5 ... % (6n:8n-1): A1, A2 fraction for 3c, e1(1) ~ e1(n), e2(1) ~ e2(n) 
                            0.88 0.75 ... % (8n:9n-1): A1+A2 fraction for DA1, E1(1) ~ E1(n)
                            0.2 ...     % (9n): A1 leak into A2 (= A2/(A1+A2))
                            0.5 0.5 0.5]';   % (9n+1:9n+3): bkg cnt (A1, A2, D)
                LUbounds=[30 1000; 30 1000; 50 1000; 50 1000;   % n0
                        0.1 10; 0.1 10; 0.1 10; 0.1 10;     % tau_d
                        0.1 10;  % k
                        0.01 0.95; % peq
                        0.01 0.8;   % species fsp
                        0.05 0.23; 0.15 0.35; 0.5 0.75; 0.45 0.63;   % DA1A2 (3c)
                        0.8 0.95; 0.7 0.85;     % DA1
                        0.15 0.3;    % A1 leak
                        0.1 5; 0.1 5; 0.1 5];  % bkg
            case 2
                % 2-state, burstML1-Kin
                numC=3;
                initparams=[100 120 ... % (1:n): n0(1) ~ n0(n) for 3c 
                            1.5 2.5 ... % (n+1:2n): tau_d(1) ~ tau_d(n) for 3c 
                            0.8 ...   % (2n+1:3n-1): rate, k(1) ~ k(n-1)
                            0.45 ...   % (3n:4n-2): Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
                            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
                            0.16 0.26 0.7 0.5 ... % (4n-1:6n-2): A1, A2 fraction for DA1A2, e1(1) ~ e1(n), e2(1) ~ e2(n) 
                            0.5 0.5 0.5]';   % (6n-1:6n+1): bkg cnt (A1, A2, D)
                LUbounds=[30 1000; 30 1000;   % n0
                        0.1 10; 0.1 10;     % tau_d
                        0.1 10;  % k
                        0.01 0.95; % peq
                        0.05 0.23; 0.15 0.35; 0.5 0.75; 0.45 0.63;   % DA1A2 (3c)
                        0.1 5; 0.1 5; 0.1 5];  % bkg
        end
                
        switch trjcolortype
            case 1
                % 2-state, burstML1
                numSp=2;
                resparams=mlhDiffNTR3cDA1CWbkg_MT(initparams,LUbounds,burstbint3r,cumindex,indexone,Diffparams{1},numC);
                nstate=(length(initparams)+2-numC-numSp)/(5+2*numSp);
                if nstate > 2
                    pfactor=cumprod([1; 1-resparams(5*nstate+1:6*nstate-2)]);
                    resparams(5*nstate+1:6*nstate-1)=resparams(5*nstate+1:6*nstate-1).*pfactor; % Convert f to population p
                end
                if numSp > 2
                    pfactor=cumprod([1; 1-resparams(6*nstate:6*nstate+numSp-3)]);
                    resparams(6*nstate:6*nstate+numSp-2)=resparams(6*nstate:6*nstate+numSp-2).*pfactor; % Convert f to population p
                end                                
                resparam(:,nn,kk)=resparams;
                resparam(:,nn,kk)
                % Error is calculated for k (sum of rates between adjacent states) and relative population p.
                [errorparams logmlh(nn,kk) bic(nn,kk)]=mlhDiffNTR3cDA1CWbkgerrorC(resparam(:,nn,kk),burstbint3r,cumindex,indexone,Diffparams,numC,drvintv,drvpnt);
                errorparam(:,nn,kk)=errorparams;                                                        

                save(['trjkin3cCWDiff_' filetext '_Elong3cDA1.mat'],'resparam','errorparam','logmlh','bic','-mat');
            case 2
                % 2-state, burstML1
                resparams=mlhDiffNTRbkg_MT(initparams,LUbounds,burstbint3r,cumindex,indexone,Diffparams{1},numC);
                nstate=(length(initparams)+2-numC)/(3+numC);
                if nstate > 1
                    pfactor=cumprod([1; 1-resparams(3*nstate:4*nstate-3)]);
                    resparams(3*nstate:4*nstate-2)=resparams(3*nstate:4*nstate-2).*pfactor; % Convert f to population p
                end
                resparam(:,nn,kk)=resparams;
                resparam(:,nn,kk)
                % Error is calculated for k (sum of rates between adjacent states) and relative population p.
                [errorparams logmlh(nn,kk) bic(nn,kk)]=mlhDiffNTRbkgerrorC(resparam(:,nn,kk),burstbint3r,cumindex,indexone,Diffparams,numC,drvintv,drvpnt);
                errorparam(:,nn,kk)=errorparams;                            

                save(['trjkin3cCWDiff_' filetext '_Elong3c.mat'],'resparam','errorparam','logmlh','bic','-mat');
        end
    end
end


%% burst selection ALEX
%$$$$$$$$$$$ ALEX analysis (i.e., Dex) $$$$$$$$$$$$$
phintvcrit = 200:-50:100; 	% maximum photon interval in microsec to be merged into the same burst
phintv = phintvcrit(2); % 150 us
phthreshold = 60:20:120;
% find bursts ------------
load(['trjkin3c_' filetext '.mat'],'bint3r','-mat');
clear burstbint3r
for kk=1:length(phthreshold)
    for nn=1:numrep
        t3r=zeros(size(bint3r{nn},1),3);
        t3r(:,1)=bint3r{nn}(:,1);
        t3r(:,3)=bint3r{nn}(:,2);

        photoninterval = diff(t3r(:,1));
        burstsepid = [0; find(photoninterval*tick*1e-3 > phintv); size(t3r,1)];
        phthreshid = find(diff(burstsepid) >= phthreshold(kk));  % number of photons >= threshold

        burstids = zeros(size(t3r,1),1);
        burstindex = 1;
        for j = 1:length(phthreshid)
            burstids(burstsepid(phthreshid(j))+1:burstsepid(phthreshid(j)+1)) = burstindex;
            burstindex = burstindex + 1;
        end

        sigphotonids = find(burstids > 0);
        burstbint3r{kk}{nn}=zeros(length(sigphotonids),5);
        burstbint3r{kk}{nn}(:,1)=burstids(sigphotonids);
        burstbint3r{kk}{nn}(:,3:end) = t3r(sigphotonids,:);
    end
    save(['trjkin3cALEXBurst_' filetext '.mat'] ,'burstbint3r','-mat');
end

% 2D plot FRET histograms for check
load(['trjkin3cALEXBurst_' filetext '.mat']);

xlabeltext={'E1','E1','E2','E1','E2'};
ylabeltext={'E2','E12','E12','S','S'};
numcol=length(xlabeltext);
numrow=length(phthreshold);
islog = true;
isFRET3c=false; % calcfret3c Returns FRET efficiency (true) or acceptor fraction (false) 

clear frhist2
for nn=1:numrep
figure;
    for kk=1:length(phthreshold)
        [frhist xhist ev frburst] = calcfret3cSim(burstbint3r{kk}{nn}, donchan, accchan, accchan2, isFRET3c);
        frhist2{1}=histcounts2(frburst(:,3),frburst(:,2),xhist,xhist);
        frhist2{2}=histcounts2(frburst(:,4),frburst(:,2),xhist,xhist);
        frhist2{3}=histcounts2(frburst(:,4),frburst(:,3),xhist,xhist);
        frhist2{4}=histcounts2(frburst(:,10),frburst(:,2),xhist,xhist);
        frhist2{5}=histcounts2(frburst(:,10),frburst(:,3),xhist,xhist);

        for ii=1:numcol
            frhist2one=frhist2{ii};
            subplot(numrow,numcol,ii+(kk-1)*numcol);
            if islog
                tempval=log10(frhist2one);tempval(find(frhist2one == 0))=0;
                [xyz phandle]=contourf(xhist(1:end-1)+diff(xhist(1:2))/2,xhist(1:end-1)+diff(xhist(1:2))/2,tempval,max(max(tempval))*(0:0.1:1));ylim([0 1]);xlim([0 1]);axis square
            else
                [xyz phandle]=contourf(xhist(1:end-1)+diff(xhist(1:2))/2,xhist(1:end-1)+diff(xhist(1:2))/2,frhist2one,max(max(frhist2one))*(0:0.1:1));ylim([0 1]);xlim([0 1]);axis square
            end
            colormap([[1 1 1]; colormap2d1(8)]);
            set(phandle,'EdgeColor','none');
            text(0.7,0.9,sprintf('T=%d', phthreshold(kk)),'fontsize',12);
            colorbar;
            xlabel(xlabeltext{ii});ylabel(ylabeltext{ii});
        end
    end
end

%% burstML analysis ALEX
% extract rates, FRET efficiencies, and equilibrium constant ---------------
drvintv=0.005;
drvpnt=5;
iscolor = true;
jmax = 30;
qmax = 4;
t_th = phintvcrit(2)/1000;  % 150 us

ww=focalvol(2)/focalvol(1);
y0=2*ww*log(ww+sqrt(ww^2-1))/sqrt(ww^2-1);
vv=2*log(ww+sqrt(ww^2-1))/log((ww+sqrt(ww^2-1))*(sqrt(ww^2+y0)-sqrt(ww^2-1))/sqrt(1+y0));
elongfactor=2*y0/(2*vv-1-sqrt(4*vv+1));     % Conversion factor of extracted diffusion time in isotropic approximation into that of elongated focal shape

trjcolortype = 1;    % 1: DA1A2 (3c) + DA1, 2: DA1A2 (3-color)
clear resparam errorparam logmlh bic
for kk=1:length(phthreshold)
    kk
    load(['trjkin3cALEXBurst_' filetext '.mat']);
    N_th = phthreshold(kk);
    Diffparams={[jmax; qmax; t_th; N_th],bkgcnt};

    burstbint3rtmp=burstbint3r{kk};
    for nn=1:numrep
        nn
        burstbint3r=[];burstphotons=[];fileidtrack=[];
        burstbint3r=burstbint3rtmp{nn};

        burstfridone=unique(burstbint3r(:,1));
        histone=histc(burstbint3r(:,1),burstfridone);
        cumindex=[0; cumsum(histone)];
        indexone=1:length(burstfridone);

        switch trjcolortype
            case 1
                % 2-state
                numC=3;
                initparams=[[200 200 300 300] ... % (1:2n): n0(1) ~ n0(n) for 3c, DA1 
                            1.5 2.5 1.5 2.5 ... % (2n+1:4n): tau_d(1) ~ tau_d(n) for 3c, DA1, 
                            0.8 ...   % (4n+1:5n-1): rate, k(1) ~ k(n-1)
                            0.45 ...   % (5n:6n-2): Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
                            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
                            0.4 ...  % (6n-1): fraction parameters of DA1A2 (same format as p above), fsp(1)
                            0.16 0.26 0.7 0.5 ... % (6n:8n-1): A1, A2 fraction for DA1A2, e1(1) ~ e1(n), e2(1) ~ e2(n) by Dex
                            0.88 0.75 ... % (8n:9n-1): A1+A2 fraction for DA1, E1(1) ~ E1(n) by Dex
                            0.8 0.6 ... % (9n:10n-1): E12(1) ~ E12(n) for DA1A2 by Aex
                            0.2 ...     % (10n): A1 leak into A2 (= A2/(A1+A2))
                            0.4 0.5 0.5 0.5 ... % (10n+1:12n): Stoichiometry
                            0.5 0.5 0.5 0.5 0.5]';   % (12n+1:12n+5): bkg cnt (A1, A2, D) by Dex + (A1, A2) by Aex
                LUbounds=[10 1000; 10 1000; 10 1000; 10 1000;   % n0
                        0.1 10; 0.1 10; 0.1 10; 0.1 10;     % tau_d
                        0.1 10;  % k
                        0.01 0.95; % peq
                        0.01 0.8;   % species fsp
                        0.05 0.23; 0.15 0.35; 0.5 0.75; 0.45 0.63;   % DA1A2
                        0.8 0.95; 0.7 0.85;     % DA1
                        0.65 0.95; 0.55 0.8;     % DA1A2, E12
                        0.15 0.35;    % A1 leak
                        0.3 0.7; 0.3 0.7; 0.3 0.7; 0.3 0.7;  % stoichiometry
                        0.1 5; 0.1 5; 0.1 5; 0.1 5; 0.1 5];  % bkg

            case 2
                % 2-state
                numC=3;
                initparams=[200 200 ... % (1:n): n0(1) ~ n0(n) for 3c 
                            1.5 2.5 ... % (n+1:2n): tau_d(1) ~ tau_d(n) for 3c 
                            1 ...   % (2n+1:3n-1): rate, k(1) ~ k(n-1)
                            0.5 ...   % (3n:4n-2): Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
                            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
                            0.16 0.26 0.7 0.5 ... % (4n-1:6n-2): A1, A2 fraction for DA1A2, e1(1) ~ e1(n), e2(1) ~ e2(n) by Dex
                            0.8 0.6 ... % (6n-1:7n-2): E12(1) ~ E12(n) for DA1A2 by Aex
                            0.45 0.45 ... % (7n-1:8n-2): Stoichiometry
                            0.5 0.5 0.5 0.5 0.5]';   % (8n-1:8n+3): bkg cnt (A1, A2, D) by Dex + (A1, A2) by Aex
                LUbounds=[50 1000; 50 1000;   % n0
                        0.1 10; 0.1 10;     % tau_d
                        0.1 10;  % k
                        0.01 0.95; % peq
                        0.05 0.23; 0.15 0.35; 0.5 0.75; 0.45 0.63;   % DA1A2
                        0.65 0.95; 0.55 0.8;     % DA1A2, E12
                        0.3 0.7; 0.3 0.7;  % stoichiometry
                        0.1 5; 0.1 5; 0.1 5; 0.1 5; 0.1 5];  % bkg
            case 3
            case 4
        end
                
        switch trjcolortype
            case 1
                numSp=2;
                resparams=mlhDiffNTR3cDA1ALEXbkg_MT(initparams,LUbounds,burstbint3r,cumindex,indexone,Diffparams{1},numC);
                nstate=(length(initparams)+3-2*numC-numSp)/(6+3*numSp);
                if nstate > 2
                    pfactor=cumprod([1; 1-resparams(5*nstate+1:6*nstate-2)]);
                    resparams(5*nstate+1:6*nstate-1)=resparams(5*nstate+1:6*nstate-1).*pfactor; % Convert f to population p
                end
                if numSp > 2
                    pfactor=cumprod([1; 1-resparams(6*nstate:6*nstate+numSp-3)]);
                    resparams(6*nstate:6*nstate+numSp-2)=resparams(6*nstate:6*nstate+numSp-2).*pfactor; % Convert f to population p
                end                                
                resparam(:,nn,kk)=resparams;
                resparam(:,nn,kk)
                % Error is calculated for k (sum of rates between adjacent states) and relative population p.
                [errorparams logmlh(nn,kk) bic(nn,kk)]=mlhDiffNTR3cDA1ALEXbkgerrorC(resparam(:,nn,kk),burstbint3r,cumindex,indexone,Diffparams,numC,drvintv,drvpnt);
                errorparam(:,nn,kk)=errorparams;                                                        

                save(['trjkin3cALEXDiff_' filetext '_Elong3cDA1.mat'],'resparam','errorparam','logmlh','bic','-mat');
            case 2
                % 2-state, burstML1-Kin
                resparams=mlhDiffNTR3cALEXbkg_MT(initparams,LUbounds,burstbint3r,cumindex,indexone,Diffparams{1},numC);
                nstate=(length(initparams)+3-2*numC)/(2+2*numC);
                if nstate > 2
                    pfactor=cumprod([1; 1-resparams(3*nstate:4*nstate-3)]);
                    resparams(3*nstate:4*nstate-2)=resparams(3*nstate:4*nstate-2).*pfactor; % Convert f to population p
                end
                resparam(:,nn,kk)=resparams;
                resparam(:,nn,kk)
                % Error is calculated for k (sum of rates between adjacent states) and relative population p.
                [errorparams logmlh(nn,kk) bic(nn,kk)]=mlhDiffNTR3cALEXbkgerrorC(resparam(:,nn,kk),burstbint3r,cumindex,indexone,Diffparams,numC,drvintv,drvpnt);
                errorparam(:,nn,kk)=errorparams;                            

                save(['trjkin3cALEXDiff_' filetext '_Elong3c.mat'],'resparam','errorparam','logmlh','bic','-mat');
            case 3
            case 4
        end
    end
end
