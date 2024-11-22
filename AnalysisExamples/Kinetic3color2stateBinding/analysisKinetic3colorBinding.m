clear all

channels=1:5;
numchans=length(channels);
donchan=3;
accchan=1;
accchan2=2;

numrep=1;

trjT = 1*1e7;  % in 100 ns
syncrate = 1; % 1: 100 ns, 2: 50 ns
Tsep = 5000; % 5 ms, time separation between trajectories of trjT
tick = 100/syncrate;    % 20 MHz, syncrate = 2;
focalvol = [0.5 1.5]; % confocal volume shape [xy z]*(1 um)

simid=1;
switch simid
    case 1
% two components, TAD-NCBD
nstate=2;
n0 = [100 200 180];     % ms-1, [Bound_3c, Unound_DA1, Bound_DA1]
D_diff = focalvol(1)^2*1.608./[1.8 1.2 1.2]/1000;     % um^2/ms, [B_3c, U_DA1, B_DA1]
eff1 = [0.08];  % acceptor 1 fraction, B_3c
eff2 = [0.7];   % acceptor 2 fraction, B_3c
leak12 = 0.2;  % A1 -> A2
eff = [0.38 0.8];   % 2-color FRET efficiency, U_DA1, B_DA1
dleak=0.05;     % Leak of donor (D) photons into A1 channel
eff1 = [eff1 eff*(1-leak12)+(1-eff)*dleak];     % combine acceptor 1 fraction of 3c and DA1 states
eff2 = [eff2 eff*leak12];   % combine acceptor 2 fraction of 3c and DA1 states, ignore dleak to A2 channel
peq0 = [0.3 0.7];    % equilibrium population, B (3c + DA1), U
p3c=0.45;    % A2 labeling efficiency, i,e., fraction of B_3c among the bound state
peq=[peq0(1)*p3c peq0(2) peq0(1)*(1-p3c)];
ratesum = 0.3;  % ms-1
rates=ratesum*([peq0(2) peq0(1)]);      % kinetic model, B_3c <-> U <-> B_DA1
Ntrj = 100;    % 3c
dt = 0.1;     % in us, simulated photon arrival times in 100 ns
rout = 4;
bkgcnt = [1.4 2.5 2.2];     % ms-1, [A1, A2, D]
filetext='TADNCBD_pB0.3k0.3p3c0.45_elong';
    case 2
end

%% burst selection
%$$$$$$$$$$$ CW analysis (i.e., Dex) $$$$$$$$$$$$$
phintvcrit = 300:-50:200; 	% maximum photon interval in microsec to be merged into the same burst
%phintv = phintvcrit(1);
phintv = phintvcrit(3);     % interphoton time threshold, 200 us
phthreshold = [30 40 50 60];
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

%% burstML analysis
% extract rates, FRET efficiencies, and equilibrium constant ---------------
drvintv=0.005;
drvpnt=5;
iscolor = true;
jmax = 30;
qmax = 4;
%t_th = phintvcrit(1)/1000;
t_th = phintvcrit(3)/1000;

ww=focalvol(2)/focalvol(1);
y0=2*ww*log(ww+sqrt(ww^2-1))/sqrt(ww^2-1);
vv=2*log(ww+sqrt(ww^2-1))/log((ww+sqrt(ww^2-1))*(sqrt(ww^2+y0)-sqrt(ww^2-1))/sqrt(1+y0));
elongfactor=2*y0/(2*vv-1-sqrt(4*vv+1));     % Conversion factor of extracted diffusion time in isotropic approximation into that of elongated focal shape

fretrngcrit = 1;     % 1: DA1A2 + DA1
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

        switch fretrngcrit
            case 1
                % 2-state, burstML1
                numC=3;
                initparams=[100 220 140 ... % (1:2n-1): n0(1) ~ n0(2n-1) for [B3c (DA1A2), U, B_DA1] 
                            1.5 1 1 ... % (2n:4n-2): tau_d(1) ~ tau_d(2n-1) for [B3c, U, B_DA1], 
                            0.5 ...   % (4n-1:5n-3): rate, k(1) ~ k(n-1)
                            0.25 ...   % (5n-2:6n-4): Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
                            ...   % f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
                            0.5 ...  % (6n-3): fraction parameters of DA1A2 (same format as p above), p3c, same for all bound species
                            0.1 0.65 ... % (6n-2:8n-5): A1, A2 fraction for DA1A2, e1(1) ~ e1(n-1), e2(1) ~ e2(n-1) 
                            0.45 0.75 ... % (8n-4:9n-5): A1+A2 fraction for DA1, [EU, E1(1) ~ E1(n-1)]
                            0.25 ...     % (9n-4): A1 leak into A2 (= A2/(A1+A2))
                            0.5 1.8 1.2]';   % (9n-3:9n-1): bkg cnt (A1, A2, D)
                LUbounds=[30 200; 50 300; 50 250;   % n0
                        0.1 3; 0.1 3; 0.1 3;     % tau_d
                        0.05 5;  % k
                        0.05 0.5; % peq
                        0.01 0.8;   % p3c
                        0.01 0.23; 0.6 0.85;   % DA1A2
                        0.3 0.55; 0.6 0.85;     % DA1
                        0.15 0.4;    % A1 leak
                        0.1 3; 0.1 5; 0.1 3];  % bkg
            case 2
        end
                
        switch fretrngcrit
            case 1
                % 2-state, burstML1
                resparams=mlhDiffNTR3cDA1CWBindbkg_MT(initparams,LUbounds,burstbint3r,cumindex,indexone,Diffparams{1},numC);
                nstate=(length(initparams)+4-numC)/9;
                if nstate > 2
                    pfactor=cumprod([1; 1-resparams(5*nstate-2:6*nstate-5)]);
                    resparams(5*nstate-2:6*nstate-4)=resparams(5*nstate-2:6*nstate-4).*pfactor; % Convert f to population p
                end
                resparam(:,nn,kk)=resparams;
                resparam(:,nn,kk)
                % Error is calculated for k (sum of rates between adjacent states) and relative population p.
                [errorparams logmlh(nn,kk) bic(nn,kk)]=mlhDiffNTR3cDA1CWBindbkgerrorC(resparam(:,nn,kk),burstbint3r,cumindex,indexone,Diffparams,numC,drvintv,drvpnt);
                errorparam(:,nn,kk)=errorparams;                                                        

                save(['trjkin3cCWDiff_' filetext '_Elong3cDA1.mat'],'resparam','errorparam','logmlh','bic','-mat');
            case 2
        end
    end
end
