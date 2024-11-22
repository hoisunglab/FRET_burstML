clear all

channels=[1 2];
numchans=length(channels);
donchan=2;
accchan=1;

numrep=1;

trjT = 1*1e7;  % in 100 ns
syncrate = 1; % 1: 100 ns, 2: 50 ns
Tsep = 5000; % 5 ms, time separation between trajectories of trjT
tick = 100/syncrate;    % 20 MHz, syncrate = 2;
focalvol = [0.5 1.5]; % confocal volume shape [xy z]*(1 um)

simid=1;
switch simid
    case 1
% two components
n0 = [180 240];     % ms-1, [F_DA1, U_DA1]
D_diff = focalvol(1)^2*1.608./[1.5 2]/1000;     % um^2/ms, [F_DA1, U_DA1]
eff = [0.8 0.4];   % DA1
peq = [0.4 0.6];    % F, U
ratesum = 1;  % ms-1
rates=ratesum*([peq(2) peq(1)]);
Ntrj = 50;    % DA1
ncolor = 1;   % DA1
dt = 0.1;     % in us, simulated photon arrival times in 100 ns
rout = 4;   % in um
bkgcnt = [2 1];     % ms-1, [A, D]
filetext='A3D_pF0.4_elong';
Erange=[0.1 0.6; 0.6 0.9];
    case 2

end

%% burst selection
phintvcrit = 300:-50:200; 	% maximum photon interval in microsec to be merged into the same burst
phintv = phintvcrit(1);
phthreshold = [30 40 50 60];
% find bursts ------------
load(['trjkin2cTR_' filetext '.mat'],'bint3r','-mat');
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
    save(['trjkin2cCWBurst_' filetext '.mat'] ,'burstbint3r','-mat');
end

% plot FRET histograms for check
load(['trjkin2cCWBurst_' filetext '.mat']);
numcol=max(numrep,4);
figure;
for kk=1:length(phthreshold)
    for nn=1:numrep
        [frhist x ev frburst] = calcfret(burstbint3r{kk}{nn}, donchan, accchan);
        subplot(length(phthreshold), numcol, numcol*(kk-1)+nn), bar(x, frhist, 1, 'edgecolor', 'none'), ylabel('counts'), ...
                xlim([-0.06 1.06]), ylim([0 max(frhist)*1.15]), xlabel('FRET efficiency'), ...
                text(0.25, max(frhist), sprintf('T=%d', phthreshold(kk)));
    end
end

% 2D duration vs. countr rate plot
plotcolor=[0 1 0; 1 0.5 0; 1 0 0];
figure;
for k=1:size(Erange,1)
%    for kk=1:length(phthreshold)
    for kk=1
        for nn=1:numrep
            [frhist x ev frburst] = calcfret(burstbint3r{kk}{nn}, donchan, accchan);
            brightout = calcburstL(burstbint3r{kk}{nn},[],syncrate(1));
            burstL=brightout(:,end);
            brightness=brightout(:,3);     % # photons per ms

            burstid=find(frburst(:,2) >= Erange(k,1) & frburst(:,2) < Erange(k,2));
            subplot(4,numcol,nn+(k-1)*numcol);
            plot(brightness(burstid),burstL(burstid),'.','color',plotcolor(k,:));
            xlim([0 300]);ylim([0.1 100]);
            set(gca,'yscale','log');
            if nn == 1, title([num2str(Erange(k,1)) ' < E < ',num2str(Erange(k,2))]);xlabel('Counts per ms'); end
            if k == size(Erange,1), ylabel('Burst duration (ms)'); end

            burstNdur(nn,1:2,k,kk)=[length(burstid) mean(burstL(burstid))];
        end
    end
end
%save(['trjkin2cCWBurstNdur_' filetext '.mat'] ,'burstNdur','-mat');

%% burstML analysis
% extract rates, FRET efficiencies, and equilibrium constant ---------------
drvintv=0.005;
drvpnt=5;
jmax = 50;
qmax = 4;
t_th = phintvcrit(1)/1000;

ww=focalvol(2)/focalvol(1);
y0=2*ww*log(ww+sqrt(ww^2-1))/sqrt(ww^2-1);
vv=2*log(ww+sqrt(ww^2-1))/log((ww+sqrt(ww^2-1))*(sqrt(ww^2+y0)-sqrt(ww^2-1))/sqrt(1+y0));
elongfactor=2*y0/(2*vv-1-sqrt(4*vv+1));     % Conversion factor of extracted diffusion time in isotropic approximation into that of elongated focal shape

clear resparam errorparam logmlh bic
for kk=1:length(phthreshold)
    kk
    load(['trjkin2cCWBurst_' filetext '.mat']);
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
    %    [frburstidsub,junk,indexone]=intersect(burstfrid,burstfridone);

        % n-state (2-state) burstML1
        numC=2;
        initparams=[150 250 ...     % (1:n): n0(1) ~ n0(n),
                    0.8 2.5 ...     % (n+1:2n): tau_d(1) ~ tau_d(n), 
                    0.8 ...         % (2n+1:3n-1): k(1) ~ k(n-1), relaxation rate, sum of forward and backward rates
                    0.5 ...         % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
                    ...             % (3n:4n-2): f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
                    0.85 0.45 ...    % (4n-1:5n-2): E(1) ~ E(n), FRET efficiency
                    2.5 0.8]';      % (5n-1:5n): bkg photon count rate, [acceptor donor]
                                    % for 3 color, (4n-1:5n-2): e1(1) ~ e1(n), (5n-1:6n-2): e2(1) ~ e2(n), (6n-1:6n+1): bkgcnt (A1, A2, D)
        LUbounds=[50 400; 100 500;
                0.1 3; 0.3 10;
                0.1 10;
                0.2 0.8;
                0.7 0.95; 0.3 0.6;
                0.5 5; 0.5 3];
                
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

        save(['trjkin2cCWDiff_' filetext '_burstML1.mat'],'resparam','errorparam','logmlh','bic','-mat');
    end
end
