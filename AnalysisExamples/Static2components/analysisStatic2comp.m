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
n0 = [150; 300];     % ms-1
eff = [0.1; 0.6];
D_diff = [0.001; 0.00025];     % um^2/ms
peq = [0.9; 0.1];
Ntrj = [36; 4];
peqeff = peq.*D_diff; peqeff=peqeff./sum(peqeff);
statecrit = [0; cumsum(peqeff(1:end-1)); 1.1];
dt = 0.1;     % in us, simulated photon arrival times in 100 ns
rout = 4;   % in um
bkgcnt = [3 1];     % ms-1
filetext='n150-300_eff0.1-0.6_D0.001-0.00025_p0.9_bkg3-1elong';
Erange=[-0.1 0.4; 0.4 0.99];
    case 2
end

%% burst selection
phintvcrit = 300:-50:100; 	% maximum photon interval in microsec to be merged into the same burst
phintv = phintvcrit(1);
phthreshold = [30 40 50 60];
% find bursts ------------
load(['trjkin_' filetext '.mat'],'bint3r','-mat');
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
    save(['trjkinBurst_' filetext '.mat'] ,'burstbint3r','-mat');
end

% plot FRET histograms for check
load(['trjkinBurst_' filetext '.mat']);
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
%save(['trjkinBurstNdur_' filetext '.mat'] ,'burstNdur','-mat');

%% burstML analysis
% extract rates, FRET efficiencies, and equilibrium constant ---------------
drvintv=0.005;
drvpnt=5;
iscolor = true;
isiptML=false;
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
    load(['trjkinBurst_' filetext '.mat']);
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

        if isiptML
            % n-state (2-state) iptML
            peq1=peq(1);    % iptML requires population of state 1:n-1
            initparamsIPT = [200 400]';     % count rates, n0(1) ~ n0(n)
            LUboundsIPT = [50 400; 150 1000];
            DiffparamsIPT = {Diffparams{1}, elongfactor*focalvol(1)^2./(D_diff*1000), peq1, Diffparams{2}};  % Diffusion time, population, and bkg need to be predetermined.
                            
            isbkg=false;    % bkg is not fitting parameter in iptML
            ispop=false;    % populateion is not determined in iptML
            iscoloript=false;   % total count rate is determined
            fixparamsIPT = [DiffparamsIPT{2}; DiffparamsIPT{3}];
            resparamsIPT = mlhipt_MT(initparamsIPT,LUboundsIPT,burstbint3r,cumindex,indexone,Diffparams{1},Diffparams{2},iscoloript,ispop,isbkg,fixparamsIPT);
            resparam(:,nn,kk)=resparamsIPT;
                        
            [errorparams logmlh(nn,kk) bic(nn,kk)]=mlhipterrorC(resparam(:,nn,kk),burstbint3r,cumindex,indexone,Diffparams,iscoloript,ispop,isbkg,fixparamsIPT,drvintv,drvpnt);
            errorparam(:,nn,kk)=errorparams;

            save(['trjkinDiff_' filetext '_iptML.mat'],'resparam','errorparam','logmlh','bic','-mat');
        else
            % n-state (2-state) burstML1
            initparams=[200 400 ...     % (1:n): n0(1) ~ n0(n),
                        0.8 3.5 ...     % (n+1:2n): tau_d(1) ~ tau_d(n), 
                        0.8 ...         % Instead of population p, a relative fraction f is used. f(i) = p(i)/(p(i)+p(i+1). Each f ranges from 0 to 1.
                        ...             % (2n+1:3n-1): f(1) ~ f(n-1), p(1) = f(1), p(2) = [1-f(1)]*f(2), p(3) = [1-f(1)]*[1-f(2)]*f(3), ... p(n) = [1-f(1)]*[1-f(2)]* ... *[1-f(n-1)]
                        0.06 0.5 ...    % (3n:4n-1): E(1) ~ E(n), FRET efficiency
                        4 1.5]';        % (4n:4n+1): bkg photon count rate, [acceptor donor]
            LUbounds=[50 400; 150 1000;
                    0.1 3; 0.3 10;
                    0.2 0.99;
                    0.01 0.2; 0.4 0.8;
                    0.5 5; 0.5 2];
                
            resparams=mlhDiffNbkg_MT(initparams,LUbounds,burstbint3r,cumindex,indexone,Diffparams{1},iscolor);
            
            if iscolor, nstate=(length(initparams)-1)/4;
            else nstate=(length(initparams)+1)/3;
            end
            if nstate > 1
                pfactor=cumprod([1; 1-resparams(2*nstate+1:3*nstate-2)]);
                resparams(2*nstate+1:3*nstate-1)=resparams(2*nstate+1:3*nstate-1).*pfactor; % Convert f to population p
            end
            resparam(:,nn,kk)=resparams;
            resparam(:,nn,kk)

            % Error is calculated for k (sum of rates between adjacent states) and relative population p.
            [errorparams logmlh(nn,kk) bic(nn,kk)]=mlhDiffNbkgerrorC(resparam(:,nn,kk),burstbint3r,cumindex,indexone,Diffparams,iscolor,drvintv,drvpnt);
            errorparam(:,nn,kk)=errorparams;
            
            save(['trjkinDiff_' filetext '_burstML1.mat'],'resparam','errorparam','logmlh','bic','-mat');
        end            
    end
end
