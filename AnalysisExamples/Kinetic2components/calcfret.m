function varargout = calcfret(inputdata, donorChans, acceptChans, varargin)
%
% calcfret(inputdata, donorChans, acceptChans)
% calcfret(inputdata, donorChans, acceptChans, gamma)
% calcfret(inputdata, donorChans, acceptChans, gamma, xAxis)
% calcfret(inputdata, donorChans, acceptChans, gamma, xAxis, donorleak)
% frethist = calcfret(....)
% [frethist xAxis] = calcfret(....)
% [frethist xAxis events] = calcfret(....)
% [frethist xAxis events fretBurst] = calcfret(....)
%
% Calculate FRET histogram for data using intensity ratios.
% inputdata may either be a 5-column matrix of burst-bin-tagged T3R data,
% or a filename containing .mat files.
% The donor and acceptor channel indices must be specified.  
% If gamma is not specified, gamma=1. 
% Optional input is the x-axis used in calculating the histogram.   
% donorleak is optional input specifying probability of donor photon leakage
% into acceptor channel.  (DEFAULT = 0).
% 
% If no output arguments are specified, then plot the FRET data. 
% Optional output are the FRET histogram data, x-axis, number of events,
% and the FRET value for each individual burst.  
% fretBurst is a n-by-4 matrix with: 
%       [Burst index, FRET values, Total counts, donor counts] 


% initialize gamma
gamma = 1;
if nargin >=4 && ~isempty(varargin{1})
    gamma = varargin{1};
end

% initialize xAxis
xAxis = [0:0.05:1]';
if nargin >=5 && ~isempty(varargin{2})
    xAxis = varargin{2};
end

% initialize donorleak
donorleak = 0;
if nargin >=6 && ~isempty(varargin{3})
    donorleak = varargin{3};
end

if isnumeric(inputdata)  % if inputdata is matrix of data, analyze data.
    [rows, columns] = size(inputdata);
    if columns ~= 5
        error('input must be 5-column matrix OR filename.');
    end
    
    [frethist, events, fretBurst] = calcfretsub(inputdata, donorChans, ...
        acceptChans, gamma, xAxis, donorleak);
    
elseif ischar(inputdata) % if inputdata is filename, load and analyze each file.
    
    datafiles = textread(inputdata, '%s', 'delimiter', '\n');
    numfiles = length(datafiles);
    frethist = zeros(length(xAxis), 1);
    events = 0;
    fretBurst = [];
    
    % loop over each file, and add results from each file to Totalfrethis
    % and Totalevents.
    for fileindex=1:numfiles
        
        data = loadt3rfile(datafiles{fileindex});

        [OneFrethist, OneEvents, OnefretBurst] = calcfretsub(data, ...
            donorChans, acceptChans, gamma, xAxis, donorleak);
        
        frethist = frethist + OneFrethist;
        events = events + OneEvents;
        fretBurst = [fretBurst; OnefretBurst];
        clear data;
    end 
else
    error('Problem with inputdata.  Must be a 5-column matrix OR filename.');
end 

% plot frethist and pack output arguments as necessary.
if nargout==0
    figure;
    bar(xAxis, frethist)
else
    k{1} = frethist;
    k{2} = xAxis;
    k{3} = events;
    k{4} = fretBurst;
    
    varargout=k(1:nargout);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end calcfret  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Calculate FRET histogram for burstbint3r data using intensity ratios.  
function [frethist, events, fretBurst] = calcfretsub(burstbint3r, donorChans, ...
    acceptChans, gamma, xAxis, donorleak)

% handling of empty burstbint3r matrices
if isempty(burstbint3r)
    frethist = zeros(length(xAxis), 1);
    events = 0;
    fretBurst = zeros(0,4);  % 4-column empty matrix.
    
    return
end

% remove marker records
markerindex = find(burstbint3r(:,5)<0);
markerrecords = burstbint3r(markerindex,:);
burstbint3r(markerindex,:) = [];

% separate donor and acceptor events in burstbint3r.  
donorData = burstbint3r(getphotonindices(burstbint3r, donorChans),:);
accptData = burstbint3r(getphotonindices(burstbint3r, acceptChans),:);
minbin = min(burstbint3r(:,1));
maxbin = max(burstbint3r(:,1));

% Calculate number of donor and acceptor counts for each burst, and
% calculate the FRET value for each burst.  
burstnumber = unique(burstbint3r(:,1));
donorBurst = histc(donorData(:,1), burstnumber);
accptBurst = histc(accptData(:,1), burstnumber);

% correct donor and acceptor bursts for donor leakage, and then calculate
% FRET efficiency with these corrected values and gamma.
donorBurstcorr =  donorBurst*(1+donorleak);
accptBurstcorr =  accptBurst - donorBurst*donorleak;

effBurst = accptBurstcorr./(donorBurstcorr*gamma + accptBurstcorr);

fretBurst = [burstnumber, effBurst, donorBurst+accptBurst, donorBurst];

% build FRET histogram, 
frethist = histc(fretBurst(:,2), xAxis);
if length(frethist(1,:)) > 1  % transpose row vector to column vector if necessary.
    frethist = frethist';
end
events = sum(frethist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end calcfretsub %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%% begin getphotonindices %%%%%%%%%%%%%%%%%%%%%%%%%
function indices = getphotonindices(data,chans)
indices = [];

% retrieve row indices for all photons coming on chans
for i=1:length(chans)
    indices = [indices; find(data(:,end) == chans(i))];
end
indices = sort(indices);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end getphotonindices %%%%%%%%%%%%%%%%%%%%%%%%%