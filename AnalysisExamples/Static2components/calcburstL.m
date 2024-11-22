function varargout = calcburstL(inputdata, varargin)
%
% bLout = calcburstL(inputdata)
% [bLout brighthist xaxis] = calcburstL(inputdata)
% [bLout brighthist] = calcburstL(inputdata, xaxis, syncrate)
%
% Calculates the brightness of a burst (or sub-burst) and burst duration.
% Brightness is defined as photons/1 ms. inputdata is a 5
% column matrix with burst-bin-tagged T3R data, or a filename containing a
% list of burstbint3r data files.  
%
% Defalt unittime (time resolution) is 100 ns (definition is different from
% calcbrightness) and can be assigned by syncrate (in 10 MHz. e.g., 2 for
% 20 MHz).
% bLout is a 5 column matrix with 
% bLout(i,:)=[burst index, photon counts, brightness, duration by bins, duration by time (ms)]
%
% optional output is a histogram of burst length in bin time unit, with an optional user
% supplied xaxis. DEFAULT = 0:maximum(burstlength)+1.


% decide if burstbint3r is a matrix or a filename.  Calculate brightness if
% it is a matrix, or read in each file and calculate brightness for each
% file if inputdata is filename.

% initiailize unittime value.
unittime = 100;     % 100 ns
if nargin >=3
    unittime = unittime / varargin{2};
end

if isnumeric(inputdata)
    
    bLout = calcburstLsub(inputdata, unittime);
    
elseif ischar(inputdata)
    
    datafiles = textread(inputdata, '%s', 'delimiter', '\n');
    numfiles = length(datafiles);
    
    bLout = [];
    
    % loop over each file, and add results from each file to frethist
    % and events.
    for fileindex=1:numfiles
        load(datafiles{fileindex}, '-mat');
        if ~exist('burstbint3r')
            error('Data in file %s must be in burstbint3r format.', datafiles{fileindex});
        end
        onebLout = calcburstLsub(burstbint3r, unittime);
        bLout = [bLout; onebLout];       
        clear burstbint3r;
    end
else
    error('Problem with inputdata.  Must be a 5-column matrix OR filename.');
end 


maxbL = max(bLout(:,4));

xaxis = 0:maxbL+1;
if nargin >= 2 && ~isempty(varargin{1})
    xaxis = varargin{1};
end

bLhist = [];
if nargout>=2
    bLhist = histc(bLout(:,4), xaxis);
end

% pack output arguments
k{1} = bLout;
k{2} = bLhist;
k{3} = xaxis;
varargout = k(1:nargout);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end calcburstL  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% subfunction to calculate burst length of a burst-bin-tagged T3R matrix.
function bLout = calcburstLsub(burstbint3r, unittime)

% exit from function if burstbint3r is empty
if isempty(burstbint3r)
    bLout = [];
    return
end

% remove marker records
markerindex = find(burstbint3r(:,end) < 0);
burstbint3r(markerindex,:) = [];

% find burst indices and counts of all photon records.
burstnumber = unique(burstbint3r(:,1));
burstcount = histc(burstbint3r(:,1), burstnumber);
burstcum = [0; cumsum(burstcount)];
numbursts = length(burstnumber);

bLout = zeros(numbursts, 5);
bLout(:,1:2) = [burstnumber burstcount];

% calculate burst length for each burst and assign values to bLout.
for i=1:numbursts
    starttime = burstbint3r(burstcum(i)+1,3);
    endtime   = burstbint3r(burstcum(i+1),3);
    duration = endtime - starttime;
    numbins = length(unique(burstbint3r(burstcum(i)+1:burstcum(i+1),2)));
    bLout(i,3:5) = [burstcount(i)/(duration*unittime*1e-6) numbins duration*unittime*1e-6];     % conversion of duration in ms
end