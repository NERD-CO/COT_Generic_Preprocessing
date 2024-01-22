% Do fractional interval firing rates, with arbitrary bin sizes and
% allowing multiple spikes inside any bin size. Going to assume all bins
% are the same size though. Sorry it's slow. Rex Tien Jan 2020

% % for testing
% for t = 1:150
%     for p = 1:50
% spiketimes = (C.spikes{t,p}');
% binedges = FRbedges1ms;
% binedges = 
% binedges = C.t.reach(t):dt(t):C.t.threshold(t);

function firates = Spikes2FIR_Arbitrary( spiketimes, binedges)

% eliminate unnecessary spikes

addfrontspike = [];
addbackspike = [];
firstspikei = find(spiketimes < binedges(1),1,'last');
if isempty(firstspikei)
    firstspikei = 1;
    addfrontspike = -Inf;
end
lastspikei = find(spiketimes > binedges(end),1,'first');
if isempty(lastspikei)
    lastspikei = length(spiketimes);
    addbackspike = Inf;
end

if (size(spiketimes(firstspikei:lastspikei),1) > size(spiketimes(firstspikei:lastspikei),2))
    stimes = [addfrontspike; spiketimes(firstspikei:lastspikei); addbackspike];
else
    stimes = [addfrontspike spiketimes(firstspikei:lastspikei) addbackspike];
end

binsize = binedges(2)-binedges(1);

nbin = length(binedges)-1;
intervals = diff(stimes);
nint = length(intervals);
lastinti = 1;

firates = nan(1,nbin);
if isempty(stimes)
    firates = zeros(1,nbin);
else
    for i = 1:nbin
        
        % For each bin, identify which intervals belong to it. First ID which
        % spikes belong to it
        spikesin = (stimes >= binedges(i)) & (stimes < binedges(i+1));
        intsin = spikesin(2:end);
        
        % carry over the last interval from the previous bin (this keeps the
        % first and last interval propagating if there are no more spikes
        intsin(lastinti) = 1;
        
        % also, the next interval is in play if you counted any spikes! But not if you're at the end
        if lastinti < nint && sum(spikesin)>0
            lastinti = find(intsin,1,'last');
            intsin(lastinti+1) = 1;
            lastinti = lastinti + 1;
        end
        
        
        % rates are inverse of the relevant intervals
        rates = 1./intervals(intsin);
        
        % calculate
        if sum(spikesin) == 0
            fracts = 1;
        else
            fracts = nan(sum(intsin),1);
            stimesin = stimes(spikesin);
            fracts(1) = (stimesin(1)-binedges(i))/binsize;
            for j = 1:sum(spikesin)-1
                fracts(j+1) = (stimesin(j+1)-stimesin(j))/binsize;
            end
            fracts(end) = (binedges(i+1)-stimesin(end))/binsize;
        end
        
        firates(i) = sum(fracts.*rates);
        if isnan(firates(i))
            ERRRRRORRRRRR
        end
    end
end


% for diagnostics
%     figure(1); clf; hold on; plot(binedges,ones(length(binedges)),'.k'); plot(stimes,1.1*ones(length(stimes)),'*b'); plot(binedges(1:end-1)+binsize/2,firates,'r');
%     k = waitforbuttonpress;
%     end
% end