function [FR] = calculate_ktime_rates(ktime, N, slidepad)

V.gkerns = [0.015, 0.05];

nneu = length(N.SpkTimes);
nkp = length(ktime);

% Do ms resolution FIR rates, add a 5 second buffer on the front and back
% to avoid edge effects when smoothing
msedge = ((ktime(1)-5):0.001:(ktime(end)+5))';
mslooktime = msedge(1:end-1)+0.0005;

msrates = nan(length(mslooktime), nneu);
for j = 1:nneu
    msrates(:,j) = Spikes2FIR_Arbitrary( N.SpkTimes{j}, msedge)';
end

for kerni = 1:length(V.gkerns)
    thisname = ['SmoothG' num2str(V.gkerns(kerni)*1000,'%0.2i')];
    FR.(thisname).krates = nan(nkp,nneu);
    for neui = 1:nneu
        FR.(thisname).krates(:,neui) = GaussSmooth_Arbitrary( mslooktime, msrates(:,neui), ktime, V.gkerns(kerni), 5);
    end

    padrates = [nan(slidepad,nneu); FR.(thisname).krates; nan(slidepad,nneu)];

    if any(any(isnan(FR.(thisname).krates)))
        OHNOO
    end

    % And do sliding renormalized and mean subtracted rates
    FR.(thisname).slidenormrates = nan(size(FR.(thisname).krates));
    FR.(thisname).slidemeanrates = nan(size(FR.(thisname).krates));
    for slidei = 1:nkp
        FR.(thisname).slidemeanrates(slidei,:) = (FR.(thisname).krates(slidei,:) - nanmean(padrates(slidei:(slidei+slidepad*2),:),1));
        FR.(thisname).slidenormrates(slidei,:) = FR.(thisname).slidemeanrates(slidei,:)./nanstd(padrates(slidei:(slidei+slidepad*2),:),0,1);
    end
end