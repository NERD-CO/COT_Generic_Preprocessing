function [gapktime, rawpos] = get_tartare_time_and_pos(K, kinname, likecutoff)
% Find the label in the DLC data
labi = find(strcmp(K.Labels,kinname));
rawpos = [K.Xpos(:,labi), K.Ypos(:,labi), K.Zpos(:,labi)];

% snip the low likelihood 
rawpos(K.Likelihood(:,labi)<Settings.likecutoff,:) = NaN;
disp(['Cut ' num2str(sum(K.Likelihood(:,labi)<likecutoff)) ' samples for low likelihood']);

gapktime = K.Time;

% Remove first two times to attempt to avoid start pauses
gapktime = gapktime(3:end);
rawpos = rawpos(3:end,:);

% Try to detect a 1s end pause and stop before it.
dgkt = diff(gapktime);
lastsamp = find(dgkt > 0.98 & dgkt < 1.02);
if ~isempty(lastsamp)
    if lastsamp > length(gapktime)*0.9
        gapktime = gapktime(1:lastsamp);
        rawpos = rawpos(1:lastsamp,:);
    end
end
disp('Detected end pause, clipping.');

