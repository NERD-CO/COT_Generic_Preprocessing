% Function to switch target indices from "1" on the right, going
% counterclockwise to "1" on the left clockwise

function newdirs = TargetIdxSwap(olddirs)
    newdirs = olddirs;
    newdirs(olddirs==1) = 5;
    newdirs(olddirs==2) = 4;
    newdirs(olddirs==4) = 2;
    newdirs(olddirs==5) = 1;
    newdirs(olddirs==6) = 8;
    newdirs(olddirs==8) = 6;
end