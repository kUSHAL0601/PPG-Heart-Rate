function [ curLoc, curBPM, Trap_count, Lssr ] = JOSS( sigWindow,prevLoc, prevBPM, Trap_count, N, Fs, Lssr )
[~, SSRsig] = JOSS_STEP_1(sigWindow, Fs, N);
Lssr = length(SSRsig);
[curLoc, curBPM, Trap_count] = JOSS_STEP_2(SSRsig,Fs, prevLoc, prevBPM, Trap_count);
end