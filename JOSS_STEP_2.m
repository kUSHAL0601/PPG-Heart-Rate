function [curLoc, curBPM, Trap_count] = JOSS_STEP_1(SSRsig, Fs,prevLoc, prevBPM, Trap_count)
delta1 = 15;
delta2 = 25;
N = length(SSRsig);

if ((prevLoc == -1 && prevBPM == -1))
curLoc = find(SSRsig == max(SSRsig));
curBPM = 60 * (curLoc - 1) / N * Fs;
else

R1 = (prevLoc - delta1) : (prevLoc + delta1);

N = length(SSRsig);
H = zeros(N,1);

H(R1(R1 > 0)) = 1;

SSRsig = SSRsig.*H;

[~,locs] = findpeaks(SSRsig);

numOfPks = length(locs);
if (numOfPks >= 1)

[~,temp_index] = min(abs(prevLoc - locs));
curLoc = locs(temp_index);
curBPM = 60 * (curLoc - 1) / N * Fs;
else

R2 = (prevLoc - delta2) : (prevLoc + delta2);

N = length(SSRsig);
H = zeros(N,1);

H(R2(R2 > 0)) = 1;

SSRsig = SSRsig.*H;

[~,locs] = findpeaks(SSRsig);

numOfPks = length(locs);
if (numOfPks >= 1)

[~, maxIndex] = max(pks);
curLoc = locs(maxIndex);
curBPM = 60 * (curLoc - 1) / N * Fs;
else

curLoc = prevLoc;
curBPM = prevBPM;

end
end
end
% validation
if (curLoc == prevLoc)
Trap_count = Trap_count + 1;
if (Trap_count > 2)
R3 = 40:200;

N = length(SSRsig);
H = zeros(N,1);

H(R3(R3 > 0)) = 1;

SSRsig = SSRsig.*H;

[pks,locs] = findpeaks(SSRsig);

[~,temp_index] = min(abs(prevLoc - locs));
curLoc = locs(temp_index);
curBPM = 60 * (curLoc - 1) / N * Fs;
end
else
Trap_count = 0;
end
end
