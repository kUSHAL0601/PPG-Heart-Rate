function [ BPMF ] = FindBPM( signal,time_left )
L=size(signal,2);
tic;
len=1;
BPM_approx=Wiener(signal);
BPMF=BPM_approx;
h = [1/2 1/2];
binomialCoeff = conv(h,h);
for n = 1:4
    binomialCoeff = conv(binomialCoeff,h);
end
temp1 = filter(binomialCoeff, 1, BPMF);
temp2=temp1(4:end);
temp2(1:3)=BPMF(1:3);
temp2=[temp2 BPMF(length(BPMF)-2:end)];
BPMF=temp2;
BPM_approx=BPMF;
d=abs(diff(BPM_approx));
d=(d>=5);
for i=3:2:floor(L/125)-7
    len=len+1;
    time_now=toc;
    if(time_now>=time_left-200)
        break;
    end
    if len>125
	break
    end
    if(d(len-1))
        [~ ,curBPM ,~ ,~] = JOSS(signal(:,(i-1)*125+1:(i+7)*125), -1,BPM_approx(len-1),0,60*125,125,0);
        if(abs(curBPM - BPM_approx(len-1)) < abs(BPM_approx(len) - BPM_approx(len-1)))
            BPMF(1,len)=curBPM;
        end
    end
end

h = [1/2 1/2];
binomialCoeff = conv(h,h);
for n = 1:4
    binomialCoeff = conv(binomialCoeff,h);
end
BPMF=BPMF(1,1:125);
end
