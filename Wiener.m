function [finalBPM] =Wiener(sig)
srate = 125;
FFTres = 125*60;   % FFT resolution - user tunable parm1
WFlength = 15;

CutoffFreqHzHP = 0.75;
CutoffFreqHzLP = 3.5;
[b,a] = butter(4, [0.4 4]/(125/2),'bandpass');

window   = 8 * srate;  % window length is 8 seconds
step     = 2 * srate;  % step size is 2 seconds

windowNb = floor((length(sig)-window)/step) + 1;  % total number of windows(estimates)

BPM_est = zeros(1,windowNb); % estimated BPM
rangeIdx = []; % range of search for the next estimates


for i =  [1 :  windowNb]
    curSegment = (i-1)*step+1 : (i-1)*step+window;
    curData = sig(:,curSegment);
    
    PPG1 = curData(1,:); PPG2 = curData(2,:);
    ACC_X = curData(3,:); ACC_Y = curData(4,:); ACC_Z = curData(5,:);
    
    % filtering
    PPG1 = filter(b,a,PPG1);
    PPG2 = filter(b,a,PPG2);
    ACC_X = filter(b,a,ACC_X);
    ACC_Y = filter(b,a,ACC_Y);
    ACC_Z = filter(b,a,ACC_Z);
    PPG_ave = 0.5*(PPG1-mean(PPG1))/(std(PPG1))+0.5*(PPG2-mean(PPG2))/(std(PPG2)); % mean overall
    
    % downsampling to 25Hz
    PPG_ave = downsample(PPG_ave,5);
    ACC_X = downsample(ACC_X,5);
    ACC_Y = downsample(ACC_Y,5);
    ACC_Z = downsample(ACC_Z,5);
    srate = 25; % new sampling rate
    
    % Periodogram
    PPG_ave_FFT = fft(PPG_ave,FFTres);
    FreqRange = linspace(0,srate,size(PPG_ave_FFT,2));
    
    % finding the indices for the range of interest
    [extra,lowR] = (min(abs(FreqRange-CutoffFreqHzHP)));
    [extra,highR] = (min(abs(FreqRange-CutoffFreqHzLP)));
    
    %  Getting rid of most spectra outside the range of interest
    FreqRange = FreqRange(lowR:highR);
    PPG_ave_FFT = PPG_ave_FFT(lowR:highR);
    ACC_X_FFT= fft(ACC_X,FFTres); ACC_X_FFT = ACC_X_FFT(lowR:highR);
    ACC_Y_FFT= fft(ACC_Y,FFTres); ACC_Y_FFT = ACC_Y_FFT(lowR:highR);
    ACC_Z_FFT= fft(ACC_Z,FFTres); ACC_Z_FFT = ACC_Z_FFT(lowR:highR);
    
    
    % phase vocoder to refine spectral estimations
    FreqRangePPG = FreqRange;
    if i>1 % start phase vocoder for current and previous frames
        for ii=1:size(FreqRangePPG,2)
            curPhase = angle(PPG_ave_FFT(ii));
            prevPhase = angle(PPG_ave_FFTpr(ii)); vocoder = zeros(1,20);
            for n = 1:20
                vocoder(n) = ((curPhase-prevPhase)+(2*pi*(n-1)))/(2*pi*2);
            end
            difference = vocoder - FreqRange(ii);
            [extra, deltaidx] = min(abs(difference));
            FreqRangePPG(ii) = vocoder(deltaidx);
        end
    end
    
    FreqRangePPG = moving(FreqRangePPG,3);
    
    PPG_ave_FFTpr = PPG_ave_FFT;
        

    WC1 = WFlength; WC2 = WFlength;
    
    W1_FFTi(i,:) = (abs(PPG_ave_FFT))/max(abs(PPG_ave_FFT));
    if i==1, W1_PPG_ave_FFT_ALL = W1_FFTi(i,:); else W1_PPG_ave_FFT_ALL = mean(W1_FFTi(max(1,i-WC1):i,:),1); end
    W1_PPG_ave_FFT_ALL_norm = (W1_PPG_ave_FFT_ALL)/max(W1_PPG_ave_FFT_ALL);
    W1_ACC_X_FFT_norm = (abs(ACC_X_FFT))/max(abs(ACC_X_FFT));
    W1_ACC_Y_FFT_norm = (abs(ACC_Y_FFT))/max(abs(ACC_Y_FFT));
    W1_ACC_Z_FFT_norm = (abs(ACC_Z_FFT))/max(abs(ACC_Z_FFT));
    WF1 = (1 - 1/3*(W1_ACC_X_FFT_norm+W1_ACC_Y_FFT_norm+W1_ACC_Z_FFT_norm)./(W1_PPG_ave_FFT_ALL_norm)); 
    WF1 (WF1<0) = -1; % limit negative -inf to -1
    W1_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*WF1;
    
    W11_FFTi(i,:) = abs(PPG_ave_FFT).^2;
    if i==1, W11_PPG_ave_FFT_ALL = W11_FFTi(i,:); else W11_PPG_ave_FFT_ALL = mean(W11_FFTi(max(1,i-WC1):i,:),1); end
    W11_PPG_ave_FFT_ALL_norm = (W11_PPG_ave_FFT_ALL)/max(W11_PPG_ave_FFT_ALL);
    W11_ACC_X_FFT_norm = (abs(ACC_X_FFT).^2)/max(abs(ACC_X_FFT.^2));
    W11_ACC_Y_FFT_norm = (abs(ACC_Y_FFT).^2)/max(abs(ACC_Y_FFT).^2);
    W11_ACC_Z_FFT_norm = (abs(ACC_Z_FFT).^2)/max(abs(ACC_Z_FFT).^2);
    WF11 = (1 - 1/3*(W11_ACC_X_FFT_norm+W11_ACC_Y_FFT_norm+W11_ACC_Z_FFT_norm)./(W11_PPG_ave_FFT_ALL_norm)); 
    W11_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*(WF11);
    
    W2_FFTi(i,:) = (abs(PPG_ave_FFT))/max(abs(PPG_ave_FFT));
    if i==1, W2_PPG_ave_FFT_ALL = W2_FFTi(i,:); else W2_PPG_ave_FFT_ALL = mean(W2_FFTi(max(1,i-WC2):i,:),1); end
    W2_PPG_ave_FFT_ALL_norm = (W2_PPG_ave_FFT_ALL)/max(W2_PPG_ave_FFT_ALL);
    W2_ACC_X_FFT_norm = (abs(ACC_X_FFT))/max(abs(ACC_X_FFT));
    W2_ACC_Y_FFT_norm = (abs(ACC_Y_FFT))/max(abs(ACC_Y_FFT));
    W2_ACC_Z_FFT_norm = (abs(ACC_Z_FFT))/max(abs(ACC_Z_FFT));
    WF2 = W2_PPG_ave_FFT_ALL_norm./(((W2_ACC_X_FFT_norm+W2_ACC_Y_FFT_norm+W2_ACC_Z_FFT_norm)/3)+W2_PPG_ave_FFT_ALL_norm); 
    W2_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*WF2;
    W2_FFTi(i,:) = (W2_PPG_ave_FFT_Clean(i,:))/max(W2_PPG_ave_FFT_Clean(i,:));
    
    W21_FFTi(i,:) = abs(PPG_ave_FFT).^2;
    if i==1, W21_PPG_ave_FFT_ALL = W21_FFTi(i,:); else W21_PPG_ave_FFT_ALL = mean(W21_FFTi(max(1,i-WC2):i,:),1); end
    W21_PPG_ave_FFT_ALL_norm = W21_PPG_ave_FFT_ALL/max(W21_PPG_ave_FFT_ALL);
    W21_ACC_X_FFT_norm = abs(ACC_X_FFT).^2/max(abs(ACC_X_FFT).^2);
    W21_ACC_Y_FFT_norm = abs(ACC_Y_FFT).^2/max(abs(ACC_Y_FFT).^2);
    W21_ACC_Z_FFT_norm = abs(ACC_Z_FFT).^2/max(abs(ACC_Z_FFT).^2);
    WF21 = W21_PPG_ave_FFT_ALL_norm./(((W21_ACC_X_FFT_norm+W21_ACC_Y_FFT_norm+W21_ACC_Z_FFT_norm)/3)+W21_PPG_ave_FFT_ALL_norm);
    W21_PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*WF21;
    W21_FFTi(i,:) = W21_PPG_ave_FFT_Clean(i,:).^2;
    
    
    W1_PPG_ave_FFT_Clean(i,:) = W1_PPG_ave_FFT_Clean(i,:)/std(W1_PPG_ave_FFT_Clean(i,:)); 
    W11_PPG_ave_FFT_Clean(i,:) = W11_PPG_ave_FFT_Clean(i,:)/std(W11_PPG_ave_FFT_Clean(i,:)); 
    W2_PPG_ave_FFT_Clean(i,:) = W2_PPG_ave_FFT_Clean(i,:)/std(W2_PPG_ave_FFT_Clean(i,:)); 
    W21_PPG_ave_FFT_Clean(i,:) = W21_PPG_ave_FFT_Clean(i,:)/std(W21_PPG_ave_FFT_Clean(i,:)); 

    
    
    PPG_ave_FFT_FIN(i,:) = W1_PPG_ave_FFT_Clean(i,:)+ W2_PPG_ave_FFT_Clean(i,:) ; % ensambling W1 & W2
    
    hist_int = 30; % We start with +- 25 BPM history tracking window
    % History tracking
    if i>30, hist_int = max(abs(diff(BPM_est(1:i-1))))+5; end; 
    
    % HR estimation
    if isempty(rangeIdx)
        [extra, idx]= max(PPG_ave_FFT_FIN(i,:));
        BPM_est(i) = FreqRangePPG(idx(1))*60; 
        rangeIdx = idx(1)-round(hist_int/((FreqRange(2)-FreqRange(1))*60)):idx(1)+round(hist_int/((FreqRange(2)-FreqRange(1))*60));
    else
        [extra, idx]= max(PPG_ave_FFT_FIN(i,rangeIdx));
        BPM_est(i) = FreqRangePPG(rangeIdx(idx(1)))*60; 
        rangeIdx = rangeIdx(idx(1))-round(hist_int/((FreqRange(2)-FreqRange(1))*60)):rangeIdx(idx(1))+round(hist_int/((FreqRange(2)-FreqRange(1))*60));
    end
    rangeIdx(rangeIdx<1) = []; rangeIdx(rangeIdx>length(FreqRange)) = [];
    
    
    % Mild smoothing with linear regression prediction
    if i>5 && abs(BPM_est(i)-BPM_est(i-1))>5
        %BPM_est(i) = BPM_est(i-1)+sign(BPM_est(i)-BPM_est(i-1))*10;
        ddd= polyfit(1:length(BPM_est(max(1,i-5):i-1)),BPM_est(max(1,i-5):i-1),1);
        BPM_est(i) = 0.8*BPM_est(i)+0.2*polyval(ddd,length(BPM_est(max(1,i-5):i-1))+1);
    end
    
    mul=0.1;
    BPM_est(i) = BPM_est(i)+sum(sign(BPM_est(max(2,i-6):i)-BPM_est(max(1,i-7):i-1))*mul);
end
finalBPM=BPM_est;
end
