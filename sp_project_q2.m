

fileList = { "File1.wav", "File2.wav", "File3.wav", "File4.wav" };

winLen  = 1024;
hopLen  = 256;
nfft    = 1024;
window  = hann(winLen, "periodic");

K = 3;

energySmoothWinFrames = 5;
minPeakDistance_sec    = 0.05;
energyThresholdFactor  = 1.0;
decayRatio = 0.1;

for f = 1:numel(fileList)
    fprintf("\n=====================================\n");
    fprintf("Processing %s\n", fileList{f});
    fprintf("=====================================\n");

   
    [audio, fs] = audioread(fileList{f});
    if size(audio,2) > 1, audio = mean(audio,2); end
    audio = audio / (max(abs(audio)) + 1e-8);

    
    [S, F, T] = spectrogram(audio, window, winLen-hopLen, nfft, fs);
    magSpec = abs(S);
    frameEnergy = sum(magSpec.^2, 1);

    smoothEnergy = movmean(frameEnergy, energySmoothWinFrames);
    smoothEnergy = smoothEnergy / (max(smoothEnergy) + 1e-8);

    thr = mean(smoothEnergy) + energyThresholdFactor * std(smoothEnergy);
    minPeakDist_frames = round(minPeakDistance_sec / (hopLen/fs));

    [~, pkLocs] = findpeaks(smoothEnergy, ...
        "MinPeakHeight", thr,...
        "MinPeakDistance", minPeakDist_frames);

    onsetFrames = pkLocs;
    onsetTimes = T(onsetFrames);
    numHits = numel(onsetTimes);
    fprintf("Detected %d drum hits.\n", numHits);

 
    durations = zeros(numHits,1);
    for i = 1:numHits
        startIdx = onsetFrames(i);
        peakVal = smoothEnergy(startIdx);
        cutoff = peakVal * decayRatio;

        j = startIdx;
        while j <= length(smoothEnergy) && smoothEnergy(j) > cutoff
            j = j+1;
        end
        endIdx = min(j, length(smoothEnergy));
        durations(i) = max(T(endIdx) - T(startIdx), 0);
    end

    
    [stftMag,freqVec,timeVec] = stft(audio,fs,"Window",window,...
        "OverlapLength",winLen-hopLen,"FFTLength",nfft);
    stftMag = abs(stftMag);

    features = zeros(numHits,6);
    for i = 1:numHits
        [~, idx] = min(abs(timeVec - onsetTimes(i)));
        idxStart = max(1, idx-1);
        idxEnd = min(size(stftMag,2), idx+1);

        patch = mean(stftMag(:,idxStart:idxEnd),2);
        patch = patch + 1e-8;

        centroid = sum(freqVec .* patch) / sum(patch);
        bandwidth = sqrt(sum(((freqVec-centroid).^2).*patch)/sum(patch));
        cumMag = cumsum(patch);
        rolloffIdx = find(cumMag>=0.85*cumMag(end),1);
        rolloffFreq = freqVec(rolloffIdx);
        logEnergy = log(sum(patch.^2));
        splitFreq = 1000;
        hi = sum(patch(freqVec>splitFreq));
        lo = sum(patch(freqVec<=splitFreq));
        hiLowRatio = hi / (lo + 1e-8);

        features(i,:) = [centroid, bandwidth, rolloffFreq, logEnergy, hiLowRatio, durations(i)];
    end

   
    if numHits >= K
        [idxCluster,~] = kmeans(features,K,"Replicates",10);
        K_used = K;
    else
        idxCluster = ones(numHits,1);
        K_used = 1;
    end

    
    fprintf("Number of instruments: %d\n", K_used);
    for k = 1:K_used
        idx = find(idxCluster == k);
        fprintf(" Cluster %d: %d hits => ", k, numel(idx));
        fprintf("%.2f ", onsetTimes(idx));
        fprintf("\n");
    end

    
    t = (0:length(audio)-1)/fs;

    figure("Name", fileList{f}, "NumberTitle","off");
    subplot(3,1,1); plot(t,audio); title("Waveform"); xlabel("Time"); ylabel("Amp"); grid on;
    subplot(3,1,2); plot(T,smoothEnergy); hold on;
    yline(thr,"--r"); stem(onsetTimes,smoothEnergy(onsetFrames),"filled");
    title("Energy + Hits"); xlabel("Time"); ylabel("Energy"); grid on;
    subplot(3,1,3);
    imagesc(T,F,20*log10(magSpec+1e-8)); axis xy; colormap jet; colorbar;
    hold on; scatter(onsetTimes, max(F)*0.9*ones(size(onsetTimes)),30,idxCluster,"filled");
    title("Spectrogram"); xlabel("Time"); ylabel("Freq (Hz)");
end

fprintf("\nDone.\n");