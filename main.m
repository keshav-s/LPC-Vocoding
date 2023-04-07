%% LPC Vocoding:
% We generate LPC coefficients and gain estimates for windows of a given
% speech signal. We then vocode the speech signal by filtering an
% excitation with the LPC coefficients/gain for a frame. We introduce
% robustness by ensuring 50% overlap between frames so our coefficients
% the output has twice the LPC coefficients used per frame.

%% Generate LPC Analysis Frames, apply LPC to each frame
audio = audioread('welcome16k.wav');
fs = 16000; % sampling frequency of audio file
nFFT = 512; % length of fft
pulse_train = zeros(length(audio),1);
pulse_train(40:160:end) = 1;
white = wgn(length(audio),1,0); %0dBW power => unity power

frame_len = fs * 20/1000; % 20 millisecond frame
h = hamming(frame_len);
overlap = frame_len/2;
num_frames = floor((length(audio) - overlap)/(overlap));
frames = zeros(num_frames,frame_len);
pframes = zeros(num_frames, frame_len);
wframes = zeros(num_frames, frame_len);
for i = 0:num_frames-1
    arange = i*overlap + 1: i*overlap+frame_len;
    frames(i+1,:) = h.*audio(arange);

    % Note that we generate the 100Hz excitation and white noise frames 
    % here in order to preserve appropriate frame alignment including
    % overlaps
    pframes(i+1,:) = pulse_train(arange);
    wframes(i+1,:) = white(arange);
end

P = 14;
as = zeros(num_frames, P+1); % LPC predictor coefficients per frame
gs = zeros(num_frames); % gain parameters per frame

for i = 1:num_frames
    [a,g] = lpc(frames(i,:), P);
    as(i,:) = a;
    gs(i) = sqrt(g);
end

figure(1);
colormap jet
% Plot spectogram with 320-sample long Hamming windows on each STFT frame
% 512 length FFT (closest power of 2 greater than 320)
% 160 sample overlap between frames (50% overlap)
spectrogram(audio,hamming(frame_len),overlap,nFFT,fs,'yaxis')

%% 100Hz Excitation
% To resynthesize, we have to take into account the 50% overlap between LPC
% frames. To resynthesize the frames, we take the LPC coeffs a(i) and gain G from
% every even frame i, apply the filter using the difference equation:
% G / (1 - a_1(i) - a_2(i) ...) to the 100Hz excitation. We do the same
% for the odd frames. We then offset the resynthesized signal from the odd
% frames by 160 samples (half a window) and add the offset signal to the
% resynthesized even signal to get our final output.
even_resynth = zeros(ceil(num_frames/2), frame_len);
odd_resynth = zeros(floor(num_frames/2), frame_len);
[y,zfe] = filter(gs(1,:), as(1,:), pframes(1,:));
even_resynth(1, :) = y;
[y,zfo] = filter(gs(2,:), as(2,:), pframes(2,:));
odd_resynth(1, :) = y;

zie = zfe;
zio = zfo;
for i = 3:2:num_frames
    [y,zfe] = filter(gs(i,:), as(i,:), pframes(i,:), zie);
    even_resynth((i-1)/2, :) = y;
    zie = zfe;
end
for i = 4:2:num_frames
    [y,zfo] = filter(gs(i,:), as(i,:), pframes(i,:), zio);
    odd_resynth(i/2, :) = y;
    zio = zfo;
end

evens = reshape(even_resynth.',1,[]);
odds = reshape(odd_resynth.',1,[]);
odds_padded = [zeros(1,overlap), odds, zeros(1,overlap)];

out100 = evens + odds_padded;
audiowrite('100HzExcite.wav', out100, fs);
figure(2);
colormap jet
spectrogram(out100,hamming(frame_len),overlap,nFFT,fs,'yaxis')

%% White Noise Excitation
even_resynth = zeros(ceil(num_frames/2), frame_len);
odd_resynth = zeros(floor(num_frames/2), frame_len);
[y,zfe] = filter(gs(1,:), as(1,:), wframes(1,:));
even_resynth(1, :) = y;
[y,zfo] = filter(gs(2,:), as(2,:), wframes(2,:));
odd_resynth(1, :) = y;

zie = zfe;
zio = zfo;
for i = 3:2:num_frames
    [y,zfe] = filter(gs(i,:), as(i,:), wframes(i,:), zie);
    even_resynth((i-1)/2, :) = y;
    zie = zfe;
end
for i = 4:2:num_frames
    [y,zfo] = filter(gs(i,:), as(i,:), wframes(i,:), zio);
    odd_resynth(i/2, :) = y;
    zio = zfo;
end

evens = reshape(even_resynth.',1,[]);
odds = reshape(odd_resynth.',1,[]);
odds_padded = [zeros(1,overlap), odds, zeros(1,overlap)];

outw = evens + odds_padded;
audiowrite('whiteNoise.wav', out100, fs);
figure(3);
colormap jet
spectrogram(outw,hamming(frame_len),overlap,nFFT,fs,'yaxis')