% Meander Analysis
% Diagnostic to forcast frequency cutoff.
% Determine width of plume to measure diffusion length and calculate mean
% wind to compute frequency cuttoff between diffusion and meandering.
% Mean Wind/width of plume = Fequency cutoff
%
% Date Author Comment
% 05-18-2011 Luna Rodriguez and George Young

% I. Housekeeping
clear all
close all
clc
% II. Load wind data
load ('output_Wind_Trial_54.mat');

% III. Fill in hard-coded variables
samplerate = 10;
Duration = 10*60;
y = WindDirCen;
y = detrend(y);
n = length(y);
DetrendLine = WindDirCen - y;
SweepRate = ((max (DetrendLine) - min(DetrendLine))/length(DetrendLine))*samplerate;%degrees per second
SweepAngle = Duration *SweepRate;
nfreq = (2^ceil(log2(n)))/2;
% III.
Y = fft(y,n);
Pyy = Y.*conj(Y)/n;
f = samplerate/n*(0:nfreq-1);
%.III.Create power spectra figure
figure
plot(log10(f),log10(Pyy(1:nfreq)))
axis([ -3.5 1.0 -5 5])
title('Wind Direction Power spectral density')
xlabel('Frequency (Hz)')
% IV.1
WindSpeed = sqrt((((U_WindsCen(:,1)).^2))+((V_WindsCen(:,1)).^2));
meanWindSpeed = mean(WindSpeed);
x = WindSpeed;
x = detrend(x);
NN = length(x);
DetrendLineWS = WindSpeed - x;
SweepRateWS = ((max (DetrendLineWS) - min(DetrendLineWS))/length(DetrendLineWS))*samplerate;%degrees per second
SweepAngleWS = Duration *SweepRateWS;
X = fft(x,NN);
Pxx = X.*conj(X)/NN;
figure
plot(log10(f),log10(Pxx(1:nfreq)))
axis([ -3.5 1.0 -6 3])
title('Wind Speed Power spectral density')
xlabel('Frequency (Hz)')
% IV.2 Calculate width of plume *** approximated***
widthPlume = 200;
% IV.3 Calculte cutoff frequnecy
cutfreq = meanWindSpeed/widthPlume;
% IV.4 Create weight matrix to apply to raw data
Fpass = 0.9 * cutfreq;
Fstop = 1.1 * cutfreq;
Dpass = 0.057501127785;
Dstop = .1;
dens = 20;
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(samplerate/2), [1,0],[Dpass, Dstop]);
b = firpm(N, Fo, Ao, W, (dens));
a = sum(b);
figure
freqz(b,a,1024); %
saveas(gcf,'freq.png')
yfilt = filtfilt(b,a,y);
xfilt = filtfilt(b,a,x);
figure
plot(y,'-b')
hold on
plot(yfilt,'-m','LineWidth',4)
title('Wind Direction')
xlabel('Time (tenths of seconds)')
ylabel('Wind Direction (degrees))')
hold off
figure
plot(x,'-b')
hold on
plot(xfilt,'-m','LineWidth',4)
title('Wind Speed')
xlabel('Time (tenths of seconds)')
ylabel('Wind Speed (m/s))')
hold off
% IV.5 Compute Statistics
stdevWindDir = std(y)+SweepAngle;
stdevWindDirFil = std(yfilt)+SweepAngle;
stdevWindSpd = std(x)+SweepAngleWS;
stdevWindSpdFil = std(xfilt)+SweepAngleWS;
save (['analysis_output_Wind_Trial_' num2str(trial) '.mat'])
