%% Frequency-domian statistical properties of a mixed signal through the higher-order spectral density
%%
% O4681 
% © 2023. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
% National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
% Department of Energy/National Nuclear Security Administration. All rights in the program are.
% reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
% Security Administration. The Government is granted for itself and others acting on its behalf a
% nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
% derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
% others to do so.

% This program is Open-Source under the BSD-3 License.
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and
% the following disclaimer.
% 
% 2.Redistributions in binary form must reproduce the above copyright notice, this list of conditions
% and the following disclaimer in the documentation and/or other materials provided with the
% distribution.
% 
% 3.Neither the name of the copyright holder nor the names of its contributors may be used to endorse
% or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
% OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
% WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
% OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


%%

%Created by Neil Loychik, Los Alamos National Laboratory

% First part, create a mixture of two signals where the low frequency is
% Gausian, and the high-frequency is non-Gaussian. 

%The general statement is a random seed makes gaussian and non-Gaussian
%samples, which are applied to the time domain, and filtered. There are 5
%distinct zones:
%<20 Hz Noise
%20-180 Hz Gaussian
%180-220 Hz Mixing of Gaussian and Lognormal
%220-2000 Hz Lognormal 
%>3000 Noise

%With this model, we can a posterori claim that this is a mixture of
%Gaussian and non-Gaussian data bound by predetermined frequencies. This
%can therefore test the validity of a higher-order spectral density where
%it should show Gaussian behavior in the low-frequency morph into
%non-Gaussian in the high-frequency.

% close all
clear all
% close all
Fs=2000;
Tr=1000; %Primary variable, need to change the file names in lines 124-128 if saving
dF=1/Tr;

t=0:1/Fs:Tr-1/Fs; t=t';
x=randn(length(t),1);
% x=x';

x=interp(x,6); %Anti aliasing results by interpolation
Fs=6*Fs;
N=Fs;
t=0:1/Fs:Tr-1/Fs;
t=t';

x2=-lognrnd(0,1,[length(x),1]);
x2=x2-mean(x2);

% x2=x2*std(x)/std(x2);

% [std(x) std(x2)]
% return

%Initial band-limited data. 
[b,a] = butter(6,20/(Fs/2),'high'); %6pole buterworth at 10 Hz - Highpass
x = filtfilt(b,a,x);
x2 = filtfilt(b,a,x2);
[b,a] = butter(6,2000/(Fs/2)); %6pole buterworth at 2000 Hz - Lowpass
x = filtfilt(b,a,x);
x2 = filtfilt(b,a,x2);
% x=x+xs;
% x2=x2+xs;
clear b a

%Hear a sound clip
xclip1=x(1:Fs*9.9+1);
xclip2=x2(1:Fs*9.9+1);


%Combined gaussian, non-gaussian
x3a=x;
[b,a] = butter(6,20/(Fs/2),'high'); %6pole buterworth at 20 Hz - Highpass
x3a = filtfilt(b,a,x3a);
[b,a] = butter(6,200/(Fs/2)); %6pole buterworth at 200 Hz - Lowpass
x3a = filtfilt(b,a,x3a);
x3b=x2;
[b,a] = butter(6,180/(Fs/2),'high'); %6pole buterworth at 200 Hz - Highpass
x3b = filtfilt(b,a,x3b);
[b,a] = butter(6,2000/(Fs/2)); %6pole buterworth at 2000 Hz - Lowpass
x3b = filtfilt(b,a,x3b);
x3a=x3a*1.05; %Scale to be 0.001 g^2/Hz
x3b=x3b*1.2; %Scale to be 0.001 g^2/Hz
x3=x3a+x3b;

[GxxW(:,1),f]=pwelch(x3a,hann(N),[],N,Fs);
[GxxW(:,2),f]=pwelch(x3b,hann(N),[],N,Fs);
[GxxW(:,3),f]=pwelch(x3,hann(N),[],N,Fs);

clear b a x3a x3b

% Engage code to hear the sound clips and notice character difference in
% noise, though pitch and volume are identical. 

% xclip3=x3(1:Fs*9.9+1);
% tclip=1:length(xclip1); tclip=((tclip-1)/Fs)'; 
% sound(xclip1/std(x)/20,Fs) %Gaussian
% pause(12)
% sound(xclip2/std(x2)/20,Fs) %Non-Gaussian
% pause(12)
% sound(xclip3/std(x3)/20,Fs) %Mix

%% Probability Densities

%Graph a clip
xclip1=x(1:Fs*10+1);
xclip2=x2(1:Fs*10+1);
xclip3=x3(1:Fs*10+1);
tclip=1:length(xclip1); tclip=((tclip-1)/Fs)'; 

[xaaa,faaa]=ksdensity(x);
[xaaa2,faaa2]=ksdensity(x2);


%% Spectral Densities

%Octspace for octave spacing as opposed to lognormal
%https://www.mathworks.com/matlabcentral/fileexchange/70448-octspace
Fbe=octspace(10,6000,12); %Octave spacing to optimize random and bias error as a function of frequency. Generally 8-16 works well in most situation. 
Fbe=Fbe.center;

[M, F, Gx, Sx, Kx] = loychikSD_NP( x3, Fs, Fbe); %Integrated spectral densities
Gxx=gradient(Gx,F); %Power Spectral Density
Sxx=gradient(Sx,F); %Skewnesss Spectral Density
Kxx=gradient(Kx,F); %Kurtosis Spectral Density

%Note: You can calculate a central moment spectral density
%Ex. M4xx = gradient(Kx.*Gx.^2,F); %= Fourth Central Moment Spectral Density

clear x x2 x3 

save KSDMixedSignal1000-1.mat %Convenient to save long recordings and re-load for graphs

%% Graphs

load KSDMixedSignal1000-1.mat

figure(1)
subplot(1,2,1)
loglog(f,GxxW(:,3),'LineWidth',1.5)
ylabel('PSD (g^2/Hz)')
xlabel('Frequency (Hz)')
title('PSD of Mixed Signals')

subplot(1,2,2)
loglog(f,GxxW(:,1:2),'LineWidth',1.5)
ylabel('PSD (g^2/Hz)')
xlabel('Frequency (Hz)')
title('PSD from Different Seeds at Different Frequencies')
legend('Gaussian Process','Lognormal Process','location','sw')
% The signal is comprised of Gaussian and non-Gaussian components, which
% are not distinguishable from the PSD.

figure(2)
clf
subplot(2,2,1:2)
hold on
plot(faaa,xaaa,'LineWidth',2,'Color',[.4 .4 .4])
plot(faaa2,xaaa2,'LineWidth',2,'Color',[.7 .7 .7])
hold off
xlim([-20 10])
set(gca,'YScale','log')

ylim([.000001 1.1])
legend('Gaussian','Non-Gaussian')
ylabel('Probability Density (log-scale)')
xlabel('Sample')
set(gca, 'fontsize', 18)

subplot(2,2,3)
plot(tclip,xclip1,'Color',[.4 .4 .4])
xlim([0 2])
ylim([-10 10])
title('Gaussian')
xlabel('Time (s)')
ylabel('Amplitude')
set(gca, 'fontsize', 18)
subplot(2,2,4)
plot(tclip,xclip2,'Color',[.7 .7 .7])
xlim([0 2])
ylim([-10 10])
title('Non-Gaussian')
xlabel('Time (s)')
ylabel('Amplitude')
set(gca, 'fontsize', 18)

% Probability densities of the seed functions and associated waveforms
% showing the Gaussian and non-Gaussian component behavior. The
% non-Gaussian signal has a negative skew and high kurtosis. 

figure(3)
loglog(F,Gxx,f,GxxW(:,3),'LineWidth',1.5)
legend('Loychik PSD','Fourier-Based (Welch) PSD','location','sw')
title('PSD')
xlabel('Frequency (Hz)')
ylabel('PSD (V^2/Hz)')
xlim([10 5000])


hf=figure(4)
semilogx(F,Sxx,'LineWidth',1.5,'Color',[0 0.4470 0.7410])
ha=annotation('textbox')
ha.Parent = hf.CurrentAxes;


title('SSD')
xlabel('Frequency (Hz)')
ylabel('SSD (-/Hz)')
xlim([10 5000])
ylim([-.004 .004])
%After stabelizing (when variance is not near-zero) at ~30 Hz, skewness
%(Sx) stabelizes to 0 and the SSD reports 0 in the Gaussian regime. When
%the non-Gaussian signal is injected, skewness decreases as expected as
%shown by a negative SSD from 300-2000 Hz. 


figure(5)
semilogx(F,Kxx,'LineWidth',1.5)
title('KSD')
xlabel('Frequency (Hz)')
ylabel('KSD (-/Hz)')
xlim([10 5000])
ylim([-.015 .015])
%After stabelizing, Gaussian shows no change in kurtosis (KSD=0) with
%frequencies from 30-200 Hz. In the non-Gaussian regime, kurtosis increases
%with frequency, inditige non-Gaussian behavior is increasing. 


%Observations:
%1) The PSD from the new theory and convention (Welch) are near-identical
%with some difference in the noise regions. This implies a derivative
%approach is valid, and substantiates the higher-order spectral densities. 

%2) The SSD and KSD have very eradic values below 23 Hz. This is because a
%normalized moment is calculated by dividing by variance. Variance in the
%noise regime is near-zero, so these values are invalid until sufficient
%variance exists in the signal. 

%3) A spectral desnity should be thought of the change in a statistic with
%frequency. Therefore, skewness and kurtosis do not change significantly in
%the Gaussian regime.

%4) Main Point
% Looking at the PSD alone, one would not know it is comprised of two
%different random processes. However, the higher-order spectral desnities
% show show decreases in skewness and increases in kurtosis with frequency.
% After the bandlimit, contributions drop to zero in the high-frequency
%noise regime. Thus, the higher-order spectral densities show where the
%non-Gaussian behavior come from as intended.




