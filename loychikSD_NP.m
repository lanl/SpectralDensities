function [M, F, Gx, Sx, Kx] = loychikSD_NP( x, Fs, Fbe )
%%
%O4680 
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


%Identical from loychikSD, except remove PDSD for faster computation.

%By: Neil Loychik
%Last Edited: 09/14/2023

%ABSTRACT: This script calculates a spectral density in accordance with 
%"Spectral Densities: Statistics and Probability in the Frequency Domain."
%This algorithm was the subject of a patent filing, and after the filing
%released as part of a conference presentation and paper. 

%As described in the patent, the novelty of this method is the ability to
%compute higher-order spectral densities. Namely, the current "state of the
%art" is to compute a power spectral density (PSD) with algorithms like
%MATLAB's pwelch script, which only calculates a spectral density of power:
%[Gxx F]=pwelch(x,...);

%In comparison, this script calculates spectral densities of power
%(variance) and higher-order statistics of skewness, kurtosis, and the
%probability density. E.g.:

%[M, F, Gx, Sx, Kx, Px] = loychikSD( x, Fs, Fbe, xbin )
%Gxx(:,ii)=gradient(Gx,F); %PSD
%Sxx(:,ii)=gradient(Sx,F); %SSD
%Kxx(:,ii)=gradient(Kx,F); %KSD
%Pxx(:,ii)=gradient(Px,F); %PDSD

%The PSD (Gxx) calculated by either Welch's of this algorithm should be 
%equivalent. 

%This script outlines a functional approach and is not optimized for
%computational efficiency. Notable efficiency improvements would be
%parrallelization, filtering, variable inputs, and optimized bandwidth
%error.

%Due to computational intensity, this script only computes a single
%array (e.g. x(:,1) ) of data instead of an matrix (e.g. x(:,:)). 




% LOYCHIKSD
%Output
%M = Mean and standard moments of the data. Useful to justify unitary
%    [mean(x); var(x), skewness(x), kurtosis(x), moment(x,1); moment(x,2);
%     moment(x,3)]
%F - Frequency for the spectral density. Bin center.
%Gx - Integrated Power Spectral Density - PSD - Results are comparable to
%     the integral of the pwelch function in Matlab or the Daniel method.
%Sx - Integrated Skewness Spectral Density - SSD
%Kx - Integrated Kurtosis Spectral Density - KSD
%Px - Integrated Probability Density Spectral Density - PDSD


%Input
% x = data
% Fs = Sampling Frequency
% Fbe = Frequency bin edges. Easier to implement by defining the edges and
% computing the bin center F, than defining F and computing edges.
%xn = bins for the probability density



%%

F=Fbe(1:end-1)+diff(Fbe); %Frequency bin centers

X=fft(x);
t=1:length(x); t=t-1; t=t/Fs; t=t';
dF=1/t(end);
Ff = (0:length(x)-1)*dF; Ff=Ff';%FFT frequency step

M=[mean(x); var(x); skewness(x); kurtosis(x); moment(x,1); moment(x,2); moment(x,3)];

for ii=1:length(F)
    [ii length(F)] %Counter to determine where you are in run
    %Begin low-pass filter by zeroing bins
    Xn=find(Ff>=Fbe(ii+1),1,'first');
    
    cutX=X(2:Xn,:);
    XI=zeros([size(X,1) size(X,2)]);
    XI(2:Xn,:)=cutX;
    XI=flipud(XI);
    XI(1:Xn-1,:)=conj(cutX); %Removed symetric tag on ifft and replaced with complex conjugate to make code more similar to Python.  8/4/2023. 
    XI=flipud(XI);
    
    xifft=ifft(XI);
    %End low-pass filter by zeroing bins
    
    Gx(ii,1)=var(xifft);
    Sx(ii,1)=skewness(xifft);
    Kx(ii,1)=kurtosis(xifft);
    M3x(ii,1)=moment(xifft,3);
    M4x(ii,1)=moment(xifft,4);
end



