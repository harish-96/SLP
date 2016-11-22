function f = fftfreq(z, Fs, Fcen)
%FFTFREQ f = fftfreq(z) OR  f = fftfreq(z, Fs, Fcen)
% f = fftfreq(z,Fs)
% f = fftfreq(z,Fs, Fcen);
% *** If single-argument input is used: ***
%	assumes that z is an FFT vector, and produces a vector of bin
%	frequencies of the same length.  Frequency is measured
%	in cycles per sampling period.  Frequencies are expressed
%	as being between -1/2 and +1/2.
% *** If two-argument input is used: ***
% Converts frequency scale to cycles per time unit (e.g., Hertz) using
% "Fs" as the sample rate. (i.e., scales the f vector by Fs). 
% *** If three-argument input is used: ***
% "Aliases" frequency representation to a band centered at "Fcen"
% (cycles/unit-time). 

M=length(z);
f=(0:(M-1))/M;
f=f-round(f);
f=f';
if nargin == 2
    f = f * Fs;
elseif nargin == 3;
    samplePeriod = 1/Fs;
    fshft = Fcen * samplePeriod;
    nfshft = round(fshft);
    f = f + nfshft;
    f( f<(fshft-0.5) ) = f( f<(fshft-0.5) ) + 1;
    f( f>=(fshft+0.5) ) = f( f>=(fshft+0.5) ) - 1;
    f = f * Fs;
end