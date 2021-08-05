close all;clc;clear all;
N=55
q=(N-1)/2
r=q+1

fs1 = 800
fps = 1200
fpe = 3000
fs2 = 3600

fsampling=10000 %sampling rate (Hz).
deltap=0.01 %passband tolerance.
deltas=0.01 %stopband tolerance.

thetas1=2*pi*fs1/fsampling
thetaps=2*pi*fps/fsampling
thetape=2*pi*fpe/fsampling
thetas2=2*pi*fs2/fsampling

density=16

L=r*density

i1 = floor(thetas1 / pi * L)
i2 = floor(thetaps / pi * L)
i3 = floor(thetape / pi * L)
i4 = floor(thetas2 / pi * L)

b1 = i1
b2 = i3-i2
b3 = L-i4



theta = pi / L
indexq = 0:q;
indexL = [(1:b1),(i2:i3),(i4:L)]';
indexM = indexL * indexq * theta;
A = cos(indexM); % get A

d=[ones(1,b1),zeros(1,b2+1),ones(1,b3+1)]';
Wpre = [ones(1,b1),ones(1,b2+2)*deltap/deltas,ones(1,b3)]
W=diag(Wpre);
% after here everything are the same
R=A'*W*A;
a=A'*W*d;
c=R\a; h=[flipud(c(2:end))/2;c(1);c(2:end)/2];
%h=h/sum(h);
k=0:N-1;
M=L;
theta=0:pi/M:pi;
Q=[0:q]'*theta;
Ha=c'*cos([0:q]'*theta);
H=h.'*exp(-1j*k'*theta); theta1=[0 1 2 3 4]*pi/4;
Ha1=h.'*exp(-1j*k'*theta1); Mag=abs(Ha1)
figure(1)
stem(k,h);grid on; xlabel('k');ylabel('h(k)');
title('Impulse Response');
figure(2);
plot(theta/pi,Ha);grid on; xlabel('\theta/\pi'); ylabel('|Ha(\theta)|'); title('Amplitude Response (Linear)');
figure(3);
plot(theta/pi,abs(H));grid on; xlabel('\theta/\pi'); ylabel('|H(\theta)|'); title('Magnitude Response (Linear)');
figure(4);
plot(theta/pi,20*log10(abs(H)));grid on;
xlabel('\theta/\pi'); ylabel('20 log_1_0(|H(\theta)|)'); title('Magnitude Response (dB)');
figure(5);
plot(theta/pi,angle(H)*180/pi);grid on; xlabel('\theta/\pi'); ylabel('\angleH(\theta) (deg)'); title('Phase Response');