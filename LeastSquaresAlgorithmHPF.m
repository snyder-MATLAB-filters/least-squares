clear all;close all;clc;
N=55
q=(N-1)/2
r=q+1
fp=2400 %passband cutoff frequencies (Hz).
fs=2000 %stopband cutoff frequencies (Hz).
fsampling=10000 %sampling rate (Hz).
deltap=0.01 %passband tolerance.
deltas=0.01 %stopband tolerance.
thetap=2*pi*fp/fsampling
thetas=2*pi*fs/fsampling
density=16
L=r*density
Df=fsampling/2/L
Np=floor(fp/Df)
Ns=ceil(fs/Df)
ip=Np:L;
is=0:Ns;
fsp=Df*ip;
fss=Df*is;
thetap=ip*pi/L;
thetas=is*pi/L;
dp=ones(1,L-Np+1); %Desired in the passband.
ds=zeros(1,Ns+1); %Desired in the stopband.
Wp=diag(ones(1,L-Np+1)); %Wp = W in the passband.
Ws=diag(ones(1,Ns+1)*deltap/deltas);%Ws = W in the stopband.
Ap=cos(thetap'*[0:q]);
As=cos(thetas'*[0:q]);
A=[Ap;As];
W=diag([diag(Wp);diag(Ws)]);
R=A'*W*A;
d=[dp,ds]'; 
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
plot(theta/pi,Ha);grid on; 
xlabel('\theta/\pi'); 
ylabel('|Ha(\theta)|'); 
title('Amplitude Response (Linear)');

figure(3);
plot(theta/pi,abs(H));grid on;
xlabel('\theta/\pi'); 
ylabel('|H(\theta)|'); 
title('Magnitude Response (Linear)');

figure(4);
plot(theta/pi,20*log10(abs(H)));grid on;
xlabel('\theta/\pi');
ylabel('20 log_1_0(|H(\theta)|)'); 
title('Magnitude Response (dB)');

figure(5);
plot(theta/pi,angle(H)*180/pi);grid on; 
xlabel('\theta/\pi'); 
ylabel('\angleH(\theta) (deg)'); 
title('Phase Response');