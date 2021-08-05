%Cosine wave plus noise
clear all;format long;
N=16 %N = Number of taps.
var=5
mu=0.0001
N1=512
k=1:N1;
f1=1000;f2=1500;fs=50000;

s=5*cos(2*pi*f1*(k-1)/fs).'+5*cos(2*pi*f2*(k-1)/fs).';
w=sqrt(var)*randn([1,N1]).';

xa=s+w;
d=s;
w=zeros(1,N)';
x=[zeros(1,N-1),xa.'].';

for i=1:N1
x1=flipud(x(i:i+N-1));
y(i)=w'*x1;
e(i)=d(i)-y(i);
w=w+mu*conj(e(i))*x1;
end

theta=0:pi/(N1-1):pi;

figure(1);
plot(k,s);grid;
xlabel('k');ylabel('s(k)');
title('Signal');

figure(2);
plot(k,xa);grid;
xlabel('k');ylabel('x(k)');
title('Signal+Noise');

figure(3);
plot(theta/pi,10*log(abs(fft(xa)).^2));grid;
xlabel('\theta/\pi');ylabel('10log|fft(x(k))^2|');
title('Mag Squared of x(k) (dB)');

figure(4);
plot(k,y);grid;
xlabel('k');ylabel('y(k)');
title('Output');

figure(5);
plot(theta/pi,10*log(abs(fft(y)).^2));grid;
xlabel('\theta/\pi');ylabel('10log|fft(y(k))^2|');
title('Mag Squared of y(k) (dB)');

figure(6);
plot(k,e);grid;
xlabel('k');ylabel('e(k)');
title('Error');