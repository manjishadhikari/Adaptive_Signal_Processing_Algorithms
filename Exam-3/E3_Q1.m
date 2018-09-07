clear;
%close all
load('V:\EECS-844\Exam-3\P1.mat');

M=19;     %Length of each snapshot
K=length(d);
N=K-M+1;

X=complex(zeros(M,N));
D=complex(zeros(M,N));
for k=1:N
  X(:,k)=flipud(x(k:k+M-1));   %Snapshot matrix
  D(:,k)=flipud(d(k:k+M-1));   %Snapshot matrix
end

%% Using Linear Prediction

R=1/N*D*D';     %Correlation Matrix

r=zeros(M,1);    %Cross Correlation Matrix
for i=M:K-1
  r=r+ flipud(d(i-M+1:i)).*conj(d(i+1));
end
r=r/N;

wf=R\r;   %Weiner filter
am=[1; -wf];   %Linear prediction error filter


%% Inverse system ID using d as input and x as output
M2=M+1;
N2=K-M2+1;
for k=1:N2
  D2(:,k)=flipud(d(k:k+M2-1));   %Snapshot matrix
  
end
R2=1/N2*D2*D2';       %Correlation Matrix

P=zeros(M2,1);     %Cross Correlation Matrix
for i=M2:N2
  P=P+ flipud(d(i-M2+1:i)).*conj(x(i));
end
P=P/N2;
w=R2\P;     %Filter using Weiner Hopf Equation

[H_lpf,ang_lpf]=freqz(am,1,512,'whole');
[H_w,ang_w]=freqz(w,1,512,'whole');

%% Plotting Filter Properties
figure(1);

subplot(2,1,1); plot([1:M+1],20*log10(abs(am)));
hold on; plot([1:M+1],20*log10(abs(w)));  legend('LPE ','Weiner')
xlabel('Filter Coefficients');ylabel('Magnitude in dB')
title('Magnitude of  Filters')
subplot(2,1,2);plot([1:M+1],angle(am)*180/pi); title(' Phase of LPE Filter');
hold on;plot([1:M+1],angle(w)*180/pi,'r');  legend('LPE ','Weiner')
title(' Phase of Filters'); xlabel('Filter Coefficients');ylabel('Phase in deg')

figure(2); subplot(2,1,1); plot(ang_lpf./pi,20*log10(abs(H_lpf)))
hold on; plot(ang_w./pi,20*log10(abs(H_w))); legend('LPE ','Weiner')
legend('LPE','Weiner');  xlabel('Normalized Frequency');ylabel('Magnitude in dB')
title('Frequency Response of Filters')
subplot(2,1,2); plot(ang_lpf./pi,angle(H_lpf)*180/pi);
hold on; plot(ang_w./pi,angle(H_w)*180/pi);
legend('LPE ','Weiner');xlabel('Normalized Frequency');ylabel('Phase in deg')
title('Phase Response of  Filters')

y=am'*D2;   %Linear Prediction Filter output

figure;freqz(y);title('Linear Prediction FIlter Output freqz resp for M=9')
