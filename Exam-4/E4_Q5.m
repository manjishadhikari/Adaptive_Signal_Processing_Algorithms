%Yule Walker Method for PSD Estimation

clear all;
% close all;
load('V:\EECS-844\Exam-4\P4.mat');
pn=[2 4 8];   %AR model order

for idx=1:length(pn)
 clearvars -except pn x idx
  p=pn(idx);
  M=p+1;    %No of AR coefficients
  K=length(x);
  N=K-M+1;
 
  for k=1:N
    X(:,k)=flipud(x(k:k+M-1));  %Snapshot matrix
  end
  R=1/N*X*X';    %Sample Correlation matrix from Data
  rxx=R(2:M,1);   %Autocorrelation 
  Rxx=R(1:p,1:p);  
  a=-inv(Rxx)*rxx;
  ar=[1;conj(a)];   %AR(p) model filter coefficients
  w=linspace(0,2*pi,1000);
  [H,W]=freqz(1,ar,1000,'whole');
  figure(1);hold on;plot(w,20*log10(abs(H)));
end
title('Yule Walker Method PSD Estimation');
xlabel('Normalized frequency rad/sample');
ylabel('Normalized Power in dB');
legend('p=2','p=4','p=8')
