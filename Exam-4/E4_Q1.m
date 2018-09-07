
clear all;
close all
load('V:\EECS-844\Exam-4\P1.mat');

M=60;         %Length of each snapshot
K=length(x);  %Length of data
N=K-M+1;      %Number of Snapshots

X=complex(zeros(M,N));

for k=1:N
  X(:,k)=flipud(x(k:k+M-1));   %Snapshot matrix
end

R=1/N*(X*X');       %Correlation Matrix

P=zeros(M,1);     %Cross Correlation Matrix
for i=M:K
  P=P+ flipud(x(i-M+1:i)).*conj(d(i));
end
P=P/N;
w_opt=inv(R)*P;      %Optimum Weiner filter

%% NLMS

mu=0.5;           %Step Size
del=0.02;         %Leakage Factor

w=zeros(M,1);
dev=zeros(K,1);
squared_error=zeros(K,1);
for n=M:K
  error=d(n)-w'*flipud(x(n-M+1:n));    %Error
  mu_normalized=mu/(del+ctranspose(x(n-M+1:n))*x(n-M+1:n));
  w=w+mu_normalized*conj(error)*flipud(x(n-M+1:n));
  dev(n)=(w-w_opt)'*(w-w_opt);
  squared_error(n)=abs(error)^2;
end
 
 figure(1);plot(20*log10(dev(M:K)));title('Squared Deviation')
 xlabel('Number of Iterations');
 ylabel('Squared deviation in dB')
 figure(2);plot(20*log10(squared_error(M:K)));title('Squared Error')
  xlabel('Number of Iterations');
 ylabel('Squared Error in dB')