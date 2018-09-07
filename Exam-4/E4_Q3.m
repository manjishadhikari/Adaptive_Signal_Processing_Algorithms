clear;
%close all
load('V:\EECS-844\Exam-4\P1.mat');

M=60;     %Length of each snapshot
K=length(x);
N=K-M+1;

X=complex(zeros(M,N));

for k=1:N
  X(:,k)=flipud(x(k:k+M-1));   %Snapshot matrix
end

R=1/N*X*X';       %Correlation Matrix

P=zeros(M,1);     %Cross Correlation Matrix
for i=M:K
  P=P+ flipud(x(i-M+1:i)).*conj(d(i));
end
P=P/N;
w_opt=inv(R)*P;      %Optimum Weiner filter

%% RLS

lambda_vect=[0.9998 0.999 0.99];

for lidx=1:length(lambda_vect)
  lambda=lambda_vect(lidx);
  P=eye(M);
w=zeros(M,1);
for n=M:K
  PI=P*(flipud(x(n-M+1:n)));
  k=PI/(lambda+(flipud(x(n-M+1:n))'*PI));
  error=d(n)-w'*(flipud(x(n-M+1:n)));
  w=w+k*conj(error);
  P=P/lambda-1/lambda*(k*(flipud(x(n-M+1:n))'*P));
  squared_error(n)=(abs(error)^2);
  dev(n)=(w-w_opt)'*(w-w_opt);
end

figure(1);hold on;plot(20*log10(dev(M:K)));
figure(2);hold on;plot(20*log10(squared_error(M:K)));
end
figure(1);title('Squared Deviation RLS');
 xlabel('Number of Iterations');
 ylabel('Squared deviation in dB')
figure(2);title('Squared Error RLS');
 xlabel('Number of Iterations');
 ylabel('Squared Error in dB')

figure(1);legend('lambda=0.9998','lambda=0.999','lambda=0.99')
figure(2);legend('lambda=0.9998','lambda=0.999','lambda=0.99')
