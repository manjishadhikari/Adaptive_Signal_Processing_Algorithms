
clear;

load('V:\EECS-844\Exam-4\P1.mat');
M=60;     %Length of each snapshot
K=length(d);
N=K-M+1;

X=complex(zeros(M,N));
D=complex(zeros(M,N));
for k=1:N
  X(:,k)=flipud(x(k:k+M-1));   %Snapshot matrix
  D(:,k)=flipud(d(k:k+M-1));   %Snapshot matrix
end

R=1/N*(X*X');       %Correlation Matrix

P=zeros(M,1);     %Cross Correlation Matrix
for i=M:K
  P=P+ flipud(x(i-M+1:i)).*conj(d(i));
end
P=P/N;
w_opt=R\P;      % Optimum Weiner filter

%% TDLMS
mu=0.5/M;
del=0.02;

[Q,D]=eig(R);
eig_vals=eig(R);
D=D+del*eye(size(D));
Z=Q'*X;            %Transformed input
w_optT=Q'*w_opt;    %Transformed Optimum weight vector

w_norm=zeros(M,1);
w_td=zeros(M,1);
for n=M:K
  e=d(n)-(w_norm'*flipud(x(n-M+1:n)));    %Error in normal domain
  z=Q'*flipud(x(n-M+1:n));
  w_td=w_td+mu*inv(D)*z*conj(e);
  w_norm=Q*w_td;          %Normal domain filter coeff
  squared_error(n)=(abs(e).^2);   %Error in normal domain
  dev(n)=(w_norm-w_opt)'*(w_norm-w_opt);   %Error in normal domain
end
figure(1);hold on;plot(20*log10(dev(M:K)))
title('Squared Deviation TDLMS');
xlabel('Number of Iterations');
ylabel('Squared deviation in dB')
figure(2);hold on; plot(20*log10(squared_error(M:K)));
title('Squared Error TDLMS')
xlabel('Number of Iterations');
ylabel('Squared Error in dB')