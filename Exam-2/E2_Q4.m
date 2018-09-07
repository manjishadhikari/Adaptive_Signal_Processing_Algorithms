load('V:\EECS-844\Exam-2\P4.mat')

filter_length=[10 20 40 80];
for idx=1:length(filter_length)
  
  %% Making snapshot matrix from x
  M=filter_length(idx);     %Length of each snapshot
  K=length(x);
  N=K-M+1;
  num_samples=50*M;
  X=complex(zeros(M,N));
  theta=linspace(-pi,pi,50*M);
  
  for k=1:N
    X(:,k)=flipud(x(k:k+M-1));
  end
  
  %% Correlation matrix
  
  clear i j
  R=1/N*X*ctranspose(X);
  R_inv=inv(R);
  mvdr_power=zeros(1,num_samples);
  for i=1:num_samples
    sv=transpose(exp(-j*theta(i)*[0:M-1]));
    mvdr_power(i)=1/(ctranspose(sv)*R_inv*sv);
  end
  figure(idx);plot(theta*180/pi,20*log10(abs(mvdr_power)));
  title(sprintf('MVDR power Spectrum, M=%d',M)); 
  xlabel('Theta in angles'); ylabel('Power in dB');
  
end
