load('V:\EECS-844\Exam-2\P1.mat')
a=x;
x=d;
d=a;

for filter_length=0:40
  
  %% Making snapshot matrix from x
  M=filter_length;     %Length of each snapshot
  K=length(x);
  N=K-M+1;
  X=complex(zeros(M,N));
  D=complex(zeros(M,N));
  
  if M==0
    w=1;
  else
    for k=1:N
      X(:,k)=flipud(x(k:k+M-1));
      D(:,k)=flipud(x(k:k+M-1));
    end
    
    %% Correlation matrix
    
    clear i j
    
    
    R=1/N*X*ctranspose(X);
    
    %Cross correlation matrix
    P=zeros(M,1);
    for i=M:size(X,2)
      P=P+ flipud(x(i-M+1:i)).*conj(d(i));
    end
    P=P/N;
    
    %w0=inv(R)*P;
    w=inv(R)*P;
    
  end
  
  y=filter(conj(w),1,x);
  error=d-y(1:length(d));
  %error_m1=error_m1/size(X,2);
  %error_m1=std(error);
  error_m1=ctranspose(error)*error/length(error);
  error_m1_dB(filter_length+1)=20*log10(error_m1);
  % error_m2=std(d).^2-(conj(w1)')*P1-(conj(P1)')*w1+(conj(w1)')*R*w1;
  if filter_length==0
    error_m2_dB(filter_length+1)=error_m1_dB(filter_length+1);
  else
    error_m2= std(d).^2-real(ctranspose(P)*inv(R)*P);
    error_m2_dB(filter_length+1)=20*log10(error_m2);
  end
end
figure;plot(error_m2_dB);title('Error performance surface')
figure;plot(error_m1_dB);title('Direct method')
figure;plot(20*log10(abs(w)));title('Magnitude of filter coefficients');
figure;plot(angle(w));title('Phase of filter coefficients');
figure;freqz(w)
[H,a]=freqz(w);