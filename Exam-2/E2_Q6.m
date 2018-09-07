load('V:\EECS-844\Exam-2\P6.mat')

close all;
[h1,w]=freqz(x1);
figure(1);plot(lp(h1));
[h2,w]= freqz(x2);
figure(1);hold on ; plot(lp(h2));

[h3]=freqz(x3);
figure(1);hold on ;plot(lp(h3));
[h4]=freqz(x4);
hold on ; plot(lp(h4))


hold on ;plot(lp(h5))
legend('1','2','3','4');

figure(2);plot(abs(x1)); hold on ;plot(abs(x2));
hold on ;plot(abs(x3)); hold on ;plot(abs(x4));
%hold on;plot(abs(x5));
legend('x1','x2','x3','x4');
xlabel('Index');ylabel('Magnitude');
title('All four signals')

%Change ip op to find error surface and convergence
x=x1;
d=x2;

for filter_length=0:30
  
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
    
    
    P=zeros(M,1);
    for i=M:size(X,2)
      P=P+ flipud(x(i-M+1:i)).*conj(d(i));
    end
    P=P/N;
    
    %w0=inv(R)*P;
    w=inv(R)*P;
    
  end
  
  y=filter(conj(w),1,x);

  error=d-y;
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
%figure;plot([0:filter_length],error_m2_dB);title('Error performance surface')
figure;plot([0:filter_length],error_m1_dB);title('Direct method')
xlabel('Filter length');ylabel('Error in bB')
figure;plot(20*log10(abs(w)));title('Magnitude of filter coefficients');
xlabel('Index');ylabel('Power')
%figure;plot(angle(w));title('Phase of filter coefficients');
%figure;freqz(w)
