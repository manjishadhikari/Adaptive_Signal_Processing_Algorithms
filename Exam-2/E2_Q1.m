load('V:\EECS-844\Exam-3\P4.mat')
close all

for filter_length=0:40
  
  %% Making snapshot matrix from x
  M=filter_length;    
  K=length(x);
  N=K-M+1;            %Snapshots
  X=complex(zeros(M,N)); %Snapshot matrix
  D=complex(zeros(M,N)); 
  
  if M==0
    w=1;
  else
    for k=1:N
      X(:,k)=flipud(x(k:k+M-1));
      D(:,k)=flipud(x(k:k+M-1));
    end
  
    clear i j 
    R=1/N*X*ctranspose(X); %% Correlation matrix
   
    P=zeros(M,1);
    for i=M:size(X,2)
      P=P+ flipud(x(i-M+1:i)).*conj(d(i)); %Cross correlation matrix
    end
    P=P/N;
    
    w=inv(R)*P;     %Weiner filter
    
  end
  
  y=filter(conj(w),1,x); %Apply the filter with weiner filter coefficients to get estimate

  error=d-y;
  error_m1=ctranspose(error)*error/length(error); %Error from direct method
  error_m1_dB(filter_length+1)=20*log10(error_m1);
  if filter_length==0
    error_m2_dB(filter_length+1)=error_m1_dB(filter_length+1);
  else
    error_m2= std(d).^2-real(ctranspose(P)*inv(R)*P); %Error from analytical method
    error_m2_dB(filter_length+1)=20*log10(error_m2);
  end
end
figure;plot(error_m2_dB);title('Error performance surface');
xlabel('Filter Length');ylabel('Error in dB');
figure;plot(error_m1_dB);title('Error from Direct method')
xlabel('Filter Length');ylabel('Error in dB');
figure;plot(abs(w));title('Magnitude of filter coefficients');
xlabel('Filter Length');ylabel('Magnitude');
figure;plot(angle(w));title('Phase of filter coefficients');
xlabel('Filter Length');ylabel('Phase in radians');
figure;freqz(w,1000,'whole')
title('frequency response of filter for M=40')