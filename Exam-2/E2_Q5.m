load('V:\EECS-844\Exam-2\P4.mat')

filter_length=40;
for idx=1:length(filter_length)
  
  %% Making snapshot matrix from x
  M=filter_length(idx);     %Length of each snapshot
  K=length(x);
  N=K-M+1;
  num_samples=20*M;
  X=complex(zeros(M,N));
  theta=linspace(-pi,pi,num_samples);
  
  for k=1:N
    X(:,k)=flipud(x(k:k+M-1));
  end
  
  %% Correlation matrix
  
  clear i j
  R=1/N*X*ctranspose(X);
  R_inv=inv(R);
  
  p=zeros(1,num_samples);
  C=transpose(exp(-j*pi/2*[0:M-1]));  %Gain constraint for theta=2*pi*0.25 or steering vector
  [U,~,~]=svd(C);
  Ca=U(:,2:end);
  g=1;
  v=(ctranspose(C)*C)\g;
  wa=((ctranspose(Ca)*R*Ca))\(ctranspose(Ca)*R*C*v);
  w_gsc=C*v-Ca*wa;
 
  for i=1:num_samples
    sv=transpose(exp(-j*theta(i)*[0:M-1]));
    gsc_spectrum(i)=ctranspose(sv)*w_gsc;
  end
  
  figure(1);plot(theta*180/pi,20*log10(abs(gsc_spectrum)));
  w_mvdr=R_inv*C/(ctranspose(C)*R_inv*C);    %C is also steering vector
  %C is sv at theta =pi/2
  mvdr_spectrum=zeros(num_samples,1);
  for i=1:num_samples
    sv=transpose(exp(-j*theta(i)*[0:M-1]));
    mvdr_spectrum(i)=ctranspose(sv)*w_mvdr;
  end
  
  figure(1); hold on; plot(theta*180/pi,20*log10(abs(mvdr_spectrum)));
  legend('GSC','MVDR')
  xlabel('Theta in angles');ylabel('Power in dB')
  
  
  
  
  
  
end
