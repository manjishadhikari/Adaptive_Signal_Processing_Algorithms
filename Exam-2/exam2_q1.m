load('V:\EECS-844\Exam-2\P1.mat')


for filter_length=1:40

%% Making snapshot matrix from x
M=filter_length;     %Length of each snapshot
K=length(x);
N=K-M+1;
X=complex(zeros(M,N));
D=complex(zeros(M,N));

% for i=1:M
%   for j=1:N
%         col=i:K-M+i;
%          X(M-i+1,j)=(x(col(j)));    %Input Snapshot matrix 
%          D(M-i+1,j)=(d(col(j)));    %Output snapsot matrix
%   end
% end
for k=1:N
    X(:,k)=flipud(x(k:k+M-1));
    D(:,k)=flipud(x(k:k+M-1));
end

%% Correlation matrix

clear i j
R=zeros(M,M);
P=zeros(M,M);
P1=zeros(M,1);

    
  for id=1:size(X,2)
      
      R=R+(X(:,id)*X(:,id)');     
  end
  R=R/size(X,2);     %Correlation Matrix 
% Rn=1/size(X,2)*X*ctranspose(X);  
%   for id=1:size(X,2)
%       
%       P=P+(X(:,id)*conj(D(:,id)'));     
%   end
%   P=P/size(X,2);     %Cross Correlation Matrix
%   P=mean(P,2);
  
  for i=M:size(X,2)
    P1=P1+ flipud(x(i-M+1:i)).*conj(d(i));
  end
  P1=P1/size(X,2);
  
  
  %w0=inv(R)*P;
  w1=inv(R)*P1;
%   d_bar=zeros(1,size(X,2));
%   for i=M:size(X,2)
%     d_bar(i)=(conj(w1))'*x(i-M+1:i);
%   end
%   error_m1=0;
%   for i=M:size(X,2)
%     error_m1=error_m1+(d(i)-d_bar(i)).^2;
%   end
  y=filter(conj(w1),1,x);
  error=d-y(1:length(d));
  %error_m1=error_m1/size(X,2);
  error_m1=std(error);
  error_m1_dB(filter_length)=20*log10(error_m1);
 % error_m2=std(d).^2-(conj(w1)')*P1-(conj(P1)')*w1+(conj(w1)')*R*w1;
 error_m2= std(d).^2+real((conj(P1)')*inv(R)*P1);
 error_m2_dB(filter_length)=20*log10(error_m2);
  
end
figure;plot(error_m2_dB);title('Error performance surface')
figure;plot(error_m1_dB);title('Direct method')