
% EECS 844 Exam-1
% Manjish Adhikari - 2870257
% Question 2-4

clear;
load('V:\EECS-844\P2.mat');

M=20;     %Length of each snapshot
K=length(x);
N=K-M+1;
X=complex(zeros(M,N));
for i=1:M
  for j=1:N
        col=i:K-M+i;
         X(20-i+1,j)=(x(col(j)));    %Snapshot matrix 
  end
end


%No diagonal loading Q. no 2

clear i j
for i=1:N
  temp=X(:,1:i);
  R=zeros(M,M);
  for id=1:size(temp,2)
      
      R=R+(temp(:,id)*temp(:,id)');     
  end
  R=R/size(temp,2);     %Correlation Matrix 
  eig_d0=eig(R);         %Eigen values

  max_eig_d0(i)=max(abs(eig_d0));     %Max eigen values
  min_eig_d0(i)=min(abs(eig_d0));       %Min eigen values
  cond_no_d0(i)=max_eig_d0(i)/min_eig_d0(i);    %Condition number
end
figure;plot(eig_d0);title('Eigen values with no diag loading Q4 for all snapshots'); xlabel('Number of snapshots'); ylabel('Eigen values in dB')
figure;plot(10*log10(abs(max_eig_d0)));title('Max eigen val with no diagonal loading Q2');xlabel('Number of snapshots'); ylabel('Max Eigen values in dB')
figure;plot(10*log10(abs(min_eig_d0)));title('Min eigen val with no diagonal loading Q2');xlabel('Number of snapshots'); ylabel('  Min Eigen values in dB')
figure;plot(cond_no_d0);title('Cond no.  with no diagonal loading Q2');xlabel('Number of snapshots'); ylabel('Condition number')



%Diagonal loading var =1
clear i j
clear temp;
var=1;
for i=1:N
  temp=X(:,1:i);
  R=zeros(M,M);
  for id=1:size(temp,2)
      
      R=R+(temp(:,id)*temp(:,id)');
  end
  R=R/size(temp,2)+var*eye(M);   %Correlation matrix with diagonal loading
  eig_d1=eig(R);
  max_eig_d1(i)=max(abs(eig_d1));
  min_eig_d1(i)=min(abs(eig_d1));
  cond_no_d1(i)=max_eig_d1(i)/min_eig_d1(i);
end
% figure;plot(eig_d1);title('Eigen values with  diag loading var =1 Q4');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
% figure;plot(10*log10(abs(max_eig_d1)));title('Max eig val with diag load var =1');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
% figure;plot(10*log10(abs(min_eig_d1)));title('Min eig val with diag load var =1');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
 figure;plot(cond_no_d1);title('Condition no with diag load var =1');xlabel('Number of snapshots'); ylabel('Condition number')


%Diagonal loading var =10
clear i j
var=10;
for i=1:N
  temp=X(:,1:i);
  R=zeros(M,M);
  for id=1:size(temp,2)
      
      R=R+(temp(:,id)*temp(:,id)');
  end
  R=R/size(temp,2)+var*eye(M);
  eig_d2=eig(R);
  max_eig_d2(i)=max(abs(eig_d2));
  min_eig_d2(i)=min(abs(eig_d2));
  cond_no_d2(i)=max_eig_d2(i)/min_eig_d2(i);
end
% figure;plot(eig_d2);title('Eigen values with  diag loading var =10');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
% figure;plot(10*log10(abs(max_eig_d2)));title('Max eig val with diag load var =10');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
% figure;plot(10*log10(abs(min_eig_d2)));title('Min eig val with diag load var =10');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
 figure;plot(cond_no_d2);title('Condition no with diag load var =10');xlabel('Number of snapshots'); ylabel('Condition number')
% 
%Diagonal loading var =100
clear i j
var=100;
for i=1:N
  temp=X(:,1:i);
  R=zeros(M,M);
  for id=1:size(temp,2)
      
      R=R+(temp(:,id)*temp(:,id)');
  end
  R=R/size(temp,2)+var*eye(M);
  eig_d3=eig(R);
  max_eig_d3(i)=max(abs(eig_d3));
  min_eig_d3(i)=min(abs(eig_d3));
  cond_no_d3(i)=max_eig_d3(i)/min_eig_d3(i);
end
% figure;plot(eig_d3);title('Eigen values with  diag loading var =100');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
% figure;plot(10*log10(abs(max_eig_d3)));title('Max eig val with diag load var =100');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
% figure;plot(10*log10(abs(min_eig_d3)));title('Min eig val with diag load var =100');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
 figure;plot(cond_no_d3);title('Condition no with diag load var =100');xlabel('Number of snapshots'); ylabel('Condition number')

 %Q. no. 3
figure(20);plot(10*log10(abs(max_eig_d0)),'r','DisplayName',' no diag load ');
figure(20);hold on;plot(10*log10(abs(max_eig_d1)),'b','DisplayName','diag loading var =1');
figure(20);hold on;plot(10*log10(abs(max_eig_d2)),'g','DisplayName','diag loading var =10');title('Max eigen value')
figure(20);hold on;plot(10*log10(abs(max_eig_d3)),'k','DisplayName',' diag loading var =100');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
legend('show')

figure(21);plot(10*log10(abs(min_eig_d0)),'r','DisplayName','no diag load ');
figure(21);hold on;plot(10*log10(abs(min_eig_d1)),'b','DisplayName',' diag load var =1');
figure(21);hold on;plot(10*log10(abs(min_eig_d2)),'g','DisplayName','diag load var =10');title('Min eigen value')
figure(21);hold on;plot(10*log10(abs(min_eig_d3)),'k','DisplayName',' diag load var =100');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
legend('show')

%Q. no 4 
figure(22);plot(eig_d0,'r','DisplayName','No diag loading');
figure(22);hold on; plot(eig_d1,'b','DisplayName',' Diag loading var=1');
figure(22);hold on; plot(eig_d2,'g','DisplayName',' Diag loading var=10');title('Eigen values with all snapshots')
figure(22);hold on; plot(eig_d3,'k','DisplayName',' Diag loading var=100');xlabel('Number of eigen values'); ylabel('Eigen values')
legend('show');

figure(23);plot(10*log10(eig_d0),'r','DisplayName','No diag loading');
figure(23);hold on; plot(10*log10(eig_d1),'b','DisplayName',' Diag loading var=1');
figure(23);hold on; plot(10*log10(eig_d2),'g','DisplayName',' Diag loading var=10');title('Eigen values')
figure(23);hold on; plot(10*log10(eig_d3),'k','DisplayName',' Diag loading var=100');xlabel('Number of snapshots'); ylabel('Eigen values in dB')
legend('show');


%figure(24);plot(cond_no_d0,'r','DisplayName','Cond no with no diag load ');
figure(24);hold on;plot(cond_no_d1,'b','DisplayName','Cond no with diag load var =1');
figure(24);hold on;plot(cond_no_d2,'g','DisplayName','Cond no with diag load var =10');title('Condition number')
figure(24);hold on;plot(cond_no_d3,'k','DisplayName','Cond no with diag load var =100');xlabel('Number of snapshots'); ylabel('Cond no')
legend('show')

