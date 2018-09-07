load('V:\EECS-844\Exam-2\P4.mat')
  
  %% Making snapshot matrix from x
  M=40;     %Length of each snapshot
  K=length(x);
  N=K-M+1;
  num_samples=20*M;
  X=complex(zeros(M,N));
  theta=linspace(-pi,pi,num_samples);
  
  for k=1:N
    X(:,k)=flipud(x(k:k+M-1));  %Snapshot matrix
  end
  
  %% Correlation matrix
  
  clear i j
  R=1/N*X*ctranspose(X);
  R_inv=inv(R);
  
  p=zeros(1,num_samples);
  C=transpose(exp(-j*pi/2*[0:M-1]));  %Gain constranit for theta=2*pi*0.25
  null_ang1=-0.8*pi;
  null_ang2=-0.4*pi;
  null_ang=linspace(null_ang1,null_ang2,800);
  Rc=zeros(M,M);
  for i=1:length(null_ang)
    sv=transpose(exp(-j*null_ang(i)*[0:M-1]));
    Rc=Rc+sv*sv';
  end
  eig_val=eig(Rc);
  [~,idx]=sort(eig_val,'descend');
  [Q,~,~]=eig(Rc);
  Q_sorted=Q(:,idx);
  null_const=[0 10 12 14];
  for m=1:length(null_const)
    g=[1 zeros(1,null_const(m))]';
   if m==1
       Q_top=[];
   else
    Q_top=Q_sorted(:,1:null_const(m));
   end
    C_new=[C Q_top];
    
    [U,~,~]=svd(C_new);
    Ca=U(:,length(g)+1:end);
    
    v=(ctranspose(C_new)*C_new)\g;
    wa=((ctranspose(Ca)*R*Ca))\(ctranspose(Ca)*R*C_new*v);
    w_gsc=C_new*v-Ca*wa;
     gsc_spectrum=zeros(num_samples,1);
    for i=1:num_samples
      sv=transpose(exp(-j*theta(i)*[0:M-1]));
      gsc_spectrum(i)=ctranspose(sv)*w_gsc;
    end
    figure(m); plot(theta,20*log10(abs(gsc_spectrum)));title(sprintf('With null constraints N=%d',null_const(m)'));
  xlabel('Theta in radians'); ylabel('Power in dB')
  end

figure;plot([1:40],10*log10(abs(eig_val)));title('Eigen values')
xlabel('Number of null constraints'); ylabel('Eigen values in dB')