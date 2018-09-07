load('V:\EECS-844\Exam-3\P4.mat')
M=length(x);
h_nmf=conj(flipud(x))./(x'*x);  %Normalized matched filter

K=2*M;
impulse_pos=1.5*M;

A=toeplitz([transpose(x) zeros(1,K-1)],[x(1) zeros(1,K-1)]);
%figure(1);imagesc(lp(A));

A2=A;
row_zero_idx=[impulse_pos-2 impulse_pos-1 impulse_pos+1 impulse_pos+2];
A2(row_zero_idx,:)=0;   % A with some rowsrows set to 0
%figure(2);imagesc(lp(A2));


em=zeros(size(A,1),1);%Elementary vector
em(impulse_pos)=1;

eig_vals=eig(A'*A);
eig_max=max(eig_vals);

h_mmf1=(inv(A'*A)*A'*em);         %h_mmf with no diag loading
h_mmf2=(inv(A2'*A2)*A2'*em);        %h_mmf with no diag loading and some rows zero

h_mmf3=(inv(A'*A+0.01*eig_max*eye(size(A'*A,1)))*A'*em);   %h_mmf with  diag loading
h_mmf4=(inv(A2'*A2+0.01*eig_max*eye(size(A'*A,1)))*A2'*em); %h_mmf with  diag loading and some rows zero

h_nmmf1=h_mmf1/(sqrt(h_mmf1'*h_mmf1)*sqrt(x'*x));
h_nmmf2=h_mmf2/(sqrt(h_mmf2'*h_mmf2)*sqrt(x'*x));
h_nmmf3=h_mmf3/(sqrt(h_mmf3'*h_mmf3)*sqrt(x'*x));
h_nmmf4=h_mmf4/(sqrt(h_mmf4'*h_mmf4)*sqrt(x'*x));

conv_with_mf=conv(h_nmf,x);  %Convolution with norm matched filter
conv_with_mmf1=conv(h_nmmf1,x);  %Convolution with norm mismatched filter1
conv_with_mmf2=conv(h_nmmf2,x);  %Convolution with norm mismatched filter2
conv_with_mmf3=conv(h_nmmf3,x);  %Convolution with norm mismatched filter3
conv_with_mmf4=conv(h_nmmf4,x);  %Convolution with norm mismatched filter4

mism_loss1=-20*log10(max(abs(conv_with_mmf1))/max(abs(conv_with_mf)));  %Mismatch loss for mmf1
mism_loss2=-20*log10(max(abs(conv_with_mmf2))/max(abs(conv_with_mf)));  %Mismatch loss for mmf1
mism_loss3=-20*log10(max(abs(conv_with_mmf3))/max(abs(conv_with_mf)));  %Mismatch loss for mmf1
mism_loss4=-20*log10(max(abs(conv_with_mmf4))/max(abs(conv_with_mf)));  %Mismatch loss for mmf1

cmf=transpose(conv_with_mf);
conv_mf=horzcat(zeros(1,50),cmf,zeros(1,50));
l=size(conv_mf,2)+1;
l2=length(conv_with_mmf1)+1;
figure(3)
subplot(5,1,1);plot(-(l/2-1):(l/2-1),20*log10(abs(conv_mf)));title('Conv with matched filter')
xlim([-149 149]);ylim([min(20*log10(abs(conv_with_mf))) max(20*log10(abs(conv_with_mf)))])
xlabel('Delay ');ylabel('Filter output (dB)');
subplot(5,1,2);plot((-(l2/2-1):l2/2-1),20*log10(abs(conv_with_mmf1)));title('Conv with mismatched filter1- Normal ')
xlabel('Delay ');ylabel('Filter Output (dB)');ylim([min(20*log10(abs(conv_with_mmf1))) max(20*log10(abs(conv_with_mmf1)))])
subplot(5,1,3);plot((-(l2/2-1):l2/2-1),20*log10(abs(conv_with_mmf2)));title('Conv with mismatched filter2- Row zeros')
xlabel('Delay ');ylabel('Filter Output (dB)');ylim([min(20*log10(abs(conv_with_mmf2))) max(20*log10(abs(conv_with_mmf2)))])
subplot(5,1,4);plot((-(l2/2-1):l2/2-1),20*log10(abs(conv_with_mmf3)));title('Conv with mismatched filter3- Diagonal loading')
xlabel('Delay ');ylabel('Filter Output (dB)');ylim([min(20*log10(abs(conv_with_mmf3))) max(20*log10(abs(conv_with_mmf3)))])
subplot(5,1,5);plot((-(l2/2-1):l2/2-1),20*log10(abs(conv_with_mmf4)));title('Conv with mismatched filter4-Row zero and diag  loading')
ylim([min(20*log10(abs(conv_with_mmf4))) max(20*log10(abs(conv_with_mmf4)))])
xlabel('Delay ');ylabel('Filter Output (dB)')

disp('============================')
disp(sprintf('mismatchedfilter1- Normal  ==> %2.4f',(mism_loss1)))
disp(sprintf('mismatchedfilter2- Rows Zero==>  %2.4f',(mism_loss2)))
disp(sprintf('mismatchedfilter3- Diagonal Loading==> %2.4f',(mism_loss3)))
disp(sprintf('mismatchedfilter4-Row zero and diagonal loading==> %2.4f',(mism_loss4)))
disp('============================')

%Deconvolution
unknwn_sys_mf=conv(y,h_nmf);

unknwn_sys_mmf1=conv(y,h_nmmf1);
unknwn_sys_mmf2=conv(y,h_nmmf2);
unknwn_sys_mmf3=conv(y,h_nmmf3);
unknwn_sys_mmf4=conv(y,h_nmmf4);

cmf=transpose(unknwn_sys_mf);
conv_mf=horzcat(zeros(1,50),cmf,zeros(1,50));
l=length(unknwn_sys_mf)+1;
l2=length(unknwn_sys_mmf1)+1;

figure(4)
subplot(5,1,1);plot(-(l/2-1):(l/2-1),20*log10(abs(unknwn_sys_mf)));title('Deconv with matched filter')
xlim([-250 250]);ylim([min(20*log10(abs(unknwn_sys_mf))) max(20*log10(abs(unknwn_sys_mf)))])
xlabel('Delay ');ylabel('Deconvoluted Output (dB)');
subplot(5,1,2);plot((-(l2/2-1):l2/2-1),20*log10(abs(unknwn_sys_mmf1)));title('Deconv with mismatched filter1- Normal ')
xlabel('Delay ');ylabel('Deconvoluted Output (dB)');%ylim([min(20*log10(abs(unknwn_sys_mmf1))) max(20*log10(abs(conv_with_mmf1)))])

subplot(5,1,3);plot((-(l2/2-1):l2/2-1),20*log10(abs(unknwn_sys_mmf2)));title('Deconv with mismatched filter2- Row zeros')
xlabel('Delay ');ylabel('Deconvoluted Output (dB)');%ylim([min(20*log10(abs(unknwn_sys_mmf2))) max(20*log10(abs(conv_with_mmf1)))])

subplot(5,1,4);plot((-(l2/2-1):l2/2-1),20*log10(abs(unknwn_sys_mmf3)));title('Deconv with mismatched filter3- Diagonal loading')
xlabel('Delay ');ylabel('Deconvoluted Output (dB)');%ylim([min(20*log10(abs(unknwn_sys_mmf3))) max(20*log10(abs(conv_with_mmf1)))])

subplot(5,1,5);plot((-(l2/2-1):l2/2-1),20*log10(abs(unknwn_sys_mmf4)));title('Deconnv with mismatched filter4-Row zero and diag  loading')
xlabel('Delay ');ylabel('Deconvoluted Output (dB)');%ylim([min(20*log10(abs(unknwn_sys_mmf4))) max(20*log10(abs(conv_with_mmf1)))])






