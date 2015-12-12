%clear all
%{ 
if ispc
    mex('-v','-largeArrayDims','eig_omp.cpp','-lmwlapack','COMPFLAGS=$COMPFLAGS /openmp','LINKFLAGS=$LINKFLAGS /nodefaultlib:vcomp "$MATLABROOT\bin\win64\libiomp5md.lib"')
else
    mex('-g','-v','-largeArrayDims','eig_omp.cpp','-lmwlapack','COMPFLAGS=$COMPFLAGS -fopenmp -pthread')
end
%}

% Parameters
%n = 50;  % Matrix size == n x n
%nb = 500;  % Number of matrices nb

% Test matrices (real nonsymmetric, complex nonhermitian, real symmetric, complex hermitian)
dge = rand(n);
zge = rand(n)+1i*rand(n);
dsy = rand(n); dsy = triu(dsy)+triu(dsy,1)';
zhe = rand(n)+1i*rand(n); zhe = triu(zhe,1)+triu(zhe,1)'+diag(real(diag(zhe)));

% Functionality tests - should not give any errors, just want to check that correct number of outputs given.
%{
[V,E]=eig_omp(zhe);
D=eig_omp(zhe);
eig_omp(zhe);
[V,E]=eig_omp(dsy);
D=eig_omp(dsy);
eig_omp(dsy);
[V,E]=eig_omp(zge);
D=eig_omp(zge);
eig_omp(zge);
[V,E]=eig_omp(dge);
D=eig_omp(dge);
eig_omp(dge);
%}

% Numerical test vs built-in eig()
% {
%disp(sprintf('\nTest vs builtin:  \teig_omp Eigenvalues\teig Eigenvalues  \tEigenvectors differences'));
%[V1,E1]=eig_omp(dsy); [V2,E2]=eig(dsy); disp(sprintf('Real symmetric   \t%10.5f\t\t%10.5f\t\t%10.5f\n',[diag(E1),diag(E2),sum(abs(V1)-abs(V2),2)]'));
%[V1,E1]=eig_omp(zhe); [V2,E2]=eig(zhe); disp(sprintf('Complex hermitian\t%10.5f\t\t%10.5f\t\t%10.5f\n',[diag(E1),diag(E2),sum(abs(V1)-abs(V2),2)]'));
%[V1,E1]=eig_omp(dge); [V2,E2]=eig(dge); disp(sprintf('Real general     \t%10.5f +i%10.5f\t%10.5f +i%10.5f\t%10.5f\n',[real(diag(E1)),imag(diag(E1)),real(diag(E2)),imag(diag(E2)),sum(abs(V1)-abs(V2),2)]'));
%[V1,E1]=eig_omp(zge); [V2,E2]=eig(zge); disp(sprintf('Complex general  \t%10.5f +i%10.5f\t%10.5f +i%10.5f\t%10.5f\n',[real(diag(E1)),imag(diag(E1)),real(diag(E2)),imag(diag(E2)),sum(abs(V1)-abs(V2),2)]'));
disp(sprintf('\nTest vs builtin:  \tEigenvalues diffs  \tEigenvectors differences'));
disp(sprintf('------------------------------------------------------------------------'));
[V1,E1]=eig_omp(dsy); [V2,E2]=eig(dsy); disp(sprintf('Real symmetric   \t% 10.5g\t\t% 10.5g',abs(sum(diag(E1)-diag(E2))),sum(sum(abs(V1)-abs(V2)))));
[V1,E1]=eig_omp(zhe); [V2,E2]=eig(zhe); disp(sprintf('Complex hermitian\t% 10.5g\t\t% 10.5g',abs(sum(diag(E1)-diag(E2))),sum(sum(abs(V1)-abs(V2)))));
[V1,E1]=eig_omp(dge); [V2,E2]=eig(dge); disp(sprintf('Real general     \t% 10.5g\t\t% 10.5g',abs(sum(diag(E1)-diag(E2))),sum(sum(abs(V1)-abs(V2)))));
[V1,E1]=eig_omp(zge); [V2,E2]=eig(zge); disp(sprintf('Complex general  \t% 10.5g\t\t% 10.5g',abs(sum(diag(E1)-diag(E2))),sum(sum(abs(V1)-abs(V2)))));
disp(sprintf('\nConsistency tests:\tV*E*inv(V)-A\tV*E*V''-A\tV*V''-I\t\tV''*V-I'))
disp(sprintf('------------------------------------------------------------------------------'));
[V,E]=eig_omp(dsy); disp(sprintf('Real symmetric   \t% 10.4g\t% 10.4g\t% 10.4g\t% 10.4g',sum(sum(V*E*inv(V)-dsy)),sum(sum(V*E*V'-dsy)),sum(sum(V*V'-eye(n))),sum(sum(V'*V-eye(n)))));
[V,E]=eig_omp(zhe); disp(sprintf('Complex hermitian\t% 10.4g\t% 10.4g\t% 10.4g\t% 10.4g',sum(sum(V*E*inv(V)-zhe)),sum(sum(V*E*V'-zhe)),sum(sum(V*V'-eye(n))),sum(sum(V'*V-eye(n)))));
[V,E]=eig_omp(dge); disp(sprintf('Real general     \t% 10.4g\t% 10.4g\t% 10.4g\t% 10.4g',sum(sum(V*E*inv(V)-dge)),sum(sum(V*E*V'-dge)),sum(sum(V*V'-eye(n))),sum(sum(V'*V-eye(n)))));
[V,E]=eig_omp(zge); disp(sprintf('Complex general  \t% 10.4g\t% 10.4g\t% 10.4g\t% 10.4g',sum(sum(V*E*inv(V)-zge)),sum(sum(V*E*V'-zge)),sum(sum(V*V'-eye(n))),sum(sum(V'*V-eye(n)))));
%}

matfn{1} = @(a,b) triu(a)+triu(a,1)';
matfn{2} = @(a,b) triu(a,1)+1i*triu(b,1)+triu(a,1)'-1i*triu(b,1)'+diag(diag(a));
matfn{3} = @(a,b) a;
matfn{4} = @(a,b) a+1i*b;

% Numerical test of sorting functionality
% {
clear V1; clear V2; clear E1; clear E2; clear B;
nt = 10; if(nt>nb); nt=nb; end
for tp=1:4
  for ii=1:nt; B(:,:,ii) = matfn{tp}(rand(n),rand(n)); end
  for ii=1:nt; 
    [vv,D]=eig(B(:,:,ii)); [~,isr]=sort(real(diag(D)),'descend'); 
    V1(:,:,ii)=vv(:,isr); E1(:,ii)=diag(D); E1(:,ii)=E1(isr,ii);
  end;
  [V2,E2]=eig_omp(B,'sort',-1);
  sortres(:,tp) = [sum(sum(sum(abs(V1)-abs(V2)))) sum(sum(abs(E1-E2)))];
end
disp(sprintf('\nSort Test:  \t\tEigenvecDiffs   EigenvalDiffs'));
disp(sprintf('-----------------------------------------------------'));
disp(sprintf('Real symmetric   \t% 12.5g\t% 12.5g',sortres(:,1)));
disp(sprintf('Complex hermitian\t% 12.5g\t% 12.5g',sortres(:,2)));
disp(sprintf('Real general     \t% 12.5g\t% 12.5g',sortres(:,3)));
disp(sprintf('Complex general  \t% 12.5g\t% 12.5g',sortres(:,4)));
%}

% Numerical test of orthogonalisation functionality
% {
clear V1; clear V2; clear E1; clear E2; clear B;
for tp=1:4
  for ii=1:nb; B(:,:,ii) = matfn{tp}(rand(n),rand(n)); end
 %for ii=1:nb;
 %  [vv,D]=eigorth(B(:,:,ii),eps);  [~,isr]=sort(real(D)); 
 %  V1(:,:,ii)=vv(:,isr); E1(:,ii) = D(isr);
 %end;
  tic; [V1,E1]=eigorth(B,eps); t1(tp)=toc;
  tic; [V2,E2]=eig_omp(B,'orth','sort'); t2(tp)=toc;
  for ii=1:nb; [~,isr]=sort(real(E1(:,ii))); V1(:,:,ii)=V1(:,isr,ii); E1(:,ii) = E1(isr,ii); end
  orthres(:,tp) = [sum(sum(sum(abs(V1)-abs(V2)))) sum(sum(abs(E1-E2)))];
end
disp(sprintf('\nOrth Test:  \t\tEigenvecDiffs   EigenvalDiffs   Time(eigorth)   Time(eig_omp)         Speedup'));
disp(sprintf('-----------------------------------------------------------------------------------------------------'));
disp(sprintf('Real symmetric   \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g',orthres(:,1),t1(1),t2(1),t1(1)/t2(1)));
disp(sprintf('Complex hermitian\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g',orthres(:,2),t1(2),t2(2),t1(2)/t2(2)));
disp(sprintf('Real general     \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g',orthres(:,3),t1(3),t2(3),t1(3)/t2(3)));
disp(sprintf('Complex general  \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g',orthres(:,4),t1(4),t2(4),t1(4)/t2(4)));
%}

% Speed tests
% {
clear E1; clear E2; clear B; clear ch; clear ch2;
for tp=1:4
  for ii=1:nb; B(:,:,ii) = matfn{tp}(rand(n),rand(n)); end
  tic; [V,E]=eig_omp(B); t1(tp)=toc; 
  for ii=1:nb; 
    ch(ii)=sum(sum(V(:,:,ii)*diag(E(:,ii))*V(:,:,ii)'-B(:,:,ii))); ch2(ii)=sum(sum(V(:,:,ii)'*V(:,:,ii)-eye(size(V(:,:,ii))))); 
  end; check1(tp,:)=[sum(ch) sum(ch2)];
  tic; for ii=1:size(E,2); [V(:,:,ii),ee]=eig(B(:,:,ii)); E2(:,ii)=diag(ee); end; t2(tp)=toc; 
  for ii=1:nb; 
    ch(ii)=sum(sum(V(:,:,ii)*diag(E(:,ii))*V(:,:,ii)'-B(:,:,ii))); ch2(ii)=sum(sum(V(:,:,ii)'*V(:,:,ii)-eye(size(V(:,:,ii))))); 
  end; check2(tp,:)=[sum(ch) sum(ch2)];
  check3(tp)=sum(sum(E-E2));
end
disp(sprintf('\nSpeed Test:  \t\tTime(eig_omp)\t   Time(eig)\t   Speedup    Eigenvalue diffs   Sum(V*E*V''-A)\tSum(V''*V-I)'));
disp(sprintf('-------------------------------------------------------------------------------------------------------------------'));
disp(sprintf('Real symmetric   \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g',t1(1),t2(1),t2(1)/t1(1),abs(check3(1)),abs(check1(1,1)),abs(check1(1,2))));
disp(sprintf('Complex hermitian\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g',t1(2),t2(2),t2(2)/t1(2),abs(check3(2)),abs(check1(2,1)),abs(check1(2,2))));
disp(sprintf('Real general     \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g',t1(3),t2(3),t2(3)/t1(3),abs(check3(3)),abs(check1(3,1)),abs(check1(3,2))));
disp(sprintf('Complex general  \t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g\t% 12.5g',t1(4),t2(4),t2(4)/t1(4),abs(check3(4)),abs(check1(4,1)),abs(check1(4,2))));
%}

