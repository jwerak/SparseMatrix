function []=make()
  chdir('/home/veva/Dropbox/SKOLA/numerika_programovani/c++/matice/matrix_solver/scilab');
  A = make_10_10();
  //  make_10_10();
  B = make25_25();
  
//  a =inv(B);
  C = A*A;
  //  disp(C);
  Path = "b_10.dat";
  [b,x]=make_b(A,Path);
//  disp(inv(A)*b);
  
  Path = "b_25.dat";
  [b,x]=make_b(B,Path);
  disp(b);
endfunction


//=============================================================
//          VYTVOR CAST b
//=============================================================
function [b,x]= make_b(A,Path)
  n = size(A);
  x=(1:1:n(1,1))';
  b=A*x
  [fd,err]=sl_openf(Path)
  
  fprintf (fd,"%i",length(b));
  for i=1:length(b)
    fprintf(fd,"%f",b(i));
  end
endfunction
 
  
  
//VRACI JENOM MATICI LDU
function [LDU]= make_10_10() 
  filename = "test_mtx_10x10.dat";
  [fd,err]=sl_openf(filename);
  LDU = [100,-2,-3,-4,-5,-6,-7,-8,-9,-10;...
  -1,200,-3,-4,-5,-6,-7,-8,-9,-10;...
  -1,-2,300,-4,-5,-6,-7,-8,-9,-10;...
  -1,-2,-3,400,-5,-6,-7,-8,-9,-10;...
  -1,-2,-3,-4,500,-6,-7,-8,-9,-10;...
  -1,-2,-3,-4,-5,600,-7,-8,-9,-10;...
  -1,-2,-3,-4,-5,-6,700,-8,-9,-10;...
  -1,-2,-3,-4,-5,-6,-7,800,-9,-10;...
  -1,-2,-3,-4,-5,-6,-7,-8,900,-10;...
  -1,-2,-3,-4,-5,-6,-7,-8,-9,1000]
//  [m,n] = size(LDU);
//  
////  disp ("Matice LDU jak by mela vypadat po rozlozeni:")
////  disp(LDU);
//  
//  
//  L = zeros(m,n);
//  D = zeros(m,n);
//  U = zeros(m,n);
//  I = eye(m,n);
//  
//  //prvni radek
//  
//  D(1,1) = LDU(1,1);
//  U(1,2:n) = LDU(1,2:n);
//  
//  //  ostatni radky
//  for i=2:m
//    L(i,1:i-1) = LDU(1,1:i-1);
//    D(i,i) = LDU(i,i);
//    U(i,i+1:n) = LDU(i,i+1:n);
//  end
////  pause;
//  L = L+I;
//  U = U+I;
//  
////  disp(L);
////  disp(D);
////  disp(U);
//  
//  A = L*D*U;
////  disp ("Matice A:")
////  disp(A);
//  print_matrix(A,fd);
endfunction


function [A]= make25_25()
  filename = "test_mtx_25x25.dat";
  [fd,err]=sl_openf(filename);
  
  M = [4,-1,0,0,0;...
  -1,5,-1,0,0;...
  0,-1,6,-1,0;...
  0,0,-1,7,-1;...
  0,0,0,-1,8];
  [m,n] = size(M);
  
  I = -eye(m,n);  
  Z = zeros(m,n);
  A = [M,I,Z,Z,Z;...
  I,M,I,Z,Z;...
  Z,I,M,I,Z;...
  Z,Z,I,M,I;...
  Z,Z,Z,I,M];
//  A = [M,I,Z;...
//  I,M,I;...
//  Z,I,M];
  print_matrix(A,fd)
endfunction



function []= print_matrix(A,fd)
  [r,c] = size(A);
  nonzero = find(A);
  
  fprintf(fd,"%i %i %i ",max(size(nonzero)),r,c);
  for i=1:r
    for j=1:c
      if A(i,j) ~= 0
        fprintf(fd,'%i %i %d ', i,j, A(i,j)); 
      end
//      fprintf(fd,'\n');
    end
  end
endfunction


//VYTVOR SOUBOR
function [fd,err]=sl_openf(Path)
  file('close',file());//close all prev opened files
  [fd,err]=file('open', Path, 'unknown');
endfunction


function spy(A)
    [i,j] = find(A~=0)
    [N,M] = size(A)
    xsetech([0,0,1,1],[1,0,M+1,N])
    xrects([j;N-i+1;ones(i);ones(i)],ones(i));
    xrect(1,N,M,N);
endfunction

make;