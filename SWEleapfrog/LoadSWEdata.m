clear all
close all

load matrix.txt


figure(1);
surf(matrix)



load matrix2.txt

[rows,cols]=size(matrix2)

 dx=0.025
    dy=dx
    num=2/dx+1

  Nt=floor(rows/num)


for ii=1:1:Nt

  pts=[1:num]+(ii-1)*num;

  figure(111);
  surf(matrix2(pts,1:num))
  axis([1,num,1,num,-1,1.5])
  title(['time=',num2str(ii)])
  pause(0.001)

  end
