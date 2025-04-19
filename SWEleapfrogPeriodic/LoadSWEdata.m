clear all
close all

##load matrix.txt
##
##
##figure(1);
##surf(matrix)



data=load('matrix2.txt');

[rows,cols]=size(data)
##L=5e6;
## dx=0.025
##    dy=dx
##    num=L/dx+1

     num=400
    L=5e6
    dx=L/(num-1)


  Nt=floor(rows/num)


for ii=1:1:Nt-5

  pts=[1:num]+(ii-1)*num;

  figure(222);
 imagesc(data(pts,1:num))

##  axis([1,num,1,num,-1,2.5])
  title(['time=',num2str(ii)])


  figure(111);
  surf(data(pts,1:num))

##  axis([1,num,1,num,-1,2.5])
  title(['time=',num2str(ii)])
  pause(0.05)

end


##figure(222);
##plot(data(:,num+1))
##title('Volume')
##
##figure(223);
##plot(data(:,num+2))
##title('Potential Energy')
##
##figure(224);
##plot(data(:,num+3))
##title('Kinetic Energy')
##
##
##
##

##figure(225);
##plot(data(:,num+3)+data(:,num+2))
##title('Total Energy')
##ii=50
## pts=[1:num]+(ii-1)*num;
##
##  figure(111);
##  surf(data(pts,1:num))
##  axis([1,num,1,num,-1,2.5])
##  title(['time=',num2str(ii)])
##  pause(0.01)






