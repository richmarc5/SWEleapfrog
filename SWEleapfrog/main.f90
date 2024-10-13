
program SWE
    implicit none
    real ,allocatable :: h(:,:),x(:),y(:),u(:,:),v(:,:)
    real, allocatable ::unew(:,:),vnew(:,:),hnew(:,:),uold(:,:),vold(:,:),hold(:,:)
    real :: f,g,hval,dx,dy,dt,gam,Lw,x2,y2
    integer :: Nt,ii,jj,num,unit_number,it,i,j,n,kk,pp


    dx=0.025
    dy=dx
    num=2/dx+1
    Lw=7**2
    Nt=4000
     gam=0.5
     dt=0.001
     hval=1;
     g=1

    allocate(x(num))
    do ii=1,num
      x(ii)=-1+dx*(ii-1)
    end do
    y=x
    allocate(h(num,num))
    allocate(u(num,num))
    allocate(v(num,num))
    allocate(unew(num,num))
   allocate(vnew(num,num))
   allocate(hnew(num,num))
   allocate(uold(num,num))
   allocate(vold(num,num))
   allocate(hold(num,num))
    h=0

! Open the file for writing
unit_number = 10
open(unit=unit_number, file='matrix.txt', status='unknown')

    do ii=1,num
      x2 = dx*(ii - 1)
        do jj=1,num
         y2 = dy*(jj - 1)
         h(ii,jj)=2*exp(-Lw*((x(ii)))**2-Lw*((y(jj)))**2)

        end do
        write(unit_number,*) h(ii,:)
    end do





!     First integration step


i=1;
j=1;

do ii=1,num
  it=ii
  if (ii==1)then
     it=ii+1
 end if
  if(ii==num)then
  it=ii-1
end if
  do jj=1,num
    j=jj
         if(jj==1)then
         j=jj+1

        end if

        if(jj==num)then
          j=jj-1
        end if

         u(ii,jj)=u(ii,jj)+dt*(-g*(h(it+1,j)-h(it,j))+0.25*f*(v(it,j)+v(it+1,j)+v(it+1,j-1)+v(it,j-1)))


       v(ii,jj)=v(ii,jj)+dt*(-g*(h(it,j+1)-h(it,j))/dy- 0.25*f*(u(it,j)+u(it,j+1)+u(it-1,j+1)+u(it-1,j)))


         h(ii,jj)=h(ii,jj)-dt*hval*((u(it,j)-u(it-1,j))/dx +(v(it,j)-v(it,j-1))/dy)


  end do

end do


uold=u;
vold=v;
hold=h;



open(unit=unit_number, file='matrix2.txt', status='unknown')


   do n=1,Nt

        do ii=1,num
             it=ii;

                 if (ii==1) then
                     it=ii+1
                  end if
                  if(ii==num) then
                  it=ii-1
                  end if
              do jj=1,num
                j=jj
                    if(jj==1)then
                    j=jj+1

                     end if

                  if(jj==num) then
                    j=jj-1
                  end if

               unew(ii,jj)=uold(ii,jj)+2*dt*(-g*(h(it+1,j)-h(it,j))/dx+0.25*f*(v(it,j)+v(it+1,j)+v(it+1,j-1)+v(it,j-1)));
               vnew(ii,jj)=vold(ii,jj)+2*dt*(-g*(h(it,j+1)-h(it,j))/dy- 0.25*f*(u(it,j)+u(it,j+1)+u(it-1,j+1)+u(it-1,j)));
               hnew(ii,jj)=hold(ii,jj)-2*dt*hval*((u(it,j)-u(it-1,j))/dx +(v(it,j)-v(it,j-1))/dy);


          end do
      end do

      do ii=1,num
        do jj=1,num
          u(ii,jj)=u(ii,jj)+gam*(uold(ii,jj)-2*u(ii,jj)+unew(ii,jj));
          v(ii,jj)=v(ii,jj)+gam*(vold(ii,jj)-2*v(ii,jj)+vnew(ii,jj));
          h(ii,jj)=h(ii,jj)+gam*(hold(ii,jj)-2*h(ii,jj)+hnew(ii,jj));
        end do
      end do

      uold=u;
      u=unew;
      vold=v;
      v=vnew;
      hold=h;
      h=hnew;
!
!      h(:,1) = h(:,2);      u(:,1) = u(:,2);       v(:,1) = -v(:,2);
!       h(:,num) = h(:,num-1);  u(:,n) = u(:,num-1);   V(:,num) = -v(:,num-1);
!       h(1,:) = h(2,:);      u(1,:) = -u(2,:);      v(1,:) = v(2,:);
!       h(num,:) = h(num-1,:);  u(num,:) = -u(num-1,:);  v(num,:) = v(num-1,:);

        do kk=1,num
         do pp=1,num

        h(kk,1) = h(kk,2)
          u(kk,1) = u(kk,2)
          v(ii,1) = -v(kk,2);
          h(kk,num) = h(kk,num-1)
          u(kk,num) = u(kk,num-1)
            v(kk,num) = -v(kk,num-1)
           h(1,pp) = h(2,pp)
            u(1,pp) = -u(2,pp)
              v(1,pp) = v(2,pp);
           h(num,pp) = h(num-1,pp)
             u(num,pp) = -u(num-1,pp)
               v(num,pp) = v(num-1,pp)

            end do
         end do


if(mod(n,10)==0)then

 do ii=1,num
  write(unit_number,*) h(ii,:)
  end do
!print*,'h= ',h
end if

    print *, 'time = ',Nt-n



    enddo






end program

