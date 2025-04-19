
program SWE
    implicit none
    real ,allocatable :: h(:,:),x(:),y(:),u(:,:),v(:,:)
    real, allocatable ::unew(:,:),vnew(:,:),hnew(:,:),uold(:,:),vold(:,:),hold(:,:)
    real :: f,g,hval,dx,dy,dt,gam,Lw,x2,y2,Stabcriteria,L
    real :: Volume,Epot,Ekine
    integer :: Nt,ii,jj,num,unit_number,it,i,j,n,kk,pp


    !==========================
    ! Stability criteria
    !   (gH)^(1/2)dt/dx<1/sqrt(8)


    num=400
    L=5e6
    dx=L/(num-1)
    !dx=0.025
    dy=dx
    !num=2/dx+1
    Lw=1e6/7  !10*7**2
    Nt=9500 ! Number of time steps
    gam=0.05! Asselin filter parameter
    dt=1.0!0.001
    hval=1000;
    g=9.81
    f=10e-4

  Stabcriteria=sqrt(g*hval)*dt/dx

  print*,'Stability= ',Stabcriteria



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
!open(unit=unit_number, file='matrix.txt', status='unknown')
!
!    do ii=1,num
!      x2 = dx*(ii - 1)
!        do jj=1,num
!         y2 = dy*(jj - 1)
!         h(ii,jj)=2*exp(-Lw*((x(ii)-0.5))**2-Lw*((y(jj)-0.5))**2)
!
!        end do
!        write(unit_number,*) h(ii,:)
!    end do

! 4 guassians
        do ii=1,num
          x2 = dx*(ii - 1)
            do jj=1,num
             y2 = dy*(jj - 1)
           !  h(ii,jj)=2*exp(-Lw*((x(ii)-0.5))**2-Lw*((y(jj)-0.5))**2)!+2*exp(-Lw*((x(ii)+0.5))**2-Lw*((y(jj)-0.5))**2)&
     !       h(ii,jj)=50000*exp(-Lw*((x(ii)-2.5e6))**2-Lw*((y(jj)-2.5e6))**2)!+2*exp(-Lw*((x(ii)+0.5))**2-Lw*((y(jj)-0.5))**2)&

 h(ii,jj)=10*exp(-((x(ii)-2e6)/Lw)**2-((y(jj)-2e6)/Lw)**2)!+2*exp(-Lw*((x(ii)+0.5))**2-Lw*((y(jj)-0.5))**2)&


            ! +2*exp(-Lw*((x(ii)-0.5))**2-Lw*((y(jj)+0.5))**2)

            end do
         !   write(unit_number,*) h(ii,:)
        end do

!h(100,:)=10



close(unit_number)
unit_number = 1

!     First integration step Euler method


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
! Save results
open(unit=unit_number, file='matrix2.txt', status='unknown')


!open(unit=unit_number, file='SWEdata.txt', status='unknown')
! Cgrid finite difference

!   do n=1,Nt
!       volume=0
!       Epot=0
!       Ekine=0
!
!!        u(1, :)     = u(num, :)
!!        u(:, 1)     = u(:, num)
!!        v(1, :)     = v(num, :)
!!        v(:, 1)     = v(:, num)
!!        h(1, :)     = h(num, :)
!!        h(:, 1)     = h(:, num)
!
!
!        do ii=2,num-1
!             it=ii;
!
!                 if (ii==2) then
!                     it=ii+1
!                  end if
!                  if(ii==num-1) then
!                  it=ii-1
!                  end if
!              do jj=2,num-1
!                j=jj
!                    if(jj==2)then
!                    j=jj+1
!
!                     end if
!
!                  if(jj==num-1) then
!                    j=jj-1
!                  end if
!
!               unew(ii,jj)=uold(ii,jj)+2*dt*(-g*(h(it+1,j)-h(it,j))/dx+0.25*f*(v(it,j)+v(it+1,j)+v(it+1,j-1)+v(it,j-1)));
!               vnew(ii,jj)=vold(ii,jj)+2*dt*(-g*(h(it,j+1)-h(it,j))/dy- 0.25*f*(u(it,j)+u(it,j+1)+u(it-1,j+1)+u(it-1,j)));
!               hnew(ii,jj)=hold(ii,jj)-2*dt*hval*((u(it,j)-u(it-1,j))/dx +(v(it,j)-v(it,j-1))/dy);
!
!               Volume=Volume+hnew(ii,jj)*dx*dy
!               Epot=Epot+0.5*g*hnew(ii,jj)**2*dx*dy
!               Ekine=Ekine+0.5*Hval*(unew(ii,jj)**2+vnew(ii,jj)**2)*dx*dy
!
!          end do
!      end do
!
!       do kk=1,num
!    do pp=1,num
!        unew(1,pp)   = unew(num-1,pp)    ! left = right-1
!        unew(num,pp) = unew(2,pp)        ! right = left+1
!        unew(kk,1)   = unew(kk,num-1)    ! bottom = top-1
!        unew(kk,num) = unew(kk,2)        ! top = bottom+1
!
!         vnew(1,pp)   = vnew(num-1,pp)    ! left = right-1
!        vnew(num,pp) = vnew(2,pp)        ! right = left+1
!        vnew(kk,1)   = vnew(kk,num-1)    ! bottom = top-1
!        vnew(kk,num) = vnew(kk,2)        ! top = bottom+1
!
!         hnew(1,pp)   = hnew(num-1,pp)    ! left = right-1
!        hnew(num,pp) = hnew(2,pp)        ! right = left+1
!        hnew(kk,1)   = hnew(kk,num-1)    ! bottom = top-1
!        hnew(kk,num) = hnew(kk,2)        ! top = bottom+1
!
!    end do
!        end do
!
!
!
!
!
!
!
!
!         ! Robert Asselin Filter
!      do ii=2,num-1
!        do jj=2,num-1
!          unew(ii,jj)=u(ii,jj)+gam*(uold(ii,jj)-2*u(ii,jj)+unew(ii,jj));
!          vnew(ii,jj)=v(ii,jj)+gam*(vold(ii,jj)-2*v(ii,jj)+vnew(ii,jj));
!          hnew(ii,jj)=h(ii,jj)+gam*(hold(ii,jj)-2*h(ii,jj)+hnew(ii,jj));
!!          u(ii,jj)=u(ii,jj)+gam*(uold(ii,jj)-2*u(ii,jj)+unew(ii,jj));
!!          v(ii,jj)=v(ii,jj)+gam*(vold(ii,jj)-2*v(ii,jj)+vnew(ii,jj));
!!          h(ii,jj)=h(ii,jj)+gam*(hold(ii,jj)-2*h(ii,jj)+hnew(ii,jj));
!        end do
!      end do
!
!
!
!
!      uold=u;
!      u=unew;
!      vold=v;
!      v=vnew;
!      hold=h;
!      h=hnew;
!
!!
!!      h(:,1) = h(:,2);      u(:,1) = u(:,2);       v(:,1) = -v(:,2);
!!       h(:,num) = h(:,num-1);  u(:,n) = u(:,num-1);   V(:,num) = -v(:,num-1);
!!       h(1,:) = h(2,:);      u(1,:) = -u(2,:);      v(1,:) = v(2,:);
!!       h(num,:) = h(num-1,:);  u(num,:) = -u(num-1,:);  v(num,:) = v(num-1,:);
!
! !! Rebounding walls
!!        do kk=1,num
!!         do pp=1,num
!!
!!          h(kk,1) = h(kk,2)
!!          u(kk,1) = u(kk,2)
!!          v(ii,1) = -v(kk,2);
!!          h(kk,num) = h(kk,num-1)
!!          u(kk,num) = u(kk,num-1)
!!          v(kk,num) = -v(kk,num-1)
!!          h(1,pp) = h(2,pp)
!!          u(1,pp) = -u(2,pp)
!!          v(1,pp) = v(2,pp);
!!          h(num,pp) = h(num-1,pp)
!!          u(num,pp) = -u(num-1,pp)
!!          v(num,pp) = v(num-1,pp)
!!
!!            end do
!!         end do
!! do kk=1,num
!!    do pp=1,num
!!   ! Periodic Boundary Conditions
!!        u(1, pp)     = u(num, pp)
!!        u(kk, 1)     = u(kk, num)
!!        v(1, pp)     = v(num, pp)
!!        v(kk, 1)     = v(kk, num)
!!        h(1, pp)     = h(num, pp)
!!        h(kk, 1)     = h(kk, num)
!!    end do
!!         end do
!
!
!!        u(1, :)     = u(num, :)
!!        u(:, 1)     = u(:, num)
!!        v(1, :)     = v(num, :)
!!        v(:, 1)     = v(:, num)
!!        h(1, :)     = h(num, :)
!!        h(:, 1)     = h(:, num)
!
!! do kk=1,num
!!    do pp=1,num
!!        u(1,pp)   = u(num-1,pp)    ! left = right-1
!!        u(num,pp) = u(2,pp)        ! right = left+1
!!        u(kk,1)   = u(kk,num-1)    ! bottom = top-1
!!        u(kk,num) = u(kk,2)        ! top = bottom+1
!!
!!         v(1,pp)   = v(num-1,pp)    ! left = right-1
!!        v(num,pp) = v(2,pp)        ! right = left+1
!!        v(kk,1)   = v(kk,num-1)    ! bottom = top-1
!!        v(kk,num) = v(kk,2)        ! top = bottom+1
!!
!!         h(1,pp)   = h(num-1,pp)    ! left = right-1
!!        h(num,pp) = h(2,pp)        ! right = left+1
!!        h(kk,1)   = h(kk,num-1)    ! bottom = top-1
!!        h(kk,num) = h(kk,2)        ! top = bottom+1
!!
!!    end do
!!        end do
!
!
!if(mod(n,100)==0)then
!
! do ii=1,num
!  write(unit_number,*) hnew(ii,:),Volume,Epot,Ekine
!  end do
!!print*,'h= ',h
!end if
!
!    print *, 'time = ',Nt-n
!
!
!enddo
!
!
!
!print*,'f=',f
!! Close the output file
!close(unit_number)





do n = 1, Nt
    volume = 0.0
    Epot = 0.0
    Ekine = 0.0

    ! Main update loop for interior points
    do ii = 1, num
        it = ii
        if (ii == 1) it = ii + 1
        if (ii == num) it = ii - 1

        do jj = 1, num
            j = jj
            if (jj == 1) j = jj + 1
            if (jj == num) j = jj - 1

!            unew(ii,jj) = uold(ii,jj) + 2*dt * (-g*(h(it+1,j) - h(it,j))/dx + 0.25*f*(v(it,j) + v(it+1,j) + v(it+1,j-1) + v(it,j-1)))
!            vnew(ii,jj) = vold(ii,jj) + 2*dt * (-g*(h(it,j+1) - h(it,j))/dy - 0.25*f*(u(it,j) + u(it,j+1) + u(it-1,j+1) + u(it-1,j)))
!            hnew(ii,jj) = hold(ii,jj) - 2*dt*hval * ((u(it,j) - u(it-1,j))/dx + (v(it,j) - v(it,j-1))/dy)

      unew(ii,jj)=uold(ii,jj)+2*dt*(-g*(h(it+1,j)-h(it,j))/dx+0.25*f*(v(it,j)+v(it+1,j)+v(it+1,j-1)+v(it,j-1)));
      vnew(ii,jj)=vold(ii,jj)+2*dt*(-g*(h(it,j+1)-h(it,j))/dy- 0.25*f*(u(it,j)+u(it,j+1)+u(it-1,j+1)+u(it-1,j)));
     hnew(ii,jj)=hold(ii,jj)-2*dt*hval*((u(it,j)-u(it-1,j))/dx +(v(it,j)-v(it,j-1))/dy);




            volume = volume + hnew(ii,jj)*dx*dy
            Epot = Epot + 0.5*g*hnew(ii,jj)**2*dx*dy
            Ekine = Ekine + 0.5*hval*(unew(ii,jj)**2 + vnew(ii,jj)**2)*dx*dy
        end do
    end do

    ! ✅ Apply periodic BCs to unew, vnew, hnew BEFORE copying to u, v, h
!    do kk = 1, num
!        do pp = 1, num
!            unew(1,pp)   = unew(num-1,pp)
!            unew(num,pp) = unew(2,pp)
!            unew(kk,1)   = unew(kk,num-1)
!            unew(kk,num) = unew(kk,2)
!
!            vnew(1,pp)   = vnew(num-1,pp)
!            vnew(num,pp) = vnew(2,pp)
!            vnew(kk,1)   = vnew(kk,num-1)
!            vnew(kk,num) = vnew(kk,2)
!
!            hnew(1,pp)   = hnew(num-1,pp)
!            hnew(num,pp) = hnew(2,pp)
!            hnew(kk,1)   = hnew(kk,num-1)
!            hnew(kk,num) = hnew(kk,2)
!
!
!
!            u(1,pp)   = u(num-1,pp)
!            u(num,pp) = u(2,pp)
!            u(kk,1)   = u(kk,num-1)
!            u(kk,num) = u(kk,2)
!
!            v(1,pp)   = v(num-1,pp)
!            v(num,pp) = v(2,pp)
!            v(kk,1)   = v(kk,num-1)
!            v(kk,num) = v(kk,2)
!
!            h(1,pp)   = h(num-1,pp)
!            h(num,pp) = h(2,pp)
!            h(kk,1)   = h(kk,num-1)
!            h(kk,num) = h(kk,2)
!
!
!
!
!
!
!        end do
!    end do

    ! Robert-Asselin Filter
    do ii = 1, num
        do jj = 1, num
            u(ii,jj) = u(ii,jj) + gam*(uold(ii,jj) - 2*u(ii,jj) + unew(ii,jj))
            v(ii,jj) = v(ii,jj) + gam*(vold(ii,jj) - 2*v(ii,jj) + vnew(ii,jj))
            h(ii,jj) = h(ii,jj) + gam*(hold(ii,jj) - 2*h(ii,jj) + hnew(ii,jj))
        end do
    end do

    u(1,:)   = u(num-1,:)    ! left = right-1
    u(num,:) = u(2,:)        ! right = left+1
    u(:,1)   = u(:,num-1)    ! bottom = top-1
    u(:,num) = u(:,2)        ! top = bottom+1

    unew(1,:)   = unew(num-1,:)    ! left = right-1
    unew(num,:) = unew(2,:)        ! right = left+1
    unew(:,1)   = unew(:,num-1)    ! bottom = top-1
    unew(:,num) = unew(:,2)        ! top = bottom+1

    v(1,:)   = v(num-1,:)    ! left = right-1
    v(num,:) = v(2,:)        ! right = left+1
    v(:,1)   = v(:,num-1)    ! bottom = top-1
    v(:,num) = v(:,2)        ! top = bottom+1

    vnew(1,:)   = vnew(num-1,:)    ! left = right-1
    vnew(num,:) = vnew(2,:)        ! right = left+1
    vnew(:,1)   = vnew(:,num-1)    ! bottom = top-1
    vnew(:,num) = vnew(:,2)        ! top = bottom+1

     h(1,:)   = h(num-1,:)    ! left = right-1
    h(num,:) = h(2,:)        ! right = left+1
    h(:,1)   = h(:,num-1)    ! bottom = top-1
    h(:,num) = h(:,2)        ! top = bottom+1

    hnew(1,:)   = hnew(num-1,:)    ! left = right-1
    hnew(num,:) = hnew(2,:)        ! right = left+1
    hnew(:,1)   = hnew(:,num-1)    ! bottom = top-1
    hnew(:,num) = hnew(:,2)        ! top = bottom+1

    ! Update fields
    uold = u
    u = unew
    vold = v
    v = vnew
    hold = h
    h = hnew

    ! Output
    if (mod(n,100) == 0) then
        do ii = 1, num
            write(unit_number,*) h(ii,:), volume, Epot, Ekine
        end do
    end if

    print *, 'time = ', Nt - n
end do

print *, 'dx = ', dx


end program
!
