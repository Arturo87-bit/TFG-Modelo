      program sir_bootstrap
      implicit none

      integer, parameter :: Npasos=600, Nboot=60
      integer b,t,iseed,t_ini

      real*8 Sar(Npasos), Iar(Npasos), Rar(Npasos)
      real*8 beta, gamma, R0, lambda
      real*8 beta_arr(Nboot), gamma_arr(Nboot), R0_arr(Nboot)
      real*8 beta_mean, gamma_mean, R0_mean
      real*8 beta_std, gamma_std, R0_std
      real*8 N0

      iseed=143920261
      call init_random_seed(iseed)

      t_ini=60

      do b=1,Nboot

         call simular_red(Sar,Iar,Rar)

         N0 = Sar(1)+Iar(1)+Rar(1)

         do t=1,Npasos
            Sar(t)=Sar(t)/N0
            Iar(t)=Iar(t)/N0
            Rar(t)=Rar(t)/N0
         end do

         call ajuste_exponencial(Iar,t_ini,lambda)
         call ajustar_sir(Sar,Iar,Rar,lambda,beta,gamma)

         if (gamma.gt.1d-12) then
            R0 = beta/gamma
         else
            R0 = 0.d0
         end if

         beta_arr(b)=beta
         gamma_arr(b)=gamma
         R0_arr(b)=R0

         write(*,*) "iter",b,"beta=",beta,"gamma=",gamma,
     &              "R0=",R0,"lambda=",lambda

      end do

      call media_std(beta_arr,Nboot,beta_mean,beta_std)
      call media_std(gamma_arr,Nboot,gamma_mean,gamma_std)
      call media_std(R0_arr,Nboot,R0_mean,R0_std)

      write(*,*)
      write(*,*) "======== RESULTADO FINAL ========"
      write(*,*) "beta  =",beta_mean," +- ",beta_std
      write(*,*) "gamma =",gamma_mean," +- ",gamma_std
      write(*,*) "R0    =",R0_mean," +- ",R0_std

      end


C========================================================
      subroutine simular_red(Sar,Iar,Rar)
      implicit none

      integer, parameter :: L=150, Npasos=600
      integer i,j,t,m,k,total
      integer nx,ny,movil,rand_pos
      integer tmp_i,tmp_j

      integer, allocatable :: poblacion(:,:),estado(:,:),estado_new(:,:)
      integer, allocatable :: virus(:,:)
      logical, allocatable :: movido(:,:)
      integer, allocatable :: idx_i(:), idx_j(:)

      integer mov(4,2)
      data mov /1,0,-1,0,0,1,0,-1/

      real*8 Sar(Npasos),Iar(Npasos),Rar(Npasos)
      real*8 d,p

      allocate(poblacion(L,L), estado(L,L), estado_new(L,L))
      allocate(virus(L,L), movido(L,L))
      allocate(idx_i(L*L), idx_j(L*L))

      p=0.043
      poblacion=0
      estado=0
      virus=0

      do i=1,L
      do j=1,L
         call random_number(d)
         if(d.le.p) then
            poblacion(i,j)=1
            estado(i,j)=1
         end if
      end do
      end do

      do j=1,L
         if(poblacion(1,j).eq.1) estado(1,j)=2
      end do

      total=L*L

      do t=1,Npasos

         estado_new=estado

! contagio
         do i=1,L
         do j=1,L
            if(estado(i,j).eq.2) then
               do m=1,4
                  nx=i+mov(m,1)
                  ny=j+mov(m,2)
                  if(nx.gt.0 .and. nx.le.L .and.
     &               ny.gt.0 .and. ny.le.L) then
                     if(estado(nx,ny).eq.1) estado_new(nx,ny)=2
                  end if
               end do
            end if
         end do
         end do

! evolucion
         do i=1,L
         do j=1,L
            if(estado(i,j).eq.2) then
               virus(i,j)=virus(i,j)+1
               if(virus(i,j).ge.100) estado_new(i,j)=3
            end if
         end do
         end do

         estado=estado_new

! movimiento
         movido=.false.

         k=0
         do i=1,L
         do j=1,L
            k=k+1
            idx_i(k)=i
            idx_j(k)=j
         end do
         end do

         do k=total,2,-1
            call random_number(d)
            rand_pos=int(d*dble(k))+1

            tmp_i=idx_i(k)
            tmp_j=idx_j(k)

            idx_i(k)=idx_i(rand_pos)
            idx_j(k)=idx_j(rand_pos)

            idx_i(rand_pos)=tmp_i
            idx_j(rand_pos)=tmp_j
         end do

         do k=1,total

            i=idx_i(k)
            j=idx_j(k)

            if(estado(i,j).ne.0 .and. .not.movido(i,j)) then

               call random_number(d)
               movil=int(d*4.d0)+1

               nx=i+mov(movil,1)
               ny=j+mov(movil,2)

               if(nx.gt.0 .and. nx.le.L .and.
     &            ny.gt.0 .and. ny.le.L) then
                  if(poblacion(nx,ny).eq.0) then

                     poblacion(nx,ny)=poblacion(i,j)
                     estado(nx,ny)=estado(i,j)
                     virus(nx,ny)=virus(i,j)

                     movido(nx,ny)=.true.

                     poblacion(i,j)=0
                     estado(i,j)=0
                     virus(i,j)=0

                  end if
               end if

            end if

         end do

         Sar(t)=dble(count(estado.eq.1))
         Iar(t)=dble(count(estado.eq.2))
         Rar(t)=dble(count(estado.eq.3))

      end do

      deallocate(poblacion,estado,estado_new,virus,movido,idx_i,idx_j)

      end


C========================================================
      subroutine ajustar_sir(Sar,Iar,Rar,lambda,beta,gamma)
      implicit none

      integer, parameter :: Npasos=600
      integer, parameter :: Ncoarse=80, Nfine=120
      real*8 Sar(Npasos), Iar(Npasos), Rar(Npasos)
      real*8 beta,gamma,lambda
      real*8 chi2,chi2_best,penal
      real*8 beta_try,gamma_try,beta_best,gamma_best
      real*8 beta_min,beta_max,gamma_min,gamma_max
      real*8 db,dg,lambda_mod,Nini
      real*8 Sm(Npasos),Im(Npasos),Rm(Npasos)
      integer t,ib,ig

      chi2_best=1.d99
      beta_best=0.d0
      gamma_best=0.d0

! primera busqueda amplia
      beta_min=0.001d0
      beta_max=2.500d0
      gamma_min=0.001d0
      gamma_max=2.000d0

      db=(beta_max-beta_min)/dble(Ncoarse-1)
      dg=(gamma_max-gamma_min)/dble(Ncoarse-1)

      do ib=1,Ncoarse
         beta_try=beta_min+dble(ib-1)*db

         do ig=1,Ncoarse
            gamma_try=gamma_min+dble(ig-1)*dg

            call integrar_sir(Sar(1),Iar(1),Rar(1),
     &           beta_try,gamma_try,Sm,Im,Rm)

            chi2=0.d0

            do t=1,Npasos
               chi2=chi2
     &         +(Sar(t)-Sm(t))**2
     &         +(Iar(t)-Im(t))**2
     &         +(Rar(t)-Rm(t))**2
            end do

            Nini=Sar(1)+Iar(1)+Rar(1)
            if(Nini.gt.0.d0) then
               lambda_mod=beta_try*Sar(1)/Nini-gamma_try
               penal=(lambda-lambda_mod)**2
               chi2=chi2+0.02d0*penal
            end if

            if(chi2.lt.chi2_best) then
               chi2_best=chi2
               beta_best=beta_try
               gamma_best=gamma_try
            end if

         end do
      end do

! segunda busqueda fina alrededor del mejor punto
      beta_min=max(0.0001d0,beta_best-0.080d0)
      beta_max=beta_best+0.080d0
      gamma_min=max(0.0001d0,gamma_best-0.080d0)
      gamma_max=gamma_best+0.080d0

      db=(beta_max-beta_min)/dble(Nfine-1)
      dg=(gamma_max-gamma_min)/dble(Nfine-1)

      do ib=1,Nfine
         beta_try=beta_min+dble(ib-1)*db

         do ig=1,Nfine
            gamma_try=gamma_min+dble(ig-1)*dg

            call integrar_sir(Sar(1),Iar(1),Rar(1),
     &           beta_try,gamma_try,Sm,Im,Rm)

            chi2=0.d0

            do t=1,Npasos
               chi2=chi2
     &         +(Sar(t)-Sm(t))**2
     &         +(Iar(t)-Im(t))**2
     &         +(Rar(t)-Rm(t))**2
            end do

            Nini=Sar(1)+Iar(1)+Rar(1)
            if(Nini.gt.0.d0) then
               lambda_mod=beta_try*Sar(1)/Nini-gamma_try
               penal=(lambda-lambda_mod)**2
               chi2=chi2+0.02d0*penal
            end if

            if(chi2.lt.chi2_best) then
               chi2_best=chi2
               beta_best=beta_try
               gamma_best=gamma_try
            end if

         end do
      end do

      beta=beta_best
      gamma=gamma_best

      end


!========================================================
      subroutine ajuste_exponencial(Iar,t_ini,lambda)
      implicit none

      integer, parameter :: Npasos=600
      integer t_ini,t,nval
      real*8 Iar(Npasos),lambda
      real*8 sumx,sumy,sumxx,sumxy,x,y,den

      sumx=0.d0
      sumy=0.d0
      sumxx=0.d0
      sumxy=0.d0
      nval=0

      do t=1,t_ini
         if(Iar(t).gt.0.d0) then
            x=dble(t)
            y=log(Iar(t))
            sumx=sumx+x
            sumy=sumy+y
            sumxx=sumxx+x*x
            sumxy=sumxy+x*y
            nval=nval+1
         end if
      end do

      den=dble(nval)*sumxx-sumx*sumx

      if(nval.gt.1 .and. abs(den).gt.1.d-14) then
         lambda=(dble(nval)*sumxy-sumx*sumy)/den
      else
         lambda=0.d0
      end if

      end


!========================================================
      subroutine integrar_sir(S0,I0,R0,beta,gamma,S,I,R)
      implicit none

      integer, parameter :: Npasos=600
      real*8 S(Npasos),I(Npasos),R(Npasos)
      real*8 S0,I0,R0,beta,gamma

      real*8 y(3),dydx(3),yout(3),x,h
      integer t

      real*8 beta_glob,gamma_glob
      common /param_sir/ beta_glob,gamma_glob

      external derivs,rk4

      beta_glob=beta
      gamma_glob=gamma

      y(1)=S0
      y(2)=I0
      y(3)=R0

      S(1)=S0
      I(1)=I0
      R(1)=R0

      do t=2,Npasos
         x=dble(t-1)
         h=1.d0

         call derivs(x,y,dydx)
         call rk4(y,dydx,3,x,h,yout,derivs)

         y=yout

         if(y(1).lt.0.d0) y(1)=0.d0
         if(y(2).lt.0.d0) y(2)=0.d0
         if(y(3).lt.0.d0) y(3)=0.d0

         S(t)=y(1)
         I(t)=y(2)
         R(t)=y(3)
      end do

      end


!========================================================
      subroutine derivs(x,y,dydx)
      implicit none

      real*8 x,y(3),dydx(3),N
      real*8 beta_glob,gamma_glob
      common /param_sir/ beta_glob,gamma_glob

      N=y(1)+y(2)+y(3)

      if(x.lt.-1.d99) print*,x

      if(N.gt.0.d0) then
         dydx(1)=-beta_glob*y(1)*y(2)/N
         dydx(2)= beta_glob*y(1)*y(2)/N-gamma_glob*y(2)
         dydx(3)= gamma_glob*y(2)
      else
         dydx(1)=0.d0
         dydx(2)=0.d0
         dydx(3)=0.d0
      end if

      end


!========================================================
      subroutine media_std(x,n,media,std)
      implicit none

      integer n,i
      real*8 x(n),media,std

      media=0.d0
      do i=1,n
         media=media+x(i)
      end do
      media=media/dble(n)

      std=0.d0
      do i=1,n
         std=std+(x(i)-media)**2
      end do

      if(n.gt.1) then
         std=sqrt(std/dble(n-1))
      else
         std=0.d0
      end if

      end


!========================================================
      subroutine rk4(y,dydx,n,x,h,yout,derivs)
      implicit none

      integer n,i
      real*8 y(n),dydx(n),yout(n)
      real*8 x,h
      external derivs

      real*8 hh,h6,xh
      real*8 dym(50),dyt(50),yt(50)

      hh=h*0.5d0
      h6=h/6.d0
      xh=x+hh

      do i=1,n
         yt(i)=y(i)+hh*dydx(i)
      end do

      call derivs(xh,yt,dyt)

      do i=1,n
         yt(i)=y(i)+hh*dyt(i)
      end do

      call derivs(xh,yt,dym)

      do i=1,n
         yt(i)=y(i)+h*dym(i)
         dym(i)=dyt(i)+dym(i)
      end do

      call derivs(x+h,yt,dyt)

      do i=1,n
         yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))
      end do

      end


!========================================================
      subroutine init_random_seed(iseed)
      implicit none

      integer iseed
      integer n,i
      integer, allocatable :: seed(:)

      call random_seed(size=n)
      allocate(seed(n))

      do i=1,n
         seed(i)=iseed+37*(i-1)
      end do

      call random_seed(put=seed)

      deallocate(seed)

      end
