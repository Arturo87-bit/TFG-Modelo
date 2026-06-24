!=======================================================================
!  SUBROUTINE LINFIT
!-----------------------------------------------------------------------
!  Realiza ajustes a una recta por minimos cuadrados; y=a*x+b
!  Si no se proporcionan errores en y (hay que poner e=1.d0), se
!  se estiman los errores de la recta a partir de chi^2.
!-----------------------------------------------------------------------
!  n    -- tamaño de los vectores de datos
!  x(n) -- vector con la coordenada independiente
!  y(n) -- vector con la coordenada y
!  e(n) -- vector con los errores en y. Si se han estimado: e=1.d0.
!  a    -- pendiente de la recta de ajuste
!  b    -- termino independiente de la recta de ajuste
!  da   -- estimacion del error en la pendiente
!  db   -- estimacion del error en b
!  r    -- coeficiente de regresion
!-----------------------------------------------------------------------
      subroutine linfit(n,x,y,e,a,b,da,db,r)
      implicit none
      integer(kind=4), intent(in) :: n
      real(kind=8),intent(in)     :: x(n),y(n),e(n)
      real(kind=8), intent(out)   :: a,b,da,db,r
      integer(kind=4) :: i
      real(kind=8)    :: sx0,sx1,sx2,sy1,sy2,sxy,xx,yy
      real(kind=8)    :: cte,den,dey,chi2,arg,errn,ay
!
      sx0=0.d0
      sx1=0.d0
      sx2=0.d0
      sy1=0.d0
      sy2=0.d0
      sxy=0.d0
      do i=1,n
        xx=x(i)
        yy=y(i)
! si e=0.0d0, no hay errores y lo ponemos a 1.d0.
        if(e(i).ne.0.d0) then
           cte=1.d0/e(i)/e(i)
        else
           cte=1.d0
        endif
        sx0=sx0+cte
        sx1=sx1+xx*cte
        sx2=sx2+xx*xx*cte
        sy1=sy1+yy*cte
        sy2=sy2+yy*yy*cte
        sxy=sxy+xx*yy*cte
      enddo
      den=sx2*sx0-sx1*sx1
      dey=sy2*sx0-sy1*sy1
      a=(sx0*sxy-sx1*sy1)/den
      ay=(sx0*sxy-sx1*sy1)/dey
      b=(sx2*sy1-sx1*sxy)/den
      da=dsqrt(sx0/den)
      db=dsqrt(sx2/den)
      r=dsqrt(a*ay)    ! coeficiente de regresion
!
      chi2=0.d0
      do i=1,n
         if(e(i).ne.0.d0) then
           arg=(y(i)-a*x(i)-b)/e(i)
         else
           arg=(y(i)-a*x(i)-b)
         endif
         chi2=chi2+arg*arg
      enddo
! Si no se proporcionan errores, se estiman a partir de chi2,
! suponiendo los errores de todos los puntos iguales.
      if(dabs(e(1)-1.d0).lt.1.d-12)then
         errn=dsqrt(chi2/dble(n-2))
         da=da*errn
         db=db*errn
      endif
      end subroutine linfit
!----------------------------------------------------------
