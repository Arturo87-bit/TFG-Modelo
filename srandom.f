!=======================================================================
!   SUBROUTINE SRANDOM(is1)
!------------------------------------------------------------------
!  subroutina para inicializar del generador de numeros aleatorios
!  intrinseco de fortran 90 utilizando solo 1 semilla: is1
!-----------------------------------------------------------------------
      subroutine srandom(is1)
      implicit none
      integer(kind=4), intent(inout) :: is1
      integer(kind=4), allocatable :: iseed(:)
      integer(kind=4) :: n,i,itop
!
      itop=2**30
      if(is1 > itop) is1=is1-itop
      call random_seed(size=n)
!      write(*,*) 'numero de seeds=',n
      allocate(iseed(n))
      call srand(is1)
      do i=1,n
         iseed(i)=irand(0)
      enddo
      call random_seed(put=iseed(1:n))
      deallocate(iseed)
!
      end subroutine srandom
