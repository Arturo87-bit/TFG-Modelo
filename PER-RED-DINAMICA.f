        module parametros
        implicit none

        integer, allocatable :: poblacion(:,:),estado(:,:)
        integer, allocatable :: estado_new(:,:),virus(:,:)
        logical, allocatable :: movido(:,:)

        integer L,Npasos
        real*8 Nper,dp,Prob

        end module parametros


        program Percolacion
        use parametros
        implicit none

        integer i,j,t,Nbos,o,m
        integer n_x,n_y,movil
        real*8 d,p
        real*8 a,b,da,db,r

        integer, dimension(4,2) :: mov
        logical llego

        integer, allocatable :: idx_i(:), idx_j(:)
        integer k,total,rand_pos,tmp_i,tmp_j

        real*8, parameter :: pmin=0.0d0,pmax=1.0d0
        integer, parameter :: npuntos=100,Nmuestra=10000
        real*8 x(npuntos),y(npuntos),e(npuntos)

        data mov / 1,0,-1,0,0,1,0,-1 /

        call srand(143920261)

        open(20,file='DATOS_AJUSTE.dat')

         do m = 1,4

        L=10*m

        allocate(poblacion(L,L),estado(L,L),estado_new(L,L),virus(L,L),
     &  movido(L,L))

          total = L*L
           allocate(idx_i(total), idx_j(total))

          write(*,*) "Simulacion para L =",L

          dp=(pmax-pmin)/real(npuntos,8)
          p=pmin-dp

          open(10,file='L'//trim(adjustl(int2str(m)))//'.dat')

          do o=1,npuntos

           p=p+dp
         Nper=0.d0

         do Nbos=1,Nmuestra

          llego=.false.

         poblacion=0
         estado=0
         virus=0

!----------------------------
! Generar poblacion
!----------------------------
         do i=1,L
         do j=1,L
         call random_number(d)
           if(d<=p) then
           poblacion(i,j)=1
           estado(i,j)=1
         end if
         end do
         end do

!----------------------------
! Infectar primera columna
!----------------------------
          do j=1,L
           if(poblacion(1,j)==1) estado(1,j)=2
             end do

            

            t = 0
            llego = .false.

            do while ( any(estado == 2) .and. .not. llego )
            t = t + 1


            estado_new=estado

!============================
! PROPAGACION DE LA INFECCION
!============================
            do i=1,L
            do j=1,L
            if(estado(i,j)==2) then
            do movil=1,4
            n_x=i+mov(movil,1)
            n_y=j+mov(movil,2)
            if(n_x>0 .and. n_x<=L .and. n_y>0 .and. n_y<=L) then
            if(estado(n_x,n_y)==1) estado_new(n_x,n_y)=2
            end if
            end do
            end if
            end do
            end do

!============================
! EVOLUCION DEL ESTADO INFECCION
!============================
            do i=1,L
            do j=1,L
            if(estado(i,j)==2) then
            virus(i,j)=virus(i,j)+1
            if(virus(i,j)>=2) estado_new(i,j)=3 !Aqui variamos tau

            end if
            end do
            end do

            estado=estado_new

!==================================
! PERCOLACION EN ESTADO INFECCIONSO
!==================================
           if(any(estado(L,:)==2)) llego=.true.

!============================
! MOVIMIENTO DE LA RED
!============================

           movido = .false.

! Construir lista
           k=0
           do i=1,L
           do j=1,L
             k=k+1
             idx_i(k)=i
             idx_j(k)=j
             end do
           end do

! Barajar, elegimos aleatoriamente la posición
          do k=total,2,-1
           call random_number(d)
           rand_pos = int(d*k)+1

           tmp_i = idx_i(k)
           tmp_j = idx_j(k)

          idx_i(k) = idx_i(rand_pos)
          idx_j(k) = idx_j(rand_pos)

          idx_i(rand_pos) = tmp_i
          idx_j(rand_pos) = tmp_j
           end do

! Movimiento
          do k=1,total

          i = idx_i(k)
          j = idx_j(k)

         if(estado(i,j)/=0 .and. .not. movido(i,j)) then

         call random_number(d)
          movil=int(d*4)+1

          n_x=i+mov(movil,1)
          n_y=j+mov(movil,2)

          if(n_x>0 .and. n_x<=L .and. n_y>0 .and. n_y<=L) then

          if(poblacion(n_x,n_y)==0) then

               poblacion(n_x,n_y)=poblacion(i,j)
               estado(n_x,n_y)=estado(i,j)
               virus(n_x,n_y)=virus(i,j)

               movido(n_x,n_y)=.true.

               poblacion(i,j)=0
               estado(i,j)=0
               virus(i,j)=0

               end if

             end if

           end if

         end do

         end do   ! pasos

         if(llego) Nper=Nper+1

         end do   ! muestras

         Prob=Nper/Nmuestra

         x(o)=p
         y(o)=Prob
        e(o)=sqrt(Prob*(1.d0-Prob)/Nmuestra)

        write(10,*) p,Prob,e(o)

         end do

        call linfit(npuntos,x,y,e,a,b,da,db,r)
        write(20, *) L,a,da,log(a),b

        deallocate(poblacion,estado,estado_new,virus,movido,idx_i,idx_j)

         close(10)

        end do

        close(20)

        contains

        function int2str(n) result(str)
         integer,intent(in)::n
         character(len=5)::str
         write(str,'(I5)')n
         end function

        end program Percolacion


        include 'srandom.f'
        include 'linfit.f'
