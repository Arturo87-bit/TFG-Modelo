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
        

        integer, dimension(4,2) :: mov
        logical llego

        integer, allocatable :: idx_i(:), idx_j(:)
        
        integer k,total,rand_pos,tmp_i,tmp_j
        integer tiempo_llegada
        real*8 suma_tiempos
        integer cuenta_llegadas
        real*8 tiempo_medio, velocidad_media


        real*8, parameter :: pmin=0.0d0,pmax=1.0d0
        integer, parameter :: npuntos=100,Nmuestra=3000, Npas=300
        real*8 suma_infectados_t(Npas)
        integer cuenta_t(Npas)
        integer infectados
        real*8 media_infectados_t(Npas)
        real*8 suma_rho_total
        real*8 rho_sim
        real*8 suma_rho_t
        real*8 rho_media
        real*8 suma_tiempos2
        real*8 error_tiempo
        real*8 error_velocidad
        real*8 suma_rho_total2
        real*8 error_rho
        
        integer ic
        real*8, allocatable :: J_t(:)

        real*8 flujo_total
        integer duracion_pulso
        real*8 J_medio

! Para promediar sobre simulaciones
        real*8 suma_J_total, suma_J_total2
        real*8 J_media, error_J
        integer cuenta_J
        real*8 J_global
        real*8 error_Prob, error_J_global
        
        

        data mov / 1,0,-1,0,0,1,0,-1 /

        call srand(143920261)

        
        
        allocate(J_t(Npas))
       
        do m=1,4

        L=10*m

        allocate(poblacion(L,L),estado(L,L),estado_new(L,L),virus(L,L),
     &  movido(L,L))
        
        ic = L/2
        

          total = L*L
           allocate(idx_i(total), idx_j(total))

          write(*,*) "Simulacion para L =",L

          dp=(pmax-pmin)/real(npuntos,8)
          p=pmin-dp

        open(10,file='L'//trim(adjustl(int2str(m)))//'.dat')
        

          do o=1,npuntos

           p=p+dp
         Nper=0.d0
         
         suma_infectados_t = 0.d0
         cuenta_t = 0
         
         suma_tiempos2 = 0.d0
         suma_rho_total  = 0.d0
         suma_rho_total2 = 0.d0
         suma_tiempos = 0.d0
         cuenta_llegadas = 0
         
         suma_J_total = 0.d0
         suma_J_total2 = 0.d0
         cuenta_J = 0

         do Nbos=1,Nmuestra

          llego=.false.

         poblacion=0
         estado=0
         virus=0
         tiempo_llegada = -1
         J_t = 0.d0
         
         

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

            Npasos=300
            suma_rho_t = 0.d0

             t = 0
            llego = .false.

            do while ( any(estado == 2) .and. .not. llego )
            t = t + 1
            
            
            infectados = count(estado == 2)
            
            J_t(t) = real(count(estado(ic,:) == 2),8) / L
            
            suma_rho_t = suma_rho_t + real(infectados,8)/(L*L)
            
      suma_infectados_t(t)=suma_infectados_t(t)+real(infectados,8)/(L*L)
         cuenta_t(t) = cuenta_t(t) + 1

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
            if(virus(i,j)>=3) estado_new(i,j)=3

            end if
            end do
            end do

            estado=estado_new

!==================================
! PERCOLACION EN ESTADO INFECCIONSO
!==================================
           if(any(estado(L,:)==2)) then
            llego = .true.
           if (tiempo_llegada == -1) tiempo_llegada = t
           end if

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
         
         if(llego .and. tiempo_llegada > 0) then

          rho_sim = suma_rho_t / tiempo_llegada

          suma_rho_total = suma_rho_total + rho_sim
          suma_rho_total2 = suma_rho_total2 + rho_sim*rho_sim

         end if

         if(llego) then
          Nper = Nper + 1
          suma_tiempos = suma_tiempos + tiempo_llegada
          cuenta_llegadas = cuenta_llegadas + 1
          suma_tiempos2 = suma_tiempos2 + tiempo_llegada**2
         end if
         
         do t=1,Npasos
          if(cuenta_t(t)>0) then
           media_infectados_t(t) = suma_infectados_t(t)/cuenta_t(t)
          else
          media_infectados_t(t) = 0.d0
          end if
          
          
      
         
         end do   ! muestras
         
         flujo_total = sum(J_t)
         duracion_pulso = count(J_t > 0.d0)

         if(duracion_pulso > 0) then
         J_medio = flujo_total / duracion_pulso
         else
          J_medio = 0.d0
         end if
         
         
         
         

         if(cuenta_llegadas > 0) then
          rho_media = suma_rho_total / cuenta_llegadas
          
          error_rho=sqrt((suma_rho_total2/cuenta_llegadas-rho_media**2) 
     & / cuenta_llegadas )
         else
         rho_media = 0.d0
         error_rho = 0.d0
         end if

         

          if(cuenta_llegadas > 0) then
          tiempo_medio = suma_tiempos / cuenta_llegadas
          velocidad_media = L / tiempo_medio
          error_tiempo = sqrt( (suma_tiempos2/cuenta_llegadas  
     & - tiempo_medio**2)/ cuenta_llegadas )
          error_velocidad = L * error_tiempo / (tiempo_medio**2)
          else
          tiempo_medio = 0.d0
          velocidad_media = 0.d0
           error_tiempo = 0.d0
           error_velocidad = 0.d0
          end if
          
          if(duracion_pulso > 0) then
          suma_J_total = suma_J_total + J_medio
          suma_J_total2 = suma_J_total2 + J_medio*J_medio
          end if
          
          if(cuenta_llegadas > 0) then
          J_media = suma_J_total / cuenta_llegadas
          error_J = sqrt( (suma_J_total2/cuenta_llegadas - J_media**2)
     & / cuenta_llegadas )
          else
          J_media = 0.d0
          error_J = 0.d0
          end if
         

          end do
         Prob=Nper/Nmuestra
         
         J_global = Prob * J_media
         
         error_Prob = sqrt(Prob*(1.d0 - Prob)/Nmuestra)

        error_J_global = sqrt((J_media*error_Prob)**2+(Prob*error_J)**2)
        

        write(10,*) p, Prob, tiempo_medio,error_tiempo,
     &  velocidad_media,error_velocidad,rho_media/(3.d0/L),
     &  error_rho/(3.d0/L),J_media, error_J, J_global,error_J_global
        

         end do

       

      

       
         
         

        deallocate(poblacion,estado,estado_new,virus,movido,idx_i,idx_j)
        

         close(10)

        end do
        
        deallocate(J_t)

        

        contains

        function int2str(n) result(str)
         integer,intent(in)::n
         character(len=5)::str
         write(str,'(I5)')n
         end function

        end program Percolacion


        include 'srandom.f'
        include 'linfit.f'
