      module fmcp_dijkstra
      use fmcp_data
      implicit none
    
      contains
    
      subroutine calculate_dijkstra()
    
      !==============================================================
      ! Calcular caminos mínimos desde cada interactuante
      ! hacia todos los residuos de nlist1
      !
      ! Además:
      !
      ! 1) calcular path_sum(v), como antes
      !
      ! 2) identificar los residuos que participan en cada camino
      !
      ! 3) para cada residuo:
      !       - contar cuántos caminos contienen al residuo
      !       - acumular el peso de esos caminos
      !
      ! Los residuos pertenecientes al grupo destino NO se cuentan
      ! como residuos intervinientes.
      !
      ! La información se acumula por:
      !
      !       cadena + grupo + residuo
      !
      !==============================================================
    
    
      !==============================================================
      ! Variables locales para el análisis de los caminos
      !==============================================================
    
      integer, allocatable :: path_c(:,:,:)
      real*8, allocatable :: path_w(:,:,:)
    
      logical, allocatable :: path_included(:)
    
      integer :: current
      integer :: istat
      integer :: target_group
      integer :: target_chain
      integer :: node
      integer :: prev
      integer :: ierr
      integer :: nsort, temp_i
      integer, dimension(nmax) :: sort_c, sort_node
 
      real*8 :: path_score,temp_w
      real*8, dimension(nmax) :: sort_w

      character*100 filename

      !==============================================================
      ! Arrays para almacenar:
      !
      ! path_c(chain,group,residue)
      !
      ! path_w(chain,group,residue)
      !
      !==============================================================
    
      allocate(path_c(nchains,ngroup,natom),stat=ierr)
    
      if (ierr /= 0) then
         write(*,*) 'ERROR: no se pudo asignar path_c'
         stop
      endif
    
      allocate(path_w(nchains,ngroup,natom),stat=ierr)
    
      if (ierr /= 0) then
         write(*,*) 'ERROR: no se pudo asignar path_w'
         stop
      endif
    
      allocate(path_included(natom),stat=ierr)
    
      if (ierr /= 0) then
         write(*,*) 'ERROR: no se pudo asignar path_included'
         stop
      endif
    

      path_c  = 0
      path_w = 0.0d0
    
    
      !==============================================================
      ! Inicializar suma de inversas de los caminos
      !==============================================================
    
      path_sum = 0.0d0
    
    
      !==============================================================
      ! Calcular caminos mínimos desde cada interactuante
      ! hacia cada target de nlist1
      !==============================================================
    
      do js = 1,nres_interacting
    
         !-----------------------------------------------------------
         ! Residuo interactuante que actúa como origen
         !-----------------------------------------------------------
    
         u = interacting_residues(js)
    
    
         !-----------------------------------------------------------
         ! Inicialización de Dijkstra
         !-----------------------------------------------------------
    
         min_dist = 1.0d15
         visited = .false.
         predecessor = -1
    
         min_dist(u) = 0.0d0
    
    
         !===========================================================
         ! Dijkstra desde el interactuante u
         !===========================================================
    
         do i = 1,natom
    
            min_value = 1.0d15
            min_idx   = -1
    
            !--------------------------------------------------------
            ! Buscar nodo no visitado con distancia mínima
            !--------------------------------------------------------
    
            do v = 1,natom
    
               if (.not.visited(v).and.min_dist(v)<min_value) then
    
                  min_value = min_dist(v)
                  min_idx   = v
    
               endif
    
            enddo
    
    
            !--------------------------------------------------------
            ! No quedan nodos alcanzables
            !--------------------------------------------------------
    
            if (min_idx == -1) exit
    
            visited(min_idx) = .true.
    
    
            !--------------------------------------------------------
            ! Relajación de aristas
            !--------------------------------------------------------
    
            do v = 1,natom
    
             if (.not.visited(v).and.dpath(min_idx,v) < 1.0d15) then
    
              if (min_dist(v) > min_dist(min_idx) + 
     &            dpath(min_idx,v)) then
    
                 min_dist(v) = min_dist(min_idx) + 
     &                         dpath(min_idx,v)
    
                 predecessor(v) = min_idx
    
              endif
    
             endif
    
            enddo
    
         enddo
    
         !===========================================================
         ! Reconstruir el camino hacia cada target de nlist1
         !===========================================================
    
         do jt = 1,ncon
    
            v = nlist1(jt)
    
    
            !--------------------------------------------------------
            ! Verificar que el target sea alcanzable
            !--------------------------------------------------------
    
            if (min_dist(v) < 1.0d15) then
    
    
               !-----------------------------------------------------
               ! Contribución de este camino al score del target
               !
               ! Es exactamente el peso utilizado en path_sum.
               !-----------------------------------------------------
    
               path_score = pinter(u) / min_dist(v)
    
    
               !-----------------------------------------------------
               ! Mantener el cálculo original del score
               !-----------------------------------------------------
    
               path_sum(v) = path_sum(v) + path_score
    
    
               !=====================================================
               ! Información del grupo al que llega este camino
               !=====================================================
    
               target_chain = chain_residue(v)
               target_group = group_residue(v)
    
    
               !-----------------------------------------------------
               ! Inicializar control de residuos incluidos.
               !
               ! Esto garantiza que un residuo se cuente solamente
               ! una vez dentro de ESTE camino.
               !-----------------------------------------------------
    
               path_included = .false.
    
    
               !=====================================================
               ! Reconstrucción del camino:
               !
               ! v
               ! |
               ! predecessor(v)
               ! |
               ! predecessor(...)
               ! |
               ! u
               !
               ! Se recorre hacia atrás hasta llegar al origen.
               !=====================================================
    
               current = v
    
               do while (current /= -1)
    
                  !--------------------------------------------------
                  ! Los residuos pertenecientes al grupo destino
                  ! NO se incluyen en la métrica de participación.
                  !--------------------------------------------------
    
                  if (group_residue(current) /= target_group .or. 
     &                chain_residue(current) /= target_chain) then
    
    
                     !-----------------------------------------------
                     ! Evitar contar dos veces el mismo residuo
                     ! dentro del mismo camino.
                     !-----------------------------------------------
    
                     if (.not.path_included(current)) then
    
                      path_included(current) = .true.
  
                      !---------------------------------------------
                      ! Contar aparición del residuo en este camino
                      !---------------------------------------------
   
                      path_c(target_chain,target_group,current) = 
     &                path_c(target_chain,target_group,current) + 1
    
                      !---------------------------------------------
                      ! Acumular el peso de este camino
                      !---------------------------------------------
    
                      path_w(target_chain,target_group,current) =
     &                path_w(target_chain,target_group,current) +
     &                path_score
    
                     endif
    
                  endif
    
                  !--------------------------------------------------
                  ! Llegamos al origen
                  !--------------------------------------------------
    
                  if (current == u) exit
    
                  !--------------------------------------------------
                  ! Avanzar hacia el nodo anterior del camino
                  !--------------------------------------------------
    
                  prev = predecessor(current)
    
                  current = prev
    
               enddo
    
            endif
    
         enddo
    
      enddo
    
    
      !=========================================================
      ! Inicializar scores de los grupos
      !=========================================================
    
      grp_score = 0.0d0
   
      !=========================================================
      ! Contar cantidad de residuos de cada grupo en cada cadena
      !=========================================================

      ngrp_member = 0

      do v = 1,natom

         ichain = chain_residue(v)
         igrp   = group_residue(v)

         if (ichain .gt. 0 .and. ichain .le. nchains .and.
     &       igrp   .gt. 0 .and. igrp .le. ngroup) then

            ngrp_member(ichain,igrp) =
     &         ngrp_member(ichain,igrp) + 1

         endif

      enddo 
    
      !=========================================================
      ! Sumar path_sum de cada residuo al grupo correspondiente
      !
      ! group_residue(v) = grupo del residuo v
      ! chain_residue(v) = cadena del residuo v
      !=========================================================

      do jt = 1,ncon

         v = nlist1(jt)

         ichain = chain_residue(v)
         igrp   = group_residue(v)

         grp_score(ichain,igrp) = grp_score(ichain,igrp) +
     &                            path_sum(v)

      enddo



      !=========================================================
      ! Normalizar el score de cada grupo por su numero
      ! de residuos y escribir
      !=========================================================

      open(26,file='score.groups.dat',status='replace')
    
      write(26,*)
      write(26,*) '# ============================================'
      write(26,*) '# SCORE POR GRUPO'
      write(26,*) '# ============================================'
      write(26,*)
      write(26,*)

      do ichain = 1,nchains

         write(26,*) '# Cadena   Grupo   Score'

         do igrp = 1,ngroup

            if (ngrp_member(ichain,igrp) .gt. 0) then

            grp_score(ichain,igrp) = grp_score(ichain,igrp) /
     &           float(ngrp_member(ichain,igrp))

            write(26,96) ichain,igrp,grp_score(ichain,igrp)

            endif

         enddo

      enddo


      close(26)
    
   96 format(2(2x,I4),2X,f12.3) 
    
!==============================================================
! Escribir participacion de residuos en los caminos
! Un archivo independiente para cada cadena y grupo
!
! Solamente se escriben residuos que participaron al menos
! en un camino.
!==============================================================
!==============================================================
!     Crear directorio para los resultados de los grupos
!==============================================================

      call execute_command_line('mkdir -p groups',exitstat=istat)

      if (istat /= 0) then
         write(*,*) 'ERROR: no se pudo crear el directorio groups'
         stop
      endif

!==============================================================
!     Escribir archivos por cadena y grupo
!     Ordenados por W_score de mayor a menor
!==============================================================

      do ichain = 1,nchains

         do igrp = 1,ngroup

            write(filename,'("groups/path.chain.",I0,
     &                       ".group.",I0,".dat")') ichain,igrp

            open(27,file=filename,status='replace')

!--------------------------------------------------------------
!           Guardar residuos que participaron en caminos
!--------------------------------------------------------------

            nsort = 0

            do node = 1,natom

               if (path_c(ichain,igrp,node) .gt. 0) then

                  nsort = nsort + 1

                  sort_node(nsort) = node
                  sort_c(nsort) = path_c(ichain,igrp,node)
                  sort_w(nsort) = path_w(ichain,igrp,node)

               endif

            enddo

!--------------------------------------------------------------
!           Ordenar por W_score de mayor a menor
!--------------------------------------------------------------

            do i = 1,nsort-1

               do j = i+1,nsort

                  if (sort_w(j) .gt. sort_w(i)) then

                     temp_w = sort_w(i)
                     sort_w(i) = sort_w(j)
                     sort_w(j) = temp_w

                     temp_i = sort_node(i)
                     sort_node(i) = sort_node(j)
                     sort_node(j) = temp_i

                     temp_i = sort_c(i)
                     sort_c(i) = sort_c(j)
                     sort_c(j) = temp_i

                  endif

               enddo

            enddo

!--------------------------------------------------------------
!           Escribir archivo
!--------------------------------------------------------------

            write(27,*)
            write(27,*) '# ======================================'
            write(27,*) '#     Residue path Participation        '
            write(27,*) '# ======================================'
            write(27,*) '# Chain : ',ichain
            write(27,*) '# Group : ',igrp
            write(27,*) '#'
            write(27,*) '# Residue   Npaths     W_score'

            do i = 1,nsort
               write(27,97) sort_node(i),sort_c(i),sort_w(i)
            enddo

            close(27)

         enddo

      enddo

   97 format(1x,I8,2x,I8,2x,f12.3) 

!==============================================================
! Liberar memoria
!==============================================================

      deallocate(path_c)
      deallocate(path_w)
      deallocate(path_included)


      end subroutine calculate_dijkstra

      end module fmcp_dijkstra
