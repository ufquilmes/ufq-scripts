      module fmcp_targets
      use fmcp_data
      implicit none

      contains


      subroutine calculate_targets()

!==============================================================
! Clustering jerarquico independiente para cada cadena.
!
! ngroup = numero de grupos solicitado por el usuario.
!
! Cada grupo final debe tener al menos 7 residuos.
!
! Si nchain_res >= 7*ngroup:
!
!       ngroup_eff = ngroup
!
! Si no:
!
!       ngroup_eff = nchain_res/7
!
! Primero se realiza clustering jerarquico por average linkage
! hasta obtener ngroup_eff clusters.
!
! Luego, si algun cluster tiene menos de 7 residuos, se
! redistribuyen residuos desde clusters con mas de 7 miembros.
!
! IMPORTANTE:
! No se impone ningun tamaño objetivo para los grupos.
! Por lo tanto los grupos pueden ser muy heterogeneos.
!
!==============================================================

      integer :: ngroup_eff
      integer :: nchain_res
      integer :: ncluster
      integer :: old_cluster
      integer :: new_cluster
      integer :: best_cluster
      integer :: best_res
      integer :: donor_cluster
      integer :: small_cluster
      integer :: min_cluster_size
      integer :: nsmall
      integer :: g, i, j, ii, jj, ifr
      integer :: candidate
      integer :: ntemp

      double precision :: best_dist
      double precision :: cluster_dist
      double precision :: dtemp2
      double precision :: min_small_dist

!==============================================================
! Construir nlist1:
!
! Se excluyen:
!   - residuos que interactuan con el ligando
!   - primeros 3 residuos de cada cadena
!   - ultimos 3 residuos de cada cadena
!==============================================================

      ncon = 0

      do ichain = 1,nchains

         do i = chain_start(ichain),chain_end(ichain)

            excluded = .false.

!---------- Primeros tres residuos -----------------------------

            if (i .le. chain_start(ichain)+2) then
               excluded = .true.
            endif

!---------- Ultimos tres residuos ------------------------------

            if (i .ge. chain_end(ichain)-2) then
               excluded = .true.
            endif

!---------- Residuos que interactuan con ligando ---------------

            do j = 1,nres_interacting

               if (i .eq. interacting_residues(j)) then
                  excluded = .true.
               endif

            enddo

            if (.not.excluded) then

               ncon = ncon + 1
               nlist1(ncon) = i

            endif

         enddo

      enddo

!==============================================================
! Numero de grupos solicitado
!==============================================================

      if (ngroup .lt. 1) then

         write(*,*) 'ERROR: ngroup debe ser >= 1'
         stop

      endif

!==============================================================
! Calcular davg y conn.
!
! ESTA PARTE SE MANTIENE IGUAL.
!==============================================================

      conn = .false.
      igrp = 1

      do i = 1,natom-1

         do j = i+1,natom

            ii = 0
            davg(i,j) = 0.0d0

            do ifr = 1,nfr(igrp)

               dx = xca(i,ifr,igrp) - xca(j,ifr,igrp)
               dy = yca(i,ifr,igrp) - yca(j,ifr,igrp)
               dz = zca(i,ifr,igrp) - zca(j,ifr,igrp)

               tempd = dsqrt(dx*dx + dy*dy + dz*dz)

               if (tempd .le. 10.d0) then
                  ii = ii + 1
               endif

               davg(i,j) = davg(i,j) + tempd

            enddo

            davg(i,j) = davg(i,j) / dble(nfr(igrp))
            davg(j,i) = davg(i,j)

            if (dble(ii)/dble(nfr(igrp)) .ge. 0.75d0) then

               conn(i,j) = .true.
               conn(j,i) = .true.

            endif

         enddo

      enddo


!==============================================================
! Inicializar informacion de cadena/grupo.
!==============================================================

      do i = 1,natom
         group_residue(i) = 0
         chain_residue(i) = 0
      enddo

!==============================================================
! Abrir groups.dat UNA SOLA VEZ
!==============================================================

      open(26,file='groups.dat',status='replace')

!==============================================================
! Clustering independiente para cada cadena
!==============================================================




      do ichain = 1,nchains

!---------- Construir lista de residuos de la cadena ------------

         nchain_res = 0

         do i = chain_start(ichain),chain_end(ichain)

            excluded = .false.

            if (i .le. chain_start(ichain)+2) then
               excluded = .true.
            endif

            if (i .ge. chain_end(ichain)-2) then
               excluded = .true.
            endif

            do j = 1,nres_interacting

               if (i .eq. interacting_residues(j)) then
                  excluded = .true.
               endif

            enddo

            if (.not.excluded) then

               nchain_res = nchain_res + 1
               chain_list(nchain_res) = i

               chain_residue(i) = ichain

            endif

         enddo

!==============================================================
! Determinar numero efectivo de grupos.
!==============================================================

         ngroup_eff = nchain_res / 7

         if (ngroup_eff .gt. ngroup) then
            ngroup_eff = ngroup
         endif

         if (ngroup_eff .lt. 1) then

            write(*,*) ' '
            write(*,*) 'Cadena ',ichain
            write(*,*) 'No hay suficientes residuos para'
            write(*,*) 'formar un grupo de al menos 7 residuos.'
            cycle

         endif

!==============================================================
! Informacion
!==============================================================

         write(*,*) ' '
         write(*,*) '----------------------------------------------'
         write(*,*) 'Cadena            ',ichain
         write(*,*) 'Inicio            = ',chain_start(ichain)
         write(*,*) 'Fin               = ',chain_end(ichain)
         write(*,*) 'Residuos incluidos = ',nchain_res
         write(*,*) 'Grupos solicitados = ',ngroup
         write(*,*) 'Grupos efectivos   = ',ngroup_eff
         write(*,*) '----------------------------------------------'

!==============================================================
! Inicializar cada residuo como un cluster individual.
!
! cluster(global_residue) contiene el identificador del cluster.
!==============================================================

         do i = 1,nchain_res

            cluster(chain_list(i)) = i
            cluster_new(i) = i

         enddo

         ncluster = nchain_res

!==============================================================
! CLUSTERING JERARQUICO
!
! Average linkage:
!
! distancia entre cluster A y B =
!
!     promedio de davg(i,j)
!
! para todos los residuos i de A y j de B.
!
! Se fusionan los dos clusters mas cercanos hasta obtener
! exactamente ngroup_eff clusters.
!==============================================================

         do while (ncluster .gt. ngroup_eff)

            best_dist = 1.0d30
            best_cluster = 0
            candidate = 0
            old_cluster = 0
            new_cluster = 0


!---------- Buscar los dos clusters mas cercanos ---------------

            do i = 1,nchain_res-1

               old_cluster = cluster(chain_list(i))

               do j = i+1,nchain_res

                  new_cluster = cluster(chain_list(j))

                  if (old_cluster .eq. new_cluster) cycle

!------------------ Average linkage ----------------------------

                  dtemp2 = 0.0d0
                  ntemp = 0

                  do ii = 1,nchain_res

                   if (cluster(chain_list(ii)).ne.old_cluster) cycle

                     do jj = 1,nchain_res

                      if (cluster(chain_list(jj)).ne.new_cluster) cycle

                       dtemp2=dtemp2+davg(chain_list(ii),chain_list(jj))
                       ntemp = ntemp + 1

                     enddo

                  enddo

                  if (ntemp .gt. 0) then

                     cluster_dist = dtemp2 / dble(ntemp)

                     if (cluster_dist .lt. best_dist) then

                        best_dist = cluster_dist
                        best_cluster = old_cluster
                        candidate = new_cluster

                     endif

                  endif

               enddo

            enddo

!---------- Fusionar los dos clusters --------------------------

            do i = 1,nchain_res

               if (cluster(chain_list(i)) .eq. candidate) then
                  cluster(chain_list(i)) = best_cluster
               endif

            enddo

            ncluster = ncluster - 1

         enddo

!==============================================================
! Renumerar los clusters finales consecutivamente:
!
! 1, 2, 3, ..., ngroup_eff
!==============================================================

         new_cluster = 0

         do i = 1,nchain_res

            old_cluster = cluster(chain_list(i))

            if (old_cluster .eq. 0) cycle

            ntemp = 0

            do j = 1,i-1

               if (cluster(chain_list(j))
     &             .eq. old_cluster) then

                  ntemp = 1
                  exit

               endif

            enddo

            if (ntemp .eq. 0) then

               new_cluster = new_cluster + 1

               do j = 1,nchain_res

                  if (cluster(chain_list(j))
     &                .eq. old_cluster) then

                     cluster(chain_list(j)) = -new_cluster

                  endif

               enddo

            endif

         enddo

         do i = 1,nchain_res

            if (cluster(chain_list(i)) .lt. 0) then

               cluster(chain_list(i)) =
     &              -cluster(chain_list(i))

            endif

         enddo

!==============================================================
! Contar miembros de cada grupo.
!==============================================================

         do g = 1,ngroup_eff
            nmember(g) = 0
         enddo

         do i = 1,nchain_res

            g = cluster(chain_list(i))

            if (g .ge. 1 .and. g .le. ngroup_eff) then
               nmember(g) = nmember(g) + 1
            endif

         enddo

!==============================================================
! REPARACION DE GRUPOS PEQUENOS
!
! Si un grupo tiene menos de 7 residuos:
!
!   - se buscan residuos en grupos con >7 miembros
!   - se mueve el residuo cuya distancia promedio al grupo
!     pequeno sea minima.
!
! No se imponen tamaños iguales.
!==============================================================

         do

            small_cluster = 0
            min_cluster_size = 1000000

!---------- Buscar el grupo mas pequeño ------------------------

            do g = 1,ngroup_eff

               if (nmember(g) .lt. 7) then

                  if (nmember(g) .lt. min_cluster_size) then

                     min_cluster_size = nmember(g)
                     small_cluster = g

                  endif

               endif

            enddo

!---------- Ya no quedan grupos pequeños ----------------------

            if (small_cluster .eq. 0) exit

!---------- Buscar el mejor residuo donante --------------------

            min_small_dist = 1.0d30
            best_res = 0
            donor_cluster = 0

            do i = 1,nchain_res

               if (cluster(chain_list(i))
     &             .eq. small_cluster) cycle

               donor_cluster = cluster(chain_list(i))

               if (nmember(donor_cluster) .le. 7) cycle

!-------------- Distancia promedio al grupo pequeño -----------

               dtemp2 = 0.0d0
               ntemp = 0

               do j = 1,nchain_res

                  if (cluster(chain_list(j)).eq. small_cluster) then

                     dtemp2 = dtemp2 + davg(chain_list(i),
     &                             chain_list(j))

                     ntemp = ntemp + 1

                  endif

               enddo

               if (ntemp .gt. 0) then

                  dtemp2 = dtemp2 / dble(ntemp)

                  if (dtemp2 .lt. min_small_dist) then
                     min_small_dist = dtemp2
                     best_res = chain_list(i)
                  endif

               endif

            enddo

!---------- Seguridad ------------------------------------------

            if (best_res .eq. 0) then

               write(*,*) ' '
               write(*,*) 'ERROR: no se pudo reparar grupo'
               write(*,*) 'pequeno en cadena ',ichain
               write(*,*) 'Grupo = ',small_cluster
               stop

            endif

!---------- Identificar cluster donante ------------------------

            donor_cluster = cluster(best_res)

!---------- Mover residuo --------------------------------------

            cluster(best_res) = small_cluster
            nmember(donor_cluster) = nmember(donor_cluster) - 1
            nmember(small_cluster) = nmember(small_cluster) + 1

         enddo

!==============================================================
! Verificacion final
!==============================================================

         nsmall = 0

         do g = 1,ngroup_eff

            if (nmember(g) .lt. 7) then
               nsmall = nsmall + 1
            endif

         enddo

         if (nsmall .ne. 0) then
            write(*,*) ' '
            write(*,*) 'ERROR: quedaron grupos con menos de 7'
            write(*,*) 'residuos en cadena ',ichain
            stop
         endif

!==============================================================
! Guardar grupo y cadena para cada residuo global.
!
! IMPORTANTE:
!
! group_residue(residuo_global) = grupo
! chain_residue(residuo_global) = cadena
!
! Esto es lo que necesita Dijkstra.
!==============================================================

         do i = 1,nchain_res

            j = chain_list(i)

            group_residue(j) = cluster(j)
            chain_residue(j) = ichain

         enddo

!==============================================================
! Imprimir grupos
!==============================================================

         write(*,*) '----------------------------------------------'
         write(*,*) 'Grupos de la cadena ',ichain
         write(*,*) '----------------------------------------------'

         do g = 1,ngroup_eff

            write(*,100) g,nmember(g)

 100        format(' Group ',i2,': ',i7,' Residues')

         enddo


         write(26,*) '----------------------------------------------'
         write(26,*) '# Cadena ',ichain
         write(26,*) '# Inicio ',chain_start(ichain)
         write(26,*) '# Fin    ',chain_end(ichain)
         write(26,*) '# Grupos ',ngroup_eff
         write(26,*) '#'
         write(26,*) '----------------------------------------------'

         do g = 1,ngroup_eff

            write(26,*) ' '
            write(26,*) '# GRUPO ',g
            write(26,*) '----------------------------------------------'

            do i = 1,nchain_res

               if (cluster(chain_list(i)) .eq. g) then
                  write(26,'(i6)',advance='no') chain_list(i)
               endif

            enddo

            write(26,*)

         enddo

         write(26,*) ' '
         write(26,*) '----------------------------------------------'
         write(26,*) '----------------------------------------------'
         write(26,*) '# Group Number    Number of residues          '
         write(26,*) '----------------------------------------------'
         write(26,*) '----------------------------------------------'

         do g = 1,ngroup_eff
            write(26,'(6X,i3,10x,i4)') g,nmember(g)
         enddo

      enddo


!==============================================================
! Cerrar groups.dat después de procesar todas las cadenas
!==============================================================

      close(26)


      return
      end subroutine calculate_targets

      end module fmcp_targets
