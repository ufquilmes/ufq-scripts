      module fmcp_interactions
      use fmcp_data
      implicit none

      contains

      subroutine calculate_interactions()

!=======================================================================
!
!     CALCULO DE INTERACCIONES
!
!=======================================================================
!
!     Solamente grupo 2.
!
!     Para cada conformacion: cada atomo pesado de la proteina
!     contra cada atomo pesado del ligando
!
!     Si d < 5 A: ninter(residuo) = ninter(residuo) + 1
!
!=======================================================================

      ninter = 0

      do ifr = 1,nfr(2)


!=======================================================================
!        Todos los atomos pesados de la proteina
!=======================================================================

         do i = 1,nph(ifr,2)

!=======================================================================
!           Indice del residuo al que pertenece el atomo
!=======================================================================

            ii = phres(i,ifr,2)

!=======================================================================
!           Todos los atomos pesados del ligando
!=======================================================================

            do j = 1,nlig

!=======================================================================
!              Diferencias de coordenadas
!=======================================================================

               dx = xph(i,ifr,2) - xlig(j,ifr)
               dy = yph(i,ifr,2) - ylig(j,ifr)
               dz = zph(i,ifr,2) - zlig(j,ifr)

!=======================================================================
!              Distancia al cuadrado
!              equivale a:  d2 < cutoff_inter2 A2
!=======================================================================

               d2 = dx*dx + dy*dy + dz*dz

               if (d2 < cutoff_inter2) then
                  ninter(ii) = ninter(ii) + 1
               endif

            enddo

         enddo

      enddo


!=======================================================================
!
!     Ordenamos los residuos de mayor a menor numero de interacciones.
!
!=======================================================================

      do i = 1,natom
         order(i) = i
      enddo


      do i = 1,natom-1

         max_idx = i

         do j = i+1,natom

            if (ninter(order(j)) > ninter(order(max_idx))) then
               max_idx = j
            endif

         enddo

         if (max_idx /= i) then
            temp = order(i)
            order(i) = order(max_idx)
            order(max_idx) = temp
         endif

      enddo

!=======================================================================
!
!     Numero de residuos que interactuan al menos una vez
!     y promedio interacciones por conformacion
!
!=======================================================================

      nres_interacting = 0

      do i = 1, natom
       pinter(i) = 0.0d0
        if (ninter(i) > 0) then
         nres_interacting = nres_interacting + 1
         ! Guardar el índice del residuo
         interacting_residues(nres_interacting) = i
         pinter(i) = float(ninter(i))/float(nfr(2))
        endif
      enddo

!=======================================================================
!
!     RESULTADOS EN PANTALLA
!
!=======================================================================

      write(*,*)
      write(*,*) '===================================================='
      write(*,*) 'INTERACCIONES PROTEINA - LIGANDO'
      write(*,*) '===================================================='
      write(*,*)

      write(*,'(A19,2X,f6.2)') 'Dist. de corte (A)=', cutoff_inter
      write(*,*) 'Conformaciones analizadas = ',nfr(2)
      write(*,*) 'Atomos pesados de proteina = ',nph_current
      write(*,*) 'Atomos pesados de ligando = ',nlig
      write(*,*) 'Residuos con interacciones = ',nres_interacting
      write(*,*)

      write(*,'(A8,2X,A8,2X,A14)') 'ResPDB','Residuo','Interacc./Conf.'
      write(*,'(A8,2X,A8,2X,A14)') '-------','-------','---------------'


!=======================================================================
!     Escribimos tambien a archivo
!=======================================================================

      open(30,file='protein_ligand_interactions.dat', 
     &                status='replace',action='write')
      write(30,'(A)') '# ResPDB  Residue  Interact./Conf.'

      do i = 1,natom

         ii = order(i)
         if (ninter(ii) > 0) then
            write(30,30) resnum_ca(ii),resname_ca(ii),pinter(ii)
         endif

      enddo

      close(30)

!=======================================================================
!     FIN DE ESTA ETAPA
!=======================================================================

      write(*,*)
      write(*,*) 'Resultados escritos en:'
      write(*,*) 'protein_ligand_interactions.dat'
      write(*,*)
      write(*,*) 'Fin de lectura y calculo de interacciones.'


30    format(I6,2X,A8,4X,f10.4)
      end subroutine calculate_interactions

      end module fmcp_interactions
