      module fmcp_vmd
      use fmcp_data
      implicit none

      contains

      subroutine write_vmd_tcl()

      integer :: node
      integer :: node_start
      integer :: node_end
      integer :: node_prev
      integer :: ichain
      integer :: igrp

      real(8) :: score_min
      real(8) :: score_max
      real(8) :: score_norm

!==============================================================
! Crear script TCL para VMD
!==============================================================

      open(28,file='score.groups.tcl',status='replace')

!--------------------------------------------------------------
! Eliminar representaciones existentes
!--------------------------------------------------------------

      write(28,*) 'set nrep [molinfo top get numreps]'
      write(28,200)
      write(28,*) '    mol delrep $rep top'
      write(28,*) '}'

!==============================================================
! Construir selección de residuos pertenecientes a grupos
!==============================================================

      write(28,*)
      write(28,*) '# ==============================================='
      write(28,*) '# Grouped residues'
      write(28,*) '# ==============================================='

      write(28,*) 'set grouped_resids "'

      do node = 1,natom

         if (group_residue(node) .gt. 0) then
            write(28,103) node
         endif

      enddo

      write(28,*) '"'

!==============================================================
! Residuos no pertenecientes a grupos
!==============================================================

      write(28,*)
      write(28,*) '# ==============================================='
      write(28,*) '# Non-target residues'
      write(28,*) '# ==============================================='

!--------------------------------------------------------------
! Interior claro y transparente de los no-target
!--------------------------------------------------------------

      write(28,*) 'mol representation Ribbons 0.47 10 10'
      write(28,*) 'mol color ColorID 8'
      write(28,*) 'mol selection "protein and not'
      write(28,*) '     (resid $grouped_resids)"'
      write(28,*) 'mol material Transparent'
      write(28,*) 'mol addrep top'
      write(28,*) 'material change opacity Transparent 0.02'

!--------------------------------------------------------------
! Borde negro transparente de los no-target
!--------------------------------------------------------------

      write(28,*) 'mol representation Ribbons 0.57 10 10'
      write(28,*) 'mol color ColorID 16'
      write(28,*) 'mol selection "protein and not'
      write(28,*) '     (resid $grouped_resids)"'
      write(28,*) 'mol material BlownGlass'
      write(28,*) 'mol addrep top'
      write(28,*) 'material change opacity Transparent 0.08'

!==============================================================
! Información de scores
!==============================================================

      write(28,*)
      write(28,*) '# ============================================='
      write(28,*) '# FMCP group scores'
      write(28,*) '# ============================================='

      do ichain = 1,nchains

         score_min = 1.0d30
         score_max = -1.0d30

         do igrp = 1,ngroup

            if (ngrp_member(ichain,igrp) .gt. 0) then

               score_min = min(score_min,
     &                         grp_score(ichain,igrp))

               score_max = max(score_max,
     &                         grp_score(ichain,igrp))

            endif

         enddo

         write(28,*)
         write(28,*) '# -----------------------------------------------'
         write(28,100) ichain,score_min,score_max
         write(28,*) '# -----------------------------------------------'

         do igrp = 1,ngroup

            if (ngrp_member(ichain,igrp) .gt. 0) then

               if (score_max .gt. score_min) then

                  score_norm =
     &            (grp_score(ichain,igrp)-score_min) /
     &            (score_max-score_min)

               else

                  score_norm = 0.5d0

               endif

               score_norm = max(0.0d0,
     &                           min(1.0d0,score_norm))

               write(28,101) ichain,igrp,
     &                      grp_score(ichain,igrp),
     &                      score_norm

            endif

         enddo

      enddo

!==============================================================
! Representaciones de los grupos FMCP
!==============================================================

      write(28,*)
      write(28,*) '# ==============================================='
      write(28,*) '# FMCP groups'
      write(28,*) '# Ribbons + Beta'
      write(28,*) '# Blue -> White -> Red'
      write(28,*) '# ==============================================='

      do ichain = 1,nchains

!--------------------------------------------------------------
! Normalización independiente por cadena
!--------------------------------------------------------------

         score_min = 1.0d30
         score_max = -1.0d30

         do igrp = 1,ngroup

            if (ngrp_member(ichain,igrp) .gt. 0) then

               score_min = min(score_min,
     &                         grp_score(ichain,igrp))

               score_max = max(score_max,
     &                         grp_score(ichain,igrp))

            endif

         enddo

!--------------------------------------------------------------
! Recorrer grupos
!--------------------------------------------------------------

         do igrp = 1,ngroup

            if (ngrp_member(ichain,igrp) .gt. 0) then

               if (score_max .gt. score_min) then

                  score_norm =
     &            (grp_score(ichain,igrp)-score_min) /
     &            (score_max-score_min)

               else

                  score_norm = 0.5d0

               endif

               score_norm = max(0.0d0,
     &                           min(1.0d0,score_norm))

               write(28,*)
               write(28,102) ichain,igrp,
     &                      grp_score(ichain,igrp),
     &                      score_norm

!--------------------------------------------------------------
! Buscar segmentos CONTIGUOS del grupo
!--------------------------------------------------------------

               do node_start = 1,natom

                  if (chain_residue(node_start) .eq. ichain
     &                .and.
     &                group_residue(node_start) .eq. igrp) then

!--------------------------------------------------------------
! Determinar si node_start es el comienzo de un segmento
!--------------------------------------------------------------

                     if (node_start .eq. 1) then

                        node_prev = 0

                     else

                        node_prev = node_start - 1

                     endif

                     if (node_prev .eq. 0) then

                        node_end = node_start

                     else if
     &                  (chain_residue(node_prev) .eq. ichain
     &                  .and.
     &                  group_residue(node_prev) .eq. igrp
     &                  .and.
     &                  node_prev .eq. node_start-1) then

                        cycle

                     else

                        node_end = node_start

                     endif

!--------------------------------------------------------------
! Extender hasta el final del segmento continuo
!--------------------------------------------------------------

                     do node = node_start+1,natom

                        if (chain_residue(node) .eq. ichain
     &                      .and.
     &                      group_residue(node) .eq. igrp
     &                      .and.
     &                      node .eq. node_end+1) then

                           node_end = node

                        else

                           exit

                        endif

                     enddo

                     write(28,*) '# Segment ',node_start,
     &                           ' - ',node_end

!==============================================================
! Borde negro del grupo
!==============================================================

                     write(28,*)
                     write(28,*) '# FMCP outline'

                     write(28,*) 'mol representation ',
     &                           'Ribbons 0.57 10 10'

                     write(28,*) 'mol color ColorID 16'

                     write(28,*) 'mol selection {resid \'

                     do node = node_start,node_end

                        write(28,'(I0,1X)',advance='no') node

                     enddo

                     write(28,*) '}'

                     write(28,*) 'mol material Opaque'
                     write(28,*) 'mol addrep top'

!==============================================================
! Grupo FMCP coloreado
!==============================================================

                     write(28,*)
                     write(28,*) '# FMCP colored group'

                     write(28,*) 'set sel [atomselect top "resid \'

                     do node = node_start,node_end

                        write(28,'(I0,1X)',advance='no') node

                     enddo

                     write(28,*) '"]'

!--------------------------------------------------------------
! Invertir Beta:
! score bajo -> azul
! score alto -> rojo
!--------------------------------------------------------------

                     write(28,104) 1.0d0-score_norm

                     write(28,*) '$sel delete'

                     write(28,*) 'mol representation ',
     &                           'Ribbons 0.47 10 10'

                     write(28,*) 'mol color Beta'

                     write(28,*) 'mol selection {resid \'

                     do node = node_start,node_end

                        write(28,'(I0,1X)',advance='no') node

                     enddo

                     write(28,*) '}'

                     write(28,*) 'mol material Opaque'
                     write(28,*) 'mol addrep top'

              write(28,*) 'set rep [expr {[molinfo top get numreps]-1}]'

                     write(28,*) 'mol scaleminmax top $rep 0.0 1.0'

                  endif

               enddo

            endif

         enddo

      enddo

!==============================================================
! Ligando
!==============================================================

      write(28,*)
      write(28,*) '# ============================================='
      write(28,*) '# Ligand'
      write(28,*) '# ============================================='

      write(28,*) 'mol representation VDW'
      write(28,*) 'mol material EdgyShiny'
      write(28,*) 'mol color Element'
      write(28,*) 'mol selection {resname ',
     &             trim(lig_resname),'}'
      write(28,*) 'mol addrep top'

      close(28)

!==============================================================
! Formatos
!==============================================================

  100 format('# Chain ',I4,
     &       ' score_min = ',F12.5,
     &       ' score_max = ',F12.5)

  101 format('# Chain ',I4,
     &       ' Group ',I4,
     &       ' Score = ',F12.5,
     &       ' Normalized = ',F8.5)

  102 format('# Chain ',I4,
     &       ' Group ',I4,
     &       ' Score = ',F12.5,
     &       ' Normalized = ',F8.5)

  103 format(I8,1X)

  104 format('$sel set beta ',F12.8)

  200 format('for {set rep [expr {$nrep-1}]} ',
     &       '{$rep >= 0} {incr rep -1} {')

      end subroutine write_vmd_tcl

      end module fmcp_vmd
