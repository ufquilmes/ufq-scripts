      module fmcp_read
      use fmcp_data
      implicit none

      contains
      subroutine read_trajectories()

!=======================================================================
!
!     Lee los nombres de las dos trayectorias y el nombre del ligando
!
!=======================================================================

      givename = iargc()

      if (givename < 4) then
         write(*,*)
         write(*,*) 'Uso:'
         write(*,*) 'fmcp traj1 traj2 ligando ngroup'
         write(*,*)
         write(*,*) 'Ejemplo:'
         write(*,*) 'fmcp apo.pdb complejo.pdb 1PR 5'
         write(*,*)
         stop
      endif

!=======================================================================
!
!     Nombre del ligando
!     Ejemplo: 1PR, ATP
!
!=======================================================================

      call getarg(3,nombre)
      read(nombre,'(A3)') lig_resname
      lig_resname = adjustl(lig_resname)

      write(*,*)
      write(*,*) 'Ligando buscado: ',trim(lig_resname)
      write(*,*)

!=======================================================================
!
!     Numero de grupos por cadena  Ejemplo: 5 
!
!=======================================================================

      call getarg(4,nombre)
      read(nombre,*) ngroup

      if (ngroup < 1) then
         write(*,*)
         write(*,*) 'ERROR: ngroup debe ser >= 1'
         write(*,*)
         stop
      endif

      write(*,*) 'Numero de grupos = ',ngroup
      write(*,*)

!=======================================================================
!
!     Inicializamos
!
!=======================================================================

      natom = 0
      nlig  = 0
      ninter = 0
      nfr = 0
      nph = 0
      nchains = 0
      chain_start = 0
      chain_end = 0

!=======================================================================
!
!     Leemos las dos trayectorias
!
!     Grupo 1:
!
!         C-alpha + atomos pesados de la proteina
!
!     Grupo 2:
!
!         C-alpha + atomos pesados de la proteina
!         + atomos pesados del ligando
!
!=======================================================================

      do igrp = 1, ngrp

!=======================================================================
!
!        Nombre del archivo
!
!=======================================================================

         call getarg(igrp,nombre)

         write(*,*)
         write(*,*) 'Leyendo ',trim(nombre)
         open(22,file=trim(nombre),status='old')
         nfr(igrp) = 0

   10    continue


!=======================================================================
!
!        Comenzamos una nueva conformacion
!
!=======================================================================

         nfr(igrp) = nfr(igrp) + 1
         iatom = 0
         nph_frame = 0
         nlig_frame = 0

         if (igrp == 1 .and. nfr(igrp) == 1) then
           new_chain = .true.
         else
           new_chain = .false.
         endif

!=======================================================================
!
!        Variables para determinar el residuo actual
!
!=======================================================================

         current_res = 0
         current_resnum = -999999
           
         if (igrp == 1 .and. nfr(igrp) == 1) then
           nchains = 0
           chain_start = 0
           chain_end = 0
         endif
!=======================================================================
!
!        Leemos la conformacion
!
!=======================================================================

         do

            read(22,'(A)',end=100) line

!=======================================================================
!           Fin de la conformacion
!=======================================================================

            if (line(1:6) == 'ENDMDL') exit


!=======================================================================
!
!           Fin de una cadena
!
!=======================================================================

           if (line(1:3) == 'TER') then

             if (nchains > 0) then
               chain_end(nchains) = current_res
             endif
   
             if (igrp == 1 .and. nfr(igrp) == 1) then
              new_chain = .true.
             endif
             cycle

           endif

!=======================================================================
!
!           Solo nos interesan lineas ATOM
!
!=======================================================================

            if (line(1:4) /= 'ATOM') cycle

!=======================================================================
!           Nombre del atomo Columnas 13-16
!=======================================================================

            read(line(13:16),'(A4)') atnm_tmp
            atnm_tmp = adjustl(atnm_tmp)


!=======================================================================
!           Nombre del residuo Columnas 18-20
!=======================================================================

            read(line(18:20),'(A3)') resname_tmp
            resname_tmp = adjustl(resname_tmp)

!=======================================================================
!           Numero del residuo Columnas 23-26
!=======================================================================

            read(line(23:26),'(I4)') resnum_tmp

!=======================================================================
!
!           Elemento quimico Columnas 77-78
!
!=======================================================================

            element_tmp = '  '
            if (len_trim(line) >= 78) then

               read(line(77:78),'(A2)') element_tmp

            endif
            element_tmp = adjustl(element_tmp)

!=======================================================================
!
!           PROTEINA Si el residuo NO es el ligando.
!
!=======================================================================

            if (trim(resname_tmp) /= trim(lig_resname)) then

!=======================================================================
!
!              Detectamos cambio de residuo.
!
!              current_res es el indice interno del residuo,
!              que coincide con el orden de los C-alpha.
!              Tambien tenemos el numero de cadenas
!
!=======================================================================

            if (resnum_tmp /= current_resnum .or. new_chain) then

               if (new_chain) then

                  if (nchains > 0) then
                     chain_end(nchains) = current_res
                  endif

                  nchains = nchains + 1

                  if (nchains > max_chain) then
                     write(*,*)
                     write(*,*) 'ERROR: #cadenas > ',max_chain
                     stop
                  endif

                  chain_start(nchains) = current_res + 1

                  new_chain = .false.

               endif

               current_resnum = resnum_tmp
               current_res = current_res + 1

            endif


!=======================================================================
!              C-alpha
!=======================================================================

               if (trim(atnm_tmp) == 'CA') then

                  iatom = iatom + 1

                  if (iatom > nmax) then

                     write(*,*)
                     write(*,*) 'ERROR: demasiados C-alpha'
                     write(*,*) 'nmax = ',nmax
                     stop

                  endif

!                 Coordenadas

                  read(line(31:38),*) xca(iatom,nfr(igrp),igrp)
                  read(line(39:46),*) yca(iatom,nfr(igrp),igrp)
                  read(line(47:54),*) zca(iatom,nfr(igrp),igrp)

!                 Informacion del residuo

                  resnum_ca(iatom) = resnum_tmp
                  resname_ca(iatom) = resname_tmp

               endif

!=======================================================================
!
!              TODOS LOS ATOMOS PESADOS DE LA PROTEINA
!
!              Incluye los C-alpha.  Todo lo que no sea H.
!
!=======================================================================

               if (trim(element_tmp) /= 'H') then

                  nph_frame = nph_frame + 1

                  if (nph_frame > nph_max) then

                     write(*,*)
                     write(*,*) 'ERROR: demasiados atomos pesados'
                     write(*,*) 'de la proteina'
                     write(*,*) 'nph_max = ',nph_max
                     stop

                  endif

!                 Coordenadas

                  read(line(31:38),*) xph(nph_frame,nfr(igrp),igrp)
                  read(line(39:46),*) yph(nph_frame,nfr(igrp),igrp)
                  read(line(47:54),*) zph(nph_frame,nfr(igrp),igrp)


!                 Guardamos el indice del residuo

                  phres(nph_frame,nfr(igrp),igrp) = current_res


               endif


!=======================================================================
!
!           LIGANDO Solamente se busca en el grupo 2.
!
!=======================================================================

            else

               if (igrp == 2) then

!=======================================================================
!
!                 Solamente atomos pesados
!
!=======================================================================

                  if (trim(element_tmp) /= 'H') then

                     nlig_frame = nlig_frame + 1

                     if (nlig_frame > nlig_max) then

                        write(*,*)
                        write(*,*) 'ERROR: too many atoms in ligand'
                        write(*,*) 'nlig_max = ',nlig_max
                        stop

                     endif

!                    Coordenadas X,Y,Z

                     read(line(31:38),*) xlig(nlig_frame,nfr(igrp))
                     read(line(39:46),*) ylig(nlig_frame,nfr(igrp))
                     read(line(47:54),*) zlig(nlig_frame,nfr(igrp))

                  endif

               endif

            endif

         enddo


!=======================================================================
!
!        Terminamos de leer una conformacion
!
!=======================================================================

         last_iatom = iatom

         if (igrp == 1 .and. nfr(igrp) == 1) then

             if (nchains > 0) then
                chain_end(nchains) = current_res
             endif

          endif

!=======================================================================
!
!        Guardamos el numero de atomos pesados de la proteina
!        de esta conformacion.
!
!=======================================================================

         nph(nfr(igrp),igrp) = nph_frame

!=======================================================================
!
!        En el grupo 2:
!
!        Guardamos el numero de atomos pesados del ligando
!        encontrado en la primera conformacion.
!
!=======================================================================

         if (igrp == 2) then

            if (nfr(igrp) == 1) then

               nlig = nlig_frame

               write(*,*)
               write(*,*) 'Atomos pesados del ligando = ',nlig
               write(*,*)


            else


               if (nlig_frame /= nlig) then

                  write(*,*)
                  write(*,*) 'ERROR: distinto numero de atomos pesados'
                  write(*,*) 'del ligando entre conformaciones'
                  write(*,*)
                  write(*,*) 'esperados  = ',nlig
                  write(*,*) 'encontrados = ',nlig_frame
                  write(*,*)

                  stop

               endif

            endif

         endif


!=======================================================================
!
!        Leemos la siguiente conformacion
!
!=======================================================================

         goto 10

!=======================================================================
!
!        Llegamos a EOF
!
!=======================================================================

  100    continue

!=======================================================================
!
!        Eliminamos la ultima entrada porque no termino
!        mediante ENDMDL.
!
!=======================================================================

         nfr(igrp) = nfr(igrp) - 1
         nat_grp = last_iatom

!=======================================================================
!
!        Verificamos que ambos grupos tengan exactamente
!        el mismo numero de C-alpha.
!
!=======================================================================

         if (igrp == 1) then

            natom = nat_grp

         else

            if (nat_grp /= natom) then

               write(*,*)
               write(*,*) 'ERROR: distinto numero de CA'
               write(*,*) 'grupo   = ',igrp
               write(*,*) 'natom   = ',natom
               write(*,*) 'nat_grp = ',nat_grp
               stop

            endif

         endif

         close(22)


      enddo


!=======================================================================
!
!     Verificamos el numero de atomos pesados de la proteina
!
!     Debe ser el mismo en todas las conformaciones del grupo 2.
!
!=======================================================================

      nph_current = nph(1,2)

      do ifr = 2,nfr(2)

         if (nph(ifr,2) /= nph_current) then

            write(*,*)
            write(*,*) 'ERROR: distinto numero de atomos pesados'
            write(*,*) 'de la proteina entre conformaciones'
            write(*,*)
            write(*,*) 'Conformacion 1 = ',nph_current
            write(*,*) 'Conformacion ',ifr,' = ',nph(ifr,2)
            write(*,*)
            stop

         endif

      enddo


!=======================================================================
!
!     INFORMACION GENERAL
!
!=======================================================================

      write(*,*)
      write(*,*) '=============================================='
      write(*,*) 'Resumen de lectura'
      write(*,*) '=============================================='
      write(*,*)


      write(*,*) 'Numero de cadenas = ',nchains
      write(*,*)


      write(*,*) 'Numero de residuos / CA = ',natom
      write(*,*) 'Conformaciones grupo 1 = ',nfr(1)
      write(*,*) 'Conformaciones grupo 2 = ',nfr(2)
      write(*,*) 'Atomos pesados proteina = ',nph_current
      write(*,*) 'Atomos pesados ligando = ',nlig
      write(*,*)

      end subroutine read_trajectories

      end module fmcp_read
