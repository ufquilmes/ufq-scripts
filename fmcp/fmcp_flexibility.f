      module fmcp_flexibility
      use fmcp_data
      implicit none

      contains

      subroutine calculate_flexibility()

      !
      ! Calculamos las fluctuaciones de cada carbono Ca,
      ! de cada grupo. para eso, primero va la estructura promedio
      ! 
      !==========================================================

      do igrp = 1, ngrp

         do i = 1, natom

            xmean(i,igrp) = 0.d0
            ymean(i,igrp) = 0.d0
            zmean(i,igrp) = 0.d0

            do ifr = 1, nfr(igrp)
              xmean(i,igrp) = xmean(i,igrp) + xca(i,ifr,igrp)
              ymean(i,igrp) = ymean(i,igrp) + yca(i,ifr,igrp)
              zmean(i,igrp) = zmean(i,igrp) + zca(i,ifr,igrp)
            enddo

            xmean(i,igrp) = xmean(i,igrp) / dble(nfr(igrp))
            ymean(i,igrp) = ymean(i,igrp) / dble(nfr(igrp))
            zmean(i,igrp) = zmean(i,igrp) / dble(nfr(igrp))

            flu(i,igrp) = 0.d0

            do ifr = 1, nfr(igrp)
              dx = xca(i,ifr,igrp) - xmean(i,igrp)
              dy = yca(i,ifr,igrp) - ymean(i,igrp)
              dz = zca(i,ifr,igrp) - zmean(i,igrp)
              flu(i,igrp) = flu(i,igrp) + dx*dx + dy*dy + dz*dz
            enddo

            flu(i,igrp) = dsqrt(flu(i,igrp)/dble(nfr(igrp)))

         enddo

      enddo


      !==================================================
      !
      ! Compute flexibility change
      !
      !==================================================

       open(24,file='rmsf.dat')
        write(24,*) "# Residues   Apo-Form   Holo-Form"
       do i = 1, natom
        drmsf(i) = flu(i,2) - flu(i,1)
        write(24,'(i6,2(3x,f12.6))') i, flu(i,1), flu(i,2)
       enddo
       close(24)

       do i = 1, natom
        do j = i, natom
         fcorr(i,j) = drmsf(i)*drmsf(j)/dble(ngrp)
        enddo
       enddo

      !============================================================
      !
      ! fill lower triangle
      !
      !============================================================

      do i = 1, natom
       do j = i+1, natom
        fcorr(j,i) = fcorr(i,j)
       enddo
      enddo

      !=============================================================
      !
      ! metodologia MCP con distancia de flexibilidad
      !
      !==============================================================

      !========================================================
      ! asigna "distancias de paso" con correlacion de fluctuacion
      !========================================================
       
       do i = 1, natom
          do j = 1, natom

             if (conn(i,j)) then
                  if (dabs(fcorr(i,j)) >= 1.d-12) then
                   dpath(i,j) = dlog( 1.0d0 + (1.0d0/dabs(fcorr(i,j))))
                   else
                   dpath(i,j) = 1.d15
                  endif
             else
                dpath(i,j) = 1.0d15
             endif

          enddo
       enddo    


      end subroutine calculate_flexibility

      end module fmcp_flexibility
