      module fmcp_data
      implicit none

!=======================================================================
!  FMCP - VARIABLES GLOBALES
!
!  Organizacion:
!    1. Parametros y dimensiones
!    2. Contadores e indices generales
!    3. Informacion de cadenas y residuos
!    4. Coordenadas de la proteina (C-alpha)
!    5. Coordenadas de atomos pesados de la proteina
!    6. Coordenadas del ligando
!    7. Interacciones proteina-ligando
!    8. Seleccion de targets y conectividad
!    9. Flexibilidad y cambio de flexibilidad
!   10. Grafo y caminos minimos (Dijkstra)
!   11. Variables auxiliares y resultados
!   12. Variables de texto y lectura PDB
!
!  IMPORTANTE:
!  No se modificaron nombres, tipos ni dimensiones respecto de la
!  version modular anterior. Solamente se reorganizaron las declaraciones.
!=======================================================================

!=======================================================================
! 1. PARAMETROS Y DIMENSIONES MAXIMAS
!=======================================================================

      integer, parameter :: nmax = 10000
      integer, parameter :: nfr_max = 10000
      integer, parameter :: ngrp = 2
      integer, parameter :: nlig_max = 500
      integer, parameter :: max_chain = 6
      integer, parameter :: ngroup_max=30

!     Numero maximo de atomos pesados de la proteina
      integer, parameter :: nph_max = 30000

!     Distancia de corte para interacciones proteina-ligando
      real(8), parameter :: cutoff_inter = 7.5d0
      real(8), parameter :: cutoff_inter2 = cutoff_inter**2

!=======================================================================
! 2. CONTADORES E INDICES GENERALES
!=======================================================================

!     Indices generales de bucles y calculos
      integer :: i, ii, iii, j, jj, k, u, v

!     Indices y contadores de caminos / Dijkstra
      integer :: min_idx, current, ntop, leg
      integer :: js, jt, path_len

!     Numero de residuos / C-alpha
      integer :: natom

!     Contadores generales
      integer :: temp, ncon, givename
      integer :: ifr, igrp, iatom, nat_grp, last_iatom, ntot
      integer, dimension(ngrp) :: nfr

!     Conformaciones y atomos del ligando
      integer :: nlig
      integer :: nlig_frame

!     Cadenas
      integer :: nchains
      integer :: chain_tmp_idx
      integer, dimension(max_chain) :: chain_start
      integer, dimension(max_chain) :: chain_end

!=======================================================================
! 3. INFORMACION DE RESIDUOS, CADENAS Y ATOMOS PESADOS
!=======================================================================

!     Variables utilizadas durante la lectura de cada conformacion
      integer :: nph_frame
      integer :: current_res
      integer :: current_resnum
      integer :: resnum_tmp
      integer :: nph_current

!     Informacion de cada C-alpha
!     ResPDB y nombre del residuo asociado a cada residuo interno
      integer, dimension(nmax) :: resnum_ca
      character(len=3), dimension(nmax) :: resname_ca

!     Numero de atomos pesados de la proteina por conformacion
!     nph(ifr,igrp)
      integer, dimension(nfr_max,ngrp) :: nph

!     Residuo al que pertenece cada atomo pesado
!     phres(i,ifr,igrp) = indice interno del residuo/C-alpha
      integer, dimension(nph_max,nfr_max,ngrp) :: phres

!=======================================================================
! 4. COORDENADAS DE LOS C-alpha DE LA PROTEINA
!=======================================================================

!     xca(i,ifr,igrp), yca(i,ifr,igrp), zca(i,ifr,igrp)
!       i    = residuo / C-alpha
!       ifr  = conformacion
!       igrp = grupo
!
!     Grupo 1 = apo
!     Grupo 2 = holo / complejo
      real(8), dimension(nmax,nfr_max,ngrp) :: xca
      real(8), dimension(nmax,nfr_max,ngrp) :: yca
      real(8), dimension(nmax,nfr_max,ngrp) :: zca

!=======================================================================
! 5. COORDENADAS DE TODOS LOS ATOMOS PESADOS DE LA PROTEINA
!=======================================================================

!     Incluye los C-alpha.
!     xph(i,ifr,igrp), yph(i,ifr,igrp), zph(i,ifr,igrp)
      real(8), dimension(nph_max,nfr_max,ngrp) :: xph
      real(8), dimension(nph_max,nfr_max,ngrp) :: yph
      real(8), dimension(nph_max,nfr_max,ngrp) :: zph

!=======================================================================
! 6. COORDENADAS DE LOS ATOMOS PESADOS DEL LIGANDO
!=======================================================================

!     El ligando se encuentra solamente en el grupo 2.
!     xlig(i,ifr), ylig(i,ifr), zlig(i,ifr)
      real(8), dimension(nlig_max,nfr_max) :: xlig
      real(8), dimension(nlig_max,nfr_max) :: ylig
      real(8), dimension(nlig_max,nfr_max) :: zlig

!=======================================================================
! 7. INTERACCIONES PROTEINA-LIGANDO
!=======================================================================

!     Numero de interacciones de cada residuo
      integer, dimension(nmax) :: ninter

!     Promedio de interacciones por conformacion
      real(8), dimension(nmax) :: pinter

!     Orden de los residuos de mayor a menor numero de interacciones
      integer, dimension(nmax) :: order

!     Lista de residuos que interactuan con el ligando
      integer, dimension(nmax) :: interacting_residues
      integer :: nres_interacting

!     Variables para calcular distancias proteina-ligando
      real(8) :: d2

!=======================================================================
! 8. SELECCION DE TARGETS Y CONECTIVIDAD DEL GRAFO
!=======================================================================

!     Residuos que pueden actuar como targets
      integer, dimension(nmax) :: nlist1

!     Variables para excluir residuos
      logical :: excluded

!     Conectividad del grafo:
!     conn(i,j) = TRUE si i y j son vecinos segun el criterio de distancia
      logical, dimension(nmax,nmax) :: conn

!     Matriz de proximidad (reservada para el calculo correspondiente)
      logical, dimension(nmax,nmax) :: prox

!     Listas auxiliares relacionadas con seleccion/ordenamiento
      integer :: max_idx

      !para el agrupamiento
      integer, dimension(nmax):: cluster, cluster_new
      integer, dimension(nmax):: nmember, chain_list
      integer :: ncluster, imin, jmin, ichain, nchain_res, ngroup 

      real(8), dimension(nmax,nmax) :: davg
      real(8), dimension(max_chain,ngroup_max) :: grp_score
      real(8) :: dmin, dtemp
 
      integer group_residue(nmax)
      integer chain_residue(nmax)

!=======================================================================
! 9. FLEXIBILIDAD Y CAMBIO DE FLEXIBILIDAD
!=======================================================================

!     Coordenadas medias de cada C-alpha para cada grupo
      real(8), dimension(nmax,ngrp) :: xmean
      real(8), dimension(nmax,ngrp) :: ymean
      real(8), dimension(nmax,ngrp) :: zmean

!     Fluctuacion / flexibilidad de cada residuo en cada grupo
      real(8), dimension(nmax,ngrp) :: flu

!     Cambio de flexibilidad entre los dos estados
      real(8), dimension(nmax) :: drmsf

!     Matriz de correlacion/cambio de flexibilidad
      real(8), dimension(nmax,nmax) :: fcorr

!=======================================================================
! 10. GRAFO, DISTANCIAS DE CAMINO Y DIJKSTRA
!=======================================================================

!     Distancias de paso del grafo
      real(8), dimension(nmax,nmax) :: dpath

!     Distancia minima acumulada desde un origen
      real(8), dimension(nmax) :: min_dist

!     Nodo anterior en el camino minimo
      integer, dimension(nmax) :: predecessor

!     Nodos ya visitados durante Dijkstra
      logical, dimension(nmax) :: visited

!     Camino reconstruido
      integer, dimension(nmax) :: path

!     Longitud / acumulacion asociada a los caminos hacia targets
      real(8), dimension(nmax) :: path_sum

!     group members
      integer, dimension(max_chain,ngroup_max) :: ngrp_member
!=======================================================================
! 11. VARIABLES AUXILIARES Y RESULTADOS
!=======================================================================

!     Variables escalares auxiliares
      real(8) :: d, min_value
      real(8) :: tempd
      real(8) :: dx, dy, dz

!=======================================================================
! 12. VARIABLES LOGICAS AUXILIARES
!=======================================================================

      logical :: close
      logical :: new_chain

!=======================================================================
! 13. VARIABLES DE TEXTO Y LECTURA DE ARCHIVOS PDB
!=======================================================================

!     Linea completa leida del PDB
      character(len=200) :: line

!     Cadenas de texto auxiliares
      character(len=200) :: str
      character(len=200) :: nombre

!     Informacion del atomo
      character(len=4) :: atnm_tmp
      character(len=2) :: element_tmp

!     Informacion de cadena
      character(len=1) :: chain_tmp
      character(len=1) :: current_chain

!     Informacion del ligando y residuos temporales
      character(len=3) :: lig_resname
      character(len=3) :: resname_tmp

!     Variables de texto auxiliares utilizadas por el programa
      character(len=2) :: exh
      character(len=5) :: suffix

!=======================================================================

      end module fmcp_data
