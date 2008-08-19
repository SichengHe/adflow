!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 2.2.4 (r2308) - 03/04/2008 10:03
!  
!  Differentiation of applyallbcforcecouplingadj in reverse (adjoint) mode:
!   gradient, with respect to input variables: winfadj padj pinfcorradj
!                wadj skadj sjadj siadj normadj
!   of linear combination of output variables: padj skadj sjadj
!                siadj
!
!      ******************************************************************
!      *                                                                *
!      * File:          applyAllBCForcesAdj.f90                         *
!      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
!      *                C.A.(Sandy) Mader                               *
!      * Starting date: 04-16-2008                                      *
!      * Last modified: 06-09-2008                                      *
!      *                                                                *
!      ******************************************************************
!
SUBROUTINE APPLYALLBCFORCECOUPLINGADJ_B(winfadj, winfadjb, pinfcorradj, &
&  pinfcorradjb, wadj, wadjb, padj, padjb, siadj, siadjb, sjadj, sjadjb&
&  , skadj, skadjb, normadj, normadjb, iibeg, iiend, jjbeg, jjend, i2beg&
&  , i2end, j2beg, j2end, secondhalo, mm)
  USE blockpointers
  USE bctypes
  USE inputdiscretization
  USE flowvarrefstate
  IMPLICIT NONE
!print *,'eulerwall applied'
!!$       ! Domain-interface boundary conditions,
!!$       ! when coupled with other solvers.
!!$       
!!$       call bcDomainInterface(secondHalo, correctForK)
!!$       
!!$       ! Supersonic inflow boundary conditions.
!!$       
!!$       call bcSupersonicInflow(secondHalo, correctForK)
!!$         enddo domains
!!$       enddo spectralLoop
!
!      ******************************************************************
!      *                                                                *
!      * applyAllBCAdj applies the possible boundary conditions for the *
!      * halo cells adjacent to the cell for which the residual needs   *
!      * to be computed.                                                *
!      *                                                                *
!      ******************************************************************
!
!, only : ie, ib, je, jb, ke, kb, nBocos, &
!         BCFaceID, BCType, BCData,p,w
!precond,choimerkle, etc...
!!$       use blockPointers
!!$       use flowVarRefState
!!$       use inputDiscretization
!!$       use inputTimeSpectral
!!$       use iteration
!!$       implicit none
!
!      Subroutine arguments.
!
  LOGICAL, INTENT(IN) :: secondhalo
  INTEGER(KIND=INTTYPE), INTENT(IN) :: iibeg, iiend, jjbeg, jjend
  INTEGER(KIND=INTTYPE), INTENT(IN) :: i2beg, i2end, j2beg, j2end
  INTEGER(KIND=INTTYPE), INTENT(IN) :: mm
!       integer(kind=intType) :: iCell, jCell, kCell
  REAL(KIND=REALTYPE), DIMENSION(0:ib, 0:jb, 0:kb, nw) :: wadj
  REAL(KIND=REALTYPE), DIMENSION(0:ib, 0:jb, 0:kb, nw) :: wadjb
  REAL(KIND=REALTYPE), DIMENSION(0:ib, 0:jb, 0:kb) :: padj
  REAL(KIND=REALTYPE), DIMENSION(0:ib, 0:jb, 0:kb) :: padjb
!real(kind=realType), dimension(-2:2,-2:2,-2:2,nw), &
!                                            intent(inout) :: wAdj
!real(kind=realType), dimension(-2:2,-2:2,-2:2),    &
!                                            intent(inout) :: pAdj
!real(kind=realType), dimension(-2:2,-2:2,-2:2,3), &
!                                            intent(in) :: siAdj, sjAdj, skAdj
!real(kind=realType), dimension(0:0,0:0,0:0), intent(in) :: volAdj
  REAL(KIND=REALTYPE) :: pinfcorradj
  REAL(KIND=REALTYPE) :: pinfcorradjb
  REAL(KIND=REALTYPE), DIMENSION(2, iibeg:iiend, jjbeg:jjend, 3) :: &
&  siadj
  REAL(KIND=REALTYPE), DIMENSION(2, iibeg:iiend, jjbeg:jjend, 3) :: &
&  siadjb
! notice the range of y dim is set 1:2 which corresponds to 1/jl
  REAL(KIND=REALTYPE), DIMENSION(iibeg:iiend, 2, jjbeg:jjend, 3) :: &
&  sjadj
  REAL(KIND=REALTYPE), DIMENSION(iibeg:iiend, 2, jjbeg:jjend, 3) :: &
&  sjadjb
! notice the range of z dim is set 1:2 which corresponds to 1/kl
  REAL(KIND=REALTYPE), DIMENSION(iibeg:iiend, jjbeg:jjend, 2, 3) :: &
&  skadj
  REAL(KIND=REALTYPE), DIMENSION(iibeg:iiend, jjbeg:jjend, 2, 3) :: &
&  skadjb
  REAL(KIND=REALTYPE), DIMENSION(iibeg:iiend, jjbeg:jjend, 3) :: normadj
  REAL(KIND=REALTYPE), DIMENSION(iibeg:iiend, jjbeg:jjend, 3) :: &
&  normadjb
!real(kind=realType), dimension(nBocos,-2:2,-2:2,3), intent(in) :: normAdj
  REAL(KIND=REALTYPE), DIMENSION(nw) :: winfadj
  REAL(KIND=REALTYPE), DIMENSION(nw) :: winfadjb
!
!      Local variables.
!
  INTEGER(KIND=INTTYPE) :: nn, sps
  INTEGER(KIND=INTTYPE) :: i, j, k, ii, jj, kk, l
  INTEGER(KIND=INTTYPE) :: istart, iend, jstart, jend, kstart, kend
  LOGICAL :: correctfork
  EXTERNAL TERMINATE
  INTEGER :: branch
  CALL PUSHREAL8ARRAY(padj, (ib+1)*(jb+1)*(kb+1))
  CALL PUSHREAL8ARRAY(wadj, (ib+1)*(jb+1)*(kb+1)*nw)
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!moved outside
!!$       ! Determine whether or not the total energy must be corrected
!!$       ! for the presence of the turbulent kinetic energy.
!!$
!!$       if( kPresent ) then
!!$         if((currentLevel <= groundLevel) .or. turbCoupled) then
!!$           correctForK = .true.
!!$         else
!!$           correctForK = .false.
!!$         endif
!!$       else
!!$         correctForK = .false.
!!$       endif
!!$       ! Loop over the number of spectral solutions.
!!$
!!$       spectralLoop: do sps=1,nTimeIntervalsSpectral
!!$
!!$         ! Loop over the number of blocks.
!!$
!!$         domains: do nn=1,nDom
!!$
!!$           ! Set the pointers for this block.
!!$
!!$           call setPointers(nn, currentLevel, sps)
!!$
!!$           ! Apply all the boundary conditions. The order is important.
! The symmetry boundary conditions.
!*************************
  CALL BCSYMMFORCECOUPLINGADJ(secondhalo, wadj, padj, normadj, mm, iibeg&
&                        , iiend, jjbeg, jjend, i2beg, i2end, j2beg, &
&                        j2end)
!**************************
!###       call bcSymmPolar(secondHalo)
!!$       ! call bcEulerWall(secondHalo, correctForK)
!!$
!!$       ! The viscous wall boundary conditions.
!!$
!!$       call bcNSWallAdiabatic( secondHalo, correctForK)
!!$       call bcNSWallIsothermal(secondHalo, correctForK)
!!$
!!$       ! The farfield is a special case, because the treatment
!!$       ! differs when preconditioning is used. Make that distinction
!!$       ! and call the appropriate routine.
!!$       
!!$!*******************************
  SELECT CASE  (precond) 
  CASE (noprecond) 
    CALL PUSHREAL8ARRAY(padj, (ib+1)*(jb+1)*(kb+1))
    CALL PUSHREAL8ARRAY(wadj, (ib+1)*(jb+1)*(kb+1)*nw)
!print *,'applying farfield'
    CALL BCFARFIELDFORCECOUPLINGADJ(secondhalo, winfadj, pinfcorradj, &
&                              wadj, padj, siadj, sjadj, skadj, normadj&
&                              , mm, iibeg, iiend, jjbeg, jjend, i2beg, &
&                              i2end, j2beg, j2end)
    CALL PUSHINTEGER4(1)
  CASE (turkel) 
    CALL PUSHINTEGER4(2)
  CASE (choimerkle) 
    CALL PUSHINTEGER4(3)
  CASE DEFAULT
    CALL PUSHINTEGER4(0)
  END SELECT
  CALL BCEULERWALLFORCECOUPLINGADJ_B(secondhalo, wadj, wadjb, padj, &
&                               padjb, siadj, siadjb, sjadj, sjadjb, &
&                               skadj, skadjb, normadj, normadjb, mm, &
&                               iibeg, iiend, jjbeg, jjend, i2beg, i2end&
&                               , j2beg, j2end)
  CALL POPINTEGER4(branch)
  IF (branch .LT. 2) THEN
    IF (branch .LT. 1) THEN
      winfadjb(:) = 0.0
      pinfcorradjb = 0.0
    ELSE
      CALL POPREAL8ARRAY(wadj, (ib+1)*(jb+1)*(kb+1)*nw)
      CALL POPREAL8ARRAY(padj, (ib+1)*(jb+1)*(kb+1))
      CALL BCFARFIELDFORCECOUPLINGADJ_B(secondhalo, winfadj, winfadjb, &
&                                  pinfcorradj, pinfcorradjb, wadj, &
&                                  wadjb, padj, padjb, siadj, sjadj, &
&                                  skadj, normadj, normadjb, mm, iibeg, &
&                                  iiend, jjbeg, jjend, i2beg, i2end, &
&                                  j2beg, j2end)
    END IF
  ELSE IF (branch .LT. 3) THEN
    winfadjb(:) = 0.0
    pinfcorradjb = 0.0
  ELSE
    winfadjb(:) = 0.0
    pinfcorradjb = 0.0
  END IF
  CALL POPREAL8ARRAY(wadj, (ib+1)*(jb+1)*(kb+1)*nw)
  CALL POPREAL8ARRAY(padj, (ib+1)*(jb+1)*(kb+1))
  CALL BCSYMMFORCECOUPLINGADJ_B(secondhalo, wadj, wadjb, padj, padjb, &
&                          normadj, normadjb, mm, iibeg, iiend, jjbeg, &
&                          jjend, i2beg, i2end, j2beg, j2end)
END SUBROUTINE APPLYALLBCFORCECOUPLINGADJ_B
