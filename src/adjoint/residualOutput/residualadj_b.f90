!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 2.2.4 (r2308) - 03/04/2008 10:03
!  
!  Differentiation of residualadj in reverse (adjoint) mode:
!   gradient, with respect to input variables: *(cdisrk) kappacoef
!                voladj padj dwadj wadj skadj sjadj siadj
!   of linear combination of output variables: dwadj
!
!      ******************************************************************
!      *                                                                *
!      * File:          residual.f90                                    *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 04-21-2008                                      *
!      * Last modified: 04-28-2008                                      *
!      *                                                                *
!      ******************************************************************
!
SUBROUTINE RESIDUALADJ_B(wadj, wadjb, padj, padjb, siadj, siadjb, sjadj&
&  , sjadjb, skadj, skadjb, voladj, voladjb, normadj, dwadj, dwadjb, &
&  icell, jcell, kcell, correctfork)
  USE inputiteration
  USE cgnsgrid
  USE blockpointers
  USE inputtimespectral
  USE inputdiscretization
  USE iteration
  USE flowvarrefstate
  IMPLICIT NONE
!* real(iblank(iCell,jCell,kCell), realType)
!
!      ******************************************************************
!      *                                                                *
!      * residual computes the residual of the mean flow equations on   *
!      * the current MG level.                                          *
!      *                                                                *
!      ******************************************************************
!
!       Subroutine Variables
  INTEGER(KIND=INTTYPE) :: icell, jcell, kcell
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2, 3) :: siadj, sjadj, &
&  skadj
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2, 3) :: siadjb, sjadjb&
&  , skadjb
  REAL(KIND=REALTYPE), DIMENSION(0:0, 0:0, 0:0) :: voladj
  REAL(KIND=REALTYPE), DIMENSION(0:0, 0:0, 0:0) :: voladjb
  REAL(KIND=REALTYPE), DIMENSION(nbocos, -2:2, -2:2, 3), INTENT(INOUT) &
&  :: normadj
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2, nw) :: wadj
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2, nw) :: wadjb
  REAL(KIND=REALTYPE), DIMENSION(nw) :: dwadj
  REAL(KIND=REALTYPE), DIMENSION(nw) :: dwadjb
  REAL(KIND=REALTYPE), DIMENSION(nw) :: dwadj2
!integer(kind=intType), intent(in) :: discr
  LOGICAL, INTENT(IN) :: correctfork
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2) :: padj
  REAL(KIND=REALTYPE), DIMENSION(-2:2, -2:2, -2:2) :: padjb
  REAL(KIND=REALTYPE), DIMENSION(nw) :: fwadj
!
!      Local variables.
!
  INTEGER(KIND=INTTYPE) :: sps, nn, discr
  INTEGER(KIND=INTTYPE) :: i, j, k, l
  LOGICAL :: finegrid
  INTEGER :: branch
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!   Come back to this later....
!!$       ! Add the source terms from the level 0 cooling model.
!!$
!!$       call level0CoolingModel
! Set the value of rFil, which controls the fraction of the old
! dissipation residual to be used. This is only for the runge-kutta
! schemes; for other smoothers rFil is simply set to 1.0.
! Note the index rkStage+1 for cdisRK. The reason is that the
! residual computation is performed before rkStage is incremented.
  IF (smoother .EQ. rungekutta) THEN
    rfil = cdisrk(rkstage+1)
  ELSE
    rfil = one
  END IF
! Initialize the local arrays to monitor the massflows to zero.
! Set the value of the discretization, depending on the grid level,
! and the logical fineGrid, which indicates whether or not this
! is the finest grid level of the current mg cycle.
  discr = spacediscrcoarse
  IF (currentlevel .EQ. 1) discr = spacediscr
  finegrid = .false.
  IF (currentlevel .EQ. groundlevel) finegrid = .true.
!moved outside...
!!$
!!$       ! Loop over the number of spectral solutions and local
!!$       ! number of blocks.
!!$
!!$       spectralLoop: do sps=1,nTimeIntervalsSpectral
!!$         domainLoop: do nn=1,nDom
!!$
!!$           ! Set the pointers to this block and compute the central
!!$           ! inviscid flux.
!!$
!!$           call setPointers(nn, currentLevel, sps)
!********************
!               call inviscidUpwindFluxAdj2(wAdj,  pAdj,  dwAdj2, &
!                                        iCell, jCell, kCell,finegrid)
!!$               call inviscidUpwindFluxAdj2(w(icell-2:icell+2,jcell-2:jcell+2,kcell-2:kcell+2,:), p(icell-2:icell+2,jcell-2:jce
!ll+2,kcell-2:kcell+2),  dwAdj2, &
!!$                                        iCell, jCell, kCell,finegrid)
!!$               
!!$               fw(:,:,:,:) = 0.0
!!$
!!$               call inviscidUpwindFlux(fineGrid)
!!$               do l=1,nwf
!!$                  do k=2,kl
!!$                     do j=2,jl
!!$                        do i=2,il
!!$                           dw(i,j,k,l) = (dw(i,j,k,l) + fw(i,j,k,l)) &
!!$                                * real(iblank(i,j,k), realType)
!!$                        enddo
!!$                     enddo
!!$                  enddo
!!$               enddo
!!$               do i = 1,nw
!!$                  !if (abs(dwAdj(i)-dwAdj2(i))>0.0) then
!!$                  if (abs(dwAdj(i)-fw(icell,jcell,kcell,i))>0.0) then
!!$                  !if (1.0>0.0) then
!!$                     print *,abs(dwAdj(i)-fw(icell,jcell,kcell,i)),'dwadjup',dwAdj(i),'up2',dwAdj2(i),i,icell,jcell,kcell,fw(i
!cell,jcell,kcell,i)
!!$                  endif
!!$               enddo
! Compute the artificial dissipation fluxes.
! This depends on the parameter discr.
  SELECT CASE  (discr) 
  CASE (upwind) 
    CALL PUSHINTEGER4(1)
  CASE DEFAULT
    CALL PUSHINTEGER4(0)
  END SELECT
  l = 0
  CALL POPINTEGER4(branch)
  IF (branch .LT. 1) THEN
    padjb(-2:2, -2:2, -2:2) = 0.0
    wadjb(-2:2, -2:2, -2:2, :) = 0.0
    skadjb(-2:2, -2:2, -2:2, :) = 0.0
    sjadjb(-2:2, -2:2, -2:2, :) = 0.0
    siadjb(-2:2, -2:2, -2:2, :) = 0.0
  ELSE
    CALL INVISCIDUPWINDFLUXADJ_B(wadj, wadjb, padj, padjb, dwadj, dwadjb&
&                           , siadj, siadjb, sjadj, sjadjb, skadj, &
&                           skadjb, icell, jcell, kcell, finegrid)
  END IF
  CALL INVISCIDCENTRALFLUXADJ_B(wadj, wadjb, padj, padjb, dwadj, dwadjb&
&                          , siadj, siadjb, sjadj, sjadjb, skadj, skadjb&
&                          , voladj, voladjb, icell, jcell, kcell)
!  cdisrkb(:) = 0.0
!  kappacoefb = 0.0
END SUBROUTINE RESIDUALADJ_B
