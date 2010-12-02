   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
   !
   !  Differentiation of extrapolate2ndhalonkpc in reverse (adjoint) mode:
   !   gradient     of useful results: padj0 padj1 padj2 wadj0 wadj1
   !                wadj2
   !   with respect to varying inputs: padj0 padj1 padj2 wadj0 wadj1
   !                wadj2
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          extrapolate2ndHaloAdj.f90                       *
   !      * Author:        Edwin van der Weide                             *
   !      *                Seongim Choi                                    *
   !      * Starting date: 03-21-2006                                      *
   !      * Last modified: 03-21-2006                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE EXTRAPOLATE2NDHALONKPC_B(nn, icbeg, icend, jcbeg, jcend, &
   &  ioffset, joffset, wadj0, wadj0b, wadj1, wadj1b, wadj2, wadj2b, padj0, &
   &  padj0b, padj1, padj1b, padj2, padj2b)
   USE CONSTANTS
   USE ITERATION
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * extrapolate2ndHaloAdj determines the states of the second      *
   !      * layer halo cells of subface nn of the block to which the       *
   !      * pointers in blockPointers currently point.                     *
   !      *                                                                *
   !      ******************************************************************
   !
   !irhoE
   !gammaInf, kPresent
   !nt1MG, nt2MG
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: nn
   INTEGER(kind=inttype), INTENT(IN) :: icbeg, icend, jcbeg, jcend
   INTEGER(kind=inttype), INTENT(IN) :: ioffset, joffset
   REAL(kind=realtype), DIMENSION(-2:2, -2:2, nw) :: wadj0, wadj1
   REAL(kind=realtype), DIMENSION(-2:2, -2:2, nw) :: wadj0b, wadj1b
   REAL(kind=realtype), DIMENSION(-2:2, -2:2, nw) :: wadj2
   REAL(kind=realtype), DIMENSION(-2:2, -2:2, nw) :: wadj2b
   REAL(kind=realtype), DIMENSION(-2:2, -2:2) :: padj0, padj1
   REAL(kind=realtype), DIMENSION(-2:2, -2:2) :: padj0b, padj1b
   REAL(kind=realtype), DIMENSION(-2:2, -2:2) :: padj2
   REAL(kind=realtype), DIMENSION(-2:2, -2:2) :: padj2b
   !
   !      Local parameter.
   !
   REAL(kind=realtype), PARAMETER :: factor=0.5_realType
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, j, l, ii, jj
   REAL(kind=realtype) :: ovgm1, gm53, factk
   REAL(kind=realtype) :: tmp
   REAL(kind=realtype) :: tmp0
   INTEGER :: branch
   INTRINSIC MAX
   REAL(kind=realtype) :: tmpb
   REAL(kind=realtype) :: tmp0b
   REAL(kind=realtype) :: tempb
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Easier storage of variables involving gamma.
   ovgm1 = one/(gammainf-one)
   gm53 = gammainf - five*third
   factk = -(ovgm1*gm53)
   ! Loop over the generic subface to set the state in the
   ! halo cells.
   DO j=jcbeg,jcend
   DO i=icbeg,icend
   CALL PUSHINTEGER4ARRAY(ii, inttype/4)
   ii = i - ioffset
   CALL PUSHINTEGER4ARRAY(jj, inttype/4)
   jj = j - joffset
   CALL PUSHREAL8ARRAY(wadj0(ii, jj, irho), realtype/8)
   ! Extrapolate the density, velocities and pressure.
   ! Make sure that a certain threshold is kept for the
   ! density and pressure.
   wadj0(ii, jj, irho) = two*wadj1(ii, jj, irho) - wadj2(ii, jj, irho&
   &        )
   IF (factor*wadj1(ii, jj, irho) .LT. wadj0(ii, jj, irho)) THEN
   CALL PUSHREAL8ARRAY(wadj0(ii, jj, irho), realtype/8)
   wadj0(ii, jj, irho) = wadj0(ii, jj, irho)
   CALL PUSHINTEGER4(1)
   ELSE
   CALL PUSHREAL8ARRAY(wadj0(ii, jj, irho), realtype/8)
   wadj0(ii, jj, irho) = factor*wadj1(ii, jj, irho)
   CALL PUSHINTEGER4(0)
   END IF
   CALL PUSHREAL8ARRAY(wadj0(ii, jj, ivx), realtype/8)
   wadj0(ii, jj, ivx) = two*wadj1(ii, jj, ivx) - wadj2(ii, jj, ivx)
   CALL PUSHREAL8ARRAY(wadj0(ii, jj, ivy), realtype/8)
   wadj0(ii, jj, ivy) = two*wadj1(ii, jj, ivy) - wadj2(ii, jj, ivy)
   CALL PUSHREAL8ARRAY(wadj0(ii, jj, ivz), realtype/8)
   wadj0(ii, jj, ivz) = two*wadj1(ii, jj, ivz) - wadj2(ii, jj, ivz)
   padj0(ii, jj) = two*padj1(ii, jj) - padj2(ii, jj)
   IF (factor*padj1(ii, jj) .LT. padj0(ii, jj)) THEN
   CALL PUSHINTEGER4(0)
   padj0(ii, jj) = padj0(ii, jj)
   ELSE
   padj0(ii, jj) = factor*padj1(ii, jj)
   CALL PUSHINTEGER4(1)
   END IF
   ! Extrapolate the turbulent variables. Use constant
   ! extrapolation.
   DO l=nt1mg,nt2mg
   CALL PUSHREAL8ARRAY(wadj0(i, j, l), realtype/8)
   wadj0(i, j, l) = wadj1(i, j, l)
   END DO
   ! Compute the total energy.
   tmp = ovgm1*padj0(ii, jj) + half*wadj0(ii, jj, irho)*(wadj0(ii, jj&
   &        , ivx)**2+wadj0(ii, jj, ivy)**2+wadj0(ii, jj, ivz)**2)
   CALL PUSHREAL8ARRAY(wadj0(ii, jj, irhoe), realtype/8)
   wadj0(ii, jj, irhoe) = tmp
   IF (kpresent) THEN
   tmp0 = wadj0(ii, jj, irhoe) - factk*wadj0(ii, jj, irho)*wadj0(ii&
   &          , jj, itu1)
   CALL PUSHREAL8ARRAY(wadj0(ii, jj, irhoe), realtype/8)
   wadj0(ii, jj, irhoe) = tmp0
   CALL PUSHINTEGER4(2)
   ELSE
   CALL PUSHINTEGER4(1)
   END IF
   END DO
   END DO
   DO j=jcend,jcbeg,-1
   DO i=icend,icbeg,-1
   CALL POPINTEGER4(branch)
   IF (.NOT.branch .LT. 2) THEN
   ii = i - ioffset
   jj = j - joffset
   CALL POPREAL8ARRAY(wadj0(ii, jj, irhoe), realtype/8)
   tmp0b = wadj0b(ii, jj, irhoe)
   wadj0b(ii, jj, irhoe) = tmp0b
   wadj0b(ii, jj, irho) = wadj0b(ii, jj, irho) - factk*wadj0(ii, jj&
   &          , itu1)*tmp0b
   wadj0b(ii, jj, itu1) = wadj0b(ii, jj, itu1) - factk*wadj0(ii, jj&
   &          , irho)*tmp0b
   END IF
   CALL POPREAL8ARRAY(wadj0(ii, jj, irhoe), realtype/8)
   tmpb = wadj0b(ii, jj, irhoe)
   wadj0b(ii, jj, irhoe) = 0.0
   tempb = half*wadj0(ii, jj, irho)*tmpb
   padj0b(ii, jj) = padj0b(ii, jj) + ovgm1*tmpb
   wadj0b(ii, jj, irho) = wadj0b(ii, jj, irho) + half*(wadj0(ii, jj, &
   &        ivx)**2+wadj0(ii, jj, ivy)**2+wadj0(ii, jj, ivz)**2)*tmpb
   wadj0b(ii, jj, ivx) = wadj0b(ii, jj, ivx) + 2*wadj0(ii, jj, ivx)*&
   &        tempb
   wadj0b(ii, jj, ivy) = wadj0b(ii, jj, ivy) + 2*wadj0(ii, jj, ivy)*&
   &        tempb
   wadj0b(ii, jj, ivz) = wadj0b(ii, jj, ivz) + 2*wadj0(ii, jj, ivz)*&
   &        tempb
   DO l=nt2mg,nt1mg,-1
   CALL POPREAL8ARRAY(wadj0(i, j, l), realtype/8)
   wadj1b(i, j, l) = wadj1b(i, j, l) + wadj0b(i, j, l)
   wadj0b(i, j, l) = 0.0
   END DO
   CALL POPINTEGER4(branch)
   IF (.NOT.branch .LT. 1) THEN
   padj1b(ii, jj) = padj1b(ii, jj) + factor*padj0b(ii, jj)
   padj0b(ii, jj) = 0.0
   END IF
   padj1b(ii, jj) = padj1b(ii, jj) + two*padj0b(ii, jj)
   padj2b(ii, jj) = padj2b(ii, jj) - padj0b(ii, jj)
   padj0b(ii, jj) = 0.0
   CALL POPREAL8ARRAY(wadj0(ii, jj, ivz), realtype/8)
   wadj1b(ii, jj, ivz) = wadj1b(ii, jj, ivz) + two*wadj0b(ii, jj, ivz&
   &        )
   wadj2b(ii, jj, ivz) = wadj2b(ii, jj, ivz) - wadj0b(ii, jj, ivz)
   wadj0b(ii, jj, ivz) = 0.0
   CALL POPREAL8ARRAY(wadj0(ii, jj, ivy), realtype/8)
   wadj1b(ii, jj, ivy) = wadj1b(ii, jj, ivy) + two*wadj0b(ii, jj, ivy&
   &        )
   wadj2b(ii, jj, ivy) = wadj2b(ii, jj, ivy) - wadj0b(ii, jj, ivy)
   wadj0b(ii, jj, ivy) = 0.0
   CALL POPREAL8ARRAY(wadj0(ii, jj, ivx), realtype/8)
   wadj1b(ii, jj, ivx) = wadj1b(ii, jj, ivx) + two*wadj0b(ii, jj, ivx&
   &        )
   wadj2b(ii, jj, ivx) = wadj2b(ii, jj, ivx) - wadj0b(ii, jj, ivx)
   wadj0b(ii, jj, ivx) = 0.0
   CALL POPINTEGER4(branch)
   IF (branch .LT. 1) THEN
   CALL POPREAL8ARRAY(wadj0(ii, jj, irho), realtype/8)
   wadj1b(ii, jj, irho) = wadj1b(ii, jj, irho) + factor*wadj0b(ii, &
   &          jj, irho)
   wadj0b(ii, jj, irho) = 0.0
   ELSE
   CALL POPREAL8ARRAY(wadj0(ii, jj, irho), realtype/8)
   END IF
   CALL POPREAL8ARRAY(wadj0(ii, jj, irho), realtype/8)
   wadj1b(ii, jj, irho) = wadj1b(ii, jj, irho) + two*wadj0b(ii, jj, &
   &        irho)
   wadj2b(ii, jj, irho) = wadj2b(ii, jj, irho) - wadj0b(ii, jj, irho)
   wadj0b(ii, jj, irho) = 0.0
   CALL POPINTEGER4ARRAY(jj, inttype/4)
   CALL POPINTEGER4ARRAY(ii, inttype/4)
   END DO
   END DO
   END SUBROUTINE EXTRAPOLATE2NDHALONKPC_B
