   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
   !
   !  Differentiation of bceulerwallforcesadj in reverse (adjoint) mode:
   !   gradient     of useful results: padj
   !   with respect to varying inputs: padj
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          bcEulerWallAdj.f90                              *
   !      * Author:        Edwin van der Weide                             *
   !      *                Seongim Choi,C.A.(Sandy) Mader                  *
   !      * Starting date: 03-21-2006                                      *
   !      * Last modified: 06-09-2008                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE BCEULERWALLFORCESADJ_B(wadj, padj, padjb)
   USE INPUTDISCRETIZATION
   USE FLOWVARREFSTATE
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * bcEulerWallAdj applies inviscid wall bcoundary condition       *
   !      * to the small set of cells wAdj (2,2,2) cube.                   *
   !      * This function is based on BCEulerWarll.f90 in src/solver       *
   !      * It only works for constant and linear pressure extrapolation   *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Subroutine Arguments
   REAL(kind=realtype), INTENT(IN) :: wadj(2, 2, 2, nw)
   REAL(kind=realtype), INTENT(INOUT) :: padj(3, 2, 2)
   REAL(kind=realtype) :: padjb(3, 2, 2)
   ! Local variables.
   INTEGER(kind=inttype) :: i, j, k
   EXTERNAL TERMINATE
   INTEGER :: branch
   INTRINSIC MAX
   SELECT CASE  (wallbctreatment) 
   CASE (constantpressure) 
   DO j=1,2
   DO i=1,2
   padj(1, i, j) = zero
   END DO
   END DO
   CALL PUSHINTEGER4(1)
   CASE (linextrapolpressure) 
   DO j=1,2
   DO i=1,2
   padj(1, i, j) = padj(3, i, j) - padj(2, i, j)
   END DO
   END DO
   CALL PUSHINTEGER4(2)
   CASE (quadextrapolpressure) 
   CALL PUSHINTEGER4(0)
   CASE (normalmomentum) 
   CALL PUSHINTEGER4(0)
   CASE DEFAULT
   CALL PUSHINTEGER4(0)
   END SELECT
   ! Determine the state in the halo cell. Again loop over
   ! the cell range for this subface.
   DO j=1,2
   DO i=1,2
   IF (zero .LT. padj(2, i, j) - padj(1, i, j)) THEN
   CALL PUSHINTEGER4(2)
   ELSE
   CALL PUSHINTEGER4(1)
   END IF
   END DO
   END DO
   DO j=2,1,-1
   DO i=2,1,-1
   CALL POPINTEGER4(branch)
   IF (branch .LT. 2) THEN
   padjb(1, i, j) = 0.0
   ELSE
   padjb(2, i, j) = padjb(2, i, j) + padjb(1, i, j)
   padjb(1, i, j) = -padjb(1, i, j)
   END IF
   END DO
   END DO
   CALL POPINTEGER4(branch)
   IF (branch .LT. 2) THEN
   IF (.NOT.branch .LT. 1) THEN
   DO j=2,1,-1
   DO i=2,1,-1
   padjb(1, i, j) = 0.0
   END DO
   END DO
   END IF
   ELSE
   DO j=2,1,-1
   DO i=2,1,-1
   padjb(3, i, j) = padjb(3, i, j) + padjb(1, i, j)
   padjb(2, i, j) = padjb(2, i, j) - padjb(1, i, j)
   padjb(1, i, j) = 0.0
   END DO
   END DO
   END IF
   END SUBROUTINE BCEULERWALLFORCESADJ_B
