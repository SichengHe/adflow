!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 2.2.4 (r2308) - 03/04/2008 10:03
!  
!  Differentiation of leftrightstateadj in reverse (adjoint) mode:
!   gradient, with respect to input variables: left right du1 du2
!                du3
!   of linear combination of output variables: left right du1 du2
!                du3
!      ==================================================================
SUBROUTINE LEFTRIGHTSTATEADJ_B(du1, du1b, du2, du2b, du3, du3b, left, &
&  leftb, right, rightb, nwint, omk, opk, factminmod, firstorderk)
  USE precision
  USE inputdiscretization
  USE constants
  IMPLICIT NONE
!
!        ****************************************************************
!        *                                                              *
!        * leftRightState computes the differences in the left and      *
!        * right state compared to the first order interpolation. For a *
!        * monotonic second order discretization the interpolations     *
!        * need to be nonlinear. The linear second order scheme can be  *
!        * stable (depending on the value of kappa), but it will have   *
!        * oscillations near discontinuities.                           *
!        *                                                              *
!        ****************************************************************
!
!
!        Local parameter.
!
  REAL(KIND=REALTYPE), PARAMETER :: epslim=1.e-10_realType
!
!        Subroutine arguments.
!
  INTEGER(KIND=INTTYPE) :: nwint
  REAL(KIND=REALTYPE) :: omk, opk, factminmod
  REAL(KIND=REALTYPE), DIMENSION(*) :: du1, du2, du3
  REAL(KIND=REALTYPE), DIMENSION(*) :: du1b, du2b, du3b
  REAL(KIND=REALTYPE), DIMENSION(*) :: left, right
  REAL(KIND=REALTYPE), DIMENSION(*) :: leftb, rightb
  LOGICAL :: firstorderk
!
!        Local variables.
!
  INTEGER(KIND=INTTYPE) :: l
  REAL(KIND=REALTYPE) :: rl1, rl2, rr1, rr2, tmp
  REAL(KIND=REALTYPE) :: rl1b, rl2b, rr1b, rr2b, tmpb
  REAL(KIND=REALTYPE) :: temp
  REAL(KIND=REALTYPE) :: temp0
  REAL(KIND=REALTYPE) :: temp1
  REAL(KIND=REALTYPE) :: temp2
  REAL(KIND=REALTYPE) :: temp3
  REAL(KIND=REALTYPE) :: temp4
  INTEGER :: branch
  REAL(KIND=REALTYPE) :: x3b
  REAL(KIND=REALTYPE) :: y1b
  REAL(KIND=REALTYPE) :: x6b
  REAL(KIND=REALTYPE) :: y4b
  REAL(KIND=REALTYPE) :: max2b
  REAL(KIND=REALTYPE) :: temp0b
  REAL(KIND=REALTYPE) :: max5b
  INTRINSIC MAX
  REAL(KIND=REALTYPE) :: x6
  INTRINSIC SIGN
  REAL(KIND=REALTYPE) :: x5
  REAL(KIND=REALTYPE) :: x4
  REAL(KIND=REALTYPE) :: x3
  REAL(KIND=REALTYPE) :: temp3b
  INTRINSIC ABS
  REAL(KIND=REALTYPE) :: x2
  REAL(KIND=REALTYPE) :: x1
  REAL(KIND=REALTYPE) :: x2b
  REAL(KIND=REALTYPE) :: temp2b3
  REAL(KIND=REALTYPE) :: temp2b2
  REAL(KIND=REALTYPE) :: temp2b1
  REAL(KIND=REALTYPE) :: temp2b0
  REAL(KIND=REALTYPE) :: x5b
  REAL(KIND=REALTYPE) :: y3b
  REAL(KIND=REALTYPE) :: max1b
  REAL(KIND=REALTYPE) :: tempb
  REAL(KIND=REALTYPE) :: max4b
  REAL(KIND=REALTYPE) :: temp2b
  REAL(KIND=REALTYPE) :: x1b
  REAL(KIND=REALTYPE) :: x4b
  REAL(KIND=REALTYPE) :: y2b
  REAL(KIND=REALTYPE) :: max3b
  INTRINSIC MIN
  REAL(KIND=REALTYPE) :: max6
  REAL(KIND=REALTYPE) :: max5
  REAL(KIND=REALTYPE) :: temp1b
  REAL(KIND=REALTYPE) :: max4
  REAL(KIND=REALTYPE) :: max6b
  REAL(KIND=REALTYPE) :: max3
  REAL(KIND=REALTYPE) :: max2
  REAL(KIND=REALTYPE) :: max1
  REAL(KIND=REALTYPE) :: y4
  REAL(KIND=REALTYPE) :: y3
  REAL(KIND=REALTYPE) :: temp4b
  REAL(KIND=REALTYPE) :: y2
  REAL(KIND=REALTYPE) :: y1
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution.                                             *
!        *                                                              *
!        ****************************************************************
!
! Linear interpolation; no limiter.
! Loop over the number of variables to be interpolated.
  SELECT CASE  (limiter) 
  CASE (nolimiter) 
    CALL PUSHINTEGER4(1)
  CASE (vanalbeda) 
!          ==============================================================
! Nonlinear interpolation using the van albeda limiter.
! Loop over the number of variables to be interpolated.
    DO l=1,nwint
      IF (du2(l) .GE. 0.) THEN
        x1 = du2(l)
        CALL PUSHINTEGER4(1)
      ELSE
        x1 = -du2(l)
        CALL PUSHINTEGER4(0)
      END IF
      IF (x1 .LT. epslim) THEN
        CALL PUSHREAL8(max1)
        max1 = epslim
        CALL PUSHINTEGER4(0)
      ELSE
        CALL PUSHREAL8(max1)
        max1 = x1
        CALL PUSHINTEGER4(1)
      END IF
      CALL PUSHREAL8(tmp)
! Compute the limiter argument rl1, rl2, rr1 and rr2.
! Note the cut off to 0.0.
      tmp = one/SIGN(max1, du2(l))
      IF (du1(l) .GE. 0.) THEN
        x3 = du1(l)
        CALL PUSHINTEGER4(1)
      ELSE
        x3 = -du1(l)
        CALL PUSHINTEGER4(0)
      END IF
      IF (x3 .LT. epslim) THEN
        CALL PUSHREAL8(max3)
        max3 = epslim
        CALL PUSHINTEGER4(0)
      ELSE
        CALL PUSHREAL8(max3)
        max3 = x3
        CALL PUSHINTEGER4(1)
      END IF
      y1 = du2(l)/SIGN(max3, du1(l))
      IF (zero .LT. y1) THEN
        CALL PUSHREAL8(rl1)
        rl1 = y1
        CALL PUSHINTEGER4(0)
      ELSE
        CALL PUSHREAL8(rl1)
        rl1 = zero
        CALL PUSHINTEGER4(1)
      END IF
      IF (zero .LT. du1(l)*tmp) THEN
        CALL PUSHREAL8(rl2)
        rl2 = du1(l)*tmp
        CALL PUSHINTEGER4(1)
      ELSE
        CALL PUSHREAL8(rl2)
        rl2 = zero
        CALL PUSHINTEGER4(0)
      END IF
      IF (zero .LT. du3(l)*tmp) THEN
        CALL PUSHREAL8(rr1)
        rr1 = du3(l)*tmp
        CALL PUSHINTEGER4(1)
      ELSE
        CALL PUSHREAL8(rr1)
        rr1 = zero
        CALL PUSHINTEGER4(0)
      END IF
      IF (du3(l) .GE. 0.) THEN
        x4 = du3(l)
        CALL PUSHINTEGER4(1)
      ELSE
        x4 = -du3(l)
        CALL PUSHINTEGER4(0)
      END IF
      IF (x4 .LT. epslim) THEN
        CALL PUSHREAL8(max4)
        max4 = epslim
        CALL PUSHINTEGER4(0)
      ELSE
        CALL PUSHREAL8(max4)
        max4 = x4
        CALL PUSHINTEGER4(1)
      END IF
      y2 = du2(l)/SIGN(max4, du3(l))
      IF (zero .LT. y2) THEN
        CALL PUSHREAL8(rr2)
        rr2 = y2
        CALL PUSHINTEGER4(0)
      ELSE
        CALL PUSHREAL8(rr2)
        rr2 = zero
        CALL PUSHINTEGER4(1)
      END IF
      CALL PUSHREAL8(rl1)
! Compute the corresponding limiter values.
      rl1 = rl1*(rl1+one)/(rl1*rl1+one)
      CALL PUSHREAL8(rl2)
      rl2 = rl2*(rl2+one)/(rl2*rl2+one)
      CALL PUSHREAL8(rr1)
      rr1 = rr1*(rr1+one)/(rr1*rr1+one)
      CALL PUSHREAL8(rr2)
      rr2 = rr2*(rr2+one)/(rr2*rr2+one)
! Compute the nonlinear corrections to the first order
! scheme.
    END DO
    CALL PUSHINTEGER4(2)
  CASE (minmod) 
!          ==============================================================
! Nonlinear interpolation using the minmod limiter.
! Loop over the number of variables to be interpolated.
    DO l=1,nwint
      IF (du2(l) .GE. 0.) THEN
        x2 = du2(l)
        CALL PUSHINTEGER4(1)
      ELSE
        x2 = -du2(l)
        CALL PUSHINTEGER4(0)
      END IF
      IF (x2 .LT. epslim) THEN
        CALL PUSHREAL8(max2)
        max2 = epslim
        CALL PUSHINTEGER4(0)
      ELSE
        CALL PUSHREAL8(max2)
        max2 = x2
        CALL PUSHINTEGER4(1)
      END IF
      CALL PUSHREAL8(tmp)
! Compute the limiter argument rl1, rl2, rr1 and rr2.
! Note the cut off to 0.0.
      tmp = one/SIGN(max2, du2(l))
      IF (du1(l) .GE. 0.) THEN
        x5 = du1(l)
        CALL PUSHINTEGER4(1)
      ELSE
        x5 = -du1(l)
        CALL PUSHINTEGER4(0)
      END IF
      IF (x5 .LT. epslim) THEN
        CALL PUSHREAL8(max5)
        max5 = epslim
        CALL PUSHINTEGER4(0)
      ELSE
        CALL PUSHREAL8(max5)
        max5 = x5
        CALL PUSHINTEGER4(1)
      END IF
      y3 = du2(l)/SIGN(max5, du1(l))
      IF (zero .LT. y3) THEN
        CALL PUSHREAL8(rl1)
        rl1 = y3
        CALL PUSHINTEGER4(0)
      ELSE
        CALL PUSHREAL8(rl1)
        rl1 = zero
        CALL PUSHINTEGER4(1)
      END IF
      IF (zero .LT. du1(l)*tmp) THEN
        CALL PUSHREAL8(rl2)
        rl2 = du1(l)*tmp
        CALL PUSHINTEGER4(1)
      ELSE
        CALL PUSHREAL8(rl2)
        rl2 = zero
        CALL PUSHINTEGER4(0)
      END IF
      IF (zero .LT. du3(l)*tmp) THEN
        CALL PUSHREAL8(rr1)
        rr1 = du3(l)*tmp
        CALL PUSHINTEGER4(1)
      ELSE
        CALL PUSHREAL8(rr1)
        rr1 = zero
        CALL PUSHINTEGER4(0)
      END IF
      IF (du3(l) .GE. 0.) THEN
        x6 = du3(l)
        CALL PUSHINTEGER4(1)
      ELSE
        x6 = -du3(l)
        CALL PUSHINTEGER4(0)
      END IF
      IF (x6 .LT. epslim) THEN
        CALL PUSHREAL8(max6)
        max6 = epslim
        CALL PUSHINTEGER4(0)
      ELSE
        CALL PUSHREAL8(max6)
        max6 = x6
        CALL PUSHINTEGER4(1)
      END IF
      y4 = du2(l)/SIGN(max6, du3(l))
      IF (zero .LT. y4) THEN
        CALL PUSHREAL8(rr2)
        rr2 = y4
        CALL PUSHINTEGER4(0)
      ELSE
        CALL PUSHREAL8(rr2)
        rr2 = zero
        CALL PUSHINTEGER4(1)
      END IF
      IF (one .GT. factminmod*rl1) THEN
        rl1 = factminmod*rl1
        CALL PUSHINTEGER4(1)
      ELSE
        rl1 = one
        CALL PUSHINTEGER4(0)
      END IF
      IF (one .GT. factminmod*rl2) THEN
        rl2 = factminmod*rl2
        CALL PUSHINTEGER4(1)
      ELSE
        rl2 = one
        CALL PUSHINTEGER4(0)
      END IF
      IF (one .GT. factminmod*rr1) THEN
        rr1 = factminmod*rr1
        CALL PUSHINTEGER4(1)
      ELSE
        rr1 = one
        CALL PUSHINTEGER4(0)
      END IF
      IF (one .GT. factminmod*rr2) THEN
        rr2 = factminmod*rr2
        CALL PUSHINTEGER4(1)
      ELSE
        rr2 = one
        CALL PUSHINTEGER4(0)
      END IF
    END DO
    CALL PUSHINTEGER4(3)
  CASE DEFAULT
    CALL PUSHINTEGER4(0)
  END SELECT
! In case only a first order scheme must be used for the
! turbulent transport equations, set the correction for the
! turbulent kinetic energy to 0.
  IF (firstorderk) THEN
    rightb(itu1) = 0.0
    leftb(itu1) = 0.0
  END IF
  CALL POPINTEGER4(branch)
  IF (branch .LT. 2) THEN
    IF (.NOT.branch .LT. 1) THEN
      DO l=nwint,1,-1
        du3b(l) = du3b(l) - omk*rightb(l)
        du2b(l) = du2b(l) + opk*leftb(l) - opk*rightb(l)
        rightb(l) = 0.0
        du1b(l) = du1b(l) + omk*leftb(l)
        leftb(l) = 0.0
      END DO
    END IF
  ELSE IF (branch .LT. 3) THEN
    DO l=nwint,1,-1
      rr1b = -(opk*du2(l)*rightb(l))
      du2b(l) = du2b(l) + opk*rl2*leftb(l) - opk*rr1*rightb(l)
      rr2b = -(omk*du3(l)*rightb(l))
      du3b(l) = du3b(l) - omk*rr2*rightb(l)
      rightb(l) = 0.0
      rl1b = omk*du1(l)*leftb(l)
      du1b(l) = du1b(l) + omk*rl1*leftb(l)
      rl2b = opk*du2(l)*leftb(l)
      leftb(l) = 0.0
      CALL POPREAL8(rr2)
      temp2b = rr2b/(one+rr2**2)
      rr2b = (2*rr2-rr2**2*(one+rr2)*2/(one+rr2**2)+one)*temp2b
      CALL POPREAL8(rr1)
      temp2b0 = rr1b/(one+rr1**2)
      rr1b = (2*rr1-rr1**2*(one+rr1)*2/(one+rr1**2)+one)*temp2b0
      CALL POPREAL8(rl2)
      temp2b1 = rl2b/(one+rl2**2)
      rl2b = (2*rl2-rl2**2*(one+rl2)*2/(one+rl2**2)+one)*temp2b1
      CALL POPREAL8(rl1)
      temp2b2 = rl1b/(one+rl1**2)
      rl1b = (2*rl1-rl1**2*(one+rl1)*2/(one+rl1**2)+one)*temp2b2
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(rr2)
        y2b = rr2b
      ELSE
        CALL POPREAL8(rr2)
        y2b = 0.0
      END IF
      temp1 = SIGN(max4, du3(l))
      temp1b = -(du2(l)*y2b/temp1**2)
      du2b(l) = du2b(l) + y2b/temp1
      max4b = SIGN(1.d0, max4*du3(l))*temp1b
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(max4)
        x4b = 0.0
      ELSE
        CALL POPREAL8(max4)
        x4b = max4b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        du3b(l) = du3b(l) - x4b
      ELSE
        du3b(l) = du3b(l) + x4b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(rr1)
        tmpb = 0.0
      ELSE
        CALL POPREAL8(rr1)
        du3b(l) = du3b(l) + tmp*rr1b
        tmpb = du3(l)*rr1b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(rl2)
      ELSE
        CALL POPREAL8(rl2)
        du1b(l) = du1b(l) + tmp*rl2b
        tmpb = tmpb + du1(l)*rl2b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(rl1)
        y1b = rl1b
      ELSE
        CALL POPREAL8(rl1)
        y1b = 0.0
      END IF
      temp0 = SIGN(max3, du1(l))
      temp0b = -(du2(l)*y1b/temp0**2)
      du2b(l) = du2b(l) + y1b/temp0
      max3b = SIGN(1.d0, max3*du1(l))*temp0b
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(max3)
        x3b = 0.0
      ELSE
        CALL POPREAL8(max3)
        x3b = max3b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        du1b(l) = du1b(l) - x3b
      ELSE
        du1b(l) = du1b(l) + x3b
      END IF
      CALL POPREAL8(tmp)
      temp = SIGN(max1, du2(l))
      tempb = -(one*tmpb/temp**2)
      max1b = SIGN(1.d0, max1*du2(l))*tempb
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(max1)
        x1b = 0.0
      ELSE
        CALL POPREAL8(max1)
        x1b = max1b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        du2b(l) = du2b(l) - x1b
      ELSE
        du2b(l) = du2b(l) + x1b
      END IF
    END DO
  ELSE
    DO l=nwint,1,-1
      rr1b = -(opk*du2(l)*rightb(l))
      du2b(l) = du2b(l) + opk*rl2*leftb(l) - opk*rr1*rightb(l)
      rr2b = -(omk*du3(l)*rightb(l))
      du3b(l) = du3b(l) - omk*rr2*rightb(l)
      rightb(l) = 0.0
      rl1b = omk*du1(l)*leftb(l)
      du1b(l) = du1b(l) + omk*rl1*leftb(l)
      rl2b = opk*du2(l)*leftb(l)
      leftb(l) = 0.0
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        rr2b = 0.0
      ELSE
        rr2b = factminmod*rr2b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        rr1b = 0.0
      ELSE
        rr1b = factminmod*rr1b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        rl2b = 0.0
      ELSE
        rl2b = factminmod*rl2b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        rl1b = 0.0
      ELSE
        rl1b = factminmod*rl1b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(rr2)
        y4b = rr2b
      ELSE
        CALL POPREAL8(rr2)
        y4b = 0.0
      END IF
      temp4 = SIGN(max6, du3(l))
      temp4b = -(du2(l)*y4b/temp4**2)
      du2b(l) = du2b(l) + y4b/temp4
      max6b = SIGN(1.d0, max6*du3(l))*temp4b
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(max6)
        x6b = 0.0
      ELSE
        CALL POPREAL8(max6)
        x6b = max6b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        du3b(l) = du3b(l) - x6b
      ELSE
        du3b(l) = du3b(l) + x6b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(rr1)
        tmpb = 0.0
      ELSE
        CALL POPREAL8(rr1)
        du3b(l) = du3b(l) + tmp*rr1b
        tmpb = du3(l)*rr1b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(rl2)
      ELSE
        CALL POPREAL8(rl2)
        du1b(l) = du1b(l) + tmp*rl2b
        tmpb = tmpb + du1(l)*rl2b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(rl1)
        y3b = rl1b
      ELSE
        CALL POPREAL8(rl1)
        y3b = 0.0
      END IF
      temp3 = SIGN(max5, du1(l))
      temp3b = -(du2(l)*y3b/temp3**2)
      du2b(l) = du2b(l) + y3b/temp3
      max5b = SIGN(1.d0, max5*du1(l))*temp3b
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(max5)
        x5b = 0.0
      ELSE
        CALL POPREAL8(max5)
        x5b = max5b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        du1b(l) = du1b(l) - x5b
      ELSE
        du1b(l) = du1b(l) + x5b
      END IF
      CALL POPREAL8(tmp)
      temp2 = SIGN(max2, du2(l))
      temp2b3 = -(one*tmpb/temp2**2)
      max2b = SIGN(1.d0, max2*du2(l))*temp2b3
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        CALL POPREAL8(max2)
        x2b = 0.0
      ELSE
        CALL POPREAL8(max2)
        x2b = max2b
      END IF
      CALL POPINTEGER4(branch)
      IF (branch .LT. 1) THEN
        du2b(l) = du2b(l) - x2b
      ELSE
        du2b(l) = du2b(l) + x2b
      END IF
    END DO
  END IF
END SUBROUTINE LEFTRIGHTSTATEADJ_B
