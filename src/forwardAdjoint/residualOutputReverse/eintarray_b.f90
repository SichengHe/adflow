   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.6 (r4159) - 21 Sep 2011 10:11
   !
   !  Differentiation of eintarray in reverse (adjoint) mode:
   !   gradient     of useful results: eint rho
   !   with respect to varying inputs: rgas tref cv0 *cptrange *cpeint
   !                *(*cptempfit.constants) *cptempfit.eint0 cvn gammaconstant
   !                k p rho
   !   Plus diff mem management of: cptrange:in cpeint:in cptempfit:in
   !      ==================================================================
   SUBROUTINE EINTARRAY_B(rho, rhob, p, pb, k, kb, eint, eintb, correctfork&
   &  , kk)
   USE CPCURVEFITS
   USE INPUTPHYSICS
   USE CONSTANTS
   USE FLOWVARREFSTATE
   USE DIFFSIZES
   !  Hint: ISIZE1OFDrfcptempfit should be the size of dimension 1 of array *cptempfit
   IMPLICIT NONE
   INTEGER(kind=inttype), INTENT(IN) :: kk
   !
   !      ******************************************************************
   !      *                                                                *
   !      * EintArray computes the internal energy per unit mass from the  *
   !      * given density and pressure (and possibly turbulent energy) for *
   !      * the given kk elements of the arrays.                           *
   !      * For a calorically and thermally perfect gas the well-known     *
   !      * expression is used; for only a thermally perfect gas, cp is a  *
   !      * function of temperature, curve fits are used and a more        *
   !      * complex expression is obtained.                                *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Subroutine arguments.
   !
   REAL(kind=realtype), DIMENSION(kk), INTENT(IN) :: rho, p, k
   REAL(kind=realtype), DIMENSION(kk) :: rhob, pb, kb
   REAL(kind=realtype), DIMENSION(kk) :: eint
   REAL(kind=realtype), DIMENSION(kk) :: eintb
   LOGICAL, INTENT(IN) :: correctfork
   !
   !      Local parameter.
   !
   REAL(kind=realtype), PARAMETER :: twothird=two*third
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: i, nn, mm, ii, start
   REAL(kind=realtype) :: ovgm1, factk, pp, t, t2, scale
   REAL(kind=realtype) :: ppb, tb, t2b
   INTEGER :: ad_count
   INTEGER :: i0
   INTEGER :: branch
   INTEGER :: ad_to
   REAL(kind=realtype) :: tempb0
   REAL(kind=realtype) :: tempb
   INTRINSIC LOG
   INTEGER :: ii1
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Determine the cp model used in the computation.
   SELECT CASE  (cpmodel) 
   CASE (cpconstant) 
   ! Abbreviate 1/(gamma -1) a bit easier.
   ovgm1 = one/(gammaconstant-one)
   ! Second step. Correct the energy in case a turbulent kinetic
   ! energy is present.
   IF (correctfork) THEN
   factk = ovgm1*(five*third-gammaconstant)
   kb = 0.0_8
   DO i=kk,1,-1
   kb(i) = kb(i) - factk*eintb(i)
   END DO
   ELSE
   kb = 0.0_8
   END IF
   pb = 0.0_8
   DO i=kk,1,-1
   tempb = ovgm1*eintb(i)/rho(i)
   pb(i) = pb(i) + tempb
   rhob(i) = rhob(i) - p(i)*tempb/rho(i)
   eintb(i) = 0.0_8
   END DO
   CASE (cptempcurvefits) 
   !        ================================================================
   ! Cp as function of the temperature is given via curve fits.
   ! Store a scale factor to compute the nonDimensional
   ! internal energy.
   scale = rgas/tref
   ! Loop over the number of elements of the array
   DO i=1,kk
   CALL PUSHREAL8(pp)
   ! Compute the dimensional temperature.
   pp = p(i)
   IF (correctfork) THEN
   pp = pp - twothird*rho(i)*k(i)
   CALL PUSHCONTROL1B(0)
   ELSE
   CALL PUSHCONTROL1B(1)
   END IF
   t = tref*pp/(rgas*rho(i))
   ! Determine the case we are having here.
   IF (t .LE. cptrange(0)) THEN
   CALL PUSHCONTROL2B(0)
   ELSE IF (t .GE. cptrange(cpnparts)) THEN
   CALL PUSHCONTROL2B(1)
   ELSE
   ! Temperature is in the curve fit range.
   ! First find the valid range.
   ii = cpnparts
   start = 1
   ad_count = 1
   100    CALL PUSHINTEGER4(nn)
   ! Next guess for the interval.
   nn = start + ii/2
   ! Determine the situation we are having here.
   IF (t .GT. cptrange(nn)) THEN
   CALL PUSHCONTROL1B(0)
   ! Temperature is larger than the upper boundary of
   ! the current interval. Update the lower boundary.
   start = nn + 1
   ii = ii - 1
   ELSE IF (t .GE. cptrange(nn-1)) THEN
   GOTO 110
   ELSE
   CALL PUSHCONTROL1B(1)
   END IF
   ! This is the correct range. Exit the do-loop.
   ! Modify ii for the next branch to search.
   ii = ii/2
   ad_count = ad_count + 1
   GOTO 100
   110    CALL PUSHINTEGER4(ad_count)
   DO ii=1,cptempfit(nn)%nterm
   IF (cptempfit(nn)%exponents(ii) .EQ. -1) THEN
   CALL PUSHCONTROL1B(1)
   ELSE
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   CALL PUSHINTEGER4(ii - 1)
   CALL PUSHCONTROL2B(2)
   END IF
   ! Add the turbulent energy if needed.
   IF (correctfork) THEN
   CALL PUSHCONTROL1B(1)
   ELSE
   CALL PUSHCONTROL1B(0)
   END IF
   END DO
   kb = 0.0_8
   pb = 0.0_8
   DO i=kk,1,-1
   CALL POPCONTROL1B(branch)
   IF (branch .NE. 0) kb(i) = kb(i) + eintb(i)
   CALL POPCONTROL2B(branch)
   IF (branch .EQ. 0) THEN
   tb = scale*cv0*eintb(i)
   eintb(i) = 0.0_8
   ELSE IF (branch .EQ. 1) THEN
   tb = scale*cvn*eintb(i)
   eintb(i) = 0.0_8
   ELSE
   eintb(i) = scale*eintb(i)
   t = tref*pp/(rgas*rho(i))
   tb = 0.0_8
   CALL POPINTEGER4(ad_to)
   DO ii=ad_to,1,-1
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   mm = cptempfit(nn)%exponents(ii) + 1
   t2b = cptempfit(nn)%constants(ii)*eintb(i)/mm
   IF (.NOT.(t .LE. 0.0 .AND. (mm .EQ. 0.0 .OR. mm .NE. INT(mm)&
   &                ))) tb = tb + mm*t**(mm-1)*t2b
   ELSE
   tb = tb + cptempfit(nn)%constants(ii)*eintb(i)/t
   END IF
   END DO
   tb = tb - eintb(i)
   eintb(i) = 0.0_8
   CALL POPINTEGER4(ad_count)
   DO i0=1,ad_count
   IF (i0 .NE. 1) CALL POPCONTROL1B(branch)
   CALL POPINTEGER4(nn)
   END DO
   END IF
   tempb0 = tref*tb/(rgas*rho(i))
   ppb = tempb0
   rhob(i) = rhob(i) - pp*tempb0/rho(i)
   CALL POPCONTROL1B(branch)
   IF (branch .EQ. 0) THEN
   rhob(i) = rhob(i) - twothird*k(i)*ppb
   kb(i) = kb(i) - twothird*rho(i)*ppb
   END IF
   CALL POPREAL8(pp)
   pb(i) = pb(i) + ppb
   END DO
   CASE DEFAULT
   kb = 0.0_8
   pb = 0.0_8
   END SELECT
   rgasb = 0.0_8
   trefb = 0.0_8
   cv0b = 0.0_8
   cptrangeb = 0.0_8
   cpeintb = 0.0_8
   DO ii1=1,ISIZE1OFDrfcptempfit
   cptempfitb(ii1)%constants = 0.0_8
   END DO
   cptempfitb%eint0 = 0.0_8
   cvnb = 0.0_8
   gammaconstantb = 0.0_8
   END SUBROUTINE EINTARRAY_B
