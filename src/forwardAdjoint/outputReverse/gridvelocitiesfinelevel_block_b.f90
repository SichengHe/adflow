   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade 3.6 (r4159) - 21 Sep 2011 10:11
   !
   !  Differentiation of gridvelocitiesfinelevel_block in reverse (adjoint) mode:
   !   gradient     of useful results: *sfacei *sfacej *sfacek *x
   !                *si *sj *sk
   !   with respect to varying inputs: *x *si *sj *sk omegafourbeta
   !                *coscoeffourmach *coefpolbeta *coscoeffouralpha
   !                rotpoint *coscoeffourbeta omegafourxrot *sincoeffourmach
   !                *coefpolalpha omegafouryrot omegafourzrot omegafouralpha
   !                omegafourmach *coefpolmach *sincoeffourbeta *sincoeffouralpha
   !                *cgnsdoms.rotcenter *cgnsdoms.rotrate deltat
   !   Plus diff mem management of: sfacei:in sfacej:in sfacek:in
   !                x:in si:in sj:in sk:in coscoeffourzrot:in sincoeffourxrot:in
   !                sincoeffouryrot:in sincoeffourzrot:in coefpolxrot:in
   !                coefpolyrot:in coefpolzrot:in coscoeffourxrot:in
   !                coscoeffouryrot:in coeftime:in cgnsdoms:in
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          gridVelocities.f90                              *
   !      * Author:        Edwin van der Weide                             *
   !      * Starting date: 02-23-2004                                      *
   !      * Last modified: 06-28-2005                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   SUBROUTINE GRIDVELOCITIESFINELEVEL_BLOCK_B(useoldcoor, t, sps)
   USE MONITOR
   USE CGNSGRID
   USE BLOCKPOINTERS_B
   USE INPUTTSSTABDERIV
   USE INPUTUNSTEADY
   USE INPUTPHYSICS
   USE COMMUNICATION
   USE ITERATION
   USE INPUTMOTION
   USE FLOWVARREFSTATE
   USE DIFFSIZES
   !  Hint: ISIZE1OFDrfcgnsdoms should be the size of dimension 1 of array *cgnsdoms
   IMPLICIT NONE
   !
   !      ******************************************************************
   !      *                                                                *
   !      * gridVelocitiesFineLevel computes the grid velocities for       *
   !      * the cell centers and the normal grid velocities for the faces  *
   !      * of moving blocks for the currently finest grid, i.e.           *
   !      * groundLevel. The velocities are computed at time t for         *
   !      * spectral mode sps. If useOldCoor is .true. the velocities      *
   !      * are determined using the unsteady time integrator in           *
   !      * combination with the old coordinates; otherwise the analytic   *
   !      * form is used.                                                  *
   !      *                                                                *
   !      ******************************************************************
   !
   !
   !      Subroutine arguments.
   !
   INTEGER(kind=inttype), INTENT(IN) :: sps
   LOGICAL, INTENT(IN) :: useoldcoor
   REAL(kind=realtype), DIMENSION(*), INTENT(IN) :: t
   !
   !      Local variables.
   !
   INTEGER(kind=inttype) :: nn, mm
   INTEGER(kind=inttype) :: i, j, k, ii, iie, jje, kke
   REAL(kind=realtype) :: oneover4dt, oneover8dt
   REAL(kind=realtype) :: velxgrid, velygrid, velzgrid, ainf
   REAL(kind=realtype) :: velxgrid0, velygrid0, velzgrid0
   REAL(kind=realtype), DIMENSION(3) :: sc, xc, xxc
   REAL(kind=realtype), DIMENSION(3) :: scb, xcb, xxcb
   REAL(kind=realtype), DIMENSION(3) :: rotcenter, rotrate
   REAL(kind=realtype), DIMENSION(3) :: rotratetemp
   REAL(kind=realtype), DIMENSION(3) :: offsetvector
   REAL(kind=realtype), DIMENSION(3, 3) :: rotratetrans
   REAL(kind=realtype), DIMENSION(3) :: rotationpoint
   REAL(kind=realtype), DIMENSION(3, 3) :: rotationmatrix, &
   &  derivrotationmatrix
   REAL(kind=realtype) :: tnew, told
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: sface
   REAL(kind=realtype), DIMENSION(:, :), POINTER :: sfaceb
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: xx, ss
   REAL(kind=realtype), DIMENSION(:, :, :), POINTER :: xxb, ssb
   REAL(kind=realtype), DIMENSION(:, :, :, :), POINTER :: xxold
   INTEGER(kind=inttype) :: liftindex
   REAL(kind=realtype) :: alpha, beta, intervalmach, alphats, &
   &  alphaincrement, betats, betaincrement
   REAL(kind=realtype), DIMENSION(3) :: veldir, liftdir, dragdir
   REAL(kind=realtype), DIMENSION(3) :: refdirection
   !Function Definitions
   REAL(kind=realtype) :: TSALPHA, TSBETA, TSMACH
   INTEGER :: ad_to
   INTEGER :: ad_to0
   INTEGER :: branch
   INTEGER :: ad_to1
   INTEGER :: ad_to2
   INTEGER :: ad_to3
   INTEGER :: ad_to4
   INTRINSIC COS
   REAL(kind=realtype) :: tempb9
   REAL(kind=realtype) :: tempb8
   REAL(kind=realtype) :: tempb7
   REAL(kind=realtype) :: tempb6
   REAL(kind=realtype) :: tempb5
   REAL(kind=realtype) :: tempb4
   REAL(kind=realtype) :: tempb3
   INTRINSIC SIN
   REAL(kind=realtype) :: tempb2
   REAL(kind=realtype) :: tempb1
   REAL(kind=realtype) :: tempb0
   REAL(kind=realtype) :: tempb10
   REAL(kind=realtype) :: tempb
   INTEGER :: ii1
   INTRINSIC SQRT
   INTERFACE 
   SUBROUTINE PUSHPOINTER4(pp)
   REAL, POINTER :: pp
   END SUBROUTINE PUSHPOINTER4
   SUBROUTINE LOOKPOINTER4(pp)
   REAL, POINTER :: pp
   END SUBROUTINE LOOKPOINTER4
   SUBROUTINE POPPOINTER4(pp)
   REAL, POINTER :: pp
   END SUBROUTINE POPPOINTER4
   END INTERFACE
      !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Compute the mesh velocity from the given mesh Mach number.
   !  aInf = sqrt(gammaInf*pInf/rhoInf)
   !  velxGrid = aInf*MachGrid(1)
   !  velyGrid = aInf*MachGrid(2)
   !  velzGrid = aInf*MachGrid(3)
   ainf = SQRT(gammainf*pinf/rhoinf)
   velxgrid0 = ainf*machgrid*(-veldirfreestream(1))
   velygrid0 = ainf*machgrid*(-veldirfreestream(2))
   velzgrid0 = ainf*machgrid*(-veldirfreestream(3))
   ! Compute the derivative of the rotation matrix and the rotation
   ! point; needed for velocity due to the rigid body rotation of
   ! the entire grid. It is assumed that the rigid body motion of
   ! the grid is only specified if there is only 1 section present.
   CALL DERIVATIVEROTMATRIXRIGID(derivrotationmatrix, rotationpoint, t&
   &                             (1))
   !compute the rotation matrix to update the velocities for the time
   !spectral stability derivative case...
   IF (tsstability) THEN
   ! Determine the time values of the old and new time level.
   ! It is assumed that the rigid body rotation of the mesh is only
   ! used when only 1 section is present.
   tnew = timeunsteady + timeunsteadyrestart
   told = tnew - t(1)
   IF ((tspmode .OR. tsqmode) .OR. tsrmode) THEN
   ! Compute the rotation matrix of the rigid body rotation as
   ! well as the rotation point; the latter may vary in time due
   ! to rigid body translation.
   CALL ROTMATRIXRIGIDBODY(tnew, told, rotationmatrix, &
   &                           rotationpoint)
   velxgrid0 = rotationmatrix(1, 1)*velxgrid0 + rotationmatrix(1, 2)*&
   &        velygrid0 + rotationmatrix(1, 3)*velzgrid0
   velygrid0 = rotationmatrix(2, 1)*velxgrid0 + rotationmatrix(2, 2)*&
   &        velygrid0 + rotationmatrix(2, 3)*velzgrid0
   velzgrid0 = rotationmatrix(3, 1)*velxgrid0 + rotationmatrix(3, 2)*&
   &        velygrid0 + rotationmatrix(3, 3)*velzgrid0
   ELSE IF (tsalphamode) THEN
   ! get the baseline alpha and determine the liftIndex
   CALL GETDIRANGLE(veldirfreestream, liftdirection, liftindex, &
   &                    alpha, beta)
   !Determine the alpha for this time instance
   alphaincrement = TSALPHA(degreepolalpha, coefpolalpha, &
   &        degreefouralpha, omegafouralpha, coscoeffouralpha, &
   &        sincoeffouralpha, t(1))
   alphats = alpha + alphaincrement
   !Determine the grid velocity for this alpha
   refdirection(:) = zero
   refdirection(1) = one
   CALL GETDIRVECTOR(refdirection, alphats, beta, veldir, &
   &                     liftindex)
   !do I need to update the lift direction and drag direction as well?
   !set the effictive grid velocity for this time interval
   velxgrid0 = ainf*machgrid*(-veldir(1))
   velygrid0 = ainf*machgrid*(-veldir(2))
   velzgrid0 = ainf*machgrid*(-veldir(3))
   ELSE IF (tsbetamode) THEN
   ! get the baseline alpha and determine the liftIndex
   CALL GETDIRANGLE(veldirfreestream, liftdirection, liftindex, &
   &                    alpha, beta)
   !Determine the alpha for this time instance
   betaincrement = TSBETA(degreepolbeta, coefpolbeta, &
   &        degreefourbeta, omegafourbeta, coscoeffourbeta, sincoeffourbeta&
   &        , t(1))
   betats = beta + betaincrement
   !Determine the grid velocity for this alpha
   refdirection(:) = zero
   refdirection(1) = one
   CALL GETDIRVECTOR(refdirection, alpha, betats, veldir, &
   &                     liftindex)
   !do I need to update the lift direction and drag direction as well?
   !set the effictive grid velocity for this time interval
   velxgrid0 = ainf*machgrid*(-veldir(1))
   velygrid0 = ainf*machgrid*(-veldir(2))
   velzgrid0 = ainf*machgrid*(-veldir(3))
   ELSE IF (tsmachmode) THEN
   !determine the mach number at this time interval
   intervalmach = TSMACH(degreepolmach, coefpolmach, &
   &        degreefourmach, omegafourmach, coscoeffourmach, sincoeffourmach&
   &        , t(1))
   !set the effective grid velocity
   velxgrid0 = ainf*(intervalmach+machgrid)*(-veldirfreestream(1))
   velygrid0 = ainf*(intervalmach+machgrid)*(-veldirfreestream(2))
   velzgrid0 = ainf*(intervalmach+machgrid)*(-veldirfreestream(3))
   END IF
   END IF
   IF (blockismoving) THEN
   ! Determine the situation we are having here.
   IF (useoldcoor) THEN
   !
   !            ************************************************************
   !            *                                                          *
   !            * The velocities must be determined via a finite           *
   !            * difference formula using the coordinates of the old      *
   !            * levels.                                                  *
   !            *                                                          *
   !            ************************************************************
   !
   ! Set the coefficients for the time integrator and store
   ! the inverse of the physical nonDimensional time step,
   ! divided by 4 and 8, a bit easier.
   CALL SETCOEFTIMEINTEGRATOR()
   oneover4dt = fourth*timeref/deltat
   !
   !            ************************************************************
   !            *                                                          *
   !            * Grid velocities of the cell centers, including the       *
   !            * 1st level halo cells.                                    *
   !            *                                                          *
   !            ************************************************************
   !
   ! Loop over the cells, including the 1st level halo's.
   DO k=1,ke
   DO j=1,je
   DO i=1,ie
   ! The velocity of the cell center is determined
   ! by a finite difference formula. First store
   ! the current coordinate, multiplied by 8 and
   ! coefTime(0) in sc.
   sc(1) = (x(i-1, j-1, k-1, 1)+x(i, j-1, k-1, 1)+x(i-1, j, k-1&
   &              , 1)+x(i, j, k-1, 1)+x(i-1, j-1, k, 1)+x(i, j-1, k, 1)+x(i&
   &              -1, j, k, 1)+x(i, j, k, 1))*coeftime(0)
   sc(2) = (x(i-1, j-1, k-1, 2)+x(i, j-1, k-1, 2)+x(i-1, j, k-1&
   &              , 2)+x(i, j, k-1, 2)+x(i-1, j-1, k, 2)+x(i, j-1, k, 2)+x(i&
   &              -1, j, k, 2)+x(i, j, k, 2))*coeftime(0)
   sc(3) = (x(i-1, j-1, k-1, 3)+x(i, j-1, k-1, 3)+x(i-1, j, k-1&
   &              , 3)+x(i, j, k-1, 3)+x(i-1, j-1, k, 3)+x(i, j-1, k, 3)+x(i&
   &              -1, j, k, 3)+x(i, j, k, 3))*coeftime(0)
   ! Loop over the older levels to complete the
   ! finite difference formula.
   DO ii=1,noldlevels
   sc(1) = sc(1) + (xold(ii, i-1, j-1, k-1, 1)+xold(ii, i, j-&
   &                1, k-1, 1)+xold(ii, i-1, j, k-1, 1)+xold(ii, i, j, k-1, &
   &                1)+xold(ii, i-1, j-1, k, 1)+xold(ii, i, j-1, k, 1)+xold(&
   &                ii, i-1, j, k, 1)+xold(ii, i, j, k, 1))*coeftime(ii)
   sc(2) = sc(2) + (xold(ii, i-1, j-1, k-1, 2)+xold(ii, i, j-&
   &                1, k-1, 2)+xold(ii, i-1, j, k-1, 2)+xold(ii, i, j, k-1, &
   &                2)+xold(ii, i-1, j-1, k, 2)+xold(ii, i, j-1, k, 2)+xold(&
   &                ii, i-1, j, k, 2)+xold(ii, i, j, k, 2))*coeftime(ii)
   sc(3) = sc(3) + (xold(ii, i-1, j-1, k-1, 3)+xold(ii, i, j-&
   &                1, k-1, 3)+xold(ii, i-1, j, k-1, 3)+xold(ii, i, j, k-1, &
   &                3)+xold(ii, i-1, j-1, k, 3)+xold(ii, i, j-1, k, 3)+xold(&
   &                ii, i-1, j, k, 3)+xold(ii, i, j, k, 3))*coeftime(ii)
   END DO
   END DO
   END DO
   END DO
   !
   !            ************************************************************
   !            *                                                          *
   !            * Normal grid velocities of the faces.                     *
   !            *                                                          *
   !            ************************************************************
   !
   ! Loop over the three directions.
   loopdir:DO mm=1,3
   ! Set the upper boundaries depending on the direction.
   SELECT CASE  (mm) 
   CASE (1_intType) 
   ! normals in i-direction
   iie = ie
   jje = je
   kke = ke
   CASE (2_intType) 
   ! normals in j-direction
   iie = je
   jje = ie
   kke = ke
   CASE (3_intType) 
   ! normals in k-direction
   iie = ke
   jje = ie
   kke = je
   END SELECT
   !
   !              **********************************************************
   !              *                                                        *
   !              * Normal grid velocities in generalized i-direction.     *
   !              * Mm == 1: i-direction                                   *
   !              * mm == 2: j-direction                                   *
   !              * mm == 3: k-direction                                   *
   !              *                                                        *
   !              **********************************************************
   !
   DO i=0,iie
   ! Set the pointers for the coordinates, normals and
   ! normal velocities for this generalized i-plane.
   ! This depends on the value of mm.
   SELECT CASE  (mm) 
   CASE (1_intType) 
   CALL PUSHPOINTER4(xxb)
   xxb => xb(i, :, :, :)
   CALL PUSHPOINTER4(xx)
   ! normals in i-direction
   xx => x(i, :, :, :)
   xxold => xold(:, i, :, :, :)
   CALL PUSHPOINTER4(ssb)
   ssb => sib(i, :, :, :)
   CALL PUSHPOINTER4(ss)
   ss => si(i, :, :, :)
   CALL PUSHPOINTER4(sfaceb)
   sfaceb => sfaceib(i, :, :)
   CALL PUSHPOINTER4(sface)
   sface => sfacei(i, :, :)
   CALL PUSHCONTROL2B(2)
   CASE (2_intType) 
   CALL PUSHPOINTER4(xxb)
   xxb => xb(:, i, :, :)
   CALL PUSHPOINTER4(xx)
   ! normals in j-direction
   xx => x(:, i, :, :)
   xxold => xold(:, :, i, :, :)
   CALL PUSHPOINTER4(ssb)
   ssb => sjb(:, i, :, :)
   CALL PUSHPOINTER4(ss)
   ss => sj(:, i, :, :)
   CALL PUSHPOINTER4(sfaceb)
   sfaceb => sfacejb(:, i, :)
   CALL PUSHPOINTER4(sface)
   sface => sfacej(:, i, :)
   CALL PUSHCONTROL2B(1)
   CASE (3_intType) 
   CALL PUSHPOINTER4(xxb)
   xxb => xb(:, :, i, :)
   CALL PUSHPOINTER4(xx)
   ! normals in k-direction
   xx => x(:, :, i, :)
   xxold => xold(:, :, :, i, :)
   CALL PUSHPOINTER4(ssb)
   ssb => skb(:, :, i, :)
   CALL PUSHPOINTER4(ss)
   ss => sk(:, :, i, :)
   CALL PUSHPOINTER4(sfaceb)
   sfaceb => sfacekb(:, :, i)
   CALL PUSHPOINTER4(sface)
   sface => sfacek(:, :, i)
   CALL PUSHCONTROL2B(0)
   CASE DEFAULT
   CALL PUSHCONTROL2B(3)
   END SELECT
   ! Loop over the k and j-direction of this
   ! generalized i-face. Note that due to the usage of
   ! the pointers xx and xxOld an offset of +1 must be
   ! used in the coordinate arrays, because x and xOld
   ! originally start at 0 for the i, j and k indices.
   DO k=1,kke
   DO j=1,jje
   CALL PUSHREAL8(sc(1))
   ! The velocity of the face center is determined
   ! by a finite difference formula. First store
   ! the current coordinate, multiplied by 4 and
   ! coefTime(0) in sc.
   sc(1) = coeftime(0)*(xx(j+1, k+1, 1)+xx(j, k+1, 1)+xx(j+1&
   &                , k, 1)+xx(j, k, 1))
   CALL PUSHREAL8(sc(2))
   sc(2) = coeftime(0)*(xx(j+1, k+1, 2)+xx(j, k+1, 2)+xx(j+1&
   &                , k, 2)+xx(j, k, 2))
   CALL PUSHREAL8(sc(3))
   sc(3) = coeftime(0)*(xx(j+1, k+1, 3)+xx(j, k+1, 3)+xx(j+1&
   &                , k, 3)+xx(j, k, 3))
   ! Loop over the older levels to complete the
   ! finite difference.
   DO ii=1,noldlevels
   CALL PUSHREAL8(sc(1))
   sc(1) = sc(1) + coeftime(ii)*(xxold(ii, j+1, k+1, 1)+&
   &                  xxold(ii, j, k+1, 1)+xxold(ii, j+1, k, 1)+xxold(ii, j&
   &                  , k, 1))
   CALL PUSHREAL8(sc(2))
   sc(2) = sc(2) + coeftime(ii)*(xxold(ii, j+1, k+1, 2)+&
   &                  xxold(ii, j, k+1, 2)+xxold(ii, j+1, k, 2)+xxold(ii, j&
   &                  , k, 2))
   CALL PUSHREAL8(sc(3))
   sc(3) = sc(3) + coeftime(ii)*(xxold(ii, j+1, k+1, 3)+&
   &                  xxold(ii, j, k+1, 3)+xxold(ii, j+1, k, 3)+xxold(ii, j&
   &                  , k, 3))
   END DO
   ! Determine the dot product of sc and the normal
   ! and divide by 4 deltaT to obtain the correct
   ! value of the normal velocity.
   sface(j, k) = sc(1)*ss(j, k, 1) + sc(2)*ss(j, k, 2) + sc(3&
   &                )*ss(j, k, 3)
   sface(j, k) = sface(j, k)*oneover4dt
   END DO
   CALL PUSHINTEGER4(j - 1)
   END DO
   CALL PUSHINTEGER4(k - 1)
   END DO
   CALL PUSHINTEGER4(i - 1)
   END DO loopdir
   scb = 0.0_8
   DO mm=3,1,-1
   CALL POPINTEGER4(ad_to1)
   DO i=ad_to1,0,-1
   CALL POPINTEGER4(ad_to0)
   DO k=ad_to0,1,-1
   CALL POPINTEGER4(ad_to)
   DO j=ad_to,1,-1
   sfaceb(j, k) = oneover4dt*sfaceb(j, k)
   scb(1) = scb(1) + ss(j, k, 1)*sfaceb(j, k)
   ssb(j, k, 1) = ssb(j, k, 1) + sc(1)*sfaceb(j, k)
   scb(2) = scb(2) + ss(j, k, 2)*sfaceb(j, k)
   ssb(j, k, 2) = ssb(j, k, 2) + sc(2)*sfaceb(j, k)
   scb(3) = scb(3) + ss(j, k, 3)*sfaceb(j, k)
   ssb(j, k, 3) = ssb(j, k, 3) + sc(3)*sfaceb(j, k)
   sfaceb(j, k) = 0.0_8
   DO ii=noldlevels,1,-1
   CALL POPREAL8(sc(3))
   CALL POPREAL8(sc(2))
   CALL POPREAL8(sc(1))
   END DO
   CALL POPREAL8(sc(3))
   tempb2 = coeftime(0)*scb(3)
   xxb(j+1, k+1, 3) = xxb(j+1, k+1, 3) + tempb2
   xxb(j, k+1, 3) = xxb(j, k+1, 3) + tempb2
   xxb(j+1, k, 3) = xxb(j+1, k, 3) + tempb2
   xxb(j, k, 3) = xxb(j, k, 3) + tempb2
   scb(3) = 0.0_8
   CALL POPREAL8(sc(2))
   tempb3 = coeftime(0)*scb(2)
   xxb(j+1, k+1, 2) = xxb(j+1, k+1, 2) + tempb3
   xxb(j, k+1, 2) = xxb(j, k+1, 2) + tempb3
   xxb(j+1, k, 2) = xxb(j+1, k, 2) + tempb3
   xxb(j, k, 2) = xxb(j, k, 2) + tempb3
   scb(2) = 0.0_8
   CALL POPREAL8(sc(1))
   tempb4 = coeftime(0)*scb(1)
   xxb(j+1, k+1, 1) = xxb(j+1, k+1, 1) + tempb4
   xxb(j, k+1, 1) = xxb(j, k+1, 1) + tempb4
   xxb(j+1, k, 1) = xxb(j+1, k, 1) + tempb4
   xxb(j, k, 1) = xxb(j, k, 1) + tempb4
   scb(1) = 0.0_8
   END DO
   END DO
   CALL POPCONTROL2B(branch)
   IF (branch .LT. 2) THEN
   IF (branch .EQ. 0) THEN
   CALL POPPOINTER4(sface)
   CALL POPPOINTER4(sfaceb)
   CALL POPPOINTER4(ss)
   CALL POPPOINTER4(ssb)
   CALL POPPOINTER4(xx)
   CALL POPPOINTER4(xxb)
   ELSE
   CALL POPPOINTER4(sface)
   CALL POPPOINTER4(sfaceb)
   CALL POPPOINTER4(ss)
   CALL POPPOINTER4(ssb)
   CALL POPPOINTER4(xx)
   CALL POPPOINTER4(xxb)
   END IF
   ELSE IF (branch .EQ. 2) THEN
   CALL POPPOINTER4(sface)
   CALL POPPOINTER4(sfaceb)
   CALL POPPOINTER4(ss)
   CALL POPPOINTER4(ssb)
   CALL POPPOINTER4(xx)
   CALL POPPOINTER4(xxb)
   END IF
   END DO
   END DO
   DO k=ke,1,-1
   DO j=je,1,-1
   DO i=ie,1,-1
   tempb = coeftime(0)*scb(3)
   xb(i-1, j-1, k-1, 3) = xb(i-1, j-1, k-1, 3) + tempb
   xb(i, j-1, k-1, 3) = xb(i, j-1, k-1, 3) + tempb
   xb(i-1, j, k-1, 3) = xb(i-1, j, k-1, 3) + tempb
   xb(i, j, k-1, 3) = xb(i, j, k-1, 3) + tempb
   xb(i-1, j-1, k, 3) = xb(i-1, j-1, k, 3) + tempb
   xb(i, j-1, k, 3) = xb(i, j-1, k, 3) + tempb
   xb(i-1, j, k, 3) = xb(i-1, j, k, 3) + tempb
   xb(i, j, k, 3) = xb(i, j, k, 3) + tempb
   scb(3) = 0.0_8
   tempb0 = coeftime(0)*scb(2)
   xb(i-1, j-1, k-1, 2) = xb(i-1, j-1, k-1, 2) + tempb0
   xb(i, j-1, k-1, 2) = xb(i, j-1, k-1, 2) + tempb0
   xb(i-1, j, k-1, 2) = xb(i-1, j, k-1, 2) + tempb0
   xb(i, j, k-1, 2) = xb(i, j, k-1, 2) + tempb0
   xb(i-1, j-1, k, 2) = xb(i-1, j-1, k, 2) + tempb0
   xb(i, j-1, k, 2) = xb(i, j-1, k, 2) + tempb0
   xb(i-1, j, k, 2) = xb(i-1, j, k, 2) + tempb0
   xb(i, j, k, 2) = xb(i, j, k, 2) + tempb0
   scb(2) = 0.0_8
   tempb1 = coeftime(0)*scb(1)
   xb(i-1, j-1, k-1, 1) = xb(i-1, j-1, k-1, 1) + tempb1
   xb(i, j-1, k-1, 1) = xb(i, j-1, k-1, 1) + tempb1
   xb(i-1, j, k-1, 1) = xb(i-1, j, k-1, 1) + tempb1
   xb(i, j, k-1, 1) = xb(i, j, k-1, 1) + tempb1
   xb(i-1, j-1, k, 1) = xb(i-1, j-1, k, 1) + tempb1
   xb(i, j-1, k, 1) = xb(i, j-1, k, 1) + tempb1
   xb(i-1, j, k, 1) = xb(i-1, j, k, 1) + tempb1
   xb(i, j, k, 1) = xb(i, j, k, 1) + tempb1
   scb(1) = 0.0_8
   END DO
   END DO
   END DO
   ELSE
   !
   !            ************************************************************
   !            *                                                          *
   !            * The velocities must be determined analytically.          *
   !            *                                                          *
   !            ************************************************************
   !
   ! Store the rotation center and determine the
   ! nonDimensional rotation rate of this block. As the
   ! reference length is 1 timeRef == 1/uRef and at the end
   ! the nonDimensional velocity is computed.
   j = nbkglobal
   rotcenter = cgnsdoms(j)%rotcenter
   offsetvector = rotcenter - rotpoint
   rotrate = timeref*cgnsdoms(j)%rotrate
   IF (usewindaxis) THEN
   !determine the current angles from the free stream velocity
   CALL GETDIRANGLE(veldirfreestream, liftdirection, liftindex, &
   &                      alpha, beta)
   IF (liftindex .EQ. 2) THEN
   ! different coordinate system for aerosurf
   ! Wing is in z- direction
   rotratetrans(1, 1) = COS(alpha)*COS(beta)
   rotratetrans(1, 2) = -SIN(alpha)
   rotratetrans(1, 3) = -(COS(alpha)*SIN(beta))
   rotratetrans(2, 1) = SIN(alpha)*COS(beta)
   rotratetrans(2, 2) = COS(alpha)
   rotratetrans(2, 3) = -(SIN(alpha)*SIN(beta))
   rotratetrans(3, 1) = SIN(beta)
   rotratetrans(3, 2) = 0.0
   rotratetrans(3, 3) = COS(beta)
   ELSE IF (liftindex .EQ. 3) THEN
   ! Wing is in y- direction
   !Rotate the rotation rate from the wind axis back to the local body axis
   rotratetrans(1, 1) = COS(alpha)*COS(beta)
   rotratetrans(1, 2) = -(COS(alpha)*SIN(beta))
   rotratetrans(1, 3) = -SIN(alpha)
   rotratetrans(2, 1) = SIN(beta)
   rotratetrans(2, 2) = COS(beta)
   rotratetrans(2, 3) = 0.0
   rotratetrans(3, 1) = SIN(alpha)*COS(beta)
   rotratetrans(3, 2) = -(SIN(alpha)*SIN(beta))
   rotratetrans(3, 3) = COS(alpha)
   END IF
   rotratetemp = rotrate
   rotrate = 0.0
   DO i=1,3
   DO j=1,3
   rotrate(i) = rotrate(i) + rotratetemp(j)*rotratetrans(i, j)
   END DO
   END DO
   END IF
   !!$             if (useWindAxis)then
   !!$                !determine the current angles from the free stream velocity
   !!$                call getDirAngle(velDirFreestream,liftDirection,liftIndex,alpha,beta)
   !!$                !Rotate the rotation rate from the wind axis back to the local body axis
   !!$                !checkt he relationship between the differnt degrees of freedom!
   !!$                rotRateTrans(1,1)=cos(alpha)*cos(beta)
   !!$                rotRateTrans(1,2)=-cos(alpha)*sin(beta)
   !!$                rotRateTrans(1,3)=-sin(alpha)
   !!$                rotRateTrans(2,1)=sin(beta)
   !!$                rotRateTrans(2,2)=cos(beta)
   !!$                rotRateTrans(2,3)=0.0
   !!$                rotRateTrans(3,1)=sin(alpha)*cos(beta)
   !!$                rotRateTrans(3,2)=-sin(alpha)*sin(beta)
   !!$                rotRateTrans(3,3)=cos(alpha)
   !!$
   !!$                rotRateTemp = rotRate
   !!$                rotRate=0.0
   !!$                do i=1,3
   !!$                   do j=1,3
   !!$                      rotRate(i)=rotRate(i)+rotRateTemp(j)*rotRateTrans(i,j)
   !!$                   end do
   !!$                end do
   !!$             end if
   !subtract off the rotational velocity of the center of the grid
   ! to account for the added overall velocity.
   !             velxGrid =velxgrid0+ 1*(rotRate(2)*rotCenter(3) - rotRate(3)*rotCenter(2))
   !             velyGrid =velygrid0+ 1*(rotRate(3)*rotCenter(1) - rotRate(1)*rotCenter(3))
   !             velzGrid =velzgrid0+ 1*(rotRate(1)*rotCenter(2) - rotRate(2)*rotCenter(1))
   velxgrid = velxgrid0 + 1*(rotrate(2)*offsetvector(3)-rotrate(3)*&
   &        offsetvector(2)) + derivrotationmatrix(1, 1)*offsetvector(1) + &
   &        derivrotationmatrix(1, 2)*offsetvector(2) + derivrotationmatrix(&
   &        1, 3)*offsetvector(3)
   velygrid = velygrid0 + 1*(rotrate(3)*offsetvector(1)-rotrate(1)*&
   &        offsetvector(3)) + derivrotationmatrix(2, 1)*offsetvector(1) + &
   &        derivrotationmatrix(2, 2)*offsetvector(2) + derivrotationmatrix(&
   &        2, 3)*offsetvector(3)
   velzgrid = velzgrid0 + 1*(rotrate(1)*offsetvector(2)-rotrate(2)*&
   &        offsetvector(1)) + derivrotationmatrix(3, 1)*offsetvector(1) + &
   &        derivrotationmatrix(3, 2)*offsetvector(2) + derivrotationmatrix(&
   &        3, 3)*offsetvector(3)
   !add in rotmatrix*rotpoint....
   !
   !            ************************************************************
   !            *                                                          *
   !            * Grid velocities of the cell centers, including the       *
   !            * 1st level halo cells.                                    *
   !            *                                                          *
   !            ************************************************************
   !
   ! Loop over the cells, including the 1st level halo's.
   DO k=1,ke
   DO j=1,je
   DO i=1,ie
   ! Determine the coordinates of the cell center,
   ! which are stored in xc.
   xc(1) = eighth*(x(i-1, j-1, k-1, 1)+x(i, j-1, k-1, 1)+x(i-1&
   &              , j, k-1, 1)+x(i, j, k-1, 1)+x(i-1, j-1, k, 1)+x(i, j-1, k&
   &              , 1)+x(i-1, j, k, 1)+x(i, j, k, 1))
   xc(2) = eighth*(x(i-1, j-1, k-1, 2)+x(i, j-1, k-1, 2)+x(i-1&
   &              , j, k-1, 2)+x(i, j, k-1, 2)+x(i-1, j-1, k, 2)+x(i, j-1, k&
   &              , 2)+x(i-1, j, k, 2)+x(i, j, k, 2))
   xc(3) = eighth*(x(i-1, j-1, k-1, 3)+x(i, j-1, k-1, 3)+x(i-1&
   &              , j, k-1, 3)+x(i, j, k-1, 3)+x(i-1, j-1, k, 3)+x(i, j-1, k&
   &              , 3)+x(i-1, j, k, 3)+x(i, j, k, 3))
   ! Determine the coordinates relative to the
   ! center of rotation.
   xxc(1) = xc(1) - rotcenter(1)
   xxc(2) = xc(2) - rotcenter(2)
   xxc(3) = xc(3) - rotcenter(3)
   ! Determine the rotation speed of the cell center,
   ! which is omega*r.
   sc(1) = rotrate(2)*xxc(3) - rotrate(3)*xxc(2)
   sc(2) = rotrate(3)*xxc(1) - rotrate(1)*xxc(3)
   sc(3) = rotrate(1)*xxc(2) - rotrate(2)*xxc(1)
   ! Determine the coordinates relative to the
   ! rigid body rotation point.
   xxc(1) = xc(1) - rotationpoint(1)
   xxc(2) = xc(2) - rotationpoint(2)
   xxc(3) = xc(3) - rotationpoint(3)
   ! Determine the total velocity of the cell center.
   ! This is a combination of rotation speed of this
   ! block and the entire rigid body rotation.
   END DO
   END DO
   END DO
   !
   !            ************************************************************
   !            *                                                          *
   !            * Normal grid velocities of the faces.                     *
   !            *                                                          *
   !            ************************************************************
   !
   ! Loop over the three directions.
   loopdirection:DO mm=1,3
   ! Set the upper boundaries depending on the direction.
   SELECT CASE  (mm) 
   CASE (1_intType) 
   ! Normals in i-direction
   iie = ie
   jje = je
   kke = ke
   CASE (2_intType) 
   ! Normals in j-direction
   iie = je
   jje = ie
   kke = ke
   CASE (3_intType) 
   ! Normals in k-direction
   iie = ke
   jje = ie
   kke = je
   END SELECT
   !
   !              **********************************************************
   !              *                                                        *
   !              * Normal grid velocities in generalized i-direction.     *
   !              * mm == 1: i-direction                                   *
   !              * mm == 2: j-direction                                   *
   !              * mm == 3: k-direction                                   *
   !              *                                                        *
   !              **********************************************************
   !
   DO i=0,iie
   ! Set the pointers for the coordinates, normals and
   ! normal velocities for this generalized i-plane.
   ! This depends on the value of mm.
   SELECT CASE  (mm) 
   CASE (1_intType) 
   CALL PUSHPOINTER4(xxb)
   xxb => xb(i, :, :, :)
   CALL PUSHPOINTER4(xx)
   ! normals in i-direction
   xx => x(i, :, :, :)
   CALL PUSHPOINTER4(ssb)
   ssb => sib(i, :, :, :)
   CALL PUSHPOINTER4(ss)
   ss => si(i, :, :, :)
   CALL PUSHPOINTER4(sfaceb)
   sfaceb => sfaceib(i, :, :)
   CALL PUSHCONTROL2B(2)
   CASE (2_intType) 
   CALL PUSHPOINTER4(xxb)
   xxb => xb(:, i, :, :)
   CALL PUSHPOINTER4(xx)
   ! normals in j-direction
   xx => x(:, i, :, :)
   CALL PUSHPOINTER4(ssb)
   ssb => sjb(:, i, :, :)
   CALL PUSHPOINTER4(ss)
   ss => sj(:, i, :, :)
   CALL PUSHPOINTER4(sfaceb)
   sfaceb => sfacejb(:, i, :)
   CALL PUSHCONTROL2B(1)
   CASE (3_intType) 
   CALL PUSHPOINTER4(xxb)
   xxb => xb(:, :, i, :)
   CALL PUSHPOINTER4(xx)
   ! normals in k-direction
   xx => x(:, :, i, :)
   CALL PUSHPOINTER4(ssb)
   ssb => skb(:, :, i, :)
   CALL PUSHPOINTER4(ss)
   ss => sk(:, :, i, :)
   CALL PUSHPOINTER4(sfaceb)
   sfaceb => sfacekb(:, :, i)
   CALL PUSHCONTROL2B(0)
   CASE DEFAULT
   CALL PUSHCONTROL2B(3)
   END SELECT
   ! Loop over the k and j-direction of this generalized
   ! i-face. Note that due to the usage of the pointer
   ! xx an offset of +1 must be used in the coordinate
   ! array, because x originally starts at 0 for the
   ! i, j and k indices.
   DO k=1,kke
   DO j=1,jje
   ! Determine the coordinates of the face center,
   ! which are stored in xc.
   xc(1) = fourth*(xx(j+1, k+1, 1)+xx(j, k+1, 1)+xx(j+1, k, 1&
   &                )+xx(j, k, 1))
   xc(2) = fourth*(xx(j+1, k+1, 2)+xx(j, k+1, 2)+xx(j+1, k, 2&
   &                )+xx(j, k, 2))
   xc(3) = fourth*(xx(j+1, k+1, 3)+xx(j, k+1, 3)+xx(j+1, k, 3&
   &                )+xx(j, k, 3))
   ! Determine the coordinates relative to the
   ! center of rotation.
   xxc(1) = xc(1) - rotcenter(1)
   xxc(2) = xc(2) - rotcenter(2)
   xxc(3) = xc(3) - rotcenter(3)
   CALL PUSHREAL8(sc(1))
   ! Determine the rotation speed of the face center,
   ! which is omega*r.
   sc(1) = rotrate(2)*xxc(3) - rotrate(3)*xxc(2)
   CALL PUSHREAL8(sc(2))
   sc(2) = rotrate(3)*xxc(1) - rotrate(1)*xxc(3)
   CALL PUSHREAL8(sc(3))
   sc(3) = rotrate(1)*xxc(2) - rotrate(2)*xxc(1)
   ! Determine the coordinates relative to the
   ! rigid body rotation point.
   xxc(1) = xc(1) - rotationpoint(1)
   xxc(2) = xc(2) - rotationpoint(2)
   xxc(3) = xc(3) - rotationpoint(3)
   CALL PUSHREAL8(sc(1))
   ! Determine the total velocity of the cell face.
   ! This is a combination of rotation speed of this
   ! block and the entire rigid body rotation.
   sc(1) = sc(1) + velxgrid + derivrotationmatrix(1, 1)*xxc(1&
   &                ) + derivrotationmatrix(1, 2)*xxc(2) + &
   &                derivrotationmatrix(1, 3)*xxc(3)
   CALL PUSHREAL8(sc(2))
   sc(2) = sc(2) + velygrid + derivrotationmatrix(2, 1)*xxc(1&
   &                ) + derivrotationmatrix(2, 2)*xxc(2) + &
   &                derivrotationmatrix(2, 3)*xxc(3)
   CALL PUSHREAL8(sc(3))
   sc(3) = sc(3) + velzgrid + derivrotationmatrix(3, 1)*xxc(1&
   &                ) + derivrotationmatrix(3, 2)*xxc(2) + &
   &                derivrotationmatrix(3, 3)*xxc(3)
   ! Store the dot product of grid velocity sc and
   ! the normal ss in sFace.
   sface(j, k) = sc(1)*ss(j, k, 1) + sc(2)*ss(j, k, 2) + sc(3&
   &                )*ss(j, k, 3)
   END DO
   CALL PUSHINTEGER4(j - 1)
   END DO
   CALL PUSHINTEGER4(k - 1)
   END DO
   CALL PUSHINTEGER4(i - 1)
   END DO loopdirection
   xcb = 0.0_8
   xxcb = 0.0_8
   scb = 0.0_8
   DO mm=3,1,-1
   CALL POPINTEGER4(ad_to4)
   DO i=ad_to4,0,-1
   CALL POPINTEGER4(ad_to3)
   DO k=ad_to3,1,-1
   CALL POPINTEGER4(ad_to2)
   DO j=ad_to2,1,-1
   scb(1) = scb(1) + ss(j, k, 1)*sfaceb(j, k)
   ssb(j, k, 1) = ssb(j, k, 1) + sc(1)*sfaceb(j, k)
   scb(2) = scb(2) + ss(j, k, 2)*sfaceb(j, k)
   ssb(j, k, 2) = ssb(j, k, 2) + sc(2)*sfaceb(j, k)
   scb(3) = scb(3) + ss(j, k, 3)*sfaceb(j, k)
   ssb(j, k, 3) = ssb(j, k, 3) + sc(3)*sfaceb(j, k)
   sfaceb(j, k) = 0.0_8
   CALL POPREAL8(sc(3))
   xxcb(1) = xxcb(1) + derivrotationmatrix(3, 1)*scb(3)
   xxcb(2) = xxcb(2) + derivrotationmatrix(3, 2)*scb(3)
   xxcb(3) = xxcb(3) + derivrotationmatrix(3, 3)*scb(3)
   CALL POPREAL8(sc(2))
   xxcb(1) = xxcb(1) + derivrotationmatrix(2, 1)*scb(2)
   xxcb(2) = xxcb(2) + derivrotationmatrix(2, 2)*scb(2)
   xxcb(3) = xxcb(3) + derivrotationmatrix(2, 3)*scb(2)
   CALL POPREAL8(sc(1))
   xxcb(1) = xxcb(1) + derivrotationmatrix(1, 1)*scb(1)
   xxcb(2) = xxcb(2) + derivrotationmatrix(1, 2)*scb(1)
   xxcb(3) = xxcb(3) + derivrotationmatrix(1, 3)*scb(1)
   xcb(3) = xcb(3) + xxcb(3)
   xxcb(3) = 0.0_8
   xcb(2) = xcb(2) + xxcb(2)
   xxcb(2) = 0.0_8
   xcb(1) = xcb(1) + xxcb(1)
   xxcb(1) = 0.0_8
   CALL POPREAL8(sc(3))
   xxcb(2) = xxcb(2) + rotrate(1)*scb(3)
   xxcb(1) = xxcb(1) - rotrate(2)*scb(3)
   scb(3) = 0.0_8
   CALL POPREAL8(sc(2))
   xxcb(1) = xxcb(1) + rotrate(3)*scb(2)
   xxcb(3) = xxcb(3) - rotrate(1)*scb(2)
   scb(2) = 0.0_8
   CALL POPREAL8(sc(1))
   xxcb(3) = xxcb(3) + rotrate(2)*scb(1)
   xxcb(2) = xxcb(2) - rotrate(3)*scb(1)
   scb(1) = 0.0_8
   xcb(3) = xcb(3) + xxcb(3)
   xxcb(3) = 0.0_8
   xcb(2) = xcb(2) + xxcb(2)
   xxcb(2) = 0.0_8
   xcb(1) = xcb(1) + xxcb(1)
   xxcb(1) = 0.0_8
   tempb8 = fourth*xcb(3)
   xxb(j+1, k+1, 3) = xxb(j+1, k+1, 3) + tempb8
   xxb(j, k+1, 3) = xxb(j, k+1, 3) + tempb8
   xxb(j+1, k, 3) = xxb(j+1, k, 3) + tempb8
   xxb(j, k, 3) = xxb(j, k, 3) + tempb8
   xcb(3) = 0.0_8
   tempb9 = fourth*xcb(2)
   xxb(j+1, k+1, 2) = xxb(j+1, k+1, 2) + tempb9
   xxb(j, k+1, 2) = xxb(j, k+1, 2) + tempb9
   xxb(j+1, k, 2) = xxb(j+1, k, 2) + tempb9
   xxb(j, k, 2) = xxb(j, k, 2) + tempb9
   xcb(2) = 0.0_8
   tempb10 = fourth*xcb(1)
   xxb(j+1, k+1, 1) = xxb(j+1, k+1, 1) + tempb10
   xxb(j, k+1, 1) = xxb(j, k+1, 1) + tempb10
   xxb(j+1, k, 1) = xxb(j+1, k, 1) + tempb10
   xxb(j, k, 1) = xxb(j, k, 1) + tempb10
   xcb(1) = 0.0_8
   END DO
   END DO
   CALL POPCONTROL2B(branch)
   IF (branch .LT. 2) THEN
   IF (branch .EQ. 0) THEN
   CALL POPPOINTER4(sfaceb)
   CALL POPPOINTER4(ss)
   CALL POPPOINTER4(ssb)
   CALL POPPOINTER4(xx)
   CALL POPPOINTER4(xxb)
   ELSE
   CALL POPPOINTER4(sfaceb)
   CALL POPPOINTER4(ss)
   CALL POPPOINTER4(ssb)
   CALL POPPOINTER4(xx)
   CALL POPPOINTER4(xxb)
   END IF
   ELSE IF (branch .EQ. 2) THEN
   CALL POPPOINTER4(sfaceb)
   CALL POPPOINTER4(ss)
   CALL POPPOINTER4(ssb)
   CALL POPPOINTER4(xx)
   CALL POPPOINTER4(xxb)
   END IF
   END DO
   END DO
   DO k=ke,1,-1
   DO j=je,1,-1
   DO i=ie,1,-1
   xcb(3) = xcb(3) + xxcb(3)
   xxcb(3) = 0.0_8
   xcb(2) = xcb(2) + xxcb(2)
   xxcb(2) = 0.0_8
   xcb(1) = xcb(1) + xxcb(1)
   xxcb(1) = 0.0_8
   xxcb(2) = xxcb(2) + rotrate(1)*scb(3)
   xxcb(1) = xxcb(1) - rotrate(2)*scb(3)
   scb(3) = 0.0_8
   xxcb(1) = xxcb(1) + rotrate(3)*scb(2)
   xxcb(3) = xxcb(3) - rotrate(1)*scb(2)
   scb(2) = 0.0_8
   xxcb(3) = xxcb(3) + rotrate(2)*scb(1)
   xxcb(2) = xxcb(2) - rotrate(3)*scb(1)
   scb(1) = 0.0_8
   xcb(3) = xcb(3) + xxcb(3)
   xxcb(3) = 0.0_8
   xcb(2) = xcb(2) + xxcb(2)
   xxcb(2) = 0.0_8
   xcb(1) = xcb(1) + xxcb(1)
   xxcb(1) = 0.0_8
   tempb5 = eighth*xcb(3)
   xb(i-1, j-1, k-1, 3) = xb(i-1, j-1, k-1, 3) + tempb5
   xb(i, j-1, k-1, 3) = xb(i, j-1, k-1, 3) + tempb5
   xb(i-1, j, k-1, 3) = xb(i-1, j, k-1, 3) + tempb5
   xb(i, j, k-1, 3) = xb(i, j, k-1, 3) + tempb5
   xb(i-1, j-1, k, 3) = xb(i-1, j-1, k, 3) + tempb5
   xb(i, j-1, k, 3) = xb(i, j-1, k, 3) + tempb5
   xb(i-1, j, k, 3) = xb(i-1, j, k, 3) + tempb5
   xb(i, j, k, 3) = xb(i, j, k, 3) + tempb5
   xcb(3) = 0.0_8
   tempb6 = eighth*xcb(2)
   xb(i-1, j-1, k-1, 2) = xb(i-1, j-1, k-1, 2) + tempb6
   xb(i, j-1, k-1, 2) = xb(i, j-1, k-1, 2) + tempb6
   xb(i-1, j, k-1, 2) = xb(i-1, j, k-1, 2) + tempb6
   xb(i, j, k-1, 2) = xb(i, j, k-1, 2) + tempb6
   xb(i-1, j-1, k, 2) = xb(i-1, j-1, k, 2) + tempb6
   xb(i, j-1, k, 2) = xb(i, j-1, k, 2) + tempb6
   xb(i-1, j, k, 2) = xb(i-1, j, k, 2) + tempb6
   xb(i, j, k, 2) = xb(i, j, k, 2) + tempb6
   xcb(2) = 0.0_8
   tempb7 = eighth*xcb(1)
   xb(i-1, j-1, k-1, 1) = xb(i-1, j-1, k-1, 1) + tempb7
   xb(i, j-1, k-1, 1) = xb(i, j-1, k-1, 1) + tempb7
   xb(i-1, j, k-1, 1) = xb(i-1, j, k-1, 1) + tempb7
   xb(i, j, k-1, 1) = xb(i, j, k-1, 1) + tempb7
   xb(i-1, j-1, k, 1) = xb(i-1, j-1, k, 1) + tempb7
   xb(i, j-1, k, 1) = xb(i, j-1, k, 1) + tempb7
   xb(i-1, j, k, 1) = xb(i-1, j, k, 1) + tempb7
   xb(i, j, k, 1) = xb(i, j, k, 1) + tempb7
   xcb(1) = 0.0_8
   END DO
   END DO
   END DO
   END IF
   END IF
   omegafourbetab0 = 0.0_8
   coscoeffourmachb0 = 0.0_8
   coefpolbetab0 = 0.0_8
   coscoeffouralphab0 = 0.0_8
   rotpointb = 0.0_8
   coscoeffourbetab0 = 0.0_8
   omegafourxrotb0 = 0.0_8
   sincoeffourmachb0 = 0.0_8
   coefpolalphab0 = 0.0_8
   omegafouryrotb0 = 0.0_8
   omegafourzrotb0 = 0.0_8
   omegafouralphab0 = 0.0_8
   omegafourmachb0 = 0.0_8
   coefpolmachb0 = 0.0_8
   sincoeffourbetab0 = 0.0_8
   sincoeffouralphab0 = 0.0_8
   DO ii1=1,ISIZE1OFDrfcgnsdoms
   cgnsdomsb(ii1)%rotcenter = 0.0_8
   END DO
   DO ii1=1,ISIZE1OFDrfcgnsdoms
   cgnsdomsb(ii1)%rotrate = 0.0_8
   END DO
   deltatb = 0.0_8
   END SUBROUTINE GRIDVELOCITIESFINELEVEL_BLOCK_B