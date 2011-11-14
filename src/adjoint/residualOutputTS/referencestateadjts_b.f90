   !        Generated by TAPENADE     (INRIA, Tropics team)
   !  Tapenade - Version 2.2 (r1239) - Wed 28 Jun 2006 04:59:55 PM CEST
   !  
   !  Differentiation of referencestateadjts in reverse (adjoint) mode:
   !   gradient, with respect to input variables: machadj
   !   of linear combination of output variables: uinfadj
   !
   !      ******************************************************************
   !      *                                                                *
   !      * File:          referenceStateAdj.f90                           *
   !      * Author:        Edwin van der Weide, Seonghyeon Hahn            *
   !      *                C.A.(Sandy) Mader                               *
   !      * Starting date: 05-29-2003                                      *
   !      * Last modified: 05-14-2008                                      *
   !      *                                                                *
   !      ******************************************************************
   !
   !       subroutine referenceStateAdj(velDirFreestreamAdj,liftDirectionAdj,&
   !            dragDirectionAdj, Machadj, MachCoefAdj,uInfAdj,prefAdj,&
   !            rhorefAdj, pinfdimAdj, rhoinfdimAdj, rhoinfAdj, pinfAdj,&
   !            murefAdj, timerefAdj)
   SUBROUTINE REFERENCESTATEADJTS_B(machadj, machadjb, machcoefadj, uinfadj&
   &  , uinfadjb, prefadj, rhorefadj, pinfdimadj, rhoinfdimadj, rhoinfadj, &
   &  pinfadj, murefadj, timerefadj)
   USE constants
   USE flowvarrefstate
   USE inputphysics
   IMPLICIT NONE
   !!$       !=================================================================
   !!$
   !!$       contains
   !!$
   !!$         !===============================================================
   !!$
   !!$         function maxValueSubface(var)
   !!$!
   !!$!        ****************************************************************
   !!$!        *                                                              *
   !!$!        * maxValueSubface determines the maximum value of var for      *
   !!$!        * currently active subface. If var is not associated this      *
   !!$!        * function returs -1.0.                                        *
   !!$!        *                                                              *
   !!$!        ****************************************************************
   !!$!
   !!$         implicit none
   !!$!
   !!$!        Function type
   !!$!
   !!$         real(kind=realType) :: maxValueSubface
   !!$!
   !!$!        Function argument.
   !!$!
   !!$         real(kind=realType), dimension(:,:), pointer :: var
   !!$!
   !!$!        Local variables.
   !!$!
   !!$         integer(kind=intType) :: i, j
   !!$!
   !!$!        ****************************************************************
   !!$!        *                                                              *
   !!$!        * Begin execution                                              *
   !!$!        *                                                              *
   !!$!        ****************************************************************
   !!$!
   !!$         ! Initialize the function to -1 and return immediately if
   !!$         ! var is not associated with data.
   !!$
   !!$         maxValueSubface = -one
   !!$         if(.not. associated(var)) return
   !!$
   !!$         ! Loop over the owned faces of the subface. As the cell range
   !!$         ! may contain halo values, the nodal range is used.
   !!$
   !!$         do j=(BCData(mm)%jnBeg+1),BCData(mm)%jnEnd
   !!$           do i=(BCData(mm)%inBeg+1),BCData(mm)%inEnd
   !!$             maxValueSubface = max(maxValueSubface,var(i,j))
   !!$           enddo
   !!$         enddo
   !!$
   !!$         end function maxValueSubface
   !!$
   REAL(KIND=REALTYPE) :: machadj, machadjb, machcoefadj, uinfadj, &
   &  uinfadjb
   REAL(KIND=REALTYPE) :: murefadj, timerefadj
   REAL(KIND=REALTYPE) :: pinfdimadj, rhoinfdimadj
   REAL(KIND=REALTYPE) :: pinfadj, rhoinfadj
   REAL(KIND=REALTYPE) :: prefadj, rhorefadj
   REAL(KIND=REALTYPE) :: dirglob(3), dirloc(3)
   INTEGER :: branch, ierr
   INTEGER(KIND=INTTYPE) :: mm, nn, sps
   REAL(KIND=REALTYPE) :: mx, my, mz, re, tinfdim, v
   REAL(KIND=REALTYPE) :: gm1, ratio, tmp
   REAL(KIND=REALTYPE) :: valglob(5), valloc(5)
   INTRINSIC SQRT
   !
   !      ******************************************************************
   !      *                                                                *
   !      * referenceState computes the reference state values in case     *
   !      * these have not been specified. A distinction is made between   *
   !      * internal and external flows. In case nothing has been          *
   !      * specified for the former a dimensional computation will be     *
   !      * made. For the latter the reference state is set to an          *
   !      * arbitrary state for an inviscid computation and computed for a *
   !      * viscous computation. Furthermore for internal flows an average *
   !      * velocity direction is computed from the boundary conditions,   *
   !      * which is used for initialization.                              *
   !      *                                                                *
   !      ******************************************************************
   !
   !!$       use BCTypes
   !!$       use block
   !!$       use communication
   !!$       use couplerParam
   !!$       use inputMotion
   !
   !      subroutine Variables
   !
   !       real(kind=realType), dimension(3), intent(inout) :: velDirFreestreamAdj
   !       real(kind=realType), dimension(3), intent(inout) :: liftDirectionAdj
   !       real(kind=realType), dimension(3), intent(inout) :: dragDirectionAdj
   !
   !      Local variables.
   !
   !!$       type(BCDataType), dimension(:), pointer :: BCData
   !!$!
   !!$!      Interfaces
   !!$!
   !!$       interface
   !!$         subroutine velMagnAndDirectionSubface(vmag, dir, &
   !!$                                               BCData, mm)
   !!$           use block
   !!$           implicit none
   !!$
   !!$           integer(kind=intType), intent(in) :: mm
   !!$           real(kind=realType), intent(out) :: vmag
   !!$           real(kind=realType), dimension(3), intent(inout) :: dir
   !!$           type(BCDataType), dimension(:), pointer :: BCData
   !!$         end subroutine velMagnAndDirectionSubface
   !!$       end interface
   !
   !      ******************************************************************
   !      *                                                                *
   !      * Begin execution                                                *
   !      *                                                                *
   !      ******************************************************************
   !
   ! Initialize the dimensional free stream temperature and pressure.
   ! From these values the density and viscosity is computed. For
   ! external viscous and internal computation this is corrected
   ! later on.
   pinfdimadj = prefadj
   IF (prefadj .LE. zero) THEN
   pinfdimadj = 101325.0_realType
   CALL PUSHINTEGER4(1)
   ELSE
   CALL PUSHINTEGER4(0)
   END IF
   tinfdim = tempfreestream
   rhoinfdimadj = pinfdimadj/(rgasdim*tinfdim)
   ! Check the flow type we are having here.
   IF (.NOT.flowtype .EQ. internalflow) THEN
   ! External flow. Compute the value of gammaInf.
   !         call computeGammaAdj(tempFreestream, gammaInf, 1_intType)
   CALL COMPUTEGAMMAADJTS(tempfreestream, gammainf, 1)
   ! In case of a viscous problem, compute the
   ! dimensional free stream density and pressure.
   !!$         if(equations == NSEquations .or. &
   !!$            equations == RANSEquations) then
   !!$
   !!$           ! Compute the x, y, and z-components of the Mach number
   !!$           ! relative to the body; i.e. the mesh velocity must be
   !!$           ! taken into account here.
   !!$
   !!$           mx = MachCoefAdj*velDirFreestreamAdj(1)
   !!$           my = MachCoefAdj*velDirFreestreamAdj(2)
   !!$           mz = MachCoefAdj*velDirFreestreamAdj(3)
   !!$
   !!$           ! Reynolds number per meter, the viscosity using sutherland's
   !!$           ! law and the free stream velocity relative to the body.
   !!$
   !!$           Re = Reynolds/ReynoldsLength
   !!$           muDim = muSuthDim                    &
   !!$                 * ((TSuthDim + SSuthDim)       &
   !!$                 /  (tempFreestream + SSuthDim)) &
   !!$                 * ((tempFreestream/TSuthDim)**1.5)
   !!$           v  = sqrt((mx*mx + my*my + mz*mz) &
   !!$              *      gammaInf*RGasDim*tempFreestream)
   !!$
   !!$           ! Compute the free stream density and pressure.
   !!$           ! Set TInfDim to tempFreestream.
   !!$
   !!$           rhoInfDimAdj = Re*muDim/v
   !!$           pInfDimAdj   = rhoInfDimAdj*RGasDim*tempFreestream
   !!$           TInfDim   = tempFreestream
   !!$
   !!$         endif
   ! In case the reference pressure, density and temperature were
   ! not specified, set them to the infinity values.
   IF (prefadj .LE. zero) THEN
   prefadj = pinfdimadj
   CALL PUSHINTEGER4(1)
   ELSE
   CALL PUSHINTEGER4(0)
   END IF
   IF (rhorefadj .LE. zero) THEN
   rhorefadj = rhoinfdimadj
   CALL PUSHINTEGER4(1)
   ELSE
   CALL PUSHINTEGER4(0)
   END IF
   ! Compute the value of muRef, such that the nonDimensional
   ! equations are identical to the dimensional ones.
   ! Note that in the non-dimensionalization of muRef there is
   ! a reference length. However this reference length is 1.0
   ! in this code, because the coordinates are converted to
   ! meters.
   ! Compute timeRef for a correct nonDimensionalization of the
   ! unsteady equations. Some story as for the reference viscosity
   ! concerning the reference length.
   ! Compute the nonDimensional pressure, density, velocity,
   ! viscosity and gas constant.
   pinfadj = pinfdim/prefadj
   rhoinfadj = rhoinfdimadj/rhorefadj
   machadjb = SQRT(gammainf*pinfadj/rhoinfadj)*uinfadjb
   CALL POPINTEGER4(branch)
   CALL POPINTEGER4(branch)
   CALL POPINTEGER4(branch)
   END IF
   END SUBROUTINE REFERENCESTATEADJTS_B
