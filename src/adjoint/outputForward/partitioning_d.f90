!        generated by tapenade     (inria, tropics team)
!  tapenade 3.10 (r5363) -  9 sep 2014 09:53
!
module partitioning_d
  implicit none
! ----------------------------------------------------------------------
!                                                                      |
!                    no tapenade routine below this line               |
!                                                                      |
! ----------------------------------------------------------------------

contains
!  differentiation of timeperiodspectral in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: *sections.timeperiod
!   with respect to varying inputs: *sections.timeperiod omegafourier
!   rw status of diff variables: *sections.timeperiod:in-out omegafourier:in
!   plus diff mem management of: sections:in
  subroutine timeperiodspectral_d()
!
!       timeperiodspectral determines the time of one period for the
!       time spectral method. it is possible that sections have
!       different periodic times.
!
    use constants
    use communication, only : myid, adflow_comm_world
    use inputmotion, only : degreefourxrot, degreefouryrot, &
&   degreefourzrot, omegafouralpha, omegafourbeta, omegafourmach, &
&   omegafourxrot, omegafouryrot, omegafourzrot, degreefourmach, &
&   degreefouralpha, degreefourbeta
    use inputphysics, only : equationmode, flowtype
    use inputtimespectral, only : omegafourier, omegafourierd
    use section, only : sections, sectionsd, nsections
    use utils_d, only : terminate
    implicit none
! periodic time could not be determined. print an error
! message and exit.
!
!      local parameter.
!
    real(kind=realtype), parameter :: tol=1.e-5_realtype
!
!      local variables.
!
    integer :: ierr
    integer(kind=inttype) :: nn
    real(kind=realtype) :: tt, omega
    real(kind=realtype) :: ttd
    real(kind=realtype) :: timeperiod
    real(kind=realtype) :: timeperiodd
    logical :: timedetermined
    real(kind=realtype) :: tmpresult
! this routine is only used for the spectral solutions. return
! immediately if a different mode is solved.
    if (equationmode .ne. timespectral) then
      return
    else
! first check if a rotational frequency has been specified.
! only for external flows.
      timedetermined = .false.
      if (flowtype .eq. externalflow) then
! aeroelastic case
        if (omegafourier .gt. 0) then
          ttd = -(two*pi*omegafourierd/omegafourier**2)
          tt = two*pi/omegafourier
! check if a time period was already determined. if so, try
! to determine a common time. otherwise just copy the data.
          if (timedetermined) then
            timeperiodd = commontimespectral_d(timeperiod, tt, ttd, &
&             tmpresult)
            timeperiod = tmpresult
          else
            timeperiodd = ttd
            timeperiod = tt
            timedetermined = .true.
          end if
        else
          timeperiodd = 0.0_8
        end if
      else
        timeperiodd = 0.0_8
      end if
!!$         ! altitude.
!!$
!!$         if(degreefouraltitude > 0) then
!!$           tt = two*pi/omegafouraltitude
!!$
!!$           ! check if a time period was already determined. if so, try
!!$           ! to determine a common time. otherwise just copy the data.
!!$
!!$           if( timedetermined ) then
!!$             timeperiod = commontimespectral(timeperiod, tt)
!!$           else
!!$             timeperiod     = tt
!!$             timedetermined = .true.
!!$           endif
!!$         endif
! if it was possible to determine the time, copy it to the
! sections and return.
      if (timedetermined) then
        do nn=1,nsections
          sectionsd(nn)%timeperiod = timeperiodd/sections(nn)%nslices
          sections(nn)%timeperiod = timeperiod/sections(nn)%nslices
        end do
!print *,'sectiontimeperiod',sections(nn)%timeperiod,nn
        return
      else
! divide the periodic time by the number of slices to get the
! characteristic time for every section.
        do nn=1,nsections
          sectionsd(nn)%timeperiod = timeperiodd/sections(nn)%nslices
          sections(nn)%timeperiod = timeperiod/sections(nn)%nslices
        end do
! return if it was possible to determine the time.
        if (timedetermined) return
      end if
    end if
  end subroutine timeperiodspectral_d
  subroutine timeperiodspectral()
!
!       timeperiodspectral determines the time of one period for the
!       time spectral method. it is possible that sections have
!       different periodic times.
!
    use constants
    use communication, only : myid, adflow_comm_world
    use inputmotion, only : degreefourxrot, degreefouryrot, &
&   degreefourzrot, omegafouralpha, omegafourbeta, omegafourmach, &
&   omegafourxrot, omegafouryrot, omegafourzrot, degreefourmach, &
&   degreefouralpha, degreefourbeta
    use inputphysics, only : equationmode, flowtype
    use inputtimespectral, only : omegafourier
    use section, only : sections, nsections
    use utils_d, only : terminate
    implicit none
! periodic time could not be determined. print an error
! message and exit.
!
!      local parameter.
!
    real(kind=realtype), parameter :: tol=1.e-5_realtype
!
!      local variables.
!
    integer :: ierr
    integer(kind=inttype) :: nn
    real(kind=realtype) :: tt, omega
    real(kind=realtype) :: timeperiod
    logical :: timedetermined
! this routine is only used for the spectral solutions. return
! immediately if a different mode is solved.
    if (equationmode .ne. timespectral) then
      return
    else
! first check if a rotational frequency has been specified.
! only for external flows.
      timedetermined = .false.
      if (flowtype .eq. externalflow) then
! aeroelastic case
        if (omegafourier .gt. 0) then
          tt = two*pi/omegafourier
! check if a time period was already determined. if so, try
! to determine a common time. otherwise just copy the data.
          if (timedetermined) then
            timeperiod = commontimespectral(timeperiod, tt)
          else
            timeperiod = tt
            timedetermined = .true.
          end if
        end if
      end if
!!$         ! altitude.
!!$
!!$         if(degreefouraltitude > 0) then
!!$           tt = two*pi/omegafouraltitude
!!$
!!$           ! check if a time period was already determined. if so, try
!!$           ! to determine a common time. otherwise just copy the data.
!!$
!!$           if( timedetermined ) then
!!$             timeperiod = commontimespectral(timeperiod, tt)
!!$           else
!!$             timeperiod     = tt
!!$             timedetermined = .true.
!!$           endif
!!$         endif
! if it was possible to determine the time, copy it to the
! sections and return.
      if (timedetermined) then
        do nn=1,nsections
          sections(nn)%timeperiod = timeperiod/sections(nn)%nslices
        end do
!print *,'sectiontimeperiod',sections(nn)%timeperiod,nn
        return
      else
! divide the periodic time by the number of slices to get the
! characteristic time for every section.
        do nn=1,nsections
          sections(nn)%timeperiod = timeperiod/sections(nn)%nslices
        end do
! return if it was possible to determine the time.
        if (timedetermined) return
      end if
    end if
  end subroutine timeperiodspectral
!  differentiation of commontimespectral in forward (tangent) mode (with options i4 dr8 r8):
!   variations   of useful results: commontimespectral
!   with respect to varying inputs: t2
!   rw status of diff variables: commontimespectral:out t2:in
  function commontimespectral_d(t1, t2, t2d, commontimespectral)
!
!       the function commontimespectral determines the smallest
!       possible common time between t1 and t2, such that
!       tcommon = n1*t1 = n2*t2 and n1, n2 integers.
!
    use constants
    use communication, only : myid, adflow_comm_world
    use utils_d, only : terminate
    implicit none
!
!      function definition
!
    real(kind=realtype) :: commontimespectral
    real(kind=realtype) :: commontimespectral_d
!
!      function arguments.
!
    real(kind=realtype), intent(in) :: t1, t2
    real(kind=realtype), intent(in) :: t2d
!
!      local parameters.
!
    integer(kind=inttype), parameter :: nmax=100
    real(kind=realtype), parameter :: tol=1.e-5_realtype
!
!      local variables.
!
    integer :: ierr
    integer(kind=inttype) :: n1, n2
    real(kind=realtype) :: tt1, tt2, tt, ratio
    real(kind=realtype) :: tt1d
    intrinsic max
    intrinsic min
    intrinsic nint
    intrinsic abs
    if (t1 .lt. t2) then
      tt1d = t2d
      tt1 = t2
    else
      tt1 = t1
      tt1d = 0.0_8
    end if
    if (t1 .gt. t2) then
      tt2 = t2
    else
      tt2 = t1
    end if
    ratio = tt1/tt2
! loop to find the smallest integer values of n1 and n2, such
! that n1*tt1 = n2*tt2. note that due to the previous definition
! n2 >= n1.
    do n1=1,nmax
      tt = n1*ratio
      n2 = nint(tt)
      if (tt - n2 .ge. 0.) then
        tt = tt - n2
      else
        tt = -(tt-n2)
      end if
      if (tt .le. tol) goto 100
    end do
! check if a common time was found
! set the common time.
 100 commontimespectral_d = n1*tt1d
    commontimespectral = n1*tt1
  end function commontimespectral_d
  function commontimespectral(t1, t2)
!
!       the function commontimespectral determines the smallest
!       possible common time between t1 and t2, such that
!       tcommon = n1*t1 = n2*t2 and n1, n2 integers.
!
    use constants
    use communication, only : myid, adflow_comm_world
    use utils_d, only : terminate
    implicit none
!
!      function definition
!
    real(kind=realtype) :: commontimespectral
!
!      function arguments.
!
    real(kind=realtype), intent(in) :: t1, t2
!
!      local parameters.
!
    integer(kind=inttype), parameter :: nmax=100
    real(kind=realtype), parameter :: tol=1.e-5_realtype
!
!      local variables.
!
    integer :: ierr
    integer(kind=inttype) :: n1, n2
    real(kind=realtype) :: tt1, tt2, tt, ratio
    intrinsic max
    intrinsic min
    intrinsic nint
    intrinsic abs
    if (t1 .lt. t2) then
      tt1 = t2
    else
      tt1 = t1
    end if
    if (t1 .gt. t2) then
      tt2 = t2
    else
      tt2 = t1
    end if
    ratio = tt1/tt2
! loop to find the smallest integer values of n1 and n2, such
! that n1*tt1 = n2*tt2. note that due to the previous definition
! n2 >= n1.
    do n1=1,nmax
      tt = n1*ratio
      n2 = nint(tt)
      if (tt - n2 .ge. 0.) then
        tt = tt - n2
      else
        tt = -(tt-n2)
      end if
      if (tt .le. tol) goto 100
    end do
! check if a common time was found
! set the common time.
 100 commontimespectral = n1*tt1
  end function commontimespectral
end module partitioning_d
