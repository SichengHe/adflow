module LST

    use constants
    ! MPI comes from constants, so we need to avoid MPIF_H in PETSc
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    Mat dRdwLST

contains

    subroutine setupLST

        ! Setup the PETSc objects for the Newton-Krylov solver.
        ! destroyLST can be used to destroy the objects created in this function.

        use constants
        use stencils, only: euler_drdw_stencil, visc_drdw_stencil, N_euler_drdw, N_visc_drdw
        use communication, only: adflow_comm_world
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use flowVarRefState, only: nw, viscous
        use ADjointVars, only: nCellsLocal
        use utils, only: EChk
        use adjointUtils, only: myMatCreate, statePreAllocation
        use amg, only: setupAMG
        implicit none

        ! Working Variables
        integer(kind=intType) :: ierr, nDimw
        integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
        integer(kind=intType) :: n_stencil
        integer(kind=intType), dimension(:, :), pointer :: stencil
        integer(kind=intType) :: level


        nDimW = nw * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

        ! Create Pre-Conditioning Matrix
        allocate (nnzDiagonal(nCellsLocal(1_intType) * nTimeIntervalsSpectral), &
                    nnzOffDiag(nCellsLocal(1_intType) * nTimeIntervalsSpectral))

        if (viscous) then
            stencil => visc_drdw_stencil
            n_stencil = N_visc_drdw
        else
            stencil => euler_drdw_stencil
            n_stencil = N_euler_drdw
        end if

        level = 1
        call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW / nw, stencil, n_stencil, &
                                level, .False.)
        call myMatCreate(dRdwLST, nw, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
                            __FILE__, __LINE__)

        call matSetOption(dRdwLST, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        deallocate (nnzDiagonal, nnzOffDiag)

    end subroutine setupLST

    subroutine FormJacobianLST

        use constants
        use utils, only: EChk
        use adjointUtils, only: setupStateResidualMatrix
        implicit none

        ! Local Variables
        integer(kind=intType) :: ierr
        logical :: useAD, usePC, useTranspose, useObjective, useCoarseMats

        ! Assemble the approximate PC (fine level, level 1)
        useAD = .True.
        usePC = .False.
        useTranspose = .False.
        useObjective = .False.
        useCoarseMats = .False.

        call setupStateResidualMatrix(dRdwLST, useAD, usePC, useTranspose, &
                                      useObjective, .False., 1_intType, useCoarseMats=useCoarseMats)
    end subroutine FormJacobianLST


    subroutine destroyLST

        ! Destroy all the PETSc objects for lst.

        use constants
        use utils, only: EChk
        use amg, only: destroyAMG
        implicit none
        integer(kind=intType) :: ierr

        call MatDestroy(dRdwLST, ierr)
        call EChk(ierr, __FILE__, __LINE__)

    end subroutine destroyLST

end module LST