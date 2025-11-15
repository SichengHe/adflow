!
!      ******************************************************************
!      *                                                                *
!      * File:          resolventAPI.F90                                *
!      * Author:        ADflow development team                         *
!      * Starting date: 2025-11-15                                      *
!      * Last modified: 2025-11-15                                      *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * resolventAPI: Module containing subroutines for resolvent     *
!      * analysis interfacing with Python.                             *
!      *                                                                *
!      * Based on "Large-Scale Flow Control Performance Optimization   *
!      * via Differentiable Resolvent Analysis" by He et al.           *
!      *                                                                *
!      ******************************************************************
!
module resolventAPI

    use constants, only: realType, intType
    use ADjointPETSc, only: dRdWT
#include <petsc/finclude/petsc.h>
    use petsc

    implicit none

contains

    subroutine setupResolventMatrix(frozenTurb)
        !
        ! Setup the Jacobian matrix (dRdW) for resolvent analysis.
        ! This assembles the state residual Jacobian at the current flow solution.
        !
        ! Parameters:
        ! -----------
        ! frozenTurb : logical, intent(in)
        !     If true, use frozen turbulence (flow variables only)
        !     If false, include turbulence variables
        !
        use adjointUtils, only: setupStateResidualMatrix
        use adjointVars, only: adjointPETScVarsAllocated
        use iteration, only: currentLevel

        implicit none

        ! Input
        logical, intent(in) :: frozenTurb

        ! Local variables
        integer(kind=intType) :: level
        logical :: useAD, usePC, useTranspose, useObjective

        ! Set parameters for Jacobian assembly
        useAD = .True.          ! Use automatic differentiation
        usePC = .False.         ! Use exact Jacobian (not preconditioner)
        useTranspose = .True.   ! dRdWT is stored (transpose)
        useObjective = .False.  ! Don't need objective
        level = 1               ! Finest level

        ! Check if adjoint PETSc variables are allocated
        if (.not. adjointPETScVarsAllocated) then
            print *, "ERROR: Adjoint PETSc variables not allocated!"
            print *, "Call setAeroProblem and ensure adjoint is initialized."
            return
        end if

        ! Assemble the Jacobian matrix
        ! Note: This assembles dRdWT (transpose), so we'll need to
        ! handle transposition when converting to Python
        call setupStateResidualMatrix(dRdWT, useAD, usePC, useTranspose, &
                                      useObjective, frozenTurb, level)

    end subroutine setupResolventMatrix


    subroutine getResolventMatrixInfo(nRows, nCols, nnz)
        !
        ! Get information about the Jacobian matrix for resolvent analysis.
        !
        ! Parameters:
        ! -----------
        ! nRows : integer, intent(out)
        !     Number of rows in the matrix
        ! nCols : integer, intent(out)
        !     Number of columns in the matrix
        ! nnz : integer, intent(out)
        !     Number of non-zeros (approximate)
        !
        implicit none

        ! Output
        integer(kind=intType), intent(out) :: nRows, nCols, nnz

        ! Local
        PetscErrorCode :: ierr
        MatInfo :: info

        ! Get matrix dimensions
        call MatGetSize(dRdWT, nRows, nCols, ierr)

        ! Get matrix info (including number of non-zeros)
        call MatGetInfo(dRdWT, MAT_GLOBAL_SUM, info, ierr)
        nnz = int(info%nz_used, kind=intType)

    end subroutine getResolventMatrixInfo


    subroutine getResolventMatrixDense(J_array, n)
        !
        ! Export the Jacobian matrix as a dense array.
        !
        ! IMPORTANT: This returns J = ∂R/∂w (NOT the transpose!)
        ! ADflow stores dRdWT = (∂R/∂w)^T, so we transpose it.
        !
        ! For resolvent analysis: R(ω) = (jω*I - J)^{-1}
        ! where J = ∂R/∂w is the Jacobian of residuals w.r.t. states
        !
        ! WARNING: Only use for small problems! Very memory intensive.
        !
        ! Parameters:
        ! -----------
        ! n : integer, intent(in)
        !     Size of the matrix (must match matrix size)
        ! J_array : real(n, n), intent(out)
        !     Dense array to store J = ∂R/∂w
        !
        use communication, only: myid

        implicit none

        ! Input/Output
        integer(kind=intType), intent(in) :: n
        real(kind=realType), dimension(n, n), intent(out) :: J_array

        ! Local variables
        PetscErrorCode :: ierr
        Mat :: J_dense
        integer(kind=intType) :: nRows, nCols, i, j
        PetscScalar, pointer :: array(:, :)

        ! Get matrix size
        call MatGetSize(dRdWT, nRows, nCols, ierr)

        if (nRows /= n .or. nCols /= n) then
            if (myid == 0) then
                print *, "ERROR: Matrix size mismatch in getResolventMatrixDense"
                print *, "  Expected:", n, "x", n
                print *, "  Got:", nRows, "x", nCols
            end if
            return
        end if

        ! Convert sparse matrix to dense
        ! Note: This is very memory intensive!
        call MatConvert(dRdWT, MATDENSE, MAT_INITIAL_MATRIX, J_dense, ierr)

        ! Get array from dense matrix
        call MatDenseGetArrayF90(J_dense, array, ierr)

        ! CRITICAL: dRdWT stores (∂R/∂w)^T, so we need to transpose
        ! to get J = ∂R/∂w for resolvent analysis
        do i = 1, n
            do j = 1, n
                J_array(i, j) = array(j, i)  ! Transpose
            end do
        end do

        ! Restore array and destroy dense matrix
        call MatDenseRestoreArrayF90(J_dense, array, ierr)
        call MatDestroy(J_dense, ierr)

    end subroutine getResolventMatrixDense


    subroutine exportResolventMatrixToFile(filename)
        !
        ! Export the Jacobian matrix to a binary file.
        ! This can be loaded in Python using PETSc.
        !
        ! Parameters:
        ! -----------
        ! filename : character(*), intent(in)
        !     Filename to save the matrix
        !
        implicit none

        ! Input
        character(len=*), intent(in) :: filename

        ! Local
        PetscErrorCode :: ierr
        PetscViewer :: viewer

        ! Create binary viewer
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), &
                                   FILE_MODE_WRITE, viewer, ierr)

        ! Write matrix (transpose to get dRdW from dRdWT)
        call MatView(dRdWT, viewer, ierr)

        ! Destroy viewer
        call PetscViewerDestroy(viewer, ierr)

    end subroutine exportResolventMatrixToFile

end module resolventAPI
