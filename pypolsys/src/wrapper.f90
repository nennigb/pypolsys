! This file is part of pypolsys, A python wrapper to `POLSYS_PLP`
! fortran90 package from Layne T. Watson, Steven M. Wise,
! Andrew J. Sommese, August, 1998.
! pypolsys is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! pypolsys is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with pypolsys.  If not, see <https://www.gnu.org/licenses/>.

! ======================================================================
!                 BASIC POLSYS_PLP WRAPPER to python
! ======================================================================
! Based on main program and user written subroutine provided with
! the POLSYS_PLP package.  Layne T. Watson, Steven M. Wise, Andrew
! J. Sommese, August, 1998.  Cosmetic changes, 10/1999.

!=======================================================================
!=======================================================================
module polystats
! This moduls contains common 'stats' on the input polynomial
! These variable are usefull in DHORNER_MV subroutine (avoid to compute
! them at each iteration)
implicit none
integer, allocatable, dimension(:,:) :: vdmax
integer, allocatable, dimension(:,:) :: coefsteps

contains

subroutine init_steps(deg, steps)
! Create the steps to walk in the coef tree for Horner scheme.
!
! Parameters
! ----------
! deg, 2D array
!    The maximal degree for each variable, for each equation (for Horner)
!
! Returns
! -------
! steps, 2D array
! 	The steps to walk in the coef tree (for Horner)
implicit none
integer, dimension(: ,:), intent(in) :: deg
integer, dimension(size(deg,1), size(deg, 2) + 1), intent(out) :: steps

integer :: i, N

N = size(deg, 2)
steps(:, N+1) = 1
do i=N, 1, -1
  steps(:, i) = steps(:, i+1) * deg(:, i)
end do

end subroutine init_steps

end module polystats

!=======================================================================
!=======================================================================
module polsyswrap
! This module will be import in python as polsys.
! All subroutines and variables here will be availlable in python.

USE POLSYS
use polystats

implicit none
! See in `solve` for more description.
! f2py works better with allocatable than with pointers.
COMPLEX (KIND=R8), allocatable, dimension(:,:) :: myROOTS
real (KIND=R8), allocatable, dimension(:) :: variables_SCALE_FACTORS
integer, allocatable, dimension(:) :: path_status, num_jac_eval
integer :: solve_status

contains

!=======================================================================
subroutine init_poly(N, n_coef_per_eq, all_coef, all_deg)
! Init the polynomial object from input.

! Parameters
! ----------
! N : int
!	The number of variable and equations
! n_coef_per_eq :  array
!   The number of terms in each polynomial equation. Sum(n_coef_per_eq) give the
!   total number of coef.
! all_coef : array 
!   Contains the 'list' of all the coefficients, sorting such poly1 coef, poly2 coef...
! all_deg : array
!   Contains the 'list' of all the degree, sorting such poly1 coef, poly2 coef...
!   This is an array such the line number is the terms number and the column contain
!   contains the degree of each variable.
!   ex : if ther 4th term is (3.12+2j) x1^2 x2^3 x3^0, this yields 
!   all_deg(4, :) = [2, 3, 0] and all_coef(4) = 3.12+2j
!
! Attributes
! ----------
! POLYNOMIAL, POLYNOMIAL
!    The polynomial 'object'
! VDMAX, array
!    The maximal degree for each variable, for each equation (useful for DHORNER_MV,
!    if used)
! coefsteps, array
! 	The steps to walk in the coef tree (usefull for DHORNER_MV, if used).
implicit none

integer, intent(in)  :: N
integer, dimension(:), intent(in) :: n_coef_per_eq
COMPLEX (KIND=8), dimension(:), intent(in) :: all_coef
integer, dimension(:,:), intent(in) :: all_deg

integer :: k, i, j, L

! Init global coef counter k
k = 1
! Allocate polynomials
CALL CLEANUP_POL
ALLOCATE(POLYNOMIAL(N))
! Init VDMAX
IF (ALLOCATED(VDMAX)) THEN
  DEALLOCATE(VDMAX)
END IF
ALLOCATE(VDMAX(N, N))
VDMAX = 0

! ALLOCATE THE SHARED POLYNOMIAL OBJECT
DO I=1,N
  POLYNOMIAL(I)%NUM_TERMS = n_coef_per_eq(I)
  ALLOCATE(POLYNOMIAL(I)%TERM(n_coef_per_eq(I)))
  DO J=1, n_coef_per_eq(I)
    ALLOCATE(POLYNOMIAL(I)%TERM(J)%DEG(N+1))
    !POLYNOMIAL(I)%TERM(J)%COEF = COEF(I,J)
    POLYNOMIAL(I)%TERM(J)%COEF = all_coef(k)  
    !POLYNOMIAL(I)%TERM(J)%DEG(1:N) = DEG(I,J,1:N)
    POLYNOMIAL(I)%TERM(J)%DEG(1:N) = all_deg(k, 1:N)
    ! Store Max degree in each variable for each equation
    ! DMAX is used only for horner alg.
    DO L=1, N
      VDMAX(I, L) = MAX(VDMAX(I, L), all_deg(k, L))
    END DO     
    ! count all the coeff
    k = k + 1
  END DO
  VDMAX(I, :) = VDMAX(I, :) + 1
END DO

! init coefsteps (used in DHORNER_MV)
IF (ALLOCATED(coefsteps)) THEN
  DEALLOCATE(coefsteps)
END IF
ALLOCATE(coefsteps(N, N+1))
call init_steps(vdmax, coefsteps)

end subroutine init_poly

!=======================================================================
subroutine init_partition(N, NUM_SETS, NUM_INDICES, INDEX)
! Init the partition object from input.

! Parameters
! ----------
! N : int
!	The number of variable and equations
! NUM_SETS :  array(N)
!   The number of set for each polynomial equation.
! NUM_INDICES : array(N,  NUM_SETS)
!   Contains the 'list' of all the coefficients, sorting such poly1 coef, poly2 coef...
! INDEX : array (N, NUM_SETS, NUM_INDICES)
!   Contains the 'list' of all the variable, in each set for each equation
!
! Attributes
! ----------
! PARTITION, PARTITION
!    The partition 'object'	
implicit none
integer, intent(in)  :: N
integer, dimension(:), intent(in) :: NUM_SETS
integer, dimension(:,:), intent(in) :: NUM_INDICES
integer, dimension(:,:,:), intent(in) :: INDEX
integer :: i, j

! Allocate storage for the system partition in PARTITION. 
CALL CLEANUP_PAR
ALLOCATE(PARTITION_SIZES(N))
PARTITION_SIZES(1:N) = NUM_SETS(1:N)
ALLOCATE(PARTITION(N))
DO I=1,N
  ALLOCATE(PARTITION(I)%SET(PARTITION_SIZES(I)))
  DO J=1,PARTITION_SIZES(I)
    PARTITION(I)%SET(J)%NUM_INDICES = NUM_INDICES(I,J)
    ALLOCATE(PARTITION(I)%SET(J)%INDEX(NUM_INDICES(I,J)))
    PARTITION(I)%SET(J)%INDEX(1:NUM_INDICES(I,J)) = &
      INDEX(I,J,1:NUM_INDICES(I,J))
  END DO 
END DO

end subroutine init_partition

!=======================================================================
subroutine show_coef(homo)
! Show the polynomial representation from fortran.

implicit none
logical, optional :: homo
integer :: I, J, N

! defaut homo=.FALSE.
if (.NOT. present(homo)) then
  homo=.FALSE.
end if

! Check if POLYNOMIAL is allocated
IF (.NOT. ALLOCATED(POLYNOMIAL)) THEN
  write(*,*) 'The polynomial has not been initialiazed. Call `init_poly` before.'
  RETURN
END IF

! Print the POLYNOMIAL coefs
N = size(POLYNOMIAL)
write(*,*) 'Polynomial system with,', N, 'equations and variables.'
DO I=1,N
  write(*,*) 'P', I, ', with ', POLYNOMIAL(I)%NUM_TERMS, ' terms.'
  DO J=1, POLYNOMIAL(I)%NUM_TERMS
    if (homo) then
      ! show also with homogenization variable degree
      write(*,*) POLYNOMIAL(I)%TERM(J)%COEF , ' x**(', POLYNOMIAL(I)%TERM(J)%DEG(1:N+1), ')'
    else
      write(*,*) POLYNOMIAL(I)%TERM(J)%COEF , ' x**(', POLYNOMIAL(I)%TERM(J)%DEG(1:N), ')'
    end if  
  END DO
END DO

end subroutine show_coef


!=======================================================================
subroutine show_partition
! Show the partition representation from fortran.

implicit none
integer :: I, J, N

! Check if PARTITION is allocated
IF (.NOT. ALLOCATED(PARTITION)) THEN
  write(*,*) 'The partition has not been initialiazed. Call `init_partition` before.'
  RETURN
END IF

! Print the PARTITION
N = size(PARTITION)
write(*,*) 'Polynomial system with,', N, 'equations and variables.'
DO I=1,N
  write(*,*) 'Eq. (', I, ') , with ', PARTITION_SIZES(I), ' sets.'
  DO J=1, PARTITION_SIZES(I)
    write(*,*) '  Set #', J,  ', with variables : {', PARTITION(I)%SET(J)%INDEX(1:PARTITION(I)%SET(J)%NUM_INDICES), '}'
  END DO
END DO

end subroutine show_partition

!=======================================================================
subroutine solve(TRACKTOL, FINALTOL, SINGTOL, DENSE, BPLP)
! Call polsys_plp solver.
!
! Parameters
! ----------
! TRACKTOL is the local error tolerance allowed the path tracker along
!   the path.  ABSERR and RELERR (of STEPNX) are set to TRACKTOL.
!
! FINALTOL is the accuracy desired for the final solution.  It is used
!   for both the absolute and relative errors in a mixed error criterion.
!
! SINGTOL is the singularity test threshold used by SINGSYS_PLP.  If
!   SINGTOL <= 0.0 on input, then SINGTOL is reset to a default value.
!
! DENSE (logical, optional) If .True., select the TARGET_SYSTEM_USER 
!   optimized for dense polynomial (horner). If .FALSE. (default) the defaut 
!   POLSYS_PLP TARGET_SYSTEM is used. The default choice is safer. Using Horner
!   requires that `utils.toDense` is called to create the polynomial system.
!
! Returns
! --------
! BPLP is the generalized Bezout number corresponding to the
!   partitioned linear product (PLP) formulation defined by the system
!   partition P.  This is the number of paths tracked and the number of
!   roots returned (counting multiplicity).

! Attributs
! ---------
! MYROOTS(1:N,I) are the complex roots (untransformed and unscaled) of
!   the target polynomial corresonding to the Ith path, for I=1,...,BPLP.
!   ROOTS(N+1,I) is the homogeneous variable of the target polynomial
!   system in complex projective space corresponding to ROOTS(1:N,I).
! path_status :  is an integer array which returns information about
!   each path tracked. Correspond to POLSYS_PLP IFLAG2.
! solve_status : is an integer containing the exist status of POLSYS_PLP.
! num_jac_eval : Array(BPLP) containing the number of Jacobian evaluation 
! variables_scale_factors : Array(N) containing the scale factor.
!   (Needed for `refine`).

implicit none
INTEGER, PARAMETER:: NN = 30
REAL (KIND=R8), INTENT(IN):: TRACKTOL, FINALTOL
REAL (KIND=R8), INTENT(INOUT):: SINGTOL
INTEGER, INTENT(OUT) :: BPLP
LOGICAL, OPTIONAL :: DENSE

INTEGER :: IFLAG1, N
REAL (KIND=R8), DIMENSION(:), POINTER:: ARCLEN, LAMBDA
INTEGER, DIMENSION(:), POINTER:: IFLAG2, NFE
REAL (KIND=R8), DIMENSION(NN):: SCALE_FACTORS
REAL (KIND=R8), DIMENSION(8):: SSPAR
COMPLEX (KIND=R8), DIMENSION(:,:), POINTER:: ROOTS

! Set SSPAR to zeros, activate default value
SSPAR(1:8) = 0.0_R8

! Optional `dense` argument
if (.NOT. present(DENSE)) then
  DENSE=.FALSE.
end if

IF (ALLOCATED(myROOTS)) THEN
  DEALLOCATE(myROOTS)
END IF

IF (.NOT. ALLOCATED(POLYNOMIAL)) THEN
  write(*,*) 'The polynomial has not been initialiazed. Call `init_poly` before.'
  RETURN
ELSE
  N = size(POLYNOMIAL)
END IF

IF (.NOT. ALLOCATED(PARTITION)) THEN
  write(*,*) 'The partition has not been initialiazed. Call `init_partition` before.'
  RETURN
END IF

! Disassociate pointers. From previous computation...
NULLIFY(IFLAG2, NFE, ARCLEN, LAMBDA, ROOTS) 

! Call POLSYS with TARGET_SYSTEM_USER if DENSE
if (DENSE) then
  CALL POLSYS_PLP(N,TRACKTOL,FINALTOL,SINGTOL,SSPAR,BPLP,IFLAG1,IFLAG2, &
        ARCLEN,LAMBDA,ROOTS,NFE,SCALE_FACTORS, USER_F_DF=DENSE)
else
  CALL POLSYS_PLP(N,TRACKTOL,FINALTOL,SINGTOL,SSPAR,BPLP,IFLAG1,IFLAG2, &
        ARCLEN,LAMBDA,ROOTS,NFE,SCALE_FACTORS)
end if

! Export results to module variable
allocate (myROOTS(size(ROOTS,1), size(ROOTS,2)))
myROOTS = ROOTS

! Map usefull parameters for reporting
! ----------------------------------------------------------------------
! Set solve_status
solve_status = IFLAG1
! Set path_status
IF (ALLOCATED(path_status)) THEN
  DEALLOCATE(path_status)
END IF
ALLOCATE(path_status(BPLP))
path_status = IFLAG2

! Set module parameter usefull for refine subroutine
! ----------------------------------------------------------------------
! Set variables_scale_factors
IF (ALLOCATED(variables_scale_factors)) THEN
  DEALLOCATE(variables_scale_factors)
END IF
ALLOCATE(variables_scale_factors(N))
variables_scale_factors = SCALE_FACTORS(1:N)
! Set num_jac_eval
IF (ALLOCATED(num_jac_eval)) THEN
  DEALLOCATE(num_jac_eval)
END IF
ALLOCATE(num_jac_eval(BPLP))
num_jac_eval = NFE

end subroutine solve


!=======================================================================
SUBROUTINE refine(path, TRACKTOL, FINALTOL, SINGTOL, DENSE, K)
! Retrack certain homotopy paths to imporove accuracy.
! A prior call to `solve` is mandatory.
!
! Parameters
! ----------
! PATH is an integer array containing the index of the paths to retrack.
!
! TRACKTOL is the local error tolerance allowed the path tracker along
!   the path.  ABSERR and RELERR (of STEPNX) are set to TRACKTOL.
!
! FINALTOL is the accuracy desired for the final solution.  It is used
!   for both the absolute and relative errors in a mixed error criterion.
!
! SINGTOL is the singularity test threshold used by SINGSYS_PLP.  If
!   SINGTOL <= 0.0 on input, then SINGTOL is reset to a default value.
!
! DENSE (logical, optional) If .True., select the TARGET_SYSTEM_USER 
!   optimized for dense polynomial (horner). If .FALSE. (default) the defaut 
!   POLSYS_PLP TARGET_SYSTEM is used. The default choice is safer. Using Horner
!   requires that `utils.toDense` is called to create the polynomial system.
!
! K is the number of path to retrack, K = size(path).

! Remarks 
! --------
! In `POLSYS_PLP` documentation, the RECALL argument
!   is used to retrack certain homotopy paths.  It's use assumes
!   BPLP contains the Bezout number (which is not recalculated),
!   SCALE_FACTORS contains the variable unscaling factors, and that
!   IFLAG2(1:BPLP) exists.  The Ith homotopy path is retracked if
!   IFLAG2(I) = -2, and skipped otherwise.

implicit none
integer, intent(in) :: K
integer, intent(in) :: path(K)
REAL (KIND=R8), INTENT(IN):: TRACKTOL, FINALTOL
REAL (KIND=R8), INTENT(INOUT):: SINGTOL
LOGICAL, OPTIONAL :: DENSE

integer :: i, bplp, N
INTEGER, DIMENSION(:), POINTER:: IFLAG2, NFE
REAL (KIND=R8), DIMENSION(8):: SSPAR
REAL (KIND=R8), DIMENSION(:), POINTER:: ARCLEN, LAMBDA
COMPLEX (KIND=R8), DIMENSION(:,:), POINTER:: ROOTS
! Optional `dense` argument
if (.NOT. present(DENSE)) then
  DENSE=.FALSE.
end if

! Basic initial checks
IF (.NOT. ALLOCATED(POLYNOMIAL)) THEN
  write(*,*) 'The polynomial has not been initialiazed. Call `init_poly` before.'
  RETURN
ELSE
  N = size(POLYNOMIAL)
END IF

IF (.NOT. ALLOCATED(PARTITION)) THEN
  write(*,*) 'The partition has not been initialiazed. Call `init_partition` before.'
  RETURN
END IF

IF ((.NOT. ALLOCATED(path_status)).OR.&
    (.NOT. ALLOCATED(num_jac_eval)).OR.&
    (.NOT. ALLOCATED(myROOTS))) THEN
  write(*,*) 'No path to retrack. Need to call `solve` before.'
  RETURN
END IF

! Setup required parameters
! ---------------------------------------------------------------------
! Set SSPAR to zeros, activate default value
SSPAR(1:8) = 0.0_R8
! Recover BPLP value
BPLP = size(path_status)
! Disassociate pointers.
NULLIFY(ARCLEN, LAMBDA, ROOTS, IFLAG2, NFE) 
! Recover variables from previous run
ALLOCATE(IFLAG2(BPLP))
ALLOCATE(NFE(BPLP))
IFLAG2 = path_status
NFE = num_jac_eval
! Need to allocate other array
ALLOCATE(ROOTS(N+1,BPLP))
ROOTS = myROOTS
ALLOCATE(LAMBDA(BPLP))
ALLOCATE(ARCLEN(BPLP))

! Change iflag2 to retrack, see POLSYS_PLP doc.
do i=1, K
  IFLAG2(path(i)) = -2
end do

! Now retrack the paths `path` using the RECALL option:
if (DENSE) then
  CALL POLSYS_PLP(N,TRACKTOL,FINALTOL,SINGTOL,SSPAR,BPLP, solve_status,IFLAG2, &
        ARCLEN,LAMBDA,ROOTS,NFE,variables_scale_factors, USER_F_DF=DENSE, RECALL=.TRUE.)
else
  CALL POLSYS_PLP(N,TRACKTOL,FINALTOL,SINGTOL,SSPAR,BPLP,solve_status,IFLAG2, &
        ARCLEN,LAMBDA,ROOTS,NFE,variables_scale_factors, RECALL=.TRUE.)
end if  

! Export refined roots
do i=1, K
  myROOTS(:, path(i)) = ROOTS(:, path(i))
  path_status(path(i)) = IFLAG2(path(i))
  num_jac_eval(path(i)) = NFE(path(i))
end do

END SUBROUTINE refine


!=======================================================================
SUBROUTINE bezout(SINGTOL, BPLP)
! Compute and returns only the generalized Bezout number BPLP 
! of the target polynomial system, based on the variable partition.
!
! Such computation is very quick, which is useful for exploring 
! alternative partitions.
!
! Parameters
! ----------
! SINGTOL is the singularity test threshold used by SINGSYS_PLP.  If
!   TOL <= 0.0 on input, TOL is reset to the default value
!   SQRT(EPSILON(1.0_R8)).
!
! GLOBAL_PLP allocatable objects POLYNOMIAL, PARTITION_SIZES, and
!   PARTITION (see GLOBAL_PLP documentation) must be allocated and
!   defined in the calling program.
!
! Returns
! --------
! BPLP is the generalized Bezout number for the target system based on
!   the variable partition P defined in the module GLOBAL_PLP.


implicit none
integer, intent(out) :: bplp
real (KIND=R8), intent(INOUT):: SINGTOL

integer :: i, N, MAXT

! Initial check
IF (.NOT. ALLOCATED(POLYNOMIAL)) THEN
  write(*,*) 'The polynomial has not been initialiazed. Call `init_poly` before.'
  RETURN
ELSE
  N = size(POLYNOMIAL)
END IF

IF (.NOT. ALLOCATED(PARTITION)) THEN
  write(*,*) 'The partition has not been initialiazed. Call `init_partition` before.'
  RETURN
END IF

! Compute MAXT is the maximum number of terms in any component of the target
!   system.  MAXT = MAX((/(NUMT(I),I=1,N)/)).
MAXT = 0
do I=1, N
  MAXT = MAX(POLYNOMIAL(I)%NUM_TERMS, MAXT)
end do
! Compute the BEsout number wusing polysys_plp
CALL BEZOUT_PLP(N, MAXT, SINGTOL, BPLP)

END SUBROUTINE bezout


!=======================================================================
SUBROUTINE report(status)
! Show Homotopy solve report with comprehensive messages.
!
! Returns
! -------
! status : integer
!   Retruns 0 if all path exit without warning.
implicit none
integer, intent(out) :: status
integer :: i
logical :: warn

warn = .FALSE.
write(*,*) 'Solve report :'
write(*,*) '--------------'
! Mapping of `POLSYS_PLP` warning messages
!   = 0  for a normal return.
!   = -1 if either POLYNOMIAL or PARTITION was improperly allocated.
!   = -2 if any POLYNOMIAL(I)%TERM(J)%DEG(K) is less than zero.
!   = -3 if F_I(X) = CONSTANT for some I.
!   = -4 if SUM_{J=1}^{PARTITION_SIZES(I)} 
!          PARTITION(I)SET(J)%NUM_INDICES /= N, for some I.
!   = -5 if UNION_{J=1}^{PARTITION_SIZES}
!               S(I,J) /= {1,2,...,N-1,N}, for some I.
!   = -6 if the optional argument RECALL was present but any of BPLP
!        or the arrays ARCLEN, IFLAG2, LAMBDA, NFE, ROOTS are
!        inconsistent with the previous call to POLSYS_PLP.
!   = -7 if the array SCALE_FACTORS is too small.

! SOLVE_STATUS Messages
! ----------------------------------------------------------------------
select case (solve_status)
  case (0) 
    print*, '  Normal return.'
  case (-1) 
    print*, '  POLSYS_PLP Warning : Either POLYNOMIAL or PARTITION was improperly allocated.'   
  case (-2) 
    print*, '  POLSYS_PLP Warning : Any POLYNOMIAL(I)%TERM(J)%DEG(K) is less than zero.'
  case (-3) 
    print*, '  POLSYS_PLP Warning : F_I(X) = CONSTANT for some I.'    
  case (-4) 
    print*, '  POLSYS_PLP Warning : SUM_{J=1}^{PARTITION_SIZES(I)} PARTITION(I)SET(J)%NUM_INDICES /= N, for some I.'   
  case (-5) 
    print*, '  POLSYS_PLP Warning :  UNION_{J=1}^{PARTITION_SIZES} S(I,J) /= {1,2,...,N-1,N}, for some I.'
  case (-6) 
    print*, '  POLSYS_PLP Warning : the optional argument RECALL was present but any of'
    print*, '    BPLP, ARCLEN, IFLAG2, LAMBDA, NFE, ROOTS are inconsistent with the previous'
    print*, '    call to POLSYS_PLP.'
  case (-7) 
    print*, '  POLSYS_PLP Warning :  the array SCALE_FACTORS is too small.'    
end select

! IFLAG2(1:BPLP) is an integer array which returns information about
!   each path tracked.  Precisely, for each path I that was tracked,
!   IFLAG2(I):
!   = 1 + 10*C, where C is the cycle number of the path, for a normal return.
!   = 2 if the specified error tolerance could not be met.  Increase
!       TRACKTOL and rerun.
!   = 3 if the maximum number of steps allowed was exceeded.  To track
!       the path further, increase NUMRR and rerun the path.
!   = 4 if the Jacobian matrix does not have full rank.  The algorithm has
!       failed (the zero curve of the homotopy map cannot be followed any
!       further).
!   = 5 if the tracking algorithm has lost the zero curve of the homotopy
!       map and is not making progress.  The error tolerances TRACKTOL and
!       FINALTOL were too lenient.  The problem should be restarted with
!       smaller error tolerances.
!   = 6 if the normal flow Newton iteration in STEPNX or ROOT_PLP failed
!       to converge.  The error error tolerances TRACKTOL or FINALTOL may
!       be too stringent.
!   = 7 if ROOT_PLP failed to find a root in 10*NUMRR iterations.

! PATH_STATUS Messages
! ----------------------------------------------------------------------
do i=1, size(path_status)
  ! Display Warnings messages
  if (path_status(i) < 10 ) then
    warn = .TRUE.
    print *, '  >In path ', i
    select case (path_status(i))
      case (2) 
        print*, '  The specified error tolerance could not be met.  Increase TRACKTOL and rerun.'   
      case (3) 
        print*, '  The maximum number of steps allowed was exceeded.  To track the path further,'
        print*, '  increase NUMRR and rerun the path.'
      case (4) 
        print*, '  The Jacobian matrix does not have full rank.  The algorithm has'
        print*, '    failed (the zero curve of the homotopy map cannot be followed any'
        print*, '    further).'
      case (5) 
        print*, '  The tracking algorithm has lost the zero curve of the homotopy'
        print*, '    map and is not making progress.  The error tolerances TRACKTOL and'
        print*, '    FINALTOL were too lenient.  The problem should be restarted with'
        print*, '    smaller error tolerances.'
      case (6)       
        print*, '   The normal flow Newton iteration in STEPNX or ROOT_PLP failed'
        print*, '    to converge.  The error error tolerances TRACKTOL or FINALTOL may'
        print*, '    be too stringent.'
      case (-7) 
        print*, '   ROOT_PLP failed to find a root in 10*NUMRR iterations.'    
    end select
  end if
end do

! If no warnings, display success message
if (.NOT.warn) then
  print*, '  All paths were normally tracked.'  
  status = 0
else
  !There were some warnings 
  status = 1
end if

END SUBROUTINE report


!=======================================================================
SUBROUTINE CLEANUP_POL
! Deallocates structure POLYNOMIAL.
implicit none
integer :: i, j

IF (.NOT. ALLOCATED(POLYNOMIAL)) RETURN
DO I=1,SIZE(POLYNOMIAL)
  DO J=1,NUMT(I)
    DEALLOCATE(POLYNOMIAL(I)%TERM(J)%DEG)
  END DO
  DEALLOCATE(POLYNOMIAL(I)%TERM)
END DO
DEALLOCATE(POLYNOMIAL)
RETURN
END SUBROUTINE CLEANUP_POL

!=======================================================================
SUBROUTINE CLEANUP_PAR
! Deallocates structure PARTITION.
implicit none
integer :: i, j

IF (.NOT. ALLOCATED(PARTITION)) RETURN
DO I=1,SIZE(PARTITION)
  DO J=1,PARTITION_SIZES(I)
    DEALLOCATE(PARTITION(I)%SET(J)%INDEX)
  END DO
  DEALLOCATE(PARTITION(I)%SET)
END DO
DEALLOCATE(PARTITION)
DEALLOCATE(PARTITION_SIZES)
RETURN
END SUBROUTINE CLEANUP_PAR


end module polsyswrap
!=======================================================================
!=======================================================================

!=======================================================================
SUBROUTINE DHORNER_MV(eq, Ncoef, Nv, x0, p, dp)
USE REAL_PRECISION
use GLOBAL_PLP
use polystats
! Evaluate the polynomial of index `eq` and its derivatives with Horner 
! scheme at x0.

! Instead of working with a matrix of coefficients, the algorithm is express in term of
! degree and coefficients list.
!  
! The Horner coefficient are stored in place and only coef[k] cummulate intermediate results.
! The `coef` vairable is walked as in a tree. Each node has known number of child 
! depending of the degree of the other variable. Each layer of the tree correspond 
! to a variable.
! Use `polystats` module to factor out intermediate variable, like coefsteps,
! computation.
!
! Parameters
! ----------
! eq, int
!   The index of the Polynomial equations to evaluate.
! Ncoef, int
! 	The Number of coeffcients in the polynomial
! Nv, int
!   The number of variable involded
! x0, 1d array
!   The point where the polynomial will be evaluated 

! Returns
! -------
! p, complex
!	  The polynamial value
! dp, array
!	  The derivative with respect to all variables
implicit none
integer, intent(in) :: eq
integer, intent(in) :: Ncoef, Nv
complex(R8), dimension(Nv), intent(in) :: x0
complex(R8), intent(out) :: p
complex(R8), dimension(Nv), intent(out) :: dp

! local variables
complex(R8), dimension(Ncoef) :: coef
complex(R8), dimension(Nv, Ncoef) :: dcoef
integer :: v, k, i, v_, ki

! Init coef with initial Polynomial coefficient
do i = 1, Ncoef
  coef(i) = POLYNOMIAL(eq)%TERM(i)%COEF
end do
dcoef = (0.0_R8, 0.0_R8)
! Horner scheme with derivatives
do v = Nv, 1, -1
    do k = coefsteps(eq, v), Ncoef, coefsteps(eq, v)
        do i = 1, vdmax(eq, v) - 1
            ! Update derivatives tree with respect to v
            dcoef(v, k) = coef(k) + dcoef(v, k)*x0(v)
            ! Update other derivatives (horner steps)
            ki = k - i*coefsteps(eq, v+1)
            do v_ = v+1, Nv
                dcoef(v_, k) = dcoef(v_, ki) + dcoef(v_, k)*x0(v)
            end do ! v_
            ! Update horner step for the polynomial tree
            coef(k) = coef(ki) + coef(k)*x0(v)      
        end do ! i
  end do ! k
end do ! v
p = coef(Ncoef)
dp = dcoef(:, Ncoef)
    
END SUBROUTINE DHORNER_MV

!=======================================================================
SUBROUTINE TARGET_SYSTEM_USER(N,PROJ_COEF,XC,F,DF)
! Provide an Horner based evaluation the polynomial system and
! it jacobian matrix. 
!
! To select this method, call `solve` with `dense=True` and be sure that
! all monomials are present (use `utils.toDense` method for instance).
!
! These method is *optimized* for speed on dense polynomial systems.
! To speed up the evaluation, evaluation is carried out only with z 
! coordinante (skip the homogenisation variable). Since XC(N+1) may
! tend to 0 for solution at infinty, it is possible that numerical error
! may increase for this cases.
! The interface must be the same as `TARGET_SYSTEM` from `POLSYS_PLP`.

! Parameters
! ----------
! N : int
!	The number of variable/equations
! XC : array (1:N+1)
! 	Is in complex projective coordinates, and the homogeneous coordinate
! 	XC(N+1) is explicitly eliminated from F(XC) and DF(XC) using the
! 	projective transformation (cf. the comments in START_POINTS_PLP).
!   XC(N+1) is the homogeneous coordinate 
!
! Returns
! -------
! F, array(N)
!	The system value at XC (complex).
! DF, 2D array NxN
!   The System jacobian matrix. DF(:,N+1) is not referenced by the
!   calling program.
!	
! Remarks
! -------
! The comments in the internal subroutine `TARGET_SYSTEM` should be read before
! attempting to write this subroutine; pay particular attention to the
! handling of the homogeneous coordinate XC(N+1).  

USE REAL_PRECISION
USE GLOBAL_PLP
IMPLICIT NONE
INTEGER, INTENT(IN):: N
COMPLEX (KIND=R8), INTENT(IN), DIMENSION(N+1):: PROJ_COEF, XC
COMPLEX (KIND=R8), INTENT(OUT), DIMENSION(N) :: F(N)
COMPLEX (KIND=R8), INTENT(OUT), DIMENSION(N, N+1) :: DF
COMPLEX (KIND=R8), DIMENSION(N) :: X, dp

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
INTEGER:: I, L
COMPLEX (KIND=R8):: gam, S, p
INTEGER, DIMENSION(N) :: DI

DF = (0.0_R8,0.0_R8)
! May be move elsewhere, call once for all
DO I=1, N
  DI(I) = SUM(POLYNOMIAL(I)%TERM(1)%DEG(:))
END DO
! Recover original variable. Possible Behave not well since XC(N+1) may tend to
! 0 for solution at infinty.
X = XC(1:N) / XC(N+1)
! May use openmp for this loop
DO I=1,N
  call DHORNER_MV(I, POLYNOMIAL(I)%NUM_TERMS, N, X, p, dp)
  ! Store Polynomial value
  F(I) = p
  ! Store Jacobian value
  DF(I,1:N) = dp
  
  ! Special treatement for J=N+1
  ! F' = w_(N+1)^di * F(z1, z2, ...) = w_(N+1)^di * F(w1/w_(N+1), w2/w_(N+1), ...)
  ! zj = wj/w_(N+1)
  S = (0.0_R8,0.0_R8)
  DO L=1,N
    S = S + DF(I, L) * XC(L)
  END DO
  ! factor XC(N+1)**(DI(I) - 1) !
  DF(I, N+1) = F(I) * DI(I) * XC(N+1)**(DI(I) - 1) - S * XC(N+1)**(DI(I)  - 2)
  
  ! Reapply transform
  gam = XC(N+1)**DI(I) 
  F(I) = F(I) * gam
  ! dF'/d wj = w_(N+1)^di * dF(z1, z2, ...) / dwj = w_(N+1)^di * dF(z1, z2, ...) /dz * dz/dwj
  DF(I, 1:N) = DF(I, 1:N) * gam / XC(N+1)  
  ! Apply chain rule
  DF(I,1:N) = DF(I,1:N) + PROJ_COEF(1:N) * DF(I,N+1)

END DO

RETURN

END SUBROUTINE TARGET_SYSTEM_USER

