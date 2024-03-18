! *********************************************************************************************************************************************************************************************
MODULE matrix_solver
! *********************************************************************************************************************************************************************************************
INTERFACE
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for tridiagonal dot product
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE SUBROUTINE dot_tridiag(a,b,c,d,x)
	! a --- Lower Diagonal (Size of n-1)
	! b --- Main  Diagonal (Size of n  )
	! c --- Upper Diagonal (Size of n-1)
	! d --- Column vector multiplied to matrix (Size of n)
	! x --- Dot product result
	IMPLICIT NONE
    
    ! Data Dicitonary --- Inteface Variables
	REAL(8), INTENT(IN), DIMENSION(2:) :: a
	REAL(8), INTENT(IN), DIMENSION(:)  :: b
	REAL(8), INTENT(IN), DIMENSION(:)  :: c
	REAL(8), INTENT(IN), DIMENSION(:)  :: d
	REAL(8), INTENT(OUT),DIMENSION(:)  :: x
    
END SUBROUTINE dot_tridiag
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for tridiagonal solver
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE SUBROUTINE solve_tridiag(a,b,c,d,x)
	! a --- Lower Diagonal (Size of n-1)
	! b --- Main  Diagonal (Size of n  )
	! c --- Upper Diagonal (Size of n-1)
	! d --- RHS            (Size of n  )
	! x --- sol
	IMPLICIT NONE
	
	! Data Dicitonary --- Inteface Variables
	REAL(8), INTENT(IN), DIMENSION(2:) :: a
	REAL(8), INTENT(IN), DIMENSION(:)  :: b
	REAL(8), INTENT(IN), DIMENSION(:)  :: c
	REAL(8), INTENT(IN), DIMENSION(:)  :: d
	REAL(8), INTENT(OUT),DIMENSION(:)  :: x

END SUBROUTINE solve_tridiag
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for GE solver
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE SUBROUTINE solve_GE(A,b,x)
	IMPLICIT NONE
	
	! Data Dictionary --- Interface Variables
	REAL(8),INTENT(IN), DIMENSION(:,:) :: A
	REAL(8),INTENT(IN), DIMENSION(:)   :: b
	REAL(8),INTENT(OUT),DIMENSION(:)   :: x

END SUBROUTINE solve_GE
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for GE solver with pivoting
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE SUBROUTINE solve_GE_piv(A,b,x)
	IMPLICIT NONE
	
	! Data Dictionary --- Interface Variables
	REAL(8),INTENT(IN), DIMENSION(:,:) :: A
	REAL(8),INTENT(IN), DIMENSION(:)   :: b
	REAL(8),INTENT(OUT),DIMENSION(:)   :: x

END SUBROUTINE solve_GE_piv
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for sawpping rows of given matrix
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE SUBROUTINE SwapRows(v,i,j,nx,ny)
    IMPLICIT NONE
    
    ! Data Dictionary --- Interface Variables
    INTEGER, INTENT(IN) :: i, j, nx, ny
    REAL(8), INTENT(INOUT), DIMENSION(nx,ny) :: v
    
END SUBROUTINE SwapRows
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for locating the maximum value (including index)
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE SUBROUTINE max_val_loc(arr,val_max,inx_max)
    IMPLICIT NONE
    
    ! Data Dictionary --- Interface Variables
    REAL(8), DIMENSION(:), INTENT(IN) :: arr
    REAL(8), INTENT(OUT) :: val_max
    INTEGER, INTENT(OUT) :: inx_max
    
END SUBROUTINE max_val_loc
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
END INTERFACE
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
END MODULE matrix_solver
! *********************************************************************************************************************************************************************************************

! *********************************************************************************************************************************************************************************************
SUBMODULE(matrix_solver) matrix_solver_exe
! *********************************************************************************************************************************************************************************************
CONTAINS
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for tridiagonal dot product
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE PROCEDURE dot_tridiag
    ! Data Dictionary --- Local Variables
	IMPLICIT NONE
    INTEGER :: i
    INTEGER :: n
    ! Get the size
    n = SIZE(b)
    ! Check conformability
	IF(n-1 /= SIZE(a) .OR. n-1 /= SIZE(c) .OR. n /= SIZE(d)) THEN
		WRITE(*,'(A)') 'Invalid Dimension (tridiag)'
		STOP
	END IF
    ! Compute the dot product
    DO i = 1,n
        IF(i == 1) THEN
            x(i) =               b(i)*d(i) + c(i)*d(i+1)
        ELSE IF(i == n) THEN
            x(n) = a(n)*d(i-1) + b(i)*d(i)
        ELSE
            x(i) = a(i)*d(i-1) + b(i)*d(i) + c(i)*d(i+1)
        END IF
    END DO
    
END PROCEDURE dot_tridiag
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for tridiagonal solver
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE PROCEDURE solve_tridiag
	! Data Dictionary --- Local Variables
	IMPLICIT NONE
	INTEGER :: i
	INTEGER :: n
	REAL(8), DIMENSION(:), ALLOCATABLE :: c_prime
	REAL(8), DIMENSION(:), ALLOCATABLE :: d_prime
	REAL(8) :: m
	! Get the Size
	n = SIZE(b)
	! Check conformability
	IF(n-1 /= SIZE(a) .OR. n-1 /= SIZE(c) .OR. n /= SIZE(d)) THEN
		WRITE(*,'(A)') 'Invalid Dimension (tridiag)'
		STOP
	END IF
	! Allocate the shape
	ALLOCATE(c_prime(n-1))
	ALLOCATE(d_prime(n))
	! Initialize c_prime & d_prime 
	c_prime(1) = c(1)/b(1) 
	d_prime(1) = d(1)/b(1)
	! Forward Substitution 
	DO i = 2, n-1
		m = b(i) - c_prime(i-1)*a(i)
		c_prime(i) = c(i)/m
		d_prime(i) = (d(i)-d_prime(i-1)*a(i))/m
	END DO
	m = b(n) - c_prime(n-1)*a(n) 
	d_prime(n) = (d(n)-d_prime(n-1)*a(n))/m
	! Initialize x 
	x(n) = d_prime(n)
	! Backward Substitution 
	DO i = n-1, 1, -1
		x(i) = d_prime(i) - c_prime(i)*x(i+1) 
	END DO

END PROCEDURE solve_tridiag 
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for GE solver
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE PROCEDURE solve_GE
	! Data Dictionary --- Local Variables
	IMPLICIT NONE
	INTEGER :: n
	INTEGER :: i, j, k
	REAL(8), DIMENSION(:,:), ALLOCATABLE :: mat_A
	REAL(8), DIMENSION(:),   ALLOCATABLE :: vec_b
	REAL(8) :: sum, mult
	! Get the size
	n = SIZE(A, DIM = 1)
	! Check conformability
	IF(n /= SIZE(b)) THEN
		WRITE(*,*) 'Invalid Dimension (GE)'
		STOP
	END IF
	! Allocate the size
	ALLOCATE(mat_A(n,n))
	ALLOCATE(vec_b(n))
	! Forward Substitution
	mat_A = A 
	vec_b = b
	DO k = 1, n-1 
		DO i = k+1, n
			mult = mat_A(i,k)/mat_A(k,k)
			mat_A(i,k) = mult
			DO j = k+1, n
				mat_A(i,j) = mat_A(i,j) - mult*mat_A(k,j)
			END DO
			vec_b(i)=vec_b(i)-mult*vec_b(k)
		END DO
	END DO
	! Initialize x
	x(n) = vec_b(n)/mat_A(n,n)
	! Backward Substitution
	DO i=n-1,1,-1 
		sum=vec_b(i)
		DO j=i+1,n 
			sum=sum-mat_A(i,j)*x(j)
		END DO
		x(i)=sum/mat_A(i,i) 
	END DO
	
END PROCEDURE solve_GE
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for GE solver with pivoting
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE PROCEDURE solve_GE_piv
	IMPLICIT NONE
	! Data Dictionary --- Local Variables
    INTEGER :: n
    INTEGER :: i,j,k,p,h
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: my_MAT, my_A
    REAL(8), DIMENSION(:)  , ALLOCATABLE :: vec_b, vec_x, vec_s
    REAL(8) :: max_A, c
    REAL(8), PARAMETER :: eps = 1.0E-6
    
    n = SIZE(A,DIM = 1)
    ALLOCATE(my_MAT(n,n)); ALLOCATE(my_A(n,n))
    ALLOCATE(vec_b(n));    ALLOCATE(vec_x(n)); ALLOCATE(vec_s(n))
    
    my_MAT = A
    vec_b  = b
    vec_x  = x
    
    my_A = my_MAT
    vec_s = 0.0
    p = 0
    
    ! Scale factor array
    DO i = 1,n
        vec_s(i) = MAXVAL(ABS(my_A(i,:)))
    END DO
    
    DO k =1, n-1
        ! Exchage rows if necessary
        CALL max_val_loc(ABS(my_A(k:n,k)/vec_s(k:n)),max_A,p)
        p = p + k - 1
        IF(max_A .LT. eps) THEN
            WRITE(*,'(A)') 'Matrix is Singluar'
            STOP
        END IF
        If(p .NE. k) THEN
            CALL SwapRows(vec_b,k,p,n,1)
            CALL SwapRows(vec_s,k,p,n,1)
            CALL SwapRows(my_A ,k,p,n,n)
        END IF
        ! Elimination
        DO i = k+1,n
            IF(my_A(i,k) .NE. 0.0) THEN
                c = my_A(i,k) / my_A(k,k)
                my_A(i,k+1:n) = my_A(i,k+1:n) - c * my_A(k,k+1:n)
                vec_b(i) = vec_b(i) - c * vec_b(k)
            END IF
        END DO
    END DO
    
    ! Back Substitution
    DO k = n,1,-1
        vec_b(k) = (vec_b(k) - DOT_PRODUCT(my_A(k,k+1:n),vec_b(k+1:n))) / my_A(k,k)
    END DO
    vec_x = vec_b
    x = vec_x
    
END PROCEDURE solve_GE_piv
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for sawpping rows of given matrix
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE PROCEDURE SwapRows
    IMPLICIT NONE
    ! Data Dictionary --- Local Variables
    INTEGER :: kk
    REAL(8), DIMENSION(ny) :: temp_row
    DO kk=1,ny
        temp_row(kk)=v(i,kk);
        v(i,kk)=v(j,kk)
        v(j,kk)=temp_row(kk)
    END DO
    
END PROCEDURE SwapRows
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE for locating the maximum value (including index)
! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
MODULE PROCEDURE max_val_loc
    IMPLICIT NONE  
    ! Data Dictionary --- Local Variables
    INTEGER :: n
    INTEGER :: i, inx
    REAL(8) :: tmp
    ! Initial condition
    n = SIZE(arr)
    inx = 1
    tmp = arr(inx)
    ! SEARCH
    DO i = 2,n
        IF(arr(i) .GT. tmp) THEN
            inx = i
            tmp = arr(i)
        END IF
    END DO
    ! Return
    val_max = tmp
    inx_max = inx
    
END PROCEDURE max_val_loc
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
END SUBMODULE matrix_solver_exe
! *********************************************************************************************************************************************************************************************