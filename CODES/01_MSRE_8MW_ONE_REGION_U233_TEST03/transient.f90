SUBROUTINE transient
    USE variables
    USE matrix_solver
    IMPLICIT NONE
    ! +++ INTERFACE
    INTERFACE
        SUBROUTINE solve_PKE_MSR(rho_in,sol_delay,sol_prev,sol_in)
            USE variables
            USE matrix_solver
            IMPLICIT NONE
            REAL(8), INTENT(IN)     :: rho_in
            REAL(8), INTENT(IN)     :: sol_delay(:)
            REAL(8), INTENT(IN)     :: sol_prev (:)
            REAL(8), INTENT(INOUT)  :: sol_in   (:)
        END SUBROUTINE solve_PKE_MSR
    END INTERFACE
    INTEGER :: inx_t,i,j
    REAL(8) :: A(2,2)
    REAL(8) :: b(2)
    REAL(8) :: x(2)
    OPEN(UNIT = 1,FILE = 'output.txt', STATUS = 'REPLACE', ACTION = 'WRITE')
    WRITE(1,'(I6,99ES15.5)') 0, 0.d0, T_fin,T_f1,T_f2,T_g,sol_PKE(1),P_RX
    time_march: DO inx_t = 1,n_step
        ! =====================================================================================
        ! [1] APPEND the previous (time-step) information
        ! =====================================================================================
        T_f1_prev           = T_f1
        T_f2_prev           = T_f2
        T_g_prev            = T_g
        T_pr_prev           = T_pr
        T_IHX_salt_avg_prev = T_IHX_salt_avg
        T_sr_prev           = T_sr
        sol_PKE_prev        = sol_PKE
        Rho_dyn_prev        = Rho_dyn
        ! =====================================================================================
        ! [2] Update inlet temperature
        ! =====================================================================================
        T_fin = T_out_IHX_bank(1)
        T_out_IHX_bank(1:bsize_Lc-1) = T_out_IHX_bank(2:bsize_Lc)
        T_out_IHX_bank(  bsize_Lc)   = 0.d0
        ! =====================================================================================
        ! [3] Update delayed PKE soluition
        ! =====================================================================================
        sol_PKE_delay = mat_sol_PKE(:,1)
        DO i = 1,bsize_L-1
            mat_sol_PKE(:,i) = mat_sol_PKE(:,i+1)
        END DO
        mat_sol_PKE(:,bsize_L) = 0.d0
        ! =====================================================================================
        ! [4] Core Calculation 
        ! =====================================================================================
        DO
            ! Temprary Storage + Restore PKE solution
            T_f1_temp = T_f1
            T_f2_temp = T_f2
            T_g_temp  = T_g
            sol_PKE_temp = sol_PKE
            sol_PKE      = sol_PKE_prev
            ! Solving PKE equation
            Rho_dyn = Rho_dyn0 + a_f*((T_f1 + T_f2)/2.d0-T_f0) + a_g*(T_g-T_g0) + Rho_ext
            CALL solve_PKE_MSR(Rho_dyn,sol_PKE_delay,sol_PKE_prev,sol_PKE)
            ! Solving Balance equation for Temperatures within the core
            P_RX = P_RX0 * sol_PKE(1)
            CALL update_core_temperature
            ! Check the convergence
            IF((ABS(sol_PKE(1) - sol_PKE_temp(1))/sol_PKE(1) < 1.0E-9) .AND. (ABS(T_f2 - T_f2_temp)/T_f2 < 1.0E-9)) THEN
                WRITE(1,'(I6,99ES15.5)') inx_t, inx_t*del_t, T_fin, T_f1, T_f2, T_g, sol_PKE(1),P_RX
                EXIT
            END IF
        END DO
        mat_sol_PKE  (:,bsize_L)  = sol_PKE
        T_out_RX_bank(2:bsize_Lh) = T_out_RX_bank(1:bsize_Lh-1)
        T_out_RX_bank(1) = T_f2
        ! =====================================================================================
        ! HEAT BALANCE AT IHX
        ! =====================================================================================
        T_IHX_in  = T_out_RX_bank(bsize_Lh)
        T_IHX_salt_in = T_out_SHX_SALT_bank(1); T_out_SHX_SALT_bank(1:bsize_r_HX-1) = T_out_SHX_SALT_bank(2:bsize_r_HX)
        
        A(1,1) = (M_IHX*cp_f/del_t + 2.0d0*m_dot_IHX*cp_f + GAMMA_IHX)
        A(1,2) = -GAMMA_IHX
        A(2,1) = -GAMMA_IHX
        A(2,2) = (M_IHX_salt*cp_HX/del_t + 2.0d0*m_dot_SHX*cp_HX + GAMMA_IHX)
        
        b(1) = (M_IHX*cp_f/del_t) * T_pr_prev + 2.0d0*m_dot_IHX*cp_f*T_IHX_in
        b(2) = (M_IHX_salt*cp_HX/del_t) * T_IHX_salt_avg_prev + 2.0d0*m_dot_SHX*cp_HX*T_IHX_salt_in
        
        CALL solve_GE_piv(A,b,x)
        
        T_pr = x(1)
        T_IHX_out = 2.0d0*T_pr - T_IHX_in
        T_IHX_salt_avg = x(2)
        T_IHX_salt_out = 2.0d0*T_IHX_salt_avg - T_IHX_salt_in
        
        T_out_IHX_bank(bsize_Lc) = T_IHX_out
        ! =====================================================================================
        ! HEAT BALANCE AT SHX (RADIATOR FOR MSRE)
        ! =====================================================================================
        T_SHX_in = T_out_IHX_SALT_bank(bsize_HX_r)
        T_out_IHX_SALT_bank(2:bsize_HX_r) = T_out_IHX_SALT_bank(1:bsize_HX_r-1)
        T_out_IHX_SALT_bank(1) = T_IHX_salt_out        
        
        T_sr     = (2.0d0*m_dot_SHX*del_t/M_SHX*T_SHX_in + T_sr_prev + GAMMA_SHX*del_t/M_SHX/cp_HX*T_sink) / &
                   (1.0d0+GAMMA_SHX*del_t/M_SHX/cp_HX    + 2.0d0*m_dot_SHX*del_t/M_SHX)
        T_SHX_out = 2.0d0*T_sr - T_SHX_in
        T_out_SHX_SALT_bank(bsize_r_HX) = T_SHX_out
        
    END DO time_march
END SUBROUTINE transient

! +++ SUBROUTINE FOR UPDATING TEMPERATURE
SUBROUTINE update_core_temperature
    USE variables
    USE matrix_solver
    IMPLICIT NONE
    REAL(8) :: A(3,3)
    REAL(8) :: b(3)
    REAL(8) :: x(3)
    ! --- LHS
    A(1,1) = 1.0d0/del_t + W_f/m_f1 + (k_1/(k_1+k_2)) * hA_fg / (m_f1*cp_f)
    A(1,2) = 0.0d0
    A(1,3) =                        - (k_1/(k_1+k_2)) * hA_fg / (m_f1*cp_f)
    A(2,1) = - W_f/m_f2
    A(2,2) = 1.0d0/del_t + W_f/m_f2 + (k_2/(k_1+k_2)) * hA_fg / (m_f2*cp_f)
    A(2,3) =                        - (k_2/(k_1+k_2)) * hA_fg / (m_f2*cp_f)
    A(3,1) =      -0.5d0 * hA_fg / (m_g*cp_g)
    A(3,2) =      -0.5d0 * hA_fg / (m_g*cp_g)
    A(3,3) = 1.0d0/del_t + hA_fg / (m_g*cp_g)
    ! --- RHS
    b(1) = T_fin * (W_f/m_f1) + T_f1_prev * 1.0d0/del_t + k_f1*P_RX/(m_f1*cp_f)
    b(2) =                      T_f2_prev * 1.0d0/del_t + k_f2*P_RX/(m_f2*cp_f)
    b(3) =                      T_g_prev  * 1.0d0/del_t + k_g *P_RX/(m_g *cp_g)
    ! --- SOLVE
    CALL solve_GE_piv(A,b,x)
    T_f1 = x(1)
    T_f2 = x(2)
    T_g  = x(3)
END SUBROUTINE update_core_temperature

! +++ SUBROUTINE FOR SOLVING FLOW INCLUDED PKE
SUBROUTINE solve_PKE_MSR(rho_in,sol_delay,sol_prev,sol_in)
    USE variables
    USE matrix_solver
    IMPLICIT NONE
    REAL(8), INTENT(IN)     :: rho_in
    REAL(8), INTENT(IN)     :: sol_delay(:)
    REAL(8), INTENT(IN)     :: sol_prev (:)
    REAL(8), INTENT(INOUT)  :: sol_in   (:)
    REAL(8), ALLOCATABLE    :: A(:,:)
    REAL(8), ALLOCATABLE    :: b(:)
    INTEGER :: i,j
    ALLOCATE(A(max_dg+1,max_dg+1))
    ALLOCATE(b(max_dg+1))
    ! --- LHS
    A(:,:) = 0.0d0
    A(1,1) = 1.0d0/del_t - (rho_in - beta)/LAMBDA
    DO i = 2,max_dg+1
        A(1,i) = -d_lbda(i-1)
        A(i,1) = -d_beta(i-1)/LAMBDA
        A(i,i) = 1.0d0/del_t + d_lbda(i-1) + 1.0d0/tau_c
    END DO
    ! --- RHS
    b(1) = sol_prev(1)/del_t
    DO i = 2,max_dg+1
        b(i) = sol_prev(i)/del_t + EXP(-d_lbda(i-1)*tau_L)/tau_c*sol_delay(i)
    END DO
    ! --- SOLVE
    CALL solve_GE_piv(A,b,sol_in)
END SUBROUTINE solve_PKE_MSR