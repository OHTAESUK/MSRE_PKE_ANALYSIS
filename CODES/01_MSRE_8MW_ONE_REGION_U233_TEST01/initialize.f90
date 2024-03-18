SUBROUTINE initialize
    USE variables
    USE matrix_solver
    IMPLICIT NONE
    INTEGER :: i
    REAL(8) :: A(3,3)
    REAL(8) :: b(3)
    REAL(8) :: x(3)
    ! *** SIMULATION TIME RELATED
    n_step    = FLOOR(time /del_t + 1.0E-9) + 1
    num_delay = FLOOR(tau_L/del_t + 1.0E-9)
    ! *** DELAYED NEUTRON RELATED
    ALLOCATE(d_lbda(max_dg)); d_lbda = [1.260E-02, 3.370E-02, 1.390E-01, 3.250E-01, 1.130E+00, 2.500E+00]
    ALLOCATE(d_beta(max_dg)); d_beta = [0.00023,0.00079,0.00067,0.00073,0.00013,0.00009]
    lbda = SUM(d_lbda)
    beta = SUM(d_beta)
    ! *** FLOW RELATED
    m_f  = W_f * tau_c
    m_f1 = m_f/2.0d0
    m_f2 = m_f/2.0d0
    m_g  = v_g*rho_g
    k_f1 = k_f/2.0d0
    k_f2 = k_f/2.0d0
    ! *** TEMPERATURE RELATED (SOLVE 3x3 matrix equation)
    A(:,:) = 0.d0
    b(:)   = 0.d0
    A(1,1) = -W_f/m_f1 - (k_1/(k_1+k_2)*hA_fg/(m_f1*cp_f))
    A(1,2) = +0.d0
    A(1,3) =           + (k_1/(k_1+k_2)*hA_fg/(m_f1*cp_f))
    A(2,1) = +W_f/m_f2
    A(2,2) = -W_f/m_f2 - (k_2/(k_1+k_2)*hA_fg/(m_f2*cp_f))
    A(2,3) =           + (k_2/(k_1+k_2)*hA_fg/(m_f2*cp_f))
    A(3,1) = 0.5d0*hA_fg/(m_g*cp_g)
    A(3,2) = 0.5d0*hA_fg/(m_g*cp_g)
    A(3,3) =      -hA_fg/(m_g*cp_g)
    b(1) = -k_f1*P_RX0/(m_f1*cp_f) - (W_f/m_f1*T_fin)
    b(2) = -k_f2*P_RX0/(m_f2*cp_f)
    b(3) = -k_g *P_RX0/(m_g*cp_g)
    CALL solve_GE_piv(A,b,x)
    T_f1 = x(1)
    T_f2 = x(2)
    T_g  = x(3)  
    
    ! (MUTED) T_f1 = T_fin + (T_f2 - T_fin)/2.0d0
    ! (MUTED) T_g  = T_f1  + (k_g*P_RX0/hA_fg);
    
    T_f  = 0.5d0 * (T_f1 + T_f2)  
    T_g0 = T_g
    T_f0 = T_f
    T_f10 = T_f1
    T_f20 = T_f2
    T_fin0 = T_fin
    dT_f  = T_f2 - T_fin
    ! >>> STORES PREVIOUS TEMPERATURE
    T_f1_prev = T_f1
    T_f2_prev = T_f2
    T_g_prev  = T_g
    ! >>> STORES T_f2 (for T_fin)    
    ALLOCATE(arr_T_f2(num_delay)); arr_T_f2 = T_f2
    ! *** REACTIVITY RELATED
    ALLOCATE(sol_PKE      (max_dg + 1))
    ALLOCATE(sol_PKE_prev (max_dg + 1))
    ALLOCATE(sol_PKE_temp (max_dg + 1))
    ALLOCATE(sol_PKE_delay(max_dg + 1))
    ALLOCATE(mat_sol_PKE  (max_dg + 1,num_delay))
    Rho_dyn0 = beta
    DO i = 1,max_dg
        Rho_dyn0 = Rho_dyn0 - d_beta(i)/(1.0d0 + (1.0d0/(d_lbda(i)*tau_c))*(1.0d0 - EXP(-d_lbda(i)*tau_L)))
    END DO
    Rho_dyn      = Rho_dyn0
    Rho_dyn_prev = Rho_dyn0 
    ! *** PKE RELATED
    sol_PKE(1) = 1.0d0
    DO i = 1,max_dg
        sol_PKE(i+1) = (d_beta(i)/LAMBDA) * (1.0d0/(d_lbda(i) - (EXP(-d_lbda(i)*tau_L) - 1.0d0)/tau_c))
    END DO
    sol_PKE_prev  = sol_PKE
    sol_PKE_temp  = sol_PKE
    sol_PKE_delay = sol_PKE
    DO i = 1,num_delay
        mat_sol_PKE(:,i) = sol_PKE
    END DO
    ! *** MISCELLENIOUS
    P_RX      = P_RX0
    P_RX_prev = P_RX0
END SUBROUTINE initialize