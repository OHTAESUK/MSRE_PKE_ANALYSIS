SUBROUTINE initialize
    USE variables
    USE matrix_solver
    IMPLICIT NONE
    INTEGER :: i,j
    REAL(8) :: A(3,3)
    REAL(8) :: b(3)
    REAL(8) :: x(3)
    REAL(8) :: tmp_x1,tmp_y1
    REAL(8) :: tmp_x2,tmp_y2
    ! *** SIMULATION TIME RELATED
    n_step    = FLOOR(time  /del_t + 1.0E-9)
    bsize_L   = FLOOR(tau_L /del_t + 1.0E-9)
    bsize_Lc  = FLOOR(tau_Lc/del_t + 1.0E-9)
    bsize_Lh  = FLOOR(tau_Lh/del_t + 1.0E-9)
    ! *** DELAYED NEUTRON RELATED
    ALLOCATE(d_lbda(max_dg)); d_lbda = [1.260E-02, 3.370E-02, 1.390E-01, 3.250E-01, 1.130E+00, 2.500E+00]
    ALLOCATE(d_beta(max_dg)); d_beta = [0.0002300, 0.0007900, 0.0006700, 0.0007300, 0.0001300, 0.0000900]
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
    ALLOCATE(T_out_RX_bank(bsize_Lh)); T_out_RX_bank = T_f2
    ALLOCATE(T_out_HX_bank(bsize_Lc)); T_out_HX_bank = T_fin
    ! *** REACTIVITY RELATED
    ALLOCATE(sol_PKE      (max_dg + 1))
    ALLOCATE(sol_PKE_prev (max_dg + 1))
    ALLOCATE(sol_PKE_temp (max_dg + 1))
    ALLOCATE(sol_PKE_delay(max_dg + 1))
    ALLOCATE(mat_sol_PKE  (max_dg + 1,bsize_L))
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
    DO i = 1,bsize_L
        mat_sol_PKE(:,i) = sol_PKE
    END DO
    ! *** MISCELLENIOUS
    P_RX      = P_RX0
    P_RX_prev = P_RX0
    ! *** INTERMEDIATE LOOP RELATED
    M_HX      = v_he_fuel*in_m*rho_f
    m_dot_HX  = W_f
    T_HX_in   = T_f2
    T_HX_out  = T_fin
    T_pr      = (T_f2 + T_fin)/2.0d0
    T_pr_prev = T_pr
    T_HX_salt = (546.0 + 579.0)/2.0d0
    GAMMA_HX  = P_RX0/(T_pr - T_HX_salt)
    ! *** UPWIND SCHEME EMPLOYED
    node_M_HX     = M_HX    /n_node_HX
    node_gamma_HX = GAMMA_HX/n_node_HX
    ALLOCATE(node_T_HX     (0:n_node_HX));      
    ALLOCATE(node_T_HX_prev(0:n_node_HX));
    ALLOCATE(node_T_pr     (1:n_node_HX));
    ALLOCATE(node_T_pr_prev(1:n_node_HX));
    node_T_HX(0) = T_f2
    !> ITERATE UNTIL node_T_HX(end) = T_HX_out (=T_fin) / FOR UPWIND SCHEME
    IF(mode_IHX == 3) THEN
        j = 0
        adjust_gamma: DO
            j = j + 1
            IF(j > 1) THEN
                tmp_x1 = tmp_x2
                tmp_y1 = tmp_y2
            END IF
            ! --- CALCULATE NODE-WISE Temperature within IHX
            DO i = 1,n_node_HX
                node_T_HX(i) = (m_dot_HX*cp_f*node_T_HX(i-1) + node_gamma_HX*T_HX_SALT) / (m_dot_HX*cp_f + node_gamma_HX)
            END DO
            ! --- USE SECANT METHOD
            tmp_x2 = node_gamma_HX
            tmp_y2 = node_T_HX(n_node_HX) - T_HX_out
            IF(j == 1) THEN
                node_gamma_HX = node_gamma_HX * 1.01
            ELSE
                node_gamma_HX = (-tmp_y1) * (tmp_x2 - tmp_x1) / (tmp_y2 - tmp_y1) + tmp_x1
            END IF
            IF(ABS(tmp_y2) < 1E-6) EXIT adjust_gamma
        END DO adjust_gamma
        node_T_HX_prev = node_T_HX
    !> INITIALIZATION FOR 1-D SCHEME (LUMPED)
    ELSE IF(mode_IHX == 4) THEN
        j = 0
        adjust_gamma2: DO
            j = j + 1
            IF(j > 1) THEN
                tmp_x1 = tmp_x2
                tmp_y1 = tmp_y2
            END IF
            ! --- CALCULATE NODE-WISE Temperature within IHX
            DO i = 1,n_node_HX
                node_T_pr(i) = (2.0d0*m_dot_HX*cp_f*node_T_HX(i-1) + node_gamma_HX*T_HX_SALT) / (2.0d0*m_dot_HX*cp_f + node_gamma_HX)
                node_T_HX(i) = 2.0d0*node_T_pr(i) - node_T_HX(i-1)
            END DO
            ! --- USE SECANT METHOD
            tmp_x2 = node_gamma_HX
            tmp_y2 = node_T_HX(n_node_HX) - T_HX_out
            IF(j == 1) THEN
                node_gamma_HX = node_gamma_HX * 1.01
            ELSE
                node_gamma_HX = (-tmp_y1) * (tmp_x2 - tmp_x1) / (tmp_y2 - tmp_y1) + tmp_x1
            END IF
            IF(ABS(tmp_y2) < 1E-6) EXIT adjust_gamma2
        END DO adjust_gamma2
        node_T_HX_prev = node_T_HX
        node_T_pr_prev = node_T_pr
    END IF
END SUBROUTINE initialize