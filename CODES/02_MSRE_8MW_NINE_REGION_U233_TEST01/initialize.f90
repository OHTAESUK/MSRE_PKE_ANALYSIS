SUBROUTINE initialize
    USE variables
    USE matrix_solver
    IMPLICIT NONE
    INTEGER :: i,j
    REAL(8) :: A(3,3)
    REAL(8) :: b(3)
    REAL(8) :: x(3)    
    REAL(8) :: tmp_val1
    REAL(8) :: tmp_x1,tmp_y1
    REAL(8) :: tmp_x2,tmp_y2
    
    ! *** SIMULATION TIME RELATED
    n_step     = FLOOR(time    /del_t + 1.0E-9)
    bsize_L    = FLOOR(tau_L   /del_t + 1.0E-9)
    bsize_Lc   = FLOOR(tau_Lc  /del_t + 1.0E-9)
    bsize_Lh   = FLOOR(tau_Lh  /del_t + 1.0E-9)
    bsize_r_HX = FLOOR(tau_r_HX/del_t + 1.0E-9) 
    bsize_HX_r = FLOOR(tau_HX_r/del_t + 1.0E-9)
    
    ! *** DELAYED NEUTRON RELATED
    ALLOCATE(d_lbda(max_dg)); d_lbda = [1.260E-02, 3.370E-02, 1.390E-01, 3.250E-01, 1.130E+00, 2.500E+00]
    ALLOCATE(d_beta(max_dg)); d_beta = [0.0002300, 0.0007900, 0.0006700, 0.0007300, 0.0001300, 0.0000900]
    lbda = SUM(d_lbda)
    beta = SUM(d_beta)
    
    ! *** FLOW RELATED
    tau_fa = [1.386d0, 2.083d0, 1.139d0, 1.424d0, 2.084d0, 1.139d0, 1.424d0, 2.371d0, 1.610d0]
    tau_fb = [1.454d0, 1.424d0, 1.139d0, 2.772d0, 1.424d0, 1.139d0, 2.774d0, 1.380d0, 2.700d0]
    
    ka = [0.01493d0, 0.02736d0, 0.04504d0, 0.05126d0, 0.03601d0, 0.06014d0, 0.06845d0, 0.06179d0, 0.09333d0]
    kb = [0.01721d0, 0.04550d0, 0.04656d0, 0.04261d0, 0.06069d0, 0.06218d0, 0.05664d0, 0.07707d0, 0.07311d0]
    kga = [0.000946d0, 0.001685d0, 0.003029d0, 0.003447d0, 0.002216d0, 0.004044d0, 0.004603d0, 0.003920d0, 0.006277d0]
    kgb = [0.001081d0, 0.003060d0, 0.003131d0, 0.002395d0, 0.004081d0, 0.004182d0, 0.003184d0, 0.005183d0, 0.004305d0]
    
    hA = [0.000392, 0.001204, 0.000900, 0.001174, 0.001977, 0.001525, 0.001985, 0.005445, 0.005360] * (9.0d0/5.0d0)
    
    mcp_fa = [0.0151, 0.0512, 0.0280, 0.0350, 0.0866, 0.0473, 0.0592, 0.2380, 0.1615] * (9.0d0/5.0d0)
    mcp_fb = [0.0158, 0.0349, 0.0280, 0.0682, 0.0592, 0.0473, 0.1152, 0.1384, 0.2710] * (9.0d0/5.0d0)
    mcp_g  = [0.0700, 0.2114, 0.1606, 0.2056, 0.3576, 0.2718, 0.3478, 0.9612, 0.9421] * (9.0d0/5.0d0)
    
    W_f1 = (mcp_fa(1)/cp_f)/tau_fa(1)
    W_f2 = (mcp_fa(2)/cp_f)/tau_fa(2)
    W_f3 = (mcp_fa(5)/cp_f)/tau_fa(5)
    W_f4 = (mcp_fa(8)/cp_f)/tau_fa(8)
    m_fm = W_f * tau_m
    
    ! +++ NORMALIZE SO THAT SUM OF {ka,kb,kga,kgb} becomes unity
    tmp_val1 = SUM(ka) + SUM(kb) + SUM(kga) + SUM(kgb)
    ka = ka/tmp_val1
    kb = kb/tmp_val1
    kga = kga/tmp_val1
    kgb = kgb/tmp_val1    
    
    ! +++ REACTIVITY FEEDBACK WEIGHTING
    Ifa = [0.02168, 0.02197, 0.07897, 0.08249, 0.02254, 0.08255, 0.08623, 0.02745, 0.06936]
    Ifb = [0.02678, 0.06519, 0.08438, 0.04124, 0.06801, 0.08823, 0.04290, 0.05529, 0.03473]
    Ig  = [0.04443, 0.08835, 0.16671, 0.12077, 0.09181, 0.17429, 0.12612, 0.08408, 0.10343]
    
    ! *** TEMPERATURE RELATED (SOLVE 3x3 matrix equation) for each region
    temp_initialize: DO i = 1,9
        A(:,:) = 0.d0
        b(:)   = 0.d0       
        !> TEMPERATURE 'a'
        A(1,1) = -1.0d0/tau_fa(i) - (kga(i)/(kga(i)+kgb(i))*hA(i)/(mcp_fa(i)))
        A(1,2) = +0.d0
        A(1,3) =                  + (kga(i)/(kga(i)+kgb(i))*hA(i)/(mcp_fa(i)))
        !> TEMPERATURE 'b'
        A(2,1) = +1.0d0/tau_fb(i)
        A(2,2) = -1.0d0/tau_fb(i) - (kgb(i)/(kga(i)+kgb(i))*hA(i)/(mcp_fb(i)))
        A(2,3) =                  + (kgb(i)/(kga(i)+kgb(i))*hA(i)/(mcp_fb(i)))
        !> TEMPERATURE 'g'
        A(3,1) = (kga(i)/(kga(i)+kgb(i)))*hA(i)/mcp_g(i)
        A(3,2) = (kgb(i)/(kga(i)+kgb(i)))*hA(i)/mcp_g(i)
        A(3,3) =                         -hA(i)/mcp_g(i)
        !> RHS        
        IF(ANY([1,2,5,8] == i) .EQ. .TRUE.) THEN
            tmp_val1 = T_fin
        ELSE
            tmp_val1 = T_fb(i-1)
        END IF
        b(1) = -ka(i)*P_RX0/mcp_fa(i) - 1.0d0/tau_fa(i)*tmp_val1
        b(2) = -kb(i)*P_RX0/mcp_fb(i)
        b(3) = -(kga(i)+kgb(i))*P_RX0/mcp_g(i)
        CALL solve_GE_piv(A,b,x)
        T_fa(i) = x(1)
        T_fb(i) = x(2)
        T_g(i)  = x(3)      
        T_f(i)  = 0.5d0 * (x(1) + x(2))  
    END DO temp_initialize
    ! >>> SUMMATION THROUGH MIXING REGION (2sec)    
    T_fm = ((W_f1/m_fm)*T_fb(1) + (W_f2/m_fm)*T_fb(4) + (W_f3/m_fm)*T_fb(7) + (W_f4/m_fm)*T_fb(9))/((W_f1+W_f2+W_f3+W_f4)/m_fm)
    ! >>> STORE THE INITIAL INFORMATION
    T_fin0 = T_fin
    T_fm0  = T_fm
    T_fa0  = T_fa
    T_fb0  = T_fb
    T_f0   = T_f
    T_g0   = T_g
    ! >>> STORE FOR PREVIOUS STEP (ITERATION RELATED)
    T_fa_prev = T_fa
    T_fb_prev = T_fb
    T_g_prev  = T_g
    ! >>> STORES T_f2 (for T_fin)    
    ALLOCATE(T_out_RX_bank (bsize_Lh)); T_out_RX_bank  = T_fm
    ALLOCATE(T_out_IHX_bank(bsize_Lc)); T_out_IHX_bank = T_fin
    
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
    M_IHX      = v_he_fuel*in_m*rho_f
    m_dot_IHX  = W_f
    T_IHX_in   = T_fm
    T_IHX_out  = T_fin
    T_pr      = (T_fm + T_fin)/2.0d0
    T_pr_prev = T_pr
    
    M_IHX_salt     = v_cool*in_m*rho_s 
    T_IHX_salt_in  = 546.0d0
    T_IHX_salt_out = 579.0d0
    T_IHX_salt_avg = (T_IHX_salt_in + T_IHX_salt_out)/2.0d0
    T_sr      = (T_IHX_salt_in + T_IHX_salt_out)/2.0d0
    T_sr_prev = T_sr
    GAMMA_IHX = P_RX0/(T_pr - T_IHX_salt_avg)
    
    M_SHX     = v_rp * rho_s                                 ! [kg]   coolant salt mass in rad (kg)
    m_dot_SHX = P_RX0/cp_HX/(T_IHX_salt_out - T_IHX_salt_in) ! [kg/s] coolant flow rate (kg/s)

    T_SHX_in  = T_IHX_salt_out
    T_SHX_out = T_IHX_salt_in
    T_SINK    = (37.78 + 148.9)/2.0d0 ! [C] air inlet & outlet temperature ORNL-TM-1647 p.2
    GAMMA_SHX = P_RX0/(T_sr - T_SINK)
    ! >>> STORES PREVIOUS TEMPERATURE
    ALLOCATE(T_out_IHX_SALT_bank(bsize_HX_r)); T_out_IHX_SALT_bank = T_IHX_salt_out
    ALLOCATE(T_out_SHX_SALT_bank(bsize_r_HX)); T_out_SHX_SALT_bank = T_IHX_salt_in
END SUBROUTINE initialize