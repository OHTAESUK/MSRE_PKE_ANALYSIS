MODULE variables
    IMPLICIT NONE
    ! *** GIVEN PARAMETERS & VARIABLES FOR MSRE 8MW U233 BENCHMARK (1GROUP)    
    REAL(8) :: del_t  = 0.020   ! [s] Time-step for calculation 
    REAL(8) :: time   = 100.0   ! [s] Total simulation time
    INTEGER :: n_step           ! [-] Number of step for transient simulation
    
    REAL(8) :: tau_L    = 16.73   ! [s] Loop 
    REAL(8) :: tau_c    =  8.46   ! [s] Core
    REAL(8) :: tau_Lc   =  8.365  ! [s] Cold leg time-delay
    REAL(8) :: tau_Lh   =  8.365  ! [s] Hot  leg time-delay
    REAL(8) :: tau_r_HX =  8.240  ! [s] Time-delay for Radiator (SHX) -> IHX
    REAL(8) :: tau_HX_r =  4.710  ! [s] Time-delay for IHX -> Radiator (SHX)
    
    INTEGER :: bsize_L    ! [-] Number of bank-size for loop (tau_L)
    INTEGER :: bsize_c    ! [-] Number of bank-size for core (tau_c)
    INTEGER :: bsize_Lc   ! [-] Number of bank-size for cold leg (tau_Lc)
    INTEGER :: bsize_Lh   ! [-] Number of bank-size for hot  leg (tau_Lh)
    INTEGER :: bsize_r_HX ! [-] Number of bank-size for SHX -> IHX time leg (tau_r_HX) 
    INTEGER :: bsize_HX_r ! [-] Number of bank-size for IHX -> SHX time leg (tau_HX_r) 
    
    INTEGER :: max_dg = 6       ! [-] Number of delayed neutron groups
    REAL(8) :: LAMBDA = 4.0E-04 ! [s] Generation time 
    
    REAL(8) :: vdot_f = 7.57080E-02           ! [m^3/s]   Flow Rate
    REAL(8) ::  rho_f = 2.14647E+03           ! [kg/m^3]  Density of fuel salt
    REAL(8) ::   cp_f = 1.96650E-03           ! [MJ/kg-K] Specific Heat of the fuel salt
    REAL(8) ::    W_f = 1.623879934566580e+02 ! [kg/s]    Fuel flow rate
    
    REAL(8) ::   v_g = 1.95386  ! [m^3]     Graphite Volume
    REAL(8) :: rho_g = 1.860E3  ! [kg/m^3]  Graphite Density
    REAL(8) ::  cp_g = 1.773E-3 ! [MJ/kg-K] Specific Heat of Graphite
    
    REAL(8) :: hA_fg = 0.02*9/5 ! [MW/K] Fuel to Graphite heat transfer coefficient
    
    REAL(8) :: k_f  = 0.93 ! [-] fraction of heat generated in fuel - that generated in external loop
    REAL(8) :: k_g  = 0.07 ! [-] fraction of total power generated in the graphite
    REAL(8) :: k_1  = 0.5  ! [-] fraction of heat transferred from graphite which goes to first  fuel lump
    REAL(8) :: k_2  = 0.5  ! [-] fraction of heat transferred from graphite which goes to second fuel lump
    
    REAL(8) :: a_f = -11.034E-5 ! [1/K] dRho/dTf
    REAL(8) :: a_g =  -5.814E-5 ! [1/K] dRho/dTg
    
    REAL(8) :: m_f   ! [kg] Fuel mass in core
    REAL(8) :: m_f1  ! [kg] Fuel mass in fuel node 1
    REAL(8) :: m_f2  ! [kg] Fuel mass in fuel node 2
    REAL(8) :: m_g   ! [kg] Graphite mass

    REAL(8) :: k_f1  ! [-] fraction of total power generated in lump f1
    REAL(8) :: k_f2  ! [-] fraction of total power generated in lump f2
    
    REAL(8) :: beta  ! [-] sum(d_beta)
    REAL(8) :: lbda  ! [-] sum(d_lbda)
    
    REAL(8), ALLOCATABLE :: d_lbda(:) ! [-] decay constant
    REAL(8), ALLOCATABLE :: d_beta(:) ! [-] delayed neutron fraction
    
    REAL(8) :: Rho_ext = 14.24E-5 ! [-] External Reactivity
    
    ! *** VALUES THAT CHANGE DURING CALCULATION
    REAL(8) :: T_fin = 6.3222E+02 ! [C]
    REAL(8) :: T_f2  = 6.5727E+02 ! [C]
    REAL(8) :: T_f1               ! [C]
    REAL(8) :: T_g                ! [C]
    REAL(8) :: T_f                ! [C] = 0.5 * (T_f1 + T_f2)
    REAL(8) :: T_fin0             ! [C] Initial Temperature
    REAL(8) :: T_g0               ! [C] Initial Temperature
    REAL(8) :: T_f0               ! [C] Initial Temperature
    REAL(8) :: T_f10              ! [C] Initial Temperature
    REAL(8) :: T_f20              ! [C] Initial Temperature
    REAL(8) :: dT_f               ! [C] Temperature Difference (fixed)
    
    REAL(8) :: T_f1_prev, T_f1_temp ! [C]
    REAL(8) :: T_f2_prev, T_f2_temp ! [C]
    REAL(8) :: T_g_prev,  T_g_temp  ! [C]
    
    REAL(8) :: P_RX           ! [MW] Thermal power
    REAL(8) :: P_RX_prev      ! [MW] Thermal power
    REAL(8) :: P_RX0 = 8.0d0  ! [MW] Thermal Power (initial)
    REAL(8) :: Rho_dyn        ! [-]  Reactivity change in going from stationary to circulating fuel
    REAL(8) :: Rho_dyn_prev   ! [-]  Reactivity change in going from stationary to circulating fuel
    REAL(8) :: Rho_dyn0       ! [-]  Reactivity change in going from stationary to circulating fuel (initial)
    
    REAL(8), ALLOCATABLE :: sol_PKE      (:)  ! [-] Solution to PKE at k   step --- micro interval
    REAL(8), ALLOCATABLE :: sol_PKE_prev (:)  ! [-] Solution to PKE at k-1 step --- micro interval
    REAL(8), ALLOCATABLE :: sol_PKE_temp (:)  ! [-] Solution to PKE / Temporary Storage
    REAL(8), ALLOCATABLE :: sol_PKE_delay(:)  ! [-] Solution to PKE / delayed tau_L
    REAL(8), ALLOCATABLE :: mat_sol_PKE(:,:)  ! [-] Stores time-step wise PKE solution -> used_as sol_PKE(t - tau_L); bsize_L
    
    ! VALUES RELATED TO INTERMIDATE LOOP
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! --- [RETRIEVED DATA FROM SIMULINK CODE; Annals of Nuclear Energy 113 (2018) 177-193]
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    REAL(8), PARAMETER :: pi   = 4.D0*DATAN(1.D0)      ! PI-VALUE
    REAL(8), PARAMETER :: in_m = 1.63871e-5            ! 1 cubic inch = 1.63871e-5 cubic meters
    REAL(8), PARAMETER :: cp_HX = 2.39E-3              ! [MJ/kg-K] Specific Heat of IHX coolant
    REAL(8), PARAMETER :: rho_s = 1.922E3              ! [kg/m^3]  coolant salt density (in IHX)
    
    REAL(8), PARAMETER :: d_he = 16.d0;                             ! (in) he diameter ORNL-TM-0728 p. 164
    REAL(8), PARAMETER :: h_he = 72.d0;                             ! (in) active height % 96; %(in) he height ORNL-TM-0728 p. 164
    REAL(8), PARAMETER :: od_tube = 0.5d0;                          ! (in) coolant tube OD ORNL-TM-0728 p. 164
    REAL(8), PARAMETER :: id_tube = od_tube - 2*0.042;              ! (in) coolant tube ID ORNL-TM-0728 p. 164
    INTEGER, PARAMETER :: n_tube = 159;                             ! (-) number of coolant tubes ORNL-TM-0728 p. 164
    REAL(8), PARAMETER :: a_tube = 254.d0*144.d0;                   ! (in^2) total area of tubes ORNL-TM-0728 p. 164
    REAL(8), PARAMETER :: l_tube = a_tube/n_tube/(pi*od_tube);      ! (in) tube length
    REAL(8), PARAMETER :: v_tube = n_tube*pi*(od_tube/2)**2*l_tube; ! (in^3) hx shell volume occupied by tubes
    REAL(8), PARAMETER :: v_cool = n_tube*pi*(id_tube/2)**2*l_tube; ! (in^3) hx volume occupied by coolant
    REAL(8), PARAMETER :: v_he = (d_he/2.0d0)**2*pi*h_he;           ! (in^3) volume of heat exchanger shell
    REAL(8), PARAMETER :: v_he_fuel = v_he-v_tube;                  ! (in^3) volume available to fuel in shell
    
    ! --- FOR RADIATOR (NEED MASS & mass-flow-rate for coolant salt)
    INTEGER, PARAMETER :: n_rtubes = 120              ! [-] number of tubes in radiator (rows times tubes per row) ORNL-TM-0728 p.296
    REAL(8), PARAMETER :: l_rtube  = 9.144            ! [m] length of tubes in radiator ORNL-TM-0728 p.296
    REAL(8), PARAMETER :: od_rad   = 0.01905          ! [m] outer diameter of tubes in radiator ORNL-TM-0728 p.296
    REAL(8), PARAMETER :: tube_wall_thick = 0.0018288 ! [m] thickness of tubes in radiator ORNL-TM-0728 p.296
    REAL(8), PARAMETER :: id_rad = od_rad-2*tube_wall_thick;                 ! [m]   used for the following calculation
    REAL(8), PARAMETER :: v_rp   = pi*(id_rad/2)**2*l_rtube*REAL(n_rtubes,8) ! [m^3] volume available to salt in radiator
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    REAL(8) :: m_dot_IHX            ! [kg/s] Mass flow rate to the IHX
    REAL(8) :: m_dot_SHX            ! [kg/s] Mass flow rate to the SHX
    REAL(8) :: M_IHX                ! [kg] Primary salt residing in IHX  
    REAL(8) :: M_IHX_SALT           ! [kg] 2ndary  salt residing in IHX (coolant in IHX)
    REAL(8) :: M_SHX                ! [kg] 2ndary  salt residing in SHX (rad; radiator)
    
    REAL(8) :: T_IHX_in            ! [C]  Inlet   temperature to the IHX
    REAL(8) :: T_IHX_out           ! [C]  Outlet  temperature to the IHX
    REAL(8) :: T_IHX_salt_in       ! [C]  Inlet   temperature to the IHX (connected to SHX); coolant salt
    REAL(8) :: T_IHX_salt_out      ! [C]  Outlet  temperature to the IHX (connected to SHX); coolant salt
    REAL(8) :: T_IHX_salt_avg      ! [C]  Average temperature of the IHX (connected to SHX); coolant salt
    REAL(8) :: T_IHX_salt_avg_prev ! [C]  Average temperature of the IHX (connected to SHX); coolant salt / (delta_t before)
    REAL(8) :: T_SHX_in            ! [C]  Inlet   temperature to the SHX 
    REAL(8) :: T_SHX_out           ! [C]  Outlet  temperature to the SHX
    
    REAL(8) :: T_pr                 ! [C]  Primary temperature in the IHX 
    REAL(8) :: T_pr_prev            ! [C]  Primary temperature in the IHX (PREVIOUS)  / T_pr(t_s - delta_t)
    REAL(8) :: T_pr_LAG             ! [C]  Primary temperature in the IHX (LAG)       / Finite time-lag considered
    
    REAL(8) :: T_sr                 ! [C]  Secondary temperature in the IHX (=average of T_IHX_salt_in & T_IHX_salt_out)
    REAL(8) :: T_sr_prev            ! [C]  Secondary temperature in the IHX (=average of T_IHX_salt_in & T_IHX_salt_out)
    
    REAL(8) :: T_SINK               ! [C]  Final temperature (fixed)
    
    REAL(8) :: GAMMA_IHX            ! [kg/s]  Flow rate of IHX (2ndary   loop)
    REAL(8) :: GAMMA_SHX            ! [kg/s]  Flow rate of SHX (tertiary loop)
    
    REAL(8), ALLOCATABLE :: T_out_RX_bank (:) ! BANKING of outlet temperature from the CORE T_IHX_in(t) = T_out_RX_BANK    (t - tau_Lh) 
    REAL(8), ALLOCATABLE :: T_out_IHX_bank(:) ! BANKING of outlet temperature from the HX   T_fin(t)    = T_out_IHX_bank   (t - tau_Lc)
    
    REAL(8), ALLOCATABLE :: T_out_IHX_SALT_bank(:) ! BANKING of outlet temperature from the IHX (SALT) T_SHX_in(t) = T_out_IHX_SALT_BANK(t - tau_HX_r)
    REAL(8), ALLOCATABLE :: T_out_SHX_SALT_bank(:) ! BANKING of outlet temperature from the SHX (SALT) T_IHX_in(t) = T_out_SHX_SALT_BANK(t - tau_r_HX)
    
    ! (muted) ! *** VARIABLES RELATED TO THE UPWIND SCHEME (1-D IHX)
    ! (muted) REAL(8)              :: node_M_HX
    ! (muted) REAL(8)              :: node_gamma_HX
    ! (muted) REAL(8), ALLOCATABLE :: node_T_HX     (:)
    ! (muted) REAL(8), ALLOCATABLE :: node_T_HX_prev(:)
    ! (muted) REAL(8), ALLOCATABLE :: node_T_pr     (:)
    ! (muted) REAL(8), ALLOCATABLE :: node_T_pr_prev(:)
    
END MODULE variables