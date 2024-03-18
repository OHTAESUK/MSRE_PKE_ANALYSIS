MODULE variables
    IMPLICIT NONE
    ! *** GIVEN PARAMETERS FOR MSRE 8MW U233 BENCHMARK (1GROUP)
    REAL(8) :: del_t  = 0.02    ! [s] Time-step for calculation 
    REAL(8) :: time   = 100.0   ! [s] Total simulation time
    INTEGER :: n_step
    
    REAL(8) :: tau_L  = 16.73   ! [s] Loop 
    REAL(8) :: tau_c  =  8.46   ! [s] Core
    REAL(8) :: LAMBDA = 4.0E-04 ! [s] Generation time 
    
    INTEGER :: max_dg = 6
    
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
    
    ! *** VALUES TO BE INITIALIZED
    REAL(8) :: m_f   ! [kg] Fuel mass in core
    REAL(8) :: m_f1  ! [kg] Fuel mass in fuel node 1
    REAL(8) :: m_f2  ! [kg] Fuel mass in fuel node 2
    REAL(8) :: m_g   ! [kg] Graphite mass

    REAL(8) :: k_f1  ! [-] fraction of total power generated in lump f1
    REAL(8) :: k_f2  ! [-] fraction of total power generated in lump f2
    
    REAL(8) :: beta
    REAL(8) :: lbda    
    REAL(8), ALLOCATABLE :: d_lbda(:)
    REAL(8), ALLOCATABLE :: d_beta(:)
    
    ! *** VALUES RELATED TO EXTERNAL REACTIVITY
    REAL(8) :: Rho_ext = 0.d0
    
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
    
    REAL(8) :: T_f1_prev, T_f1_temp
    REAL(8) :: T_f2_prev, T_f2_temp
    REAL(8) :: T_g_prev,  T_g_temp
    
    REAL(8) :: P_RX           ! [MW] Thermal power
    REAL(8) :: P_RX_prev      ! [MW] Thermal power
    REAL(8) :: P_RX0 = 8.0d0  ! [MW] Thermal Power (initial)
    REAL(8) :: Rho_dyn        ! [-]  Reactivity change in going from stationary to circulating fuel
    REAL(8) :: Rho_dyn_prev   ! [-]  Reactivity change in going from stationary to circulating fuel
    REAL(8) :: Rho_dyn0       ! [-]  Reactivity change in going from stationary to circulating fuel (initial)
    
    REAL(8), DIMENSION(:), ALLOCATABLE :: sol_PKE        ! Solution to PKE at k   step --- micro interval
    REAL(8), DIMENSION(:), ALLOCATABLE :: sol_PKE_prev   ! Solution to PKE at k-1 step --- micro interval
    REAL(8), DIMENSION(:), ALLOCATABLE :: sol_PKE_temp   ! Solution to PKE / Temporary Storage
    REAL(8), DIMENSION(:), ALLOCATABLE :: sol_PKE_delay  ! Solution to PKE / delayed tau_L

    ! STORING VARIABLES TO CONSIDER TIME-LAG (tau_L)
    INTEGER              :: num_delay        ! [-] Number of elements to be stored for delayed response (tau_L)
    REAL(8), ALLOCATABLE :: arr_T_f2   (:)   ! [C] Stores time-step wise T_f2 value   -> used as T_fin(t) = T_f2(t - tau_L)
    REAL(8), ALLOCATABLE :: mat_sol_PKE(:,:) ! [-] Stores time-step wise PKE solution -> used_as         sol_PKE(t - tau_L)
    
END MODULE variables