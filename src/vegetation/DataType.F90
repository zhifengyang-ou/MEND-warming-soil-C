module DataType
  ! use
  
  implicit none
  
  !public:: 
  
  !----------------- constants ---------------------------------------
  !===============constants===============
  logical, public, parameter :: read_from_parameter_file = .TRUE.
  integer, public, parameter :: days_per_year  = 365
  integer, public, parameter :: hours_per_year = 365 * 24  ! 8760
  real,    public, parameter :: seconds_per_year = 365. * 24. * 3600.
  real,    public, parameter :: seconds_per_day = 24. * 3600.  
  ! Physical constants
  real, public, parameter :: TFREEZE = 273.16
  real, public, parameter :: Rugas = 8.314472 ! universal gas constant, J K-1 mol-1
  real, public, parameter :: mol_C = 12.0e-3 ! molar mass of carbon, kg
  real, public, parameter :: mol_air = 28.96440e-3 ! molar mass of air, kg
  real, public, parameter :: mol_CO2 = 44.00995e-3 ! molar mass of CO2,kg
  real, public, parameter :: mol_h2o = 18.0e-3 ! molar mass of water, kg
  real, public, parameter :: cpair = 1010.
  real, public, parameter :: H2OLv0=2.501e6   !latent heat H2O (J/kg)
  real, public, parameter :: p_sea = 101325.  ! atmospheric pressure  (Pa)
  real, public, parameter :: DENS_H2O = 1000. ! kg m-3
  real, public, parameter :: PI = 3.1415926   
 
  ! Vegetation and soil types
  integer, public, parameter :: MSPECIES = 15
  integer, public, parameter :: LEAF_ON = 1
  ! -------------------------------------------------------------------

  !===============data types ==============================
  
  !Forcing data type
  type :: climate_data_type
     integer :: year          ! Year
     integer :: doy           ! day of the year
     real :: hod           ! hour of the day
     real :: PAR           ! umol m-2 s-1
     real :: radiation     ! W/m2
     real :: Tair          ! air temperature,  K
     real :: Tsoil         ! soil temperature, K
     real :: RH            ! relative humidity
     real :: rain          ! kgH2O m-2 s-1
     real :: windU         ! wind velocity (m s-1)
     real :: P_air         ! pa
     real :: CO2           ! ppm
     real :: soilwater     ! soil moisture, vol/vol
  end type climate_data_type
  
  !-----------PFT data type----------------
  type spec_data_type
    integer :: lifeform  ! 0 for grasses, 1 for trees
    integer :: phenotype ! phenology type: 0 for deciduous, 1 for evergreen
    integer :: pt        ! photosynthetic physiology of species
    ! leaf traits
    real :: LMA          ! leaf mass per unit area, kg C/m2
    real :: leafLS       ! leaf life span
    real :: alpha_L      ! leaf turn over rate
    real :: LNA          ! leaf Nitrogen per unit area, kg N/m2
    real :: LNbase       ! basal leaf Nitrogen per unit area, kg N/m2, (Rubisco)
    real :: CN0leafST    ! leaf structural tissues, 175
    real :: leaf_size    ! characteristic leaf size
    real :: leafTK       ! leaf thickness, m
    real :: rho_leaf     ! leaf mass density (kgC/m3)
    real :: alpha_ps     ! photosynthesis efficiency
    real :: m_cond       ! factor of stomatal conductance
    real :: Vmax         ! max rubisco rate, mol m-2 s-1
    real :: Vannual      ! annual productivity per unit area at full fun (kgC m-2 yr-1)
    real :: gamma_L      ! leaf respiration coeficient (per yr)
    real :: gamma_LN     ! leaf respiration coeficient per unit N
    real :: ps_wet       ! wet leaf photosynthesis down-regulation
    real :: LFR_rate     ! Leaf filling rate per day
    ! root traits
    real :: rho_FR       ! material density of fine roots (kgC m-3)
    real :: root_r       ! radius of the fine roots, m
    real :: root_zeta    ! e-folding parameter of root vertical distribution (m)
    real :: root_frac(soil_L)    ! root fraction
    real :: SRA          ! speific fine root area, m2/kg C
    real :: SRL          ! specific root lenght
    real :: gamma_FR     ! Fine root respiration rate, kgC kgC-1 yr-1
    real :: alpha_FR     ! Turnover rate of Fine roots, fraction yr-1
    real :: Kw_root      ! fine root water conductivity mol m m-2 s−1 Pa−1 !
    real :: root_perm
    !  real :: rho_N_up0   ! maximum N uptake rate
    !  real :: N_roots0    ! root biomass at half of max. N-uptake rate
    real :: NfixRate0    ! Reference N fixation rate (kgN kgC-1 root)
    real :: NfixCost0    ! Carbon cost of N fixation (kgC kgN-1)
    ! wood traits
    real :: rho_wood     ! woody density, kg C m-3 wood
    real :: gamma_SW     ! sapwood respiration rate, kgC m-2 Acambium yr-1
    real :: f_taper

    ! Plant hydraulics
    real :: kx0  ! xylem conductivity, (mm/s)/(Mpa/m)
    real :: WTC0 ! xylem water transfer capacity, m/lifetime
    real :: CR_Leaf ! leaf compression ratio per MPa
    real :: CR_Wood ! Wood compression ratio per MPa
    real :: w0L_max ! leaf maximum water/carbon ratio
    real :: w0S_max ! stem maximum water/carbon ratio
    real :: w0L_min ! leaf minimum water/carbon ratio
    real :: w0S_min ! stem minimum water/carbon ratio
    real :: psi0_LF ! minimum leaf water potential
    real :: psi0_WD ! minimum stem wood potential
    real :: psi50_WD !wood potential at which 50% conductivity lost, MPa
    real :: Kexp_WD  ! exponent of the PLC curve
    real :: f_supply ! fraction of stem water available for leaves per hour

    ! Allometry
    real :: alphaHT, thetaHT ! height = alphaHT * DBH ** thetaHT
    real :: alphaCA, thetaCA ! crown area = alphaCA * DBH ** thetaCA
    real :: alphaBM, thetaBM ! biomass = alphaBM * DBH ** thetaBM
    real :: phiRL            ! ratio of fine root to leaf area
    real :: phiCSA           ! ratio of sapwood CSA to target leaf area
    real :: tauNSC           ! residence time of C in NSC (to define storage capacity)
    real :: fNSNmax          ! multiplier for NSNmax
    real :: f_N_add
    real :: transT           ! Structural transitional time for canopy layer trees
    ! Default C/N ratios
    real :: CNleaf0
    real :: CNroot0
    real :: CNsw0
    real :: CNwood0
    real :: CNseed0
    ! phenology
    real :: tc_crit_off     ! K, for turning OFF a growth season
    real :: tc_crit_on      ! K, for turning ON a growth season
    real :: gdd_crit        ! K, critical value of GDD5 for turning ON growth season
    real :: gdd_par1
    real :: gdd_par2
    real :: gdd_par3
    real :: betaON         ! Critical soil moisture for phenology ON
    real :: betaOFF        ! Critical soil moisture for phenology OFF
    !  vital rates
    real :: AgeRepro       ! the age that can reproduce
    real :: v_seed           ! fracton of G_SF to G_F
    real :: s0_plant     ! size of the seedlings, kgC/indiv
    real :: prob_g,prob_e    ! germination and establishment probabilities
    real :: r0mort_c     ! yearly mortality rate in canopy
    real :: D0mu         ! Reference diameter for size-dependent mortality
    real :: A_un         ! Parameter for understory mortality affected by layers
    real :: A_sd         ! Max multiplier for seedling mortality
    real :: B_sd         ! Mortality sensitivity for seedlings
    real :: A_D          ! Sensitivity to dbh
    real :: s_hu         ! hydraulic mortality sensitivity
    ! Population level variables
    real :: LAImax    ! max. LAI
    real :: LAImax_u  ! max. LAI understorey
    real :: LAI_light ! light controlled maximum LAI
    integer :: n_cc   ! for calculating LAImax via cc%LAImax derived from cc%NSN
    real :: f_cGap    ! fraction of internal gaps in the canopy
    ! "internal" gaps are the gaps that are created within the canopy by the
    ! branch fall processes.
  end type

  !----------cohort-----------------
  type :: cohort_type
    ! for climate-vegetation type
    integer :: phenotype    ! phenology type: 0 for deciduous, 1 for evergreen
    logical :: firstday     ! First day of a growing season
    integer :: pt           ! photosynthetic physiology of species

    ! ---- biological prognostic variables
    integer :: ccID   = 0   ! cohort ID
    integer :: species= 0   ! vegetation species
    real :: gdd       = 0.0   ! for phenology
    real :: ALT       = 0.0  ! growing season accumulative cold temperature
    integer :: Ngd    = 0   ! growing days
    integer :: Ndm    = 0   ! dormant days
    integer :: Ncd    = 0   ! number of cold days in non-growing season
    integer :: status = 0   ! growth status of plant: 1 for ON, 0 for OFF
    integer :: layer  = 1   ! the layer of this cohort (numbered from top, top layer=1)
    real :: layerfrac = 0.0 ! fraction of layer area occupied by this cohort
    real :: leafage   = 0.0 ! leaf age (year)

    ! for populatin structure
    real :: nindivs= 1.0 ! density of vegetation, individuals/m2
    real :: mu     = 0.02 ! Cohort mortality rate
    real :: age    = 0.0 ! age of cohort, years
    real :: dbh    = 0.0 ! diameter at breast height, m
    real :: height = 0.0 ! vegetation height, m
    real :: Acrown = 1.0 ! crown area, m2/individual
    real :: Aleaf  = 0.0 ! total area of leaves, m2/individual
    real :: lai    = 0.0 ! crown leaf area index, m2/m2
    real :: D_bark = 0.0 ! thickness of bark
    ! carbon pools
    real :: bl     = 0.0 ! biomass of leaves, kg C/individual
    real :: br     = 0.0 ! biomass of fine roots, kg C/individual
    real :: bsw    = 0.0 ! biomass of sapwood, kg C/individual
    real :: bHW    = 0.0 ! biomass of heartwood, kg C/individual
    real :: seedC  = 0.0 ! biomass put aside for future progeny, kg C/individual
    real :: nsc    = 0.0 ! non-structural carbon, kg C/individual

    ! ----- carbon fluxes
    real :: gpp  = 0.0 ! gross primary productivity kg C/step
    real :: npp  = 0.0 ! net primary productivity kg C/step
    real :: resp = 0.0 ! plant respiration
    real :: resl = 0.0 ! leaf respiration
    real :: resr = 0.0 ! root respiration
    real :: resg = 0.0 ! growth respiration
    real :: NPPleaf,NPProot,NPPwood ! to record C allocated to leaf, root, and wood

    ! for hydraulics-mortality
    integer :: Nrings = 1
    real :: psi_s0   ! Equilibrium stem base water potential (soil-stem flux=0)
    real :: psi_leaf ! MPa, leaf water potential
    real :: psi_stem ! MPa, stem water potential
    real :: H_leaf ! Leaf capacitance, kgH2O MPa-1 (per tree)
    real :: H_stem ! Stem capacitance, kgH2O MPa-1 (per tree)
    real :: W_leaf ! Leaf water content, kgH2O (per tree)
    real :: W_stem ! Stem water content, kgH2O (per tree)
    real :: W_dead ! water storage in heartwood, just for balance counting.
    real :: Wmax_l ! Leaf max water content, kgH2O (per tree)
    real :: Wmax_s ! Stem max water content, kgH2O (per tree)
    real :: Wmin_l ! Leaf min water content, kgH2O (per tree)
    real :: Wmin_s ! Stem min water content, kgH2O (per tree)
    real :: V_stem ! Volumn of stems (including trunk)
    real :: V_leaf ! Volumn of leaves
    real :: Q_stem ! water flux from soil to stems (kg/tree/step)
    real :: Q_leaf ! water flux from stems to leaves (kg/tree/step)

    real :: Ktrunk ! trunk water conductance, m/(s MPa)
    real :: Asap ! Functional cross sectional area
    real :: Atrunk ! Sum of all rings
    real :: treeHU ! total water transported by the functional sapwood, m^3
    real :: treeW0 ! total WTC0 of the sapwood, m^3
    real :: Kx(Ysw_max) = 0.0 ! Initial conductivity of the woody generated in each year
    real :: WTC0(Ysw_max) = 0.0 ! lifetime water transfer capacity
    real :: accH(Ysw_max) = 0.0 ! m, total water transport for functional conduits
    real :: farea(Ysw_max) = 0.0 ! fraction of functional area, 1.0/(exp(r_DF*(1.0-accH[j]/W0[j]))+1.0)
    real :: Rring(Ysw_max) = 0.0 ! Radius to the outer edge
    real :: Lring(Ysw_max) = 0.0 ! Length of xylem conduits
    real :: Aring(Ysw_max) = 0.0 ! Area of each ring
    real :: Kring(Ysw_max) = 0.0 ! Conductance of each ring

    ! for diagnostics
    real :: Aleafmax = 0.0  ! Yearly maximum leaf area
    real :: dailyTrsp
    real :: dailyGPP   ! kgC/tree day-1
    real :: dailyNPP
    real :: dailyResp
    real :: dailyNup
    real :: annualTrsp
    real :: annualGPP ! C flux/tree
    real :: annualNPP
    real :: annualResp

    ! ---- Nitrogen model related parameters
    real :: NSNmax = 0.
    real :: NSN = 0.    ! non-structural N pool
    real :: leafN = 0.
    real :: sapwN= 0.
    real :: woodN = 0. ! N of heart wood
    real :: rootN = 0. ! N of fine roots
    real :: seedN = 0. !
    real :: N_uptake = 0.
    real :: annualNup  = 0.0
    real :: fixedN ! fixed N at each stem per tree
    real :: dailyfixedN
    real :: annualfixedN = 0.0 ! annual N fixation per unit crown area

    real :: bl_max  = 0.0 ! Max. leaf biomass, kg C/individual
    real :: br_max  = 0.0 ! Max. fine root biomass, kg C/individual
    real :: CSAsw   = 0.0
    real :: topyear = 0.0 ! the years that a plant in top layer
    real :: DBH_ys        ! DBH at the begining of a year (growing season)

    ! ---- water uptake-related variables
    real :: root_length(soil_L) ! m
    real :: rootarea ! total fine root area per tree
    real :: rootdepth  ! maximum depth of fine roots
    real :: ArootL(soil_L) = 0.0 ! Root area per layer
    real :: WupL(soil_L) = 0.0 ! normalized vertical distribution of uptake
    real :: Q_soil(soil_L) = 0.0 ! Soil to roots water flux (kg H2O/tree/step)
    real :: W_supply  ! potential water uptake rate per unit time per tree
    real :: transp   ! transpiration rate per tree per hour
    real :: uptake_frac(soil_L) ! for LM3 soil water uptake, Weng, 2017-10-28
    real :: K_r,r_r
    real :: root_zeta
    ! for photosynthesis
    real :: An_op = 0.0 ! mol C/(m2 of leaf per year)
    real :: An_cl = 0.0 ! mol C/(m2 of leaf per year)
    real :: w_scale =-9999
    real :: C_growth = 0.0 ! carbon gain since last growth, kg C/individual
    real :: N_growth = 0.0 ! Nitrogen used for plant tissue growth
    real :: extinct = 0.75     ! light extinction coefficient in the canopy for photosynthesis

  end type cohort_type

  !---------------------------
  type :: vegn_tile_type
    integer :: tag ! kind of the tile
    integer :: landuse = LU_NTRL
    integer :: n_cohorts = 0
    integer :: n_initialCC = 0
    integer :: n_years   = 0
    integer :: n_canopycc = 0
    type(cohort_type), pointer :: cohorts(:)=>NULL()
    type(cohort_type), pointer :: initialCC(:)=>NULL()
    real :: area      ! m2
    real :: age = 0.0 ! tile age
    real :: LAI  ! leaf area index
    real :: CAI  ! crown area index
    real :: LAIlayer(9) = 0.0 ! LAI of each crown layer, max. 9
    real :: f_gap(9)    = 0.0 ! gap fraction of each crown layer
    real :: kp(9)       = 0.0 ! light extinction coefficient for each layer
    ! uptake-related variables
    real :: root_distance(soil_L) ! characteristic half-distance between fine roots, m
    real :: ArootL(soil_L) = 0.0 ! Root are per layer
    ! averaged quantities for PPA phenology
    real :: tc_daily = 0.0
    real :: tc_pheno = 0.0 ! smoothed canopy air temperature for phenology

    ! litter and soil carbon pools
    real :: litter = 0.0 ! litter flux
    real :: SOC(5) = 0. ! metabolicL, structuralL, microbial, fastSOM, slowSOM
    real :: SON(5) = 0.

    !!  Nitrogen pools, Weng 2014-08-08
    real :: mineralN= 0.  ! Mineral nitrogen pool, (kg N/m2)
    real :: totN    = 0.
    real :: N_input = 0.      ! annual N input (kgN m-2 yr-1)
    real :: N_uptake= 0.0  ! kg N m-2 hour-1
    real :: fixedN  = 0.0  ! kg N/step
    real :: annualN = 0.0  ! annual available N in a year
    real :: Nloss_yr= 0.0  ! annual N loss
    real :: N_P2S_yr= 0.0  ! annual N from plants to soil
    real :: previousN      ! an weighted annual available N
    real :: initialN0

    ! Soil water
    integer :: soiltype    ! lookup table for soil hydrologic parameters
    real :: FLDCAP  ! soil field capacity
    real :: WILTPT  ! soil wilting point (0.xx)
    real :: evap           ! kg m-2 per unit fast time step (mm/hour)
    real :: transp         ! kg m-2 hour-1
    real :: runoff        ! Water runoff of the veg tile, unit?
    real :: thetaS     ! moisture index (ws - wiltpt)/(fldcap - wiltpt)
    real :: wcl(soil_L)   ! volumetric soil water content for each layer
    real :: freewater(soil_L) ! Available water in each layer
    real :: psi_soil(soil_L) ! MPa
    real :: K_soil(soil_L)   ! Kg H2O/(m2 s MPa)
    real :: soilWater      ! kg m-2 in root zone

    ! Vegetation water
    real :: W_leaf
    real :: W_stem
    real :: W_dead

    ! water uptake-related variables
    real :: RAI ! root area index
    real :: RAIL(soil_L) = 0.0 ! Root length per layer, m of root/m
    real :: W_uptake  ! water uptake rate per unit time per m2

    !  Carbon fluxes
    real :: gpp =0 ! gross primary production, kgC m-2 yr-1
    real :: npp =0 ! net primary productivity
    real :: resp = 0 ! auto-respiration of plants
    real :: rh  =0 ! soil carbon lost to the atmosphere

    !  fire disturbance
    real :: C_combusted = 0.0 ! Carbon released to atmosphere via fire
    real :: treecover   = 0.0 ! tree CAI in the top layer, for fire spread
    real :: grasscover  = 0.0 ! grass CAI, for the initial fire severity
    ! daily diagnostics
    real :: dailyGPP
    real :: dailyNPP
    real :: dailyResp
    real :: dailyRh
    real :: dailyNup
    real :: dailyfixedN
    ! for annual diagnostics
    real :: dailyPrcp=0.0, annualPrcp = 0.0 ! mm m-2 yr-1
    real :: dailyTrsp=0.0,dailyEvap=0.0, dailyRoff=0.0 ! mm m-2 yr-1
    real :: annualTrsp=0.0,annualEvap=0.0, annualRoff=0.0 ! mm m-2 yr-1
    real :: annualGPP = 0.0 ! kgC m-2 ground yr-1
    real :: annualNPP = 0.0
    real :: annualResp = 0.0
    real :: annualRh   = 0.0
    real :: annualNup  = 0.0   ! accumulated N uptake kgN m-2 yr-1
    real :: annualfixedN = 0.  ! fixe N in a tile
    ! for annual reporting at tile level
    real :: NSC, SeedC, leafC, rootC, SapwoodC, WoodC
    real :: NSN, SeedN, leafN, rootN, SapwoodN, WoodN
    real :: totSeedC,totSeedN
    ! for cohort plant types (climate-vegetation relationship, Biome, LM3)
    real :: t_ann  = 0.0 ! annual mean T, degK
    real :: t_cold = 0.0 ! average temperature of the coldest month, degK
    real :: p_ann  = 0.0 ! annual mean precip
    real :: ncm    = 0.0 ! number of cold months
  end type vegn_tile_type

  ! PFT-specific parameters
  type(spec_data_type), save :: spdata(0:MSPECIES) ! define PFTs
  
  
contains

!============================ Subroutines =================================

 !===========================
 FUNCTION esat(T) ! pressure, Pa
   IMPLICIT NONE
   REAL :: esat
   REAL, INTENT(IN) :: T ! degC
   esat=610.78*exp(17.27*T/(T+237.3))
 END FUNCTION esat

  
  
  
end module DataType