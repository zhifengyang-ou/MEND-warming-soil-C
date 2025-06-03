module PhotosynthesisMod
  !use
  use DataType
  
  implicit none
  
 ! public:: 
  
 ! private::

contains
  =============== Plant physiology =======================================
  ! Weng 2017-10-18:compute stomatal conductance, photosynthesis and respiration
  ! updates cc%An_op and cc%An_cl, from LM3
  subroutine vegn_photosynthesis (forcing, vegn)
    !arguments
    type(climate_data_type),intent(in):: forcing
    type(vegn_tile_type), intent(inout) :: vegn

    !----- local var --------------
    type(cohort_type),pointer :: cc
    real :: rad_top  ! downward radiation at the top of the canopy, W/m2
    real :: rad_net  ! net radiation absorbed by the canopy, W/m2
    real :: Tair, TairK     ! air temperature, degC and degK
    real :: cana_q   ! specific humidity in canopy air space, kg/kg
    real :: cana_co2 ! co2 concentration in canopy air space, mol CO2/mol dry air
    real :: p_surf   ! surface pressure, Pa
    real :: water_supply ! water supply per m2 of leaves
    real :: fw, fs ! wet and snow-covered fraction of leaves
    real :: psyn   ! net photosynthesis, mol C/(m2 of leaves s)
    real :: resp   ! leaf respiration, mol C/(m2 of leaves s)
    real :: tempLAI,w_scale2, transp ! mol H20 per m2 of leaf per second
    real :: kappa  ! light extinction coefficient of crown layers
    real :: f_light(10)=0.0      ! light fraction of each layer
    integer :: i, layer

    !! Light supply for photosynthesis
    vegn%kp = 0.0 ! light extinction coefficient for each of the layers
    do i = 1, vegn%n_cohorts
       cc => vegn%cohorts(i)
       layer = Max (1, Min(cc%layer,9))
       ! Calculate kappa according to sun zenith angle ! kappa = cc%extinct/max(cosz,0.01) !
       vegn%kp(layer) = vegn%kp(layer)  &  ! -0.75
                      + cc%extinct * cc%Acrown * cc%nindivs
    enddo

    ! Light fraction
    f_light = 0.0
    f_light(1) = 1.0
    do i =2, layer !MIN(int(vegn%CAI+1.0),9)
        f_light(i) = f_light(i-1)  &
             * (exp(0.-vegn%kp(i-1)*vegn%LAIlayer(i-1)) + vegn%f_gap(i-1)) !flag: make sure LAIlayer is accumulated or not; zheng
    enddo

    ! Photosynthesis
    do i = 1, vegn%n_cohorts
       cc => vegn%cohorts(i)
       associate ( sp => spdata(cc%species) ) ! spdata(:) (type(spec_data_type)) from use DataType 
         if(cc%status == LEAF_ON .and. cc%Aleaf > 1.0E-4) then ! growth status of plant: 1 for ON, 0 for OFF ! Aleaf: total area of leaves, m2/individual
           ! Convert forcing data
            layer = Max (1, Min(cc%layer,9)) !make sure less or equal 9 layers
            rad_top = f_light(layer) * forcing%radiation ! downward radiation at the top of the canopy, W/m2 ! flag: is this top or bottom of the canopy; zheng
            rad_net = f_light(layer) * forcing%radiation * 0.9 ! net radiation absorbed by the canopy, W/m2
            p_surf  = forcing%P_air  ! Pa
            TairK   = forcing%Tair ! K
            Tair    = forcing%Tair - 273.16 ! degC
            cana_q  = (esat(Tair)*forcing%RH*mol_h2o)/(p_surf*mol_air)  ! air specific humidity, kg/kg
            cana_co2= forcing%CO2 ! co2 concentration in canopy air space, mol CO2/mol dry air
           ! recalculate the water supply to mol H20 per m2 of leaf per second
            water_supply = cc%W_supply/(cc%Aleaf*step_seconds*mol_h2o) ! mol m-2 leafarea s-1

           fw = 0.0; fs = 0.0
           call gs_Leuning(rad_top, rad_net, TairK, cana_q, cc%lai, &
                           p_surf, water_supply, cc%species, sp%pt, &
                           cana_co2, cc%extinct, fs+fw, cc%layer,   &
                           psyn, resp,w_scale2,transp ) ! output
           ! put into cohort data structure for future use in growth
           cc%An_op  = psyn  ! molC s-1 m-2 of leaves
           cc%An_cl  = -resp  ! molC s-1 m-2 of leaves
           cc%w_scale  = w_scale2
           cc%transp = transp * mol_h2o * cc%Aleaf * step_seconds  ! Transpiration (kgH2O/(tree step), Weng, 2017-10-16
           cc%gpp  = (psyn-resp) * mol_C * cc%Aleaf * step_seconds ! kgC step-1 tree-1

           !if(isnan(cc%gpp))cc%gpp=0.0
           if(isnan(cc%gpp))stop '"gpp" is a NaN'
           if(isnan(cc%transp))then
              write(*,*)'w_scale2,transp,lai',w_scale2,transp,cc%lai
              stop '"transp" is a NaN'
           endif
         else
           ! no leaves means no photosynthesis and no stomatal conductance either
           cc%An_op  = 0.0
           cc%An_cl  = 0.0
           cc%gpp    = 0.0
           cc%transp = 0.0
           cc%w_scale  = -9999
         endif
       end associate
    enddo ! vegn, go through all cohorts
  end subroutine vegn_photosynthesis



  !============================================================================
  subroutine gs_Leuning(rad_top, rad_net, tl, ea, lai, &
                     p_surf, ws, pft, pt, ca, kappa, f_w, layer, &
                     apot, acl,w_scale2, transp)
    real,    intent(in)   :: rad_top ! PAR dn on top of the canopy, w/m2
    real,    intent(in)   :: rad_net ! PAR net on top of the canopy, w/m2
    real,    intent(in)   :: tl   ! leaf temperature, degK
    real,    intent(in)   :: ea   ! specific humidity in the canopy air, kg/kg
    real,    intent(in)   :: lai  ! leaf area index
    !real,    intent(in)   :: leafage ! age of leaf since budburst (deciduos), days
    real,    intent(in)   :: p_surf ! surface pressure, Pa
    real,    intent(in)   :: ws   ! water supply, mol H20/(m2 of leaf s)
    integer, intent(in)   :: pft  ! species
    integer, intent(in)   :: pt   ! physiology type (C3 or C4)
    real,    intent(in)   :: ca   ! concentartion of CO2 in the canopy air space, mol CO2/mol dry air
    real,    intent(in)   :: kappa! canopy extinction coefficient (move inside f(pft))
    real,    intent(in)   :: f_w ! fraction of leaf that's wet or snow-covered
    integer, intent(in)   :: layer  ! the layer of this canopy
    ! note that the output is per area of leaf; to get the quantities per area of
    ! land, multiply them by LAI
    !real,    intent(out)   :: gs   ! stomatal conductance, m/s
    real,    intent(out)   :: apot ! net photosynthesis, mol C/(m2 s)
    real,    intent(out)   :: acl  ! leaf respiration, mol C/(m2 s)
    real,    intent(out)   :: w_scale2,transp  ! transpiration, mol H20/(m2 of leaf s)

    ! ---- local vars
    ! photosynthesis
    real :: vm
    real :: kc,ko  ! Michaelis-Menten constants for CO2 and O2, respectively
    real :: ci
    real :: capgam ! CO2 compensation point
    real :: f2,f3
    real :: coef0,coef1
    real :: Resp

    ! conductance related
    real :: gs
    real :: b
    real :: ds  ! humidity deficit, kg/kg
    real :: hl  ! saturated specific humidity at the leaf temperature, kg/kg
    real :: do1

    ! misceleneous
    real :: dum2
    real, parameter :: light_crit = 0
    real, parameter :: gs_lim = 0.25
    real, parameter :: Rgas = 8.314 ! J mol-1 K-1, universal gas constant
    ! new average computations
    real :: lai_eq;
    real, parameter :: rad_phot = 0.0000046 ! PAR conversion factor of J -> mol of quanta
    real :: light_top
    real :: par_net
    real :: Ag
    real :: An
    real :: Ag_l
    real :: Ag_rb
    real :: anbar
    real :: gsbar
    real :: w_scale
    real, parameter :: p_sea = 1.0e5 ! sea level pressure, Pa
    ! soil water stress
    real :: Ed, an_w, gs_w

    b=0.01
    do1=0.09 ! kg/kg ! flag: what is this ?  zheng vapor pressure deficit? probably the parameter D0 in equation 12.20
    if (pft < 2) do1=0.15 ! flag: pft<2 are C3, C4 GRASSES? zheng

    ! Convert Solar influx from W/(m^2s) to mol_of_quanta/(m^2s) PAR,
    ! empirical relationship from McCree is light=rn*0.0000046
    light_top = rad_top*rad_phot
    par_net   = rad_net*rad_phot

    ! Humidity deficit, kg/kg
    call qscomp(tl, p_surf, hl)
    ds = max(hl - ea,0.0)

    associate ( sp => spdata(pft) )
      !  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*p_sea/p_surf
      !  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*p_sea/p_surf
      !  vm=sp%Vmax*exp(3000.0*(1.0/288.2-1.0/tl))

      ! Weng, 2013-01-10 
      ko=0.248    * exp(35948/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10 !ko25 in unit mol/mol
      kc=0.000404 * exp(59356/Rgas*(1.0/298.2-1.0/tl))*p_sea/p_surf ! Weng, 2013-01-10 !kc25 in unit mol/mol
      vm=sp%Vmax*exp(24920/Rgas*(1.0/298.2-1.0/tl)) ! / ((layer-1)*1.0+1.0) ! Ea = 33920

      !decrease Vmax due to aging of temperate deciduous leaves
      !(based on Wilson, Baldocchi and Hanson (2001)."Plant,Cell, and Environment", vol 24, 571-583)
      !! Turned off by Weng, 2013-02-01, since we can't trace new leaves
      !  if (sp%leafage_tau>0 .and. leafage>sp%leafage_onset) then
      !     vm=vm*exp(-(leafage-sp%leafage_onset)/sp%leafage_tau)
      !  endif

      ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986
      capgam=0.5*kc/ko*0.21*0.209 ! Farquhar & Caemmerer 1982 ! CO2 compensation point

      ! Find respiration for the whole canopy layer
      !  Resp=sp%gamma_resp*vm*lai /((layer-1)*1.0+1.0)  ! Weng, 2013-01-17 add '/ ((layer-1)*1.0+1.0)'

      ! 2014-09-03, for Nitrogen model: resp = D*(A + B*LMA)
      ! (A+B*LMA) = LNA, D=Vmax/LNA = 25E-6/0.0012 = 0.02 for a standard deciduous species
      !! Leaf resp as a function of nitrogen
      !  Resp=sp%gamma_resp*0.04*sp%LNA  & ! basal rate, mol m-2 s-1
      !       * exp(24920/Rgas*(1.0/298.2-1.0/tl))         & ! temperature scaled
      !       * lai                                        & ! whole canopy
      !       /((layer-1)*1.0+1.0)                         !
      !! as a function of LMA
      !  Resp=(sp%gamma_LNbase*sp%LNbase+sp%gamma_LMA*sp%LMA)  & ! basal rate, mol m-2 s-1
      !  Resp=sp%gamma_LNbase*(2.5*sp%LNA-1.5*sp%LNbase)     & ! basal rate, mol m-2 s-1
      Resp= sp%gamma_LN/seconds_per_year          & ! per seconds,  mol m-2 s-1
           * sp%LNA * lai / mol_c                 & ! whole canopy, mol m-2 s-1
           * exp(24920/Rgas*(1.0/298.2-1.0/tl))     ! temperature scaled

      ! Temperature effects : looks like fth in Bonan (thermal deactivation; needs to check)
      Resp = Resp / ((1.0+exp(0.4*(5.0-tl+TFREEZE))) &
                   * (1.0+exp(0.4*(tl-45.0-TFREEZE))))

      ! ignore the difference in [CO2] near the leaf and in the canopy air, rb=0.
      Ag_l  = 0.
      Ag_rb = 0.
      Ag    = 0.
      anbar = -Resp/lai
      gsbar = b
      ! find the LAI level at which gross photosynthesis rates are equal. flag: why carry out this step? zheng
      !flag: why not directly derived light fraction and shaded fraction zheng
      ! only if PAR is positive
      if(light_top > light_crit)then
         if(pt==PT_C4) then ! C4 species pt: photosynthesis pathway
            coef0=(1+ds/do1)/sp%m_cond;
            ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
            if (ci>capgam) then
               f2=vm
               f3=18000.0*vm*ci ! 18000 or 1800?
               dum2=min(f2,f3)

               ! find LAI level at which rubisco limited rate is equal to light limited rate
               lai_eq = -log(dum2/(kappa*sp%alpha_ps*light_top))/kappa
               lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

               ! gross photosynthesis for light-limited part of the canopy 
               Ag_l   = sp%alpha_ps * par_net     &
                      * (exp(-lai_eq*kappa)-exp(-lai*kappa)) &
                      / (1-exp(-lai*kappa))

               ! gross photosynthesis for rubisco-limited part of the canopy (flag: light fraction? even so, one still needs to calculate light-limited part which maybe bigger than this
               Ag_rb  = dum2*lai_eq
               Ag=(Ag_l+Ag_rb)/ &
                 ((1.0+exp(0.4*(5.0-tl+TFREEZE))) &
                 *(1.0+exp(0.4*(tl-45.0-TFREEZE))))
               An=Ag-Resp
               anbar=An/lai

               if(anbar>0.0) then
                   gsbar=anbar/(ci-capgam)/coef0;
               endif
            endif ! ci>capgam
         else ! C3 species
            coef0=(1+ds/do1)/sp%m_cond;
            coef1=kc*(1.0+0.209/ko);
            ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
            f2=vm*(ci-capgam)/(ci+coef1);
            f3=vm/2.;
            dum2=min(f2,f3);
            if (ci>capgam) then
               ! find LAI level at which rubisco limited rate is equal to light limited rate
               lai_eq=-log(dum2*(ci+2.*capgam)/(ci-capgam)/ &
                           (sp%alpha_ps*light_top*kappa))/kappa
               lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

               ! gross photosynthesis for light-limited part of the canopy
               Ag_l   = sp%alpha_ps              &
                    * (ci-capgam)/(ci+2.*capgam) * par_net   &
                    * (exp(-lai_eq*kappa)-exp(-lai*kappa))  &
                    / (1.0-exp(-lai*kappa))
               ! gross photosynthesis for rubisco-limited part of the canopy
               Ag_rb  = dum2*lai_eq
               Ag = (Ag_l+Ag_rb) /((1.0+exp(0.4*(5.0-tl+TFREEZE))) &
                  * (1.0+exp(0.4*(tl-45.0-TFREEZE))))
               An = Ag - Resp
               anbar = An/lai
               if(anbar>0.0) then
                 gsbar=anbar/(ci-capgam)/coef0
               endif
            endif ! ci>capgam
         endif
      endif ! light is available for photosynthesis

      an_w=anbar
      if (an_w > 0.) then
         an_w=an_w*(1-sp%ps_wet*f_w)
      endif
      gs_w = 1.56 * gsbar *(1-sp%ps_wet*f_w) !Weng: 1.56 for H2O?
      if (gs_w > gs_lim) then
          if(an_w > 0.) an_w = an_w*gs_lim/gs_w
          gs_w = gs_lim
      endif
    end associate
    ! find water availability diagnostic demand
    Ed = gs_w * ds*mol_air/mol_h2o ! ds*mol_air/mol_h2o is the humidity deficit in [mol_h2o/mol_air]
    ! the factor mol_air/mol_h2o makes units of gs_w and humidity deficit ds compatible:
    if (Ed>ws) then
       w_scale=ws/Ed
       gs_w=w_scale*gs_w
       if(an_w > 0.0) an_w = an_w*w_scale
       if(an_w < 0.0.and.gs_w >b) gs_w=b
    endif
    gs   = gs_w
    apot = an_w
    acl  = -Resp/lai
    transp = min(ws,Ed) ! mol H20/(m2 of leaf s)

    ! Convert units of stomatal conductance to m/s from mol/(m2 s) by
    ! multiplying it by a volume of a mole of gas
    gs = gs * Rugas * Tl / p_surf

    ! for reporting
    w_scale2=min(1.0,ws/Ed)

    ! Error check
    if(isnan(transp))then
      write(*,*)'ws,ed',ws,ed
      stop '"transp" is a NaN'
    endif
  end subroutine gs_Leuning






end module PhotosynthesisMod