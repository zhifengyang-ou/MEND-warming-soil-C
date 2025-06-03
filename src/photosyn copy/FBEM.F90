module photosyn
      
        implicit none
        
        private
        
        public  :: can_photosyn
        private :: Arrhenius
        
contains

    subroutine can_photosyn(gpp)
      ! will use name list for key parameters
        real(8),parameter::Ox=0.21,Ca=410,Rgas=8.314 !Ca in ppm
        real(8):: Ta=20.0,swc=0.35,I=1500.0,RH=0.6,LAI=0.0 !driving force;I: absorbed PAR
        real(8):: Rad
        real(8):: f_Ci,Ci,Tk,alpha_q
        real(8):: Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc,Ea_Ko,Kc_25,Ko_25
        real(8):: r_JmVm
        real(8):: g1,D0,kn
        real(8):: Reco0,Q10,a1
        real(8):: Vm,Gamma_star,Kc,Ko,Jm,Jc,Je,A
        real(8):: es,D,Gs,An,Ac,Reco,NEE,gppx
        character:: header
        integer:: istat3,day,hour,nday,ihour
        !real :: gpp(366*24)
        real(8), allocatable :: gpp(:)

        namelist /main_para/ f_Ci,alpha_q,Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc, &
                             Ea_Ko,Kc_25,Ko_25,r_JmVm,g1,D0,kn,Reco0,Q10,a1

        open(11,file='./src/photosyn/FBEM_namelist.nml')
        read(11,nml=main_para)
        close(11)

        ! read in climate and biologicla data
        open(12,file='./src/photosyn/tair_rh_rad_hourly.in')!,status='old',access ='sequential',form='formatted')
        open(13,file='./src/photosyn/df_swc_daily.in')
        open(14,file='./src/photosyn/df_lai_daily.in')
        read(12,*) header
        read(13,*) header
        read(14,*) header

        !write outputs header
        open(15,file='./src/photosyn/cflux.out',status='unknown')
        write(15,*)'An,Ac,gppx'
        
        nday=365
        allocate(gpp(24*nday))
        
        ihour=0
        do day = 1, nday
        read(13,*,IOSTAT=istat3)swc
        read(14,*,IOSTAT=istat3)LAI 
        if(LAI<0) LAI=0.0
        if(istat3<0)exit

        do hour = 1, 24
        ihour = ihour+1
        read(12,*,IOSTAT=istat3)Ta,RH,Rad
        RH=RH/100.0
        I = Rad*2.0
        Ci = f_Ci * Ca
        Tk = Ta+273.15
        Vm = Vm_25 * Arrhenius(Ea_vm,Rgas,Tk)
        Gamma_star = Gamma_star25*Arrhenius(Ea_gamma,Rgas,Tk)
        Kc = Kc_25*Arrhenius(Ea_Kc,Rgas,Tk)
        Ko = Ko_25*Arrhenius(Ea_Ko,Rgas,Tk)
        Jm = r_JmVm*Vm

        Jc = Vm*(Ci-Gamma_star)/(Ci+Kc*(1.+Ox/Ko))
        Je = (alpha_q*I*Jm/(sqrt(Jm**2+alpha_q**2*I**2)))*((Ci-Gamma_star)/(4*(Ci+2*Gamma_star)))  
        A = min(Jc,Je)

        es = exp(21.382-5347.5/Tk)
        D = 0.1*es*(1-RH)

        Gs = g1*A/((Ca-Gamma_star)*(1+D/D0))
        if(LAI<=0.0) then
            An=0.0
        else
            An = Gs*(Ca-Ci)! top layer canopy photosynthesis
        end if
        Ac = An*(1.0-exp(-kn*LAI))/kn

        Reco = Reco0*Q10**(Ta/10)*(swc/(swc+a1)) ! Ta in degree C;a1:moisture coefficient at which respiration is half the maximum
        NEE = Reco-Ac
        
        gpp(ihour)=Ac/((10.0**7.0)/(12.0*24.0*3600.0))
        !gppx=Ac/((10.0**7.0)/(12.0*3600.0)) ! (12.0*24.0*3600.0)

        !write outputs 
        write(15,*)An,',',Ac,',',gppx

        end do
        end do

        close(12)
        close(13)
        close(14)
        close(15)

    end subroutine can_photosyn


!      function Jc(Vmx,Ci,Gamma_star,Kc,Ox,Ko)
!             real:: Vmx,Ci,Gamma_star,Kc,Ox,Ko
!              real:: Jc
!              Jc=Vmx*(Ci-Gamma_star)/(Ci+Kc*(1.+Ox/Ko))
!      end function Jc

    function Arrhenius(Ea,Rgas,Tk)
           real(8):: Ea,Tk
           real(8):: Rgas 
           real(8):: Arrhenius
           Arrhenius=exp((Ea/(Rgas*Tk)*(Tk/298.15-1.)))
    end function Arrhenius

end module photosyn
