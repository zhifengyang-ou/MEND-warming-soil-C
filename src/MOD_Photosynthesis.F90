!! Ref:Variation of parameters in a Flux-Based Ecosystem Model across 12 sites of terrestrial ecosystems in the conterminous USA
module photosyn
      
        implicit none
        
        private
        ! public  :: can_photosyn_treatment
        public  :: can_photosyn
        public  :: can_photosyn_pft
        public  :: NLimitationGPP
        private :: Arrhenius
        
contains
    subroutine can_photosyn(gpp,treatment)
        real(8),allocatable,intent(inout):: gpp(:)
        character(len=*),intent(in)      :: treatment 
        real(8),allocatable              :: gpp_C3(:)
        real(8),allocatable              :: gpp_C4(:)
        real(8),allocatable              :: gpp_Legume(:)
        integer                          :: iiterm
        print*,'treatment: ',treatment
        call can_photosyn_pft(gpp_C3,'C3',treatment)
        call can_photosyn_pft(gpp_C4,'C4',treatment)
        call can_photosyn_pft(gpp_Legume,'Legume',treatment)
        gpp = gpp_C3 + gpp_C4 + gpp_Legume

        !write outputs header
        open(116,file='./src/photosyn/gpp_pft.out',status='unknown')
        write(116,'(2A7,A11,A3)') "gpp_C3,","gpp_C4,","gpp_legume,","gpp"
        do iiterm=1,size(gpp)
            write(116,*) gpp_C3(iiterm),',',gpp_C4(iiterm),',',gpp_legume(iiterm),',',gpp(iiterm)
        end do
        close(116)

    end subroutine !can_photosyn

    subroutine can_photosyn_pft(gpp,pft,treatment)
      ! will use name list for key parameters
        real(8),parameter::Ox=0.21,Ca=410,Rgas=8.314 !Ca in ppm
        ! real(8):: Ta=20.0,swc=0.35,I=1500.0,RH=0.6,LAI=0.0 !driving force for test;I: absorbed PAR
        real(8):: Ta,swc,I,RH,LAI !driving force for test;I: absorbed PAR
        real(8):: Rad
        real(8):: f_Ci,Ci,Tk,alpha_q
        real(8):: Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc,Ea_Ko,Kc_25,Ko_25
        real(8):: r_JmVm
        real(8):: g1,D0,kn
        real(8):: Reco0,Q10,a1
        real(8):: Vm,Gamma_star,Kc,Ko,Jm,Jc,Je,A
        real(8):: es,D,Gs,An,Ac,Reco,NEE,gppx
        character:: header
        character(len=100) :: fp, filename_full ! folder path
        integer:: istat3,day,hour,nday,ihour,eof
        !real :: gpp(366*24)
        real(8), allocatable :: gpp(:),pft_fraction
        real(8):: LAIscalar,fLAI_C3,fLAI_C4,fLAI_legume
        character(len=*)::pft,treatment 

        ! namelist /legume_para/ f_Ci,alpha_q,Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc, &
        !     Ea_Ko,Kc_25,Ko_25,r_JmVm,g1,D0,kn,Reco0,Q10,a1

        !     open(11,file='./src/photosyn/FBEM_namelist.nml')
        !     read(11,nml=legume_para,iostat=eof)
        !     if(eof/=0) then
        !         print*,"Namelist Error",eof
        !         stop
        !     end if
        !     close(11)
        
        namelist /C3_para/ f_Ci,alpha_q,Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc, &
        Ea_Ko,Kc_25,Ko_25,r_JmVm,g1,D0,kn,Reco0,Q10,a1

        namelist /C4_para/ f_Ci,alpha_q,Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc, &
                            Ea_Ko,Kc_25,Ko_25,r_JmVm,g1,D0,kn,Reco0,Q10,a1

        namelist /legume_para/ f_Ci,alpha_q,Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc, &
                                Ea_Ko,Kc_25,Ko_25,r_JmVm,g1,D0,kn,Reco0,Q10,a1

        namelist /aCaN/ LAIscalar,fLAI_C3,fLAI_C4,fLAI_legume,fp
        namelist /aCeN/ LAIscalar,fLAI_C3,fLAI_C4,fLAI_legume
        namelist /eCaN/ LAIscalar,fLAI_C3,fLAI_C4,fLAI_legume
        namelist /eCeN/ LAIscalar,fLAI_C3,fLAI_C4,fLAI_legume
  
        open(111,file='FBEM_namelist.nml')

        ! pft_fraction = 1.0
        select case(treatment)
        case('aCaN')
        read(111,nml=aCaN)
        close(111)
        
        case('aCeN')
        read(111,nml=aCeN)
        close(111)

        case('eCaN')
        read(111,nml=eCaN)
        close(111)

        case('eCeN')
        read(111,nml=eCeN)
        close(111)
        end select

        open(111,file='FBEM_namelist.nml')! either open it again or rewind(11) 
        select case(pft)
        case('C3')
        read(111,nml=C3_para)!,iostat=istat3)
        pft_fraction = fLAI_C3
        close(111)
        
        case('C4')
        read(111,nml=C4_para)!,iostat=istat3)
        pft_fraction = fLAI_C4
        close(111)

        case('Legume')
        read(111,nml=legume_para)!,iostat=istat3)
        pft_fraction = fLAI_legume
        close(111)
        end select
        

        ! read in climate and biological data
        filename_full = trim(fp)//"tair_rh_rad_hourly.in"
        open(112,file=filename_full)!,status='old',access ='sequential',form='formatted')
		filename_full = trim(fp)//"SWC_day.dat"
        open(113,file=filename_full)
        filename_full = trim(fp)//"df_lai_daily.in"
        open(114,file=filename_full)
        read(112,*) header
        read(113,*) header
        read(114,*) header

        !write outputs header
        open(115,file='./src/photosyn/cflux.out',status='unknown')
        write(115,'(2A3,A4)')"An,","Ac,","gpp"
        
        nday=4140 ! hard wired here (flag) for days from 20010101-20091231
        allocate(gpp(24*nday))
        
        ihour=0
        do day = 1, nday
            read(113,*,IOSTAT=istat3)swc
			swc=swc/100 ! transform percentage to ratio 100% to 1
            !swc = 0.35
            read(114,*,IOSTAT=istat3)LAI 
            LAI = LAI * LAIscalar * pft_fraction! respond to treatment
            if(LAI<0) LAI=0.0
            if(istat3<0)exit

            do hour = 1, 24
                ihour = ihour+1
                read(112,*,IOSTAT=istat3)Ta,RH,Rad
                RH=RH/100.0
                I = Rad*4.6/2.0 ! convert radiation to photon (4.6) and fraction (1/2) used for photosynthesis
                Ci = f_Ci * Ca  ! Ci is intercellular CO2 concentration, Ca is atmosphere CO2 concentration
                Tk = Ta+273.15
                Vm = Vm_25 * Arrhenius(Ea_vm,Rgas,Tk)
                Gamma_star = Gamma_star25*Arrhenius(Ea_gamma,Rgas,Tk)
                Kc = Kc_25*Arrhenius(Ea_Kc,Rgas,Tk)
                Ko = Ko_25*Arrhenius(Ea_Ko,Rgas,Tk)
                Jm = r_JmVm*Vm

                Jc = Vm*(Ci-Gamma_star)/(Ci+Kc*(1.+Ox/Ko)) !  Gamma_star is the photorespiratory compensation point in the absence of dark respiration
                Je = (alpha_q*I*Jm/(sqrt(Jm**2+alpha_q**2*I**2)))*((Ci-Gamma_star)/(4*(Ci+2*Gamma_star)))  
                A = min(Jc,Je) ! minimum of a light-limited and a RuBisCO-limited assimilation rate

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

                gpp(ihour)=Ac/(10.0**7.0)*(12.0*3600.0) ! convert umolCO2m-2s-1 to mgCcm-2h-1 for direct use in MEND soil
               ! gppx=Ac/((10.0**7.0)/(12.0*24.0*3600.0))

                !write outputs 
                write(115,*)An,',',Ac,',',gpp(ihour)

            end do
        end do

        close(112)
        close(113)
        close(114)
        close(115)
        print*,pft_fraction
        print*,'LAIscalar: ',LAIscalar
    end subroutine !can_photosyn_pft


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
    
    function NLimitationGPP(mineralN) result(ans)
        real(8), intent(in):: mineralN
        real(8) :: ans
        
        ans = exp(mineralN)/(1+exp(mineralN))
    
    end function NLimitationGPP

end module !photosyn
