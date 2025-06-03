## added one big leaf canopy GPP model (FBEM.F90) to MEND through MEND_IN.F90

1. the inputs and control variables (.nml) are in folder photosyn

2. Plotting figures with python script and need to fix by removing the space in df[' An']

3. run with long-term climate, lai data from Biocon

    a. copied from Genie Folder
    
    b. fix auto-compile and link (see runScript.txt)
    
    c. **due to no soil water data, I treat it as 0.35 and comment out read(13...), all the above code are ineffective** (water is used for ecosystem respiration)
    
    d. uploaded it to github under ecosys_mend; will create another branch (nfeedback) for nitrogen feedback to photosynthesis
    
            1) borrow functions from soil water limitation on gpp from teco (not so good)
            2) try exp(ax)/(1+exp(ax))

4. overhaul the code for better integration between GPP model and MEND
   1. match MEND weather data (19980101-20091231) with photosynthesis weather data (2001.01.01-2009.12.31) <DONE>
   2. The modeled GPP is not comparable with GPP data (sINI): 0.1 (1.0) <it is now due to the difference in time scale: hourly vs daily>
   3. the other modeled variables (soil c, flux) are comparable with original MEND
   4. now add plant functional groups (C3,C4 and legume) <on it> 
      1. select namelist for C3, C4, legume not working; one error is that an extra line is needed at the end (why? not sure) <done; namelist is like datatypes which need to be defined at the begining of the program>
      2. adding namelists for the four treatment combinations to control LAI scaling and functional proportions <added in namelist(FBEM_namelist.nml)> now the code is running
      3. check GPP proportion <done:output the three functional group GPPs in subroutine can_photosyn>
   5. add nitrogren fixation by symbiotic N fixation (lines 1585-1591 MOD_MEND.F90)
      1. open and write data: line 468,469 in MOD_MEND.F90; close: line 471