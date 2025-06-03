Spinoff BiomeE by Ensheng Weng (https://doi.org/10.5194/gmd-15-8153-2022; https://github.com/wengensheng/BiomeESS/) 

Questions are flagged as: flag 

going through vegn_photosynthesis, branch out to gs_leuning()\
going through gs_leuning()

going through the details is not a good idea! It will take forever and won't get me anywhere.

The wise thing to do is probably to comb the model logic and try to run it as wanted\
and then integrate it with MEND. 

1) name list file (./para_files/input.nml) is created in runBiomeE.x with the content of fparameter transferred to file input.nml
```
fparameter='./para_files/parameters_ORNL_test.nml'
echo $fparameter
cat $fparameter > ./para_files/input.nml
```

## the model BiomeE logic flow
- BiomeE_initialization() in BiomeE module was called first
  - in BiomeE_initialization() name list file was called first (./para_files/input.nml)\
  in the file (the parameters_ORNL_test.nml was transferred into input.nml).\
  soil_data_nml,vegn_parameters_nml,initial_state_nml; note that these namelists have been defined and made public in module datatypes.F90 
  - once the variables in namelists are assigned values from parameters_ORNL_test.nml, these values are used to initiate soilpars and species data (spdata)


Please check all the flags once done