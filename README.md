# HGBS_pipeline
This pipeline is meant to serve as a uniform core analysis tool for the Herschel Gould Belt Survey first-look core catalog papers.  The pipeline consists of a single Python function (get_core_masses.core_mass_fits()) that performs SED fits to HGBS getsources-detected dense cores, determines their masses, creates a core catalog that is formatted in a similar manner to Konyves et al (2015), and creates a bevy of plots to assist the analysis of the observed region.  The function can be used after core extractions have been completed on a region using both getsources and CSAR.  

If you encounter any problems with the pipeline, please let me know.  Also, if you feel you can improve it in any way, please add a pull request.

Follow the steps below when using the get_core_masses.core_mass_fits() function:

1) Perform a getsources core and protostar extraction on your region in the HGBS prescribed method.
2) Perform a CSAR extraction on your region in the HGBS prescribed method.  
3) 
