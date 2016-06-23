# HGBS_pipeline
This pipeline is meant to serve as a uniform core analysis tool for the Herschel Gould Belt Survey first-look core catalog papers.  The pipeline consists of a single Python function (get_core_masses.core_mass_fits()) that performs SED fits to HGBS getsources-detected dense cores, determines their masses, creates a core catalog that is formatted in a similar manner to Konyves et al (2015), and creates a bevy of plots to assist the analysis of the observed region.  The function can be used after core extractions have been completed on a region using both getsources and CSAR.  

If you encounter any problems with the pipeline, please let me know.  Also, if you feel you can improve it in any way, please add a pull request.

Follow the steps below when using the get_core_masses.core_mass_fits() function:

1) Perform a getsources core and protostar extraction on your region in the HGBS prescribed method.

2) Perform a CSAR extraction on your region in the HGBS prescribed method.  

3) Make sure the following Python packages are installed: numpy, pylab, scipy.optimize, scipy.stats, matplotlib, astropy, aplpy
    (I recommend using the Ureka Python environment, which is what I used to test the scripts)

4) open a terminal and move to the directory containing the HGBS_pipline scripts

5) open iPython session and import get_core_masses

6) insert the proper variable names that the function requires and run the function (example below using L1172 region):
      
      get_core_masses.core_mass_fits(
          
          region_name = 'L1172', # Name of the region observed 
          
          distance = 288., # Distance to region (pc) 
          
          getsources_core_catalog = 'L1172.sw.final.reliable.ok.cat', # getsources core extraction final catalog
          
          getsources_additional_catalog = 'L1172.sw.final.reliable.add.ok.cat', # getsources core extraction "additional" catalog
          
          YSO_catalog = 'YSO_catalog.cat', # getsources protostar extraction final catalog
          
          YSO_additional_catalog = 'YSO_catalog_add.cat', # getsources protostar extraction "additional" catalog
          
          high_res_coldens_image = 'cep1172_255_mu.image.resamp.fits', # High resolution column density image produced by prepareobs
          
          SED_figure_directory = '/Users/jkeown/Desktop/HGBS_pipeline-master/L1172_core_SED/', # Path to the directory where the figures created by the function will be stored
          
          CSAR_catalog = 'CEPl1172_CSAR.dat' # CSAR extraction catalog
          
          ***MAKE SURE FIRST TWO COLUMNS OF "CSAR_catalog" FILE ARE X_POSITION AND Y_POSITION OF SOURCE IN DECIMAL DEGREES*** 
          )
          
After the function runs, it will produce PDF files of the SEDs for each core identified by getsources that passes the HGBS detection criteria established in Konyves et al (2015).  The cores that have "_NA" appearing in the file title have unreliable SED fits and their masses were estimated based on the flux from the longest significant wavelength.  

The function also produces a CMF for the region (region_name_CMF.png), core SED temperature histogram (region_name_cores_vs_temp.png), core radius histogram (region_name_cores_vs_radius.png), a core mass/radius diagram (mass_vs_radius_region_name.png), etc.

DS9 region files for the starless, protostellar, and prestellar core populations are also created.  

The final core catalogs are saved as "region_name_core_catalog1.dat", which is similar to Table A.1 in Konyves et al (2015), and "region_name_core_catalog1.dat", which is similar to Table A.2 in Konyves et al (2015).
          
