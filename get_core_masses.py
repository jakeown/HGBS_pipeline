import numpy
from pylab import *
from scipy.optimize import curve_fit
from scipy.stats import mstats
import sys
import matplotlib.pyplot as plt
from astropy import wcs
from astropy import coordinates
from astropy import units as u

import good_cores_getsources
import good_protostars_getsources
import proto_core_cross_match
import make_catalog
import make_CMF
import make_DS9_region
import make_plots

def core_mass_fits(region_name = 'L1157', cloud_name = 'Cepheus', distance = 325, getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo/+catalogs/L1157.sw.final.reliable.ok.cat', getsources_additional_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo/+catalogs/L1157.sw.final.reliable.add.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo_proto/+catalogs/L1157.sw.final.reliable.ok.cat', YSO_additional_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo_proto/+catalogs/L1157.sw.final.reliable.add.ok.cat', high_res_coldens_image = '/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1157/080615/cep1157_255_mu.image.resamp.fits', SED_figure_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1157/L1157_core_SED/', CSAR_catalog = '/mnt/scratch-lustre/jkeown/DS9_regions/L1157/CSAR/CEPl1157_CSAR.dat', Dunham_YSOs_file = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/Dunham_YSOs.dat'):
	# These are the values in each column of "cores_array" and "protostar_array"
	#NO,    XCO_P,   YCO_P,  WCS_ACOOR,  WCS_DCOOR,  SIG_GLOB,  FG,  GOOD, SIG_MONO01, FM01, FXP_BEST01, FXP_ERRO01, FXT_BEST01, FXT_ERRO01, AFWH01, BFWH01, THEP01, SIG_MONO02, FM02, FXP_BEST02, FXP_ERRO02, FXT_BEST02, FXT_ERRO02, AFWH02, BFWH02, THEP02, SIG_MONO03, FM03, FXP_BEST03, FXP_ERRO03, FXT_BEST03, FXT_ERRO03, AFWH03, BFWH03, THEP03, SIG_MONO04, FM04, FXP_BEST04, FXP_ERRO04, FXT_BEST04, FXT_ERRO04, AFWH04, BFWH04, THEP04, SIG_MONO05, FM05, FXP_BEST05, FXP_ERRO05, FXT_BEST05, FXT_ERRO05, AFWH05, BFWH05, THEP05, SIG_MONO06, FM06, FXP_BEST06, FXP_ERRO06, FXT_BEST06, FXT_ERRO06, AFWH06, BFWH06, THEP06, SIG_MONO07, FM07, FXP_BEST07, FXP_ERRO07, FXT_BEST07, FXT_ERRO07, AFWH07, BFWH07, THEP07

	# These are the values in each column of the "additional" "cores_array" and "protostar_array"
	# NO    XCO_P   YCO_P PEAK_SRC01 PEAK_BGF01 CONV_SRC01 CONV_BGF01 PEAK_SRC02 PEAK_BGF02 CONV_SRC02 CONV_BGF02 PEAK_SRC03 PEAK_BGF03 CONV_SRC03 CONV_BGF03 PEAK_SRC04 PEAK_BGF04 CONV_SRC04 CONV_BGF04 PEAK_SRC05 PEAK_BGF05 CONV_SRC05 CONV_BGF05 PEAK_SRC06 PEAK_BGF06 CONV_SRC06 CONV_BGF06 PEAK_SRC07 PEAK_BGF07 CONV_SRC07 CONV_BGF07

	# Import the raw getsources core catalog (includes "bad" sources that don't pass HGBS selection criteria)
	cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	# Find the indices of the "good" cores that pass HGBS selection criteria
	good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
	# Create new array of only the "good" cores to be used for our analysis
	cores_array = cores_array1[numpy.array(good_core_indices)]

	# Import the raw "additional" getsources catalog
	additional_cores_array1 = numpy.loadtxt(getsources_additional_catalog,comments='!')
	# Create another new array of only the "good" cores from the getsources "additional" catalog
	additional_cores_array = additional_cores_array1[numpy.array(good_core_indices)]
	
	# Import the raw getsources protostar catalog
	protostar_array1 = numpy.loadtxt(YSO_catalog,comments='!', unpack=True)
	# Find the indices of the "good" protostars that pass HGBS selection criteria
	good_proto_indices = good_protostars_getsources.get_good_protostars(YSO_catalog)
	# Create new array of only the "good" protostars to be used for our analysis
	# Not sure why the transpose below is needed, but it seems the protostar_array1
	# is imported transposed.  This doesn't seem to be the case for the core_array1
	protostar_array = protostar_array1.T[numpy.array(good_proto_indices)]
	# Pull the 70micron flux and error columns to be used later
	flux_70um = protostar_array[:,12]
	flux_70um_err = protostar_array[:,13]

	# Import the raw getsources protostar catalog
	additional_protostar_array1 = numpy.loadtxt(YSO_additional_catalog,comments='!', unpack=True)
	# Create another new array of only the "good" YSOs from the getsources "additional" catalog
	additional_protostar_array = additional_protostar_array1.T[numpy.array(good_proto_indices)]

	# Cross-match the good cores and protostars arrays to find protostars that fall within a core's ellipse
	print "Cross-matching Core and Protostar Catalogs:"
	protostellar_core_indices = proto_core_cross_match.cross_match(getsources_core_catalog, YSO_catalog, high_res_coldens_image)
	# The first column of protostellar_core_indices is the core's index in cores_array
	cross_matched_core_indices = numpy.array(protostellar_core_indices[:,0])	
	# The second column of protostellar_core_indices is the protostar's index in protostar_array
	cross_matched_proto_indices = numpy.array(protostellar_core_indices[:,1])

	# Replace the 70 micron flux measurements of protostellar cores with the protostar extraction measurements
	if len(cross_matched_core_indices)>0:
		for i,j in zip(cross_matched_core_indices, cross_matched_proto_indices):
			cores_array[i][10]=protostar_array[j][10]
			cores_array[i][11]=protostar_array[j][11]
			cores_array[i][12]=protostar_array[j][12]
			cores_array[i][13]=protostar_array[j][13]

			additional_cores_array[i][3]=additional_protostar_array[j][3]
			additional_cores_array[i][4]=additional_protostar_array[j][4]
			additional_cores_array[i][5]=additional_protostar_array[j][5]
			additional_cores_array[i][6]=additional_protostar_array[j][6]	
 
	# Calculate the deconvolved core radii
	AFWH05 = cores_array[:,50]
	BFWH05 = cores_array[:,51]
	A = numpy.float64(((((AFWH05)/60.)/60.)*numpy.pi)/180.) #radians
	A1 = numpy.float64(numpy.tan(A/2.)*2.*distance) #pc
	B = numpy.float64(((((BFWH05)/60.)/60.)*numpy.pi)/180.) #radians
	B1 = numpy.float64(numpy.tan(B/2.)*2.*distance) #pc
	FWHM_mean = mstats.gmean([A1,B1])
	HPBW = numpy.float64(((((18.2)/60.)/60.)*numpy.pi)/180.) #radians
	HPBW1 = numpy.float64(numpy.tan(HPBW/2.)*2.*distance) #pc
	R_deconv = ((FWHM_mean**2.0) - (HPBW1**2.0))**0.5 #pc

	R_deconv = numpy.where(((FWHM_mean**2.0) - (HPBW1**2.0))<=0., FWHM_mean, R_deconv)

	# Calculate the Bonnor-Ebert masses of each core based on their R_deconvolved  
	c_s = 0.2 #km/s at 10K
	G = 4.302*10**-3 #pc/M_solar (km/s)^2
	M_BE = (2.4*R_deconv*(c_s**2.))/G #M_solar
	
	M_BE = numpy.where(((FWHM_mean**2.0) - (HPBW1**2.0))<0, 9999., M_BE)

	# Define a function that produces a Flux given wavelength, Temp, and Mass
	# We will input wavelength then find T and M using least squares minimization below
	def core_mass(wavelength, T, M):
		#wavelength input in microns, Temp in Kelvin, Mass in M_solar
		#returns S_v (i.e., Flux) in units of Jy  
		D = distance #parsecs to cloud
		wavelength_mm = numpy.array(wavelength)*10.**-3.
		exponent = 1.439*(wavelength_mm**-1)*((T/10.)**-1)
		aaa = (0.12*(numpy.exp(exponent)-1.0))**-1.0
		bbb = (0.1*((numpy.array(wavelength)/300.)**-2.0))/0.01
		ccc = (D/100.)**-2
		ddd = wavelength_mm**-3.
		return M*aaa*bbb*ccc*ddd

	# Define another function that calculates Mass directly from wavelength, Temp, and Flux
	# This will be used to find the Mass of cores that don't have reliable least-squares fits
	def core_mass_from_flux(wavelength, T, S_v):
		#wavelength input in microns, Temp in Kelvin, Mass in M_solar
		#returns S_v (i.e., Flux) in units of Jy  
		D = distance #parsecs to cloud
		wavelength_mm = wavelength*10.**-3.
		exponent = 1.439*(wavelength_mm**-1)*((T/10.)**-1)
		aaa = 0.12*(numpy.exp(exponent)-1.0)
		bbb = ((0.1*((wavelength/300.)**-2.0))/0.01)**-1.0
		ccc = (D/100.)**2.0
		ddd = wavelength_mm**3.0
		return S_v*aaa*bbb*ccc*ddd

	# Define another function that calculates Mass uncertainty due to temp
	def core_mass_err_dT(wavelength, T, S_v, dT):
		#wavelength input in microns, Temp in Kelvin, Mass in M_solar
		#returns S_v (i.e., Flux) in units of Jy  
		D = distance #parsecs to cloud
		wavelength_mm = wavelength*10.**-3.
		exponent = 1.439*(wavelength_mm**-1)*((T/10.)**-1)
		aaa = 0.12*(numpy.exp(exponent))
		bbb = ((0.1*((wavelength/300.)**-2.0))/0.01)**-1.0
		ccc = (D/100.)**2.0
		ddd = wavelength_mm**3.0
		eee = 1.439*10*(wavelength_mm**-1)*(T**-2.)
		return S_v*aaa*bbb*ccc*ddd*dT*eee

	# Define another function that calculates Mass uncertainty due to flux
	def core_mass_err_dS_v(wavelength, T, S_v, dS_v):
		#wavelength input in microns, Temp in Kelvin, Mass in M_solar
		#returns S_v (i.e., Flux) in units of Jy  
		D = distance #parsecs to cloud
		wavelength_mm = wavelength*10.**-3.
		exponent = 1.439*(wavelength_mm**-1)*((T/10.)**-1)
		aaa = 0.12*(numpy.exp(exponent)-1.0)
		bbb = ((0.1*((wavelength/300.)**-2.0))/0.01)**-1.0
		ccc = (D/100.)**2.0
		ddd = wavelength_mm**3.0
		return aaa*bbb*ccc*ddd*dS_v

	# Create some empty arrays to which we will append accepted values
	Temps = []
	Masses = []
	Temps_err = []
	Masses_err = []
	counter=0
	not_accepted_counter = []
	protostellar_core_counter = 0
	mean_dust_Temp = 10.0
	mean_dust_Temp_err = 4.0
	# Designate the initial Temp and Mass guess for the least-squares fitting below
	guess = (20.0, 0.7)
	
	print "Begin SED-fitting:"
	# Loop through all the "good" cores
	for NO,    XCO_P,   YCO_P,  WCS_ACOOR,  WCS_DCOOR,  SIG_GLOB,  FG,  GOOD, SIG_MONO01, FM01, FXP_BEST01, FXP_ERRO01, FXT_BEST01, FXT_ERRO01, AFWH01, BFWH01, THEP01, SIG_MONO02, FM02, FXP_BEST02, FXP_ERRO02, FXT_BEST02, FXT_ERRO02, AFWH02, BFWH02, THEP02, SIG_MONO03, FM03, FXP_BEST03, FXP_ERRO03, FXT_BEST03, FXT_ERRO03, AFWH03, BFWH03, THEP03, SIG_MONO04, FM04, FXP_BEST04, FXP_ERRO04, FXT_BEST04, FXT_ERRO04, AFWH04, BFWH04, THEP04, SIG_MONO05, FM05, FXP_BEST05, FXP_ERRO05, FXT_BEST05, FXT_ERRO05, AFWH05, BFWH05, THEP05, SIG_MONO06, FM06, FXP_BEST06, FXP_ERRO06, FXT_BEST06, FXT_ERRO06, AFWH06, BFWH06, THEP06, SIG_MONO07, FM07, FXP_BEST07, FXP_ERRO07, FXT_BEST07, FXT_ERRO07, AFWH07, BFWH07, THEP07 in cores_array:
		not_accepted = False
		N_SED_counter = 0
	
		fig = plt.figure()
			
		# Proceed with starless core SED-fitting procedure
		# Determine how many bands the core has significant flux measurements
		if SIG_MONO01 > 5.0:
			N_SED_counter+=1
		if SIG_MONO02 > 5.0:
			N_SED_counter+=1
		if SIG_MONO04 > 5.0:
			N_SED_counter+=1
		if SIG_MONO06 > 5.0:
			N_SED_counter+=1
		if SIG_MONO07 > 5.0:
			N_SED_counter+=1
		# If the core has more than three bands in which it is significant, 
		# and the 350um Flux is higher than the 500um flux, fit the SED
		# over the 70-500um bands
		if N_SED_counter>=3 and FXT_BEST06>FXT_BEST07:
			flux_run1 = [FXT_BEST01, FXT_BEST02, FXT_BEST04, FXT_BEST06, FXT_BEST07]
			flux_err_run1 = [FXT_ERRO01, FXT_ERRO02, FXT_ERRO04, FXT_ERRO06, FXT_ERRO07]
			
			flux_run2 = [FXT_BEST02, FXT_BEST04, FXT_BEST06, FXT_BEST07]
			flux_err_run2 = [FXT_ERRO02, FXT_ERRO04, FXT_ERRO06, FXT_ERRO07]

			wavelength_run1=[70., 160.,250.,350.,500.]
			wavelength_run2=[160.,250.,350.,500.]
			try:
				popt,pcov = curve_fit(core_mass, wavelength_run1, flux_run1, p0=guess, sigma=flux_err_run1)
			except RuntimeError:
				popt = [-9999., -9999.]
			# Perform a second round of least squares fitting (without 70um point)
			try:
				popt2, pcov2 = curve_fit(core_mass, wavelength_run2, flux_run2, p0=guess, sigma=flux_err_run2)
			except RuntimeError:
				popt2 = [-1., -1.]
			# If the best fit Mass from the two fits varies by more than a factor of 2,
			# calculate the mass from the flux of the longest wavelength with significant flux
			if (popt2[1]/popt[1]) > 2.0 or (popt2[1]/popt[1]) < 0.5:
				not_accepted = True
				not_accepted_counter.append("no_SED_fit")
				# Store fillers for Mass and T
				# Will re-calculate these values below using a median core dust temp
				# from the reliable SED fits 
				Masses.append(9999)
				Temps.append(9999)
				Temps_err.append(9999)
				Masses_err.append(9999)
			else:
				not_accepted_counter.append(" ")
				# This means the SED fit is reliable
				# Store the best-fit T and M (from second run) to corresponding arrays
				Temps.append(popt2[0])
				Masses.append(popt2[1])
				# Find the uncertainty on T and M from the square root of the diagonal 
				# terms in the covariance matrix and store in arrays
				Temps_err.append(pcov2[0][0]**0.5)
				Masses_err.append(pcov2[1][1]**0.5)
				
		else:
			# This means the SED is not reliable and we will instead get the Mass
			# directly from the longest significant wavelength's flux measurement
			not_accepted = True
			not_accepted_counter.append("no_SED_fit")
			
			# Store fillers for Mass and T
			# Will re-calculate these values below using a median core dust temp
			# from the reliable SED fits 
			Masses.append(9999)
			Temps.append(9999)
			Temps_err.append(9999)
			Masses_err.append(9999)
		# Plot the fluxes of all Herschel bands regardless of significance	
		wavelength3=[70., 160.,250.,350.,500.]
		flux3 = [FXT_BEST01, FXT_BEST02, FXT_BEST04, FXT_BEST06, FXT_BEST07]
		flux_err3 = [FXT_ERRO01, FXT_ERRO02, FXT_ERRO04, FXT_ERRO06, FXT_ERRO07]
		plt.errorbar(wavelength3,flux3,yerr=flux_err3,ecolor='r',fmt='o')
		plt.yscale('log')
		plt.xscale('log')
		plt.xlabel("Wavelength ($\mu$m)")
		plt.ylabel("Flux density (Jy)")

		# If the SED fit is reliable, also plot the best fit line to the SED points
		wavelength2=[160.,250.,350.,500.]
		if not_accepted==False:
			plt.plot(np.linspace(wavelength2[0]-40,wavelength2[3]+100, 50),core_mass(np.linspace(wavelength2[0]-40,wavelength2[3]+100, 50), popt2[0], popt2[1]), color="green")
			plt.title('Core ' + str(counter+1) + ' \n' + 'T$_{dust}$ (K) = ' + str("{0:.2f}".format(round(popt2[0],2))) + ' $\pm$ ' + str("{0:.2f}".format(round(pcov2[0][0],3))) + ', Mass (M$_\odot$) = ' + str("{0:.2f}".format(round(popt2[1],2))) + ' $\pm$ ' + str("{0:.2f}".format(round(pcov2[1][1],3))))
			# Save the plots as a PDF in specified SED_figure_directory
			# The core number (in terms of the "good" sources) is added to file name
			fig.savefig(SED_figure_directory + 'core' + str(counter+1) + '.pdf')
		
		plt.close(fig)
		# Keep track of the core's that have been fit
		counter=counter+1
		print "Core Number " + str(counter) + " of " + str(len(cores_array))

	redo_indices = numpy.where(numpy.array(Temps_err)==9999)
	print "Need to re-calculate Mass and Temp for " + str(len(redo_indices[0])) + " Cores"
	Temps = numpy.array(Temps)
	
	mean_dust_Temp = numpy.median(numpy.delete(Temps, redo_indices))
	mean_dust_Temp_err = numpy.std(numpy.delete(Temps, redo_indices))
	print "We will assume a core Temp of " + str(mean_dust_Temp) + " +/- " + str(mean_dust_Temp_err) + " for those cores"
	
	redo_indices = redo_indices[0]
	redo_counter = 0
	for NO,    XCO_P,   YCO_P,  WCS_ACOOR,  WCS_DCOOR,  SIG_GLOB,  FG,  GOOD, SIG_MONO01, FM01, FXP_BEST01, FXP_ERRO01, FXT_BEST01, FXT_ERRO01, AFWH01, BFWH01, THEP01, SIG_MONO02, FM02, FXP_BEST02, FXP_ERRO02, FXT_BEST02, FXT_ERRO02, AFWH02, BFWH02, THEP02, SIG_MONO03, FM03, FXP_BEST03, FXP_ERRO03, FXT_BEST03, FXT_ERRO03, AFWH03, BFWH03, THEP03, SIG_MONO04, FM04, FXP_BEST04, FXP_ERRO04, FXT_BEST04, FXT_ERRO04, AFWH04, BFWH04, THEP04, SIG_MONO05, FM05, FXP_BEST05, FXP_ERRO05, FXT_BEST05, FXT_ERRO05, AFWH05, BFWH05, THEP05, SIG_MONO06, FM06, FXP_BEST06, FXP_ERRO06, FXT_BEST06, FXT_ERRO06, AFWH06, BFWH06, THEP06, SIG_MONO07, FM07, FXP_BEST07, FXP_ERRO07, FXT_BEST07, FXT_ERRO07, AFWH07, BFWH07, THEP07 in cores_array[redo_indices]:
		
		# Find the longest significant wavelength and corresponding flux
		if SIG_MONO07>5.0:
			wave = 500.				
			flux_fit = FXT_BEST07
			flux_fit_err = FXT_ERRO07
		elif SIG_MONO06>5.0:
			wave = 350.
			flux_fit = FXT_BEST06
			flux_fit_err = FXT_ERRO06
		elif SIG_MONO04>5.0:
			wave = 250.
			flux_fit = FXT_BEST04
			flux_fit_err = FXT_ERRO04
		elif SIG_MONO02>5.0:
			wave = 160.
			flux_fit = FXT_BEST02
			flux_fit_err = FXT_ERRO02
		# Find the mass corresponding to that flux measurement
		# ***This uses the median of the best-fit Temps from the cores with 
		# reliable SED fits (i.e., those that pass the test above)
		Mass_fit = core_mass_from_flux(wave, mean_dust_Temp, flux_fit)
		# Can add more uncertainties (e.g., calibration, etc.) below
		Mass_error = (core_mass_err_dT(wave, mean_dust_Temp, flux_fit, mean_dust_Temp_err) + 
					core_mass_err_dS_v(wave, mean_dust_Temp, flux_fit, flux_fit_err)**2.0)**0.5
		# Store the Mass and T with uncertainties
		# Need to perform a more in-depth error analysis 
		Masses[redo_indices[redo_counter]] = Mass_fit
		Temps[redo_indices[redo_counter]] = mean_dust_Temp
		Temps_err[redo_indices[redo_counter]] = mean_dust_Temp_err
		Masses_err[redo_indices[redo_counter]] = Mass_error
		
		fig = plt.figure()
		wavelength3=[70., 160.,250.,350.,500.]
		flux3 = [FXT_BEST01, FXT_BEST02, FXT_BEST04, FXT_BEST06, FXT_BEST07]
		flux_err3 = [FXT_ERRO01, FXT_ERRO02, FXT_ERRO04, FXT_ERRO06, FXT_ERRO07]
		plt.errorbar(wavelength3,flux3,yerr=flux_err3,ecolor='r',fmt='o')
		plt.yscale('log')
		plt.xscale('log')
		plt.xlabel("Wavelength ($\mu$m)")
		plt.ylabel("Flux density (Jy)")
		plt.title('Core ' + str(redo_indices[redo_counter]+1) + ' **Not Accepted** ' + ' \n' + 'T$_{dust}$ (K) = ' + str("{0:.2f}".format(round(mean_dust_Temp,2))) + "$\pm$" + str("{0:.2f}".format(round(mean_dust_Temp_err,2))) + ', Mass (M$_\odot$) = ' + str("{0:.2f}".format(round(Mass_fit,2))) + "$\pm$" + str("{0:.2f}".format(round(Mass_error,2))))
		fig.savefig(SED_figure_directory + 'core' + str(redo_indices[redo_counter]+1) + '_NA.pdf')
		plt.close(fig)

		redo_counter+=1
		print "Re-calculating Core " + str(redo_counter) + " of " + str(len(redo_indices))

	unreliable = float(len(numpy.where(numpy.array(not_accepted_counter)=="no_SED_fit")[0]))
	print 'Fraction of Unreliable SED fits: ' + str(round((unreliable/float(len(cores_array))),3))
	not_accepted_counter = numpy.where(numpy.array(not_accepted_counter) == "no_SED_fit", not_accepted_counter, "None")
	
	#Replace nans if they exist in the M_BE array
	where_are_nans = numpy.isnan(M_BE)
	M_BE[where_are_nans] = 9999

	# Calculate the alpha_BE ratio to determine prestellar cores
	alpha_BE = numpy.array(M_BE)/numpy.array(Masses)

	alpha_BE = numpy.where(((FWHM_mean**2.0) - (HPBW1**2.0))<0, 9999., alpha_BE)
	
	# Create an array indicating a core as candidate(1)/robust(2) prestellar
	candidate_array = numpy.where(alpha_BE<=5.0, 1, 0)
	robust_candidate_array = numpy.where(alpha_BE<=2.0, 2, candidate_array)
	# Identify protostars in the array with the number (3)
	if len(cross_matched_core_indices)>0:
		robust_candidate_array[cross_matched_core_indices]=3

	# Remove protostars from the alpha_BE array and find the indices of the remaining 
	# candidate/robust prestellar cores
	robust_prestellar_indices = numpy.where(numpy.delete(alpha_BE,cross_matched_core_indices)<=2.0)
	candidate_prestellar_indices =  numpy.where(numpy.delete(alpha_BE,cross_matched_core_indices)<=5.0)

	# Remove protostars from the Mass and Radius arrays 
	# (we only want to plot starless cores in the Mass vs. Radius plot)
	R_deconv_minus_protos=numpy.delete(R_deconv, cross_matched_core_indices)
	Masses_minus_protos=numpy.delete(Masses, cross_matched_core_indices)

	# Find the final list of prestellar candidate/robust Masses with protostars removed
	prestellar_candidates = Masses_minus_protos[numpy.array(candidate_prestellar_indices[0])]
	prestellar_robust = Masses_minus_protos[numpy.array(robust_prestellar_indices[0])]
	print 'prestellar candidates: ' + str(len(prestellar_candidates))
	print 'robust prestellar candidates: ' + str(len(prestellar_robust))
	
	# Plot Mass versus Radius and save the figure
	fig = plt.figure()
	plt.scatter(R_deconv_minus_protos,Masses_minus_protos, label='starless')
	plt.scatter(R_deconv_minus_protos[numpy.array(candidate_prestellar_indices[0])],prestellar_candidates, color='red', label='candidate')
	plt.scatter(R_deconv_minus_protos[numpy.array(robust_prestellar_indices[0])],prestellar_robust, color='green', label='robust')
	plt.yscale('log')
	plt.xscale('log')
	#plt.legend()
	plt.title(region_name + ' Cores')
	plt.ylabel("Mass, M (M$_\odot$)")
	plt.xlabel("Deconvolved FWHM size, R (pc)")
	plt.xlim([10**-3, 2*10**-1])
	plt.ylim([10**-3, 10**2])
	fig.savefig(SED_figure_directory + 'mass_vs_radius_' + region_name + '.png')

	# Append the Radius, Mass, Temperature, alpha_BE, etc. arrays as columns 
	# onto the "good cores" array and save as a .dat file
	numpy.savetxt(SED_figure_directory + region_name +'_good_sources.dat', numpy.column_stack((cores_array,numpy.array(R_deconv),numpy.array(FWHM_mean),numpy.array(Masses), numpy.array(Masses_err),numpy.array(Temps), numpy.array(Temps_err),numpy.array(alpha_BE),numpy.array(robust_candidate_array))))

	# Save a text file with RA and Dec of the "good" cores
	# This is needed for the SIMBAD cross-match
	numpy.savetxt(SED_figure_directory + region_name +'_SIMBAD_RA_DEC.dat', zip(cores_array[:,3],cores_array[:,4]))
	
	# Create the catalog of good cores; includes flux measurments, positions, etc.
	make_catalog.make_catalog(region_name=region_name, cloud_name=cloud_name, distance=distance, additional_cores_array=additional_cores_array, good_cores_array = cores_array, cross_matched_core_indices=cross_matched_core_indices, cross_matched_proto_indices=cross_matched_proto_indices, alpha_BE=alpha_BE, getsources_core_catalog = getsources_core_catalog, R_deconv=R_deconv, FWHM_mean = FWHM_mean, Masses=Masses, Masses_err = Masses_err, Temps=Temps, Temps_err=Temps_err, not_accepted_counter = not_accepted_counter, CSAR_catalog = CSAR_catalog, high_res_coldens_image = high_res_coldens_image, SED_figure_directory=SED_figure_directory, Dunham_YSOs_file=Dunham_YSOs_file)

	make_CMF.CMF_plotter(region=region_name, Masses_minus_protos=Masses_minus_protos, prestellar_candidates=prestellar_candidates, prestellar_robust=prestellar_robust, SED_figure_directory=SED_figure_directory)

	# Create plot of column density PDF
	#make_plots.coldense_vs_cores(region_name=region_name, SED_figure_directory=SED_figure_directory, high_res_coldens_image=high_res_coldens_image)

	# Create histogram of background column densities for the prestellar cores 
	make_plots.bg_coldense_plotter(region_name=region_name, SED_figure_directory=SED_figure_directory)

	# Create histogram of core dust temperatures from the cores with reliable SED fits
	make_plots.core_temp_plotter(region_name=region_name, SED_figure_directory=SED_figure_directory)

	# Create histogram of core radii for the starless core population
	make_plots.core_size_plotter(region_name=region_name, SED_figure_directory=SED_figure_directory)

	# Create DS9 region files for the good core and proto catalogs at all wavelengths
	make_DS9_region.make_DS9_regions_good_cores(getsources_core_catalog=getsources_core_catalog, YSO_catalog = YSO_catalog, DS9_region_directory = SED_figure_directory, catalog_type='all_cores', cross_matched_core_indices=cross_matched_core_indices, candidate_prestellar_indices=candidate_prestellar_indices, robust_prestellar_indices=robust_prestellar_indices)
	make_DS9_region.make_DS9_regions_good_cores(getsources_core_catalog=getsources_core_catalog, YSO_catalog = YSO_catalog, DS9_region_directory = SED_figure_directory, catalog_type='proto', cross_matched_core_indices=cross_matched_core_indices, candidate_prestellar_indices=candidate_prestellar_indices, robust_prestellar_indices=robust_prestellar_indices)
	make_DS9_region.make_DS9_regions_good_cores(getsources_core_catalog=getsources_core_catalog, YSO_catalog = YSO_catalog, DS9_region_directory = SED_figure_directory, catalog_type='prestellar_candidates', cross_matched_core_indices=cross_matched_core_indices, candidate_prestellar_indices=candidate_prestellar_indices, robust_prestellar_indices=robust_prestellar_indices)
	make_DS9_region.make_DS9_regions_good_cores(getsources_core_catalog=getsources_core_catalog, YSO_catalog = YSO_catalog, DS9_region_directory = SED_figure_directory, catalog_type='prestellar_robust', cross_matched_core_indices=cross_matched_core_indices, candidate_prestellar_indices=candidate_prestellar_indices, robust_prestellar_indices=robust_prestellar_indices)
	if len(cross_matched_core_indices)>0:
		make_DS9_region.make_DS9_regions_good_cores(getsources_core_catalog=getsources_core_catalog, YSO_catalog = YSO_catalog, DS9_region_directory = SED_figure_directory, catalog_type='proto_cores', cross_matched_core_indices=cross_matched_core_indices, candidate_prestellar_indices=candidate_prestellar_indices, robust_prestellar_indices=robust_prestellar_indices)
		make_DS9_region.make_DS9_regions_good_cores(getsources_core_catalog=getsources_core_catalog, YSO_catalog = YSO_catalog, DS9_region_directory = SED_figure_directory, catalog_type='starless_cores', cross_matched_core_indices=cross_matched_core_indices, candidate_prestellar_indices=candidate_prestellar_indices, robust_prestellar_indices=robust_prestellar_indices)
	make_DS9_region.make_DS9_regions_CSAR(CSAR_catalog = CSAR_catalog, DS9_region_directory=SED_figure_directory)
	
	#plt.show()

#core_mass_fits()
#core_mass_fits(region_name ='L1172', distance = 288, getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1172/120115_flat/combo/+catalogs/L1172.sw.final.reliable.ok.cat', getsources_additional_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1172/120115_flat/combo/+catalogs/L1172.sw.final.reliable.add.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1172/120115_flat/combo_proto/+catalogs/L1172.sw.final.reliable.ok.cat', YSO_additional_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1172/120115_flat/combo_proto/+catalogs/L1172.sw.final.reliable.add.ok.cat', high_res_coldens_image = '/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1172/082315/cep1172_255_mu.image.resamp.fits', SED_figure_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1172/L1172_core_SED/', CSAR_catalog = '/mnt/scratch-lustre/jkeown/DS9_regions/L1172/CSAR/CEPl1172_CSAR.dat')
#core_mass_fits(region_name ='L1228', distance = 200, getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1228/120115_flat/combo/+catalogs/L1228.sw.final.reliable.ok.cat', getsources_additional_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1228/120115_flat/combo/+catalogs/L1228.sw.final.reliable.add.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1228/120115_flat/combo_proto/+catalogs/L1228.sw.final.reliable.ok.cat',YSO_additional_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1228/120115_flat/combo_proto/+catalogs/L1228.sw.final.reliable.add.ok.cat', high_res_coldens_image = '/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1228/082315/cep1228_255_mu.image.resamp.fits', SED_figure_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1228/L1228_core_SED/', CSAR_catalog = '/mnt/scratch-lustre/jkeown/DS9_regions/L1228/CSAR/CEPl1228_CSAR.dat')
#core_mass_fits(region_name ='L1241', distance = 300, getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1241/120115_flat/combo/+catalogs/L1241.sw.final.reliable.ok.cat', getsources_additional_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1241/120115_flat/combo/+catalogs/L1241.sw.final.reliable.add.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1241/120115_flat/combo_proto/+catalogs/L1241.sw.final.reliable.ok.cat', YSO_additional_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1241/120115_flat/combo_proto/+catalogs/L1241.sw.final.reliable.add.ok.cat', high_res_coldens_image = '/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1241/071415/cep1241_255_mu.image.resamp.fits', SED_figure_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1241/L1241_core_SED/', CSAR_catalog = '/mnt/scratch-lustre/jkeown/DS9_regions/L1241/CSAR/CEPl1241_CSAR.dat')
#core_mass_fits(region_name ='L1251', distance = 300, getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1251/120115_flat/combo/+catalogs/L1251.sw.final.reliable.ok.cat', getsources_additional_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1251/120115_flat/combo/+catalogs/L1251.sw.final.reliable.add.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1251/120115_flat/combo_proto/+catalogs/L1251.sw.final.reliable.ok.cat', YSO_additional_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1251/120115_flat/combo_proto/+catalogs/L1251.sw.final.reliable.add.ok.cat', high_res_coldens_image = '/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1251/082315/cep1251_255_mu.image.resamp.fits', SED_figure_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1251/L1251_core_SED/', CSAR_catalog = '/mnt/scratch-lustre/jkeown/DS9_regions/L1251/CSAR/CEPl1251_CSAR.dat')
