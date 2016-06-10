import aplpy
import matplotlib.pyplot as plt
import pyfits
import numpy
import pyregion
import astropy
from astropy.io import fits
from astropy import coordinates
import matplotlib.patheffects as Patheffects
from astropy import units as u

import good_cores_getsources
import good_protostars_getsources

def cross_match(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo/+catalogs/L1157.sw.final.reliable.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo_proto/+catalogs/L1157.sw.final.reliable.ok.cat', high_res_coldens_image = '/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1157/080615/cep1157_255_mu.image.resamp.fits', proto_core_indices='L1157_matched_protostellar_cores.dat'):
	# Import getsources "good core" data table
	cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
	good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
	cores_array = cores_array1[numpy.array(good_core_indices)]

	NO,    XCO_P,   YCO_P,  WCS_ACOOR,  WCS_DCOOR,  SIG_GLOB,  FG,  GOOD, SIG_MONO01, FM01, FXP_BEST01, FXP_ERRO01, FXT_BEST01, FXT_ERRO01, AFWH01, BFWH01, THEP01, SIG_MONO02, FM02, FXP_BEST02, FXP_ERRO02, FXT_BEST02, FXT_ERRO02, AFWH02, BFWH02, THEP02, SIG_MONO03, FM03, FXP_BEST03, FXP_ERRO03, FXT_BEST03, FXT_ERRO03, AFWH03, BFWH03, THEP03, SIG_MONO04, FM04, FXP_BEST04, FXP_ERRO04, FXT_BEST04, FXT_ERRO04, AFWH04, BFWH04, THEP04, SIG_MONO05, FM05, FXP_BEST05, FXP_ERRO05, FXT_BEST05, FXT_ERRO05, AFWH05, BFWH05, THEP05, SIG_MONO06, FM06, FXP_BEST06, FXP_ERRO06, FXT_BEST06, FXT_ERRO06, AFWH06, BFWH06, THEP06, SIG_MONO07, FM07, FXP_BEST07, FXP_ERRO07, FXT_BEST07, FXT_ERRO07, AFWH07, BFWH07, THEP07 = numpy.loadtxt(YSO_catalog,comments='!', unpack=True)

	good_proto_indices = good_protostars_getsources.get_good_protostars(YSO_catalog)

	#print len(good_proto_indices)

	### Loop through entire core catalog to find cores that are close to YSOs.
	### This step significantly reduces the time it would take to run through 
	### entire catalog.

	potential_matches = []
	count = 0
	for line in cores_array:
		match_counter=0
		for i in good_proto_indices:
			distance = ((line[3]-WCS_ACOOR[i])**2 + (line[4]-WCS_DCOOR[i])**2)**0.5
			if distance < 200.0/3600. and match_counter==0:
				# matched_counter prevents counting indices twice 
				# if two YSO candidates fall within getsources ellipse
				potential_matches.append(count)
				match_counter+=1
		count += 1
	#print len(potential_matches)

	### Loop through the potential matched cores identified in the step above.
	matched_cores = []
	matched_cores_proto_index = []
	for value in potential_matches:	
		line = cores_array[value]
	
		x_coor = str(line[3])
		y_coor = str(line[4])
	
		### Create a DS9 region string for the core's getsources ellipse, 
		### from which a mask will be created. 
		region = ('fk5;ellipse(' + x_coor + ', ' + y_coor + ', ' + str((line[50]/2.0)/3600.) + ', ' + 
		str((line[51]/2.0)/3600.) + ', ' + str(line[52]+90.0)+')')
	
		r = pyregion.parse(region)
		f=fits.open(high_res_coldens_image)
		mymask = r.get_mask(hdu=f[0])
		f.close()
		newmask=mymask
		### Set all values outside the core's ellipse to zero, 
		### all valus inside the ellipse are set to one.
		newmask=numpy.where(newmask==0,0,1)

		### Loop through the 70 micron point sources
		### If any fall within a core's ellipse, store that core's index
		match_counter=0
		good_protostar_index = 0
		for i in good_proto_indices:
			if newmask[int(round(YCO_P[i],0))-1][int(round(XCO_P[i],0))-1]!=0 and match_counter==0:
				matched_cores.append(value)
				matched_cores_proto_index.append(good_protostar_index)
				# matched_counter prevents counting indices twice 
				# if two YSO candidates fall within getsources ellipse
				match_counter+=1
				good_protostar_index+=1
			else: 
				good_protostar_index+=1
		#print len(matched_cores)

	#print matched_cores
	#print matched_cores_proto_index

	### Save the protostellar core indices to a file
	#numpy.savetxt(proto_core_indices, matched_cores, fmt='%i')
	# Return indices of matched cores and protostars from the "good" lists of each
	# i.e., The indices are the index from the "good" lists and not the unprocessed original, total getsources list
	return numpy.column_stack((numpy.array(matched_cores), numpy.array(matched_cores_proto_index)))

#cross_match()

#cross_match(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/DS9_regions/L1172/core_SED/L1172_good_sources.dat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1172/120115_flat/combo_proto/+catalogs/L1172.sw.final.reliable.ok.cat', high_res_coldens_image = '/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1172/082315/cep1172_255_mu.image.resamp.fits', proto_core_indices='L1172_matched_protostellar_cores.dat')

#cross_match(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/DS9_regions/L1228/core_SED/L1228_good_sources.dat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1228/120115_flat/combo_proto/+catalogs/L1228.sw.final.reliable.ok.cat', high_res_coldens_image = '/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1228/082315/cep1228_255_mu.image.resamp.fits', proto_core_indices='L1228_matched_protostellar_cores.dat')

#cross_match(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/DS9_regions/L1241/core_SED/L1241_good_sources.dat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1241/120115_flat/combo_proto/+catalogs/L1241.sw.final.reliable.ok.cat', high_res_coldens_image = '/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1241/071415/cep1241_255_mu.image.resamp.fits', proto_core_indices='L1241_matched_protostellar_cores.dat')

#cross_match(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/DS9_regions/L1251/core_SED/L1251_good_sources.dat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1251/120115_flat/combo_proto/+catalogs/L1251.sw.final.reliable.ok.cat', high_res_coldens_image = '/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1251/082315/cep1251_255_mu.image.resamp.fits', proto_core_indices='L1251_matched_protostellar_cores.dat')
