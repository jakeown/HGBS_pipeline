import numpy

import good_cores_getsources
import good_protostars_getsources

#NO,    XCO_P,   YCO_P,  WCS_ACOOR,  WCS_DCOOR,  SIG_GLOB,  FG,  GOOD, SIG_MONO01, FM01, FXP_BEST01, FXP_ERRO01, FXT_BEST01, FXT_ERRO01, AFWH01, BFWH01, THEP01, SIG_MONO02, FM02, FXP_BEST02, FXP_ERRO02, FXT_BEST02, FXT_ERRO02, AFWH02, BFWH02, THEP02, SIG_MONO03, FM03, FXP_BEST03, FXP_ERRO03, FXT_BEST03, FXT_ERRO03, AFWH03, BFWH03, THEP03, SIG_MONO04, FM04, FXP_BEST04, FXP_ERRO04, FXT_BEST04, FXT_ERRO04, AFWH04, BFWH04, THEP04, SIG_MONO05, FM05, FXP_BEST05, FXP_ERRO05, FXT_BEST05, FXT_ERRO05, AFWH05, BFWH05, THEP05, SIG_MONO06, FM06, FXP_BEST06, FXP_ERRO06, FXT_BEST06, FXT_ERRO06, AFWH06, BFWH06, THEP06, SIG_MONO07, FM07, FXP_BEST07, FXP_ERRO07, FXT_BEST07, FXT_ERRO07, AFWH07, BFWH07, THEP07 

def make_DS9_regions_good_cores(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo/+catalogs/L1157.sw.final.reliable.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo_proto/+catalogs/L1157.sw.final.reliable.ok.cat', DS9_region_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1157/', catalog_type=False, cross_matched_core_indices=False, candidate_prestellar_indices=False, robust_prestellar_indices=False):
	if catalog_type=='all_cores':
		# Import getsources "good core" data table
		cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
		good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
		cores_array = cores_array1[numpy.array(good_core_indices)]
		save_file_name = '_all_cores_good.reg'
		new_line = ' #color=green width=2\n'
	elif catalog_type=='proto':
		# Import the raw getsources protostar catalog
		protostar_array1 = numpy.loadtxt(YSO_catalog,comments='!', unpack=True)
		# Find the indices of the "good" protostars that pass HGBS selection criteria
		good_proto_indices = good_protostars_getsources.get_good_protostars(YSO_catalog)
		# Create new array of only the "good" protostars to be used for our analysis
		# Not sure why the transpose below is needed, but it seems the protostar_array1
		# is imported transposed.  This doesn't seem to be the case for the core_array1
		cores_array = protostar_array1.T[numpy.array(good_proto_indices)]
		save_file_name = '_proto_good.reg'
		new_line = ' #color=red width=2\n'
	elif catalog_type=='proto_cores':
		cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
		good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
		cores_array = cores_array1[numpy.array(good_core_indices)]
		cores_array = cores_array[numpy.array(cross_matched_core_indices)]
		save_file_name = '_protostellar_cores.reg'
		new_line = ' #color=yellow width=2\n'
	elif catalog_type=='starless_cores':
		cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
		good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
		cores_array = cores_array1[numpy.array(good_core_indices)]
		indices = numpy.arange(len(cores_array[:,0]))
		for i in cross_matched_core_indices:
			indices_minus_protostellar_cores = filter(lambda a: a!=i,indices)
			indices = indices_minus_protostellar_cores
		
		cores_array = cores_array[numpy.array(indices)]
		save_file_name = '_starless_cores.reg'
		new_line = ' #color=magenta width=2\n'
	elif catalog_type=='prestellar_candidates':
		cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
		good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
		cores_array = cores_array1[numpy.array(good_core_indices)]
		cores_array = cores_array[numpy.array(candidate_prestellar_indices[0])]
		save_file_name = '_prestellar_candidates.reg'
		new_line = ' #color=brown width=2\n'
	elif catalog_type=='prestellar_robust':
		cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
		good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
		cores_array = cores_array1[numpy.array(good_core_indices)]
		cores_array = cores_array[numpy.array(robust_prestellar_indices[0])]
		save_file_name = '_prestellar_robust.reg'
		new_line = ' #color=blue width=2\n'

	header = 'Region file format: DS9 version 4.1 \nglobal color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \nfk5'

	ellipse = []
	for i in cores_array[:,0]:
		ellipse.append('ellipse')

	number = [1, 2, 3, 4, 5, 6, 7]
	A = [cores_array[:,14], cores_array[:,23], cores_array[:,32], cores_array[:,41], cores_array[:,50], cores_array[:,59], cores_array[:,68]]
	B = [cores_array[:,15], cores_array[:,24], cores_array[:,33], cores_array[:,42], cores_array[:,51], cores_array[:,60], cores_array[:,69]]
	T = [cores_array[:,16], cores_array[:,25], cores_array[:,34], cores_array[:,43], cores_array[:,52], cores_array[:,61], cores_array[:,70]]
	f = 0
	for wavelength in number:
		string = (numpy.array(A[f])/2.)/3600. #Divide by two because DS9 takes radius rather than FWHM (or diameter) &&& convert arcseconds to degrees
		string2 = []
		for i in string:
			string2.append(str(i))
		#string3 = []
		#for i in string2:
		#	string3.append(i+'"')
		AFWHM = string2

		string = (numpy.array(B[f])/2.)/3600.
		string2 = []
		for i in string:
			string2.append(str(i))
		#string3 = []
		#for i in string2:
		#	string3.append(i+'"')
		BFWHM = string2

		# theta is given as degrees E of N in getsources, but degrees ccw from the +ve x-axis in ds9.
	
		THEP_adjusted = numpy.array(T[f]) + 90.0

		# Must divide FWHM by 2 because DS9 region files use a radius as input rather than diameter or "full width"
		data = numpy.column_stack((ellipse, cores_array[:,3], cores_array[:,4], AFWHM, BFWHM, THEP_adjusted))

		wave = ['070', '160', '165', '250', '255', '350', '500']
		
		numpy.savetxt(DS9_region_directory + wave[f]+save_file_name, data, delimiter=' ', fmt = '%s', header=header, newline=new_line)
		f = f + 1

	# Write over the third line in the created files
	# If this is not added, the "fk5" string remains commented out in the header
	files = ['070' + save_file_name,'160' + save_file_name,'165' + save_file_name, '250' + save_file_name,'255' + save_file_name,'350' + save_file_name,'500' + save_file_name]

	for i in files:
		with open(DS9_region_directory + i, 'r') as file:
			data = file.readlines()
		data[2] = "fk5 \n"
		with open(DS9_region_directory + i, 'w') as file:			
			file.writelines(data)

def make_DS9_regions_CSAR(CSAR_catalog = '/mnt/scratch-lustre/jkeown/DS9_regions/L1157/CSAR/CEPl1157_CSAR.dat', DS9_region_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1157/'):
	### Import the CSAR catalog 
	### ***MAKE SURE FIRST TWO COLUMNS OF "CSAR_catalog" FILE ARE X_POSITION AND Y_POSITION OF SOURCE IN DECIMAL DEGREES***
	CSAR_array = numpy.loadtxt(CSAR_catalog,comments='#')
	x_position = CSAR_array[:,0]
	y_position = CSAR_array[:,1]

	header = 'Region file format: DS9 version 4.1 \nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \nfk5'

	ellipse = []
	for i in x_position:
		ellipse.append('point')

	# Must divide FWHM by 2 because DS9 region files use a radius as input rather than diameter or "full width"
	data = numpy.column_stack((ellipse, x_position, y_position))
	numpy.savetxt(DS9_region_directory+'CSAR_sources.reg', data, delimiter=' ', fmt = '%s', header=header, newline=' #point=x color=green width=2\n')

	# Write over the third line in the created files
	# If this is not added, the "fk5" string remains commented out in the header
	files = [DS9_region_directory+'CSAR_sources.reg']
	
	for i in files:
		with open(i, 'r') as file:
			data = file.readlines()
		data[2] = "fk5 \n"
		with open(i, 'w') as file:			
			file.writelines(data)

def make_DS9_regions_points(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo/+catalogs/L1157.sw.final.reliable.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo_proto/+catalogs/L1157.sw.final.reliable.ok.cat', DS9_region_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1157/', catalog_type=False, cross_matched_core_indices=False, candidate_prestellar_indices=False, robust_prestellar_indices=False):
	if catalog_type=='all_cores':
		# Import getsources "good core" data table
		cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
		good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
		cores_array = cores_array1[numpy.array(good_core_indices)]
		save_file_name = '_all_cores_good.reg'
		new_line = ' #color=green width=2\n'
	elif catalog_type=='proto':
		# Import the raw getsources protostar catalog
		protostar_array1 = numpy.loadtxt(YSO_catalog,comments='!', unpack=True)
		# Find the indices of the "good" protostars that pass HGBS selection criteria
		good_proto_indices = good_protostars_getsources.get_good_protostars(YSO_catalog)
		# Create new array of only the "good" protostars to be used for our analysis
		# Not sure why the transpose below is needed, but it seems the protostar_array1
		# is imported transposed.  This doesn't seem to be the case for the core_array1
		cores_array = protostar_array1.T[numpy.array(good_proto_indices)]
		save_file_name = '_proto_good.reg'
		new_line = ' #color=red width=2\n'
	elif catalog_type=='proto_cores':
		cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
		good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
		cores_array = cores_array1[numpy.array(good_core_indices)]
		cores_array = cores_array[numpy.array(cross_matched_core_indices)]
		save_file_name = '_protostellar_cores.reg'
		new_line = ' #color=yellow width=2\n'
	elif catalog_type=='starless_cores':
		cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
		good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
		cores_array = cores_array1[numpy.array(good_core_indices)]
		indices = numpy.arange(len(cores_array[:,0]))
		for i in cross_matched_core_indices:
			indices_minus_protostellar_cores = filter(lambda a: a!=i,indices)
			indices = indices_minus_protostellar_cores
		
		cores_array = cores_array[numpy.array(indices)]
		save_file_name = '_starless_cores.reg'
		new_line = ' #color=magenta width=2\n'
	elif catalog_type=='prestellar_candidates':
		cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
		good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
		cores_array = cores_array1[numpy.array(good_core_indices)]
		cores_array = cores_array[numpy.array(candidate_prestellar_indices[0])]
		save_file_name = '_prestellar_candidates.reg'
		new_line = ' #color=red width=2 point=diamond\n'
	elif catalog_type=='prestellar_robust':
		cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	
		good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
		cores_array = cores_array1[numpy.array(good_core_indices)]
		cores_array = cores_array[numpy.array(robust_prestellar_indices[0])]
		save_file_name = '_prestellar_robust.reg'
		new_line = ' #color=blue width=2 point=diamond\n'

	header = 'Region file format: DS9 version 4.1 \nglobal color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \nfk5'

	ellipse = []
	for i in cores_array[:,0]:
		ellipse.append('point')

	number = [1, 2, 3, 4, 5, 6, 7]
	A = [cores_array[:,14], cores_array[:,23], cores_array[:,32], cores_array[:,41], cores_array[:,50], cores_array[:,59], cores_array[:,68]]
	B = [cores_array[:,15], cores_array[:,24], cores_array[:,33], cores_array[:,42], cores_array[:,51], cores_array[:,60], cores_array[:,69]]
	T = [cores_array[:,16], cores_array[:,25], cores_array[:,34], cores_array[:,43], cores_array[:,52], cores_array[:,61], cores_array[:,70]]
	f = 0
	for wavelength in number:
		string = (numpy.array(A[f])/2.)/3600. #Divide by two because DS9 takes radius rather than FWHM (or diameter) &&& convert arcseconds to degrees
		string2 = []
		for i in string:
			string2.append(str(i))
		#string3 = []
		#for i in string2:
		#	string3.append(i+'"')
		AFWHM = string2

		string = (numpy.array(B[f])/2.)/3600.
		string2 = []
		for i in string:
			string2.append(str(i))
		#string3 = []
		#for i in string2:
		#	string3.append(i+'"')
		BFWHM = string2

		# theta is given as degrees E of N in getsources, but degrees ccw from the +ve x-axis in ds9.
	
		THEP_adjusted = numpy.array(T[f]) + 90.0

		# Must divide FWHM by 2 because DS9 region files use a radius as input rather than diameter or "full width"
		data = numpy.column_stack((ellipse, cores_array[:,3], cores_array[:,4]))

		wave = ['070', '160', '165', '250', '255', '350', '500']
		
		numpy.savetxt(DS9_region_directory + wave[f]+save_file_name, data, delimiter=' ', fmt = '%s', header=header, newline=new_line)
		f = f + 1

	# Write over the third line in the created files
	# If this is not added, the "fk5" string remains commented out in the header
	files = ['070' + save_file_name,'160' + save_file_name,'165' + save_file_name, '250' + save_file_name,'255' + save_file_name,'350' + save_file_name,'500' + save_file_name]

	for i in files:
		with open(DS9_region_directory + i, 'r') as file:
			data = file.readlines()
		data[2] = "fk5 \n"
		with open(DS9_region_directory + i, 'w') as file:			
			file.writelines(data)

