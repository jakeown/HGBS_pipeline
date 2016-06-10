import numpy

def get_good_protostars(getsources_proto_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo_proto/+catalogs/L1157.sw.final.reliable.ok.cat', good_proto_indices_file = 'L1157_good_protostars_indices.dat'):

	NO,    XCO_P,   YCO_P,  WCS_ACOOR,  WCS_DCOOR,  SIG_GLOB,  FG,  GOOD, SIG_MONO01, FM01, FXP_BEST01, FXP_ERRO01, FXT_BEST01, FXT_ERRO01, AFWH01, BFWH01, THEP01, SIG_MONO02, FM02, FXP_BEST02, FXP_ERRO02, FXT_BEST02, FXT_ERRO02, AFWH02, BFWH02, THEP02, SIG_MONO03, FM03, FXP_BEST03, FXP_ERRO03, FXT_BEST03, FXT_ERRO03, AFWH03, BFWH03, THEP03, SIG_MONO04, FM04, FXP_BEST04, FXP_ERRO04, FXT_BEST04, FXT_ERRO04, AFWH04, BFWH04, THEP04, SIG_MONO05, FM05, FXP_BEST05, FXP_ERRO05, FXT_BEST05, FXT_ERRO05, AFWH05, BFWH05, THEP05, SIG_MONO06, FM06, FXP_BEST06, FXP_ERRO06, FXT_BEST06, FXT_ERRO06, AFWH06, BFWH06, THEP06, SIG_MONO07, FM07, FXP_BEST07, FXP_ERRO07, FXT_BEST07, FXT_ERRO07, AFWH07, BFWH07, THEP07 = numpy.loadtxt(getsources_proto_catalog,comments='!', unpack=True)

	#Remove bad sources using the Herschel selection criteria
	index = 0
	good_proto_indices = []
	for i in range(len(NO)):
		if SIG_MONO01[i]>5.0 and FXP_BEST01[i]>0.0 and FXT_BEST01[i]>0.0 and SIG_GLOB[i]>=1 and FXP_BEST01[i]/FXP_ERRO01[i]>1.5 and AFWH01[i]<12.6 and BFWH01[i]<12.6 and AFWH01[i]/BFWH01[i]<1.3:
			good_proto_indices.append(index)
		index += 1

	#print len(good_proto_indices)

	### Save the "good protostar" indices to a file
	#numpy.savetxt(good_proto_indices_file, good_proto_indices, fmt='%i')
	return good_proto_indices

#get_good_protostars()
#get_good_protostars(getsources_proto_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1172/120115_flat/combo_proto/+catalogs/L1172.sw.final.reliable.ok.cat', good_proto_indices_file = 'L1172_good_protostars_indices.dat')
#get_good_protostars(getsources_proto_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1228/120115_flat/combo_proto/+catalogs/L1228.sw.final.reliable.ok.cat', good_proto_indices_file = 'L1228_good_protostars_indices.dat')
#get_good_protostars(getsources_proto_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1241/120115_flat/combo_proto/+catalogs/L1241.sw.final.reliable.ok.cat', good_proto_indices_file = 'L1241_good_protostars_indices.dat')
#get_good_protostars(getsources_proto_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1251/120115_flat/combo_proto/+catalogs/L1251.sw.final.reliable.ok.cat', good_proto_indices_file = 'L1251_good_protostars_indices.dat')
