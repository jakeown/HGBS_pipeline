import numpy
from pylab import *
from scipy.optimize import curve_fit
from scipy.stats.mstats import mode
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
import pyregion
from numpy import genfromtxt
from astropy import wcs
from astropy import coordinates
import matplotlib.patheffects as Patheffects
from astropy import units as u
import aplpy

def coldense_vs_cores(region_name, SED_figure_directory, high_res_coldens_image):
	
	data = fits.getdata(high_res_coldens_image)
	#new_mask = numpy.where(data==0, 0, data)
	
	fig = plt.figure()
	histogram, bins = numpy.histogram(data, bins=numpy.logspace(numpy.log10(3.e20),numpy.log10(100.e21),51))
	root_N_err = histogram**0.5
	centers = (bins[:-1]+bins[1:])/2
	plt.errorbar(centers, histogram, yerr = root_N_err, color='blue', linestyle='None', marker='None')
	plt.plot(centers, histogram, color='blue', drawstyle='steps-mid')
	plt.semilogy()
	plt.semilogx()
	#plt.xlim([3.e20, 100.e21])
	#plt.ylim([10.,10.**7])
	plt.ylabel("Number of pixels per bin: $\Delta$N/$\Delta$logN$_{H_{2}}$")	
	plt.xlabel("Column Density, N$_{H_{2}}$ [cm$^{-2}$]")
	fig.savefig(SED_figure_directory + region_name +"_coldens_hist.png")
	plt.close()

def bg_coldense_plotter(region_name, SED_figure_directory):
	
	L1157_bg, L1157_type = numpy.loadtxt(SED_figure_directory+region_name+'_core_catalog1.dat', usecols=(57,63), dtype=[('bg',float),('type','S20')], unpack=True)

	prestellar_indices_L1157 = numpy.where(L1157_type=='prestellar')
	
	prestellar_cores_bg = L1157_bg[prestellar_indices_L1157]
	counter=0
	fig = plt.figure()
	histogram, bins = numpy.histogram(prestellar_cores_bg/10., bins=numpy.linspace(0., 30., 20))
	root_N_err = histogram**0.5
	centers = (bins[:-1]+bins[1:])/2
	plt.errorbar(centers, histogram, yerr = root_N_err, color='green', linestyle='None', markerfacecolor='none', marker='^')
	plt.plot(centers, histogram, color='green', drawstyle='steps-mid', label=region_name)

	Av_range = numpy.arange(0.,100.)
	Av_7 = numpy.zeros(len(Av_range))+(7.*0.94)
	plt.plot(Av_7, Av_range, color='black', linestyle='--')		

	plt.ylabel("Number of cores per bin: $\Delta$N/$\Delta$N$_{H_2}$")
	plt.xlabel("Background Column Density, N$_{H_2}$ [10$^{21}$ cm$^{-2}$]")
	plt.legend()
	fig.savefig(SED_figure_directory + region_name + '_cores_vs_bg_coldens.png')
	#plt.show()
	plt.close()
	#print float(len(numpy.where(prestellar_cores_bg/10.>7.*0.94)[0]))/float(len(prestellar_cores_bg))

def core_temp_plotter(region_name, SED_figure_directory):
	
	L1157_temp, L1157_type, L1157_type2 = numpy.loadtxt(SED_figure_directory+region_name+'_core_catalog2.dat', usecols=(8,17,18), dtype=[('bg',float),('type','S20'),('type2','S20')], unpack=True)
	
	starless_indices_L1157 = numpy.where(L1157_type!='protostellar')

	starless_cores_temps=[]
	for i,j in zip(L1157_temp[starless_indices_L1157], L1157_type2[starless_indices_L1157]):
		if j!="no_SED_fit":
		# Only accept the temperature if it was a reliable SED fit
			starless_cores_temps.append(i)
	
	fig = plt.figure()
	histogram, bins = numpy.histogram(starless_cores_temps, bins=numpy.linspace(7., 30., 30))
	root_N_err = histogram**0.5
	centers = (bins[:-1]+bins[1:])/2
	plt.errorbar(centers, histogram, yerr = root_N_err, color='green', linestyle='None', markerfacecolor='none', marker='^')
	plt.plot(centers, histogram, color='green', drawstyle='steps-mid', label=region_name)

	plt.ylabel("Number of cores per bin: $\Delta$N/$\Delta$T")
	plt.xlabel("Dust Temperature [K]")
	plt.legend()
	fig.savefig(SED_figure_directory + region_name + '_cores_vs_temp.png')
	#plt.show()
	plt.close()

def core_size_plotter(region_name, SED_figure_directory):
	
	L1157_radius, L1157_type, L1157_type2 = numpy.loadtxt(SED_figure_directory+region_name+'_core_catalog2.dat', usecols=(4,17,18), dtype=[('bg',float),('type','S20'),('type2','S20')], unpack=True)

	starless_indices_L1157 = numpy.where(L1157_type!='protostellar')

	prestellar_indices_L1157 = numpy.where(L1157_type=='prestellar')

	#L1157=[]
	#for i,j in zip(L1157_radius[starless_indices_L1157], L1157_type2[starless_indices_L1157]):
	#	if j!="no_SED_fit":
	#		L1157.append(i)

	fig = plt.figure()
	histogram, bins = numpy.histogram(L1157_radius[starless_indices_L1157], bins=numpy.linspace(0., 0.1, 30))
	root_N_err = histogram**0.5
	centers = (bins[:-1]+bins[1:])/2
	plt.errorbar(centers, histogram, yerr = root_N_err, color='green', linestyle='None', markerfacecolor='none', marker='^')
	plt.plot(centers, histogram, color='green', drawstyle='steps-mid', label='Cepheus')

	plt.ylabel("Number of cores per bin: $\Delta$N/$\Delta$R")
	plt.xlabel("Core Radius [pc]")
	plt.legend()
	fig.savefig(SED_figure_directory + region_name + '_cores_vs_radius.png')
	#plt.show()
	plt.close()

