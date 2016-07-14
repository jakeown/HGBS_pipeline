import aplpy
import matplotlib.pyplot as plt
import numpy
import pyregion
import astropy
from astropy.io import fits
import matplotlib.patheffects as Patheffects

import good_cores_getsources

#NO,    XCO_P,   YCO_P,  WCS_ACOOR,  WCS_DCOOR,  SIG_GLOB,  FG,  GOOD, SIG_MONO01, FM01, FXP_BEST01, FXP_ERRO01, FXT_BEST01, FXT_ERRO01, AFWH01, BFWH01, THEP01, SIG_MONO02, FM02, FXP_BEST02, FXP_ERRO02, FXT_BEST02, FXT_ERRO02, AFWH02, BFWH02, THEP02, SIG_MONO03, FM03, FXP_BEST03, FXP_ERRO03, FXT_BEST03, FXT_ERRO03, AFWH03, BFWH03, THEP03, SIG_MONO04, FM04, FXP_BEST04, FXP_ERRO04, FXT_BEST04, FXT_ERRO04, AFWH04, BFWH04, THEP04, SIG_MONO05, FM05, FXP_BEST05, FXP_ERRO05, FXT_BEST05, FXT_ERRO05, AFWH05, BFWH05, THEP05, SIG_MONO06, FM06, FXP_BEST06, FXP_ERRO06, FXT_BEST06, FXT_ERRO06, AFWH06, BFWH06, THEP06, SIG_MONO07, FM07, FXP_BEST07, FXP_ERRO07, FXT_BEST07, FXT_ERRO07, AFWH07, BFWH07, THEP07

def auto_zoomed_cores(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo/+catalogs/L1157.sw.final.reliable.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1157/120115_flat/combo_proto/+catalogs/L1157.sw.final.reliable.ok.cat', core_figure_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1157/core_figures/', Herschel_70um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1157/080615/cepL1157_070_mu_emission.image.resamp.fits', Herschel_160um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1157/080615/cepL1157_160_mu_emission.image.resamp.fits', Herschel_250um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1157/080615/cep1157_250_mu.image.resamp.fits', Herschel_350um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1157/080615/cep1157_350_mu.image.resamp.fits', Herschel_500um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1157/080615/cep1157_500_mu.image.resamp.fits', Herschel_coldens_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1157/080615/cep1157_255_mu.image.resamp.fits', distance=325.):

	# Import getsources "good core" data table
	cores_array1 = numpy.loadtxt(getsources_core_catalog,comments='!')
	good_core_indices = good_cores_getsources.get_good_cores(getsources_core_catalog)
	cores_array = cores_array1[numpy.array(good_core_indices)]

	core_center_RA = cores_array[:,3]
	core_center_Dec = cores_array[:,4]
	
	core_index = 1
	for line in cores_array:

		x_coor = str(line[3])
		y_coor = str(line[4])
		AFWHM_array = [line[14], line[23], line[41], line[59], line[68], line[50]] # 70,160,250,350,500,coldens
		BFWHM_array = [line[15], line[24], line[42], line[60], line[69], line[51]]
		Theta_array = [line[16], line[25], line[43], line[61], line[70], line[52]]
		SIG_array = [line[8], line[17], line[35], line[53], line[62], line[44]]
		images_array = [Herschel_70um_image, Herschel_160um_image, Herschel_250um_image, 
				Herschel_350um_image, Herschel_500um_image, Herschel_coldens_image]
		wavelengths = ['070', '160', '250', '350', '500', '255']
		maximums = numpy.zeros(len(AFWHM_array))
		minimums = numpy.zeros(len(AFWHM_array))
		counter = 0
		for i in wavelengths:
			### Create a DS9 region string for the core's getsources ellipse, 
			### from which a mask will be created. 
			region = ('fk5;box(' + x_coor + ', ' + y_coor + ', ' + str(0.05) + ', ' + 
			str(0.05) + ', ' + str(0.0)+')')
			r = pyregion.parse(region)

			f = fits.open(images_array[counter])
			header_primary = fits.getheader(images_array[counter])
			data = fits.getdata(images_array[counter])

			mymask = r.get_mask(hdu=f[0])
			newmask=numpy.where(mymask!=0)
			maximums[counter] = max(data[newmask])
			minimums[counter] = min(data[newmask])
			f.close()
			newmask2=numpy.where(mymask==0, 0, data)
			fits.writeto(i + '_contour_mask.fits', newmask2, header_primary, clobber=True)
			counter+=1
	
		microns = ['mu_070', 'mu_160', 'mu_250', 'mu_350', 'mu_500', 'mu_255']
		v_max = maximums
		v_min = minimums
		contour_levels = [(maximums[0]*0.3,maximums[0]*0.5,maximums[0]*0.7,maximums[0]*0.9),
		(maximums[1]*0.3,maximums[1]*0.5,maximums[1]*0.7,maximums[1]*0.9),
		(maximums[2]*0.3,maximums[2]*0.5,maximums[2]*0.7,maximums[2]*0.9),
		(maximums[3]*0.3,maximums[3]*0.5,maximums[3]*0.7,maximums[3]*0.9),
		(maximums[4]*0.3,maximums[4]*0.5,maximums[4]*0.7,maximums[4]*0.9),
		(maximums[5]*0.3,maximums[5]*0.5,maximums[5]*0.7,maximums[5]*0.9)]

		fig = plt.figure(figsize=(8,5.5))

		f=0
		for i in wavelengths:
			microns[f]=aplpy.FITSFigure(i + '_contour_mask.fits', figure=fig, subplot=(2,3,f+1))
			microns[f].recenter(line[3],line[4],0.025)
			microns[f].show_contour(i + '_contour_mask.fits', colors=('white', 'white', 'white','grey'), levels=contour_levels[f], overlap=True)
			microns[f].show_colorscale(cmap='gist_stern', vmin=v_min[f], vmax=v_max[f])
			#microns[f].show_regions('/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1157/L1157_core_SED/255_all_cores_good.reg')
			microns[f].add_colorbar()
			if f==2 or f==5:
				microns[f].colorbar.set_axis_label_text('MJy/sr')
				
			scale = str('{:.2f}'.format(2*(distance*numpy.tan(numpy.radians(0.5*(1.0/60.))))))
			if f==0:
				microns[f].add_scalebar((1./60.), path_effects=[Patheffects.withStroke(linewidth=3, foreground='white')])
				microns[f].scalebar.show(1./60.)  # length in degrees (one arcminute)
				microns[f].scalebar.set_corner('top left')
				microns[f].scalebar.set_label('1' + "'" + ' = ' + scale + ' pc')
			if SIG_array[f]>5.:
				line_type = 'solid'
			else:
				line_type = 'dashed'	
			microns[f].show_ellipses(line[3], line[4], AFWHM_array[f]/3600., BFWHM_array[f]/3600., Theta_array[f]+90., color='#00FF00', zorder=10, linewidth=1.0, linestyle=line_type)
			#microns[f].show_markers(line[3], line[4], c='#00FF00', marker='x', zorder=10, linewidths=1.5,s=20)
			#microns[f].show_markers(line[3], line[4], c='black', marker='x', zorder=9, linewidths=2.0,s=30)
			microns[f].tick_labels.hide() 
			microns[f].axis_labels.hide()
			#microns[f].set_title('L' + name)
			if i=='255':
				microns[f].add_label(0.3,0.1,'N$_{H2}$',relative=True, path_effects=[Patheffects.withStroke(linewidth=3, foreground='white')])
			else:
				microns[f].add_label(0.3,0.1,i + ' $\mu$' + 'm',relative=True, path_effects=[Patheffects.withStroke(linewidth=3, foreground='white')])
			f=f+1
		fig.subplots_adjust(hspace=0.05, wspace=0.35)
		fig.canvas.draw()
		fig.suptitle('Core Number ' + str(core_index) + ' - RA: ' + str(line[3]) + ' Dec: ' + str(line[4]))
		fig.savefig(core_figure_directory + 'core' +str(core_index) + '.pdf')
		print 'Core Number = ' + str(core_index)
		for i in range(6):
			microns[i].close()
		core_index = core_index + 1
		#plt.show()

#auto_zoomed_cores()
#auto_zoomed_cores(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1172/120115_flat/combo/+catalogs/L1172.sw.final.reliable.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1172/120115_flat/combo_proto/+catalogs/L1172.sw.final.reliable.ok.cat', core_figure_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1172/core_figures/', Herschel_70um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1172/082315/cepL1172_070_mu_emission.image.resamp.fits', Herschel_160um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1172/082315/cepL1172_160_mu_emission.image.resamp.fits', Herschel_250um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1172/082315/cep1172_250_mu.image.resamp.fits', Herschel_350um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1172/082315/cep1172_350_mu.image.resamp.fits', Herschel_500um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1172/082315/cep1172_500_mu.image.resamp.fits', Herschel_coldens_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1172/082315/cep1172_255_mu.image.resamp.fits', distance=288.)
#auto_zoomed_cores(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1228/120115_flat/combo/+catalogs/L1228.sw.final.reliable.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1228/120115_flat/combo_proto/+catalogs/L1228.sw.final.reliable.ok.cat', core_figure_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1228/core_figures/', Herschel_70um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1228/082315/cepL1228_070_mu_emission.image.resamp.fits', Herschel_160um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1228/082315/cepL1228_160_mu_emission.image.resamp.fits', Herschel_250um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1228/082315/cep1228_250_mu.image.resamp.fits', Herschel_350um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1228/082315/cep1228_350_mu.image.resamp.fits', Herschel_500um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1228/082315/cep1228_500_mu.image.resamp.fits', Herschel_coldens_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1228/082315/cep1228_255_mu.image.resamp.fits', distance=200.)
#auto_zoomed_cores(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1241/120115_flat/combo/+catalogs/L1241.sw.final.reliable.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1241/120115_flat/combo_proto/+catalogs/L1241.sw.final.reliable.ok.cat', core_figure_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1241/core_figures/', Herschel_70um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1241/071415/cepL1241_070_mu_emission.image.resamp.fits', Herschel_160um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1241/071415/cepL1241_160_mu_emission.image.resamp.fits', Herschel_250um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1241/071415/cep1241_250_mu.image.resamp.fits', Herschel_350um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1241/071415/cep1241_350_mu.image.resamp.fits', Herschel_500um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1241/071415/cep1241_500_mu.image.resamp.fits', Herschel_coldens_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1241/071415/cep1241_255_mu.image.resamp.fits', distance=300.)
#auto_zoomed_cores(getsources_core_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1251/120115_flat/combo/+catalogs/L1251.sw.final.reliable.ok.cat', YSO_catalog = '/mnt/scratch-lustre/jkeown/Getsources/Extract/cep1251/120115_flat/combo_proto/+catalogs/L1251.sw.final.reliable.ok.cat', core_figure_directory = '/mnt/scratch-lustre/jkeown/DS9_regions/HGBS_pipeline/L1251/core_figures/', Herschel_70um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1251/082315/cepL1251_070_mu_emission.image.resamp.fits', Herschel_160um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1251/082315/cepL1251_160_mu_emission.image.resamp.fits', Herschel_250um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1251/082315/cep1251_250_mu.image.resamp.fits', Herschel_350um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1251/082315/cep1251_350_mu.image.resamp.fits', Herschel_500um_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1251/082315/cep1251_500_mu.image.resamp.fits', Herschel_coldens_image='/mnt/scratch-lustre/jkeown/Getsources/Prepare/Images/cep1251/082315/cep1251_255_mu.image.resamp.fits', distance=300.)
