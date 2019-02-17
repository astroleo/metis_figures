from astropy.io import ascii
from matplotlib import pyplot as plt
import numpy as np


def plot_filters_band(band_name,color):
	all_filters = ascii.read("../data/filter_sensitivity.txt")
	filters = all_filters["col1"][all_filters["col2"]==band_name]
	sensitivities = all_filters["col3"][all_filters["col2"]==band_name]
	ylim = [0,np.max(sensitivities)*1.2]
	label_offset = ylim[1]*0.05

	for filter_name, sens in zip(filters,sensitivities):
		edges = get_filter_edges(filter_name)
		center = np.sum(edges)/2
		plt.plot(center,sens,color=color,marker="o")
		plt.plot(edges,[sens,sens],color=color)
		plt.text(center,sens-label_offset,filter_name,ha="center")
	
	return(ylim)

def get_filter_edges(filter_name):
	filter_curve = "../data/filter_curves/TC_filter_"+filter_name+".dat"
	flt_trans = ascii.read(filter_curve)
	# for the wavelength edges of the filter we use the standard definition,
	#    i.e. where the transmission is 0.5 * max.
	w = flt_trans["col1"]
	t = flt_trans["col2"]
	inband = w[np.where(t>0.5*max(t))]
	edges = [inband[0],inband[-1]]
	return(edges)
	
skycalc = ascii.read("../data/transmission_comp_median.dat") ## transmission components for our median conditions, i.e. T_tel = 282 K, airmass = 1.3, pwv = 2.5 mm, telescope emissivity (all mirrors combined) = 15 %, 3-20 micron, R=2000

wave = skycalc["col1"] ## unit: micron
mol_absorption = skycalc["col2"]
ozone = skycalc["col3"]
rayleigh = skycalc["col4"]
aerosol = skycalc["col5"]
trans_total = mol_absorption * ozone * rayleigh * aerosol


def plot_trans_atmo_filters(band_name,band_pass):
	fig, ax1 = plt.subplots()
	axcolor = "blue"
	ax1.set_xlim(band_pass)
	ax1.set_xlabel('Wavelength [micron]')
	ylim = plot_filters_band(band_name,axcolor)
	ax1.set_ylim(ylim)
	ax1.set_ylabel(r'Point source sensitivity [$\mu{}$Jy/5$\sigma{}$/1h]', color=axcolor)
	ax1.tick_params('y', colors=axcolor)

	ax2 = ax1.twinx()
	axcolor="grey"
	ax2.set_xlim(band_pass)
	ax2.plot(wave,trans_total,color=axcolor,alpha=0.5)
	ax2.set_ylabel('Atmospheric transmission', color=axcolor)
	ax2.tick_params('y', colors=axcolor)

	plt.title(band_name + " band")
	plt.tight_layout()
	plt.savefig("trans_atmo_filters_"+band_name+".pdf")
	plt.clf()

Lband = [3.0,4.2]
Mband = [4.4,5.5]
Nband = [7.5,13.5]
Qband = [16,20]


plot_trans_atmo_filters("L",Lband)
plot_trans_atmo_filters("M",Mband)
plot_trans_atmo_filters("N",Nband)
plot_trans_atmo_filters("Q",Qband)