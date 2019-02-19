from astropy.io import ascii
from matplotlib import pyplot as plt
import numpy as np
from astropy.modeling.blackbody import blackbody_lambda
from astropy import units as u
from astropy import constants as c

skycalc = ascii.read("../data/radiance_comp_median.dat") ## radiance components for our median conditions, i.e. T_tel = 282 K, airmass = 1.3, pwv = 2.5 mm, telescope emissivity (all mirrors combined) = 15 %, 3-20 micron, R=2000

wave = skycalc["col1"] ## unit: micron

atm = skycalc["col5"] ## unit: ph/s/m^2/micron/arcsec2 ("lower atmosphere", "upper atmosphere" is 0)
area_ELT = np.pi*((38)**2-(11.1)**2)/4 ## unit: m^2
atm_int = atm * area_ELT ## unit: ph/s/micron/arcsec2

tel = skycalc["col8"] ## unit: ph/s/m^2/micron/arcsec2
spider_contrib = 0.1 ## fraction of flux that the spiders contribute to the total telescope background
tel_int = (1+spider_contrib) * tel * area_ELT ## unit: ph/s/micron/arcsec2

### Emissivity of METIS entrance window
T_window = 282 * u.K
f_ratio = 17.7 ## ELT Nasmyth focus
solid_angle = (1/f_ratio)**2 * u.sr
plate_scale = 3.316 *u.mm/u.arcsec ## according to ELT ICD Instruments

window_emissivity = ascii.read("../data/TC_window_METIS.dat")
we_l = window_emissivity["col1"].data
we_e = window_emissivity["col3"].data
window_radiance = we_e * blackbody_lambda(we_l*u.micron, T_window)
window_flux = (window_radiance * solid_angle * plate_scale**2).decompose()
energy_per_photon = c.h*c.c/(we_l*u.micron)
window_flux_photon = (window_flux/energy_per_photon).to("1/(micron s arcsec**2)")
## interpolate to other wavelength scale
window = np.interp(wave,we_l,window_flux_photon)

LM = [3.0,5.5]
NQ = [7.5,20.0]

LM_ix = np.where((wave>LM[0]) & (wave<LM[1]))
NQ_ix = np.where((wave>NQ[0]) & (wave<NQ[1]))

plt.subplot(211)
plt.plot(wave[LM_ix],atm_int[LM_ix],label="Sky",color="blue",linewidth=0.5)
plt.plot(wave[LM_ix],tel_int[LM_ix],label="Telescope + spiders",color="green",linewidth=0.5)
plt.plot(wave[LM_ix],window[LM_ix],label="Window",color="red",linewidth=0.5)
plt.plot(wave[LM_ix],atm_int[LM_ix]+tel_int[LM_ix]+window[LM_ix],label="Total background",color="black",linewidth=1)
plt.title("L/M bands")
plt.xlabel("Wavelength [micron]")
plt.ylabel("[ph/s/micron/arcsec^2]")
plt.yscale("log")
plt.legend(prop={'size': 6})

plt.subplot(212)
plt.plot(wave[NQ_ix],atm_int[NQ_ix],label="Sky",color="blue",linewidth=0.5)
plt.plot(wave[NQ_ix],tel_int[NQ_ix],label="Telescope + spiders",color="green",linewidth=0.5)
plt.plot(wave[NQ_ix],window[NQ_ix],label="Window",color="red",linewidth=0.5)
plt.plot(wave[NQ_ix],atm_int[NQ_ix]+tel_int[NQ_ix]+window[NQ_ix],label="Total background",color="black",linewidth=1)
plt.title("N band")
plt.xlabel("Wavelength [micron]")
plt.ylabel("[ph/s/micron/arcsec^2]")
plt.yscale("log")

plt.tight_layout()
plt.savefig("metis_background_components.pdf")
plt.clf()
