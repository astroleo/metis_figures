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
tel_int = tel * area_ELT ## unit: ph/s/micron/arcsec2

window_emissivity = ascii.read("../data/TC_window_METIS.dat")
window_radiance = window_emissivity["col3"].data * blackbody_lambda(window_emissivity["col1"].data*u.micron,282*u.K)
f_ratio = 17.7
solid_angle = (1/f_ratio)**2 * u.sr
plate_scale = 3.316 *u.mm/u.arcsec
window_flux = window_radiance * solid_angle * plate_scale**2
energy_per_photon = c.h*c.c/(window_emissivity["col1"].data*u.micron)
window_flux_photon = (window_flux/energy_per_photon).to("1/(micron s arcsec**2)") ## unit: ph/s/micron/arcsec2
window_flux_photon_interp = np.interp(wave,window_emissivity["col1"],window_flux_photon)

plt.plot(wave,atm_int,label="Sky")
plt.plot(wave,tel_int,label="Telescope")
plt.plot(wave,window_flux_photon_interp,label="Window")
plt.plot(wave,atm_int+tel_int+window_flux_photon_interp,label="Total background")
plt.xlim([3.0,4.2])
plt.ylim([0,6e10])
plt.legend()
plt.title("L band")
plt.xlabel("Wavelength [micron]")
plt.ylabel("Specific intensity [ph/s/micron/arcsec^2]")
plt.savefig("L.png")
plt.clf()

plt.plot(wave,atm_int,label="Sky")
plt.plot(wave,tel_int,label="Telescope")
plt.plot(wave,window_flux_photon_interp,label="Window")
plt.plot(wave,atm_int+tel_int+window_flux_photon_interp,label="Total background")
plt.xlim([4.4,5.5])
plt.ylim([0,2e12])
plt.legend()
plt.title("M band")
plt.xlabel("Wavelength [micron]")
plt.ylabel("Specific intensity [ph/s/micron/arcsec^2]")
plt.savefig("M.png")
plt.clf()

plt.plot(wave,atm_int,label="Sky")
plt.plot(wave,tel_int,label="Telescope")
plt.plot(wave,window_flux_photon_interp,label="Window")
plt.plot(wave,atm_int+tel_int+window_flux_photon_interp,label="Total background")
plt.xlim([7.5,13.5])
plt.ylim([0,1e13])
plt.legend()
plt.title("N band")
plt.xlabel("Wavelength [micron]")
plt.ylabel("Specific intensity [ph/s/micron/arcsec^2]")
plt.savefig("N.png")
plt.clf()

plt.plot(wave,atm_int,label="Sky")
plt.plot(wave,tel_int,label="Telescope")
plt.plot(wave,window_flux_photon_interp,label="Window")
plt.plot(wave,atm_int+tel_int+window_flux_photon_interp,label="Total background")
plt.xlim([15,20])
plt.ylim([0,3e13])
plt.legend()
plt.title("Q band")
plt.xlabel("Wavelength [micron]")
plt.ylabel("Specific intensity [ph/s/micron/arcsec^2]")
plt.savefig("Q.png")
plt.clf()
