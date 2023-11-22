from atlast_sc.calculator import Calculator
from atlast_sc.data import Data

from astropy import constants as const
from astropy import units as u

from astropy.io import fits

import numpy as np
import os

# Weather octiles
octiles = {i: 12.5*i for i in range(1,8)}
octile  = octiles[1]

# For the moment, to use the finetune option, you need to install the AtLAST Sensitivity Calculator
# from this fork: https://github.com/lucadimascolo/AtLAST_sensitivity_calculator
calc = Calculator(finetune=True)
calc.weather = octile
calc.elevation = 45*u.deg

# Generate a linear set of frequencies spanning the whole range covered by AtLAST [35-950 GHz]
nulist = np.linspace(Data.obs_frequency.lower_value,
                     Data.obs_frequency.upper_value,1000)

nuband = np.zeros(nulist.shape)
nusens = np.zeros(nulist.shape)

nufwhm = np.zeros(nulist.shape)

# Target spectral resolution
R = 1000

for ni, nu in enumerate(nulist):
  print(ni+1,nulist.shape[0],nu)
  
# Update the sensitivity calculator
  calc.obs_freq  =  nu*u.GHz
  calc.bandwidth = (nu/R)*u.GHz

  freq = calc.user_input.obs_freq.value
  dish = calc.instrument_setup.dish_radius.value

# Compute AtLAST beam at frequency nu  
  fwhm = ((1.20*const.c/(freq*(2.00*dish))).to(u.dimensionless_unscaled)*u.radian).to(u.deg)
  area = np.pi*(fwhm**2)/4.00/np.log(2.00)

# Calculate beam-level sensitivity
  sens = calc.calculate_sensitivity(1*u.h)
  
  nusens[ni] = sens.to(u.Jy).value
  nuband[ni] = (nu/R)

  nufwhm[ni] = fwhm.to(u.deg).value

os.system('mkdir -p output')

outfits = fits.BinTableHDU.from_columns([fits.Column('nu',         format = 'D', unit = 'GHz', array = nulist),
                                         fits.Column('bandwidth',  format = 'D', unit = 'GHz', array = nuband),
                                         fits.Column('RMS_noise',  format = 'D', unit =  'Jy', array = nusens),
                                         fits.Column('FWHM_beam',  format = 'D', unit = 'deg', array = nufwhm)])
outfits.header['octile'] =  octile
outfits.header['elev']   = (calc.elevation,'degree')
outfits.header['R']      = (R,'spectral resolution')

outfits.writeto('output/AtLAST_sensitivity_table_R-{0:04d}_octile-{1}_el-{2}.fits'.format(R,octile,calc.elevation.to(u.deg).value),overwrite=True)
