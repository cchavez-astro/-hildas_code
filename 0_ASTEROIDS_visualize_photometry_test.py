#%%

import matplotlib.pyplot as plt
import numpy as np
import time as t
from astropy.io import fits
from photutils import aperture_photometry, CircularAperture, CircularAnnulus
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

#Title
print()
print("Checking code for Aperture Photometry variables")
print()
KernelX=45
KernelY=45
delta=25
year=t.strftime("%Y")
month=t.strftime("%B")
day=t.strftime("%d")
print("Date is :", day, month, year)

#Setting directories, Larissa and Schubart:
direc0L = 'C:/Users/Cristian/Documents/ASTRO/astrobservations/'

#Selection (INPUT FOR RUNNING THE CODE -l OR S- )
folder=input("Name of folder? ")
round_number=input("Number of the round? ")

direc=direc0L+folder+'/1.fits'

#Opening the FITS file
f = fits.open(direc)
image = f[0].data
image = image[200:300,200:320]
headr = f[0].header

#Wavelength, Calibration factor and date from fits header
wavelength = headr['WAVELNTH'] 
date = headr['DATE-OBS']
factor = headr['CALFCTR']  

if wavelength == 11.1: 
   r=3.6
    
if wavelength == 19.7: 
   r=3.8
        
if wavelength == 31.5: 
   r=4.3
       
if wavelength == 34.8: 
   r=4.5
    
print("print(image.shape) (y,x)", image.shape)
    
limX=image.shape[1]
limY=image.shape[0]
         
print("Limit x, y values are", limX, limY)
    
#redefining size of image array
kernel = image[KernelY:KernelY+delta,KernelX:KernelX+delta]
limXk=kernel.shape[1]
limYk=kernel.shape[0]

#Input parameters for Aperture Photometry    
f1=float(input("Input factor for r_in (r_in x factor): "))
f2=float(input("Input factor for r_out (r_out x factor): "))
r_in=f1*r
r_out=f2*r


#METHOD 1 FOR FINDING THE MAIN SOURCE ****************   
peak  = np.max(kernel)
peak_loc = np.argmax(kernel[np.where(kernel != 0)])
peakpos = np.unravel_index(peak_loc, kernel.shape, order='C')

yloc = (peakpos[0])+KernelY
xloc = (peakpos[1])+KernelX
    
print()
print("Peak is ", peak/factor)
print("peak_loc is", peak_loc)   
print("Coordinates are ", [xloc],[yloc])

#Setting the support image for Aperture Photometry Method 1
plt.title('Larissa FITS APhot. Method 1 indiv (date vs lambda)')
plt.xlabel(date)  
plt.ylabel(wavelength) 
plt.imshow(image, cmap='gray')
   
plt.plot([xloc],[yloc] , 'ro')

#Aperture Photometry variables and functions
positions = (xloc, yloc)
apertures = CircularAperture(positions, r)
annulus_apertures = CircularAnnulus(positions, r_in, r_out)
    
apertures.plot(color='blue', lw=2.0, alpha=0.5)
annulus_apertures.plot(color='blue', lw=2.0, alpha=0.5)
    

plt.savefig(direc0+'draw/'+name+im_number+'_'+round_number+'_'+day+month+year+"_"+'_M1_manual'+'.png',format='png',dpi=300)
plt.clf()


#METHOD 2 FOR FINDING THE MAIN SOURCE ****************   
factor_threshold=float(input("Input factor_threshold "))
mean, median, std = sigma_clipped_stats(image, sigma=3, iters=5)    
print("sigma_clipped_stats values: Mean, median, std", mean, median, std)  
fwhm_var=5
threshold_var=factor_threshold*std
daofinder = DAOStarFinder(fwhm=fwhm_var, threshold=threshold_var)    
sources = daofinder(image) 

print(sources)    
positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r)
norm = ImageNormalize(stretch=SqrtStretch())

#Setting the support image for Aperture Photometry Method 2
plt.title('Larissa FITS APhot. Method 2 indiv (std vs fact_threshold)')
plt.xlabel(std)   
plt.ylabel(factor_threshold) 
plt.imshow(image, cmap='seismic')
apertures.plot(color='orange', lw=5, alpha=5)

plt.savefig(direc+'draw/'+folder+month+year+"_"+'_M2_manual'+'.png',format='png',dpi=300)

#Last picture
plt.matshow(image)

#%%