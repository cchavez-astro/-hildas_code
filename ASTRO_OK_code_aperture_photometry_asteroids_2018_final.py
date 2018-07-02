#%%
#!/usr/bin/python
#Procedure to get photometry of Larissa from SOFIA/FORCAST images

from astropy.io import ascii
from astropy.io import fits
import numpy as np
import time as t
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from astropy.table import Table, Column, MaskedColumn
from astropy.time import Time
from photutils import aperture_photometry, CircularAperture, CircularAnnulus

#Setting directories, Larissa and Schubart:
direc0L = 'C:/Users/Cristian/Documents/ASTRO/astrobservations/larissa/L/'
direc0S = 'C:/Users/Cristian/Documents/ASTRO/astrobservations/schubart/S/'

#Selection (INPUT SOURCE FOR RUNNING THE CODE -l OR S- )
op=input("Aperture Photometry for Larissa (L) or Schubart (S)? ")

if op=="L":
    direc=direc0L
    name="larissa"
else:
    direc=direc0S
    name="schubart"


#Date intro
year=t.strftime("%Y")
month=t.strftime("%B")
day=t.strftime("%d")
print("Date of today is :", day, month, year)

#parameters
mpl.rcParams['legend.numpoints'] = 1

#list of filenames
asteroid_table = ascii.read(direc+name+'_imagesLC.txt',data_start=1)
files = asteroid_table["names"].data

#Initializing lists
lambdaL=[]
peakL=[]
peak_backL=[]
xlocL=[]
ylocL=[]
dateL=[]
FactL=[]
errorL=[]
error2L=[]

#Lightcurve:
plt.figure(1)
plt.figure(figsize=(12,9))
#plt.axis((2455720.9,2455721.0,0.0,0.8))

op2=input("Customized Labelplot? (Y) or (N)? ")
if op2=="Y":
    Labelplot=input("Put the name to the Graph")
else:
    round_number=input("Number of the round? ")
    Labelplot='Lightcurve of '+name+' round_'+(round_number)
    
plt.title(Labelplot)
plt.xlabel('Time (JD)')
plt.ylabel('Integrated Flux (Jy)')  

#Preparation for the Peak calculus
KernelX=45
KernelY=45
delta=25

for i in range(len(files)):
#for i in range(3):
    print('Examining image: '+files[i])
    f = fits.open(direc+'data/'+files[i])
    LC = asteroid_table[i][1]
    print('LC is: ', LC)
    if LC == 3:
        continue
    image = f[0].data
    image = image[200:300,200:320]
    headr = f[0].header
    #Wavelength, Calibration factor and date from fits header
    wavelength = headr['WAVELNTH'] 
    factor = headr['CALFCTR']  
    errcalf = headr['ERRCALF']
    date = headr['DATE-OBS']
    t = Time(date, format='isot', scale='utc')
    LCtimes=t.jd
    
    #Kernel esta entre 250 y 290 ...recuerda poner o sacar el 2 solamente
    kernel = image[KernelY:KernelY+delta,KernelX:KernelX+delta]
    
    #Peak function
    limXk=kernel.shape[1]
    limYk=kernel.shape[0]
    peak_0  = np.max(kernel)
    print("peak_0 ", peak_0)
    peak=peak_0/factor
    print("peak ", peak)
    peak_loc = np.argmax(kernel[np.where(kernel != 0)])
    peakpos = np.unravel_index(peak_loc, kernel.shape, order='C')
    #peakpos = np.unravel_index(np.max(image), (200,200))
    yloc = (peakpos[0])+KernelY
    xloc = (peakpos[1])+KernelX 
    print("peakpos is (x,y)", xloc, yloc)
     
#UPDATED FACTORS!!! (2018)
    if wavelength == 11.1: 
        xwid = 2.8/0.4/2.355
        ywid = 2.8/0.4/2.355
        C=1
        K=1.009
        r=3.6
        err_table=0.09
        edgcolor = "red"
        dotcolor = "red"
        if LC == 1: 
            marker_var='8'
        else:
            marker_var=10

    if wavelength == 19.7: 
        xwid = 2.4/0.4/2.355
        ywid = 2.4/0.4/2.355
        C=1
        K=1.0025
        r=3.8
        err_table=0.065
        edgcolor = "red"
        dotcolor = "yellow"
        if LC == 1: 
            marker_var='D'
        else:
            marker_var=9

    if wavelength == 31.5: 
        xwid = 2.8/0.4/2.355
        ywid = 2.8/0.4/2.355
        C=1
        K=1.0005
        r=4.3
        err_table=0.1
        edgcolor = "darkgreen"
        dotcolor = "forestgreen"
        if LC == 1: 
            marker_var='*'
        else:
            marker_var=11

    if wavelength == 34.8: 
        xwid = 3.1/0.4/2.355
        ywid = 3.1/0.4/2.355
        C=1
        K=1.0015
        r=4.5
        err_table=0.09
        edgcolor = "mediumblue"
        dotcolor = "dodgerblue"
        if LC == 1: 
            marker_var='P'
        else:
            marker_var=11

        
#APERTURE PHOTOMETRY PART
 #    input("APERTURE Photometry now ...")
    print("For the image ",i+1)
    AP1=input("Do you want to enter the source location? (Y) or (N) ")
    if AP1=='Y':
        xloc=int(input("Xloc?"))
        yloc=int(input("Yloc?"))
    AP2=input("Do you want to enter the annulus aperture values? (Y) or (N) ")
    if AP2=='Y':
        r_in=float(input("Inner radius factor? "))
        r_out=float(input("Outter radius factor? "))
    else:    
        r_in=r*3
        r_out=r*4.5

    positions = (xloc, yloc)
    apertures = CircularAperture(positions, r)
    annulus_apertures = CircularAnnulus(positions, r_in, r_out)
    apers = [apertures, annulus_apertures]
    
    phot_table = aperture_photometry(image, apers)
    print("BB Phot table 2018 ")
    print(phot_table)
    
    print()    
    print("Apertures Area is / Annulus Area is", apertures.area(), annulus_apertures.area())
    bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
    bkg_sum = bkg_mean * apertures.area()

    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    phot_table['Aperture_final'] = final_sum
    print("Final sum is ", final_sum, phot_table[0][5])
    
    print(phot_table)
      
    flux_ap = (phot_table[0][5])/factor
    Fcc_ap=flux_ap/K
    
    if LC == 1: 
        Fact_ap=Fcc_ap/C
    
    #ERROR CALCULATION
    #Setting the zone of annulus
    region = image*0
    nx=99
    ny=99
    for i in range(nx): 
      for j in range(ny):
        x = i - nx/2
        y = j - ny/2
        region[i,j] = np.sqrt(x**2 + y**2)
    
    std_dev = np.std(image[(region > r_in) & (region < r_out)])/factor

    if LC == 2: 
        if Fcc_ap > 0. :
            Fact_ap=Fcc_ap/C+3*std_dev
        else:
            Fact_ap=np.median(image) + 3.*std_dev
             
    error = np.sqrt((std_dev)**2+(errcalf/factor*Fact_ap)**2+(err_table*Fact_ap)**2)
    
    print("std_dev, error", std_dev, error)
    input("STOP 1 \n")
        
    #First output of the code
    print("Lambda, factor //", wavelength, factor)
    print("flux_ap, Fcc_ap, Fact_ap //", flux_ap, Fcc_ap, Fact_ap)

    print("error, calibration factor for error ", error, errcalf)
    input("STOP 2 \n")
          
    #Graph
    plt.errorbar(LCtimes, Fact_ap, yerr=error,ecolor=edgcolor)
    plt.plot(LCtimes,Fact_ap,marker=marker_var,markeredgecolor=edgcolor,markerfacecolor=dotcolor)
    
    #List appending data
    lambdaL.append(wavelength)
    peakL.append(peak)
    xlocL.append(xloc)
    ylocL.append(yloc)
    dateL.append(LCtimes)
    FactL.append(Fact_ap)
    errorL.append(error)
    error2L.append(std_dev)
    

#Third output of the code: a table of data
#data = Table([dateL, lambdaL, peakL, FactL, xlocL, ylocL], names=['DATE (JD)', 'LAMBDA', 'PEAK', 'FLUX (Jy)', 'Xloc','Yloc'])

data = Table([dateL, lambdaL, FactL, errorL, error2L], names=['DATE (JD)','lambda',  'FLUX (Jy)', 'Flux error1', 'Std dev annulus'])

ascii.write(data, direc+'draw/'+name+'_'+'_'+day+month+year+"_"+'_Lightcurve'+'_Graph_'+round_number+'.dat', overwrite=True)
print()


#Add legend
edgcolor = "red"
dotcolor = "red"
bandA = mlines.Line2D([], [], color='none', marker='8',markeredgecolor=edgcolor,markerfacecolor=dotcolor,
                          markersize=10, label='11.1 $\mu$m')
edgcolor = "red"
dotcolor = "yellow"
bandB = mlines.Line2D([], [], color='none', marker='D',markeredgecolor=edgcolor,markerfacecolor=dotcolor,
                          markersize=10, label='19.7 $\mu$m')
edgcolor = "darkgreen"
dotcolor = "forestgreen"
bandC = mlines.Line2D([], [], color='none', marker='*',markeredgecolor=edgcolor,markerfacecolor=dotcolor,
                          markersize=10, label='31.4 $\mu$m')

edgcolor = "mediumblue"
dotcolor = "dodgerblue"
bandD = mlines.Line2D([], [], color='none', marker='P',markeredgecolor=edgcolor,markerfacecolor=dotcolor,
                          markersize=10, label='34.8 $\mu$m')

edgcolor = "red"
dotcolor = "red"
bandE = mlines.Line2D([], [], color='none', marker=10,markeredgecolor=edgcolor,markerfacecolor=dotcolor,
                          markersize=10, label='11.1 $\mu$m Upper Limit')
edgcolor = "mediumblue"
dotcolor = "dodgerblue"
bandF = mlines.Line2D([], [], color='none', marker=11,markeredgecolor=edgcolor,markerfacecolor=dotcolor,
                          markersize=10, label='34.8 $\mu$m Upper Limit')

plt.legend(handles=[bandA,bandB,bandC,bandD,bandE,bandF])

#Output plot to a file
plt.savefig(direc+'draw/'+name+'_'+'_'+day+month+year+"_"+'_Lightcurve'+'_Graph_'+round_number+'.png',format='png',dpi=300)



#print(p_fit)
    
#%%