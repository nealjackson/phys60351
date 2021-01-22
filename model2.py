import numpy as np,matplotlib, astropy, os, sys
from astropy.io.fits import getheader, getdata
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord; from astropy import units as u
infile, firstfile = 'c2.txt', 'first_simple.txt'  # names of coord file, FIRST

def tohms (a): # routine to convert deg -> string hour, min, sec
    ahr = int(a); a = 60.0*(a-ahr); amin = int(a); asec = 60.*(a-amin)
    return str(ahr), str(amin), str(asec)

# from email about wget (using the os.system, otherwise use the wget package)
def first_download (ra,dec,outfile,gif=0,fits=1,imsize=2.0,\
                    imserver='third.ucllnl.org'):
    rahr, ramin, rasec = tohms (ra/15.)   # convert ra (hr=deg*15) to dms
    decdeg, decmin, decsec = tohms (dec)  # ditto for declination
    command = ('wget -O %s "http://%s/cgi-bin/firstimage?RA=%s%%20%s%%20%s%%20%s%%20%s%%20%s&Dec=&Equinox=J2000&ImageSize=%.1f&MaxInt=10&GIF=%d&FITS=%d&Download=1"'%(outfile,imserver,str(rahr),str(ramin),str(rasec),str(decdeg),str(decmin),str(decsec),imsize,gif,fits))
    os.system(command)

# part 1 - read data file and FIRST catalogues
first, data = np.loadtxt(firstfile), np.loadtxt(infile)
fb = first[first[:,3]>50]
if data.shape[1] not in [2,6]:
    print('Must have 2 or 6 columns'); exit()
if data.shape[1]==6:      # convert hmsdms -> deg
    data[:,0]=data[:,0]*15.+data[:,1]*15./60.+data[:,2]*15./3600
    data[:,1]=data[:,3]+data[:,4]/60.+data[:,5]/3600.
    data = data[:,:2]

# part 2/3/4 loop
c = SkyCoord (data, unit=u.deg)        #skycoord object for input file
fcat = SkyCoord (fb[:,:2],unit=u.deg)  #skycoord object for bright FIRST
ccorr = c.match_to_catalog_sky(fcat)   #produce set of matches
for i,imin in enumerate(ccorr[0]):     #iterate over the match indices
    first_download(fb[imin,0],fb[imin,1],'temp.fits')
    try:
        nside = int(np.sqrt(len(data)))+1   #sqrt(n)+1 plots on a side
        plt.subplot(nside,nside,i+1,xticks=[],yticks=[]) # part 4
        a,h = getdata('temp.fits').squeeze(),getheader('temp.fits') # part 3
        noise = np.median(abs(a)) # find noise level or use np.std() around edge
        plt.contour(a,np.array([3.,6.,12.,24.,48.,96.,192.])*noise)
        plt.plot([2.,2.],[2.,2+10./(3600*h['CDELT2'])])   # 10" bar on plot
    except:              # fits file didn't download properly
        pass

plt.savefig('model2.png')
